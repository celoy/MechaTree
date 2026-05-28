#include "momentum.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

inline double clip01(double v) {
    if (v < 0.0) return 0.0;
    if (v > 1.0) return 1.0;
    return v;
}

// searchsorted(bounds, val, side="right") - 1, clamped to [0, n_cells-1].
// Matches the NumPy binning in _momentum_wind_kernel.compute_momentum_wind.
inline int bin_index(const double* bounds, std::size_t nb, double val, int n_cells) {
    // upper_bound returns first element strictly greater than val; its offset
    // equals searchsorted(side="right").
    const double* hi = std::upper_bound(bounds, bounds + nb, val);
    int idx = static_cast<int>(hi - bounds) - 1;
    if (idx < 0) idx = 0;
    if (idx > n_cells - 1) idx = n_cells - 1;
    return idx;
}

// Implicit 4-neighbour edge-clamped cross-stream diffusion over one (Nz, Ny)
// slice. Mirrors _diffuse_slice_into. `out` may alias `U_prev`: nb_sum and
// w_buf are fully computed (reading U_prev) before any out write, and the
// final write out[c] reads only U_prev[c] (its own cell), so aliasing is safe.
void diffuse_slice(const double* U_prev,
                   double* out,
                   double* nb_sum,
                   double* w_buf,
                   int Nz,
                   int Ny,
                   double nu_diff,
                   double grid_size) {
    for (int k = 0; k < Nz; ++k) {
        const int km = (k > 0) ? k - 1 : 0;
        const int kp = (k < Nz - 1) ? k + 1 : Nz - 1;
        for (int j = 0; j < Ny; ++j) {
            const int jm = (j > 0) ? j - 1 : 0;
            const int jp = (j < Ny - 1) ? j + 1 : Ny - 1;
            const double y_part = U_prev[k * Ny + jm] + U_prev[k * Ny + jp];
            const double z_part = U_prev[km * Ny + j] + U_prev[kp * Ny + j];
            nb_sum[k * Ny + j] = y_part + z_part;
            // w = nu_diff / (|U_prev| * grid_size + 1e-30)
            w_buf[k * Ny + j] = nu_diff / (std::fabs(U_prev[k * Ny + j]) * grid_size + 1e-30);
        }
    }
    const int Nzy = Nz * Ny;
    for (int c = 0; c < Nzy; ++c) {
        const double w = w_buf[c];
        out[c] = (U_prev[c] + w * nb_sum[c]) / (1.0 + 4.0 * w);
    }
}

}  // namespace

void momentum_solve(const double* start,
                    const double* axis,
                    const double* D,
                    const double* L,
                    std::size_t n,
                    const double* cell_bounds_x,
                    std::size_t nbx,
                    const double* cell_bounds_y,
                    std::size_t nby,
                    const double* cell_bounds_z,
                    std::size_t nbz,
                    double grid_size,
                    const double* U_infty,
                    double C_D,
                    double nu_diff,
                    int diffusion_per_line,
                    double* U_out,
                    double* U_in_grid,
                    double* F_D_cell_grid,
                    double* U_branch,
                    double* F_N_branch,
                    double* F_D_branch,
                    double* F_vec_branch,
                    double cos_theta,
                    double sin_theta,
                    double* canopy_mean_out) {
    const int Nx = static_cast<int>(nbx) - 1;
    const int Ny = static_cast<int>(nby) - 1;
    const int Nz = static_cast<int>(nbz) - 1;
    const int Nzy = Nz * Ny;
    const std::size_t NxNzy = static_cast<std::size_t>(Nx) * static_cast<std::size_t>(Nzy);
    const double grid_size2 = grid_size * grid_size;

    // ---- Per-branch geometry + binning ----
    std::vector<int> i_idx(n), j_idx(n), k_idx(n);
    std::vector<double> cos_I(n), sin_I(n), drag_geom(n), U_loc_sq(n, 0.0);
    for (std::size_t b = 0; b < n; ++b) {
        const double sx = start[3 * b + 0];
        const double sy = start[3 * b + 1];
        const double sz = start[3 * b + 2];
        const double ax = axis[3 * b + 0];
        const double ay = axis[3 * b + 1];
        const double az = axis[3 * b + 2];
        const double half_len = 0.5 * L[b];
        const double mx = sx + half_len * ax;
        const double my = sy + half_len * ay;
        const double mz = sz + half_len * az;
        i_idx[b] = bin_index(cell_bounds_x, nbx, mx, Nx);
        j_idx[b] = bin_index(cell_bounds_y, nby, my, Ny);
        k_idx[b] = bin_index(cell_bounds_z, nbz, mz, Nz);

        const double c = ax;  // cos i = t · û with û = +x
        cos_I[b] = c;
        const double s = std::sqrt(clip01(1.0 - c * c));
        sin_I[b] = s;
        const double half_DL = 0.5 * C_D * D[b] * L[b];  // ρ = 1
        drag_geom[b] = half_DL * (s * s * s);  // streamwise scalar factor
    }

    // ---- Stable grouping by x-column (counting sort, original order kept) ----
    std::vector<int> col_start(Nx + 1, 0);
    for (std::size_t b = 0; b < n; ++b) col_start[i_idx[b] + 1]++;
    for (int i = 0; i < Nx; ++i) col_start[i + 1] += col_start[i];
    std::vector<int> col_order(n);
    {
        std::vector<int> cursor(col_start.begin(), col_start.end());
        for (std::size_t b = 0; b < n; ++b) col_order[cursor[i_idx[b]]++] = static_cast<int>(b);
    }

    // ---- Scratch slices (size Nzy), reused per column ----
    std::vector<double> U_in(Nzy), U_slice(Nzy), F_D_cell(Nzy);
    std::vector<double> nb_sum(Nzy), w_buf(Nzy), prev_slice(Nzy);

    // Inflow at column 0: U_in[k][j] = U_infty[k].
    for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j) U_in[k * Ny + j] = U_infty[k];

    for (int i = 0; i < Nx; ++i) {
        if (i > 0) {
            // Gather previous column U_out[:, :, i-1] into prev_slice.
            for (int k = 0; k < Nz; ++k)
                for (int j = 0; j < Ny; ++j)
                    prev_slice[k * Ny + j] = U_out[(static_cast<std::size_t>(k) * Ny + j) * Nx + (i - 1)];
            if (diffusion_per_line) {
                diffuse_slice(prev_slice.data(), U_in.data(), nb_sum.data(), w_buf.data(), Nz, Ny,
                              nu_diff, grid_size);
            } else {
                std::copy(prev_slice.begin(), prev_slice.end(), U_in.begin());
            }
        }
        // Scatter U_in into U_in_grid[:, :, i] (skipped when not requested).
        if (U_in_grid != nullptr)
            for (int k = 0; k < Nz; ++k)
                for (int j = 0; j < Ny; ++j)
                    U_in_grid[(static_cast<std::size_t>(k) * Ny + j) * Nx + i] = U_in[k * Ny + j];

        std::fill(F_D_cell.begin(), F_D_cell.end(), 0.0);
        const int s = col_start[i];
        const int e = col_start[i + 1];
        if (e > s) {
            // Per-cell drag (original order within column → matches np.bincount).
            for (int idx = s; idx < e; ++idx) {
                const int b = col_order[idx];
                const int cell = k_idx[b] * Ny + j_idx[b];
                const double U_loc = U_in[cell];
                const double usq = U_loc * U_loc;
                const double F_D_b = drag_geom[b] * usq;
                F_D_cell[cell] += F_D_b;
                U_branch[b] = U_loc;
                F_D_branch[b] = F_D_b;
                U_loc_sq[b] = usq;
            }
            // Actuator-disk update (vectorised over the slice in NumPy).
            for (int c = 0; c < Nzy; ++c) {
                const double uin = U_in[c];
                const double denom = grid_size2 * uin * uin + 1e-30;
                const double disc = clip01(1.0 - 4.0 * F_D_cell[c] / denom);
                U_slice[c] = uin * (0.5 * (1.0 + std::sqrt(disc)));
            }
        } else {
            // Empty column → slice = inflow.
            std::copy(U_in.begin(), U_in.end(), U_slice.begin());
        }
        // Scatter F_D_cell into F_D_cell_grid[:, :, i] (skipped when not requested).
        if (F_D_cell_grid != nullptr)
            for (int k = 0; k < Nz; ++k)
                for (int j = 0; j < Ny; ++j)
                    F_D_cell_grid[(static_cast<std::size_t>(k) * Ny + j) * Nx + i] = F_D_cell[k * Ny + j];

        // Per-line cross-stream diffusion (in place; aliasing-safe).
        if (diffusion_per_line) {
            diffuse_slice(U_slice.data(), U_slice.data(), nb_sum.data(), w_buf.data(), Nz, Ny,
                          nu_diff, grid_size);
        }

        // BCs: top z + lateral y → free stream.
        for (int j = 0; j < Ny; ++j) U_slice[(Nz - 1) * Ny + j] = U_infty[Nz - 1];
        for (int k = 0; k < Nz; ++k) {
            U_slice[k * Ny + 0] = U_infty[k];
            U_slice[k * Ny + (Ny - 1)] = U_infty[k];
        }

        // Scatter U_slice into U_out[:, :, i].
        for (int k = 0; k < Nz; ++k)
            for (int j = 0; j < Ny; ++j)
                U_out[(static_cast<std::size_t>(k) * Ny + j) * Nx + i] = U_slice[k * Ny + j];
    }
    (void)NxNzy;

    if (!diffusion_per_line) {
        // Single global diffusion pass over columns 1..Nx-1.
        for (int i = 1; i < Nx; ++i) {
            for (int k = 0; k < Nz; ++k)
                for (int j = 0; j < Ny; ++j)
                    prev_slice[k * Ny + j] = U_out[(static_cast<std::size_t>(k) * Ny + j) * Nx + i];
            diffuse_slice(prev_slice.data(), U_slice.data(), nb_sum.data(), w_buf.data(), Nz, Ny,
                          nu_diff, grid_size);
            for (int k = 0; k < Nz; ++k)
                for (int j = 0; j < Ny; ++j)
                    U_out[(static_cast<std::size_t>(k) * Ny + j) * Nx + i] = U_slice[k * Ny + j];
        }
    }

    // ---- Per-branch force outputs (F_vec horizontals rotated to world by
    // (cos_theta, sin_theta); (1, 0) leaves them in the solve frame). ----
    for (std::size_t b = 0; b < n; ++b) {
        const double c = cos_I[b];
        const double s = sin_I[b];
        const double ay = axis[3 * b + 1];
        const double az = axis[3 * b + 2];
        const double up0 = 1.0 - c * c;  // = sin²i (matches NumPy u_perp[:,0])
        const double up1 = -c * ay;
        const double up2 = -c * az;
        const double half_CDDL = 0.5 * C_D * D[b] * L[b];
        F_N_branch[b] = half_CDDL * U_loc_sq[b] * (s * s);
        const double coef = half_CDDL * s * U_loc_sq[b];
        const double fx = coef * up0;  // solve-frame force horizontals
        const double fy = coef * up1;
        F_vec_branch[3 * b + 0] = cos_theta * fx - sin_theta * fy;
        F_vec_branch[3 * b + 1] = sin_theta * fx + cos_theta * fy;
        F_vec_branch[3 * b + 2] = coef * up2;
    }

    // ---- Canopy mean of U_out over the cells that contain branches (the ε
    // convergence thermometer; the march already binned every branch, so this
    // replaces the Python searchsorted). ----
    if (canopy_mean_out != nullptr) {
        double acc = 0.0;
        for (std::size_t b = 0; b < n; ++b) {
            const std::size_t cell =
                (static_cast<std::size_t>(k_idx[b]) * Ny + j_idx[b]) * Nx + i_idx[b];
            acc += U_out[cell];
        }
        canopy_mean_out[0] = (n > 0) ? acc / static_cast<double>(n) : 0.0;
    }
}

// ---------------------------------------------------------------------------
// Step 26f: consolidated sensing entry — rotation + grid build + inflow + solve
// + world-frame force, all in C++. See momentum.h.
// ---------------------------------------------------------------------------

namespace {

// Replicate np.arange(lo, stop, step): length = ceil((stop - lo) / step) values
// lo + i*step. The grid build in MomentumWindBridge uses stop = hi + grid_size.
std::vector<double> arange(double lo, double stop, double step) {
    double raw = std::ceil((stop - lo) / step);
    int len = (raw <= 0.0) ? 0 : static_cast<int>(raw);
    std::vector<double> out(static_cast<std::size_t>(len));
    for (int i = 0; i < len; ++i) out[static_cast<std::size_t>(i)] = lo + i * step;
    return out;
}

}  // namespace

void momentum_solve_world(const double* start,
                          const double* axis,
                          const double* D,
                          const double* L,
                          std::size_t n,
                          double theta,
                          double grid_size,
                          double pad_x,
                          double pad_y,
                          double pad_z,
                          double U_uniform,
                          double ua,
                          double z0,
                          double kappa,
                          double amp,
                          double C_D,
                          double nu_diff,
                          int diffusion_per_line,
                          double* F_world,
                          double* w_world) {
    if (n == 0) return;

    const double cos_t = std::cos(theta);
    const double sin_t = std::sin(theta);
    // Rotate horizontals into the wind frame (storm → +x): matches the Python
    // bridge's ct = cos(-theta), st = sin(-theta); x_r = ct*x - st*y, etc.
    const double ct = cos_t;   // cos(-theta)
    const double st = -sin_t;  // sin(-theta)

    std::vector<double> rs(3 * n), ra(3 * n);
    double x_lo = 0.0, x_hi = 0.0, y_lo = 0.0, y_hi = 0.0, z_hi = 0.0;
    for (std::size_t b = 0; b < n; ++b) {
        const double sx = start[3 * b + 0];
        const double sy = start[3 * b + 1];
        const double sz = start[3 * b + 2];
        const double ax = axis[3 * b + 0];
        const double ay = axis[3 * b + 1];
        const double az = axis[3 * b + 2];
        const double rx = ct * sx - st * sy;
        const double ry = st * sx + ct * sy;
        const double rax = ct * ax - st * ay;
        const double ray = st * ax + ct * ay;
        rs[3 * b + 0] = rx;
        rs[3 * b + 1] = ry;
        rs[3 * b + 2] = sz;
        ra[3 * b + 0] = rax;
        ra[3 * b + 1] = ray;
        ra[3 * b + 2] = az;
        const double top = sz + L[b] * az;
        if (b == 0) {
            x_lo = x_hi = rx;
            y_lo = y_hi = ry;
            z_hi = top;
        } else {
            if (rx < x_lo) x_lo = rx;
            if (rx > x_hi) x_hi = rx;
            if (ry < y_lo) y_lo = ry;
            if (ry > y_hi) y_hi = ry;
            if (top > z_hi) z_hi = top;
        }
    }
    x_lo -= pad_x;
    x_hi += pad_x;
    y_lo -= pad_y;
    y_hi += pad_y;
    z_hi += pad_z;

    std::vector<double> cbx = arange(x_lo, x_hi + grid_size, grid_size);
    std::vector<double> cby = arange(y_lo, y_hi + grid_size, grid_size);
    std::vector<double> cbz = arange(0.0, z_hi + grid_size, grid_size);
    const int Nz = static_cast<int>(cbz.size()) - 1;

    // Inflow profile: uniform if U_uniform >= 0, else the log-law; scaled by amp.
    std::vector<double> U_infty(static_cast<std::size_t>(Nz));
    for (int k = 0; k < Nz; ++k) {
        const double zc = 0.5 * (cbz[static_cast<std::size_t>(k)] + cbz[static_cast<std::size_t>(k) + 1]);
        double u;
        if (U_uniform >= 0.0) {
            u = U_uniform;
        } else {
            const double zz = (zc > z0) ? zc : z0;
            u = (ua / kappa) * std::log(zz / z0);
        }
        U_infty[static_cast<std::size_t>(k)] = (amp != 1.0) ? u * amp : u;
    }

    const int Ny = static_cast<int>(cby.size()) - 1;
    const int Nx = static_cast<int>(cbx.size()) - 1;
    const std::size_t grid_cells =
        static_cast<std::size_t>(Nx) * static_cast<std::size_t>(Ny) * static_cast<std::size_t>(Nz);

    // Scratch: U_out (the march needs it) + per-branch buffers. U_in_grid /
    // F_D_cell_grid skipped (nullptr) — sensing doesn't need the grids.
    std::vector<double> U_out(grid_cells);
    std::vector<double> U_branch(n), F_N(n), F_D(n);

    // F_world receives the per-branch force already rotated to the world frame.
    momentum_solve(rs.data(), ra.data(), D, L, n, cbx.data(), cbx.size(), cby.data(), cby.size(),
                   cbz.data(), cbz.size(), grid_size, U_infty.data(), C_D, nu_diff,
                   diffusion_per_line, U_out.data(), nullptr, nullptr, U_branch.data(), F_N.data(),
                   F_D.data(), F_world, cos_t, sin_t, nullptr);

    // w_world = U_branch * (cos theta, sin theta, 0).
    for (std::size_t b = 0; b < n; ++b) {
        w_world[3 * b + 0] = U_branch[b] * cos_t;
        w_world[3 * b + 1] = U_branch[b] * sin_t;
        w_world[3 * b + 2] = 0.0;
    }
}
