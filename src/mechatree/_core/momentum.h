/*
 * Momentum-wind — GIL-free C++ port of the per-solve column march (Step 26e).
 *
 * Mirrors the pure-NumPy reference in
 * src/mechatree/wind/_momentum_wind_kernel.py (compute_momentum_wind). The
 * reference stays the readable spec + the equivalence-test oracle; this kernel
 * is the hot path. Numerically equivalent to atol 1e-10 (FMA contraction is
 * the only source of last-bit drift; the per-cell drag is accumulated in the
 * same original-within-column order as np.bincount so summation order matches).
 *
 * The Cython binding (mechatree._core._core.momentum_solve_kernel) declares
 * this `nogil`, so the independent n_sensing_angles solves can run on a thread
 * pool without contending for the GIL — the lever Step 26d measured we needed
 * (a thread pool over the NumPy kernel ran at 0.67x because the Python column
 * march held the GIL).
 *
 * One call = one forward solve. Storm direction is assumed +x in the input
 * frame; the Python bridge rotates the canopy before calling.
 */

#ifndef MECHATREE_MOMENTUM_H_
#define MECHATREE_MOMENTUM_H_

#include <cstddef>

// Run one forward momentum-wind solve.
//
// Per-branch geometry (row-major):
//   start  — (n, 3) branch base points
//   axis   — (n, 3) unit axis vectors
//   D, L   — (n,)   diameters, lengths
// Grid (uniform spacing == grid_size, strictly monotone bounds):
//   cell_bounds_x/y/z — (Nx+1)/(Ny+1)/(Nz+1) cell boundaries
//   U_infty           — (Nz,) per-layer inflow speed (inlet + top + lateral BC)
// Params:
//   grid_size, C_D, nu_diff, diffusion_per_line (0/1)
// Outputs (caller-allocated; grids row-major (Nz, Ny, Nx)):
//   U_out         — (Nz*Ny*Nx) post-canopy wind magnitude
//   U_in_grid     — (Nz*Ny*Nx) inflow into each cell, OR nullptr to skip
//   F_D_cell_grid — (Nz*Ny*Nx) per-cell streamwise drag (sum), OR nullptr to skip
//   U_branch      — (n,)   per-branch local wind magnitude
//   F_N_branch    — (n,)   per-branch normal-force magnitude
//   F_D_branch    — (n,)   per-branch streamwise drag
//   F_vec_branch  — (n, 3) per-branch full force vector (Nat Comms F_seg),
//                   horizontal components rotated by (cos_theta, sin_theta)
//   canopy_mean_out — (1,) mean of U_out over the cells containing branches,
//                     OR nullptr to skip (Step 26f — replaces a Python searchsorted)
//
// cos_theta / sin_theta rotate the per-branch F_vec horizontals to the world
// frame ((1, 0) leaves them in the solve frame — the storm-frame oracle).
//
// Pre-condition: max(L) <= grid_size (centroid binning is exact).
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
                    double* canopy_mean_out);

// Consolidated sensing-path entry (Step 26f): everything in C++ so the bridge
// just pools the *unrotated* geometry and passes the storm angle. Rotates the
// canopy by -theta for binning, builds the grid internally from the rotated
// bounding box (replicating np.arange(lo, hi+grid_size, grid_size)), builds the
// inflow profile, runs the solve, and rotates the per-branch force back to the
// world frame by +theta.
//
//   start, axis, D, L, n — UNROTATED pooled per-branch geometry
//   theta                — storm direction (world frame)
//   grid_size, pad_x/y/z — grid cell size + bounding-box padding
//   U_uniform            — uniform inflow speed if >= 0; else use the log-law
//   ua, z0, kappa        — log-law params (used when U_uniform < 0)
//   amp                  — storm amplitude multiplier on the inflow
//   C_D, nu_diff, diffusion_per_line — kernel params
// Outputs (caller-allocated):
//   F_world — (n, 3) per-branch force in the world frame
//   w_world — (n, 3) per-branch local wind = U_branch * (cos theta, sin theta, 0)
//
// No grid arrays are produced (built as internal scratch); the sensing sweep
// only needs the per-branch fields.
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
                          double* w_world);

#endif  // MECHATREE_MOMENTUM_H_
