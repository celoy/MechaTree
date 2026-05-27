"""Native vectorised momentum-wind wind kernel (Step 25b).

The simplest 3-D canopy wind solver that does the job:

1. **Per-branch force** (same kernel as Eloy et al. Nat Commun 2017):
   ``F = ½ ρ C_D D L U² sin² i`` is the magnitude of the normal force
   on a branch at inclination ``i`` to the wind. The streamwise
   component is ``F_D = F · sin i = ½ ρ C_D D L U² sin³ i``.
2. **Per-cell drag**: bin branches by centroid, sum per cell.
3. **Actuator-disk update**: ``U_out = ½ U_in [1 + √(1 - 4 F_D /
   (ρ grid_size² U_in²))]`` — the +-root of the momentum-balance quadratic
   ``U_out² - U_in U_out + F_D/(ρ grid_size²) = 0`` (see
   ``docs/momentum_wind_derivation.md``).
4. **Cross-stream diffusion**: implicit 4-neighbour smoothing in
   (y, z), one knob ``nu_diff``.

Occam's razor cuts that the legacy k-ε reference had and we don't:

- No cylinder subdivision (assume ``max(L) <= grid_size`` — true for unit-length
  MechaTree twigs at the default ``grid_size = 1``).
- No equivalent-obstacle / sub-cell sheltering correction (the
  ``F_D - F_D²/(ρH²U²)`` second-order term).
- No skin-friction ``4/√Re sin^(3/2)`` correction (negligible at the
  high Re MechaTree runs at).
- No lift forces (only the streamwise component matters for U_out).
- No k-ε closure: a single ``nu_diff`` replaces the
  ``I_turb · l_mix · C_μ`` chain. See the derivation doc §5.

Used by :class:`mechatree.wind.momentum_wind.MomentumWindBridge`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MomentumWindResult:
    """Output bundle of :func:`compute_momentum_wind`.

    Per-cell grid fields (``U_out``, ``U_in_grid``, ``F_D_cell``) for
    diagnostics + viz. Per-branch arrays (``U_branch``, ``F_N``,
    ``F_D_branch``, ``F_vec_branch``) for feeding back into the
    tree's pruning (option B in the wind/mechanics discussion).
    """

    U_out: np.ndarray  # (Nz, Ny, Nx) post-canopy wind magnitude
    U_in_grid: np.ndarray  # (Nz, Ny, Nx) inflow into each cell
    F_D_cell: np.ndarray  # (Nz, Ny, Nx) per-cell streamwise drag (sum)
    U_branch: np.ndarray  # (N,) per-branch local wind magnitude (the |U_in| at the branch's cell)
    F_N_branch: np.ndarray  # (N,) per-branch normal-force magnitude = ½ ρ D L U² sin²i
    F_D_branch: np.ndarray  # (N,) per-branch streamwise drag = F_N · sin i
    F_vec_branch: np.ndarray  # (N, 3) per-branch full force vector (Nat Comms F_seg)
    cell_bounds_x: np.ndarray
    cell_bounds_y: np.ndarray
    cell_bounds_z: np.ndarray
    x_centers: np.ndarray
    y_centers: np.ndarray
    z_centers: np.ndarray
    wind_direction: tuple[float, float]

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.U_out.shape  # (Nz, Ny, Nx)


def compute_momentum_wind(
    start: np.ndarray,
    axis: np.ndarray,
    D: np.ndarray,
    L: np.ndarray,
    *,
    cell_bounds_x: np.ndarray,
    cell_bounds_y: np.ndarray,
    cell_bounds_z: np.ndarray,
    grid_size: float,
    U_infty: np.ndarray,
    C_D: float = 1.0,
    nu_diff: float = 0.03,
    diffusion_per_line: bool = True,
    wind_direction: tuple[float, float] = (1.0, 0.0),
) -> MomentumWindResult:
    """Run one forward momentum-wind solve.

    Storm direction is assumed ``+x`` in the input frame; the caller
    is responsible for rotating the canopy if the world-frame storm
    angle is non-zero.

    Parameters
    ----------
    start, axis, D, L
        Per-branch geometry: ``(N, 3)`` start points, ``(N, 3)`` unit
        axis vectors, ``(N,)`` diameters and lengths. Caller must
        ensure ``max(L) <= grid_size`` so centroid-binning is exact (each
        branch belongs to exactly one cell).
    cell_bounds_x / y / z
        1-D arrays of cell boundaries, strictly monotone. Spacing
        must equal ``grid_size``.
    grid_size
        Uniform cell size.
    U_infty
        ``(Nz,)`` per-layer inflow speeds; sets the inlet + top +
        lateral boundary conditions.
    C_D
        Branch drag coefficient (default 1.0).
    nu_diff
        Dimensionless cross-stream diffusion coefficient. Per-cell
        diffusion weight is ``w = nu_diff / (|U| · grid_size)``.
    diffusion_per_line
        ``True`` interleaves the (y, z) diffusion after every x-column
        update. ``False`` defers it to a single global pass.
    wind_direction
        Stored on the result for the caller; the kernel itself solves
        with the storm pointing ``+x``.
    """
    Nx = int(cell_bounds_x.size - 1)
    Ny = int(cell_bounds_y.size - 1)
    Nz = int(cell_bounds_z.size - 1)
    n = int(start.shape[0])

    # ---- Centroid binning ----
    mid = start + 0.5 * L[:, None] * axis
    i_idx = np.clip(np.searchsorted(cell_bounds_x, mid[:, 0], side="right") - 1, 0, Nx - 1)
    j_idx = np.clip(np.searchsorted(cell_bounds_y, mid[:, 1], side="right") - 1, 0, Ny - 1)
    k_idx = np.clip(np.searchsorted(cell_bounds_z, mid[:, 2], side="right") - 1, 0, Nz - 1)

    # Nat Comms vector force kernel (Eloy et al. 2017):
    #   F_seg = ½ ρ U² d L ||t × u||² (n × t)
    # With wind = +x in this frame: cos i = t·u = t_x, and the
    # in-plane direction m = n × t = (u - cos i · t).
    # |m| = sin i, so:
    #   F_seg = ½ ρ U² d L sin²i · m̂        (vector form)
    #         = ½ ρ U² d L · sin i · m       (multiplied form, |m| = sin i)
    # Streamwise: F_D = F_seg · û = ½ ρ U² d L sin³ i
    #
    # Compute the U²-independent factor once; U_loc² scales per cell.
    cos_I = axis[:, 0]  # = t · û since û = +x
    u_perp = np.zeros_like(axis)  # (N, 3), m * sin i
    u_perp[:, 0] = 1.0 - cos_I * cos_I  # u - cos_I · t, with u = (1,0,0)
    u_perp[:, 1] = -cos_I * axis[:, 1]
    u_perp[:, 2] = -cos_I * axis[:, 2]
    # |u_perp| = sin i. Force geometric factor: ½ ρ D L · u_perp · sin i,
    # so that |F| = ½ ρ U² D L sin²i.
    sin_I = np.sqrt(np.clip(1.0 - cos_I * cos_I, 0.0, 1.0))
    half_DL = 0.5 * C_D * D * L  # ρ = 1
    drag_geom = half_DL * (sin_I * sin_I * sin_I)  # for the streamwise scalar

    # Sort by i_idx for fast column slicing.
    order = np.argsort(i_idx, kind="stable")
    i_sorted = i_idx[order]
    j_sorted = j_idx[order]
    k_sorted = k_idx[order]
    drag_geom_sorted = drag_geom[order]
    column_starts = np.searchsorted(i_sorted, np.arange(Nx + 1), side="left")

    # ---- x-march (sequential; vectorised within each (Nz, Ny) slice) ----
    U_out = np.empty((Nz, Ny, Nx))
    U_in_grid = np.empty((Nz, Ny, Nx))
    F_D_grid = np.zeros((Nz, Ny, Nx))
    U_branch_sorted = np.zeros(n)
    F_D_branch_sorted = np.zeros(n)
    U_loc_sq_sorted = np.zeros(n)  # U_loc² for later F_vec scaling
    Nzy = Nz * Ny

    # Pre-allocated scratch buffers, reused every column iteration to
    # avoid per-slice allocation churn (~10-20% win on big grids).
    U_in = np.empty((Nz, Ny))
    U_slice = np.empty((Nz, Ny))
    F_D_cell = np.empty((Nz, Ny))
    neighbour_sum = np.empty((Nz, Ny))
    w_buf = np.empty((Nz, Ny))
    grid_size2 = grid_size * grid_size

    # Inflow at column 0.
    U_in[:] = U_infty[:, None]

    for i in range(Nx):
        if i > 0:
            if diffusion_per_line:
                _diffuse_slice_into(
                    U_out[:, :, i - 1], nu_diff, grid_size, U_in, neighbour_sum, w_buf
                )
            else:
                np.copyto(U_in, U_out[:, :, i - 1])
        U_in_grid[:, :, i] = U_in

        # 2. Per-cell drag via bincount. Skip cleanly if column is
        # empty (no branches in this x-slice).
        s, e = column_starts[i], column_starts[i + 1]
        if e > s:
            j_h = j_sorted[s:e]
            k_h = k_sorted[s:e]
            U_loc = U_in[k_h, j_h]
            U_loc_sq = U_loc * U_loc
            F_D_b = drag_geom_sorted[s:e] * U_loc_sq
            U_branch_sorted[s:e] = U_loc
            F_D_branch_sorted[s:e] = F_D_b
            U_loc_sq_sorted[s:e] = U_loc_sq
            lin = k_h * Ny + j_h
            np.copyto(
                F_D_cell,
                np.bincount(lin, weights=F_D_b, minlength=Nzy).reshape(Nz, Ny),
            )
            F_D_grid[:, :, i] = F_D_cell

            # 3. Actuator-disk update (vectorised; ρ = 1):
            #    U_out = ½ U_in [1 + √(1 - 4 F_D / (grid_size² U_in²))]
            denom = grid_size2 * U_in * U_in + 1e-30
            disc = np.clip(1.0 - 4.0 * F_D_cell / denom, 0.0, 1.0)
            np.multiply(U_in, 0.5 * (1.0 + np.sqrt(disc)), out=U_slice)
        else:
            # Empty column → no drag, slice = inflow. Skip the
            # momentum-wind arithmetic entirely.
            F_D_grid[:, :, i] = 0.0
            np.copyto(U_slice, U_in)

        # 4. Optional cross-stream diffusion (in place into U_slice).
        if diffusion_per_line:
            _diffuse_slice_into(U_slice, nu_diff, grid_size, U_slice, neighbour_sum, w_buf)

        # 5. BCs: top z + lateral y → free stream.
        U_slice[-1, :] = U_infty[-1]
        U_slice[:, 0] = U_infty
        U_slice[:, -1] = U_infty

        U_out[:, :, i] = U_slice

    if not diffusion_per_line:
        tmp = np.empty((Nz, Ny))
        for i in range(1, Nx):
            np.copyto(tmp, U_out[:, :, i])
            _diffuse_slice_into(tmp, nu_diff, grid_size, U_out[:, :, i], neighbour_sum, w_buf)

    x_centers = 0.5 * (cell_bounds_x[:-1] + cell_bounds_x[1:])
    y_centers = 0.5 * (cell_bounds_y[:-1] + cell_bounds_y[1:])
    z_centers = 0.5 * (cell_bounds_z[:-1] + cell_bounds_z[1:])

    # Un-sort per-branch arrays back to caller's input order.
    inverse_order = np.empty(n, dtype=np.int64)
    inverse_order[order] = np.arange(n)
    U_branch = U_branch_sorted[inverse_order]
    F_D_branch = F_D_branch_sorted[inverse_order]
    U_loc_sq = U_loc_sq_sorted[inverse_order]
    # F_N magnitude per branch (Nat Comms |F_seg| = ½ ρ D L U² sin²i).
    F_N_branch = 0.5 * C_D * D * L * U_loc_sq * (sin_I * sin_I)
    # F_seg vector per branch = ½ ρ D L U² sin i · u_perp (with |u_perp| = sin i).
    F_vec_branch = (0.5 * C_D * D * L * sin_I * U_loc_sq)[:, None] * u_perp

    return MomentumWindResult(
        U_out=U_out,
        U_in_grid=U_in_grid,
        F_D_cell=F_D_grid,
        U_branch=U_branch,
        F_N_branch=F_N_branch,
        F_D_branch=F_D_branch,
        F_vec_branch=F_vec_branch,
        cell_bounds_x=cell_bounds_x,
        cell_bounds_y=cell_bounds_y,
        cell_bounds_z=cell_bounds_z,
        x_centers=x_centers,
        y_centers=y_centers,
        z_centers=z_centers,
        wind_direction=wind_direction,
    )


def _diffuse_slice_into(
    U_prev: np.ndarray,
    nu_diff: float,
    grid_size: float,
    out: np.ndarray,
    neighbour_sum: np.ndarray,
    w_buf: np.ndarray,
) -> None:
    """Implicit 4-neighbour edge-clamped cross-stream diffusion, in place.

    Per cell averages with its (y, z) neighbours, weighted by
    ``w = nu_diff / (|U| · grid_size)``. Edges are clamped (cells at the
    boundary take their own value as the off-grid neighbour).

    Writes into ``out`` (which may alias ``U_prev`` — internal
    slicing is index-careful enough). ``neighbour_sum`` and ``w_buf``
    are pre-allocated scratch buffers shape ``(Nz, Ny)``, kept by
    the caller so this hot path doesn't allocate.

    Manual indexing instead of ``np.pad`` — eliminates the per-call
    array allocations that dominated profiling
    (~80% of wall-clock in the previous implementation).
    """
    # Sum of 4 edge-clamped neighbours, computed without np.pad.
    # y-direction: left neighbour at column j-1 (clamped to j=0),
    # right neighbour at j+1 (clamped to j=Ny-1).
    neighbour_sum[:, 0] = U_prev[:, 0] + U_prev[:, 1]
    neighbour_sum[:, -1] = U_prev[:, -2] + U_prev[:, -1]
    neighbour_sum[:, 1:-1] = U_prev[:, :-2] + U_prev[:, 2:]
    # z-direction: bottom neighbour at row k-1 (clamped to k=0),
    # top neighbour at k+1 (clamped to k=Nz-1).
    neighbour_sum[0, :] += U_prev[0, :] + U_prev[1, :]
    neighbour_sum[-1, :] += U_prev[-2, :] + U_prev[-1, :]
    neighbour_sum[1:-1, :] += U_prev[:-2, :] + U_prev[2:, :]

    # w = nu_diff / (|U| · grid_size); ε guard against zero-velocity cells.
    np.abs(U_prev, out=w_buf)
    w_buf *= grid_size
    w_buf += 1e-30
    np.divide(nu_diff, w_buf, out=w_buf)

    # out = (U_prev + w · Σ_neighbours) / (1 + 4w)
    np.multiply(w_buf, neighbour_sum, out=neighbour_sum)
    neighbour_sum += U_prev
    np.multiply(4.0, w_buf, out=w_buf)
    w_buf += 1.0
    np.divide(neighbour_sum, w_buf, out=out)


__all__ = ["MomentumWindResult", "compute_momentum_wind"]
