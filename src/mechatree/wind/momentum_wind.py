"""Momentum-wind CFD bridge (Step 25b; renamed from actuator-disk in 26b).

A canopy-aware ``WindFn`` backed by
:func:`mechatree.wind._momentum_wind_kernel.compute_momentum_wind`.
Builds a structured grid from the forest bounding box, runs the
single-pass momentum-wind solve, and returns the canopy-mean wind in
the world frame.

The 3-D wind field ``U_out(x, y, z)`` is exposed via
``last_result.U_out`` for diagnostics + the notebook colourmap.

The kernel only needs ``(start, axis, D,
L)`` arrays which come straight from :meth:`PyTree.get_branch_data_batch`.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np

from mechatree.wind._momentum_wind_kernel import (
    MomentumWindResult,
    compute_momentum_wind,
)

if TYPE_CHECKING:
    from mechatree._core import PyTree
    from mechatree.forest import Forest


class MomentumWindBridge:
    """Stateful ``WindFn`` backed by the momentum-wind kernel.

    Use the 3-arg ``WindFn`` shape ``(generation, rng, context)``.
    The bridge:

    1. Pools the forest's PyTrees into ``(start, axis, D, L)`` arrays.
    2. Optionally rotates the geometry into the wind frame so the
       kernel (which solves with the storm pointing ``+x``) sees the
       storm direction as ``+x``.
    3. Builds a structured grid from the bounding box with the
       configured cell size ``grid_size`` and ``(pad_x, pad_y, pad_z)``
       padding.
    4. Builds the inflow profile. By default a log-law
       ``U_∞(z) = (u_a / κ) log(max(z, z₀)/z₀)``; if ``U_uniform`` is
       set, a height-independent constant ``U_∞(z) = U_uniform``
       (``u_a / z₀ / κ`` are then ignored).
    5. Runs :func:`compute_momentum_wind` once.
    6. Reduces ``U_out`` to a canopy-mean (mean over cells with
       branches in them) and rotates back to the world frame.

    ``last_result`` carries the most recent :class:`MomentumWindResult`
    for diagnostics. ``last_wind_direction`` stores the storm direction
    used.
    """

    def __init__(
        self,
        *,
        grid_size: float = 2.0,
        nu_diff: float = 0.03,
        pad_x: float = 12.0,
        pad_y: float = 2.0,
        pad_z: float = 3.0,
        ua: float = 0.4,
        z0: float = 0.1,
        kappa: float = 0.41,
        U_uniform: float | None = None,
        C_D: float = 1.0,
        diffusion_per_line: bool = True,
        angle_sampler: Callable[[np.random.Generator, int], np.ndarray] | None = None,
        amplitude_sampler: Callable[[np.random.Generator, int], np.ndarray] | None = None,
    ) -> None:
        self.grid_size = float(grid_size)
        self.nu_diff = float(nu_diff)
        self.pad_x = float(pad_x)
        self.pad_y = float(pad_y)
        self.pad_z = float(pad_z)
        self.ua = float(ua)
        self.z0 = float(z0)
        self.kappa = float(kappa)
        self.U_uniform = None if U_uniform is None else float(U_uniform)
        self.C_D = float(C_D)
        self.diffusion_per_line = bool(diffusion_per_line)
        # Step 26a: per-branch forces are now the *only* momentum behaviour.
        # The bridge always writes each branch's CFD force + local wind back
        # onto the trees so the prune loop reads them via
        # ``prune_with_stored_forces`` (a 3-D solve collapsed to one canopy-
        # mean scalar would just be a constant wind). ``writes_segment_forces``
        # is the capability flag the prune loops dispatch on via ``getattr``.
        self.writes_segment_forces = True
        self.angle_sampler = angle_sampler
        self.amplitude_sampler = amplitude_sampler
        self.last_result: MomentumWindResult | None = None
        self.last_wind_direction: tuple[float, float] = (1.0, 0.0)
        # Step 26a: the storm (θ, amplitude) is sampled once per generation
        # and held across the fixed-point pruning iterations. The Step-24 loop
        # calls this bridge once per inner iteration with the same
        # ``generation``; re-sampling each time would make a single gen's
        # "storm" a sequence of *different* storms (inflating the iteration
        # count + biasing big-storm pruning). Cache keyed on ``generation``.
        self._storm_gen: int | None = None
        self._storm_theta: float = 0.0
        self._storm_amp: float = 1.0

    @staticmethod
    def _trees_of(context: PyTree | Forest) -> list:
        # Duck-typed: anything with a `.trees` attribute is forest-like;
        # otherwise the context itself is the single PyTree.
        return context.trees if hasattr(context, "trees") else [context]

    @staticmethod
    def _pool(trees: list):
        """Concatenate every tree's ``(start, axis, D, L)`` plus a per-tree
        branch-count list (for splitting the per-branch forces back)."""
        starts, axes, Ds, Ls, counts = [], [], [], [], []
        for t in trees:
            s, a, d, ell = t.get_branch_data_batch()
            starts.append(s)
            axes.append(a)
            Ds.append(d)
            Ls.append(ell)
            counts.append(s.shape[0])
        return (
            np.concatenate(starts),
            np.concatenate(axes),
            np.concatenate(Ds),
            np.concatenate(Ls),
            counts,
        )

    def _solve_and_store(
        self,
        trees: list,
        counts: list,
        start: np.ndarray,
        axis: np.ndarray,
        D: np.ndarray,
        L: np.ndarray,
        *,
        theta: float,
        U_uniform: float | None,
        amp: float,
    ) -> float:
        """Solve the momentum field for one direction ``theta`` and write the
        per-branch screened force + local wind back onto each tree. Shared by
        the pruning call (storm direction + configured inflow) and the Step-26c
        sensing sweep (explicit angle + uniform inlet). Returns the canopy-mean
        wind magnitude (the ε convergence thermometer)."""
        cx, cy = math.cos(theta), math.sin(theta)
        self.last_wind_direction = (cx, cy)

        # Rotate horizontal position + axis into the wind frame (storm → +x).
        if theta != 0.0:
            ct, st = math.cos(-theta), math.sin(-theta)
            x_r = ct * start[:, 0] - st * start[:, 1]
            y_r = st * start[:, 0] + ct * start[:, 1]
            nx_r = ct * axis[:, 0] - st * axis[:, 1]
            ny_r = st * axis[:, 0] + ct * axis[:, 1]
            start = np.column_stack([x_r, y_r, start[:, 2]])
            axis = np.column_stack([nx_r, ny_r, axis[:, 2]])

        grid_size = self.grid_size
        x_lo = float(start[:, 0].min()) - self.pad_x
        x_hi = float(start[:, 0].max()) + self.pad_x
        y_lo = float(start[:, 1].min()) - self.pad_y
        y_hi = float(start[:, 1].max()) + self.pad_y
        z_hi = float((start[:, 2] + L * axis[:, 2]).max()) + self.pad_z

        cell_bounds_x = np.arange(x_lo, x_hi + grid_size, grid_size)
        cell_bounds_y = np.arange(y_lo, y_hi + grid_size, grid_size)
        cell_bounds_z = np.arange(0.0, z_hi + grid_size, grid_size)
        z_centers = 0.5 * (cell_bounds_z[:-1] + cell_bounds_z[1:])

        # Inflow profile: uniform if ``U_uniform`` is set, else the log-law.
        if U_uniform is not None:
            U_infty = np.full(z_centers.size, U_uniform)
        else:
            U_infty = (self.ua / self.kappa) * np.log(np.maximum(z_centers, self.z0) / self.z0)
        if amp != 1.0:
            U_infty = U_infty * amp

        result = compute_momentum_wind(
            start,
            axis,
            D,
            L,
            cell_bounds_x=cell_bounds_x,
            cell_bounds_y=cell_bounds_y,
            cell_bounds_z=cell_bounds_z,
            grid_size=grid_size,
            U_infty=U_infty,
            C_D=self.C_D,
            nu_diff=self.nu_diff,
            diffusion_per_line=self.diffusion_per_line,
            wind_direction=(cx, cy),
        )
        self.last_result = result

        # Plumb the per-branch force + local wind back onto each tree (rotate
        # F_vec from the +x solve frame to the world frame by +θ; the local
        # wind is along the storm direction at per-branch magnitude U_branch).
        Fx, Fy = result.F_vec_branch[:, 0], result.F_vec_branch[:, 1]
        F_world = np.column_stack([cx * Fx - cy * Fy, cy * Fx + cx * Fy, result.F_vec_branch[:, 2]])
        w_world = result.U_branch[:, None] * np.array([cx, cy, 0.0])
        split_at = np.cumsum(counts)[:-1]
        F_per_tree = np.split(F_world, split_at)
        w_per_tree = np.split(w_world, split_at)
        for t, Ft, wt in zip(trees, F_per_tree, w_per_tree, strict=True):
            t.set_segment_forces_batch(Ft)
            t.set_segment_winds_batch(wt)

        # Canopy-mean over occupied cells (ε convergence thermometer only —
        # pruning reads the per-branch forces above, not this scalar).
        Nx = len(cell_bounds_x) - 1
        Ny = len(cell_bounds_y) - 1
        Nz = len(cell_bounds_z) - 1
        mid = start + 0.5 * L[:, None] * axis
        i_idx = np.clip(np.searchsorted(cell_bounds_x, mid[:, 0], side="right") - 1, 0, Nx - 1)
        j_idx = np.clip(np.searchsorted(cell_bounds_y, mid[:, 1], side="right") - 1, 0, Ny - 1)
        k_idx = np.clip(np.searchsorted(cell_bounds_z, mid[:, 2], side="right") - 1, 0, Nz - 1)
        return float(np.mean(result.U_out[k_idx, j_idx, i_idx]))

    def __call__(
        self,
        generation: int,
        rng: np.random.Generator,
        context: PyTree | Forest,
    ) -> tuple[float, float, float]:
        trees = self._trees_of(context)
        if not trees or all(t.get_number_of_branches() == 0 for t in trees):
            return (0.0, 0.0, 0.0)
        start, axis, D, L, counts = self._pool(trees)

        # Sample the storm (θ, amplitude) once per generation (Step 26a) and
        # reuse it across the fixed-point loop's inner iterations.
        if generation != self._storm_gen:
            self._storm_gen = generation
            self._storm_theta = (
                float(self.angle_sampler(rng, 1)[0]) if self.angle_sampler is not None else 0.0
            )
            self._storm_amp = (
                float(self.amplitude_sampler(rng, 1)[0])
                if self.amplitude_sampler is not None
                else 1.0
            )
        u_mag = self._solve_and_store(
            trees,
            counts,
            start,
            axis,
            D,
            L,
            theta=self._storm_theta,
            U_uniform=self.U_uniform,
            amp=self._storm_amp,
        )
        cx, cy = self.last_wind_direction
        return (u_mag * cx, u_mag * cy, 0.0)

    def sensing_angles(self, rng: np.random.Generator, n: int) -> list[float]:
        """``n`` sensing directions for the Step-26c stress sweep, drawn from
        the storm angle distribution (uniform ``[0, 2π)`` when no ``angle_cdf``
        is configured)."""
        if self.angle_sampler is not None:
            return [float(a) for a in self.angle_sampler(rng, n)]
        return [float(a) for a in rng.uniform(0.0, 2.0 * math.pi, n)]

    def sense(self, context: PyTree | Forest, theta: float) -> None:
        """Step 26c: solve the momentum field for sensing direction ``theta``
        at a **uniform inlet of U = 1** (the wind scale) and write each
        branch's local screened force + wind onto the trees. The kernel still
        screens U = 1 down through the canopy, so a sheltered branch is sensed
        against the weaker wind it actually feels. No storm sampling — ``theta``
        is explicit, so this never touches the per-generation storm cache."""
        trees = self._trees_of(context)
        if not trees or all(t.get_number_of_branches() == 0 for t in trees):
            return
        start, axis, D, L, counts = self._pool(trees)
        self._solve_and_store(trees, counts, start, axis, D, L, theta=theta, U_uniform=1.0, amp=1.0)


def make_momentum_wind_fn(
    *,
    grid_size: float = 2.0,
    nu_diff: float = 0.03,
    pad_x: float = 12.0,
    pad_y: float = 2.0,
    pad_z: float = 3.0,
    ua: float = 0.4,
    z0: float = 0.1,
    kappa: float = 0.41,
    U_uniform: float | None = None,
    C_D: float = 1.0,
    diffusion_per_line: bool = True,
    angle_sampler: Callable[[np.random.Generator, int], np.ndarray] | None = None,
    amplitude_sampler: Callable[[np.random.Generator, int], np.ndarray] | None = None,
) -> MomentumWindBridge:
    """Construct an :class:`MomentumWindBridge`.

    Inflow is a log-law (``ua / z0 / kappa``) by default; pass
    ``U_uniform`` for a height-independent constant inflow instead.

    The bridge always writes each branch's CFD force + local wind back onto
    the trees so the prune loop scores branches against their own local wind
    (Step 26a — per-branch is the only momentum behaviour; the canopy-mean is
    just the convergence thermometer).
    """
    return MomentumWindBridge(
        grid_size=grid_size,
        nu_diff=nu_diff,
        pad_x=pad_x,
        pad_y=pad_y,
        pad_z=pad_z,
        ua=ua,
        z0=z0,
        kappa=kappa,
        U_uniform=U_uniform,
        C_D=C_D,
        diffusion_per_line=diffusion_per_line,
        angle_sampler=angle_sampler,
        amplitude_sampler=amplitude_sampler,
    )


__all__ = ["MomentumWindBridge", "make_momentum_wind_fn"]
