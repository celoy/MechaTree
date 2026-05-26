"""Native bulk-thinning wind model for MechaTree (Step 25).

A self-contained, PyTree-native port of DendroFlow's
``BulkThinningBranchWindModel`` ([DendroFlow.wind.branch][df]). Removes
the cross-repo dependency for the common case (one wind model in the
canopy, no advanced k-ε features), and adds **rotation of the canopy
to face the storm direction** so a tunable storm-angle distribution
(:mod:`mechatree.wind.distributions`) actually changes which side of
the tree the wind hits — the original DendroFlow bridge always solves
in the world ``+x`` frame and ignores the storm angle.

Mirrors the light-interception pattern: just as
:func:`mechatree.light.interception.intercept` rotates every leaf into
the sun frame before binning shadows, this module conceptually rotates
every branch axis into a wind frame where ``+x`` is downwind. In
practice the rotation collapses to a single dot product with the wind
direction (``cos I = axis_horizontal · wind_direction``), so we don't
materialise rotated coordinates.

Algorithm
---------
For each branch:

1. Per-branch centre height ``z_c = start.z + 0.5 * L * axis.z`` (invariant
   under horizontal rotation — the only rotation the wind direction
   induces).
2. Per-branch projected area perpendicular to wind:
   ``area = D * L * |sin I|``, where ``cos I = axis_x * wind.x + axis_y * wind.y``
   is the cosine between the branch axis and the (unit) wind direction.
3. Histogram ``area`` onto vertical z-layers (cells of height ``H``,
   centres ``z_centers``).
4. Actuator-disk thinning per z-layer:
   ``U_canopy(z) = 0.5 * U_inf(z) * (1 + sqrt(max(0, 1 - 4 * F_drag(z) / (H^2 * U_inf(z)^2))))``,
   matching the v8 formula at
   [dendroflow.wind.branch.BulkThinningBranchWindModel.compute][df].
5. Look up per-branch ``U_canopy`` at its z; return both per-branch
   wind magnitudes AND the canopy-mean wind vector (rotated back into
   the world frame as ``(Ū cos θ, Ū sin θ, 0)``).

The model defaults to a **uniform free stream** ``U_infty(z) = 1`` over
``z ∈ [0, 50]`` with ``H = 0.5``: the cheapest possible "wind is there"
setting, no boundary layer, no opt-in needed. Realistic boundary-layer
profiles (log-law, power-law) are configured via the YAML ``wind:``
block.

[df]: ../DendroFlow/src/dendroflow/wind/branch.py
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from mechatree._core import PyTree
    from mechatree.forest import Forest


@dataclass(frozen=True)
class BulkThinningResult:
    """Per-branch and per-layer outputs of one ``compute`` call.

    Captures what we need both for the next pruning sweep
    (``canopy_mean``, the single wind vector handed to MechaTree's
    ``WindFn``) and for diagnostics + the 3D wind visualisation
    (``U_branch``).
    """

    U_branch: np.ndarray  # (n_branches,) — wind magnitude at each branch
    z_branch: np.ndarray  # (n_branches,) — branch centre height
    branch_offsets: np.ndarray  # (n_trees + 1,) — slice boundaries by tree
    U_canopy: np.ndarray  # (Nz,) — wind magnitude per z-layer (rotated frame)
    z_centers: np.ndarray  # (Nz,) — z-layer centres
    wind_direction: tuple[float, float]  # (cos θ, sin θ) — world-frame direction

    @property
    def n_branches(self) -> int:
        return int(self.U_branch.shape[0])

    @property
    def canopy_mean(self) -> tuple[float, float, float]:
        """Mean horizontal wind vector ``(wx, wy, 0)`` across all branches.

        Magnitude is the mean of ``U_branch``; direction is the storm
        direction passed into ``compute``. The vertical component is
        zero by construction (this is a horizontal canopy-mean model).
        """
        if self.U_branch.size == 0:
            return (0.0, 0.0, 0.0)
        mag = float(np.mean(self.U_branch))
        cx, cy = self.wind_direction
        return (mag * cx, mag * cy, 0.0)


@dataclass(frozen=True)
class BulkThinningParams:
    """Static wind-model parameters.

    Parameters
    ----------
    U_infty
        Free-stream wind magnitude per z-layer, ``(Nz,)``. Defaults to
        ``np.ones(100)`` so the out-of-the-box behaviour is a uniform
        wind with no boundary layer.
    z_centers
        Layer mid-heights, ``(Nz,)``, strictly monotone-increasing.
        Default ``np.arange(0.25, 50.0, 0.5)`` (100 layers of height
        ``H = 0.5`` from ``z ∈ [0, 50]``).
    H
        Layer height. ``z_centers[0] - H/2`` must be ``<= 0`` so the
        trunk base lies inside the lowest cell.
    C_D
        Cylinder drag coefficient.
    eps
        Divide-by-zero guard inside the actuator-disk formula.
    """

    U_infty: np.ndarray
    z_centers: np.ndarray
    H: float = 0.5
    C_D: float = 1.0
    eps: float = 1e-30

    def __post_init__(self) -> None:
        u = np.asarray(self.U_infty, dtype=float)
        z = np.asarray(self.z_centers, dtype=float)
        if u.ndim != 1 or z.ndim != 1:
            raise ValueError("BulkThinningParams: U_infty and z_centers must be 1-D")
        if u.shape != z.shape:
            raise ValueError(
                f"BulkThinningParams: U_infty and z_centers must share shape "
                f"({u.shape} vs {z.shape})"
            )
        if u.size < 2:
            raise ValueError("BulkThinningParams: need at least 2 z-layers")
        if not np.all(np.diff(z) > 0):
            raise ValueError("BulkThinningParams.z_centers must be strictly monotone increasing")
        if self.H <= 0.0:
            raise ValueError(f"BulkThinningParams.H must be positive, got {self.H}")
        if self.C_D <= 0.0:
            raise ValueError(f"BulkThinningParams.C_D must be positive, got {self.C_D}")
        if float(z[0]) - 0.5 * self.H > 0.0:
            raise ValueError(
                "BulkThinningParams: z_centers[0] - H/2 must be <= 0 so the trunk base "
                f"is in-range; got z_centers[0]={float(z[0])} H={self.H}"
            )
        object.__setattr__(self, "U_infty", u)
        object.__setattr__(self, "z_centers", z)

    @classmethod
    def uniform(cls, U: float = 1.0, *, z_max: float = 50.0, H: float = 0.5) -> BulkThinningParams:
        """Default constructor: uniform free stream, no boundary layer.

        ``U_infty(z) = U`` for all ``z`` from 0 to ``z_max``. Matches
        the simplest sensible default ("wind is there, same everywhere")
        and reproduces the Fortran-default behaviour when ``U = 1`` so
        existing reproductions don't shift.
        """
        z_edges = np.arange(0.0, z_max + 1e-9, H)
        z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
        return cls(U_infty=np.full_like(z_centers, U, dtype=float), z_centers=z_centers, H=H)


def _edges_from_centers(centers: np.ndarray) -> np.ndarray:
    """Reconstruct z-cell edges from monotone centres.

    Port of ``dendroflow.wind.branch._edges_from_centers``.
    """
    if centers.size < 2:
        raise ValueError("Cannot reconstruct edges from fewer than 2 centers.")
    mid = 0.5 * (centers[:-1] + centers[1:])
    dx0 = centers[1] - centers[0]
    dx_end = centers[-1] - centers[-2]
    return np.concatenate([[centers[0] - 0.5 * dx0], mid, [centers[-1] + 0.5 * dx_end]])


def _stack_trees(
    trees: Sequence,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Pull (start, axis, D, L) arrays from every tree and concatenate.

    Mirrors ``mechatree.wind.dendroflow.forest_to_cylinders`` but skips
    the DendroFlow ``Cylinders`` round-trip — we don't need its frame.
    Returns ``(start, axis, D, L, offsets)`` where ``offsets[i]`` is
    the first branch index belonging to tree ``i`` and
    ``offsets[-1] == total_branches``.
    """
    starts: list[np.ndarray] = []
    axes: list[np.ndarray] = []
    Ds: list[np.ndarray] = []
    Ls: list[np.ndarray] = []
    sizes: list[int] = []
    for tree in trees:
        s, a, d, ell = tree.get_branch_data_batch()
        starts.append(s)
        axes.append(a)
        Ds.append(d)
        Ls.append(ell)
        sizes.append(int(s.shape[0]))
    if not starts:
        return (
            np.empty((0, 3), dtype=float),
            np.empty((0, 3), dtype=float),
            np.empty(0, dtype=float),
            np.empty(0, dtype=float),
            np.zeros(1, dtype=np.int64),
        )
    start = np.concatenate(starts)
    axis = np.concatenate(axes)
    D = np.concatenate(Ds)
    L = np.concatenate(Ls)
    offsets = np.zeros(len(trees) + 1, dtype=np.int64)
    np.cumsum(sizes, out=offsets[1:])
    return start, axis, D, L, offsets


def compute_bulk_thinning(
    start: np.ndarray,
    axis: np.ndarray,
    D: np.ndarray,
    L: np.ndarray,
    *,
    branch_offsets: np.ndarray,
    params: BulkThinningParams,
    wind_direction: tuple[float, float],
) -> BulkThinningResult:
    """Run one bulk-thinning solve on stacked branch arrays.

    ``wind_direction`` is the unit horizontal vector ``(cos θ, sin θ)``
    pointing downwind in the world frame. The canopy is conceptually
    rotated so that this direction becomes ``+x``; in practice the
    rotation collapses to a single dot product, so we don't allocate
    rotated coordinate arrays.

    Returns a :class:`BulkThinningResult` with per-branch wind magnitudes
    and the per-layer canopy profile (in the rotated frame).
    """
    n = start.shape[0]
    if n == 0:
        return BulkThinningResult(
            U_branch=np.empty(0, dtype=float),
            z_branch=np.empty(0, dtype=float),
            branch_offsets=np.asarray(branch_offsets, dtype=np.int64),
            U_canopy=params.U_infty.copy(),
            z_centers=params.z_centers.copy(),
            wind_direction=(float(wind_direction[0]), float(wind_direction[1])),
        )

    cx, cy = float(wind_direction[0]), float(wind_direction[1])
    # 1. Per-branch centre height (z is invariant under horizontal rotation).
    z_c = start[:, 2] + 0.5 * L * axis[:, 2]
    # 2. cos I = axis_horizontal · wind_direction.  sin I = sqrt(1 - cos^2 I).
    cos_I = axis[:, 0] * cx + axis[:, 1] * cy
    cos_I = np.clip(cos_I, -1.0, 1.0)
    sin_I = np.sqrt(1.0 - cos_I * cos_I)
    area = D * L * sin_I

    # 3. Bin by z-layer.
    z_edges = _edges_from_centers(params.z_centers)
    bin_area, _ = np.histogram(z_c, bins=z_edges, weights=area)

    # 4. Actuator-disk thinning per layer.
    U_inf = params.U_infty
    U_in_sq = U_inf * U_inf
    F_drag_z = 0.5 * params.C_D * bin_area * U_in_sq
    disc = 1.0 - 4.0 * F_drag_z / (params.H * params.H * U_in_sq + params.eps)
    disc = np.clip(disc, 0.0, 1.0)
    U_canopy = 0.5 * U_inf * (1.0 + np.sqrt(disc))

    # 5. Per-branch lookup. The "+1" branch is to handle z exactly on an
    # edge — searchsorted with side='right' returns Nz at the top edge,
    # which we then clip into [0, Nz-1].
    k_idx = np.clip(np.searchsorted(z_edges, z_c, side="right") - 1, 0, params.z_centers.size - 1)
    U_branch = U_canopy[k_idx]

    return BulkThinningResult(
        U_branch=U_branch,
        z_branch=z_c,
        branch_offsets=np.asarray(branch_offsets, dtype=np.int64),
        U_canopy=U_canopy,
        z_centers=params.z_centers.copy(),
        wind_direction=(cx, cy),
    )


class BulkThinningWindBridge:
    """Stateful ``WindFn`` backed by the native bulk-thinning model.

    Use the 3-arg ``WindFn`` shape ``(generation, rng, context)``. The
    bridge samples a fresh storm direction per call from the configured
    angle distribution (or returns ``+x`` if no distribution is set),
    rotates the canopy into the wind frame, runs the bulk-thinning
    solve, and returns the canopy-mean wind vector in the world frame.

    ``last_result`` carries the most recent
    :class:`BulkThinningResult` for diagnostics and the 3D wind viz.
    """

    def __init__(
        self,
        params: BulkThinningParams | None = None,
        *,
        angle_sampler=None,
        amplitude_sampler=None,
    ) -> None:
        self.params = params if params is not None else BulkThinningParams.uniform()
        self.angle_sampler = angle_sampler
        self.amplitude_sampler = amplitude_sampler
        self.last_result: BulkThinningResult | None = None

    def __call__(
        self,
        generation: int,
        rng: np.random.Generator,
        context: PyTree | Forest,
    ) -> tuple[float, float, float]:
        from mechatree.forest import Forest as _Forest

        trees = context.trees if isinstance(context, _Forest) else [context]
        if not trees:
            return (float(self.params.U_infty[0]), 0.0, 0.0)

        # Sample storm direction (default +x).
        theta = float(self.angle_sampler(rng, 1)[0]) if self.angle_sampler is not None else 0.0
        cx, cy = math.cos(theta), math.sin(theta)

        # Optional amplitude scaling (a multiplier on U_infty for this gen).
        if self.amplitude_sampler is not None:
            amp = float(self.amplitude_sampler(rng, 1)[0])
            scaled_params = BulkThinningParams(
                U_infty=self.params.U_infty * amp,
                z_centers=self.params.z_centers,
                H=self.params.H,
                C_D=self.params.C_D,
                eps=self.params.eps,
            )
        else:
            scaled_params = self.params

        start, axis, D, L, offsets = _stack_trees(trees)
        if start.shape[0] == 0:
            return (float(scaled_params.U_infty[0]) * cx, float(scaled_params.U_infty[0]) * cy, 0.0)

        result = compute_bulk_thinning(
            start,
            axis,
            D,
            L,
            branch_offsets=offsets,
            params=scaled_params,
            wind_direction=(cx, cy),
        )
        self.last_result = result
        return result.canopy_mean


def make_bulk_thinning_wind_fn(
    *,
    U_infty: np.ndarray | Sequence[float] | None = None,
    z_centers: np.ndarray | Sequence[float] | None = None,
    H: float = 0.5,
    C_D: float = 1.0,
    angle_sampler=None,
    amplitude_sampler=None,
) -> BulkThinningWindBridge:
    """Build a :class:`BulkThinningWindBridge`. The default constructor
    yields the simplest sensible setting: uniform free stream
    ``U = 1`` from ``z = 0`` to ``z = 50`` with ``H = 0.5`` (matches
    :meth:`BulkThinningParams.uniform`).
    """
    if U_infty is None and z_centers is None:
        params = BulkThinningParams.uniform(U=1.0)
    else:
        if U_infty is None or z_centers is None:
            raise ValueError(
                "make_bulk_thinning_wind_fn: pass both U_infty and z_centers, or neither"
            )
        params = BulkThinningParams(
            U_infty=np.asarray(U_infty, dtype=float),
            z_centers=np.asarray(z_centers, dtype=float),
            H=float(H),
            C_D=float(C_D),
        )
    return BulkThinningWindBridge(
        params,
        angle_sampler=angle_sampler,
        amplitude_sampler=amplitude_sampler,
    )


__all__ = [
    "BulkThinningParams",
    "BulkThinningResult",
    "BulkThinningWindBridge",
    "compute_bulk_thinning",
    "make_bulk_thinning_wind_fn",
]
