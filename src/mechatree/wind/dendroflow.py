"""DendroFlow bridge — wrap ``BulkThinningBranchWindModel`` as a ``WindFn``.

This is the MechaTree side of DendroFlow's milestone M6 (cf. that project's
``CLAUDE.md``). Per generation, the bridge:

1. Snapshots the live geometry of a ``PyTree`` (or all trees in a ``Forest``)
   into a DendroFlow ``Cylinders`` via ``dendroflow.from_arrays``.
2. Calls ``BulkThinningBranchWindModel.compute(...)`` — a single-pass 1-D
   vertical thinning of the inflow profile by the canopy's projected
   cross-section per z-layer. Cost: < 1 ms at 10 k cylinders.
3. Returns the result's ``canopy_mean`` — a streamwise-only ``(Ū, 0, 0)``
   triple that plugs straight into MechaTree's pruning step.

The bridge is constructed once and reused across generations (the DendroFlow
model instance is cheap but worth keeping around).

DendroFlow is an *optional* dependency. Importing this module without it
raises a friendly ImportError with the install hint.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

try:
    from dendroflow import from_arrays as _df_from_arrays
    from dendroflow.wind import BulkThinningBranchWindModel
except ImportError as _err:  # pragma: no cover - exercised by docs only
    raise ImportError(
        "mechatree.wind.dendroflow requires the DendroFlow package. "
        "Install with: pip install 'mechatree[dendroflow]' "
        "(DendroFlow is not yet on PyPI; for local dev, "
        "`uv pip install -e ../DendroFlow`)."
    ) from _err

if TYPE_CHECKING:
    from collections.abc import Sequence

    from dendroflow import Cylinders
    from dendroflow.wind import BranchWindResult

    from mechatree._core import PyTree
    from mechatree.forest import Forest


_Z_MODES = ("mean", "max", "base")


@dataclass(frozen=True)
class DendroFlowWindParams:
    """Static inputs consumed by ``BulkThinningBranchWindModel.compute``.

    Parameters
    ----------
    U_infty
        ``(Nz,)`` inflow wind magnitude at each z-layer. Same units as the
        result.
    z_centers
        ``(Nz,)`` cell-center heights, strictly monotone increasing. Must
        cover ``z = 0`` (i.e. ``z_centers[0] - H/2 <= 0``) so the trunk
        base isn't dropped by the histogram.
    H
        z-cell height.
    C_D
        Drag coefficient applied to every cylinder cross-section.
    """

    U_infty: np.ndarray
    z_centers: np.ndarray
    H: float = 0.5
    C_D: float = 1.0

    def __post_init__(self) -> None:
        u = np.asarray(self.U_infty, dtype=float)
        z = np.asarray(self.z_centers, dtype=float)
        if u.ndim != 1 or z.ndim != 1:
            raise ValueError("U_infty and z_centers must be 1-D arrays")
        if u.shape != z.shape:
            raise ValueError(f"U_infty and z_centers must share shape; got {u.shape} vs {z.shape}")
        if u.size < 2:
            raise ValueError("U_infty / z_centers must have at least 2 elements")
        if not np.all(np.diff(z) > 0):
            raise ValueError("z_centers must be strictly monotone increasing")
        if self.H <= 0.0:
            raise ValueError(f"H must be positive, got {self.H}")
        if self.C_D <= 0.0:
            raise ValueError(f"C_D must be positive, got {self.C_D}")
        if z[0] - 0.5 * self.H > 0.0:
            raise ValueError(
                "z_centers[0] - H/2 must be <= 0 so the trunk base is in-range; "
                f"got z_centers[0]={z[0]} H={self.H}"
            )
        # Replace with normalised numpy arrays (frozen dataclass requires
        # object.__setattr__ to mutate).
        object.__setattr__(self, "U_infty", u)
        object.__setattr__(self, "z_centers", z)

    def to_namespace(self) -> SimpleNamespace:
        """Pack into the SimpleNamespace shape DendroFlow expects."""
        return SimpleNamespace(
            U_infty=self.U_infty,
            z_centers=self.z_centers,
            H=self.H,
            C_D=self.C_D,
        )


def pytree_to_cylinders(tree: PyTree, *, tree_id: float = 1.0) -> Cylinders:
    """Build a DendroFlow ``Cylinders`` from a single MechaTree ``PyTree``.

    Each branch becomes one cylinder. ``start`` is the branch base location;
    ``axis`` is the unit-t direction (DendroFlow normalises internally).
    DendroFlow requires ``tree_id`` to be a float in its schema.
    """
    n = tree.get_number_of_branches()
    start = np.empty((n, 3), dtype=float)
    axis = np.empty((n, 3), dtype=float)
    D = np.empty(n, dtype=float)
    L = np.empty(n, dtype=float)
    for i in range(n):
        start[i] = tree.get_location(i)
        axis[i] = tree.get_unit_t(i)
        D[i] = tree.get_diameter(i)
        L[i] = tree.get_length(i)
    return _df_from_arrays(start=start, axis=axis, D=D, L=L, tree_id=float(tree_id))


def forest_to_cylinders(trees: Sequence[PyTree]) -> Cylinders:
    """Pool many ``PyTree`` instances into one ``Cylinders``.

    ``tree_id`` is the tree's index in the sequence (as float). When
    ``trees`` is empty, the caller should short-circuit upstream — building
    a zero-row ``Cylinders`` is not useful here and DendroFlow would
    happily produce one but ``compute`` would then divide by zero.
    """
    if len(trees) == 0:
        raise ValueError("forest_to_cylinders requires at least one tree")

    sizes = [t.get_number_of_branches() for t in trees]
    total = sum(sizes)
    start = np.empty((total, 3), dtype=float)
    axis = np.empty((total, 3), dtype=float)
    D = np.empty(total, dtype=float)
    L = np.empty(total, dtype=float)
    tree_id = np.empty(total, dtype=float)

    offset = 0
    for ti, (tree, n) in enumerate(zip(trees, sizes, strict=True)):
        for i in range(n):
            start[offset + i] = tree.get_location(i)
            axis[offset + i] = tree.get_unit_t(i)
            D[offset + i] = tree.get_diameter(i)
            L[offset + i] = tree.get_length(i)
        tree_id[offset : offset + n] = float(ti)
        offset += n

    return _df_from_arrays(start=start, axis=axis, D=D, L=L, tree_id=tree_id)


class BranchWindBridge:
    """Stateful ``WindFn`` that delegates to DendroFlow's lean wind model.

    Use the 3-arg ``WindFn`` shape ``(generation, rng, context)``, where
    ``context`` is the live ``PyTree`` (in ``grow_tree``) or the ``Forest``
    (in ``Forest.step``). ``rng`` is unused by the deterministic
    ``BulkThinningBranchWindModel`` but accepted for signature compatibility.

    The most recent ``BranchWindResult`` is stashed on ``self.last_result``
    after each call — useful for diagnostics and per-generation plotting.
    """

    def __init__(
        self,
        params: DendroFlowWindParams,
        *,
        z_representative: Literal["mean", "max", "base"] = "mean",
    ) -> None:
        if z_representative not in _Z_MODES:
            raise ValueError(
                f"z_representative must be one of {_Z_MODES}; got {z_representative!r}"
            )
        self.params = params
        self.z_representative = z_representative
        self._params_ns = params.to_namespace()
        self._model = BulkThinningBranchWindModel(z_representative=z_representative)
        self.last_result: BranchWindResult | None = None

    def __call__(
        self,
        generation: int,
        rng: np.random.Generator,
        context: PyTree | Forest,
    ) -> tuple[float, float, float]:
        # Lazy import — avoids a circular import (forest -> simulate -> wind).
        from mechatree.forest import Forest as _Forest

        if isinstance(context, _Forest):
            trees = context.trees
            if not trees:
                return (float(self.params.U_infty[0]), 0.0, 0.0)
            cylinders = forest_to_cylinders(trees)
        else:
            cylinders = pytree_to_cylinders(context)

        if cylinders.n == 0:
            return (float(self.params.U_infty[0]), 0.0, 0.0)

        result = self._model.compute(cylinders, wind_params=self._params_ns)
        self.last_result = result
        return result.canopy_mean


def make_dendroflow_wind_fn(
    *,
    U_infty: Any,
    z_centers: Any,
    H: float = 0.5,
    C_D: float = 1.0,
    z_representative: Literal["mean", "max", "base"] = "mean",
) -> BranchWindBridge:
    """Build a :class:`BranchWindBridge` from raw kwargs.

    This is the entry point used by YAML wiring in
    :class:`mechatree.config.WindConfig` and the most convenient one for
    Python callers.
    """
    params = DendroFlowWindParams(
        U_infty=np.asarray(U_infty, dtype=float),
        z_centers=np.asarray(z_centers, dtype=float),
        H=float(H),
        C_D=float(C_D),
    )
    return BranchWindBridge(params, z_representative=z_representative)


__all__ = [
    "BranchWindBridge",
    "DendroFlowWindParams",
    "forest_to_cylinders",
    "make_dendroflow_wind_fn",
    "pytree_to_cylinders",
]
