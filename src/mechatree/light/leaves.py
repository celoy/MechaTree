"""``Leaves`` — a struct-of-arrays of every leaf in a tree or a forest.

The hot path of light interception runs over thousands of leaves and tens
of sun directions every generation, so a Python ``list[Leaf]`` is the
wrong shape. The arrays here are the natural input for the NumPy /
rasterisation work in ``mechatree.light.interception``.

For Step 10 the per-leaf attributes are intentionally minimal — only what
the Fortran reference reads. ``diameter`` and ``transparency`` named in
the design doc are deferred until there is a concrete model that uses
them.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np

from mechatree._core import PyTree


@dataclass
class Leaves:
    """Struct-of-arrays for every leaf in a tree or a forest.

    Attributes
    ----------
    location:
        ``(n, 3)`` ``float64``. World-space coordinates of each leaf
        (branch tip = ``branch.location + branch.length * branch.unit_t``).
    branch_index:
        ``(n,)`` ``int32``. Index of the leaf-bearing branch within its
        owning tree (the value to pass back to ``PyTree.set_light``).
    tree_index:
        ``(n,)`` ``int32``. Which tree in the input list owns this leaf.
        Always zero for a single-tree call.
    light_per_direction:
        ``(n, n_dir)`` ``float64``. Populated in place by ``intercept``;
        used by ``aggregate_onto_trees`` to write the averaged scalar
        ``branch.light`` back onto each tree.
    """

    location: np.ndarray
    branch_index: np.ndarray
    tree_index: np.ndarray
    light_per_direction: np.ndarray

    @property
    def n_leaves(self) -> int:
        return int(self.location.shape[0])

    @property
    def n_directions(self) -> int:
        return int(self.light_per_direction.shape[1]) if self.n_leaves else 0


def extract_leaves(trees: Iterable[PyTree], n_directions: int = 0) -> Leaves:
    """Build a ``Leaves`` snapshot of every leaf in ``trees``.

    ``trees`` is always an iterable of ``PyTree`` — for a single tree,
    pass ``[tree]``. The output's ``light_per_direction`` is allocated as
    a zero ``(n, n_directions)`` array; ``intercept`` fills it in. Passing
    ``n_directions=0`` (the default) is fine for callers that only need
    the geometry.
    """
    trees_list = list(trees)

    locs: list[np.ndarray] = []
    branch_idx: list[int] = []
    tree_idx: list[int] = []

    for ti, tree in enumerate(trees_list):
        # We need the trunk's location, length, unit_t to compute the leaf
        # tip per the Fortran convention. PyTree gives us those one branch
        # at a time — there's no batched accessor in PR1's API. That's
        # fine; this is called once per generation.
        for bi in tree.leaf_indices():
            base = np.asarray(tree.get_location(bi), dtype=np.float64)
            unit_t = np.asarray(tree.get_unit_t(bi), dtype=np.float64)
            length = float(tree.get_length(bi))
            tip = base + length * unit_t
            locs.append(tip)
            branch_idx.append(bi)
            tree_idx.append(ti)

    n = len(locs)
    location = np.asarray(locs, dtype=np.float64) if n else np.zeros((0, 3), dtype=np.float64)
    return Leaves(
        location=location,
        branch_index=np.asarray(branch_idx, dtype=np.int32),
        tree_index=np.asarray(tree_idx, dtype=np.int32),
        light_per_direction=np.zeros((n, n_directions), dtype=np.float64),
    )
