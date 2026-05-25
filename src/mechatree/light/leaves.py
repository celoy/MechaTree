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

    Uses the C++ batched accessor :meth:`PyTree.get_leaf_tips_batch` so
    every leaf in a tree is read in one Cython call — required for
    island-scale (Step 21b) where the per-leaf Python loop dominated
    light-interception cost.
    """
    trees_list = list(trees)
    if not trees_list:
        return Leaves(
            location=np.zeros((0, 3), dtype=np.float64),
            branch_index=np.zeros(0, dtype=np.int32),
            tree_index=np.zeros(0, dtype=np.int32),
            light_per_direction=np.zeros((0, n_directions), dtype=np.float64),
        )

    per_tree_tips: list[np.ndarray] = []
    per_tree_branches: list[np.ndarray] = []
    per_tree_indices: list[np.ndarray] = []
    for ti, tree in enumerate(trees_list):
        tips, branches = tree.get_leaf_tips_batch()
        if tips.shape[0] == 0:
            continue
        per_tree_tips.append(tips)
        per_tree_branches.append(branches)
        per_tree_indices.append(np.full(tips.shape[0], ti, dtype=np.int32))

    if not per_tree_tips:
        return Leaves(
            location=np.zeros((0, 3), dtype=np.float64),
            branch_index=np.zeros(0, dtype=np.int32),
            tree_index=np.zeros(0, dtype=np.int32),
            light_per_direction=np.zeros((0, n_directions), dtype=np.float64),
        )

    location = np.concatenate(per_tree_tips, axis=0)
    branch_index = np.concatenate(per_tree_branches)
    tree_index = np.concatenate(per_tree_indices)
    return Leaves(
        location=location,
        branch_index=branch_index,
        tree_index=tree_index,
        light_per_direction=np.zeros((location.shape[0], n_directions), dtype=np.float64),
    )
