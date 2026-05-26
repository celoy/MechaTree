"""Light interception and aggregation onto trees.

``intercept`` ports ``legacy/fortran/mod_tree.f90:219`` (light_interception)
across every direction in a ``Sun``. Per direction:

1. Rotate every leaf so the sun is along the new Z axis (``Xp/Yp/Zp``).
2. Bin each leaf's ``(Xp, Yp)`` into an integer cell of size ``size_leaf``.
3. Each leaf's light is ``tau ** depth``, where ``depth`` is its 0-indexed
   rank in Z'-descending order within its cell (the topmost gets
   ``tau^0 = 1``). ``tau = 0`` recovers the Fortran's binary
   topmost-wins regime; ``tau = 1`` makes leaves fully transparent; the
   default ``tau = 0.5`` matches Eloy et al. (Nat Commun 2017).

The Fortran reference walks leaves in Z'-descending order through a
shadow counter. We instead use ``np.lexsort + np.unique(return_index=...)``
to identify the in-cell depth of every leaf in one vectorised pass —
~5× faster on realistic leaf counts and equivalent up to tie-breaks
(ties broken by leaf index ascending, same as the Fortran's stable
sort by original index).

``aggregate_onto_trees`` ports ``light_on_trees`` (:273): each leaf's
scalar light is the mean of its per-direction lights, written back onto
the source ``Branch`` via ``PyTree.set_light``.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from mechatree._core import PyTree
from mechatree._core._core import light_intercept_kernel
from mechatree.light.leaves import Leaves
from mechatree.light.sun import Sun


def intercept(leaves: Leaves, sun: Sun, leaf_transparency: float = 0.5) -> None:
    """Populate ``leaves.light_per_direction`` in place.

    ``leaf_transparency`` (``tau``, in [0, 1]) sets how much light passes
    through a leaf: the i-th leaf from the top of a shadow cell receives
    ``tau**i``. ``tau = 0`` is the Fortran binary regime, ``tau = 1`` is
    fully transparent, default ``tau = 0.5`` matches the value used in
    Eloy et al. (Nat Commun 2017).

    Allocates a fresh ``light_per_direction`` array if its shape doesn't
    match ``(n_leaves, sun.n_directions)`` — the caller can avoid the
    realloc by passing a Leaves built with
    ``extract_leaves(..., n_directions=sun.n_directions)``.

    Body is in C++ (:func:`mechatree._core._core.light_intercept_kernel`)
    so the per-direction rotate + bin + sort loop runs without per-leaf
    NumPy allocations. At island scale (R = 200 L, ~430k leaves, 32 sun
    directions) the C++ kernel is ~5× faster than the equivalent
    NumPy / lexsort path the function used until Step 21b.
    """
    if not (0.0 <= leaf_transparency <= 1.0):
        raise ValueError(f"intercept: leaf_transparency must be in [0, 1], got {leaf_transparency}")

    n = leaves.n_leaves
    n_dir = sun.n_directions

    if leaves.light_per_direction.shape != (n, n_dir):
        leaves.light_per_direction = np.zeros((n, n_dir), dtype=np.float64)
    else:
        leaves.light_per_direction.fill(0.0)

    if n == 0 or n_dir == 0:
        return

    # Make sure the C++ kernel gets contiguous float64 views. ``Sun`` stores
    # elev/azim as float64 ndarrays today; cast defensively.
    loc = np.ascontiguousarray(leaves.location, dtype=np.float64)
    elev = np.ascontiguousarray(sun.elev, dtype=np.float64)
    azim = np.ascontiguousarray(sun.azim, dtype=np.float64)
    if not leaves.light_per_direction.flags["C_CONTIGUOUS"]:
        leaves.light_per_direction = np.ascontiguousarray(leaves.light_per_direction)

    light_intercept_kernel(
        loc, elev, azim, sun.size_leaf, leaf_transparency, leaves.light_per_direction
    )


def aggregate_onto_trees(leaves: Leaves, trees: Iterable[PyTree]) -> None:
    """Write the mean per-direction light back onto each tree's branches.

    Mirrors ``light_on_trees`` (mod_tree.f90:273). For each leaf,
    ``branch.light = mean(light_per_direction[i, :])``. Branches that own
    no leaf (internal nodes) are not touched.

    Uses the C++ batched setter :meth:`PyTree.set_lights_batch` — one
    Cython call per tree instead of one per leaf. At island scale this
    cuts ``aggregate_onto_trees`` wall-clock by ~10×.
    """
    trees_list = list(trees)
    n_dir = leaves.n_directions
    if leaves.n_leaves == 0 or n_dir == 0:
        return

    # mean across directions — same as Fortran's `sum / (Nazim*Nelev)`.
    mean_light = leaves.light_per_direction.mean(axis=1)
    branch_idx = leaves.branch_index
    tree_idx = leaves.tree_index

    # Group leaves by tree_index and dispatch one batched setter per tree.
    # The Leaves builder concatenates trees in order, so contiguous slices
    # work — `np.unique(return_index=True)` on the sorted tree_idx gives
    # us each tree's slice in O(N).
    if trees_list:
        unique_trees, slice_starts = np.unique(tree_idx, return_index=True)
        slice_ends = np.append(slice_starts[1:], leaves.n_leaves)
        for ti, start, end in zip(unique_trees, slice_starts, slice_ends, strict=True):
            trees_list[int(ti)].set_lights_batch(branch_idx[start:end], mean_light[start:end])
