"""Pruning — wind-driven stochastic branch removal."""

from mechatree._core._core import PyTree


def prune(tree: PyTree, wind, leaf_drag_S0: float, cauchy: float) -> int:
    """Run one pruning pass under wind direction ``wind``.

    Returns the number of branches removed (including every descendant of
    any directly-cut branch). The trunk is never cut. Reproducible from a
    single seed via ``tree.set_seed(...)``.
    """
    return tree.prune(wind, leaf_drag_S0, cauchy)


def prune_with_stored_forces(tree: PyTree, leaf_drag_S0: float, cauchy: float) -> int:
    """Run one pruning pass from the per-branch forces pre-stored on the tree.

    Step 25c (option B): instead of recomputing each branch's woody drag from
    a single canopy-mean wind, this reads ``Branch::segment_force_`` (the CFD
    per-branch force) and ``Branch::segment_wind_`` (its local wind, for the
    leaf-cluster drag term). The momentum-wind bridge must have populated
    those fields for the current branch set first (via
    :meth:`PyTree.set_segment_forces_batch` / :meth:`set_segment_winds_batch`).

    Returns the number of branches removed (including every descendant of any
    directly-cut branch). The trunk is never cut.
    """
    return tree.prune_with_stored_forces(leaf_drag_S0, cauchy)


__all__ = ["prune", "prune_with_stored_forces"]
