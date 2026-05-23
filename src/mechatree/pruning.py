"""Pruning — wind-driven stochastic branch removal."""

from mechatree._core._core import PyTree


def prune(tree: PyTree, wind, leaf_drag_S0: float, cauchy: float) -> int:
    """Run one pruning pass under wind direction ``wind``.

    Returns the number of branches removed (including every descendant of
    any directly-cut branch). The trunk is never cut. Reproducible from a
    single seed via ``tree.set_seed(...)``.
    """
    return tree.prune(wind, leaf_drag_S0, cauchy)


__all__ = ["prune"]
