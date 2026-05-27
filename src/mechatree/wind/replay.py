"""Storm replay — capture each pass of the Step-24 fixed-point loop.

A targeted diagnostic for understanding what the coupled wind ↔ pruning
loop is doing under a single storm. Replays the loop on a live
:class:`~mechatree.forest.Forest`, snapshotting per-iteration cuts so a
notebook can show "the tree before the storm" → "after iter 1, these
branches fell" → "after iter 2, these too" etc., typically rendering
cut branches in black on top of the pre-storm canopy.

Typical use::

    forest = Forest(cfg, seed=0)
    for g in range(100):
        forest.step(g)  # grow under whatever wind
    storm_wind = make_momentum_wind_fn(grid_size=2.0, U_uniform=1.6)
    pre, snapshots = run_storm_replay(
        forest, storm_wind, generation=100, leaf_drag_S0=..., cauchy=...
    )
    plot_storm_replay(pre, snapshots)

The replay does NOT run growth/light/seed phases — it isolates the
storm. The forest is mutated in place; if you want to keep the
pre-storm forest around for comparison, deep-copy or re-seed.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from mechatree.pruning import prune

if TYPE_CHECKING:
    from mechatree.forest import Forest


@dataclass(frozen=True)
class StormPreSnapshot:
    """Per-tree pre-storm canopy geometry.

    Captured once before the storm starts. The replay then references
    every branch by its index in these arrays, since pruning shifts
    branch indices inside the live tree but the pre-storm arrays are
    immutable.

    Attributes
    ----------
    per_tree_geometry
        Length-``n_trees`` list of dicts with keys ``start``, ``axis``,
        ``D``, ``L`` — each a ``(n_branches_in_tree, …)`` numpy array.
    """

    per_tree_geometry: list[dict[str, np.ndarray]]

    @property
    def n_trees(self) -> int:
        return len(self.per_tree_geometry)

    @property
    def n_branches_total(self) -> int:
        return sum(g["start"].shape[0] for g in self.per_tree_geometry)


@dataclass(frozen=True)
class StormSnapshot:
    """One frame of the storm replay.

    ``iteration == 0`` is the pre-storm state (all alive); ``iteration
    == 1, 2, …`` are after each iteration of the fixed-point loop.

    Attributes
    ----------
    iteration
        0 = pre-storm; 1, 2, … = after the matching prune sweep.
    n_pruned_this_iter
        Total branches cut across the whole forest this iteration.
    n_pruned_cumulative
        Running total since the start of the storm.
    n_trees_affected_this_iter
        How many trees lost at least one branch in this iter — the
        user-visible "how many trees got hit on this pass".
    wind_used
        The ``(x, y, z)`` wind vector this iteration's prune saw. For
        iteration 0 the value is ``(0, 0, 0)`` (no wind applied yet).
    alive_mask_per_tree
        Length-``n_trees`` list of ``(n_branches_in_pre_storm_tree,)``
        boolean arrays. ``True`` means the matching pre-storm branch is
        still alive after this iteration.
    per_tree_n_alive
        Shape-``(n_trees,)`` int. Cached ``alive_mask.sum()`` per tree.
    per_tree_n_pruned_this_iter
        Shape-``(n_trees,)`` int. Cuts contributed by each tree on this
        iteration.
    """

    iteration: int
    n_pruned_this_iter: int
    n_pruned_cumulative: int
    n_trees_affected_this_iter: int
    wind_used: tuple[float, float, float]
    alive_mask_per_tree: list[np.ndarray]
    per_tree_n_alive: np.ndarray
    per_tree_n_pruned_this_iter: np.ndarray


def _start_lookup(starts: np.ndarray) -> set[tuple[float, float, float]]:
    """Build a set of branch base positions for O(1) membership lookup.

    Pruning never moves a surviving branch, so the start position is a
    stable key from pre-storm to post-storm.
    """
    return {tuple(row) for row in starts}


def _alive_mask(pre_starts: np.ndarray, current_starts: np.ndarray) -> np.ndarray:
    """For each pre-storm branch, is its start position still present?"""
    current_set = _start_lookup(current_starts)
    return np.array(
        [tuple(row) in current_set for row in pre_starts],
        dtype=bool,
    )


def run_storm_replay(
    forest: Forest,
    wind_fn: Callable,
    *,
    generation: int = 0,
    max_iterations: int = 8,
    eps_rel: float = 0.01,
    leaf_drag_S0: float,
    cauchy: float,
    rng: np.random.Generator | None = None,
) -> tuple[StormPreSnapshot, list[StormSnapshot]]:
    """Run one storm through the fixed-point loop, capturing each pass.

    Mirrors the loop in :meth:`mechatree.forest.Forest.step` step-3, but
    extracted so we can snapshot iter-by-iter. Mutates ``forest.trees``
    in place — the trees end the call in the post-storm state.

    Parameters
    ----------
    forest
        The forest to storm. Mutated in place.
    wind_fn
        A 3-arg ``WindFn(generation, rng, context)`` (canopy-aware).
        A 2-arg fn works too (context ignored) but the iter-2 ε check
        will always trigger because the wind doesn't depend on the
        canopy — only iter 1 will run, which is correct (single-pass
        semantics).
    generation
        Generation index passed to ``wind_fn``. Default 0.
    max_iterations, eps_rel
        Match the Step-24 fixed-point cap and ε early-exit.
    leaf_drag_S0, cauchy
        Forwarded to :func:`mechatree.pruning.prune` exactly as the
        normal ``Forest.step`` would.
    rng
        Optional ``numpy.random.Generator``; defaults to a fresh
        ``default_rng(0)`` so replay is reproducible.

    Returns
    -------
    ``(pre, snapshots)``
        ``pre`` carries the per-tree geometry as captured *before* the
        storm. ``snapshots`` is a list starting with the iteration-0
        (all-alive) frame, then one entry per loop iteration.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    # 1. Snapshot pre-storm geometry per tree.
    per_tree_geometry: list[dict[str, np.ndarray]] = []
    pre_start_lookup: list[set[tuple[float, float, float]]] = []
    for tree in forest.trees:
        s, a, d, L = tree.get_branch_data_batch()
        per_tree_geometry.append({"start": s, "axis": a, "D": d, "L": L})
        pre_start_lookup.append(_start_lookup(s))
    pre = StormPreSnapshot(per_tree_geometry=per_tree_geometry)

    # 2. iteration-0 snapshot (all alive).
    snapshots: list[StormSnapshot] = []
    alive_masks0 = [np.ones(g["start"].shape[0], dtype=bool) for g in per_tree_geometry]
    snapshots.append(
        StormSnapshot(
            iteration=0,
            n_pruned_this_iter=0,
            n_pruned_cumulative=0,
            n_trees_affected_this_iter=0,
            wind_used=(0.0, 0.0, 0.0),
            alive_mask_per_tree=alive_masks0,
            per_tree_n_alive=np.array([m.size for m in alive_masks0], dtype=int),
            per_tree_n_pruned_this_iter=np.zeros(len(alive_masks0), dtype=int),
        )
    )

    # 3. Fixed-point loop. Mirrors mechatree.forest.Forest.step lines
    #    ~213-260 but with per-iter snapshotting.
    n_pruned_cumulative = 0
    prev_wind: tuple[float, float, float] | None = None
    wind_arity = _arity(wind_fn)
    for n_iter in range(1, max_iterations + 1):
        wind = wind_fn(generation, rng, forest) if wind_arity >= 3 else wind_fn(generation, rng)

        # ε early-exit (only from iter 2 onwards).
        if prev_wind is not None and eps_rel > 0.0:
            delta = math.hypot(wind[0] - prev_wind[0], wind[1] - prev_wind[1])
            ref = max(math.hypot(prev_wind[0], prev_wind[1]), 1e-6)
            if delta / ref < eps_rel:
                break

        per_tree_cuts = np.zeros(len(forest.trees), dtype=int)
        for ti, tree in enumerate(forest.trees):
            n_cut = prune(tree, wind=wind, leaf_drag_S0=leaf_drag_S0, cauchy=cauchy)
            tree.reorder()
            per_tree_cuts[ti] = n_cut
        n_cut_this_iter = int(per_tree_cuts.sum())
        n_pruned_cumulative += n_cut_this_iter

        # Recompute alive masks vs pre-storm geometry.
        alive_per_tree: list[np.ndarray] = []
        per_tree_n_alive: list[int] = []
        for tree, pre_starts in zip(
            forest.trees, [g["start"] for g in per_tree_geometry], strict=True
        ):
            s_now, _, _, _ = tree.get_branch_data_batch()
            mask = _alive_mask(pre_starts, s_now)
            alive_per_tree.append(mask)
            per_tree_n_alive.append(int(mask.sum()))

        snapshots.append(
            StormSnapshot(
                iteration=n_iter,
                n_pruned_this_iter=n_cut_this_iter,
                n_pruned_cumulative=n_pruned_cumulative,
                n_trees_affected_this_iter=int(np.sum(per_tree_cuts > 0)),
                wind_used=(float(wind[0]), float(wind[1]), float(wind[2])),
                alive_mask_per_tree=alive_per_tree,
                per_tree_n_alive=np.array(per_tree_n_alive, dtype=int),
                per_tree_n_pruned_this_iter=per_tree_cuts,
            )
        )

        prev_wind = wind
        if n_cut_this_iter == 0:
            break

    return pre, snapshots


def _arity(fn: Callable) -> int:
    """Best-effort positional arg count; matches simulate._callback_arity
    semantics so a 2-arg ``wind_fn`` is treated as such."""
    import inspect

    try:
        sig = inspect.signature(fn)
    except (TypeError, ValueError):
        return 2
    return sum(
        1
        for p in sig.parameters.values()
        if p.kind in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD)
    )


__all__ = [
    "StormPreSnapshot",
    "StormSnapshot",
    "run_storm_replay",
]
