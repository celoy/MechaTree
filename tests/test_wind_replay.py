"""Tests for the Step-25 storm-replay helper + plot."""

from __future__ import annotations

import mechatree as mt
from mechatree.config import Config, ForestConfig, TreeConfig
from mechatree.wind.momentum_wind import make_momentum_wind_fn
from mechatree.wind.replay import StormPreSnapshot, run_storm_replay


def _build_warmed_forest(seed: int = 0, n_init: int = 5, gens: int = 25) -> mt.Forest:
    """Grow a small forest with zero wind so the canopy gets non-trivial
    without being pre-thinned. Used as the pre-storm state for the
    replay tests."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=15.0, n_trees_init=n_init, n_trees_max=50),
    )
    forest = mt.Forest(cfg, seed=seed, wind_fn=lambda g, r: (0.0, 0.0, 0.0))
    for g in range(gens):
        forest.step(g)
    return forest


def _momentum_storm(u: float):
    """A strong momentum-wind storm at uniform inlet ``u`` (the screening is
    resolved through the canopy by the kernel)."""
    return make_momentum_wind_fn(grid_size=2.0, pad_x=10.0, U_uniform=u)


def test_pre_snapshot_captures_full_canopy():
    f = _build_warmed_forest()
    n_branches_pre = sum(t.get_number_of_branches() for t in f.trees)
    pre, _ = run_storm_replay(
        f,
        lambda g, r, ctx: (0.0, 0.0, 0.0),
        leaf_drag_S0=f.config.tree.leaf_surface,
        cauchy=f.config.tree.cauchy,
        max_iterations=1,
    )
    assert isinstance(pre, StormPreSnapshot)
    assert pre.n_trees == len(f.trees)
    assert pre.n_branches_total == n_branches_pre


def test_calm_wind_replays_zero_iterations_after_initial():
    """Zero wind → iter 1 cuts nothing → loop exits → snapshots = [iter0, iter1]."""
    f = _build_warmed_forest()
    pre, snaps = run_storm_replay(
        f,
        lambda g, r, ctx: (0.0, 0.0, 0.0),
        leaf_drag_S0=f.config.tree.leaf_surface,
        cauchy=f.config.tree.cauchy,
        max_iterations=6,
    )
    assert len(snaps) == 2  # iter 0 (all alive) + iter 1 (zero cuts)
    assert snaps[0].iteration == 0
    assert snaps[0].n_pruned_this_iter == 0
    assert all(m.all() for m in snaps[0].alive_mask_per_tree)
    assert snaps[1].iteration == 1
    assert snaps[1].n_pruned_this_iter == 0


def test_storm_replay_with_momentum_wind():
    """Strong storm via the momentum-wind bridge: at least one iteration
    should cut branches, and per-tree counts should sum to the forest-wide
    total."""
    f = _build_warmed_forest(n_init=6, gens=30)
    pre, snaps = run_storm_replay(
        f,
        _momentum_storm(5.0),
        leaf_drag_S0=f.config.tree.leaf_surface,
        cauchy=f.config.tree.cauchy,
        max_iterations=4,
        eps_rel=0.01,
    )
    storm_snaps = [s for s in snaps if s.iteration >= 1]
    assert any(s.n_pruned_this_iter > 0 for s in storm_snaps)
    for s in storm_snaps:
        assert s.per_tree_n_pruned_this_iter.sum() == s.n_pruned_this_iter
        # alive_mask per tree gives the right per-tree alive count.
        for mask, n_alive in zip(s.alive_mask_per_tree, s.per_tree_n_alive, strict=True):
            assert int(mask.sum()) == int(n_alive)


def test_storm_replay_cumulative_pruned_monotone():
    """Cumulative cut count must be non-decreasing across snapshots."""
    f = _build_warmed_forest(n_init=4, gens=25)
    pre, snaps = run_storm_replay(
        f,
        _momentum_storm(4.0),
        leaf_drag_S0=f.config.tree.leaf_surface,
        cauchy=f.config.tree.cauchy,
        max_iterations=4,
        eps_rel=0.0,
    )
    cum = [s.n_pruned_cumulative for s in snaps]
    assert all(cum[i + 1] >= cum[i] for i in range(len(cum) - 1))


def test_plot_storm_replay_builds_figure():
    """End-to-end: the viz helper builds a multi-panel plotly figure."""
    f = _build_warmed_forest(n_init=3, gens=20)
    pre, snaps = run_storm_replay(
        f,
        _momentum_storm(4.0),
        leaf_drag_S0=f.config.tree.leaf_surface,
        cauchy=f.config.tree.cauchy,
        max_iterations=3,
    )
    fig = mt.plot_storm_replay(pre, snaps, n_cols=2)
    # Each snapshot contributes 1 or 2 traces (alive + optional cut),
    # so number of traces is between len(snaps) and 2*len(snaps).
    assert len(fig.data) >= len(snaps)
    assert len(fig.data) <= 2 * len(snaps)


def test_storm_replay_with_two_arg_wind_works():
    """A 2-arg wind callable (no canopy context) should also work — the
    replay simply doesn't iterate past iter 1 because the wind never
    changes between calls."""
    f = _build_warmed_forest(n_init=3, gens=20)
    pre, snaps = run_storm_replay(
        f,
        lambda g, r: (5.0, 0.0, 0.0),
        leaf_drag_S0=f.config.tree.leaf_surface,
        cauchy=f.config.tree.cauchy,
        max_iterations=4,
        eps_rel=0.01,
    )
    # Iter 1 runs; iter 2 sees identical wind → ε early-exit.
    assert len(snaps) >= 2
