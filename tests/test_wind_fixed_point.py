"""Tests for the Step 24 wind ↔ pruning fixed-point loop.

The loop fires only when the wind callback is 3-arg (canopy-aware,
e.g. the momentum-wind bridge). The default 2-arg wind path is
byte-identical to before Step 24.

A note on the loop's direction: for the bulk-thinning canopy model
(a bulk-thinning canopy model, and the synthetic
``canopy_intensifying_wind`` below), thinning the canopy raises the
mean wind on the survivors, so the loop iterates toward MORE cuts than
the single-pass would produce, not fewer. That's the physically
correct settled state — the single-pass under-prunes because it scores
the canopy against the wind the pre-pruning forest produces, which is
weaker than what the post-pruning forest produces.
"""

from __future__ import annotations

import time
import warnings

import numpy as np
import pytest

from mechatree.config import Config, ForestConfig, TreeConfig, WindConfig
from mechatree.forest import Forest
from mechatree.simulate import (
    TreeStats,
    _prune_to_fixed_point,
    default_wind_fn,
    grow_tree,
    make_seed_tree,
)

# ---------------------------------------------------------------------------
# Helpers.
# ---------------------------------------------------------------------------


class CanopyIntensifyingWind:
    """A toy 3-arg ``WindFn`` whose magnitude grows as the canopy thins.

    Mirrors the qualitative behaviour of a bulk-thinning canopy
    (verified empirically: dropping branches raises ``canopy_mean``)
    while staying fully deterministic and self-contained. Use
    for tests that need a controlled canopy-aware wind without
    without any optional extras.

    The wind on a tree/forest with ``n`` branches is
    ``(base * (1 + slope * (n0 - n) / n0), 0, 0)``: equal to ``base`` at
    the initial size ``n0`` and rising as ``n`` falls.
    """

    def __init__(self, base: float, n0: int, slope: float = 1.5) -> None:
        self.base = base
        self.n0 = max(1, int(n0))
        self.slope = slope
        self.calls: list[int] = []  # tracks size of canopy at each call

    def __call__(self, gen, rng, context) -> tuple[float, float, float]:
        if isinstance(context, Forest):
            n = sum(t.get_number_of_branches() for t in context.trees)
        else:
            n = context.get_number_of_branches()
        self.calls.append(n)
        scale = 1.0 + self.slope * max(0.0, (self.n0 - n) / self.n0)
        return (self.base * scale, 0.0, 0.0)


# ---------------------------------------------------------------------------
# 1. Default 2-arg wind path: n_wind_iterations stays at 1 every step.
# ---------------------------------------------------------------------------


def test_default_wind_path_single_iteration():
    """With the 2-arg default wind (canopy doesn't feed back), the
    fixed-point loop is a single pass — ``n_wind_iterations == 1`` on
    every step regardless of how many branches are cut."""
    stats_by_gen: list[TreeStats] = []

    def on_step(gen, tree, stats):
        stats_by_gen.append(stats)

    grow_tree(TreeConfig(), n_generations=40, seed=42, on_step=on_step)
    assert all(s.n_wind_iterations == 1 for s in stats_by_gen)


def test_default_wind_forest_single_iteration():
    """Same as above for the Forest path."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=30.0, n_trees_init=10, n_trees_max=200),
    )
    forest = Forest(cfg, seed=42)
    stats = [forest.step(g) for g in range(20)]
    assert all(s.n_wind_iterations == 1 for s in stats)


# ---------------------------------------------------------------------------
# 2. Calm canopy-aware wind converges in one iteration (nothing to cut).
# ---------------------------------------------------------------------------


def test_canopy_aware_calm_converges_in_one_iter():
    """A canopy-aware wind that's tiny in magnitude cuts nothing, so the
    loop exits after a single iteration."""
    seed_tree = make_seed_tree(TreeConfig())
    wind_fn = CanopyIntensifyingWind(base=1e-6, n0=200)

    stats_by_gen: list[TreeStats] = []

    def on_step(gen, tree, stats):
        stats_by_gen.append(stats)

    grow_tree(TreeConfig(), n_generations=15, seed=0, wind_fn=wind_fn, on_step=on_step)
    assert all(s.n_wind_iterations == 1 for s in stats_by_gen)
    assert all(s.n_pruned == 0 for s in stats_by_gen)
    del seed_tree  # quieten unused-var warnings while keeping imports honest


# ---------------------------------------------------------------------------
# 3. Storm with canopy-aware wind: loop iterates and converges within cap.
# ---------------------------------------------------------------------------


def test_grow_tree_storm_iterates_within_cap():
    """A strong canopy-aware wind triggers the loop. Every step's
    ``n_wind_iterations`` must stay ``<= max_pruning_iterations``, and
    at least one generation should exceed 1."""
    max_iter = 8
    cfg = Config(
        tree=TreeConfig(),
        wind=WindConfig(max_pruning_iterations=max_iter, wind_convergence_eps_rel=0.0),
        n_generations=40,
    )
    # Choose ``base`` high enough that pruning fires for several generations.
    wind_fn = CanopyIntensifyingWind(base=8.0, n0=200, slope=0.5)

    stats_by_gen: list[TreeStats] = []

    def on_step(gen, tree, stats):
        stats_by_gen.append(stats)

    grow_tree(cfg, seed=1, wind_fn=wind_fn, on_step=on_step)
    iters = [s.n_wind_iterations for s in stats_by_gen]
    pruned = [s.n_pruned for s in stats_by_gen]
    assert max(iters) >= 2, (
        "expected at least one storm gen with >1 inner iteration; "
        f"got iters={iters}, pruned={pruned}"
    )
    assert max(iters) <= max_iter
    # On gens where the loop iterated, at least one branch was cut.
    for it, p in zip(iters, pruned, strict=True):
        if it > 1:
            assert p > 0


def test_forest_storm_iterates_within_cap():
    """Forest variant of the same check."""
    max_iter = 8
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=30.0, n_trees_init=10, n_trees_max=200),
        wind=WindConfig(max_pruning_iterations=max_iter, wind_convergence_eps_rel=0.0),
    )
    forest = Forest(cfg, seed=1, wind_fn=CanopyIntensifyingWind(base=6.0, n0=500))
    stats = [forest.step(g) for g in range(30)]
    iters = [s.n_wind_iterations for s in stats]
    assert max(iters) >= 2
    assert max(iters) <= max_iter


# ---------------------------------------------------------------------------
# 4. cap=1 + eps=0 recovers the old single-pass behaviour.
# ---------------------------------------------------------------------------


def test_cap_one_recovers_single_pass_behaviour():
    """With max_pruning_iterations=1, every step does exactly one
    wind/prune pass even under a canopy-aware wind. Equivalent to the
    pre-Step-24 behaviour."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=30.0, n_trees_init=10, n_trees_max=200),
        wind=WindConfig(max_pruning_iterations=1, wind_convergence_eps_rel=0.0),
    )
    forest = Forest(cfg, seed=1, wind_fn=CanopyIntensifyingWind(base=6.0, n0=500))
    stats = [forest.step(g) for g in range(20)]
    assert all(s.n_wind_iterations == 1 for s in stats)


# ---------------------------------------------------------------------------
# 5. Higher cap cuts more (canopy-intensifying wind, fixed point < single pass).
# ---------------------------------------------------------------------------


def test_higher_cap_cuts_more_under_canopy_intensifying_wind():
    """Under a wind that intensifies as the canopy thins (bulk-thinning-like),
    the fixed-point loop settles at a smaller canopy than the single-pass
    because the surviving branches feel a higher wind than the original
    pre-pruning estimate. Verify the direction is consistent."""
    base = 6.0

    def _run(cap):
        cfg = Config(
            tree=TreeConfig(),
            forest=ForestConfig(size=30.0, n_trees_init=10, n_trees_max=200),
            wind=WindConfig(max_pruning_iterations=cap, wind_convergence_eps_rel=0.0),
        )
        forest = Forest(cfg, seed=2, wind_fn=CanopyIntensifyingWind(base=base, n0=500))
        cuts = 0
        for g in range(25):
            stats = forest.step(g)
            cuts += stats.n_pruned_total
        return cuts, forest

    cuts_single, _ = _run(cap=1)
    cuts_fixed, _ = _run(cap=8)
    assert cuts_fixed >= cuts_single, (
        f"expected the fixed-point loop to cut at least as many branches as "
        f"the single pass under intensifying wind, got fixed={cuts_fixed} vs "
        f"single={cuts_single}"
    )


# ---------------------------------------------------------------------------
# 6. cap-hit warning fires once per run().
# ---------------------------------------------------------------------------


class _MonotonicallyEscalatingWind:
    """Returns a wind that grows on every call within a generation.

    Each iteration's wind is strictly stronger than the previous one's,
    so each iteration's prune sweep is essentially guaranteed to cut
    something the previous sweep didn't — letting us exercise the
    cap-hit warning machinery deterministically.
    """

    def __init__(self, base: float, ramp: float = 2.0) -> None:
        self.base = base
        self.ramp = ramp
        self.last_gen = -1
        self.in_gen_calls = 0

    def __call__(self, gen, rng, context) -> tuple[float, float, float]:
        if gen != self.last_gen:
            self.last_gen = gen
            self.in_gen_calls = 0
        self.in_gen_calls += 1
        amp = self.base * (self.ramp ** (self.in_gen_calls - 1))
        return (amp, 0.0, 0.0)


def test_cap_hit_warns_once_per_run():
    """Force the cap to be hit by setting max_pruning_iterations very low
    against a monotonically escalating per-iteration wind, so every
    iteration is guaranteed to cut more. Only one RuntimeWarning per
    Forest.run call."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=30.0, n_trees_init=20, n_trees_max=400),
        wind=WindConfig(max_pruning_iterations=2, wind_convergence_eps_rel=0.0),
    )
    # Warm the canopy with default wind so there are many branches to
    # prune when the storm hits — without this the tiny seedlings cut
    # to extinction in iter 1 and iter 2 cuts nothing (no warning).
    warm = Forest(cfg, seed=3)
    for g in range(40):
        warm.step(g)
    # Drop in the escalating wind and run a bit more.
    warm.wind_fn = _MonotonicallyEscalatingWind(base=2.0, ramp=2.0)
    warm._wind_arity = 3  # tell the loop this wind is canopy-aware
    warm._cap_hit_warned = False
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        # Use run() so the cap-hit gate resets and we exercise the
        # once-per-run semantics.
        warm.run(10)
    cap_warnings = [w for w in caught if "fixed-point loop hit cap" in str(w.message)]
    assert len(cap_warnings) == 1, (
        f"expected exactly one cap warning per run(), got {len(cap_warnings)}"
    )


def test_cap_eq_one_never_warns():
    """Setting max_pruning_iterations=1 is the explicit single-pass knob,
    not a cap violation — must never emit the cap-hit warning even
    when iter 1 cuts a lot."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=30.0, n_trees_init=10, n_trees_max=200),
        wind=WindConfig(max_pruning_iterations=1, wind_convergence_eps_rel=0.0),
    )
    forest = Forest(cfg, seed=4, wind_fn=CanopyIntensifyingWind(base=6.0, n0=500))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        forest.run(20)
    cap_warnings = [w for w in caught if "fixed-point loop hit cap" in str(w.message)]
    assert cap_warnings == []


# ---------------------------------------------------------------------------
# 7. ε-tolerance early-exit triggers when wind change is small.
# ---------------------------------------------------------------------------


class _StrongThenEpsTrigger:
    """First call within a generation returns a strong wind that prunes;
    subsequent calls return a tiny perturbation so the ε early-exit fires
    iteration 2 before the second prune sweep.

    Resets ``in_gen_calls`` whenever the generation index advances —
    state is per-gen, not across the whole run, so the trigger pattern
    repeats every generation.
    """

    def __init__(self, base: float) -> None:
        self.base = base
        self.last_gen = -1
        self.in_gen_calls = 0
        self.total_calls = 0

    def __call__(self, gen, rng, context) -> tuple[float, float, float]:
        if gen != self.last_gen:
            self.last_gen = gen
            self.in_gen_calls = 0
        self.in_gen_calls += 1
        self.total_calls += 1
        if self.in_gen_calls == 1:
            return (self.base, 0.0, 0.0)
        return (self.base + 0.0001, 0.0, 0.0)


def test_eps_early_exit_when_wind_change_small():
    """Inject a per-gen wind whose iter-2 value barely differs from
    iter-1. Once the tree has branches and pruning fires, the loop
    must exit at iteration 2 (the ε check) without running iter-2's
    prune sweep."""
    wind_fn = _StrongThenEpsTrigger(base=12.0)
    cfg = Config(
        tree=TreeConfig(),
        wind=WindConfig(max_pruning_iterations=8, wind_convergence_eps_rel=0.01),
        n_generations=20,
    )

    seen: list[TreeStats] = []

    def on_step(gen, tree, stats):
        seen.append(stats)

    grow_tree(cfg, seed=4, wind_fn=wind_fn, on_step=on_step)
    # Find every gen where pruning happened — those are the gens the
    # loop should have iterated. With the ε trigger, each such gen
    # exits at exactly 2 iterations.
    storm_gens = [s for s in seen if s.n_pruned > 0]
    assert storm_gens, "expected at least one generation with pruning under base=12"
    assert all(s.n_wind_iterations == 2 for s in storm_gens), (
        f"expected n_wind_iterations==2 on every storm gen, got "
        f"{[(s.generation, s.n_wind_iterations, s.n_pruned) for s in storm_gens]}"
    )


def test_eps_zero_never_short_circuits_above_zero_cuts():
    """With ``eps_rel = 0`` the loop never short-circuits on wind delta;
    it only exits on ``n_cut == 0`` or the cap. Compare against the
    same wind with the default ε early-exit: ε-off must use at least
    as many iterations on every gen, and a fully-stable-wind scenario
    is enough to demonstrate the difference."""

    def _run(eps):
        wind_fn = _StrongThenEpsTrigger(base=12.0)
        cfg = Config(
            tree=TreeConfig(),
            wind=WindConfig(max_pruning_iterations=8, wind_convergence_eps_rel=eps),
            n_generations=20,
        )
        seen: list[TreeStats] = []

        def on_step(gen, tree, stats):
            seen.append(stats)

        grow_tree(cfg, seed=4, wind_fn=wind_fn, on_step=on_step)
        return seen

    with_eps = _run(0.01)
    without_eps = _run(0.0)
    # ε-on must never iterate more than ε-off (this is the contract).
    # The reverse — ε-off iterating *strictly more* on some gen — is
    # not guaranteed because under nearly-constant wind iter 2's prune
    # often cuts zero by Weibull, exiting on the n_cut==0 branch.
    pairs = list(zip(with_eps, without_eps, strict=True))
    assert all(a.n_wind_iterations <= b.n_wind_iterations for a, b in pairs), (
        "ε-on (0.01) must never iterate more than ε-off (0.0); got "
        f"{[(a.n_wind_iterations, b.n_wind_iterations) for a, b in pairs]}"
    )


# ---------------------------------------------------------------------------
# 8. Direct unit test of the helper.
# ---------------------------------------------------------------------------


def test_prune_to_fixed_point_default_wind_returns_one_iter():
    """``_prune_to_fixed_point`` short-circuits for a 2-arg wind_fn."""
    tree = make_seed_tree(TreeConfig())
    # Grow a tiny canopy via primary growth so prune has something to act on.
    grow_tree(TreeConfig(), n_generations=10, seed=0)
    rng = np.random.default_rng(0)

    wind, n_pruned, n_iter, cap_hit = _prune_to_fixed_point(
        tree,
        default_wind_fn,
        wind_arity=2,
        gen=0,
        rng=rng,
        leaf_drag_S0=TreeConfig().leaf_surface,
        cauchy=TreeConfig().cauchy,
        max_iterations=8,
        eps_rel=0.01,
    )
    assert n_iter == 1
    assert cap_hit is False
    assert isinstance(wind, tuple) and len(wind) == 3
    assert n_pruned >= 0


# ---------------------------------------------------------------------------
# 9. Wall-clock overhead guard (soft) at 10-tree scale.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_step_overhead_at_10_trees_within_budget():
    """Sanity check: at the small forest scale the fixed-point loop
    should not blow up step time by more than a factor of ``3``
    vs. the cap=1 baseline. Loose threshold so this passes consistently
    across machines but catches accidental algorithmic regressions
    (e.g. dropping the eps early-exit, or pulling forest_to_cylinders
    back into the per-branch Python loop).
    """

    def _time_run(cap, eps):
        cfg = Config(
            tree=TreeConfig(),
            forest=ForestConfig(size=30.0, n_trees_init=10, n_trees_max=200),
            wind=WindConfig(max_pruning_iterations=cap, wind_convergence_eps_rel=eps),
        )
        forest = Forest(cfg, seed=4, wind_fn=CanopyIntensifyingWind(base=6.0, n0=500))
        # Warm up to a representative steady-state size before timing.
        for g in range(10):
            forest.step(g)
        t0 = time.perf_counter()
        for g in range(10, 40):
            forest.step(g)
        return time.perf_counter() - t0

    baseline = _time_run(cap=1, eps=0.0)
    loop = _time_run(cap=8, eps=0.01)
    ratio = loop / max(baseline, 1e-6)
    assert ratio < 3.0, (
        f"overhead ratio too high: {ratio:.2f}x (baseline={baseline:.3f}s, loop={loop:.3f}s)"
    )
