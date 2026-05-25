"""Tests for the Forest container (Step 12)."""

import pytest

from mechatree.config import Config, ForestConfig, TreeConfig
from mechatree.forest import Forest, ForestStats
from mechatree.genome import ConstantAllocation, ConstantSafety
from mechatree.light import Sun


def _small_config(
    n_trees_init: int = 4,
    n_trees_max: int = 50,
    size: float = 20.0,
    max_age: int = 1000,
    min_age_for_undersize: int = 5,
    min_branches: int = 11,
) -> Config:
    return Config(
        tree=TreeConfig(),
        forest=ForestConfig(
            size=size,
            n_trees_init=n_trees_init,
            n_trees_max=n_trees_max,
            max_age=max_age,
            min_age_for_undersize=min_age_for_undersize,
            min_branches=min_branches,
        ),
    )


# ---------------------------------------------------------------------------
# ForestConfig
# ---------------------------------------------------------------------------


def test_forest_config_defaults():
    fc = ForestConfig()
    assert fc.size == 100.0
    assert fc.n_trees_init == 100
    assert fc.n_trees_max == 10000
    # Fortran death rule defaults:
    assert fc.min_branches == 11
    assert fc.min_age_for_undersize == 5
    assert fc.max_age == 1000


def test_forest_config_validation():
    with pytest.raises(ValueError):
        ForestConfig(size=0.0)
    with pytest.raises(ValueError):
        ForestConfig(n_trees_init=0)
    with pytest.raises(ValueError):
        ForestConfig(n_trees_init=200, n_trees_max=100)
    with pytest.raises(ValueError):
        ForestConfig(max_age=0)


# ---------------------------------------------------------------------------
# Initial population
# ---------------------------------------------------------------------------


def test_initial_population_has_correct_size():
    cfg = _small_config(n_trees_init=12)
    forest = Forest(cfg, seed=42)
    assert len(forest.trees) == 12
    assert len(forest.ages) == 12
    assert all(a == 0 for a in forest.ages)


def test_initial_population_inside_disk():
    """Every initial trunk lands at radius < forest.size."""
    cfg = _small_config(n_trees_init=50, size=10.0)
    forest = Forest(cfg, seed=7)
    radius_sq = cfg.forest.size**2
    for tree in forest.trees:
        x, y, z = tree.get_location(0)
        assert x * x + y * y <= radius_sq
        assert z == 0.0


def test_initial_population_has_per_tree_orientation():
    """Tree unit_b should NOT all be (1, 0, 0); each tree picks its own angle."""
    cfg = _small_config(n_trees_init=10)
    forest = Forest(cfg, seed=11)
    unit_bs = [tree.get_unit_b(0) for tree in forest.trees]
    # With 10 random angles, we expect different unit_b vectors.
    distinct = {tuple(round(c, 4) for c in ub) for ub in unit_bs}
    assert len(distinct) >= 5  # plenty of variety


def test_seed_reproducibility():
    """Same seed -> same initial trunk positions."""
    cfg = _small_config(n_trees_init=8)
    f1 = Forest(cfg, seed=42)
    f2 = Forest(cfg, seed=42)
    locs1 = sorted(f1.trees[i].get_location(0) for i in range(len(f1.trees)))
    locs2 = sorted(f2.trees[i].get_location(0) for i in range(len(f2.trees)))
    assert locs1 == locs2


# ---------------------------------------------------------------------------
# Single step + death + birth
# ---------------------------------------------------------------------------


def test_step_returns_stats():
    cfg = _small_config(n_trees_init=3)
    forest = Forest(cfg, seed=42)
    stats = forest.step(0)
    assert isinstance(stats, ForestStats)
    assert stats.generation == 0
    assert stats.n_trees == 3  # nothing should die at gen 0
    assert stats.n_branches_total >= 3
    assert stats.biomass_total > 0.0


def test_step_increments_ages():
    cfg = _small_config(n_trees_init=3)
    forest = Forest(cfg, seed=42)
    forest.step(0)
    assert all(a == 1 for a in forest.ages)
    forest.step(1)
    assert all(a == 2 for a in forest.ages)


def test_old_trees_die():
    """If max_age is very small, every initial tree should die quickly."""
    cfg = _small_config(n_trees_init=5, max_age=2, n_trees_max=200)
    forest = Forest(cfg, seed=11)
    # After 3 steps (age=3 > max_age=2), the initial trees die.
    # Any new seedlings born in the meantime may survive (their age is younger).
    n_initial_alive_after = []
    initial_ids = [id(t) for t in forest.trees]
    for gen in range(3):
        forest.step(gen)
        alive_initial = sum(1 for t in forest.trees if id(t) in initial_ids)
        n_initial_alive_after.append(alive_initial)
    assert n_initial_alive_after[-1] == 0  # all initial trees died of old age


def test_undersized_trees_die_after_age_threshold():
    """A tree below min_branches that ages past min_age_for_undersize dies."""
    # Use very small min_branches so initial single-trunk trees qualify.
    cfg = _small_config(n_trees_init=3, min_branches=100, min_age_for_undersize=2, max_age=1000)
    forest = Forest(cfg, seed=0)
    for gen in range(8):
        forest.step(gen)
    # At very high min_branches, no tree gets big enough; all die when
    # their age crosses min_age_for_undersize (with possible replenishment
    # via birth — but newborns are also undersized, so churn continues).
    # The forest may not be empty due to ongoing births, but the original
    # initials are definitely gone.
    assert all(a <= cfg.forest.min_age_for_undersize for a in forest.ages)


def test_no_growth_no_births():
    """With ConstantAllocation(p_seeds=0, p_leaves=0), primary_growth makes
    no new branches AND no seedlings. Tree count stays at the initial."""
    cfg = _small_config(n_trees_init=5, max_age=1000, min_branches=1)
    forest = Forest(
        cfg,
        seed=42,
        safety=ConstantSafety(0.0),
        allocation=ConstantAllocation(p_seeds=0.0, p_leaves=0.0, phototropism=0.0),
    )
    for gen in range(5):
        forest.step(gen)
    assert len(forest.trees) == 5


def test_seedlings_land_inside_disk():
    """Every seedling's trunk position must be inside the forest disk."""
    cfg = _small_config(n_trees_init=5, size=15.0, n_trees_max=200)
    forest = Forest(cfg, seed=42)
    radius_sq = cfg.forest.size**2
    for gen in range(8):
        forest.step(gen)
    for tree in forest.trees:
        x, y, _ = tree.get_location(0)
        assert x * x + y * y <= radius_sq + 1e-9


def test_n_trees_max_respected():
    """No matter how exuberant the birthing, the forest never exceeds n_trees_max."""
    cfg = _small_config(n_trees_init=3, n_trees_max=10, size=30.0)
    forest = Forest(cfg, seed=42)
    for gen in range(20):
        forest.step(gen)
        assert len(forest.trees) <= cfg.forest.n_trees_max


# ---------------------------------------------------------------------------
# Light competition
# ---------------------------------------------------------------------------


def test_close_trees_shade_each_other():
    """Two trees right next to each other should sum-aggregate less light than
    two trees far apart. We measure the SUM of trunk light values."""
    sun = Sun.from_arrays(elev=[0.0], azim=[0.0])  # vertical sun
    safety = ConstantSafety(0.0)
    allocation = ConstantAllocation(0.0, 0.0, 0.0)  # no growth

    # Close — within a single shadow cell (size_leaf default = 1.0).
    cfg_close = _small_config(n_trees_init=2, size=2.0)
    f_close = Forest(cfg_close, seed=0, safety=safety, allocation=allocation, sun=sun)
    # Manually pin both trees to (0, 0) to force shading.
    f_close.trees[0].set_location(0, (0.0, 0.0, 0.0))
    f_close.trees[1].set_location(0, (0.0, 0.0, 0.0))
    f_close.step(0)
    light_close = sum(t.get_light(0) for t in f_close.trees)

    # Far apart — far enough that they fall into separate cells.
    f_far = Forest(cfg_close, seed=0, safety=safety, allocation=allocation, sun=sun)
    f_far.trees[0].set_location(0, (-50.0, 0.0, 0.0))
    f_far.trees[1].set_location(0, (50.0, 0.0, 0.0))
    f_far.step(0)
    light_far = sum(t.get_light(0) for t in f_far.trees)

    assert light_far > light_close
    assert light_far == pytest.approx(2.0)  # both fully lit
    assert light_close < 2.0  # one shadows the other


# ---------------------------------------------------------------------------
# run() + on_step callback
# ---------------------------------------------------------------------------


def test_run_calls_callback_each_generation():
    cfg = _small_config(n_trees_init=3)
    forest = Forest(cfg, seed=42)
    history = []

    def cb(gen, f, stats):
        history.append((gen, stats.n_trees))

    forest.run(5, on_step=cb)
    assert [h[0] for h in history] == [0, 1, 2, 3, 4]


def test_run_progresses_state():
    """After ``run(N)``, the forest is in the state it would be after N steps."""
    cfg = _small_config(n_trees_init=3)
    f_a = Forest(cfg, seed=42)
    f_b = Forest(cfg, seed=42)
    f_a.run(5)
    for gen in range(5):
        f_b.step(gen)
    assert len(f_a.trees) == len(f_b.trees)
    locs_a = sorted(t.get_location(0) for t in f_a.trees)
    locs_b = sorted(t.get_location(0) for t in f_b.trees)
    assert locs_a == locs_b


# ---------------------------------------------------------------------------
# Wind direction
# ---------------------------------------------------------------------------


def test_zero_wind_no_pruning():
    """Zero wind ⇒ stress = 0 on every branch ⇒ no pruning *ever*.

    The previous version compared total branches across two stochastic
    runs and asserted a direction. Both runs share ``seed=42`` but
    diverge because the default storm ``wind_fn`` consumes ``rng.random``
    calls that ``zero_wind`` doesn't — so the rest of the simulation
    (births, dispersal, deaths) takes different RNG paths. With a small
    forest the resulting signal was smaller than cross-platform
    floating-point drift in the C++ random core, flipping the
    comparison on Linux / Windows builds. We now check the direct
    invariant via ``ForestStats.n_pruned_total``.
    """
    cfg = _small_config(n_trees_init=8, n_trees_max=100, size=30.0, max_age=1000)

    def zero_wind(gen, rng):
        return (0.0, 0.0, 0.0)

    f_quiet = Forest(cfg, seed=42, wind_fn=zero_wind)
    quiet_pruned = 0
    for gen in range(20):
        stats = f_quiet.step(gen)
        quiet_pruned += stats.n_pruned_total
    assert quiet_pruned == 0, f"expected no pruning under zero wind, got {quiet_pruned}"


def test_strong_wind_prunes():
    """Sanity counterpart to :func:`test_zero_wind_no_pruning`: a strong
    deterministic wind must prune *something* within a short run.

    Uses a fixed-amplitude (5.0) westerly so the test doesn't depend on
    the stochastic default storm wind — that path can leave a small
    forest untouched for ~30 gens by chance.
    """
    cfg = _small_config(n_trees_init=8, n_trees_max=100, size=30.0, max_age=1000)

    def strong_wind(gen, rng):
        return (5.0, 0.0, 0.0)

    f = Forest(cfg, seed=42, wind_fn=strong_wind)
    pruned_total = 0
    for gen in range(20):
        pruned_total += f.step(gen).n_pruned_total
    assert pruned_total > 0, f"expected strong wind to prune something, got {pruned_total}"


# ---------------------------------------------------------------------------
# Empty forest gracefully degenerates
# ---------------------------------------------------------------------------


def test_empty_forest_step_no_crash():
    """If every tree dies, subsequent steps still run without crashing."""
    # max_age=1 + zero allocation means every tree dies before reproducing.
    cfg = _small_config(n_trees_init=2, max_age=1)
    forest = Forest(
        cfg,
        seed=42,
        safety=ConstantSafety(0.0),
        allocation=ConstantAllocation(0.0, 0.0, 0.0),
    )
    # Step 0: ages become 1, not > 1 yet, all survive.
    forest.step(0)
    # Step 1: ages become 2 > 1, all die. No births (zero allocation).
    forest.step(1)
    assert len(forest.trees) == 0
    # Subsequent step on an empty forest must not crash.
    stats = forest.step(2)
    assert stats.n_trees == 0
    assert stats.biomass_total == 0.0
