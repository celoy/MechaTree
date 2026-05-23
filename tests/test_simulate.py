"""Tests for the single-tree orchestrator (Step 11)."""

import math

import numpy as np
import pytest

from mechatree.config import Config, TreeConfig
from mechatree.genome import ConstantAllocation, ConstantSafety
from mechatree.light import Sun
from mechatree.simulate import default_wind_fn, grow_tree, make_seed_tree


def test_make_seed_tree_initial_state():
    """The seed tree is a single trunk at the origin pointing up, with
    reserve = 2 * volume_twig per the Fortran reference."""
    cfg = TreeConfig(twig_length=1.0, twig_diameter=0.1)
    tree = make_seed_tree(cfg)
    assert tree.get_number_of_branches() == 1
    assert tree.get_length(0) == pytest.approx(1.0)
    assert tree.get_diameter(0) == pytest.approx(0.1)
    assert tree.get_unit_t(0) == pytest.approx((0.0, 0.0, 1.0))
    assert tree.get_unit_b(0) == pytest.approx((1.0, 0.0, 0.0))
    assert tree.get_location(0) == pytest.approx((0.0, 0.0, 0.0))
    assert tree.get_reserve() == pytest.approx(2.0 * cfg.volume_twig)


def test_default_wind_fn_shape_and_amplitude_range():
    """The default wind has z=0, rotates with generation, and amplitude is
    drawn from the Fortran ``0.835 - log(u)/6`` distribution."""
    rng = np.random.default_rng(0)
    wind = default_wind_fn(5, rng)
    assert wind[2] == 0.0
    amplitude = math.hypot(wind[0], wind[1])
    # 0.835 + (positive heavy tail). Hard lower bound 0.835 - log(1)/6 = 0.835.
    assert amplitude >= 0.835
    # Verify rotation: same amplitude at different generations gives different
    # x/y splits.
    rng2 = np.random.default_rng(0)
    wind5 = default_wind_fn(5, rng2)
    rng3 = np.random.default_rng(0)
    wind6 = default_wind_fn(6, rng3)
    assert wind5 != wind6  # different generation -> different (cos, sin)


def test_grow_tree_runs_to_completion_short():
    """A short run completes without crashing, on a default config."""
    cfg = Config()
    tree = grow_tree(cfg, n_generations=5, seed=42)
    assert tree.get_number_of_branches() >= 1


def test_grow_tree_reproducible_with_same_seed():
    """Same seed + same config => byte-identical final state on supported
    invariants (branch count, leaf count, reserve)."""
    cfg = Config()

    def run(seed):
        t = grow_tree(cfg, n_generations=10, seed=seed)
        return (
            t.get_number_of_branches(),
            t.get_total_leaves(),
            round(t.get_reserve(), 6),
        )

    assert run(seed=7) == run(seed=7)
    assert run(seed=7) != run(seed=13)


def test_grow_tree_callback_fires_every_generation():
    cfg = Config()
    seen = []

    def cb(gen, tree):
        seen.append((gen, tree.get_number_of_branches()))

    grow_tree(cfg, n_generations=5, seed=0, on_step=cb)
    assert [s[0] for s in seen] == [0, 1, 2, 3, 4]


def test_grow_tree_tree_config_form_requires_n_generations():
    """Calling with a bare TreeConfig demands ``n_generations`` keyword."""
    with pytest.raises(TypeError):
        grow_tree(TreeConfig(), seed=0)


def test_grow_tree_accepts_tree_config_with_explicit_n_generations():
    tree = grow_tree(TreeConfig(), n_generations=3, seed=0)
    assert tree.get_number_of_branches() >= 1


def test_grow_tree_custom_wind_fn():
    """User-supplied wind function is respected — verify by injecting a
    no-wind function and observing zero pruning (more leaves survive)."""
    cfg = Config()

    def zero_wind(gen, rng):
        return (0.0, 0.0, 0.0)

    tree_quiet = grow_tree(cfg, n_generations=8, seed=42, wind_fn=zero_wind)
    tree_stormy = grow_tree(cfg, n_generations=8, seed=42)  # default wind
    # Under zero wind, no pruning ever happens; the storm-driven default
    # culls many branches occasionally.
    assert tree_quiet.get_number_of_branches() >= tree_stormy.get_number_of_branches()


def test_grow_tree_custom_genome_models():
    """Custom safety / allocation models are accepted and reach the tree."""
    cfg = Config()
    # Safety = 0 means no growth requested above maintenance — the tree
    # should grow MUCH more slowly than with Safety=1.
    tree_growing = grow_tree(
        cfg,
        n_generations=8,
        seed=42,
        safety=ConstantSafety(1.0),
        allocation=ConstantAllocation(p_seeds=0.0, p_leaves=1.0, phototropism=0.0),
    )
    tree_idle = grow_tree(
        cfg,
        n_generations=8,
        seed=42,
        safety=ConstantSafety(0.0),
        allocation=ConstantAllocation(p_seeds=0.0, p_leaves=0.0, phototropism=0.0),
    )
    assert tree_growing.get_number_of_branches() > tree_idle.get_number_of_branches()


def test_grow_tree_custom_sun_is_used():
    """A coarse sun (2 directions) should still drive growth."""
    cfg = Config()
    sun = Sun.from_arrays(elev=[0.1, 0.5], azim=[0.0, math.pi])
    tree = grow_tree(cfg, n_generations=5, seed=42, sun=sun)
    assert tree.get_number_of_branches() >= 1
