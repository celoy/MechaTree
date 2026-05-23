"""Growth tests — requested, secondary and primary growth on small trees."""

import math

import pytest

from mechatree import PyTree
from mechatree.genome import ConstantAllocation, ConstantSafety
from mechatree.growth import primary_growth, requested_growth, secondary_growth


def _vertical_trunk(length=1.0, diameter=0.1):
    t = PyTree({})
    t.set_length(0, length)
    t.set_diameter(0, diameter)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    t.set_location(0, (0.0, 0.0, 0.0))
    return t


def test_requested_growth_hand_derived_single_branch():
    """One trunk, ConstantSafety(1.0), max_stress=8.

    vol_actual = pi/4 * d² * L = pi * 0.01 / 4 ≈ 7.854e-3
    maintenance_vol = pi * d * h = pi * 0.1 * 0.02 ≈ 6.283e-3
    vol_wished = 1.0 * vol_actual * 8^(2/3) = vol_actual * 4
    vol_growth = (vol_wished - vol_actual) + maintenance_vol = 3*vol_actual + 6.283e-3
    """
    t = _vertical_trunk()
    t.set_max_stress(0, 8.0)
    t.reorder()  # ensures nb_leaves=1

    requested_growth(t, ConstantSafety(1.0), maintenance_h=0.02)

    vol_actual = 0.25 * math.pi * 0.1**2 * 1.0
    maintenance_vol = math.pi * 0.1 * 0.02
    expected_vol_growth = 3 * vol_actual + maintenance_vol

    assert t.get_maintenance_vol(0) == pytest.approx(maintenance_vol)
    assert t.get_vol_growth(0) == pytest.approx(expected_vol_growth)
    assert t.get_vol_summed(0) == pytest.approx(expected_vol_growth)


def test_requested_growth_zero_stress_only_maintenance():
    """max_stress = 0 ⇒ vol_wished = 0 ⇒ vol_growth = maintenance only."""
    t = _vertical_trunk()
    t.set_max_stress(0, 0.0)
    t.reorder()
    requested_growth(t, ConstantSafety(1.0), maintenance_h=0.02)
    assert t.get_vol_growth(0) == pytest.approx(math.pi * 0.1 * 0.02)


def test_secondary_growth_grows_diameter_under_light():
    """A single leaf with light=1 and a populated vol_summed grows its
    diameter from 0.1 toward a hand-derivable value."""
    t = _vertical_trunk()
    t.set_max_stress(0, 8.0)
    t.reorder()
    requested_growth(t, ConstantSafety(1.0), maintenance_h=0.02)
    t.set_light(0, 1.0)

    d0 = t.get_diameter(0)
    secondary_growth(t, volume_per_leaf=0.5)
    d1 = t.get_diameter(0)
    assert d1 > d0  # photosynthate flowed into diameter


def test_secondary_growth_excess_feeds_reserve():
    """When light * vol_per_leaf > vol_summed, the excess goes to reserve."""
    t = _vertical_trunk()
    t.set_vol_summed(0, 0.01)
    t.set_vol_growth(0, 0.01)
    t.set_maintenance_vol(0, 0.0)
    t.set_light(0, 1.0)
    t.set_nb_leaves(0, 1)
    t.set_reserve(0.0)

    # photosynth = 1.0 * 0.5 = 0.5; vol_growth_branches = min(0.5, 0.01) = 0.01
    # reserve gains 0.5 - 0.01 = 0.49
    secondary_growth(t, volume_per_leaf=0.5)
    assert t.get_reserve() == pytest.approx(0.49)


def test_primary_growth_spawns_two_children():
    """Single trunk with ample reserve and p_leaves=1.0 grows two daughters."""
    t = _vertical_trunk()
    t.set_light(0, 1.0)
    t.set_reserve(10.0)
    t.set_nb_leaves(0, 1)
    t.set_seed(42)

    created = primary_growth(
        t,
        ConstantAllocation(p_seeds=0.0, p_leaves=1.0, phototropism=0.0),
        twig_length=0.5,
        twig_diameter=0.05,
        theta1=0.4,
        theta2=-0.4,
        gamma1=0.0,
        gamma2=math.pi,
        generation=0,
    )
    assert created == 2
    assert t.get_number_of_branches() == 3

    # Both children sit at the trunk's tip — (0,0,1).
    for child_idx in (1, 2):
        assert t.get_location(child_idx) == pytest.approx((0.0, 0.0, 1.0))
        assert t.get_length(child_idx) == pytest.approx(0.5)
        assert t.get_diameter(child_idx) == pytest.approx(0.05)

    # Mother's light is cleared since it's no longer a leaf.
    assert t.get_light(0) == pytest.approx(0.0)
    # Reserve depleted by 2 * VolumeTwig.
    volume_twig = 0.25 * math.pi * 0.5 * 0.05**2
    assert t.get_reserve() == pytest.approx(10.0 - 2 * volume_twig)


def test_primary_growth_zero_reserve_no_spawn():
    """No reserve ⇒ nothing happens."""
    t = _vertical_trunk()
    t.set_reserve(0.0)
    t.set_light(0, 1.0)
    t.set_nb_leaves(0, 1)
    t.set_seed(0)
    created = primary_growth(
        t,
        ConstantAllocation(p_seeds=0.0, p_leaves=1.0, phototropism=0.0),
        twig_length=0.5,
        twig_diameter=0.05,
        theta1=0.4,
        theta2=-0.4,
        gamma1=0.0,
        gamma2=math.pi,
        generation=0,
    )
    assert created == 0
    assert t.get_number_of_branches() == 1


def test_primary_growth_underground_leaf_no_spawn():
    """A leaf whose tip is below ground (z constraint) does not spawn."""
    t = PyTree({})
    t.set_length(0, 1.0)
    t.set_diameter(0, 0.1)
    t.set_location(0, (0.0, 0.0, 0.0))
    t.set_unit_t(0, (0.0, 0.0, -1.0))  # pointing DOWN — tip_z = -1
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    t.set_light(0, 1.0)
    t.set_reserve(10.0)
    t.set_nb_leaves(0, 1)
    t.set_seed(42)

    created = primary_growth(
        t,
        ConstantAllocation(p_seeds=0.0, p_leaves=1.0, phototropism=0.0),
        twig_length=0.5,
        twig_diameter=0.05,
        theta1=0.4,
        theta2=-0.4,
        gamma1=0.0,
        gamma2=math.pi,
        generation=0,
    )
    assert created == 0


def test_primary_growth_seed_reproducibility():
    """Same seed ⇒ identical child geometry."""

    def grow(seed):
        t = _vertical_trunk()
        t.set_light(0, 1.0)
        t.set_reserve(10.0)
        t.set_nb_leaves(0, 1)
        t.set_seed(seed)
        primary_growth(
            t,
            ConstantAllocation(p_seeds=0.0, p_leaves=1.0, phototropism=0.0),
            twig_length=0.5,
            twig_diameter=0.05,
            theta1=0.3,
            theta2=-0.3,
            gamma1=0.0,
            gamma2=math.pi,
            generation=0,
        )
        return [t.get_unit_t(i) for i in range(1, t.get_number_of_branches())]

    assert grow(seed=123) == grow(seed=123)
    assert grow(seed=123) != grow(seed=456)
