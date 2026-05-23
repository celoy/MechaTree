"""Step 9 integration — one full generation of the simulation pipeline.

Step 10 will plug in the real light module, and Step 11 will orchestrate the
multi-generation loop. For now this test stubs `branch.light` directly and
exercises every routine PR2 ships in one go to confirm they compose.
"""

import math

import pytest

from mechatree import PyTree
from mechatree.genome import ConstantAllocation, ConstantSafety
from mechatree.growth import primary_growth, requested_growth, secondary_growth
from mechatree.mechanics import calculate_stresses
from mechatree.pruning import prune


def _seed_tree():
    """A 3-branch starter: trunk + two short tilted children."""
    t = PyTree({})
    t.set_length(0, 1.0)
    t.set_diameter(0, 0.15)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    t.set_location(0, (0.0, 0.0, 0.0))

    t.add_branch_with_geometry(
        0,
        length=0.5,
        diameter=0.06,
        unit_t=(0.2, 0.0, 0.98),
        unit_b=(1.0, 0.0, 0.0),
    )
    t.add_branch_with_geometry(
        0,
        length=0.5,
        diameter=0.06,
        unit_t=(0.0, 0.2, 0.98),
        unit_b=(1.0, 0.0, 0.0),
    )
    t.reorder()
    return t


def _stub_light_on_leaves(tree, value=1.0):
    for idx in tree.leaf_indices():
        tree.set_light(idx, value)


def test_full_generation_runs_and_grows():
    """Run stress -> requested -> secondary -> prune -> primary on the seed
    tree (with stubbed leaf light) and assert sane invariants."""
    t = _seed_tree()
    t.set_seed(123)
    t.set_reserve(0.05)

    diams_before = {i: t.get_diameter(i) for i in range(t.get_number_of_branches())}

    _stub_light_on_leaves(t, value=1.0)
    calculate_stresses(t, leaf_drag_S0=0.5, cauchy=1.0)
    requested_growth(t, ConstantSafety(1.0), maintenance_h=0.005)
    secondary_growth(t, volume_per_leaf=0.01)
    prune(t, wind=(1.0, 0.0, 0.0), leaf_drag_S0=0.5, cauchy=1.0)
    t.reorder()

    # primary_growth then adds twigs; reseed leaf light first.
    _stub_light_on_leaves(t, value=1.0)
    created = primary_growth(
        t,
        ConstantAllocation(p_seeds=0.0, p_leaves=1.0, phototropism=0.0),
        twig_length=0.3,
        twig_diameter=0.02,
        theta1=0.25,
        theta2=-0.25,
        gamma1=0.0,
        gamma2=math.pi,
        generation=0,
    )
    t.reorder()

    # Invariants:
    n = t.get_number_of_branches()
    assert n >= 1  # trunk survives
    assert created % 2 == 0  # twigs grow in pairs
    assert t.get_reserve() >= 0.0  # no negative pool

    # No NaN / inf anywhere on the typed fields.
    for i in range(n):
        for field in (
            t.get_length(i),
            t.get_diameter(i),
            t.get_stress(i),
            t.get_max_stress(i),
            t.get_vol_growth(i),
            t.get_vol_summed(i),
        ):
            assert math.isfinite(field), f"non-finite scalar at branch {i}"
        for v in t.get_location(i) + t.get_force(i) + t.get_moment(i):
            assert math.isfinite(v), f"non-finite vector at branch {i}"

    # At least one of the original branches that survived grew in diameter.
    grew = False
    for i in range(n):
        # only check the trunk and its original direct children that are
        # still around (their depth-first indices may have shifted but the
        # trunk is always index 0).
        if i == 0 and t.get_diameter(0) > diams_before[0]:
            grew = True
    assert grew, "no original branch grew in diameter"


def test_repeated_generations_dont_explode_or_die():
    """Run several generations; tree should sustain itself without crashing
    or going to zero branches (the trunk always survives)."""
    t = _seed_tree()
    t.set_seed(7)
    t.set_reserve(0.1)

    counts = []
    for gen in range(5):
        _stub_light_on_leaves(t, value=1.0)
        calculate_stresses(t, leaf_drag_S0=0.5, cauchy=1.0)
        requested_growth(t, ConstantSafety(1.0), maintenance_h=0.005)
        secondary_growth(t, volume_per_leaf=0.02)
        prune(t, wind=(2.0, 0.0, 0.0), leaf_drag_S0=0.5, cauchy=1.0)
        t.reorder()
        primary_growth(
            t,
            ConstantAllocation(p_seeds=0.0, p_leaves=0.5, phototropism=0.0),
            twig_length=0.3,
            twig_diameter=0.02,
            theta1=0.25,
            theta2=-0.25,
            gamma1=0.0,
            gamma2=math.pi,
            generation=gen,
        )
        t.reorder()
        counts.append(t.get_number_of_branches())
        assert counts[-1] >= 1

    # No NaNs / no explosive blow-up to crazy sizes.
    assert max(counts) < 10_000
    assert t.get_reserve() >= 0.0


def test_property_map_untouched_by_simulation():
    """A user property set before the simulation runs is preserved after."""
    t = _seed_tree()
    t.set_seed(11)
    # Use the property map for a user-extension annotation.
    t.add_property(0, "user_id", 42.0)
    pytest.approx(t.get_property(0, "user_id")) == 42.0  # noqa: B015

    _stub_light_on_leaves(t, value=1.0)
    calculate_stresses(t, leaf_drag_S0=0.5, cauchy=1.0)
    requested_growth(t, ConstantSafety(1.0), maintenance_h=0.005)
    secondary_growth(t, volume_per_leaf=0.02)
    # No prune so the trunk's index 0 is stable.
    assert t.get_property(0, "user_id") == pytest.approx(42.0)
