"""Tests for the light interception module (Step 10)."""

import math

import numpy as np
import pytest

from mechatree import PyTree
from mechatree.light import Leaves, Sun, aggregate_onto_trees, extract_leaves, intercept

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _vertical_trunk(length=1.0, diameter=0.1, location=(0.0, 0.0, 0.0)):
    t = PyTree({})
    t.set_length(0, length)
    t.set_diameter(0, diameter)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    t.set_location(0, location)
    return t


# ---------------------------------------------------------------------------
# Sun
# ---------------------------------------------------------------------------


def test_sun_default_dimensions():
    sun = Sun()
    assert sun.n_elevations == 4
    assert sun.n_azimuths == 8
    assert sun.n_directions == 32
    assert sun.elev.shape == (32,)
    assert sun.azim.shape == (32,)


def test_sun_elev_matches_fortran_formula():
    """elev(k) = acos((i - 0.5) / N) for i = 1..N — cos-uniform on [0, 1]."""
    sun = Sun(n_elevations=4, n_azimuths=1)
    expected = np.arccos((np.arange(1, 5) - 0.5) / 4.0)
    assert sun.elev == pytest.approx(expected)


def test_sun_azim_uniform_around_circle():
    sun = Sun(n_elevations=1, n_azimuths=8)
    # 8 bins uniformly spaced over [0, 2pi)
    expected = 2.0 * np.pi * np.arange(8) / 8.0
    assert sun.azim == pytest.approx(expected)


def test_sun_validation():
    with pytest.raises(ValueError):
        Sun(n_elevations=0)
    with pytest.raises(ValueError):
        Sun(n_azimuths=-1)
    with pytest.raises(ValueError):
        Sun(size_leaf=0.0)


# ---------------------------------------------------------------------------
# extract_leaves
# ---------------------------------------------------------------------------


def test_extract_single_leaf_at_trunk_tip():
    """A bare trunk has one leaf at its tip (location + length * unit_t)."""
    t = _vertical_trunk(length=2.0)
    L = extract_leaves([t], n_directions=3)
    assert L.n_leaves == 1
    assert L.location.shape == (1, 3)
    assert L.location[0] == pytest.approx((0.0, 0.0, 2.0))
    assert L.branch_index.tolist() == [0]
    assert L.tree_index.tolist() == [0]
    assert L.light_per_direction.shape == (1, 3)
    assert L.light_per_direction.sum() == 0.0  # zero-initialised


def test_extract_only_leaf_branches_not_internals():
    """Internal branches don't contribute leaves; only childless ones do."""
    t = _vertical_trunk()
    t.add_branch_with_geometry(0, length=0.5, diameter=0.05, unit_t=(0, 0, 1), unit_b=(1, 0, 0))
    t.add_branch_with_geometry(0, length=0.5, diameter=0.05, unit_t=(1, 0, 0), unit_b=(0, 1, 0))
    t.reorder()
    L = extract_leaves([t])
    assert L.n_leaves == 2  # trunk now has children; only the two children are leaves
    # Both child bases are at (0, 0, 1); tips differ by unit_t * 0.5.
    tips = sorted(L.location.tolist())
    assert tips[0] == pytest.approx([0.0, 0.0, 1.5])  # straight-up child
    assert tips[1] == pytest.approx([0.5, 0.0, 1.0])  # x-pointing child


def test_extract_multi_tree_tags_tree_index():
    t1 = _vertical_trunk(length=1.0)
    t2 = _vertical_trunk(length=2.0, location=(5.0, 0.0, 0.0))
    L = extract_leaves([t1, t2])
    assert L.n_leaves == 2
    assert L.tree_index.tolist() == [0, 1]
    assert L.location[0] == pytest.approx((0.0, 0.0, 1.0))
    assert L.location[1] == pytest.approx((5.0, 0.0, 2.0))


def test_extract_empty_tree_list_returns_empty_leaves():
    L = extract_leaves([], n_directions=4)
    assert L.n_leaves == 0
    assert L.location.shape == (0, 3)
    assert L.light_per_direction.shape == (0, 4)


# ---------------------------------------------------------------------------
# intercept
# ---------------------------------------------------------------------------


def test_intercept_single_leaf_gets_full_light():
    """One leaf in isolation can't be shaded — light = 1 in every direction."""
    t = _vertical_trunk()
    L = extract_leaves([t])
    sun = Sun(n_elevations=2, n_azimuths=4)
    intercept(L, sun)
    assert L.light_per_direction.shape == (1, sun.n_directions)
    assert np.allclose(L.light_per_direction, 1.0)


def test_intercept_top_leaf_shadows_one_directly_below():
    """Two leaves at the same (X, Y) but different Z under a vertical sun:
    the higher one sees light=1, the lower one is fully shaded."""
    locations = np.array([[0.0, 0.0, 2.0], [0.0, 0.0, 1.0]])
    L = Leaves(
        location=locations,
        branch_index=np.array([0, 1], dtype=np.int32),
        tree_index=np.array([0, 0], dtype=np.int32),
        light_per_direction=np.zeros((2, 0)),
    )
    # Vertical sun: elev=0 ⇒ rotation is identity ⇒ shadow grid is plain (X, Y).
    sun = Sun.from_arrays(elev=[0.0], azim=[0.0])
    intercept(L, sun)
    assert L.light_per_direction[0, 0] == 1.0
    assert L.light_per_direction[1, 0] == 0.0


def test_intercept_two_well_separated_leaves_both_lit():
    """Two leaves far apart (more than size_leaf) never shade each other."""
    locations = np.array([[0.0, 0.0, 1.0], [100.0, 100.0, 1.0]])
    L = Leaves(
        location=locations,
        branch_index=np.array([0, 1], dtype=np.int32),
        tree_index=np.array([0, 0], dtype=np.int32),
        light_per_direction=np.zeros((2, 0)),
    )
    sun = Sun(n_elevations=4, n_azimuths=8, size_leaf=1.0)
    intercept(L, sun)
    assert np.allclose(L.light_per_direction, 1.0)


def test_intercept_reallocates_when_n_dir_mismatches():
    """If the input ``light_per_direction`` has the wrong shape, ``intercept``
    swaps in a correctly-sized buffer (so users don't have to thread
    ``n_directions`` through ``extract_leaves``)."""
    t = _vertical_trunk()
    L = extract_leaves([t], n_directions=0)
    sun = Sun(n_elevations=2, n_azimuths=4)
    intercept(L, sun)
    assert L.light_per_direction.shape == (1, 8)


def test_intercept_zero_leaves_no_crash():
    L = extract_leaves([])
    sun = Sun()
    intercept(L, sun)
    assert L.light_per_direction.shape == (0, sun.n_directions)


# ---------------------------------------------------------------------------
# aggregate_onto_trees
# ---------------------------------------------------------------------------


def test_aggregate_writes_branch_light_back():
    """A standalone leaf gets mean(per_direction) = 1.0 back onto its branch."""
    t = _vertical_trunk()
    L = extract_leaves([t])
    sun = Sun()
    intercept(L, sun)
    # Pre-condition: branch.light defaults to 0.
    assert t.get_light(0) == pytest.approx(0.0)
    aggregate_onto_trees(L, [t])
    assert t.get_light(0) == pytest.approx(1.0)


def test_aggregate_zero_directions_is_a_noop():
    t = _vertical_trunk()
    L = extract_leaves([t], n_directions=0)
    aggregate_onto_trees(L, [t])
    # Nothing to write — branch.light stays at its default.
    assert t.get_light(0) == pytest.approx(0.0)


def test_aggregate_routes_to_right_tree_per_index():
    """Two trees, one leaf each — aggregate writes each leaf's mean light
    onto the correct tree."""
    t1 = _vertical_trunk()
    t2 = _vertical_trunk(location=(100.0, 100.0, 0.0))
    L = extract_leaves([t1, t2])
    sun = Sun()
    intercept(L, sun)
    aggregate_onto_trees(L, [t1, t2])
    assert t1.get_light(0) == pytest.approx(1.0)
    assert t2.get_light(0) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Integration with mechanics + growth
# ---------------------------------------------------------------------------


def test_light_pipeline_drives_secondary_growth():
    """End-to-end: extract leaves -> intercept -> aggregate -> secondary_growth.

    Confirms the light module composes with the Step-9 pipeline without
    any hand-stubbing of branch.light values.
    """
    from mechatree.genome import ConstantSafety
    from mechatree.growth import requested_growth, secondary_growth
    from mechatree.mechanics import calculate_stresses

    t = _vertical_trunk()
    t.add_branch_with_geometry(0, length=0.5, diameter=0.05, unit_t=(0, 0, 1), unit_b=(1, 0, 0))
    t.add_branch_with_geometry(
        0, length=0.5, diameter=0.05, unit_t=(0.3, 0, 0.95), unit_b=(1, 0, 0)
    )
    t.reorder()

    # Light from the real module — not stubs.
    L = extract_leaves([t])
    sun = Sun()
    intercept(L, sun)
    aggregate_onto_trees(L, [t])

    # Some light reached at least one leaf.
    leaf_lights = [t.get_light(i) for i in t.leaf_indices()]
    assert any(light > 0 for light in leaf_lights)

    # Now drive one growth round.
    calculate_stresses(t, leaf_drag_S0=0.5, cauchy=1.0)
    requested_growth(t, ConstantSafety(1.0), maintenance_h=0.005)
    d_before = [t.get_diameter(i) for i in range(t.get_number_of_branches())]
    secondary_growth(t, volume_per_leaf=0.05)
    d_after = [t.get_diameter(i) for i in range(t.get_number_of_branches())]
    # At least one branch grew.
    assert any(da > db for da, db in zip(d_after, d_before, strict=False))


def test_intercept_top_view_decimates_with_n_leaves():
    """With many leaves at the same Z, only the topmost (per cell) gets
    light=1. Average light over many leaves goes down with leaf count.

    This is the basic "self-shading" sanity check: more leaves crowded
    into the same footprint ⇒ less light per leaf on average.
    """
    # 100 leaves all at height z=1.0, scattered in a small (X, Y) box.
    # size_leaf = 1.0, so they all land in roughly the same cell.
    rng = np.random.default_rng(0)
    n = 100
    locations = np.zeros((n, 3))
    locations[:, 0] = rng.uniform(-0.5, 0.5, n)
    locations[:, 1] = rng.uniform(-0.5, 0.5, n)
    locations[:, 2] = 1.0
    L = Leaves(
        location=locations,
        branch_index=np.arange(n, dtype=np.int32),
        tree_index=np.zeros(n, dtype=np.int32),
        light_per_direction=np.zeros((n, 0)),
    )
    # Vertical sun keeps the test intuition simple — under a tilted sun, two
    # leaves at the same (X, Y) but different Z don't project to the same
    # cell, so we wouldn't get self-shading from this geometry.
    sun = Sun.from_arrays(elev=[0.0], azim=[0.0], size_leaf=1.0)
    intercept(L, sun)
    # Most leaves in the same cell as a higher one ⇒ majority get 0.
    lit_fraction = L.light_per_direction.mean()
    assert lit_fraction < 0.5

    # Spreading the same leaves over a wider box (size 100) ⇒ each leaf is
    # in its own cell ⇒ all get light=1.
    locations[:, 0] = rng.uniform(-50.0, 50.0, n)
    locations[:, 1] = rng.uniform(-50.0, 50.0, n)
    L.location = locations
    intercept(L, sun)
    assert L.light_per_direction.mean() == pytest.approx(1.0)


def test_intercept_height_ordering_under_vertical_sun():
    """Three leaves vertically stacked — under a vertical sun the top is
    fully lit, the others are shaded."""
    locations = np.array(
        [
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 2.0],
            [0.0, 0.0, 3.0],
        ]
    )
    L = Leaves(
        location=locations,
        branch_index=np.array([0, 1, 2], dtype=np.int32),
        tree_index=np.zeros(3, dtype=np.int32),
        light_per_direction=np.zeros((3, 0)),
    )
    sun = Sun.from_arrays(elev=[0.0], azim=[0.0], size_leaf=1.0)
    intercept(L, sun)
    assert L.light_per_direction[2, 0] == 1.0  # top
    assert L.light_per_direction[1, 0] == 0.0
    assert L.light_per_direction[0, 0] == 0.0  # bottom


def test_extract_no_leaves_with_empty_n_directions():
    L = extract_leaves([])
    assert L.n_directions == 0
    sun = Sun()
    # No crash even with no leaves and no pre-allocated buffer.
    intercept(L, sun)
    assert L.light_per_direction.shape == (0, sun.n_directions)


# ---------------------------------------------------------------------------
# A Fortran-flavoured edge case: leaves shifted in (X, Y) by less than size_leaf
# should still fall into the same cell when nint() rounds them in.
# ---------------------------------------------------------------------------


def test_intercept_subcell_shift_is_shadowed():
    """Two leaves vertically stacked but offset in X by 0.2 (less than
    half a cell at size_leaf=1.0) still hit the same cell under a vertical
    sun and shade."""
    locations = np.array(
        [
            [0.0, 0.0, 2.0],
            [0.2, 0.0, 1.0],
        ]
    )
    L = Leaves(
        location=locations,
        branch_index=np.array([0, 1], dtype=np.int32),
        tree_index=np.zeros(2, dtype=np.int32),
        light_per_direction=np.zeros((2, 0)),
    )
    sun = Sun.from_arrays(elev=[0.0], azim=[0.0], size_leaf=1.0)
    intercept(L, sun)
    assert L.light_per_direction[0, 0] == 1.0
    assert L.light_per_direction[1, 0] == 0.0


def test_aggregate_skips_internal_branches():
    """A Y-tree with 1 trunk (internal) + 2 leaves: aggregate writes light
    only on the leaves, not on the trunk."""
    t = _vertical_trunk()
    t.add_branch_with_geometry(0, 0.5, 0.05, (0, 0, 1), (1, 0, 0))
    t.add_branch_with_geometry(0, 0.5, 0.05, (0.5, 0, 0.87), (1, 0, 0))
    t.reorder()
    # Pre-set trunk's light to a sentinel — aggregate must not touch it.
    t.set_light(0, math.pi)

    L = extract_leaves([t])
    intercept(L, Sun())
    aggregate_onto_trees(L, [t])

    assert t.get_light(0) == pytest.approx(math.pi)  # trunk untouched
    for leaf_idx in t.leaf_indices():
        # Each leaf got some non-NaN value written.
        assert math.isfinite(t.get_light(leaf_idx))
