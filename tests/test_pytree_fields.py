"""Tests for the typed mechanics fields added to Branch in Step 9 (PR1).

These fields (length, diameter, light, stress, ..., location, unit_t, ...) are
the hot path for the mechanics, growth and pruning routines that land in PR2.
They live alongside — but independent of — the existing `properties` map.
"""

import pytest

from mechatree import PyTree

SCALAR_FIELDS = [
    ("length", "get_length", "set_length", 1.5),
    ("diameter", "get_diameter", "set_diameter", 0.25),
    ("light", "get_light", "set_light", 0.8),
    ("stress", "get_stress", "set_stress", 12.7),
    ("max_stress", "get_max_stress", "set_max_stress", 99.0),
    ("vol_growth", "get_vol_growth", "set_vol_growth", 0.03),
    ("vol_summed", "get_vol_summed", "set_vol_summed", 4.5),
    ("maintenance_vol", "get_maintenance_vol", "set_maintenance_vol", 0.01),
]


VECTOR_FIELDS = [
    ("location", "get_location", "set_location", (1.0, 2.0, 3.0)),
    ("unit_t", "get_unit_t", "set_unit_t", (0.0, 1.0, 0.0)),
    ("unit_b", "get_unit_b", "set_unit_b", (0.0, 0.0, 1.0)),
    ("force", "get_force", "set_force", (-1.0, 0.5, 2.5)),
    ("moment", "get_moment", "set_moment", (3.0, -4.0, 5.0)),
]


@pytest.fixture
def trunk():
    return PyTree({"length": 1.0, "radius": 0.1})


def test_typed_scalar_defaults_are_zero(trunk):
    """Fresh branches have all typed scalars at zero, independent of any
    dict the user passed at construction. ``length`` is the one exception —
    it defaults to 1.0 (the simulator-wide ``twig_length``)."""
    assert trunk.get_length(0) == pytest.approx(1.0)
    assert trunk.get_diameter(0) == pytest.approx(0.0)
    assert trunk.get_light(0) == pytest.approx(0.0)
    assert trunk.get_stress(0) == pytest.approx(0.0)
    assert trunk.get_max_stress(0) == pytest.approx(0.0)
    assert trunk.get_vol_growth(0) == pytest.approx(0.0)
    assert trunk.get_vol_summed(0) == pytest.approx(0.0)
    assert trunk.get_maintenance_vol(0) == pytest.approx(0.0)
    assert trunk.get_nb_leaves(0) == 0


def test_typed_vector_defaults(trunk):
    """`location`, `force`, `moment` default to (0,0,0); `unit_t` is the
    upward unit vector and `unit_b` is along x — sensible defaults for the
    trunk before the orchestrator overrides them."""
    assert trunk.get_location(0) == pytest.approx((0.0, 0.0, 0.0))
    assert trunk.get_unit_t(0) == pytest.approx((0.0, 0.0, 1.0))
    assert trunk.get_unit_b(0) == pytest.approx((1.0, 0.0, 0.0))
    assert trunk.get_force(0) == pytest.approx((0.0, 0.0, 0.0))
    assert trunk.get_moment(0) == pytest.approx((0.0, 0.0, 0.0))


@pytest.mark.parametrize("name, getter, setter, value", SCALAR_FIELDS)
def test_scalar_field_roundtrip(trunk, name, getter, setter, value):
    getattr(trunk, setter)(0, value)
    assert getattr(trunk, getter)(0) == pytest.approx(value)


def test_nb_leaves_roundtrip(trunk):
    trunk.set_nb_leaves(0, 7)
    assert trunk.get_nb_leaves(0) == 7


@pytest.mark.parametrize("name, getter, setter, value", VECTOR_FIELDS)
def test_vector_field_roundtrip(trunk, name, getter, setter, value):
    getattr(trunk, setter)(0, value)
    assert getattr(trunk, getter)(0) == pytest.approx(value)


def test_add_branch_with_geometry_derives_child_location(trunk):
    """Trunk at origin pointing up, length 1 ⇒ child sits at (0, 0, 1)."""
    trunk.set_length(0, 1.0)
    trunk.set_unit_t(0, (0.0, 0.0, 1.0))
    trunk.set_unit_b(0, (1.0, 0.0, 0.0))
    idx = trunk.add_branch_with_geometry(
        0,
        length=0.5,
        diameter=0.05,
        unit_t=(0.0, 0.0, 1.0),
        unit_b=(1.0, 0.0, 0.0),
    )
    assert idx == 1
    assert trunk.get_location(idx) == pytest.approx((0.0, 0.0, 1.0))
    assert trunk.get_length(idx) == pytest.approx(0.5)
    assert trunk.get_diameter(idx) == pytest.approx(0.05)
    assert trunk.get_unit_t(idx) == pytest.approx((0.0, 0.0, 1.0))
    assert trunk.get_parent_index(idx) == 0


def test_add_branch_with_geometry_oblique_parent():
    """Parent pointing along x with length 2 ⇒ child sits at (2, 0, 0)."""
    t = PyTree({"length": 1.0})
    t.set_length(0, 2.0)
    t.set_unit_t(0, (1.0, 0.0, 0.0))
    idx = t.add_branch_with_geometry(
        0,
        length=0.3,
        diameter=0.02,
        unit_t=(1.0, 0.0, 0.0),
        unit_b=(0.0, 1.0, 0.0),
    )
    assert t.get_location(idx) == pytest.approx((2.0, 0.0, 0.0))


def test_property_map_and_typed_field_are_independent(trunk):
    """A trunk built with `{"length": 7.0}` exposes that value via the
    property map, but the typed `length` field still holds its default
    (1.0). The two storages do not auto-mirror."""
    t = PyTree({"length": 7.0})
    assert t.get_property(0, "length") == pytest.approx(7.0)
    assert t.get_length(0) == pytest.approx(1.0)

    # Writing through the typed setter does not touch the property map.
    t.set_length(0, 3.0)
    assert t.get_length(0) == pytest.approx(3.0)
    assert t.get_property(0, "length") == pytest.approx(7.0)


def test_reserve_roundtrip(trunk):
    assert trunk.get_reserve() == pytest.approx(0.0)
    trunk.set_reserve(2.5)
    assert trunk.get_reserve() == pytest.approx(2.5)


def test_setters_out_of_range_raise(trunk):
    with pytest.raises(IndexError):
        trunk.set_length(999, 1.0)
    with pytest.raises(IndexError):
        trunk.set_location(999, (1.0, 2.0, 3.0))
    with pytest.raises(IndexError):
        trunk.get_diameter(-1)


def test_vector_setter_wrong_length_raises(trunk):
    with pytest.raises((TypeError, ValueError)):
        trunk.set_location(0, (1.0, 2.0))
    with pytest.raises((TypeError, ValueError)):
        trunk.set_unit_t(0, (1.0, 2.0, 3.0, 4.0))


def test_add_branch_with_geometry_invalid_parent_raises(trunk):
    with pytest.raises(IndexError):
        trunk.add_branch_with_geometry(
            999,
            length=0.5,
            diameter=0.05,
            unit_t=(0.0, 0.0, 1.0),
            unit_b=(1.0, 0.0, 0.0),
        )


def test_dict_built_branch_preserves_property_map(trunk):
    """Existing dict-driven `add_branch` keeps populating the property map
    and does not touch the typed fields."""
    trunk.add_branch(0, {"length": 0.5, "radius": 0.05})
    assert trunk.get_property(1, "length") == pytest.approx(0.5)
    assert trunk.get_property(1, "radius") == pytest.approx(0.05)
    # Typed `length` default is 1.0; the property-map "length" is independent.
    assert trunk.get_length(1) == pytest.approx(1.0)
    assert trunk.get_diameter(1) == pytest.approx(0.0)
