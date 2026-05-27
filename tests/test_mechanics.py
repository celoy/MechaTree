"""Mechanics tests — wind force on a branch and 4-angle stress propagation.

These cases lean on hand-derivable physics on small trees rather than
Fortran reference numbers (which arrive in Step 11).
"""

import math

import pytest

from mechatree import PyTree
from mechatree.mechanics import calculate_stresses, wind_force


def _trunk(length=1.0, diameter=0.1, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)):
    t = PyTree({})
    t.set_length(0, length)
    t.set_diameter(0, diameter)
    t.set_unit_t(0, unit_t)
    t.set_unit_b(0, unit_b)
    return t


def test_wind_force_perpendicular_unit_branch():
    """Wind (1,0,0) on a trunk pointing up: force = (½·V²·d·L, 0, 0).

    The ½ is the ½ρU² dynamic-pressure factor restored to ``wind_force``
    (the Fortran omitted it and folded it into Cauchy); with d=0.1, L=1,
    V=1 that gives 0.05.
    """
    t = _trunk()
    force, moment = wind_force(t, 0, (1.0, 0.0, 0.0))
    assert force == pytest.approx((0.05, 0.0, 0.0))
    # moment = (0,0,L/2) × force = (0, L·½·V²·d/2, 0)
    assert moment == pytest.approx((0.0, 0.025, 0.0))


def test_wind_force_wind_magnitude_squared():
    """Doubling wind speed quadruples force magnitude."""
    t = _trunk()
    f1, _ = wind_force(t, 0, (1.0, 0.0, 0.0))
    f2, _ = wind_force(t, 0, (2.0, 0.0, 0.0))
    assert math.hypot(*f2) == pytest.approx(4.0 * math.hypot(*f1))


def test_wind_force_parallel_is_zero():
    """Wind aligned with the branch ⇒ no projected area ⇒ zero force/moment.

    Guards against div-by-zero on the Nn normalisation in `wind_force`.
    """
    t = _trunk()
    force, moment = wind_force(t, 0, (0.0, 0.0, 5.0))
    assert force == pytest.approx((0.0, 0.0, 0.0))
    assert moment == pytest.approx((0.0, 0.0, 0.0))


def test_wind_force_zero_wind_is_zero():
    t = _trunk()
    force, moment = wind_force(t, 0, (0.0, 0.0, 0.0))
    assert force == pytest.approx((0.0, 0.0, 0.0))
    assert moment == pytest.approx((0.0, 0.0, 0.0))


def test_calculate_stresses_single_branch_is_positive():
    """A single vertical branch under the four-angle sweep gets a non-zero
    max_stress, regardless of which of the four winds is "worst"."""
    t = _trunk(length=1.0, diameter=0.1)
    t.reorder()
    calculate_stresses(t, leaf_drag_S0=1.0, cauchy=1.0)
    assert t.get_max_stress(0) > 0.0
    # Stress on a thicker branch should be lower (d³ in the denominator).
    t2 = _trunk(length=1.0, diameter=0.2)
    t2.reorder()
    calculate_stresses(t2, leaf_drag_S0=1.0, cauchy=1.0)
    assert t2.get_max_stress(0) < t.get_max_stress(0)


def test_calculate_stresses_zero_drag_zero_stiffness_is_zero():
    """No leaf drag, no material stiffness ⇒ no stress."""
    t = _trunk()
    t.reorder()
    calculate_stresses(t, leaf_drag_S0=0.0, cauchy=0.0)
    assert t.get_max_stress(0) == pytest.approx(0.0)


def test_calculate_stresses_trunk_accumulates_children():
    """Y-tree: the trunk's force vector under any one angle is the sum of its
    children's plus its own wind drag. After 4 angles, trunk max_stress > 0
    and the trunk's force magnitude exceeds each child's."""
    t = _trunk(length=1.0, diameter=0.2)
    # Two children: same length and diameter, one slightly tilted in x, one in y.
    t.add_branch_with_geometry(
        0,
        length=0.5,
        diameter=0.1,
        unit_t=(0.3, 0.0, 0.95),
        unit_b=(1.0, 0.0, 0.0),
    )
    t.add_branch_with_geometry(
        0,
        length=0.5,
        diameter=0.1,
        unit_t=(0.0, 0.3, 0.95),
        unit_b=(1.0, 0.0, 0.0),
    )
    t.reorder()
    calculate_stresses(t, leaf_drag_S0=1.0, cauchy=1.0)

    trunk_force = t.get_force(0)
    child1_force = t.get_force(1)
    child2_force = t.get_force(2)

    trunk_mag = math.hypot(*trunk_force)
    c1_mag = math.hypot(*child1_force)
    c2_mag = math.hypot(*child2_force)
    # Trunk carries both children's forces plus its own — strictly greater than
    # either child individually.
    assert trunk_mag > c1_mag
    assert trunk_mag > c2_mag
    assert t.get_max_stress(0) > 0.0
