import math

import pytest

from mechatree import PyTree
from mechatree.geometry import distance_test


@pytest.fixture
def trunk():
    # theta = 1e-6 mirrors the archive: distance_test divides by tan(theta).
    return PyTree({"x": 0.0, "y": 0.0, "theta": 1e-6, "L": 1.0, "grow": 1})


def test_single_trunk_can_grow(trunk):
    assert distance_test(trunk, 0, influence_radius=0.1) == 1


def test_parent_does_not_block_child(trunk):
    # Child is directly attached to trunk's tip; only parent/brother/self are
    # excluded from the check, so two-branch tree => trivially clear.
    trunk.add_branch(0, {"x": 0.0, "y": 1.0, "theta": 1e-6, "L": 1.0, "grow": 1})
    assert distance_test(trunk, 1, influence_radius=0.5) == 1


def test_far_neighbour_does_not_block():
    # Linear chain: trunk -> branch1 -> branch2. Branch 2 is compared against
    # the trunk only (parent=1, brother=none, both excluded). Place branch 2
    # well away from the trunk.
    tree = PyTree({"x": 0.0, "y": 0.0, "theta": 1e-6, "L": 1.0, "grow": 1})
    tree.add_branch(0, {"x": 10.0, "y": 10.0, "theta": 1e-6, "L": 1.0, "grow": 1})
    tree.add_branch(1, {"x": 20.0, "y": 20.0, "theta": 1e-6, "L": 1.0, "grow": 1})
    assert distance_test(tree, 2, influence_radius=0.1) == 1


def test_close_parallel_neighbour_blocks():
    # Linear chain again. Branch 2's tip is placed so it lands on the trunk's
    # line, well inside the influence radius.
    tree = PyTree({"x": 0.0, "y": 0.0, "theta": 1e-6, "L": 1.0, "grow": 1})
    tree.add_branch(0, {"x": 10.0, "y": 10.0, "theta": 1e-6, "L": 1.0, "grow": 1})
    # Branch 2 base at (0, -0.5), theta=1e-6, L=1 -> tip at ~(0, 0.5) which
    # lies on the trunk's body (trunk runs from (0,0) up to ~(0,1)).
    tree.add_branch(1, {"x": 0.0, "y": -0.5, "theta": 1e-6, "L": 1.0, "grow": 1})
    assert distance_test(tree, 2, influence_radius=0.1) == 0


def test_close_neighbour_above_radius_does_not_block():
    # Same chain, but branch 2's tip is offset 0.5 horizontally — outside the
    # 0.1 influence radius around the trunk's line.
    tree = PyTree({"x": 0.0, "y": 0.0, "theta": 1e-6, "L": 1.0, "grow": 1})
    tree.add_branch(0, {"x": 10.0, "y": 10.0, "theta": 1e-6, "L": 1.0, "grow": 1})
    tree.add_branch(1, {"x": 0.5, "y": -0.5, "theta": 1e-6, "L": 1.0, "grow": 1})
    assert distance_test(tree, 2, influence_radius=0.1) == 1


def test_return_type_is_int(trunk):
    # Important: the archive contract is int (0/1), not bool.
    result = distance_test(trunk, 0, 0.1)
    assert isinstance(result, int)
    assert result in (0, 1)


def test_horizontal_branch_does_not_crash():
    # theta=pi/2: 1/tan(pi/2) ~ 1e-16 (not exactly 0). Smoke check that the
    # formula evaluates without raising.
    tree = PyTree({"x": 0.0, "y": 0.0, "theta": 1e-6, "L": 1.0, "grow": 1})
    tree.add_branch(0, {"x": 10.0, "y": 10.0, "theta": math.pi / 2, "L": 1.0, "grow": 1})
    # Just exercise the path.
    distance_test(tree, 1, 0.1)
