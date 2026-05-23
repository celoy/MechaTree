import pytest

from mechatree import PyTree


@pytest.fixture
def trunk():
    return PyTree({"length": 1.0, "radius": 0.1})


def test_add_branch_increases_count(trunk):
    trunk.add_branch(0, {"length": 0.5, "radius": 0.05})
    assert trunk.get_number_of_branches() == 2


def test_children_index_returns_list(trunk):
    trunk.add_branch(0, {"length": 0.5, "radius": 0.05})
    children = trunk.get_children_index(0)
    assert list(children) == [1]


def test_parent_index_of_trunk_is_minus_one(trunk):
    assert trunk.get_parent_index(0) == -1


def test_property_roundtrip(trunk):
    assert trunk.get_property(0, "length") == pytest.approx(1.0)
    trunk.set_property(0, "length", 2.5)
    assert trunk.get_property(0, "length") == pytest.approx(2.5)


def test_strahler_set_and_get(trunk):
    trunk.add_branch(0, {"length": 0.5, "radius": 0.05})
    trunk.add_branch(0, {"length": 0.5, "radius": 0.05})
    trunk.set_strahler()
    assert trunk.get_strahler(0) >= 1
