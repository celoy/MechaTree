import pytest

from mechatree import PyTree


def _balanced_binary_tree() -> PyTree:
    """trunk -> [A, B], each with 2 leaves; 7 branches total.

    Tree::addBranch inserts each new branch at parent_index+1, which shifts
    later branches' indices. To get predictable final indices we build the
    right subtree first, then the left, depth-first.

    Final layout: trunk=0, A=1, AA=2, AB=3, B=4, BA=5, BB=6.
    """
    t = PyTree({"length": 1.0})
    # Right subtree first (ends up at idx 4..6)
    t.add_branch(0, {"length": 0.5})  # B
    t.add_branch(1, {"length": 0.25})  # BB
    t.add_branch(1, {"length": 0.25})  # BA (BB shifts to 3)
    # Then left subtree (B and its subtree shift right by 3 idx)
    t.add_branch(0, {"length": 0.5})  # A
    t.add_branch(1, {"length": 0.25})  # AB
    t.add_branch(1, {"length": 0.25})  # AA (AB shifts to 3)
    return t


def test_strahler_leaves_are_order_one():
    t = _balanced_binary_tree()
    t.set_strahler()
    for leaf in (2, 3, 5, 6):
        assert t.get_strahler(leaf) == 1


def test_strahler_trunk_of_balanced_binary_is_three():
    t = _balanced_binary_tree()
    t.set_strahler()
    assert t.get_strahler(0) == 3


def test_strahler_distribution_balanced_binary():
    t = _balanced_binary_tree()
    t.set_strahler()
    dist = t.get_strahler_distribution()
    # 4 leaves (order 1), 2 internal-of-two-leaves (order 2), 1 trunk (order 3)
    assert dist == {1: 4, 2: 2, 3: 1}


def test_strahler_set_is_idempotent():
    t = _balanced_binary_tree()
    t.set_strahler()
    dist_once = t.get_strahler_distribution()
    t.set_strahler()
    dist_twice = t.get_strahler_distribution()
    assert dist_once == dist_twice


def test_horton_trunk_matches_strahler():
    t = _balanced_binary_tree()
    t.set_horton()
    assert t.get_horton(0) == t.get_strahler(0)


def test_horton_distribution_balanced_binary():
    t = _balanced_binary_tree()
    t.set_horton()
    dist = t.get_horton_distribution()
    # Balanced binary: one long path inherits its parent's Horton order,
    # branching off at each level produces one new Horton-1 branch per leaf.
    assert sum(dist.values()) >= 1
    assert max(dist) == t.get_horton(0)


def test_horton_set_is_idempotent():
    t = _balanced_binary_tree()
    t.set_horton()
    dist_once = t.get_horton_distribution()
    t.set_horton()
    dist_twice = t.get_horton_distribution()
    assert dist_once == dist_twice


def test_mean_agg_prop_s_roundtrip():
    t = _balanced_binary_tree()
    t.set_strahler()
    means = t.mean_agg_prop_s("length")
    # 4 leaves of length 0.25 -> mean 0.25 at order 1
    assert means[1] == pytest.approx(0.25)
    # 2 mid branches of length 0.5 -> mean 0.5 at order 2
    assert means[2] == pytest.approx(0.5)
    # 1 trunk of length 1.0 -> mean 1.0 at order 3
    assert means[3] == pytest.approx(1.0)


def test_mean_agg_prop_h_roundtrip():
    t = _balanced_binary_tree()
    t.set_horton()
    means = t.mean_agg_prop_h("length")
    assert set(means.keys()) == set(t.get_horton_distribution().keys())
    for value in means.values():
        assert value > 0
