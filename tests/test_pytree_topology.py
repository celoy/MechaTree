import pytest

from mechatree import PyTree


def _chain(n_branches: int) -> PyTree:
    """Linear chain of `n_branches` branches: trunk -> c1 -> c2 -> ... -> c_{n-1}."""
    t = PyTree({"length": 1.0})
    for i in range(n_branches - 1):
        t.add_branch(i, {"length": 1.0})
    return t


def _star(n_leaves: int) -> PyTree:
    """trunk with `n_leaves` direct children (no grandchildren)."""
    t = PyTree({"length": 1.0})
    for _ in range(n_leaves):
        t.add_branch(0, {"length": 0.5})
    return t


def test_chain_construction_counts():
    t = _chain(100)
    assert t.get_number_of_branches() == 100


def test_chain_parent_walks_back_to_trunk():
    n = 100
    t = _chain(n)
    idx = n - 1
    steps = 0
    while idx > 0:
        idx = t.get_parent_index(idx)
        steps += 1
    assert idx == 0
    assert steps == n - 1


def test_star_children_returns_full_list():
    n = 20
    t = _star(n)
    children = t.get_children_index(0)
    assert len(children) == n


def test_get_parent_of_trunk_is_minus_one():
    t = _chain(5)
    assert t.get_parent_index(0) == -1


@pytest.mark.slow
def test_large_tree_construction():
    """10 000-branch chain should construct without crashing."""
    t = _chain(10_000)
    assert t.get_number_of_branches() == 10_000


@pytest.mark.slow
def test_large_tree_parent_walk():
    """Walking the parent chain on a large tree is O(N^2) today; this test
    exercises (and pins) the hot path. Phase 4 should make it dramatically
    faster — the test itself just asserts correctness, not timing."""
    n = 2_000
    t = _chain(n)
    idx = n - 1
    while idx > 0:
        idx = t.get_parent_index(idx)
    assert idx == 0
