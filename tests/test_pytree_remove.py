from mechatree import PyTree


def _trunk_with_two_children() -> PyTree:
    t = PyTree({"length": 1.0})
    t.add_branch(0, {"length": 0.5})
    t.add_branch(0, {"length": 0.5})
    return t


def test_remove_leaf_decrements_count():
    t = _trunk_with_two_children()
    assert t.get_number_of_branches() == 3
    t.remove_branch(1)
    assert t.get_number_of_branches() == 2


def test_remove_leaf_unlinks_from_parent():
    t = _trunk_with_two_children()
    children_before = t.get_children_index(0)
    assert len(children_before) == 2
    t.remove_branch(children_before[0])
    children_after = t.get_children_index(0)
    assert len(children_after) == 1


def test_remove_internal_node_drops_full_subtree():
    t = PyTree({"length": 1.0})
    t.add_branch(0, {"length": 0.5})
    t.add_branch(1, {"length": 0.25})
    t.add_branch(1, {"length": 0.25})
    assert t.get_number_of_branches() == 4
    t.remove_branch(1)
    assert t.get_number_of_branches() == 1


def test_remove_only_child_clears_parents_children_list():
    t = PyTree({"length": 1.0})
    t.add_branch(0, {"length": 0.5})
    t.remove_branch(1)
    assert t.get_number_of_branches() == 1
    assert t.has_parent(0) == 0


def test_remove_then_add_keeps_consistency():
    t = _trunk_with_two_children()
    t.remove_branch(1)
    t.add_branch(0, {"length": 0.7})
    assert t.get_number_of_branches() == 3
    children = t.get_children_index(0)
    assert len(children) == 2
