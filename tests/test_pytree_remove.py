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


def test_remove_with_subtree_when_parent_loses_all_children():
    """Regression for the order-of-operations bug in ``Tree::removeBranch``.

    A non-leaf branch whose removal would empty its parent's children list
    must still drop its entire subtree from ``tree_branches``. The bug:
    ``last_descendant_index`` was being computed AFTER the parent's
    ``removeChild`` call, which (when this branch was the parent's only
    child) collapsed the descendant search to the parent's own (smaller)
    index, making the deletion loop a no-op.
    """
    t = PyTree({"length": 1.0})
    t.add_branch(0, {"length": 0.5})  # child A at index 1 — trunk's only child
    t.add_branch(1, {"length": 0.25})  # grandchild B at index 2 under A
    t.add_branch(1, {"length": 0.25})  # grandchild C at index 2 under A
    # tree_branches = [trunk, A, C, B] in depth-first order. trunk's only
    # child is A; A has children [B, C].
    assert t.get_number_of_branches() == 4

    # Remove A — this empties trunk's children list AND must take B, C with it.
    t.remove_branch(1)

    assert t.get_number_of_branches() == 1
    assert t.get_children_index(0) == []
    # Re-adding works (no stale state in branch_to_index):
    t.add_branch(0, {"length": 0.3})
    assert t.get_number_of_branches() == 2
    assert t.get_parent_index(1) == 0
