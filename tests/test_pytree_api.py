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


def test_get_leaf_tips_batch_matches_per_leaf_path():
    """Batched accessor must match the per-leaf get_location + length * unit_t.

    Used by ``mechatree.light.extract_leaves`` to skip the per-leaf
    Python loop. Validate on a grown tree where leaves carry non-trivial
    geometry."""
    import numpy as np

    from mechatree.config import Config
    from mechatree.simulate import grow_tree

    tree = grow_tree(Config(), n_generations=15, seed=0)
    tips, branch_idx = tree.get_leaf_tips_batch()
    assert tips.dtype == np.float64
    assert branch_idx.dtype == np.int32
    assert tips.shape == (branch_idx.size, 3)
    assert list(branch_idx) == tree.leaf_indices()

    # Reconstruct via per-leaf path and compare.
    ref = np.empty_like(tips)
    for i, bi in enumerate(branch_idx):
        base = np.asarray(tree.get_location(int(bi)))
        ut = np.asarray(tree.get_unit_t(int(bi)))
        ref[i] = base + tree.get_length(int(bi)) * ut
    np.testing.assert_allclose(tips, ref, atol=1e-12)


def test_get_leaf_tips_batch_empty_tree():
    tree = PyTree({})
    tips, branch_idx = tree.get_leaf_tips_batch()
    assert tips.shape == (1, 3)  # the trunk alone is a leaf
    assert branch_idx.shape == (1,)


def test_set_lights_batch_round_trip():
    """Writing via the batched setter must show up in subsequent get_light."""
    import numpy as np

    from mechatree.config import Config
    from mechatree.simulate import grow_tree

    tree = grow_tree(Config(), n_generations=10, seed=0)
    _tips, branch_idx = tree.get_leaf_tips_batch()
    values = np.linspace(0.1, 0.9, branch_idx.size)
    tree.set_lights_batch(branch_idx, values)
    read_back = np.array([tree.get_light(int(b)) for b in branch_idx])
    np.testing.assert_allclose(read_back, values, atol=1e-12)


def test_set_lights_batch_length_mismatch_raises():
    import numpy as np

    tree = PyTree({})
    with pytest.raises(ValueError):
        tree.set_lights_batch(np.array([0], dtype=np.int32), np.array([0.0, 1.0]))
