"""Pruning tests — wind-driven stochastic branch removal."""

from mechatree import PyTree
from mechatree.pruning import prune


def _vertical_trunk(length=1.0, diameter=0.1):
    t = PyTree({})
    t.set_length(0, length)
    t.set_diameter(0, diameter)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    t.set_location(0, (0.0, 0.0, 0.0))
    return t


def _build_y_tree():
    """Trunk plus two leaves, all geometry set."""
    t = _vertical_trunk(length=1.0, diameter=0.2)
    t.add_branch_with_geometry(
        0,
        length=0.5,
        diameter=0.05,
        unit_t=(0.3, 0.0, 0.95),
        unit_b=(1.0, 0.0, 0.0),
    )
    t.add_branch_with_geometry(
        0,
        length=0.5,
        diameter=0.05,
        unit_t=(0.0, 0.3, 0.95),
        unit_b=(1.0, 0.0, 0.0),
    )
    t.reorder()
    return t


def test_prune_zero_wind_keeps_all():
    """Zero wind ⇒ zero stress ⇒ P_fail = 0 ⇒ no cuts."""
    t = _build_y_tree()
    t.set_seed(7)
    n_before = t.get_number_of_branches()
    cut = prune(t, wind=(0.0, 0.0, 0.0), leaf_drag_S0=1.0, cauchy=1.0)
    assert cut == 0
    assert t.get_number_of_branches() == n_before


def test_prune_thin_branches_cut_under_strong_wind():
    """Hurricane wind on thin children — they go; trunk stays."""
    t = _build_y_tree()
    t.set_seed(7)
    n_before = t.get_number_of_branches()
    cut = prune(t, wind=(100.0, 0.0, 0.0), leaf_drag_S0=1.0, cauchy=1.0)
    assert cut > 0
    assert t.get_number_of_branches() == n_before - cut
    # Trunk always survives.
    assert t.get_number_of_branches() >= 1


def test_prune_trunk_never_cut_single_branch_tree():
    """A one-branch tree can never lose anything, even under huge winds."""
    t = _vertical_trunk()
    t.set_seed(0)
    cut = prune(t, wind=(100.0, 0.0, 0.0), leaf_drag_S0=1.0, cauchy=1.0)
    assert cut == 0
    assert t.get_number_of_branches() == 1


def test_prune_is_reproducible_with_seed():
    """Same seed + same inputs ⇒ same number of branches cut."""

    def run(seed):
        t = _build_y_tree()
        t.set_seed(seed)
        return prune(t, wind=(50.0, 0.0, 0.0), leaf_drag_S0=1.0, cauchy=1.0)

    assert run(seed=42) == run(seed=42)


def test_prune_removes_subtree():
    """Cutting an internal branch removes its entire subtree."""
    t = _vertical_trunk(length=1.0, diameter=0.2)
    # internal child (a bit stronger), then a grandchild (also thin)
    t.add_branch_with_geometry(
        0,
        length=0.5,
        diameter=0.05,
        unit_t=(0.3, 0.0, 0.95),
        unit_b=(1.0, 0.0, 0.0),
    )
    t.add_branch_with_geometry(
        1,
        length=0.3,
        diameter=0.02,
        unit_t=(0.5, 0.0, 0.85),
        unit_b=(1.0, 0.0, 0.0),
    )
    t.reorder()
    t.set_seed(1)
    n_before = t.get_number_of_branches()
    cut = prune(t, wind=(100.0, 0.0, 0.0), leaf_drag_S0=1.0, cauchy=1.0)
    assert cut > 0
    # If branch 1 was cut, branch 2 (the grandchild) is also gone.
    assert t.get_number_of_branches() == n_before - cut
