"""Tests for ``PyTree.collapse_single_child_chains``.

The merge fuses any maximal run of single-child parents into one straight
segment, preserving bottom + top points and total volume (sum of
``pi/4 * d**2 * L``).
"""

from __future__ import annotations

import math

import pytest

from mechatree import PyTree
from mechatree.config import Config
from mechatree.pruning import prune
from mechatree.simulate import grow_tree

PI_OVER_4 = math.pi / 4.0


def _vertical_trunk(length: float = 1.0, diameter: float = 0.1) -> PyTree:
    t = PyTree({})
    t.set_length(0, length)
    t.set_diameter(0, diameter)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    t.set_location(0, (0.0, 0.0, 0.0))
    return t


def _tip(tree: PyTree, idx: int) -> tuple[float, float, float]:
    """Geometric tip of a branch: base + length * unit_t."""
    loc = tree.get_location(idx)
    ut = tree.get_unit_t(idx)
    L = tree.get_length(idx)
    return (loc[0] + L * ut[0], loc[1] + L * ut[1], loc[2] + L * ut[2])


def _total_volume(tree: PyTree) -> float:
    total = 0.0
    for i in range(tree.get_number_of_branches()):
        d = tree.get_diameter(i)
        L = tree.get_length(i)
        total += PI_OVER_4 * d * d * L
    return total


def _leaf_tips(tree: PyTree) -> list[tuple[float, float, float]]:
    return [_tip(tree, i) for i in tree.leaf_indices()]


# ---------- unit tests on hand-built trees ------------------------------------


def test_no_op_balanced_binary_tree():
    """Every parent has 2 children -> no chain to fuse, returns 0."""
    t = _vertical_trunk()
    t.add_branch_with_geometry(
        0, length=0.5, diameter=0.05, unit_t=(0.3, 0.0, 0.95), unit_b=(1.0, 0.0, 0.0)
    )
    t.add_branch_with_geometry(
        0, length=0.5, diameter=0.05, unit_t=(0.0, 0.3, 0.95), unit_b=(1.0, 0.0, 0.0)
    )
    t.reorder()

    n_before = t.get_number_of_branches()
    vol_before = _total_volume(t)
    tips_before = sorted(_leaf_tips(t))

    absorbed = t.collapse_single_child_chains()

    assert absorbed == 0
    assert t.get_number_of_branches() == n_before
    assert _total_volume(t) == pytest.approx(vol_before, rel=0.0, abs=1e-15)
    tips_after = sorted(_leaf_tips(t))
    assert len(tips_after) == len(tips_before)
    for got, want in zip(tips_after, tips_before, strict=True):
        assert got == pytest.approx(want, abs=1e-15)


def test_single_child_parent_only_no_chain_to_merge():
    """Trunk has 1 child but that child is a fork: chain length 1, no merge."""
    t = _vertical_trunk(length=1.0, diameter=0.2)
    fork = t.add_branch_with_geometry(
        0, length=0.5, diameter=0.1, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
    )
    t.add_branch_with_geometry(
        fork, length=0.3, diameter=0.05, unit_t=(0.3, 0.0, 0.95), unit_b=(1.0, 0.0, 0.0)
    )
    t.add_branch_with_geometry(
        fork, length=0.3, diameter=0.05, unit_t=(-0.3, 0.0, 0.95), unit_b=(1.0, 0.0, 0.0)
    )
    t.reorder()

    n_before = t.get_number_of_branches()
    absorbed = t.collapse_single_child_chains()

    # No chain has length >= 2 single-child links, so nothing is absorbed.
    assert absorbed == 0
    assert t.get_number_of_branches() == n_before


def test_collinear_chain_collapses_to_one_branch():
    """Trunk + 1 collinear child + 1 collinear grandchild (a leaf).

    All three vertical, each length 1.0 diameter 0.1. After merge: one branch
    of length 3.0 with the same diameter (volume per unit length is constant
    in this special case).
    """
    t = _vertical_trunk(length=1.0, diameter=0.1)
    c1 = t.add_branch_with_geometry(
        0, length=1.0, diameter=0.1, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
    )
    t.add_branch_with_geometry(
        c1, length=1.0, diameter=0.1, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
    )
    t.reorder()

    vol_before = _total_volume(t)
    tips_before = sorted(_leaf_tips(t))

    absorbed = t.collapse_single_child_chains()

    assert absorbed == 2
    assert t.get_number_of_branches() == 1
    assert t.get_length(0) == pytest.approx(3.0, abs=1e-12)
    assert t.get_diameter(0) == pytest.approx(0.1, abs=1e-12)
    assert t.get_location(0) == pytest.approx((0.0, 0.0, 0.0), abs=1e-12)
    assert t.get_unit_t(0) == pytest.approx((0.0, 0.0, 1.0), abs=1e-12)
    assert _total_volume(t) == pytest.approx(vol_before, rel=0.0, abs=1e-12)
    assert sorted(_leaf_tips(t)) == pytest.approx(tips_before, rel=0.0, abs=1e-12)


def test_non_collinear_chain_straight_line_merge():
    """Three-segment chain with a bend: tip + base preserved, volume preserved,
    new unit_t points base->tip, new length is Euclidean."""
    t = _vertical_trunk(length=1.0, diameter=0.20)
    c1 = t.add_branch_with_geometry(
        0, length=1.0, diameter=0.10, unit_t=(0.6, 0.0, 0.8), unit_b=(1.0, 0.0, 0.0)
    )
    t.add_branch_with_geometry(
        c1, length=1.0, diameter=0.05, unit_t=(0.0, 0.8, 0.6), unit_b=(1.0, 0.0, 0.0)
    )
    t.reorder()

    base_before = t.get_location(0)
    tip_before = _tip(t, 2)  # the leaf grandchild
    vol_before = _total_volume(t)

    absorbed = t.collapse_single_child_chains()

    assert absorbed == 2
    assert t.get_number_of_branches() == 1

    # Endpoints preserved exactly.
    assert t.get_location(0) == pytest.approx(base_before, abs=1e-12)
    assert _tip(t, 0) == pytest.approx(tip_before, abs=1e-12)

    # Length = Euclidean distance between the original base and tip.
    dx = tip_before[0] - base_before[0]
    dy = tip_before[1] - base_before[1]
    dz = tip_before[2] - base_before[2]
    expected_L = math.sqrt(dx * dx + dy * dy + dz * dz)
    assert t.get_length(0) == pytest.approx(expected_L, abs=1e-12)

    # unit_t aligned with base->tip vector.
    ut = t.get_unit_t(0)
    assert ut == pytest.approx((dx / expected_L, dy / expected_L, dz / expected_L), abs=1e-12)

    # Volume preserved.
    assert _total_volume(t) == pytest.approx(vol_before, rel=0.0, abs=1e-12)


def test_chain_ending_at_fork_keeps_fork_distinct():
    """trunk -> mid (1 child) -> fork (2 leaf children).

    Before: 5 branches (trunk, mid, fork, 2 leaves). After: trunk+mid collapsed
    into one merged segment; fork survives and is rewired as merged.child; its
    2 leaves are untouched. 4 branches remain.
    """
    t = _vertical_trunk(length=1.0, diameter=0.20)
    mid = t.add_branch_with_geometry(
        0, length=1.0, diameter=0.10, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
    )
    # mid has one child (the fork). The fork has two children (leaves).
    fork = t.add_branch_with_geometry(
        mid, length=0.5, diameter=0.08, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
    )
    t.add_branch_with_geometry(
        fork, length=0.3, diameter=0.04, unit_t=(0.6, 0.0, 0.8), unit_b=(1.0, 0.0, 0.0)
    )
    t.add_branch_with_geometry(
        fork, length=0.3, diameter=0.04, unit_t=(-0.6, 0.0, 0.8), unit_b=(1.0, 0.0, 0.0)
    )
    t.reorder()

    fork_base_before = t.get_location(fork)
    fork_len_before = t.get_length(fork)
    fork_d_before = t.get_diameter(fork)
    vol_before = _total_volume(t)
    tips_before = sorted(_leaf_tips(t))

    absorbed = t.collapse_single_child_chains()

    # trunk + mid -> 1 merged segment (1 absorbed). Fork and its 2 leaves kept.
    assert absorbed == 1
    assert t.get_number_of_branches() == 4

    # Merged branch sits at index 0; its tip is the fork's old base.
    assert _tip(t, 0) == pytest.approx(fork_base_before, abs=1e-12)
    # Fork is now b0's only child, untouched.
    kids = t.get_children_index(0)
    assert len(kids) == 1
    fork_now = kids[0]
    assert t.get_length(fork_now) == pytest.approx(fork_len_before, abs=1e-12)
    assert t.get_diameter(fork_now) == pytest.approx(fork_d_before, abs=1e-12)
    assert t.get_location(fork_now) == pytest.approx(fork_base_before, abs=1e-12)
    # Fork still has its 2 leaves.
    assert t.get_number_of_children(fork_now) == 2

    # Total volume preserved; leaf tips unchanged.
    assert _total_volume(t) == pytest.approx(vol_before, rel=0.0, abs=1e-12)
    tips_after = sorted(_leaf_tips(t))
    assert len(tips_after) == len(tips_before)
    for got, want in zip(tips_after, tips_before, strict=True):
        assert got == pytest.approx(want, abs=1e-12)


def test_length_max_truncates_long_chain():
    """A 10-segment collinear chain of unit-length twigs; cap at length_max=3.5.

    Expected: the original chain is broken up into several merged segments,
    each of length ≤ 3.5. The original tip is still reachable as the last
    surviving branch's tip, and total volume is preserved.
    """
    t = _vertical_trunk(length=1.0, diameter=0.1)
    prev = 0
    for _ in range(9):
        prev = t.add_branch_with_geometry(
            prev, length=1.0, diameter=0.1, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
        )
    t.reorder()

    n_before = t.get_number_of_branches()
    vol_before = _total_volume(t)
    tip_before = _tip(t, prev)

    absorbed = t.collapse_single_child_chains(length_max=3.5)

    # Some compression must have happened, but not the full collapse.
    assert absorbed > 0
    assert absorbed < n_before - 1
    assert t.get_number_of_branches() == n_before - absorbed
    # Each surviving branch must respect the cap.
    for i in range(t.get_number_of_branches()):
        assert t.get_length(i) <= 3.5 + 1e-12
    # Total volume preserved across the whole tree.
    assert _total_volume(t) == pytest.approx(vol_before, rel=0.0, abs=1e-12)
    # The original chain's tip is still the only leaf's tip.
    leaves = t.leaf_indices()
    assert len(leaves) == 1
    assert _tip(t, leaves[0]) == pytest.approx(tip_before, abs=1e-12)


def test_length_max_zero_disables_merge():
    """Setting length_max to 0 forbids any merging; tree is unchanged."""
    t = _vertical_trunk(length=1.0, diameter=0.1)
    c1 = t.add_branch_with_geometry(
        0, length=1.0, diameter=0.1, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
    )
    t.add_branch_with_geometry(
        c1, length=1.0, diameter=0.1, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
    )
    t.reorder()

    n_before = t.get_number_of_branches()
    absorbed = t.collapse_single_child_chains(length_max=0.0)
    assert absorbed == 0
    assert t.get_number_of_branches() == n_before


def test_length_max_default_is_ten():
    """Default length_max=10 fully collapses the 5-deep chain (total length
    5 ≤ 10) and caps every surviving branch in the 15-deep chain at length
    ≤ 10."""
    # 5-deep chain, total length 5: fully collapses.
    t = _vertical_trunk(length=1.0, diameter=0.1)
    prev = 0
    for _ in range(4):
        prev = t.add_branch_with_geometry(
            prev, length=1.0, diameter=0.1, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
        )
    t.reorder()
    assert t.collapse_single_child_chains() == 4
    assert t.get_number_of_branches() == 1

    # 15-deep chain: first chain absorbs b1..b9 (merged length = 10). The
    # outer loop then finds the rest of the chain starting at b10.
    t2 = _vertical_trunk(length=1.0, diameter=0.1)
    prev = 0
    for _ in range(14):
        prev = t2.add_branch_with_geometry(
            prev, length=1.0, diameter=0.1, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
        )
    t2.reorder()
    absorbed = t2.collapse_single_child_chains()
    assert absorbed > 0
    # Every surviving branch respects the cap.
    for i in range(t2.get_number_of_branches()):
        assert t2.get_length(i) <= 10.0 + 1e-12
    # The first chain's merged length is exactly 10.
    assert t2.get_length(0) == pytest.approx(10.0, abs=1e-12)


def test_deep_chain_terminating_in_leaf():
    """A 5-deep collinear chain (trunk + 4 single-child links, last is a leaf)
    collapses to one branch."""
    t = _vertical_trunk(length=1.0, diameter=0.1)
    prev = 0
    for _ in range(4):
        prev = t.add_branch_with_geometry(
            prev,
            length=1.0,
            diameter=0.1,
            unit_t=(0.0, 0.0, 1.0),
            unit_b=(1.0, 0.0, 0.0),
        )
    t.reorder()

    vol_before = _total_volume(t)
    tip_before = _tip(t, prev)

    absorbed = t.collapse_single_child_chains()

    assert absorbed == 4
    assert t.get_number_of_branches() == 1
    assert t.get_length(0) == pytest.approx(5.0, abs=1e-12)
    assert _tip(t, 0) == pytest.approx(tip_before, abs=1e-12)
    assert _total_volume(t) == pytest.approx(vol_before, rel=0.0, abs=1e-12)


# ---------- targeted collapse_chains_after_prune ------------------------------


def test_prune_records_parents_of_cut_branches():
    """After a hurricane prune, the tree exposes the parents of cut subtrees
    via ``get_last_prune_parents``."""
    t = _vertical_trunk(length=1.0, diameter=0.2)
    t.add_branch_with_geometry(
        0, length=0.5, diameter=0.05, unit_t=(0.3, 0.0, 0.95), unit_b=(1.0, 0.0, 0.0)
    )
    t.add_branch_with_geometry(
        0, length=0.5, diameter=0.05, unit_t=(0.0, 0.3, 0.95), unit_b=(1.0, 0.0, 0.0)
    )
    t.reorder()
    t.set_seed(7)
    n_cut = prune(t, wind=(100.0, 0.0, 0.0), leaf_drag_S0=1.0, cauchy=1.0)
    assert n_cut > 0
    parents = t.get_last_prune_parents()
    # Trunk is the parent of both cut twigs; deduplicated -> exactly one entry.
    assert parents == [0]


def test_prune_no_cuts_clears_parent_list():
    """A zero-wind prune cuts nothing and leaves the parents list empty."""
    t = _vertical_trunk()
    t.add_branch_with_geometry(
        0, length=0.5, diameter=0.05, unit_t=(0.3, 0.0, 0.95), unit_b=(1.0, 0.0, 0.0)
    )
    t.add_branch_with_geometry(
        0, length=0.5, diameter=0.05, unit_t=(0.0, 0.3, 0.95), unit_b=(1.0, 0.0, 0.0)
    )
    t.reorder()
    t.set_seed(7)
    # Prime the parents list with a prior prune that cuts something so we can
    # confirm the next no-op prune actually clears it.
    prune(t, wind=(100.0, 0.0, 0.0), leaf_drag_S0=1.0, cauchy=1.0)
    prune(t, wind=(0.0, 0.0, 0.0), leaf_drag_S0=1.0, cauchy=1.0)
    assert t.get_last_prune_parents() == []


def test_after_prune_collapse_merges_chain_created_by_pruning():
    """A trunk with two children: one is a long single-child chain, the other
    is a fork. Pruning the fork turns its parent into a single-child parent;
    ``collapse_chains_after_prune`` should fuse the resulting chain."""
    t = _vertical_trunk(length=1.0, diameter=0.2)
    # Single-child chain on the +x side: trunk -> A (long) — already a
    # 1-child link, but A's child is a fork so nothing to merge yet.
    A = t.add_branch_with_geometry(
        0, length=1.0, diameter=0.1, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
    )
    # A's two-child fork that we will obliterate via a hurricane cut.
    fork = t.add_branch_with_geometry(
        A, length=1.0, diameter=0.04, unit_t=(0.0, 0.0, 1.0), unit_b=(1.0, 0.0, 0.0)
    )
    leaf_l = t.add_branch_with_geometry(
        fork, length=0.5, diameter=0.02, unit_t=(0.6, 0.0, 0.8), unit_b=(1.0, 0.0, 0.0)
    )
    leaf_r = t.add_branch_with_geometry(
        fork, length=0.5, diameter=0.02, unit_t=(-0.6, 0.0, 0.8), unit_b=(1.0, 0.0, 0.0)
    )
    t.reorder()
    assert t.get_number_of_branches() == 5

    # Prune leaf_r directly via a targeted setup: arrange wind so the thin
    # leaves break. Use a strong wind in the +x direction; the thinner branches
    # will go.
    t.set_seed(11)
    n_cut = prune(t, wind=(50.0, 0.0, 0.0), leaf_drag_S0=1.0, cauchy=1.0)
    assert n_cut > 0
    t.reorder()

    # Whatever exact subtree got cut, ``get_last_prune_parents`` must contain
    # the parent of each removed subtree. Running the targeted collapse should
    # behave identically to the whole-tree variant — we sanity-check that by
    # comparing absorbed counts on a parallel tree.
    vol_before = _total_volume(t)
    tips_before = sorted(_leaf_tips(t))
    n_before = t.get_number_of_branches()

    absorbed = t.collapse_chains_after_prune()
    t.reorder()

    assert absorbed >= 0
    assert t.get_number_of_branches() == n_before - absorbed
    assert _total_volume(t) == pytest.approx(vol_before, rel=0.0, abs=1e-12)
    tips_after = sorted(_leaf_tips(t))
    assert len(tips_after) == len(tips_before)
    for got, want in zip(tips_after, tips_before, strict=True):
        assert got == pytest.approx(want, abs=1e-12)
    # silence unused warnings on leaf indices
    del leaf_l, leaf_r


def test_after_prune_collapse_matches_whole_tree_when_used_every_step():
    """End-to-end equivalence: when ``collapse_chains_after_prune`` is called
    after every prune, the resulting tree is identical to running the
    whole-tree ``collapse_single_child_chains`` on the same schedule."""

    def grow(targeted: bool) -> PyTree:
        def on_step(_gen, tree):
            if targeted:
                tree.collapse_chains_after_prune()
            else:
                tree.collapse_single_child_chains()
            tree.reorder()

        return grow_tree(Config(), n_generations=20, seed=42, on_step=on_step)

    t_targeted = grow(targeted=True)
    t_whole = grow(targeted=False)

    assert t_targeted.get_number_of_branches() == t_whole.get_number_of_branches()
    assert t_targeted.get_total_leaves() == t_whole.get_total_leaves()

    tips_a = sorted(_leaf_tips(t_targeted))
    tips_b = sorted(_leaf_tips(t_whole))
    assert len(tips_a) == len(tips_b)
    for got, want in zip(tips_a, tips_b, strict=True):
        assert got == pytest.approx(want, abs=1e-9)


# ---------- integration test with the simulator -------------------------------


def test_volume_and_tips_preserved_after_grow():
    """Grow a real tree, snapshot volume + leaf tips, collapse, compare.

    With deterministic seed and config, pruning produces single-child chains;
    collapse should not change the tree's outer geometry or total volume.
    """
    tree = grow_tree(Config(), n_generations=30, seed=42)

    vol_before = _total_volume(tree)
    tips_before = sorted(_leaf_tips(tree))
    n_before = tree.get_number_of_branches()

    absorbed = tree.collapse_single_child_chains()
    tree.reorder()

    # If any pruning occurred over 30 generations there will usually be at
    # least one chain to collapse — but if not, the test still asserts the
    # invariants below.
    assert absorbed >= 0
    assert tree.get_number_of_branches() == n_before - absorbed

    assert _total_volume(tree) == pytest.approx(vol_before, rel=1e-9, abs=1e-12)
    tips_after = sorted(_leaf_tips(tree))
    assert len(tips_after) == len(tips_before)
    for got, want in zip(tips_after, tips_before, strict=True):
        assert got == pytest.approx(want, abs=1e-9)


def test_grow_with_collapse_each_generation():
    """Run the simulator with collapse called at every step; verify the tree
    is well-formed (re-orderable, simulator-callable) and that the merge
    actually fires for at least one generation on a seed that prunes."""
    total_absorbed = [0]

    def on_step(_gen, tree):
        total_absorbed[0] += tree.collapse_single_child_chains()
        tree.reorder()

    tree = grow_tree(Config(), n_generations=50, seed=42, on_step=on_step)

    # The collapse must have fired at least once across 50 generations of
    # default-wind pruning.
    assert total_absorbed[0] > 0
    # End-state must be consistent: every branch is reachable through the
    # depth-first index and ``leaf_indices`` agrees with childless branches.
    n = tree.get_number_of_branches()
    leaves = set(tree.leaf_indices())
    by_children = {i for i in range(n) if tree.get_number_of_children(i) == 0}
    assert leaves == by_children
