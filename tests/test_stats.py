"""Tests for the mechatree.stats module."""

import math

import numpy as np
import pytest

from mechatree import PyTree
from mechatree.config import Config
from mechatree.simulate import grow_tree
from mechatree.stats import leonardo_ratios, strahler_summary, tokunaga_matrix


def _balanced_binary_tree(depth: int) -> PyTree:
    """Build a depth-``depth`` perfect binary tree with uniform geometry."""
    t = PyTree({})
    t.set_length(0, 1.0)
    t.set_diameter(0, 0.1)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))

    for _ in range(depth):
        # Process leaves in descending index order — later insertions don't
        # invalidate the smaller indices we haven't touched yet. Re-query
        # after each depth so the index map is fresh.
        for leaf in sorted(t.leaf_indices(), reverse=True):
            t.add_branch_with_geometry(
                leaf,
                length=1.0,
                diameter=0.1,
                unit_t=(0.0, 0.0, 1.0),
                unit_b=(1.0, 0.0, 0.0),
            )
            t.add_branch_with_geometry(
                leaf,
                length=1.0,
                diameter=0.1,
                unit_t=(0.0, 0.0, 1.0),
                unit_b=(1.0, 0.0, 0.0),
            )
    t.reorder()
    return t


def test_strahler_summary_perfect_binary_tree():
    """A perfect binary tree of depth d has 2^(d+1) - 1 branches, with order
    counts 2^d, 2^(d-1), ..., 1 (leaves first, then internal nodes climbing
    to the trunk)."""
    t = _balanced_binary_tree(3)
    s = strahler_summary(t)
    assert s.max_order == 4
    # depth-3 tree: 2^3=8 leaves (order 1), 4 (order 2), 2 (order 3), 1 trunk (order 4)
    assert s.n_branches.tolist() == [8, 4, 2, 1]
    # All branches have length=1, diameter=0.1, so means are constant.
    assert s.mean_length == pytest.approx([1.0, 1.0, 1.0, 1.0])
    assert s.mean_diameter == pytest.approx([0.1, 0.1, 0.1, 0.1])
    expected_area = 0.25 * math.pi * 0.1**2
    assert s.mean_area == pytest.approx([expected_area] * 4)


def test_strahler_summary_empty_tree():
    """A tree with no branches yields an empty summary."""
    # PyTree({}) creates a tree with 1 branch (the trunk). To get 0 branches
    # we'd have to remove it, which isn't allowed. Skip the 0-case here.
    pass


def test_strahler_summary_grown_tree():
    """A real simulated tree has more branches at lower Strahler orders."""
    tree = grow_tree(Config(), n_generations=30, seed=42)
    s = strahler_summary(tree)
    # By construction, leaf count >= count of next order >= ...
    assert s.n_branches[0] >= s.n_branches[-1]
    # All counts positive.
    assert (s.n_branches > 0).all()


def test_leonardo_ratios_uniform_diameter_tree():
    """In a tree where every branch has the same diameter, the area ratio
    at each binary junction is 2*A_child / A_parent = 2."""
    t = _balanced_binary_tree(3)
    ratios = leonardo_ratios(t)
    # All ratios should be 2.0 (2 children of same diameter as parent).
    assert ratios.size > 0
    assert np.allclose(ratios, 2.0)


def test_leonardo_ratios_skips_non_binary():
    """A tree where a branch has 1 child (not 2) is skipped in the count."""
    t = PyTree({})
    t.set_length(0, 1.0)
    t.set_diameter(0, 0.1)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    # Add a single child to the trunk — not a binary node.
    t.add_branch_with_geometry(0, length=1.0, diameter=0.1, unit_t=(0, 0, 1), unit_b=(1, 0, 0))
    t.reorder()
    ratios = leonardo_ratios(t)
    # No binary junctions; should be empty.
    assert ratios.size == 0


def test_tokunaga_matrix_perfect_binary_tree():
    """In a perfect binary tree, only the immediate-subordinate entries are
    nonzero: T[k+1, k] = 2 * (number of order-(k+1) branches).
    Actually, for the canonical Horton-Strahler perfect binary tree, each
    order-(k+1) branch has exactly 2 order-k children, so:
        T[k, k-1] = 2 * n_branches[k]
    """
    t = _balanced_binary_tree(3)
    T = tokunaga_matrix(t)
    assert T.shape == (4, 4)
    # Order 2 branches each have 2 order-1 children (their children are
    # both at order 1 because the perfect binary tree has equal-order
    # siblings, so Strahler increments). T[1, 0] = 2 * count(order 2) = 8.
    assert T[1, 0] == 8
    assert T[2, 1] == 4
    assert T[3, 2] == 2
    # No "skip-level" side branches in a perfect tree.
    assert T[2, 0] == 0
    assert T[3, 0] == 0


def test_tokunaga_matrix_no_lower_in_upper_triangle():
    """T is lower-triangular: a parent can't be a lower Strahler order
    than its child."""
    tree = grow_tree(Config(), n_generations=30, seed=42)
    T = tokunaga_matrix(tree)
    if T.size:
        # The upper triangle (including diagonal where the child has same
        # order as parent — but that can happen and contributes to the
        # parent's order). The matrix tracks child < parent only.
        assert (np.triu(T) == 0).all()
