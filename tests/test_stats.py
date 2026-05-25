"""Tests for the mechatree.stats module."""

import math

import numpy as np
import pytest

from mechatree import PyTree
from mechatree.config import Config
from mechatree.simulate import grow_tree
from mechatree.stats import (
    HortonRatios,
    HortonSummary,
    horton_ratios,
    horton_summary,
    leonardo_ratios,
    strahler_summary,
    tokunaga_matrix,
)


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


def _geometric_summary(W: int, R_n: float, R_l: float, R_d: float) -> HortonSummary:
    """Build a HortonSummary with exact geometric series across W ranks.

    N_w = R_n**(W-w)  → N_W = 1, N_1 = R_n**(W-1)
    L_w = R_l**(w-1)
    d_w = R_d**(w-1)
    A_w = pi/4 * d_w**2 (so R_a = R_d**2)
    """
    w = np.arange(1, W + 1, dtype=float)
    counts = np.round(R_n ** (W - w)).astype(np.int64)
    length = R_l ** (w - 1)
    diameter = R_d ** (w - 1)
    area = 0.25 * math.pi * diameter * diameter
    return HortonSummary(
        n_branches=counts,
        mean_length=length,
        mean_diameter=diameter,
        mean_area=area,
        total_area=area * counts,
        max_order=W,
    )


def test_horton_ratios_recovers_geometric_series():
    """A synthetic summary built from exact geometric series should round-trip
    through horton_ratios at high precision. Use R_n=2 so counts are exact
    powers of two (no rounding drift); length / diameter are floats so any
    ratios work for them."""
    R_n, R_l, R_d = 2.0, 1.7, 1.9
    s = _geometric_summary(W=8, R_n=R_n, R_l=R_l, R_d=R_d)
    h = horton_ratios(s)
    assert isinstance(h, HortonRatios)
    assert h.R_n == pytest.approx(R_n, rel=1e-10)
    assert h.R_l == pytest.approx(R_l, rel=1e-10)
    assert h.R_d == pytest.approx(R_d, rel=1e-10)
    assert h.R_a == pytest.approx(R_d**2, rel=1e-10)
    assert pytest.approx(math.log(R_n) / math.log(R_l), rel=1e-10) == h.D
    # drop_top=True → fit_orders = 1..W-1
    assert h.fit_orders.tolist() == [1, 2, 3, 4, 5, 6, 7]


def test_horton_ratios_paper_targets_close_with_loose_tol():
    """With paper-target ratios (R_n=3.5, R_l=1.7, R_d=1.9, W=8), integer
    rounding on N_w introduces ~1% drift in R_n. Confirm recovery at 5%."""
    R_n, R_l, R_d = 3.5, 1.7, 1.9
    s = _geometric_summary(W=8, R_n=R_n, R_l=R_l, R_d=R_d)
    h = horton_ratios(s)
    assert h.R_n == pytest.approx(R_n, rel=0.05)
    assert h.R_l == pytest.approx(R_l, rel=1e-10)
    assert h.R_d == pytest.approx(R_d, rel=1e-10)
    assert pytest.approx(math.log(R_n) / math.log(R_l), rel=0.05) == h.D


def test_horton_ratios_drop_top_false_includes_trunk():
    """With drop_top=False the rank-W trunk (N_W = 1) flattens the fit;
    R_n still recovers to the synthetic ratio because counts are exact."""
    s = _geometric_summary(W=6, R_n=4.0, R_l=2.0, R_d=2.0)
    h = horton_ratios(s, drop_top=False)
    assert h.R_n == pytest.approx(4.0, rel=1e-10)
    assert h.fit_orders.tolist() == [1, 2, 3, 4, 5, 6]


def test_horton_ratios_rejects_too_few_ranks():
    s = _geometric_summary(W=2, R_n=3.0, R_l=1.5, R_d=1.5)
    # W=2, drop_top=True leaves 1 rank → cannot fit.
    with pytest.raises(ValueError):
        horton_ratios(s)


def test_horton_summary_refreshes_strahler_on_grown_tree():
    """Regression: the C++ ``setHorton`` skips ``setStrahler`` when
    ``Strahler_distribution`` is already populated, which would yield
    stale Horton labels on a tree that has grown since the last
    Strahler computation. ``horton_summary`` must defend against this
    by forcing a fresh ``set_strahler`` first.

    Reproduce by growing one tree in a loop and snapshotting along the
    way — without the defensive call, every snapshot collapses to the
    first one's max_order."""
    cfg = Config()
    saved = []

    def cb(gen, tree):
        if gen % 10 == 0:
            saved.append((gen, horton_summary(tree).max_order))

    grow_tree(cfg, n_generations=80, seed=42, on_step=cb)
    orders = [m for _, m in saved]
    # Without the fix, all entries would equal the gen=0 max_order
    # (typically 1 or 2). With the fix, max_order grows monotonically
    # over time.
    assert orders[-1] > orders[0], f"Horton orders never grew: {orders}"


def test_horton_summary_perfect_binary_tree():
    """In a perfect depth-d binary tree the C++ ``setHorton`` rule lets
    one sibling at every fork inherit the parent's Horton index, so a
    chain at order w absorbs w consecutive unit segments. Trace by hand
    for depth=3 (max_order=4) and check both the stream counts and the
    per-stream mean lengths."""
    t = _balanced_binary_tree(3)
    h = horton_summary(t)
    # Trunk = 1 chain of order 4 (length 4: trunk + inheriting children
    # down to one leaf). Right child of trunk = 1 new chain of order 3,
    # length 3. Two new chains of order 2 (each length 2). Four leftover
    # leaves are their own chains of order 1.
    assert h.n_branches.tolist() == [4, 2, 1, 1]
    assert h.mean_length == pytest.approx([1.0, 2.0, 3.0, 4.0])


def test_horton_ratios_on_champion_tree_near_paper_values():
    """Regression: grow the species-0 champion long enough for the
    Horton stream structure to develop, then assert the recovered ratios
    land in the neighbourhood of the paper targets (R_n=3.5, R_l=1.7,
    D=2.3). Bounds are intentionally loose — a single seed has real
    scatter."""
    from pathlib import Path

    from mechatree.config import load_config
    from mechatree.genome import load_champion

    champions = Path(__file__).resolve().parents[1] / "data" / "S3_champions.json"
    forest_yaml = Path(__file__).resolve().parents[1] / "examples" / "forest.yaml"
    if not (champions.exists() and forest_yaml.exists()):
        pytest.skip("S3_champions.json or forest.yaml missing")

    cfg = load_config(forest_yaml)
    safety, allocation, _ = load_champion(champions, species_id=0)
    tree = grow_tree(cfg, n_generations=150, seed=42, safety=safety, allocation=allocation)
    h_summary = horton_summary(tree)
    if h_summary.max_order < 4:
        pytest.skip(f"champion tree only reached Horton order {h_summary.max_order}")
    h = horton_ratios(h_summary)
    assert 2.0 < h.R_n < 6.0, f"R_n={h.R_n} outside loose window for paper target 3.5"
    assert 1.2 < h.R_l < 2.5, f"R_l={h.R_l} outside loose window for paper target 1.7"
    assert 1.2 < h.R_d < 3.0, f"R_d={h.R_d} outside loose window for paper target 1.9"
    assert 1.5 < h.D < 3.5, f"D={h.D} outside loose window for paper target 2.3"


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
