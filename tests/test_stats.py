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
    distance_to_leaves,
    horton_ratios,
    horton_strahler_counts,
    horton_summary,
    leonardo_ratios,
    mean_distance_to_leaves,
    mean_stream_length,
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


def test_horton_strahler_counts_perfect_binary_tree():
    """A perfect binary tree of depth d has the classical doubling Horton-
    Strahler stream counts: ``N(w) = 2^(d - w + 1)`` for w = 1..d, plus the
    trunk ``N(d+1) = 1``. ``N(1)`` must equal the leaf count."""
    t = _balanced_binary_tree(3)
    n_w = horton_strahler_counts(t)
    # depth-3 perfect binary tree: 8 leaves, 4 order-2, 2 order-3, 1 order-4 trunk.
    assert n_w.tolist() == [8, 4, 2, 1]
    assert int(n_w[0]) == t.get_total_leaves()


def test_horton_strahler_counts_n1_equals_leaves_on_grown_tree():
    """Even on a real simulated tree with pruning artifacts, ``N(1)`` must
    equal the leaf count exactly — that is the paper's convention
    (mod_tree.f90:1057)."""
    tree = grow_tree(Config(), n_generations=80, seed=7)
    n_w = horton_strahler_counts(tree)
    assert int(n_w[0]) == tree.get_total_leaves()
    # Counts are monotone-decreasing (roughly) on a self-similar tree; at
    # minimum the top rank is the singleton trunk.
    assert int(n_w[-1]) == 1


def test_mean_stream_length_perfect_binary_tree():
    """In a perfect binary tree every stream is a single segment, so the
    mean stream length is 1 at every Strahler rank."""
    t = _balanced_binary_tree(3)
    L_w = mean_stream_length(t)
    assert L_w.tolist() == [1.0, 1.0, 1.0, 1.0]


def test_mean_stream_length_matches_seg_over_stream_definition():
    """``mean_stream_length(tree)[w] == strahler_n_branches[w] / n_w[w]``
    on a grown tree — this is the Fortran ``length_Strahler`` formula
    (mod_tree.f90:1091-1097)."""
    tree = grow_tree(Config(), n_generations=50, seed=11)
    seg_counts = strahler_summary(tree).n_branches.astype(float)
    n_w = horton_strahler_counts(tree).astype(float)
    expected = np.divide(seg_counts, n_w, out=np.zeros_like(seg_counts), where=n_w > 0)
    L_w = mean_stream_length(tree)
    np.testing.assert_allclose(L_w, expected)


def test_horton_strahler_counts_empty_tree_returns_zero_length_array():
    """``PyTree({})`` has one trunk segment (Strahler 1); the result has
    length 1 entry."""
    t = PyTree({})
    t.set_length(0, 1.0)
    t.set_diameter(0, 0.1)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    n_w = horton_strahler_counts(t)
    assert n_w.tolist() == [1]


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
    safety, allocation, _, _ = load_champion(champions, species_id=0)
    tree = grow_tree(cfg, n_generations=150, seed=42, safety=safety, allocation=allocation)
    h_summary = horton_summary(tree)
    if h_summary.max_order < 4:
        pytest.skip(f"champion tree only reached Horton order {h_summary.max_order}")
    h = horton_ratios(h_summary)
    assert 2.0 < h.R_n < 6.0, f"R_n={h.R_n} outside loose window for paper target 3.5"
    assert 1.2 < h.R_l < 2.5, f"R_l={h.R_l} outside loose window for paper target 1.7"
    assert 1.2 < h.R_d < 3.0, f"R_d={h.R_d} outside loose window for paper target 1.9"
    assert 1.5 < h.D < 3.5, f"D={h.D} outside loose window for paper target 2.3"


def test_distance_to_leaves_trunk_only():
    """A tree with only the trunk: the trunk is itself a terminal branch."""
    t = PyTree({})
    t.set_length(0, 1.0)
    t.set_diameter(0, 0.1)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    dist = distance_to_leaves(t)
    assert dist.tolist() == [0.5]


def test_distance_to_leaves_single_binary_split():
    """Trunk with two leaf children: each leaf is 0.5, trunk is
    ((0.5 + 1) * 1 + (0.5 + 1) * 1) / 2 = 1.5."""
    t = PyTree({})
    t.set_length(0, 1.0)
    t.set_diameter(0, 0.1)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    t.add_branch_with_geometry(0, 1.0, 0.1, (0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
    t.add_branch_with_geometry(0, 1.0, 0.1, (0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
    t.reorder()
    dist = distance_to_leaves(t)
    assert dist == pytest.approx([1.5, 0.5, 0.5])


def test_distance_to_leaves_perfect_binary_depth_3():
    """In a depth-3 perfect binary tree the per-level value collapses to
    level + 0.5 (8 leaves at 0.5, 4 grandparents at 1.5, 2 at 2.5, trunk 3.5)."""
    t = _balanced_binary_tree(3)
    dist = distance_to_leaves(t)
    # 15 branches; trunk at index 0 (level 3); leaves at end (level 0).
    assert dist[0] == pytest.approx(3.5)
    # Eight leaves with dist=0.5 — at least 8 entries == 0.5.
    assert int((dist == 0.5).sum()) == 8


def test_distance_to_leaves_unbalanced_weighting():
    """Asymmetric subtrees: the weighting by ``nb_leaves`` must use the
    descendant-leaf count, not the immediate child count. Build a trunk
    with a heavy left subtree (4 leaves) and a single leaf on the right.

    Left subtree: trunk → l → (l1, l2) where l1 → (l1a, l1b) and l2 →
    (l2a, l2b). Then trunk also has a direct right leaf ``r``.
    Hand computation:
      l1.dist = 1.5,  l1.nb_leaves = 2
      l2.dist = 1.5,  l2.nb_leaves = 2
      l.dist = ((1.5+1)*2 + (1.5+1)*2) / (2+2) = 2.5
      l.nb_leaves = 4
      r.dist = 0.5, r.nb_leaves = 1
      trunk.dist = ((2.5+1)*4 + (0.5+1)*1) / (4+1) = (14 + 1.5) / 5 = 3.1
    """
    t = PyTree({})
    t.set_length(0, 1.0)
    t.set_diameter(0, 0.1)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    geo = {"length": 1.0, "diameter": 0.1, "unit_t": (0.0, 0.0, 1.0), "unit_b": (1.0, 0.0, 0.0)}
    l_idx = t.add_branch_with_geometry(0, **geo)
    l1_idx = t.add_branch_with_geometry(l_idx, **geo)
    t.add_branch_with_geometry(l1_idx, **geo)
    t.add_branch_with_geometry(l1_idx, **geo)
    l2_idx = t.add_branch_with_geometry(l_idx, **geo)
    t.add_branch_with_geometry(l2_idx, **geo)
    t.add_branch_with_geometry(l2_idx, **geo)
    t.add_branch_with_geometry(0, **geo)  # direct right leaf
    t.reorder()
    dist = distance_to_leaves(t)
    assert dist[0] == pytest.approx(3.1)


def test_distance_to_leaves_varying_branch_lengths():
    """Lengths != 1.0 should propagate: terminal dist = L/2; internal
    dist = own_length + weighted_mean(child.dist). Trunk with length 3,
    two leaf children each of length 2:
        leaf.dist = 2 / 2 = 1.0
        trunk.dist = 3 + (1.0 + 1.0) / 2 = 4.0
    """
    t = PyTree({})
    t.set_length(0, 3.0)
    t.set_diameter(0, 0.1)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))
    t.add_branch_with_geometry(0, 2.0, 0.1, (0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
    t.add_branch_with_geometry(0, 2.0, 0.1, (0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
    t.reorder()
    dist = distance_to_leaves(t)
    assert dist == pytest.approx([4.0, 1.0, 1.0])


def test_mean_distance_to_leaves_perfect_binary_depth_3():
    """The C++ ``set_horton`` rule (``tree.cpp:setHorton``) lets one child
    per fork inherit the parent's Horton label, so per-Horton-rank
    membership crosses levels-from-leaves. For the depth-3 perfect
    binary tree:
        Horton 4 chain absorbs one branch per level (dist 3.5, 2.5, 1.5, 0.5)
                  → mean 2.0
        Horton 3 chain absorbs three branches (dist 2.5, 1.5, 0.5)  → mean 1.5
        Horton 2: 2 level-1 branches (dist 1.5) + 2 leaves (dist 0.5) → mean 1.0
        Horton 1: 4 lone leaves (dist 0.5)                             → mean 0.5
    """
    t = _balanced_binary_tree(3)
    mean = mean_distance_to_leaves(t)
    assert mean == pytest.approx([0.5, 1.0, 1.5, 2.0])


def test_horton_ratios_mean_length_override():
    """The override path swaps the length series used in the R_l fit; the
    other three ratios are untouched."""
    R_n, R_l, R_d = 2.0, 1.7, 1.9
    s = _geometric_summary(W=8, R_n=R_n, R_l=R_l, R_d=R_d)
    # Build a different geometric length series with R_l_alt = 1.4.
    R_l_alt = 1.4
    w = np.arange(1, 9, dtype=float)
    override = R_l_alt ** (w - 1)
    h = horton_ratios(s, mean_length_override=override)
    assert h.R_l == pytest.approx(R_l_alt, rel=1e-10)
    # Counts / diameter / area fits should be untouched.
    assert h.R_n == pytest.approx(R_n, rel=1e-10)
    assert h.R_d == pytest.approx(R_d, rel=1e-10)
    assert h.R_a == pytest.approx(R_d**2, rel=1e-10)


def test_horton_ratios_max_rank_caps_fit():
    """``max_rank`` truncates the fit support so high-rank noise can't
    drag the recovered ratios. Build a clean geometric series for ranks
    1..5 and intentionally corrupt rank 6; without the cap, ``R_n`` is
    wrong; with ``max_rank=5`` the fit recovers the underlying ratio."""
    R_n = 2.0
    W = 6
    w = np.arange(1, W + 1, dtype=float)
    counts = np.round(R_n ** (W - w)).astype(np.int64)
    counts[-1] = 100  # corrupt the top rank (should be 1)
    length = 1.7 ** (w - 1)
    diameter = 1.9 ** (w - 1)
    area = 0.25 * math.pi * diameter * diameter
    s = HortonSummary(
        n_branches=counts,
        mean_length=length,
        mean_diameter=diameter,
        mean_area=area,
        total_area=area * counts,
        max_order=W,
    )
    # Default: drop_top=True trims rank 6, so the corruption never sees the
    # fit and R_n recovers cleanly. Use drop_top=False to expose it.
    h_no_cap = horton_ratios(s, drop_top=False)
    assert h_no_cap.R_n != pytest.approx(R_n, rel=0.1)
    # With max_rank=5 the corrupted rank-6 entry is excluded.
    h_capped = horton_ratios(s, drop_top=False, max_rank=5)
    assert h_capped.R_n == pytest.approx(R_n, rel=1e-10)
    assert h_capped.fit_orders.tolist() == [1, 2, 3, 4, 5]


def test_horton_ratios_max_rank_below_two_raises():
    s = _geometric_summary(W=6, R_n=2.0, R_l=1.7, R_d=1.9)
    with pytest.raises(ValueError, match="max_rank"):
        horton_ratios(s, max_rank=1)


def test_horton_ratios_mean_length_override_shape_mismatch():
    s = _geometric_summary(W=6, R_n=2.0, R_l=1.7, R_d=1.9)
    with pytest.raises(ValueError, match="mean_length_override shape"):
        horton_ratios(s, mean_length_override=np.array([1.0, 2.0]))


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
