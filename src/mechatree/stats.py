"""Per-tree statistical diagnostics — Python port of the MATLAB plots in
``../Eloy2017_NatComm_archive/``:

- ``strahler_summary`` mirrors ``plot_stat_single_tree.m`` and the Strahler-
  order tables emitted by ``mod_tools.f90`` ``save_statistics``.
- ``horton_ratios`` mirrors ``Fractal_dim.m`` — log-linear fit across
  Strahler ranks recovering the bifurcation / length / diameter / area
  ratios and the fractal dimension ``D = log R_n / log R_l``.
- ``leonardo_ratios`` mirrors ``plot_area_preservation_1tree.m`` — the
  child-area / parent-area ratio at each binary branching node, used to
  test Leonardo da Vinci's rule that a parent branch's cross-section is
  preserved across its children's cross-sections.
- ``tokunaga_matrix`` mirrors the ``Z_Tokunaga_*.dat`` file: for each
  pair of Strahler orders ``(i, j)`` with ``j < i``, the count of order-``j``
  side branches attached to order-``i`` branches.

All diagnostics operate on a ``PyTree`` in memory. ``set_strahler()``
must have been called first (or the helpers will call it for you).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from mechatree._core import PyTree


@dataclass
class StrahlerSummary:
    """Per-Strahler-order summary for one tree.

    Each NumPy array is indexed 0..max_order-1; entry ``k`` corresponds to
    Strahler order ``k + 1``.
    """

    n_branches: np.ndarray  # (max_order,) int — count per order
    mean_length: np.ndarray  # (max_order,) — mean length per order
    mean_diameter: np.ndarray  # (max_order,)
    mean_area: np.ndarray  # (max_order,) — pi/4 * d^2 averaged
    total_area: np.ndarray  # (max_order,) — sum of cross-sections
    max_order: int


def strahler_summary(tree: PyTree) -> StrahlerSummary:
    """Compute per-Strahler-order length / diameter / area / count tables."""
    tree.set_strahler()
    n = tree.get_number_of_branches()
    if n == 0:
        empty = np.zeros((0,))
        return StrahlerSummary(
            n_branches=empty.astype(np.int64),
            mean_length=empty,
            mean_diameter=empty,
            mean_area=empty,
            total_area=empty,
            max_order=0,
        )

    strahler = np.empty(n, dtype=np.int64)
    length = np.empty(n)
    diameter = np.empty(n)
    for i in range(n):
        strahler[i] = tree.get_strahler(i)
        length[i] = tree.get_length(i)
        diameter[i] = tree.get_diameter(i)

    max_order = int(strahler.max())
    counts = np.zeros(max_order, dtype=np.int64)
    sum_len = np.zeros(max_order)
    sum_diam = np.zeros(max_order)
    sum_area = np.zeros(max_order)
    area = 0.25 * math.pi * diameter * diameter
    for k in range(1, max_order + 1):
        mask = strahler == k
        counts[k - 1] = int(mask.sum())
        sum_len[k - 1] = float(length[mask].sum())
        sum_diam[k - 1] = float(diameter[mask].sum())
        sum_area[k - 1] = float(area[mask].sum())

    mean_len = np.divide(sum_len, counts, out=np.zeros_like(sum_len), where=counts > 0)
    mean_diam = np.divide(sum_diam, counts, out=np.zeros_like(sum_diam), where=counts > 0)
    mean_area = np.divide(sum_area, counts, out=np.zeros_like(sum_area), where=counts > 0)

    return StrahlerSummary(
        n_branches=counts,
        mean_length=mean_len,
        mean_diameter=mean_diam,
        mean_area=mean_area,
        total_area=sum_area,
        max_order=max_order,
    )


@dataclass
class HortonSummary:
    """Per-Horton-order stream summary for one tree.

    A Horton **stream** of order ``w`` is a maximal chain of consecutive
    branches that all carry Horton index ``w`` (i.e. an unbranched
    sequence of branches that share an order). Each NumPy array is
    indexed 0..max_order-1; entry ``k`` corresponds to Horton order
    ``k + 1``.

    Field semantics:

    * ``n_branches`` — number of **streams** of order w (NOT the number
      of segments). Equals ``Horton_distribution[w]`` from the C++ core.
    * ``mean_length`` — mean stream length per order: sum of segment
      lengths over segments with horton == w, divided by the number of
      streams of order w. This is the Horton/Strahler length the paper
      reports — in MechaTree every segment is a unit twig, so the
      mean *segment* length is always 1 and only the per-stream mean
      varies with order.
    * ``mean_diameter``, ``mean_area`` — per-segment mean over segments
      belonging to streams of order w (an average over chain members,
      not over chains).
    """

    n_branches: np.ndarray  # (max_order,) int — stream count per order
    mean_length: np.ndarray  # (max_order,) — mean stream length per order
    mean_diameter: np.ndarray  # (max_order,)
    mean_area: np.ndarray  # (max_order,)
    total_area: np.ndarray  # (max_order,) — sum of cross-sections over segments
    max_order: int


def horton_summary(tree: PyTree) -> HortonSummary:
    """Per-Horton-order stream statistics.

    Use this for the architectural-ratio analysis (R_n, R_l, R_d, R_a, D)
    that reproduces SI Fig. S8 of Eloy et al. 2017. The complementary
    :func:`strahler_summary` reports per-segment means — useful for
    histograms of segment properties but not for Horton's length ratio.
    """
    # set_horton() only recomputes Strahler when Strahler_distribution is
    # empty — so on a tree that has grown since the last set_strahler the
    # Horton labels would be derived from stale Strahler values. Force a
    # refresh before letting set_horton run.
    tree.set_strahler()
    tree.set_horton()
    n = tree.get_number_of_branches()
    if n == 0:
        empty = np.zeros((0,))
        return HortonSummary(
            n_branches=empty.astype(np.int64),
            mean_length=empty,
            mean_diameter=empty,
            mean_area=empty,
            total_area=empty,
            max_order=0,
        )

    horton = np.empty(n, dtype=np.int64)
    length = np.empty(n)
    diameter = np.empty(n)
    for i in range(n):
        horton[i] = tree.get_horton(i)
        length[i] = tree.get_length(i)
        diameter[i] = tree.get_diameter(i)

    max_order = int(horton.max())
    dist = tree.get_horton_distribution()
    n_streams = np.zeros(max_order, dtype=np.int64)
    for k, v in dist.items():
        if 1 <= k <= max_order:
            n_streams[k - 1] = int(v)

    sum_len = np.zeros(max_order)
    sum_diam = np.zeros(max_order)
    sum_area = np.zeros(max_order)
    seg_count = np.zeros(max_order, dtype=np.int64)
    area = 0.25 * math.pi * diameter * diameter
    for k in range(1, max_order + 1):
        mask = horton == k
        seg_count[k - 1] = int(mask.sum())
        sum_len[k - 1] = float(length[mask].sum())
        sum_diam[k - 1] = float(diameter[mask].sum())
        sum_area[k - 1] = float(area[mask].sum())

    mean_len = np.divide(sum_len, n_streams, out=np.zeros_like(sum_len), where=n_streams > 0)
    mean_diam = np.divide(sum_diam, seg_count, out=np.zeros_like(sum_diam), where=seg_count > 0)
    mean_area = np.divide(sum_area, seg_count, out=np.zeros_like(sum_area), where=seg_count > 0)

    return HortonSummary(
        n_branches=n_streams,
        mean_length=mean_len,
        mean_diameter=mean_diam,
        mean_area=mean_area,
        total_area=sum_area,
        max_order=max_order,
    )


@dataclass
class HortonRatios:
    """Geometric ratios across consecutive Strahler ranks for one tree.

    ``R_n`` is the bifurcation ratio ``N_w / N_{w+1}``; ``R_l``, ``R_d``,
    ``R_a`` are the length / diameter / area ratios ``X_{w+1} / X_w``. All
    four are > 1 for real trees. ``D`` is the Horton-Strahler fractal
    dimension ``log R_n / log R_l``. ``fit_orders`` records the Strahler
    ranks that contributed to the log-linear fit.
    """

    R_n: float
    R_l: float
    R_d: float
    R_a: float
    D: float
    fit_orders: np.ndarray


def mean_stream_length(tree: PyTree) -> np.ndarray:
    """Mean Horton-Strahler stream length per Strahler rank, in unit twigs.

    Reproduces the Fortran ``length_Strahler`` array
    (``mod_tree.f90:1091-1097``, written to ``Z_length_*.dat``), which the
    MATLAB ``Fractal_dim.m`` plots as the per-rank ``<l>`` trace and feeds
    to the ``R_l`` fit in SI Fig. S8(b).

    For each Strahler rank ``w``::

        L(w) = (# segments with Strahler == w) / (# streams of order w)

    Numerator from :func:`strahler_summary`; denominator from
    :func:`horton_strahler_counts`. In a clean binary tree every leaf is a
    single-segment order-1 stream so ``L(1) = 1``; in MechaTree
    chain-merger artifacts let ``L(1)`` drift slightly above 1, and higher
    ranks recover the geometric growth ``L(w+1) / L(w) = R_l ≈ 1.7``
    that the paper reports.

    This is the right length series for the ``R_l`` fit (the Fortran's own
    output). The radial-distance metric :func:`mean_distance_to_leaves`
    biases the slope ~8% high in MechaTree and pushes the recovered
    fractal dimension ``D = log R_n / log R_l`` ~10% low.

    Returns a 1D ``np.float64`` array of length ``max_strahler``; entry
    ``k`` is the value for rank ``k + 1``.
    """
    n_w = horton_strahler_counts(tree)
    if n_w.size == 0:
        return np.zeros(0)
    seg_counts = strahler_summary(tree).n_branches.astype(np.float64)
    # Both arrays are indexed by Strahler rank starting at 1; pad the
    # shorter one if a rank is missing from either (shouldn't happen since
    # both call set_strahler, but guard anyway).
    W = max(seg_counts.size, n_w.size)
    if seg_counts.size < W:
        seg_counts = np.concatenate([seg_counts, np.zeros(W - seg_counts.size)])
    if n_w.size < W:
        n_w = np.concatenate([n_w, np.zeros(W - n_w.size, dtype=n_w.dtype)])
    return np.divide(
        seg_counts, n_w.astype(np.float64), out=np.zeros(W, dtype=np.float64), where=n_w > 0
    )


def horton_strahler_counts(tree: PyTree) -> np.ndarray:
    """Per-Horton-Strahler-stream count N(w), matching the paper.

    Reproduces the Fortran formula in ``mod_tree.f90:1052-1070`` (writer of
    ``Z_Nsegments_*.dat``, the file the paper's SI Fig. S8(b) plots as
    "number of branches per Strahler order"):

    - ``N(1) = number of leaves`` (terminal branches).
    - ``N(w+1) = number of internal branches whose two children have
      identical Strahler order w``.

    This is the classical Horton-Strahler **stream** count, not the
    per-segment Strahler count. In a unit-twig representation (MechaTree
    and the Fortran reference) the per-segment Strahler count at high
    ranks is dominated by long chains of equal-order segments along the
    main trunk and overcounts the topology — e.g. the trunk's order-W
    main path may contain dozens of unit segments yet only one bifurcation
    actually creates a new order-W stream.

    Returns a 1D ``np.int64`` array of length ``max_strahler``; entry
    ``k`` is the count for rank ``k + 1``.
    """
    tree.set_strahler()
    n = tree.get_number_of_branches()
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    strahler = np.empty(n, dtype=np.int64)
    for i in range(n):
        strahler[i] = tree.get_strahler(i)
    max_order = int(strahler.max())
    counts = np.zeros(max_order, dtype=np.int64)
    n_leaves = 0
    for i in range(n):
        kids = tree.get_children_index(i)
        if not kids:
            n_leaves += 1
            continue
        if len(kids) != 2:
            continue
        s1 = int(strahler[kids[0]])
        s2 = int(strahler[kids[1]])
        if s1 == s2:
            # New stream of order s1 + 1 starts at this bifurcation.
            counts[s1] += 1  # 0-indexed entry for rank s1 + 1
    counts[0] = int(n_leaves)
    return counts


def horton_ratios(
    summary: HortonSummary | StrahlerSummary,
    *,
    drop_top: bool = True,
    max_rank: int | None = None,
    mean_length_override: np.ndarray | None = None,
    n_branches_override: np.ndarray | None = None,
) -> HortonRatios:
    """Log-linear fit of per-rank quantities → bifurcation ratios + fractal D.

    Mirrors ``Fractal_dim.m`` from the Eloy et al. 2017 MATLAB archive.
    For each per-rank series ``{n_branches, mean_length, mean_diameter,
    mean_area}`` the log of the quantity is regressed on rank ``w``; the
    matching geometric ratio is recovered as ``10**|slope|``. The fractal
    dimension follows as ``D = log(R_n) / log(R_l)``.

    The figure-S8 ratios in Eloy et al. 2017 are computed from a
    :func:`horton_summary` (per-Horton-stream view: ``n_branches`` is the
    stream count, ``mean_length`` is the per-stream length). A
    :class:`StrahlerSummary` will also work via duck typing but its
    ``mean_length`` is per-segment, which is constant in MechaTree, so
    the recovered ``R_l`` collapses to 1.

    ``drop_top=True`` (default) excludes the highest rank because
    ``N_W = 1`` by construction (the single root stream), which flattens
    the bifurcation slope; the means at rank W are also noisy (one data
    point) so dropping the top rank keeps the four fits anchored to the
    same support.

    ``max_rank`` caps the fit support at ``min(W - drop_top, max_rank)``.
    Eloy et al. 2017's SI Fig. S12 fits only ranks 1..7 even on trees
    that grow to max Strahler 8–9: high-rank means are noisy because
    they sit on a handful of branches, and capping the fit keeps the
    recovered ratios stable across runs of differing depth.

    ``mean_length_override`` substitutes a different per-rank length
    array (typically :func:`mean_distance_to_leaves`) for the length fit,
    leaving the count / diameter / area fits unchanged. Must have the
    same length as ``summary.n_branches``. Reproduces the paper's
    ``R_l`` from the recursive distance-to-leaves rather than the
    per-stream chain length that ``HortonSummary.mean_length`` reports.

    ``n_branches_override`` substitutes a different per-rank count array
    for the ``R_n`` fit (typically :func:`horton_strahler_counts`, the
    Fortran reference). Required when the markers in the figure are
    drawn from a different series than ``summary.n_branches``, so the
    fit line passes through the markers rather than a parallel series.

    Raises ``ValueError`` if fewer than two ranks remain after the
    drop_top exclusion and masking of zero counts.
    """
    W = summary.max_order
    if W < 2:
        raise ValueError(f"need at least 2 Strahler ranks, got {W}")
    last = W - 1 if drop_top else W
    if max_rank is not None:
        if max_rank < 2:
            raise ValueError(f"max_rank must be >= 2, got {max_rank}")
        last = min(last, max_rank)
    if last < 2:
        raise ValueError(f"need at least 2 usable ranks, got {last} (try drop_top=False)")
    w_all = np.arange(1, last + 1, dtype=float)
    if n_branches_override is not None:
        if n_branches_override.shape != summary.n_branches.shape:
            raise ValueError(
                "n_branches_override shape "
                f"{n_branches_override.shape} does not match summary "
                f"{summary.n_branches.shape}"
            )
        counts = n_branches_override[:last].astype(float)
    else:
        counts = summary.n_branches[:last].astype(float)
    if mean_length_override is not None:
        if mean_length_override.shape != summary.n_branches.shape:
            raise ValueError(
                "mean_length_override shape "
                f"{mean_length_override.shape} does not match summary "
                f"{summary.n_branches.shape}"
            )
        length = mean_length_override[:last].astype(float)
    else:
        length = summary.mean_length[:last].astype(float)
    diameter = summary.mean_diameter[:last].astype(float)
    area = summary.mean_area[:last].astype(float)
    mask = (counts > 0) & (length > 0) & (diameter > 0) & (area > 0)
    if int(mask.sum()) < 2:
        raise ValueError(f"need at least 2 positive-count ranks, got {int(mask.sum())}")
    w = w_all[mask]

    def slope(y: np.ndarray) -> float:
        m, _ = np.polyfit(w, np.log10(y[mask]), 1)
        return float(m)

    R_n = 10.0 ** (-slope(counts))
    R_l = 10.0 ** slope(length)
    R_d = 10.0 ** slope(diameter)
    R_a = 10.0 ** slope(area)
    D = math.log(R_n) / math.log(R_l)
    return HortonRatios(R_n=R_n, R_l=R_l, R_d=R_d, R_a=R_a, D=D, fit_orders=w)


def distance_to_leaves(tree: PyTree) -> np.ndarray:
    """Per-branch recursive distance from each branch to its descendant twigs.

    Terminal branches (no children) get ``length / 2``. Internal branches
    get ``length + weighted_mean(child.distance_to_leaves)`` across their
    children, weighted by each child subtree's ``nb_leaves``. The reference
    point is the **base** of the branch, and the metric is the average
    arc-length distance along the tree skeleton to the midpoint of a
    descendant terminal twig.

    Mirrors ``b%distance_leaves`` in the legacy Fortran ``save_area``
    (``legacy_fortran/mod_tree.f90:1174-1203``); the Fortran constants
    ``0.5`` and ``+1.0`` are the unit-twig instantiation (``length = 1``)
    of this length-aware form. Carrying ``b.length`` lets the metric stay
    correct when branches are merged after pruning and grow past unit length.

    Used to reproduce the per-rank ``<l>`` plotted in Eloy et al. 2017
    SI Fig. S8(b); aggregate with :func:`mean_distance_to_leaves` to
    obtain the per-Horton-rank series that feeds :func:`horton_ratios`.

    The MechaTree depth-first ordering guarantees every descendant has a
    higher index than its parent (see ``mechanics.cpp:98-100``), so a single
    reverse-index pass suffices.
    """
    tree.reorder()  # refresh nb_leaves in case the tree hasn't been stepped
    n = tree.get_number_of_branches()
    if n == 0:
        return np.zeros(0)
    dist = np.empty(n)
    for i in range(n - 1, -1, -1):
        L = tree.get_length(i)
        kids = tree.get_children_index(i)
        if len(kids) == 0:
            dist[i] = 0.5 * L
            continue
        num = 0.0
        denom = 0
        for k in kids:
            w = tree.get_nb_leaves(k)
            num += dist[k] * w
            denom += w
        dist[i] = L + (num / denom if denom > 0 else 0.5 * L)
    return dist


def mean_distance_to_leaves(tree: PyTree) -> np.ndarray:
    """Per-Horton-rank mean of :func:`distance_to_leaves`.

    For each Horton order ``w`` (1-indexed), returns the **mean over all
    branches** with ``horton == w`` of their per-branch
    ``distance_to_leaves`` value. The output is a ``(max_order,)`` array
    indexed ``0..max_order-1`` so it can be dropped in as
    ``mean_length_override`` for :func:`horton_ratios` alongside a
    :func:`horton_summary` of the same tree.

    Matches the per-rank ``<l>`` reported in Eloy et al. 2017 SI Fig. S8(b)
    — the Fortran pipeline writes per-branch ``distance_leaves`` paired
    with Strahler order and the MATLAB ``Fractal_dim.m`` script then
    averages over branches per rank.
    """
    tree.set_strahler()
    tree.set_horton()
    n = tree.get_number_of_branches()
    if n == 0:
        return np.zeros(0)
    dist = distance_to_leaves(tree)
    horton = np.empty(n, dtype=np.int64)
    for i in range(n):
        horton[i] = tree.get_horton(i)
    max_order = int(horton.max())
    sums = np.zeros(max_order)
    counts = np.zeros(max_order, dtype=np.int64)
    for w in range(1, max_order + 1):
        mask = horton == w
        counts[w - 1] = int(mask.sum())
        sums[w - 1] = float(dist[mask].sum())
    return np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0)


def leonardo_ratios(tree: PyTree) -> np.ndarray:
    """Area-preservation ratios at each binary branching node.

    Returns an array of ``(child1_area + child2_area) / parent_area``,
    one entry per branch that has exactly two children. For
    cross-section-preserving branching (Leonardo's rule) the ratio is
    ~1; for area-preservation under the WBE / pipe model it's > 1.

    Branches with 0, 1, or >2 children are skipped.
    """
    ratios: list[float] = []
    n = tree.get_number_of_branches()
    for i in range(n):
        kids = tree.get_children_index(i)
        if len(kids) != 2:
            continue
        d_p = tree.get_diameter(i)
        if d_p <= 0.0:
            continue
        d_a = tree.get_diameter(kids[0])
        d_b = tree.get_diameter(kids[1])
        a_p = d_p * d_p
        a_ab = d_a * d_a + d_b * d_b
        ratios.append(a_ab / a_p)
    return np.asarray(ratios, dtype=np.float64)


def tokunaga_matrix(tree: PyTree) -> np.ndarray:
    """Tokunaga side-branching matrix.

    Entry ``[i-1, j-1]`` (``i > j >= 1``) counts how many order-``j``
    branches sit as direct children of order-``i`` branches across the
    whole tree. Diagonal and upper-triangle are zero.

    Returns an ``(max_order, max_order)`` int64 array. For a perfectly
    self-similar Horton-Strahler tree, ``T_{i,j}`` follows the geometric
    progression ``T_k * R^k`` characteristic of the Tokunaga generator
    (Tokunaga 1978, Pelletier & Turcotte 2000).
    """
    tree.set_strahler()
    n = tree.get_number_of_branches()
    if n == 0:
        return np.zeros((0, 0), dtype=np.int64)

    max_order = int(max(tree.get_strahler(i) for i in range(n)))
    T = np.zeros((max_order, max_order), dtype=np.int64)
    for i in range(n):
        s_i = tree.get_strahler(i)
        for child_idx in tree.get_children_index(i):
            s_c = tree.get_strahler(child_idx)
            if s_c < s_i:
                T[s_i - 1, s_c - 1] += 1
    return T


__all__ = [
    "HortonRatios",
    "HortonSummary",
    "StrahlerSummary",
    "distance_to_leaves",
    "horton_ratios",
    "horton_strahler_counts",
    "horton_summary",
    "leonardo_ratios",
    "mean_distance_to_leaves",
    "mean_stream_length",
    "strahler_summary",
    "tokunaga_matrix",
]
