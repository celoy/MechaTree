"""Per-tree statistical diagnostics — Python port of the MATLAB plots in
``../Eloy2017_NatComm_archive/``:

- ``strahler_summary`` mirrors ``plot_stat_single_tree.m`` and the Strahler-
  order tables emitted by ``mod_tools.f90`` ``save_statistics``.
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


__all__ = ["StrahlerSummary", "leonardo_ratios", "strahler_summary", "tokunaga_matrix"]
