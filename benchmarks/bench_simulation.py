"""Per-phase timings for one simulation generation.

Builds a balanced binary tree of N branches, then times each phase of the
single-generation pipeline (stress -> requested_growth -> secondary_growth
-> prune -> primary_growth). The shape isn't biologically representative —
it's a fixed-size benchmark so per-phase timings are comparable across runs.

Step 9 ships this as a data-collection script, not a CI gate. Compare
against a Fortran reference (legacy/fortran/tree.f90 on a matched seed)
when one is available.
"""

import math
import time

from mechatree import PyTree
from mechatree.genome import ConstantAllocation, ConstantSafety
from mechatree.growth import primary_growth, requested_growth, secondary_growth
from mechatree.mechanics import calculate_stresses
from mechatree.pruning import prune


def _build_balanced_tree(target_branches: int) -> PyTree:
    """Construct a balanced binary tree of approximately `target_branches`."""
    t = PyTree({})
    t.set_length(0, 1.0)
    t.set_diameter(0, 0.2)
    t.set_unit_t(0, (0.0, 0.0, 1.0))
    t.set_unit_b(0, (1.0, 0.0, 0.0))

    # Each generation doubles the leaf count. Stop when total ≥ target.
    while t.get_number_of_branches() < target_branches:
        # Snapshot leaves before mutating (add_branch_with_geometry shifts
        # indices and we don't want to grow off the new twigs).
        leaves = t.leaf_indices()
        for leaf_idx in sorted(leaves, reverse=True):
            if t.get_number_of_branches() >= target_branches:
                break
            t.add_branch_with_geometry(
                leaf_idx,
                length=0.5,
                diameter=0.05,
                unit_t=(0.3, 0.0, 0.95),
                unit_b=(1.0, 0.0, 0.0),
            )
            if t.get_number_of_branches() >= target_branches:
                break
            t.add_branch_with_geometry(
                leaf_idx,
                length=0.5,
                diameter=0.05,
                unit_t=(0.0, 0.3, 0.95),
                unit_b=(1.0, 0.0, 0.0),
            )
    t.reorder()
    return t


def _seed_state(tree: PyTree) -> None:
    """Set per-branch light = 1.0 on every leaf and a non-zero reserve."""
    for idx in tree.leaf_indices():
        tree.set_light(idx, 1.0)
    tree.set_reserve(0.1)
    tree.set_seed(42)


def _time(fn, repeats: int = 5) -> float:
    """Best-of-`repeats` wall time, in milliseconds."""
    best = math.inf
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        dt = (time.perf_counter() - t0) * 1e3
        best = min(best, dt)
    return best


def main() -> None:
    sizes = [127, 511, 2047, 4095]  # 2^k - 1 — close to balanced binary sizes
    print(
        f"{'N':>6}  {'stress':>10}  {'req_growth':>10}  {'sec_growth':>10}  "
        f"{'prune':>10}  {'primary':>10}  {'total':>10}"
    )
    print(
        f"{'':6}  {'(ms)':>10}  {'(ms)':>10}  {'(ms)':>10}  "
        f"{'(ms)':>10}  {'(ms)':>10}  {'(ms)':>10}"
    )

    safety = ConstantSafety(1.0)
    alloc = ConstantAllocation(p_seeds=0.0, p_leaves=0.5, phototropism=0.0)

    for n in sizes:

        def fresh(n=n):
            t = _build_balanced_tree(n)
            _seed_state(t)
            return t

        # Non-mutating phases: best-of-5 on a single tree.
        t = fresh()
        t_stress = _time(lambda t=t: calculate_stresses(t, 0.5, 1.0))

        t = fresh()
        calculate_stresses(t, 0.5, 1.0)
        t_req = _time(lambda t=t: requested_growth(t, safety, 0.005))

        t = fresh()
        calculate_stresses(t, 0.5, 1.0)
        requested_growth(t, safety, 0.005)
        t_sec = _time(lambda t=t: secondary_growth(t, 0.01))

        # Mutating phases: pre-build a batch and time the phase across them.
        # Build cost stays out of the measurement; we report best-of-batch.
        prepped = []
        for _ in range(3):
            tt = fresh()
            calculate_stresses(tt, 0.5, 1.0)
            prepped.append(tt)
        best = math.inf
        for tt in prepped:
            t0 = time.perf_counter()
            prune(tt, (1.0, 0.0, 0.0), 0.5, 1.0)
            best = min(best, (time.perf_counter() - t0) * 1e3)
        t_prune = best

        prepped = [fresh() for _ in range(3)]
        best = math.inf
        for tt in prepped:
            t0 = time.perf_counter()
            primary_growth(tt, alloc, 0.3, 0.02, 0.25, -0.25, 0.0, math.pi, 0)
            best = min(best, (time.perf_counter() - t0) * 1e3)
        t_prim = best

        total = t_stress + t_req + t_sec + t_prune + t_prim
        print(
            f"{n:>6}  {t_stress:>10.3f}  {t_req:>10.3f}  {t_sec:>10.3f}  "
            f"{t_prune:>10.3f}  {t_prim:>10.3f}  {total:>10.3f}"
        )


if __name__ == "__main__":
    main()
