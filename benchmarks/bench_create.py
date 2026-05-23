"""Measure tree-creation throughput.

Builds chain trees of varying N branches and reports wall time per branch.
Expected scaling: O(N) — anything worse indicates a hot path that's
quadratic in tree size.
"""

import time

from mechatree import PyTree


def build_chain(n: int) -> None:
    t = PyTree({"length": 1.0, "radius": 0.1})
    for i in range(n - 1):
        t.add_branch(i, {"length": 1.0, "radius": 0.05})
    # touch a getter so the optimizer can't elide the work
    assert t.get_number_of_branches() == n


def main() -> None:
    sizes = [100, 1_000, 10_000, 100_000]
    print(f"{'N':>10}  {'total (s)':>12}  {'per branch (us)':>18}")
    for n in sizes:
        start = time.perf_counter()
        build_chain(n)
        elapsed = time.perf_counter() - start
        per_branch_us = elapsed * 1e6 / n
        print(f"{n:>10}  {elapsed:>12.4f}  {per_branch_us:>18.3f}")


if __name__ == "__main__":
    main()
