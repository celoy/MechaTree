"""Measure subtree removal throughput.

Builds a chain tree of N branches, then times `remove_branch(N//2)`
which removes the second half of the chain (~N/2 branches).
"""

import time

from mechatree import PyTree


def build_chain(n: int) -> PyTree:
    t = PyTree({"length": 1.0, "radius": 0.1})
    for i in range(n - 1):
        t.add_branch(i, {"length": 1.0, "radius": 0.05})
    return t


def main() -> None:
    sizes = [100, 1_000, 10_000, 50_000]
    print(f"{'N':>10}  {'remove (s)':>12}  {'per removed (us)':>18}")
    for n in sizes:
        t = build_chain(n)
        cut = n // 2
        remaining_before = t.get_number_of_branches()
        start = time.perf_counter()
        t.remove_branch(cut)
        elapsed = time.perf_counter() - start
        removed = remaining_before - t.get_number_of_branches()
        per_removed_us = elapsed * 1e6 / max(removed, 1)
        print(f"{n:>10}  {elapsed:>12.4f}  {per_removed_us:>18.3f}")


if __name__ == "__main__":
    main()
