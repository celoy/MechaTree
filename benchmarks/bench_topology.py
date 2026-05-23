"""Measure the leaf -> trunk parent-chain walk.

This exercises the O(N) Tree::getIndex linear scan inside every
Tree::getParentIndex call, making the whole walk O(N^2). Phase 4 of the
modernization adds a Branch*->index lookup map; this benchmark should drop
by orders of magnitude after.
"""

import time

from mechatree import PyTree


def build_chain(n: int) -> PyTree:
    t = PyTree({"length": 1.0, "radius": 0.1})
    for i in range(n - 1):
        t.add_branch(i, {"length": 1.0, "radius": 0.05})
    return t


def parent_chain_walk(t: PyTree, leaf_index: int) -> int:
    idx = leaf_index
    steps = 0
    while idx > 0:
        idx = t.get_parent_index(idx)
        steps += 1
    return steps


def main() -> None:
    sizes = [100, 500, 1_000, 2_000, 5_000]
    print(f"{'N':>8}  {'walk (s)':>12}  {'us/step':>10}")
    for n in sizes:
        t = build_chain(n)
        start = time.perf_counter()
        steps = parent_chain_walk(t, n - 1)
        elapsed = time.perf_counter() - start
        us_per_step = elapsed * 1e6 / max(steps, 1)
        print(f"{n:>8}  {elapsed:>12.4f}  {us_per_step:>10.3f}")


if __name__ == "__main__":
    main()
