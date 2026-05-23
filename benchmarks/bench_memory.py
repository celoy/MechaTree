"""Measure peak RSS growth across create+drop iterations.

Pre-Phase-2: Tree / Branch lack destructors and PyTree has no __dealloc__,
so each iteration permanently retains the C++ tree on the heap. ru_maxrss
grows monotonically with iteration count.

Post-Phase-2: RSS should plateau at one tree's footprint.
"""

import gc
import platform
import resource

from mechatree import PyTree


def maxrss_bytes() -> int:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux/BSD report kilobytes.
    if platform.system() == "Darwin":
        return rss
    return rss * 1024


def build_chain(n: int) -> PyTree:
    t = PyTree({"length": 1.0, "radius": 0.1})
    for i in range(n - 1):
        t.add_branch(i, {"length": 1.0, "radius": 0.05})
    return t


def main() -> None:
    n_branches = 1_000
    n_iter = 100

    # Warm-up so allocator pages are mapped.
    for _ in range(3):
        build_chain(n_branches)
    gc.collect()

    baseline = maxrss_bytes()
    print(f"baseline maxrss: {baseline / 1024 / 1024:.2f} MB")
    print(f"{'iter':>6}  {'maxrss (MB)':>14}  {'delta (MB)':>14}")

    for i in range(1, n_iter + 1):
        t = build_chain(n_branches)
        del t
        gc.collect()
        if i % (n_iter // 10) == 0:
            current = maxrss_bytes()
            delta = (current - baseline) / 1024 / 1024
            print(f"{i:>6}  {current / 1024 / 1024:>14.2f}  {delta:>14.2f}")


if __name__ == "__main__":
    main()
