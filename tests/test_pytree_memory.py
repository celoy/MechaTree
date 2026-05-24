"""Regression test for the C++ memory leak in `_core`.

Tree allocates each Branch with `new` but neither Tree nor Branch defines a
destructor, and PyTree has no __dealloc__ — so every PyTree leaks its entire
branch graph on garbage collection (tree.cpp:867-876, _core.pyx:20-41).

This test loops create+drop on many trees and checks the process's peak RSS
does not balloon. Marked xfail strict today; the xfail is removed in Phase 2
once the destructors land.
"""

import gc
import platform

import pytest

# `resource` is a Unix-only stdlib module — skip the whole file on Windows.
resource = pytest.importorskip("resource")

from mechatree import PyTree  # noqa: E402


def _maxrss_bytes() -> int:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports ru_maxrss in bytes; Linux/BSD report in kilobytes.
    if platform.system() == "Darwin":
        return rss
    return rss * 1024


def _build_tree(n_branches: int) -> PyTree:
    t = PyTree({"length": 1.0, "radius": 0.1})
    for i in range(n_branches - 1):
        t.add_branch(i, {"length": 1.0, "radius": 0.05})
    return t


@pytest.mark.slow
def test_repeated_create_does_not_grow_rss():
    n_branches = 1_000
    n_iter = 50

    # Warm up the allocator so the first iteration's overhead doesn't skew
    # the baseline.
    for _ in range(3):
        _build_tree(n_branches)
    gc.collect()

    baseline = _maxrss_bytes()

    for _ in range(n_iter):
        t = _build_tree(n_branches)
        del t
        gc.collect()

    after = _maxrss_bytes()
    growth_bytes = after - baseline

    # Per-branch overhead is ~200 bytes (Branch + vectors + property map),
    # so n_iter * n_branches * 200 ~= 10 MB if everything leaks. Threshold
    # at 4 MB so the test cleanly distinguishes "leaking" from "not".
    assert growth_bytes < 4 * 1024 * 1024, (
        f"RSS grew by {growth_bytes / 1024 / 1024:.1f} MB across "
        f"{n_iter} create+drop cycles of {n_branches}-branch trees — "
        f"memory is leaking."
    )
