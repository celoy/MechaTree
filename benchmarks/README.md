# MechaTree micro-benchmarks

Standalone runnable scripts (not pytest) for measuring `_core` performance
and memory behavior. Use them to capture before/after numbers when changing
the C++/Cython core.

## Running

```bash
uv run python benchmarks/bench_create.py
uv run python benchmarks/bench_remove.py
uv run python benchmarks/bench_topology.py
uv run python benchmarks/bench_memory.py
```

Each script prints a small table to stdout. Capture numbers in
[baseline.md](baseline.md) when measuring a change.

## What each script measures

| Script | Measures | Hot path it exercises |
| --- | --- | --- |
| `bench_create.py` | Wall time for building chain trees of N branches | `Tree::addBranch`, `Branch` allocation, property map |
| `bench_remove.py` | Wall time for `remove_branch` on a full subtree | `Tree::removeBranch`, `getLastDescendantIndex` |
| `bench_topology.py` | Wall time for the leaf→root parent-chain walk | `Tree::getParentIndex` → `Tree::getIndex` (O(N) scan) |
| `bench_memory.py` | Process RSS growth over create+drop iterations | Tree / Branch destructors (or lack thereof) |

## Reading the numbers

- `bench_create` and `bench_remove` should be roughly **linear in N** if
  the implementation is well-behaved.
- `bench_topology` is **O(N²)** in the current code (every
  `get_parent_index` does a linear scan of `tree_branches`); Phase 4 of
  the modernization plan adds an O(1) `Branch*→index` map and should
  bring it to O(N).
- `bench_memory` exposes the destructor leaks: pre-Phase-2, RSS grows
  monotonically across iterations; post-Phase-2, it plateaus.
