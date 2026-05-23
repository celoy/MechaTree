# Benchmark numbers

Run `uv run python benchmarks/bench_*.py` to reproduce. Numbers below were
captured on an Apple Silicon Mac (Darwin 25.4, Python 3.12.6) — your
mileage will vary in absolute terms, but the **scaling** is what matters.

Phase 0 = unmodified port from `archive/PyTreeLib/`.
Phase 4 = modernized core with O(1) `branch_to_index` lookup, single-pointer
parent links, and `unordered_map` for branch properties.

## bench_create — chain tree of N branches

| N       | Phase 0 (s) | Phase 0 (µs/branch) | Phase 4 (s) | Phase 4 (µs/branch) |
| ------: | ----------: | ------------------: | ----------: | ------------------: |
|     100 |      0.0001 |               0.760 |      0.0001 |               0.646 |
|   1 000 |      0.0006 |               0.580 |      0.0005 |               0.524 |
|  10 000 |      0.0059 |               0.589 |      0.0060 |               0.604 |
| 100 000 |      0.0536 |               0.536 |      0.0688 |               0.688 |

Roughly flat per-branch in both columns (chain tree → append at end →
no shift cost). The slight Phase-4 uptick at 100k comes from maintaining
the `branch_to_index` map on every insert.

## bench_remove — remove subtree at N/2 of an N-branch chain

| N      | Phase 0 (s) | Phase 0 (µs/removed) | Phase 4 (s) | Phase 4 (µs/removed) |
| -----: | ----------: | -------------------: | ----------: | -------------------: |
|    100 |      0.0000 |                2.167 |      0.0000 |                1.167 |
|  1 000 |      0.0000 |                8.875 |      0.0000 |                5.375 |
| 10 000 |      0.0001 |              116.375 |      0.0001 |               67.875 |
| 50 000 |      0.0012 |             1179.875 |      0.0004 |              367.958 |

Phase 0 total wall time at N=50 000 was **3.0× slower** than Phase 4
(1.2 ms vs 0.4 ms). The remaining linear-in-N cost in Phase 4 is the
O(subtree) `shift_indices` walk that keeps `branch_to_index` in sync
after the erase.

## bench_topology — parent-chain walk leaf → trunk on an N-branch chain

| N     | Phase 0 (s) | Phase 0 (µs/step) | Phase 4 (s) | Phase 4 (µs/step) |
| ----: | ----------: | ----------------: | ----------: | ----------------: |
|   100 |      0.0000 |             0.070 |      0.0000 |             0.044 |
|   500 |      0.0001 |             0.132 |      0.0000 |             0.045 |
| 1 000 |      0.0002 |             0.213 |      0.0001 |             0.068 |
| 2 000 |      0.0007 |             0.368 |      0.0001 |             0.053 |
| 5 000 |      0.0040 |             0.803 |      0.0003 |             0.056 |

The headline win: **µs/step is now constant** instead of growing linearly
with N — `Tree::getIndex` is O(1) average via the `branch_to_index`
unordered_map, so `getParentIndex` is too. Walking 5 000 steps dropped
from 4.0 ms to 0.3 ms (**~14× faster**); the gap widens with N.

## bench_memory — RSS over 100 create+drop iterations of 1 000-branch trees

| iter | Phase 0 maxrss (MB) | Phase 0 Δ (MB) | Phase 4 Δ (MB) |
| ---: | ------------------: | -------------: | -------------: |
|   10 |               17.70 |           2.47 |           0.05 |
|   20 |               20.09 |           4.86 |           0.25 |
|   30 |               22.50 |           7.27 |           0.25 |
|   40 |               24.92 |           9.69 |           0.25 |
|   50 |               27.31 |          12.08 |           0.25 |
|   60 |               29.72 |          14.48 |           0.25 |
|   70 |               32.12 |          16.89 |           0.25 |
|   80 |               34.53 |          19.30 |           0.25 |
|   90 |               36.92 |          21.69 |           0.25 |
|  100 |               39.34 |          24.11 |           0.25 |

Phase 0 leaked ~241 KB per 1 000-branch tree (perfectly linear growth).
Phase 4 plateaus at +0.25 MB across 100 iterations — that's allocator
slack from the warm-up phase, not a leak. Roughly **100× less** total
memory growth.
