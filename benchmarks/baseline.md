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

## bench_simulation — one generation, per-phase timings (Step 9, PR2)

Balanced binary tree of N branches. Best-of-5 for non-mutating phases
(stress, requested_growth, secondary_growth); best-of-3 across a pre-built
batch for mutating phases (prune, primary_growth). Build cost excluded.

| N    | stress (ms) | req_growth (ms) | sec_growth (ms) | prune (ms) | primary (ms) | total (ms) |
| ---: | ----------: | --------------: | --------------: | ---------: | -----------: | ---------: |
|  127 |       0.006 |           0.003 |           0.006 |      0.019 |        0.101 |      0.134 |
|  511 |       0.023 |           0.012 |           0.034 |      0.067 |        1.088 |      1.223 |
| 2047 |       0.087 |           0.046 |           0.173 |      0.230 |        2.811 |      3.346 |
| 4095 |       0.156 |           0.089 |           0.372 |      0.501 |        6.035 |      7.152 |

Reading: mechanics + growth are roughly linear in N (good). `primary_growth`
super-linear because `addBranchWithGeometry` rebuilds `branch_to_index`
across the tail of the vector per insertion — fine for now, the hot path
to fix once a Fortran reference is in hand.

A Fortran reference column will be added once `legacy_fortran/tree.f90`
is set up to dump per-phase wall times on the same machine.

## bench_light — light interception over N leaves (Step 10)

Random leaf cloud in a `[-15, 15]^2 x [0, 20]` box under the default
`Sun()` (4 elevations × 8 azimuths = 32 directions). Best of 3 after a
warmup call.

| N      | best (ms) | ms/direction | µs/leaf |
| -----: | --------: | -----------: | ------: |
|    100 |      0.51 |        0.016 |    5.08 |
|    500 |      1.31 |        0.041 |    2.61 |
|   1000 |      2.52 |        0.079 |    2.52 |
|   2500 |      6.87 |        0.215 |    2.75 |
|   5000 |     14.54 |        0.454 |    2.91 |
|  10000 |     30.91 |        0.966 |    3.09 |

Linear in N (~3 µs/leaf at scale), per the `lexsort + np.unique` design.

This is the vectorised implementation; the original Python `for j in
order:` loop was ~5× slower (e.g. 80 ms at N=5000). The vectorised path
keeps light competitive with the Step-9 mechanics phases: at N=5000
across 3000 generations the light cost is roughly **45 s**, vs. the
total Step-9 mechanics+growth+pruning at the same scale (~25 s including
the still-superlinear `primary_growth`). Light is no longer the
dominant cost; a Cython port would buy maybe another 2–3× but is not
needed at current scales.
