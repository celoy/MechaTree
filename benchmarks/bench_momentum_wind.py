"""Benchmark the native momentum-wind wind kernel.

Measures wall-clock vs (grid size, forest size) on a warmed-up
forest. Use this to gate optimisations: run the benchmark, apply a
change to ``mechatree.wind._momentum_wind_kernel``, re-run, compare.

Usage
-----
::

    uv run python benchmarks/bench_momentum.py            # standard sweep
    uv run python benchmarks/bench_momentum.py --profile  # cProfile output
    uv run python benchmarks/bench_momentum.py --quick    # tiny sweep

The standard sweep sizes match what Notebook 07 uses (~5k branches,
~58k grid cells) plus a 10× and 100× scale.
"""

from __future__ import annotations

import argparse
import time
from typing import NamedTuple

import numpy as np

import mechatree as mt
from mechatree.wind._momentum_wind_kernel import compute_momentum_wind


class _Probe(NamedTuple):
    """One bench point: forest setup + grid setup."""

    label: str
    n_trees_init: int
    n_trees_max: int
    forest_size: float
    n_gens: int
    grid_size: float
    pad_x: float


PROBES = [
    _Probe(
        label="notebook (5k branches, 58k cells)",
        n_trees_init=8,
        n_trees_max=80,
        forest_size=15.0,
        n_gens=30,
        grid_size=1.0,
        pad_x=12.0,
    ),
    _Probe(
        label="10x denser (~30k branches, 80k cells)",
        n_trees_init=20,
        n_trees_max=200,
        forest_size=15.0,
        n_gens=60,
        grid_size=1.0,
        pad_x=12.0,
    ),
    _Probe(
        label="fine grid (~5k branches, 350k cells, grid_size=0.5)",
        n_trees_init=8,
        n_trees_max=80,
        forest_size=15.0,
        n_gens=30,
        grid_size=0.5,
        pad_x=12.0,
    ),
]


def _build_forest(p: _Probe) -> mt.Forest:
    cfg = mt.Config(
        tree=mt.TreeConfig(),
        forest=mt.ForestConfig(
            size=p.forest_size, n_trees_init=p.n_trees_init, n_trees_max=p.n_trees_max
        ),
        n_generations=p.n_gens,
    )
    f = mt.Forest(cfg, seed=0)
    for g in range(p.n_gens):
        f.step(g)
    return f


def _build_inputs(forest: mt.Forest, p: _Probe):
    """Pool per-branch geometry + build the grid for the bench probe."""
    starts, axes, Ds, Ls = [], [], [], []
    for t in forest.trees:
        s, a, d, ell = t.get_branch_data_batch()
        starts.append(s)
        axes.append(a)
        Ds.append(d)
        Ls.append(ell)
    start = np.concatenate(starts)
    axis = np.concatenate(axes)
    D = np.concatenate(Ds)
    L = np.concatenate(Ls)

    x_lo = start[:, 0].min() - p.pad_x
    x_hi = start[:, 0].max() + p.pad_x
    y_lo = start[:, 1].min() - 2.0
    y_hi = start[:, 1].max() + 2.0
    z_hi = (start[:, 2] + L * axis[:, 2]).max() + 3.0

    grid_size = p.grid_size
    bx = np.arange(x_lo, x_hi + grid_size, grid_size)
    by = np.arange(y_lo, y_hi + grid_size, grid_size)
    bz = np.arange(0.0, z_hi + grid_size, grid_size)
    z_centers = 0.5 * (bz[:-1] + bz[1:])
    ua, z0, kappa = 0.4, 0.1, 0.41
    U_infty = (ua / kappa) * np.log(np.maximum(z_centers, z0) / z0)
    return start, axis, D, L, bx, by, bz, U_infty


def _time_once(start, axis, D, L, bx, by, bz, U_infty, grid_size, repeats: int = 5):
    """Time `compute_momentum_wind` and return (best, mean) seconds."""
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        compute_momentum_wind(
            start,
            axis,
            D,
            L,
            cell_bounds_x=bx,
            cell_bounds_y=by,
            cell_bounds_z=bz,
            grid_size=grid_size,
            U_infty=U_infty,
        )
        times.append(time.perf_counter() - t0)
    return min(times), sum(times) / len(times)


def run_probe(p: _Probe, repeats: int = 5) -> None:
    forest = _build_forest(p)
    start, axis, D, L, bx, by, bz, U_infty = _build_inputs(forest, p)
    Nx, Ny, Nz = bx.size - 1, by.size - 1, bz.size - 1
    n_branches = start.shape[0]
    best, mean = _time_once(start, axis, D, L, bx, by, bz, U_infty, p.grid_size, repeats)
    print(
        f"  {p.label}: branches={n_branches}, "
        f"grid=({Nx}x{Ny}x{Nz})={Nx * Ny * Nz}, "
        f"best={best * 1000:.1f}ms, mean={mean * 1000:.1f}ms"
    )


def run_all(repeats: int = 5) -> None:
    print(f"momentum-wind bench — repeats={repeats}")
    for p in PROBES:
        run_probe(p, repeats=repeats)


def run_profile() -> None:
    import cProfile
    import pstats

    p = PROBES[1]  # bigger probe → more signal
    forest = _build_forest(p)
    start, axis, D, L, bx, by, bz, U_infty = _build_inputs(forest, p)

    pr = cProfile.Profile()
    pr.enable()
    for _ in range(3):
        compute_momentum_wind(
            start,
            axis,
            D,
            L,
            cell_bounds_x=bx,
            cell_bounds_y=by,
            cell_bounds_z=bz,
            grid_size=p.grid_size,
            U_infty=U_infty,
        )
    pr.disable()
    s = pstats.Stats(pr).sort_stats("cumulative")
    s.print_stats(20)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick", action="store_true", help="single small probe")
    ap.add_argument("--profile", action="store_true", help="cProfile output")
    ap.add_argument("--repeats", type=int, default=5, help="timer repeats per probe")
    args = ap.parse_args()

    if args.profile:
        run_profile()
        return
    if args.quick:
        run_probe(PROBES[0], repeats=args.repeats)
        return
    run_all(repeats=args.repeats)


if __name__ == "__main__":
    main()
