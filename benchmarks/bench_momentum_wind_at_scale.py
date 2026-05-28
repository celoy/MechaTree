"""Forest-scale benchmark for the canopy-aware momentum-wind wind path.

Answers two operational questions:

1. Where does time go in a real ``Forest.step`` loop when the wind
   model is ``momentum`` (so the Step-24 fixed-point iterates
   each generation)?
2. How many wind iterations does each generation actually need —
   average + max + histogram?

Default probe: 1000 trees on a disk of radius R = 100, 100 generations.
Override via CLI flags.

Usage
-----
::

    uv run python benchmarks/bench_momentum_wind_at_scale.py
    uv run python benchmarks/bench_momentum_wind_at_scale.py --R 50 --n-trees 200 --n-gens 50

The script monkey-patches the per-phase functions in
:mod:`mechatree.forest` with a wall-clock decorator, so per-phase
totals come out automatically — same trick
:mod:`benchmarks.bench_forest_scale` uses.

Also runs a small ``(nu_diff, grid_size)`` sweep on the post-warmup forest
so the operator can pick parameters with eyes open.
"""

from __future__ import annotations

import argparse
import time
from collections import Counter, defaultdict
from typing import Any

import numpy as np

import mechatree as mt
from mechatree.wind._momentum_wind_kernel import compute_momentum_wind

# ---------------------------------------------------------------------------
# Per-phase timing instrumentation (mirrors bench_forest_scale style)
# ---------------------------------------------------------------------------

_PHASE_TOTALS: dict[str, float] = defaultdict(float)
_PHASE_COUNTS: dict[str, int] = defaultdict(int)


def _wrap_phase(module: Any, attr: str, label: str) -> None:
    """Wrap ``module.attr`` so each call adds to ``_PHASE_TOTALS[label]``."""
    original = getattr(module, attr)

    def wrapped(*args: Any, **kwargs: Any) -> Any:
        t0 = time.perf_counter()
        try:
            return original(*args, **kwargs)
        finally:
            _PHASE_TOTALS[label] += time.perf_counter() - t0
            _PHASE_COUNTS[label] += 1

    setattr(module, attr, wrapped)


def _install_instrumentation() -> None:
    """Wrap the Forest.step hot-path callees so wall-clock is attributed.

    Step 26c-aware: the momentum model senses via ``_sense_canopy`` (the
    CFD solves are timed on the bridge by :class:`_WindTracker`) and prunes
    via ``prune_with_stored_forces``; the legacy ``calculate_stresses`` is no
    longer imported into ``mechatree.forest``."""
    import mechatree.forest as F
    import mechatree.simulate as S

    _wrap_phase(F, "extract_leaves", "light.extract_leaves")
    _wrap_phase(F, "intercept", "light.intercept")
    _wrap_phase(F, "aggregate_onto_trees", "light.aggregate_onto_trees")
    _wrap_phase(F, "requested_growth", "growth.requested_growth")
    _wrap_phase(F, "secondary_growth", "growth.secondary_growth")
    _wrap_phase(F, "primary_growth", "growth.primary_growth")
    # Per-branch C++ aggregation passes (cheap; one per sensing angle / prune
    # sweep). The expensive CFD solves are timed separately on the bridge.
    _wrap_phase(S, "calculate_stresses_from_stored_forces", "mechanics.stress_agg(C++)")
    _wrap_phase(F, "prune_with_stored_forces", "pruning.cut(C++)")
    _wrap_phase(F, "prune", "pruning.prune(default)")  # only fires for model=default


# ---------------------------------------------------------------------------
# Wind-solve tracking: proxy the bridge, timing sensing vs pruning solves
# ---------------------------------------------------------------------------


class _WindTracker:
    """Proxy over the momentum bridge. Times the sensing (``sense``) and
    pruning (``__call__``) CFD solves separately, counts pruning fixed-point
    iterations per generation, and forwards everything else (``sensing_angles``,
    ``writes_segment_forces``, ``last_result``, …) to the wrapped bridge.

    (The old ``_IterTracker`` replaced ``forest.wind_fn`` with a plain
    callable, which broke the Step-26c sensing path that calls
    ``wind_fn.sense`` / ``wind_fn.sensing_angles``.)"""

    def __init__(self, bridge: Any) -> None:
        self._b = bridge
        self.per_gen: dict[int, int] = Counter()  # pruning fixed-point iters / gen
        self.sense_time = 0.0
        self.sense_calls = 0
        self.prune_time = 0.0
        self.prune_calls = 0

    def __call__(self, generation: int, rng: Any, context: Any) -> Any:
        t0 = time.perf_counter()
        out = self._b(generation, rng, context)
        self.prune_time += time.perf_counter() - t0
        self.prune_calls += 1
        self.per_gen[generation] += 1
        return out

    def sense(self, context: Any, theta: float) -> Any:
        t0 = time.perf_counter()
        out = self._b.sense(context, theta)
        self.sense_time += time.perf_counter() - t0
        self.sense_calls += 1
        return out

    def solve_directions(self, context: Any, thetas: Any) -> Any:
        # Step 26e: sensing now fans the n_sensing_angles solves out here (the
        # GIL-free kernel lets a thread pool overlap them). Time the whole
        # parallel sweep; count one "sense call" per angle so the per-solve
        # average stays comparable to the old per-angle ``sense`` timing.
        thetas = list(thetas)
        t0 = time.perf_counter()
        out = self._b.solve_directions(context, thetas)
        self.sense_time += time.perf_counter() - t0
        self.sense_calls += len(thetas)
        return out

    def sensing_angles(self, rng: Any, n: int) -> Any:
        return self._b.sensing_angles(rng, n)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._b, name)


def _summarise_iter_counts(per_gen: dict[int, int]) -> None:
    if not per_gen:
        print("  (no wind-bridge calls)")
        return
    counts = np.array(list(per_gen.values()), dtype=int)
    print(f"  generations with at least 1 iter:    {len(counts)}")
    print(f"  total wind-bridge calls:             {int(counts.sum())}")
    print(f"  mean iters / gen:                    {counts.mean():.2f}")
    print(f"  median iters / gen:                  {int(np.median(counts))}")
    print(f"  max iters / gen:                     {counts.max()}")
    print(
        "  iter-count histogram (iters: gens):  "
        + ", ".join(f"{k}: {v}" for k, v in sorted(Counter(counts.tolist()).items()))
    )


# ---------------------------------------------------------------------------
# Main bench
# ---------------------------------------------------------------------------


def run_forest_bench(R: float, n_trees: int, n_gens: int) -> None:
    print(f"\n=== Forest perf bench: R={R}, n_trees_init={n_trees}, n_gens={n_gens} ===")
    print("    wind model = momentum (canopy-aware ⇒ Step-24 loop)")

    cfg = mt.Config(
        tree=mt.TreeConfig(),
        forest=mt.ForestConfig(
            size=R,
            n_trees_init=n_trees,
            n_trees_max=int(n_trees * 1.5),
        ),
        wind=mt.WindConfig(
            model="momentum",
            grid_size=2.0,
            momentum_nu_diff=0.03,
            momentum_pad_x=12.0,
            momentum_pad_y=4.0,
            momentum_pad_z=3.0,
            n_sensing_angles=2,
            max_pruning_iterations=8,
            wind_convergence_eps_rel=0.01,
        ),
        n_generations=n_gens,
    )

    _install_instrumentation()
    forest = mt.Forest(cfg, seed=0)
    # Proxy the bridge so sensing vs pruning CFD solves are timed separately.
    tracker = _WindTracker(forest.wind_fn)
    forest.wind_fn = tracker

    t_start = time.perf_counter()
    for g in range(n_gens):
        forest.step(g)
    total = time.perf_counter() - t_start

    n_branches_final = sum(t.get_number_of_branches() for t in forest.trees)
    print(f"  final: {len(forest.trees)} trees, {n_branches_final} branches (after {n_gens} gens)")
    print(f"  wall-clock total: {total:.1f}s  ({total / n_gens * 1000:.0f} ms/gen mean)")
    print()
    print("  Momentum CFD solves (the dominant cost):")
    for label, t_sec, calls in (
        ("sensing solves (n_sensing_angles)", tracker.sense_time, tracker.sense_calls),
        ("pruning solves (fixed-point)", tracker.prune_time, tracker.prune_calls),
    ):
        pct = 100 * t_sec / total if total else 0.0
        per = t_sec / calls * 1000 if calls else 0.0
        print(
            f"    {label:36s} {t_sec * 1000:>9.0f} ms  ({pct:>4.1f}%)  "
            f"{calls:>6d} solves  ({per:.1f} ms/solve, {calls / n_gens:.1f}/gen)"
        )
    solve_total = tracker.sense_time + tracker.prune_time
    print(
        f"    {'→ all wind solves':36s} {solve_total * 1000:>9.0f} ms  "
        f"({100 * solve_total / total:>4.1f}% of wall-clock)"
    )
    print()
    print("  Per-phase wall-clock totals (non-solve):")
    rows = sorted(_PHASE_TOTALS.items(), key=lambda kv: -kv[1])
    for label, secs in rows:
        pct = 100 * secs / total
        calls = _PHASE_COUNTS[label]
        print(
            f"    {label:42s} {secs * 1000:>9.0f} ms  ({pct:>4.1f}%)  "
            f"{calls:>7d} calls  ({secs / calls * 1e6:.0f} µs/call)"
        )
    print()
    print("  Pruning fixed-point iteration distribution:")
    _summarise_iter_counts(tracker.per_gen)


def run_nu_diff_H_sweep(R: float, n_trees: int, warmup_gens: int) -> None:
    print(
        f"\n=== (nu_diff, grid_size) parameter sweep on R={R} forest "
        f"after {warmup_gens} gens of zero-wind warmup ==="
    )

    cfg = mt.Config(
        tree=mt.TreeConfig(),
        forest=mt.ForestConfig(
            size=R,
            n_trees_init=n_trees,
            n_trees_max=int(n_trees * 1.5),
        ),
        n_generations=warmup_gens,
    )

    def zero_wind(_g: int, _rng: Any) -> tuple[float, float, float]:
        return (0.0, 0.0, 0.0)

    forest = mt.Forest(cfg, seed=0, wind_fn=zero_wind)
    for g in range(warmup_gens):
        forest.step(g)

    # Pool geometry.
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
    print(f"  warmup state: {len(forest.trees)} trees, {start.shape[0]} branches")

    nu_diffs = [0.0, 0.01, 0.03, 0.1, 0.3]
    Hs = [0.5, 1.0, 2.0]
    print(
        f"  sweeping nu_diff ∈ {nu_diffs} × grid_size ∈ {Hs}; reporting "
        f"<U_branch> / <U_infty>, wake-spread (FWHM of deficit at downwind y-slice), "
        f"and wall-clock per call."
    )
    print(
        f"  {'nu_diff':>8s} {'grid_size':>5s} {'(Nx, Ny, Nz)':>16s} "
        f"{'⟨U_br/U_inf⟩':>13s} {'min U_br/U_inf':>15s} "
        f"{'wake FWHM (y/grid_size)':>16s} {'time (ms)':>10s}"
    )
    z_canopy_top = float((start[:, 2] + L * axis[:, 2]).max())
    for grid_size in Hs:
        x_lo = start[:, 0].min() - 6.0
        x_hi = start[:, 0].max() + 12.0
        y_lo = start[:, 1].min() - 4.0
        y_hi = start[:, 1].max() + 4.0
        z_hi = z_canopy_top + 3.0
        bx = np.arange(x_lo, x_hi + grid_size, grid_size)
        by = np.arange(y_lo, y_hi + grid_size, grid_size)
        bz = np.arange(0.0, z_hi + grid_size, grid_size)
        z_centers = 0.5 * (bz[:-1] + bz[1:])
        ua, z0, kappa = 0.4, 0.1, 0.41
        U_infty = (ua / kappa) * np.log(np.maximum(z_centers, z0) / z0)
        for nu_diff in nu_diffs:
            # Warm a JIT cache and measure best of 3.
            best = float("inf")
            res = None
            for _ in range(3):
                t0 = time.perf_counter()
                res = compute_momentum_wind(
                    start,
                    axis,
                    D,
                    L,
                    cell_bounds_x=bx,
                    cell_bounds_y=by,
                    cell_bounds_z=bz,
                    grid_size=grid_size,
                    U_infty=U_infty,
                    nu_diff=nu_diff,
                )
                best = min(best, time.perf_counter() - t0)
            # Per-branch ratio.
            mid_z = start[:, 2] + 0.5 * L * axis[:, 2]
            k_idx = np.clip(np.searchsorted(bz, mid_z, side="right") - 1, 0, bz.size - 2)
            U_inf_branch = U_infty[k_idx]
            ratio = res.U_branch / np.maximum(U_inf_branch, 1e-6)
            # Wake FWHM at downwind y-slice (last x cell, at canopy-top z).
            k_top = int(np.argmin(np.abs(z_centers - z_canopy_top * 0.8)))
            U_wake = res.U_out[k_top, :, -1] / max(U_infty[k_top], 1e-6)
            deficit = 1.0 - U_wake
            peak = deficit.max()
            if peak > 0.05:
                half = peak / 2
                idxs = np.where(deficit >= half)[0]
                fwhm_cells = idxs.max() - idxs.min() + 1 if len(idxs) else 0
                fwhm_y = fwhm_cells * grid_size / grid_size  # in units of grid_size
            else:
                fwhm_y = 0.0
            print(
                f"  {nu_diff:>8.2f} {grid_size:>5.1f} "
                f"{f'({bx.size - 1},{by.size - 1},{bz.size - 1})':>16s} "
                f"{ratio.mean():>13.3f} {ratio.min():>15.3f} "
                f"{fwhm_y:>16.1f} {best * 1000:>10.2f}"
            )


def run_ab_compare(R: float, n_trees: int, n_gens: int, U_uniform: float) -> None:
    """A vs B: same seeded forest, A under spatially-uniform ``default`` wind,
    B under the ``momentum`` per-branch screened field. Asks "does resolving
    the screening change the trees?" — final structure + cumulative pruning.

    Provisional framing (Step 26a): the two arms do *not* yet feel an
    identical storm sequence (``default`` uses the Fortran long-tail formula;
    the momentum arm a fixed ``U_uniform``). Step 26d will match the storms
    so the only difference is screening."""
    print(
        f"\n=== A/B: uniform default vs momentum per-branch "
        f"(R={R}, n_trees={n_trees}, n_gens={n_gens}, U_uniform={U_uniform}) ==="
    )

    cfg_a = mt.Config(
        tree=mt.TreeConfig(),
        forest=mt.ForestConfig(size=R, n_trees_init=n_trees, n_trees_max=int(n_trees * 1.5)),
        wind=mt.WindConfig(model="default"),
        n_generations=n_gens,
    )
    cfg_b = mt.Config(
        tree=mt.TreeConfig(),
        forest=mt.ForestConfig(size=R, n_trees_init=n_trees, n_trees_max=int(n_trees * 1.5)),
        wind=mt.WindConfig(
            model="momentum",
            grid_size=1.0,
            momentum_pad_x=12.0,
            momentum_pad_y=4.0,
            momentum_pad_z=3.0,
            momentum_U_uniform=U_uniform,
            max_pruning_iterations=8,
            wind_convergence_eps_rel=0.01,
        ),
        n_generations=n_gens,
    )

    def _run(cfg: Any) -> dict[str, Any]:
        forest = mt.Forest(cfg, seed=0)
        total_pruned = 0
        t0 = time.perf_counter()
        for g in range(n_gens):
            total_pruned += forest.step(g).n_pruned_total
        wall = time.perf_counter() - t0
        lr = getattr(forest.wind_fn, "last_result", None)
        return {
            "wall": wall,
            "trees": len(forest.trees),
            "branches": sum(t.get_number_of_branches() for t in forest.trees),
            "pruned": total_pruned,
            "U_min": float(lr.U_branch.min()) if lr is not None else float("nan"),
            "U_mean": float(lr.U_branch.mean()) if lr is not None else float("nan"),
            "U_max": float(lr.U_branch.max()) if lr is not None else float("nan"),
        }

    a = _run(cfg_a)
    b = _run(cfg_b)
    print(f"  {'metric':22s} {'A: uniform default':>20s} {'B: momentum per-branch':>24s}")
    print(f"  {'wall-clock (s)':22s} {a['wall']:>20.2f} {b['wall']:>24.2f}")
    print(
        f"  {'ms/gen mean':22s} {a['wall'] / n_gens * 1000:>20.0f} "
        f"{b['wall'] / n_gens * 1000:>24.0f}"
    )
    print(f"  {'final trees':22s} {a['trees']:>20d} {b['trees']:>24d}")
    print(f"  {'final branches':22s} {a['branches']:>20d} {b['branches']:>24d}")
    print(f"  {'cumulative pruned':22s} {a['pruned']:>20d} {b['pruned']:>24d}")
    print(
        f"  {'B last U_branch m/μ/M':22s} "
        f"{b['U_min']:.3f} / {b['U_mean']:.3f} / {b['U_max']:.3f}  "
        f"(per-branch wind spread the screening resolves)"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--R", type=float, default=100.0, help="forest disk radius")
    ap.add_argument("--n-trees", type=int, default=1000, help="initial tree count")
    ap.add_argument("--n-gens", type=int, default=100, help="generations to run")
    ap.add_argument(
        "--skip-perf", action="store_true", help="skip the forest perf bench (sweep only)"
    )
    ap.add_argument(
        "--skip-sweep", action="store_true", help="skip the (nu_diff, grid_size) sweep (perf only)"
    )
    ap.add_argument(
        "--skip-ab",
        action="store_true",
        help="skip the Step-25c A/B (canopy-mean vs per-branch) compare",
    )
    ap.add_argument(
        "--ab-u-uniform",
        type=float,
        default=2.0,
        help="uniform inflow for the A/B compare (storm strength)",
    )
    args = ap.parse_args()

    if not args.skip_perf:
        run_forest_bench(R=args.R, n_trees=args.n_trees, n_gens=args.n_gens)
    if not args.skip_ab:
        # Smaller scale by default — the A/B runs the forest twice.
        ab_R = min(args.R, 50.0)
        ab_n = min(args.n_trees, 200)
        ab_gens = min(args.n_gens, 60)
        run_ab_compare(R=ab_R, n_trees=ab_n, n_gens=ab_gens, U_uniform=args.ab_u_uniform)
    if not args.skip_sweep:
        # Smaller / shorter warmup for the sweep — it's a parameter
        # exploration, not a paper-scale stress test.
        sweep_R = min(args.R, 50.0)
        sweep_n = min(args.n_trees, 200)
        sweep_gens = min(args.n_gens, 30)
        run_nu_diff_H_sweep(R=sweep_R, n_trees=sweep_n, warmup_gens=sweep_gens)


if __name__ == "__main__":
    main()
