"""Phase-0 sizing benchmark for Step 24 (coupled wind ↔ pruning fixed-point loop).

Measures the per-iteration cost of the inner loop's two ingredients
(``forest_to_cylinders`` + ``BulkThinningBranchWindModel.compute``) so we
know whether the forest-wide fixed-point design is viable at scale before
wiring it into the orchestrator. Per the user's directive: "may be super
slow... so benchmarking and test are necessary on small forests of ~10
trees before implementation".

For each forest scale we:

1. Grow a ``Forest`` to a representative steady state under the
   DendroFlow bridge.
2. Time one ``forest_to_cylinders`` call and one
   ``BulkThinningBranchWindModel.compute`` call in isolation, repeated
   ``--repeats`` times (default 20) for a stable mean.
3. Time ``Forest.step`` itself at the same point so we can express the
   inner-iteration cost as a percentage of the per-step cost.
4. Project per-step overhead at 1 / 2 / 4 / 8 inner iterations.
5. Run an ε-tolerance probe: drop ~0.5 % of branches from random trees
   (a "sparse storm" iter-1 outcome), measure the relative change in
   pooled canopy mean, and report whether ε_rel = 0.01 would
   short-circuit iteration 2.

Usage::

    uv run python benchmarks/bench_wind_fixed_point.py
    uv run python benchmarks/bench_wind_fixed_point.py --full       # adds 20k-tree scale
    uv run python benchmarks/bench_wind_fixed_point.py --scales 10 100

Requires the ``dendroflow`` extra (see CLAUDE.md Step 17).
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass

import numpy as np

from mechatree.config import Config, ForestConfig, TreeConfig, WindConfig
from mechatree.forest import Forest
from mechatree.wind.dendroflow import (
    DendroFlowWindParams,
    forest_to_cylinders,
    make_dendroflow_wind_fn,
)

# Default wind profile — calmer than examples/dendroflow_wind.yaml so the
# small-scale forests don't go extinct in the warmup. Mean ~1.0, which is
# in the same ballpark as default_wind_fn's typical amplitude (~0.835).
# Module-level so the CLI ``--u-scale`` flag can rescale before any
# Forest is built.
U_INFTY = np.array([0.5, 0.7, 0.9, 1.05, 1.15, 1.25, 1.3, 1.35, 1.4, 1.45])
Z_CENTERS = np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75])

DEFAULT_SCALES = [
    (10, 30, 30.0),  # mini
    (100, 80, 50.0),  # small
    (1000, 120, 80.0),  # medium
]
FULL_EXTRA = [(20000, 60, 200.0)]


@dataclass
class ScaleResult:
    n_init: int
    n_trees_alive: int
    n_branches_total: int
    ms_per_step_baseline: float
    ms_forest_to_cylinders: float
    ms_compute: float
    eps_relative_delta: float | None  # None if no trees / no branches to drop


def _build_forest(n_init: int, size: float, gens: int, seed: int) -> Forest:
    """Grow a Forest under DendroFlow wind for ``gens`` generations."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(
            size=size,
            n_trees_init=n_init,
            n_trees_max=max(n_init * 2, 25_000),
            max_age=1000,
            min_age_for_undersize=5,
            min_branches=11,
        ),
        wind=WindConfig(
            model="dendroflow",
            U_infty=tuple(U_INFTY.tolist()),
            z_centers=tuple(Z_CENTERS.tolist()),
            H=0.5,
            C_D=1.0,
            z_representative="mean",
        ),
    )
    forest = Forest(cfg, seed=seed)
    for gen in range(gens):
        forest.step(gen)
    return forest


def _time_one_step(forest: Forest, gen: int, repeats: int = 5) -> float:
    """Time a single ``Forest.step`` (averaged over ``repeats`` runs).

    Each repeat advances the forest by one generation, so the average
    cost is meaningful only as an order-of-magnitude estimate at the
    current size. Good enough for the decision gate.
    """
    ts: list[float] = []
    for k in range(repeats):
        t0 = time.perf_counter()
        forest.step(gen + k)
        ts.append(time.perf_counter() - t0)
    return float(np.mean(ts)) * 1000.0


def _time_inner_iter(forest: Forest, repeats: int = 20) -> tuple[float, float]:
    """Time the two ingredients of the inner loop separately.

    Returns ``(ms_pool, ms_compute)``.
    """
    params = DendroFlowWindParams(
        U_infty=U_INFTY,
        z_centers=Z_CENTERS,
        H=0.5,
        C_D=1.0,
    )
    bridge = make_dendroflow_wind_fn(
        U_infty=U_INFTY,
        z_centers=Z_CENTERS,
        H=0.5,
        C_D=1.0,
        z_representative="mean",
    )
    # Warm up (first call allocates).
    cyl = forest_to_cylinders(forest.trees)
    _ = bridge._model.compute(cyl, wind_params=params.to_namespace())

    pool_times: list[float] = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        cyl = forest_to_cylinders(forest.trees)
        pool_times.append(time.perf_counter() - t0)

    compute_times: list[float] = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        _ = bridge._model.compute(cyl, wind_params=params.to_namespace())
        compute_times.append(time.perf_counter() - t0)

    return float(np.mean(pool_times)) * 1000.0, float(np.mean(compute_times)) * 1000.0


def _build_arrays(trees):
    """Rebuild the (start, axis, D, L, tree_id) arrays that
    ``forest_to_cylinders`` would feed to DendroFlow. Lets us slice rows
    out of the canopy without going through the DataFrame round-trip.
    """
    sizes = [t.get_number_of_branches() for t in trees]
    total = sum(sizes)
    start = np.empty((total, 3), dtype=float)
    axis = np.empty((total, 3), dtype=float)
    D = np.empty(total, dtype=float)
    L = np.empty(total, dtype=float)
    tree_id = np.empty(total, dtype=float)
    off = 0
    for ti, (tree, n) in enumerate(zip(trees, sizes, strict=True)):
        for i in range(n):
            start[off + i] = tree.get_location(i)
            axis[off + i] = tree.get_unit_t(i)
            D[off + i] = tree.get_diameter(i)
            L[off + i] = tree.get_length(i)
        tree_id[off : off + n] = float(ti)
        off += n
    return start, axis, D, L, tree_id


def _eps_probe(forest: Forest, drop_fraction: float, rng: np.random.Generator) -> float | None:
    """Drop ``drop_fraction`` of branches at random across the canopy, then
    measure the relative change in the pooled canopy-mean wind. Returns
    ``|wind_after - wind_before| / |wind_before|`` (horizontal magnitude),
    or ``None`` if the probe is degenerate (no trees, no branches).

    Models the "sparse storm" case: iteration 1 cut ~0.5 % of the forest's
    branches; iteration 2 would now re-pool + recompute wind — by how much
    does that change the canopy mean, and would ``eps_rel = 0.01`` skip it?
    """
    if not forest.trees:
        return None
    from dendroflow import from_arrays as _df_from_arrays

    params = DendroFlowWindParams(
        U_infty=U_INFTY,
        z_centers=Z_CENTERS,
        H=0.5,
        C_D=1.0,
    )
    bridge = make_dendroflow_wind_fn(
        U_infty=U_INFTY,
        z_centers=Z_CENTERS,
        H=0.5,
        C_D=1.0,
        z_representative="mean",
    )

    start, axis, D, L, tree_id = _build_arrays(forest.trees)
    total = start.shape[0]
    if total == 0:
        return None
    n_to_drop = max(1, int(round(drop_fraction * total)))
    if n_to_drop >= total:
        return None

    cyl_before = _df_from_arrays(start=start, axis=axis, D=D, L=L, tree_id=tree_id)
    wind_before = bridge._model.compute(cyl_before, wind_params=params.to_namespace()).canopy_mean

    drop_idxs = rng.choice(total, size=n_to_drop, replace=False)
    keep = np.ones(total, dtype=bool)
    keep[drop_idxs] = False
    cyl_after = _df_from_arrays(
        start=start[keep], axis=axis[keep], D=D[keep], L=L[keep], tree_id=tree_id[keep]
    )
    wind_after = bridge._model.compute(cyl_after, wind_params=params.to_namespace()).canopy_mean

    dwx = wind_after[0] - wind_before[0]
    dwy = wind_after[1] - wind_before[1]
    delta = math.hypot(dwx, dwy)
    ref = max(math.hypot(wind_before[0], wind_before[1]), 1e-6)
    return delta / ref


def _bench_scale(n_init: int, size: float, gens: int, repeats: int, seed: int) -> ScaleResult:
    print()
    print(f"=== scale: n_init={n_init}, size={size:g} L, grown {gens} gens ===")
    t0 = time.perf_counter()
    forest = _build_forest(n_init=n_init, size=size, gens=gens, seed=seed)
    grow_s = time.perf_counter() - t0
    n_alive = len(forest.trees)
    n_branches = sum(t.get_number_of_branches() for t in forest.trees)
    print(f"  warmup: {grow_s:.2f} s   trees alive: {n_alive}   branches total: {n_branches}")

    if n_alive == 0:
        print("  forest went extinct under DendroFlow wind — skipping inner-loop measurement")
        return ScaleResult(
            n_init=n_init,
            n_trees_alive=0,
            n_branches_total=0,
            ms_per_step_baseline=float("nan"),
            ms_forest_to_cylinders=float("nan"),
            ms_compute=float("nan"),
            eps_relative_delta=None,
        )

    ms_step = _time_one_step(forest, gen=gens, repeats=5)
    ms_pool, ms_compute = _time_inner_iter(forest, repeats=repeats)

    rng = np.random.default_rng(seed + 1)
    eps_delta = _eps_probe(forest, drop_fraction=0.005, rng=rng)

    return ScaleResult(
        n_init=n_init,
        n_trees_alive=n_alive,
        n_branches_total=n_branches,
        ms_per_step_baseline=ms_step,
        ms_forest_to_cylinders=ms_pool,
        ms_compute=ms_compute,
        eps_relative_delta=eps_delta,
    )


def _print_summary(results: list[ScaleResult], eps_rel: float) -> None:
    print()
    print("=== summary ===")
    print(
        f"  {'n_init':>7s} {'alive':>6s} {'branches':>9s} "
        f"{'step ms':>9s} {'pool ms':>9s} {'compute ms':>11s} {'inner ms':>9s} "
        f"{'1×':>6s} {'2×':>6s} {'4×':>6s} {'8×':>6s} "
        f"{'εδ@0.5%':>9s} {'skip?':>6s}"
    )
    for r in results:
        if not math.isfinite(r.ms_per_step_baseline):
            print(f"  {r.n_init:>7d} {'extinct':>6s}")
            continue
        inner = r.ms_forest_to_cylinders + r.ms_compute
        baseline = r.ms_per_step_baseline
        if baseline > 0:
            pct1 = 100 * inner / baseline
            pct2 = 100 * 2 * inner / baseline
            pct4 = 100 * 4 * inner / baseline
            pct8 = 100 * 8 * inner / baseline
        else:
            pct1 = pct2 = pct4 = pct8 = float("nan")
        if r.eps_relative_delta is None:
            eps_str = "—"
            skip = "—"
        else:
            eps_str = f"{r.eps_relative_delta * 100:.3f}%"
            skip = "yes" if r.eps_relative_delta < eps_rel else "no"
        print(
            f"  {r.n_init:>7d} {r.n_trees_alive:>6d} {r.n_branches_total:>9d} "
            f"{baseline:>9.2f} {r.ms_forest_to_cylinders:>9.3f} {r.ms_compute:>11.3f} "
            f"{inner:>9.3f} "
            f"{pct1:>5.1f}% {pct2:>5.1f}% {pct4:>5.1f}% {pct8:>5.1f}% "
            f"{eps_str:>9s} {skip:>6s}"
        )
    print()
    print(f"  ε early-exit threshold used: ε_rel = {eps_rel:g}")
    print("  εδ@0.5% column: pooled |Δwind|/|wind| after dropping 0.5% of branches.")
    print("                  If εδ < ε_rel, iteration 2 of a sparse storm would short-circuit.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scales",
        nargs="+",
        type=int,
        default=None,
        help="Override the default scale list (e.g. --scales 10 100).",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Add the 20k-tree island scale (slow; for the record).",
    )
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--eps-rel", type=float, default=0.01)
    parser.add_argument(
        "--u-scale",
        type=float,
        default=1.0,
        help=(
            "Scale the default U_infty profile (mean ~1.0) by this factor. "
            "u-scale=3 reproduces the strong wind from "
            "examples/dendroflow_wind.yaml."
        ),
    )
    args = parser.parse_args()

    # Rescale the global U_infty in-place so the helper functions see it.
    global U_INFTY  # noqa: PLW0603
    U_INFTY = U_INFTY * args.u_scale

    scales = list(DEFAULT_SCALES)
    if args.full:
        scales.extend(FULL_EXTRA)
    if args.scales is not None:
        # Use the user-specified n_init values; pair them with sane defaults
        # for (gens, size) by matching against the default table when possible.
        defaults_by_n = {n: (g, s) for (n, g, s) in DEFAULT_SCALES + FULL_EXTRA}
        scales = []
        for n in args.scales:
            g, s = defaults_by_n.get(n, (80, max(30.0, math.sqrt(n))))
            scales.append((n, g, s))

    print("=== Phase 0: forest-wide wind-prune sizing benchmark ===")
    print(
        f"  default U_infty: shape {U_INFTY.shape}, range "
        f"[{U_INFTY.min():.2f}, {U_INFTY.max():.2f}]"
    )
    print(f"  repeats: {args.repeats}   seed: {args.seed}")

    results: list[ScaleResult] = []
    for n_init, gens, size in scales:
        results.append(
            _bench_scale(n_init=n_init, size=size, gens=gens, repeats=args.repeats, seed=args.seed)
        )

    _print_summary(results, eps_rel=args.eps_rel)


if __name__ == "__main__":
    main()
