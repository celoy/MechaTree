"""Scale benchmark for a single ``Forest.step`` loop (Step 21b prep).

Times each phase of ``Forest.step`` by monkey-patching the underlying
functions with a tiny wall-clock decorator. Reports per-phase totals so we
can see where the Python port spends its time at island-scale (Eloy et
al., Nat Commun 2017: R = 200 L, N_init = 20 000 random genomes).

Defaults are intentionally modest (n_init=1000, gens=100, size=63) so the
first invocation comes back in ~minutes; pass ``--full`` for the paper-
scale config. The script also estimates the wall-clock for an n_init=20k,
gens=1000 run by extrapolation.

Usage::

    uv run python benchmarks/bench_forest_scale.py
    uv run python benchmarks/bench_forest_scale.py --n-init 2000 --gens 50
    uv run python benchmarks/bench_forest_scale.py --full
    uv run python benchmarks/bench_forest_scale.py --profile

When ``--profile`` is given, cProfile dumps stats to ``out/forest_scale.prof``
which can be inspected with ``snakeviz`` or ``python -m pstats``.
"""

from __future__ import annotations

import argparse
import cProfile
import math
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

import mechatree.forest as forest_mod
import mechatree.growth as growth_mod
import mechatree.light as light_mod
import mechatree.mechanics as mechanics_mod
import mechatree.pruning as pruning_mod
from mechatree.config import Config, ForestConfig, TreeConfig
from mechatree.forest import Forest

# Phases tracked. Each entry is (label, module, attr-name).
_PHASES = [
    ("extract_leaves", light_mod, "extract_leaves"),
    ("intercept", light_mod, "intercept"),
    ("aggregate_onto_trees", light_mod, "aggregate_onto_trees"),
    ("calculate_stresses", mechanics_mod, "calculate_stresses"),
    ("requested_growth", growth_mod, "requested_growth"),
    ("secondary_growth", growth_mod, "secondary_growth"),
    ("prune", pruning_mod, "prune"),
    ("primary_growth", growth_mod, "primary_growth"),
]

# Forest also imports these names directly into its own module — that's the
# version actually called at runtime, so we have to patch there too.
_FOREST_PATCH_TARGETS = [
    "extract_leaves",
    "intercept",
    "aggregate_onto_trees",
    "calculate_stresses",
    "requested_growth",
    "secondary_growth",
    "prune",
    "primary_growth",
]


def _install_timers(totals: dict[str, float]) -> list[tuple[object, str, object]]:
    """Wrap every phase function with a wall-clock decorator and return
    the originals so they can be restored.

    Patches both the canonical module attribute *and* the local rebinding
    inside :mod:`mechatree.forest` (which imported each name with a bare
    ``from … import …``)."""
    saved: list[tuple[object, str, object]] = []
    for label, mod, attr in _PHASES:
        original = getattr(mod, attr)

        def timed(*args, _label=label, _original=original, **kw):
            t0 = time.perf_counter()
            out = _original(*args, **kw)
            totals[_label] += time.perf_counter() - t0
            return out

        saved.append((mod, attr, original))
        setattr(mod, attr, timed)
    # Re-bind inside forest.py too.
    for attr in _FOREST_PATCH_TARGETS:
        original = getattr(forest_mod, attr)
        saved.append((forest_mod, attr, original))
        setattr(forest_mod, attr, getattr(_module_for(attr), attr))
    return saved


def _module_for(attr: str):
    """Map a phase function name back to the module that owns the timed
    version (we just patched it there a moment ago)."""
    if attr in {"extract_leaves", "intercept", "aggregate_onto_trees"}:
        return light_mod
    if attr in {"calculate_stresses"}:
        return mechanics_mod
    if attr in {"requested_growth", "secondary_growth", "primary_growth"}:
        return growth_mod
    if attr in {"prune"}:
        return pruning_mod
    raise KeyError(attr)


def _restore(saved: list[tuple[object, str, object]]) -> None:
    for mod, attr, original in saved:
        setattr(mod, attr, original)


def _bench(cfg: Config, n_generations: int, seed: int) -> dict[str, float]:
    totals: dict[str, float] = defaultdict(float)
    saved = _install_timers(totals)
    try:
        forest = Forest(cfg, seed=seed)
        per_gen_ms: list[float] = []
        n_trees_alive: list[int] = []
        n_branches_alive: list[int] = []
        t_total = time.perf_counter()
        for gen in range(n_generations):
            t_g = time.perf_counter()
            stats = forest.step(gen)
            per_gen_ms.append((time.perf_counter() - t_g) * 1000)
            n_trees_alive.append(stats.n_trees)
            n_branches_alive.append(stats.n_branches_total)
        wall = time.perf_counter() - t_total
    finally:
        _restore(saved)

    print()
    print("=== per-phase wall-clock (totals, % of step time) ===")
    in_step = wall  # we want a denominator that includes the un-instrumented bits
    measured = sum(totals.values())
    for label, _mod, _attr in _PHASES:
        t = totals[label]
        print(f"  {label:>22s}: {t * 1000:>9.1f} ms  ({100 * t / in_step:>5.2f} %)")
    print(f"  {'(measured / total)':>22s}: {measured * 1000:>9.1f} / {wall * 1000:.1f} ms")
    print()
    print("=== per-generation summary ===")
    print(f"  n_generations:       {n_generations}")
    print(f"  total wall:          {wall * 1000:>9.1f} ms ({wall:.2f} s)")
    print(f"  mean ms/gen:         {1000 * wall / n_generations:>9.2f}")
    print(f"  median ms/gen:       {float(np.median(per_gen_ms)):>9.2f}")
    print(f"  p95 ms/gen:          {float(np.percentile(per_gen_ms, 95)):>9.2f}")
    print(f"  max ms/gen:          {max(per_gen_ms):>9.2f}")
    print()
    print(f"  trees @ start / end:    {n_trees_alive[0]} → {n_trees_alive[-1]}")
    print(f"  branches @ start / end: {n_branches_alive[0]} → {n_branches_alive[-1]}")

    return {
        "wall_s": wall,
        "ms_per_gen": 1000 * wall / n_generations,
        "n_trees_final": n_trees_alive[-1],
        "n_branches_final": n_branches_alive[-1],
        **{f"phase_{k}_s": v for k, v in totals.items()},
    }


def _extrapolate(
    result: dict[str, float], target_n_init: int, target_gens: int, source_n_init: int
) -> None:
    """Project the measured per-gen time to the paper-scale config.

    Two scaling assumptions, both deliberately conservative:

    * Per-gen cost grows ~linearly with N_trees (each Python loop in
      :meth:`Forest.step` iterates `for tree in self.trees`, so wall time
      is bounded below by O(N) regardless of what C++ does internally).
    * Light competition `extract_leaves`/`intercept` has been observed to
      grow roughly N²/N_directions for dense stands. We multiply by an
      extra factor of (N_target / N_source) on top of the linear scale
      to flag that concern.
    """
    ratio = target_n_init / source_n_init
    linear_projected_s = result["wall_s"] * ratio * (target_gens / max(1, _gens_used(result)))
    print()
    print(f"=== extrapolation to paper-scale (N_init={target_n_init}, {target_gens} gens) ===")
    print(f"  measured at n_init={source_n_init}: {result['ms_per_gen']:.2f} ms/gen")
    print(
        f"  linear projection:           {linear_projected_s / 60:.1f} min "
        f"({linear_projected_s:.0f} s)"
    )
    print(
        f"  quadratic-light upper bound: "
        f"{linear_projected_s * ratio / 3600:.1f} h ({linear_projected_s * ratio:.0f} s)"
    )


def _gens_used(result: dict[str, float]) -> int:
    return max(1, int(round(result["wall_s"] * 1000 / result["ms_per_gen"])))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-init", type=int, default=1000, help="Initial founders.")
    parser.add_argument(
        "--n-max",
        type=int,
        default=None,
        help="Carrying cap (default: max(n_init, 25k) so births can happen).",
    )
    parser.add_argument(
        "--size",
        type=float,
        default=63.0,
        help="Forest radius in twig-length units. Paper uses 200.",
    )
    parser.add_argument("--gens", type=int, default=100, help="Number of generations.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--full",
        action="store_true",
        help="Paper-scale config (n_init=20000, size=200, gens=1000).",
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Run cProfile and dump stats to out/forest_scale.prof.",
    )
    parser.add_argument("--profile-out", type=Path, default=Path("out/forest_scale.prof"))
    args = parser.parse_args()

    if args.full:
        n_init = 20000
        size = 200.0
        gens = 1000
    else:
        n_init = args.n_init
        size = args.size
        gens = args.gens

    n_max = args.n_max if args.n_max is not None else max(n_init, 25000)

    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(
            size=size,
            n_trees_init=n_init,
            n_trees_max=n_max,
            max_age=1000,
            min_age_for_undersize=5,
            min_branches=11,
        ),
    )

    print("=== Forest-scale benchmark ===")
    print(f"config: N_init={n_init}, N_max={n_max}, R={size:g} L, gens={gens}, seed={args.seed}")
    # Density check: the paper uses 20000 trees on R=200, density ~ 0.16 trees/L².
    area = math.pi * size * size
    print(f"  area = π·R² = {area:.0f} L²  →  density = {n_init / area:.3f} trees / L²")
    print(f"  (paper density: {20000 / (math.pi * 200**2):.3f} trees / L²)")

    if args.profile:
        args.profile_out.parent.mkdir(parents=True, exist_ok=True)
        prof = cProfile.Profile()
        prof.enable()
        result = _bench(cfg, gens, args.seed)
        prof.disable()
        prof.dump_stats(str(args.profile_out))
        print()
        print(f"cProfile stats → {args.profile_out}")
        print("inspect with:  uv run python -m pstats out/forest_scale.prof")
        print("      or:      uv run snakeviz out/forest_scale.prof  (requires snakeviz)")
    else:
        result = _bench(cfg, gens, args.seed)

    if not args.full:
        _extrapolate(result, target_n_init=20000, target_gens=1000, source_n_init=n_init)


if __name__ == "__main__":
    main()
