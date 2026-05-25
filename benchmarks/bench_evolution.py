"""Per-generation timings for an evolving Forest (Step 21).

Drives a small Darwinian-island tournament (n_init=20, max=100, 50 gens)
and reports the wall-clock per generation plus a coarse phase breakdown:

* total per-gen time
* per-tree mutation overhead (only the ``Genome.mutate`` calls during
  ``_grow_and_disperse``)
* time spent in ``Genome.to_models`` materialising fresh C++ NN pairs

Mutation overhead should be a negligible fraction of total (a sanity
check that evolution adds ~no cost on top of the Step-12 Forest).

Not a CI gate; results go into ``baseline.md`` for tracking.
"""

from __future__ import annotations

import time

from mechatree.config import Config, ForestConfig, TreeConfig
from mechatree.evolution import Genome, run_tournament


def main() -> None:
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(
            size=30.0,
            n_trees_init=20,
            n_trees_max=100,
            max_age=1000,
            min_age_for_undersize=5,
            min_branches=11,
        ),
    )
    n_gens = 50

    # Warm up so the first generation's JIT / numpy import overhead doesn't
    # poison the average.
    _ = run_tournament(cfg, n_generations=2, seed=99)

    # Instrument Genome.mutate + to_models so we can measure their share.
    mutate_total = 0.0
    to_models_total = 0.0
    orig_mutate = Genome.mutate
    orig_to_models = Genome.to_models

    def timed_mutate(self, rng, **kw):  # type: ignore[no-untyped-def]
        nonlocal mutate_total
        t0 = time.perf_counter()
        out = orig_mutate(self, rng, **kw)
        mutate_total += time.perf_counter() - t0
        return out

    def timed_to_models(self):  # type: ignore[no-untyped-def]
        nonlocal to_models_total
        t0 = time.perf_counter()
        out = orig_to_models(self)
        to_models_total += time.perf_counter() - t0
        return out

    Genome.mutate = timed_mutate  # type: ignore[method-assign]
    Genome.to_models = timed_to_models  # type: ignore[method-assign]
    try:
        t_start = time.perf_counter()
        result = run_tournament(cfg, n_generations=n_gens, seed=0)
        total = time.perf_counter() - t_start
    finally:
        Genome.mutate = orig_mutate  # type: ignore[method-assign]
        Genome.to_models = orig_to_models  # type: ignore[method-assign]

    n_branches_total = sum(s.n_branches_total for s in result.history) / n_gens
    n_trees_final = result.history[-1].n_trees
    n_born_total = sum(s.n_born for s in result.history)
    n_died_total = sum(s.n_died for s in result.history)
    n_lineages_alive = result.history[-1].n_lineages_alive

    print("=== mechatree.evolution micro-benchmark ===")
    print(f"config: n_init=20, n_max=100, size=30, {n_gens} gens, seed=0")
    print(
        f"final: {n_trees_final} trees, {n_lineages_alive} lineages alive, "
        f"{n_born_total} births, {n_died_total} deaths"
    )
    print(f"mean n_branches over the run: {n_branches_total:.0f}")
    print()
    print(f"total wall time:    {total * 1000:>8.1f} ms")
    print(f"per generation:     {total * 1000 / n_gens:>8.1f} ms/gen")
    print(
        f"mutate total:       {mutate_total * 1000:>8.1f} ms ({100 * mutate_total / total:>5.2f} %)"
    )
    print(
        f"to_models total:    {to_models_total * 1000:>8.1f} ms "
        f"({100 * to_models_total / total:>5.2f} %)"
    )
    print(
        f"evolution overhead: {(mutate_total + to_models_total) * 1000:>8.1f} ms "
        f"({100 * (mutate_total + to_models_total) / total:>5.2f} %)"
    )

    # Quick sanity probe: rerun the same config without evolution and
    # report the relative cost.
    from mechatree.forest import Forest

    f = Forest(cfg, seed=0)
    t0 = time.perf_counter()
    for gen in range(n_gens):
        f.step(gen)
    baseline = time.perf_counter() - t0
    print()
    print(f"shared-model Forest baseline (no genomes): {baseline * 1000:>8.1f} ms total")
    overhead = (total - baseline) / baseline * 100.0
    print(f"evolution overhead vs Step-12 baseline:    {overhead:>+8.1f} %")


if __name__ == "__main__":
    main()
