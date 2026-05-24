"""Branch-merge benchmark: wall time and branch count, merge OFF vs ON.

The growth loop produces single-child parent-child chains every time pruning
deletes one of a branching pair. `PyTree.collapse_single_child_chains` fuses
those chains into one straight segment, preserving the chain's two endpoints
and its total volume. ``collapse_chains_after_prune`` is the targeted
variant that only walks chains seeded by the most recent pruning pass.

Three experiments:

  (A) End-to-end single tree. Run `grow_tree(seed=42)` twice — once without
      the merge, once with `collapse_chains_after_prune` called at the end
      of every generation — and report wall time, branches, and leaves.

  (B) Isolated per-step compute. Grow once (no merge) to ~10 k branches.
      Snapshot per-step cost. Apply one whole-tree collapse pass. Snapshot
      per-step cost again. Measures the pure savings on the same tree.

  (C) Long forest run. Many trees, many generations. The forest has more
      pruning per step (one wind hits all trees), so chains accumulate
      faster. Compares wall time + biomass + tree count between the
      untouched run and the targeted-collapse run.
"""

import math
import time
from dataclasses import replace

import numpy as np

from mechatree.config import Config
from mechatree.forest import Forest
from mechatree.genome import ConstantAllocation, ConstantSafety
from mechatree.growth import primary_growth, requested_growth, secondary_growth
from mechatree.light import Sun, aggregate_onto_trees, extract_leaves, intercept
from mechatree.mechanics import calculate_stresses
from mechatree.pruning import prune
from mechatree.simulate import default_wind_fn, grow_tree


def _run(n_generations: int, *, with_merge: bool, repeats: int = 3):
    """Return (best_seconds, final_branches, final_leaves) over `repeats` runs."""
    best = math.inf
    final_n = 0
    final_leaves = 0
    for _ in range(repeats):
        if with_merge:

            def on_step(_gen, tree):
                tree.collapse_chains_after_prune()
                tree.reorder()
        else:
            on_step = None

        t0 = time.perf_counter()
        tree = grow_tree(Config(), n_generations=n_generations, seed=42, on_step=on_step)
        dt = time.perf_counter() - t0
        best = min(best, dt)
        final_n = tree.get_number_of_branches()
        final_leaves = tree.get_total_leaves()
    return best, final_n, final_leaves


def _experiment_end_to_end() -> None:
    print("(A) End-to-end grow_tree, merge OFF vs ON at every generation:")
    header = (
        f"{'n_gen':>6}  {'n_off':>7}  {'lv_off':>7}  {'n_on':>7}  {'lv_on':>7}  "
        f"{'t_off':>8}  {'t_on':>8}  {'speedup':>8}"
    )
    units = (
        f"{'':6}  {'(brn)':>7}  {'(lvs)':>7}  {'(brn)':>7}  {'(lvs)':>7}  "
        f"{'(s)':>8}  {'(s)':>8}  {'x':>8}"
    )
    print(header)
    print(units)

    # Stopping at 100 generations: that's already ~11k branches (the user's
    # target), and going further runs into rare extreme gusts in the wind's
    # long-tailed amplitude — those wipe the merge-ON tree disproportionately
    # because its merged segments have larger moment arms.
    for n_gen in (25, 50, 100):
        t_off, n_off, lv_off = _run(n_gen, with_merge=False)
        t_on, n_on, lv_on = _run(n_gen, with_merge=True)
        speedup = t_off / t_on if t_on > 0 else float("inf")
        print(
            f"{n_gen:>6}  {n_off:>7d}  {lv_off:>7d}  {n_on:>7d}  {lv_on:>7d}  "
            f"{t_off:>8.3f}  {t_on:>8.3f}  {speedup:>8.2f}"
        )
    print("  note: merging perturbs the simulation — the two columns are not the same tree.")


def _time_one_full_step(tree, cfg, generation: int, rng) -> float:
    """Mirror simulate.py's per-generation pipeline; return seconds."""
    sun = Sun()
    safety = ConstantSafety(cfg.tree.safety if hasattr(cfg.tree, "safety") else 3.0)
    alloc = ConstantAllocation(p_seeds=0.1, p_leaves=0.5, phototropism=0.5)

    t0 = time.perf_counter()
    leaves = extract_leaves([tree], n_directions=sun.n_directions)
    intercept(leaves, sun)
    aggregate_onto_trees(leaves, [tree])
    calculate_stresses(tree, leaf_drag_S0=cfg.tree.leaf_surface, cauchy=cfg.tree.cauchy)
    requested_growth(tree, safety, maintenance_h=cfg.tree.maintenance_h)
    secondary_growth(tree, volume_per_leaf=cfg.tree.volume_per_leaf)
    wind = default_wind_fn(generation, rng)
    prune(tree, wind=wind, leaf_drag_S0=cfg.tree.leaf_surface, cauchy=cfg.tree.cauchy)
    tree.reorder()
    primary_growth(
        tree,
        alloc,
        twig_length=cfg.tree.twig_length,
        twig_diameter=cfg.tree.twig_diameter,
        theta1=cfg.tree.theta1,
        theta2=cfg.tree.theta2,
        gamma1=cfg.tree.gamma1,
        gamma2=cfg.tree.gamma2,
        generation=generation,
    )
    tree.reorder()
    return time.perf_counter() - t0


def _experiment_isolated() -> None:
    """Pure compute comparison: same physical tree, one step before and one
    after a single ``collapse_single_child_chains`` call."""
    print()
    print("(B) Isolated per-step compute on the SAME tree:")
    header = (
        f"{'grew':>6}  {'pre':>8}  {'absorb':>7}  "
        f"{'post':>8}  {'t_pre':>9}  {'t_post':>9}  {'speedup':>8}"
    )
    units = f"{'':6}  {'(brn)':>8}  {'(brn)':>7}  {'(brn)':>8}  {'(s)':>9}  {'(s)':>9}  {'x':>8}"
    print(header)
    print(units)

    for n_grow in (50, 100):
        cfg = Config()
        # Grow once, no merge, to set up the test bed.
        tree = grow_tree(cfg, n_generations=n_grow, seed=42)
        n_pre = tree.get_number_of_branches()

        # Time one more step on the (unchanged) tree.
        rng = np.random.default_rng(99)
        # Best-of-3 to damp jitter (re-run on a freshly grown copy each time).
        best_pre = math.inf
        for _ in range(3):
            t = grow_tree(cfg, n_generations=n_grow, seed=42)
            best_pre = min(best_pre, _time_one_full_step(t, cfg, n_grow, rng))

        # Now collapse and time one more step on the (smaller) tree.
        absorbed = tree.collapse_single_child_chains()
        tree.reorder()
        n_post = tree.get_number_of_branches()

        best_post = math.inf
        for _ in range(3):
            t = grow_tree(cfg, n_generations=n_grow, seed=42)
            t.collapse_single_child_chains()
            t.reorder()
            best_post = min(best_post, _time_one_full_step(t, cfg, n_grow, rng))

        speedup = best_pre / best_post if best_post > 0 else float("inf")
        print(
            f"{n_grow:>6}  {n_pre:>8d}  {absorbed:>7d}  "
            f"{n_post:>8d}  {best_pre:>9.4f}  {best_post:>9.4f}  {speedup:>8.2f}"
        )


def _forest_config(n_trees_init: int, size: float) -> Config:
    """Default config with a smaller forest plot for the benchmark."""
    base = Config()
    return replace(
        base,
        forest=replace(base.forest, n_trees_init=n_trees_init, size=size),
    )


def _run_forest(
    n_trees_init: int,
    size: float,
    n_generations: int,
    *,
    with_merge: bool,
    seed: int = 42,
    repeats: int = 2,
):
    """Return (best_seconds, n_trees_final, n_branches_total, n_leaves_total,
    biomass_total) over `repeats` runs."""
    cfg = _forest_config(n_trees_init, size)

    best = math.inf
    n_trees_final = 0
    n_branches_total = 0
    n_leaves_total = 0
    biomass_total = 0.0
    pi_over_4 = math.pi / 4.0

    for _ in range(repeats):
        forest = Forest(cfg, seed=seed)
        if with_merge:

            def on_step(_gen, _forest, _stats):
                for tree in _forest.trees:
                    tree.collapse_chains_after_prune()
                    tree.reorder()
        else:
            on_step = None

        t0 = time.perf_counter()
        forest.run(n_generations, on_step=on_step)
        dt = time.perf_counter() - t0

        best = min(best, dt)

        # Final snapshot.
        n_trees_final = len(forest.trees)
        n_branches_total = sum(t.get_number_of_branches() for t in forest.trees)
        n_leaves_total = sum(t.get_total_leaves() for t in forest.trees)
        biomass = 0.0
        for tree in forest.trees:
            for i in range(tree.get_number_of_branches()):
                d = tree.get_diameter(i)
                L = tree.get_length(i)
                biomass += pi_over_4 * d * d * L
        biomass_total = biomass

    return best, n_trees_final, n_branches_total, n_leaves_total, biomass_total


def _experiment_forest() -> None:
    print()
    print("(C) Forest run, targeted collapse OFF vs ON every generation:")
    header = (
        f"{'n_gen':>5}  {'trees':>5}  "
        f"{'br_off':>8}  {'br_on':>8}  {'lv_off':>7}  {'lv_on':>7}  "
        f"{'V_off':>8}  {'V_on':>8}  {'t_off':>8}  {'t_on':>8}  {'speedup':>8}"
    )
    units = (
        f"{'':5}  {'':>5}  "
        f"{'(brn)':>8}  {'(brn)':>8}  {'(lvs)':>7}  {'(lvs)':>7}  "
        f"{'(vol)':>8}  {'(vol)':>8}  {'(s)':>8}  {'(s)':>8}  {'x':>8}"
    )
    print(header)
    print(units)

    # Modest forest so we can afford long runs.
    n_trees_init = 15
    size = 15.0
    for n_gen in (50, 100):
        t_off, n_t_off, br_off, lv_off, V_off = _run_forest(
            n_trees_init, size, n_gen, with_merge=False
        )
        t_on, n_t_on, br_on, lv_on, V_on = _run_forest(n_trees_init, size, n_gen, with_merge=True)
        speedup = t_off / t_on if t_on > 0 else float("inf")
        # n_trees_off / n_trees_on are usually equal for matched seeds since
        # death depends on per-tree state; we print the OFF number.
        print(
            f"{n_gen:>5}  {n_t_off:>5d}  "
            f"{br_off:>8d}  {br_on:>8d}  {lv_off:>7d}  {lv_on:>7d}  "
            f"{V_off:>8.1f}  {V_on:>8.1f}  {t_off:>8.2f}  {t_on:>8.2f}  {speedup:>8.2f}"
        )
    print(
        f"  config: n_trees_init={n_trees_init}, plot size={size}, seed=42; "
        "merging perturbs the simulation slightly."
    )


def main() -> None:
    _experiment_end_to_end()
    _experiment_isolated()
    _experiment_forest()


if __name__ == "__main__":
    main()
