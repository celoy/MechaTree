#!/usr/bin/env python3
"""Plug in your own growth law / wind / sun.

Demonstrates the customization surface around
:func:`mechatree.simulate.grow_tree`:

1. **Safety / allocation** — pass non-default ``ConstantSafety`` and
   ``ConstantAllocation`` to bias the tree toward survival vs reproduction.
2. **Wind function** — a Python callable ``(generation, rng) -> (x, y, z)``.
   Override the Fortran storm with steady wind from the west, or a sudden
   storm at a chosen generation.
3. **Sun model** — ``Sun.from_arrays(elev, azim)`` for a single overhead
   direction, or an arbitrary integration scheme.

A non-constant safety / allocation model is **not** yet pluggable from
Python alone — the C++ side dispatches through a virtual ``compute()`` and
needs a concrete subclass. Adding a NeuralSafety subclass is a self-
contained later step; until then, your customization knob is the constant
value plus the wind / sun functions above.

Run with::

    uv run python examples/custom_simulation.py
    uv run python examples/custom_simulation.py --no-show
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from mechatree.config import load_config
from mechatree.genome import ConstantAllocation, ConstantSafety
from mechatree.light import Sun
from mechatree.plotting import plot_tree_3d
from mechatree.simulate import grow_tree


def steady_west_wind(generation: int, rng: np.random.Generator) -> tuple[float, float, float]:
    """Always (1, 0, 0) — no gusts, no rotation. A "trade-wind" stand-in."""
    return (1.0, 0.0, 0.0)


def storm_at_gen_50(generation: int, rng: np.random.Generator) -> tuple[float, float, float]:
    """Calm for the first 50 generations, then a fixed storm from the south."""
    if generation < 50:
        return (0.0, 0.0, 0.0)
    return (0.0, 3.0, 0.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "forest.yaml",
        help="YAML config file.",
    )
    parser.add_argument("--iterations", type=int, default=80, help="Number of generations.")
    parser.add_argument("--seed", type=int, default=42, help="Master RNG seed.")
    parser.add_argument("--no-show", action="store_true", help="Skip opening the plotly figure.")
    args = parser.parse_args()

    cfg = load_config(args.config)
    print(f"Running 3 simulations side-by-side, {args.iterations} gens each.\n")

    # 1. Default — Fortran storm wind, default sun, default genome.
    t_default = grow_tree(cfg, n_generations=args.iterations, seed=args.seed)

    # 2. Steady west wind + overhead sun + "save your reserves" genome.
    overhead_sun = Sun.from_arrays(elev=[0.1, 0.4], azim=[0.0, math.pi])
    cautious_safety = ConstantSafety(2.0)  # more growth headroom per stress unit
    save_reserves = ConstantAllocation(p_seeds=0.0, p_leaves=0.3, phototropism=0.5)
    t_steady = grow_tree(
        cfg,
        n_generations=args.iterations,
        seed=args.seed,
        safety=cautious_safety,
        allocation=save_reserves,
        wind_fn=steady_west_wind,
        sun=overhead_sun,
    )

    # 3. Calm-then-storm — same default genome, but a delayed wind event.
    t_late_storm = grow_tree(
        cfg, n_generations=args.iterations, seed=args.seed, wind_fn=storm_at_gen_50
    )

    print(f"{'scenario':>20}  {'branches':>10}  {'leaves':>8}  {'trunk d':>8}")
    for label, tree in [
        ("default", t_default),
        ("steady-west + overhead", t_steady),
        ("calm-then-storm", t_late_storm),
    ]:
        print(
            f"{label:>20}  "
            f"{tree.get_number_of_branches():>10}  "
            f"{tree.get_total_leaves():>8}  "
            f"{tree.get_diameter(0):>8.4f}"
        )

    if not args.no_show:
        from plotly.subplots import make_subplots

        scenarios = [
            ("default storm", t_default),
            ("steady-west + overhead", t_steady),
            ("calm-then-storm", t_late_storm),
        ]
        fig = make_subplots(
            rows=1,
            cols=3,
            specs=[[{"type": "scene"} for _ in scenarios]],
            subplot_titles=[label for label, _ in scenarios],
        )
        for col, (_, tree) in enumerate(scenarios, start=1):
            for trace in plot_tree_3d(tree).data:
                fig.add_trace(trace, row=1, col=col)
        fig.update_layout(width=1500, height=550, paper_bgcolor="white")
        fig.show()


if __name__ == "__main__":
    main()
