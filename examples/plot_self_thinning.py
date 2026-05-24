#!/usr/bin/env python3
"""Drive a forest and plot self-thinning (N vs M).

CLI demo for :func:`mechatree.plotting.plot_self_thinning`. The reusable
plot helper lives in the library; this script is a thin wrapper that
simulates a forest and feeds its per-step stats history to the plot.

Python port of ``../Eloy2017_NatComm_archive/self_thinning.m``.

Run with::

    uv run python examples/plot_self_thinning.py --iterations 200 --seed 42
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

from mechatree.config import load_config
from mechatree.forest import Forest
from mechatree.plotting import plot_self_thinning


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "forest.yaml",
        help="YAML config file.",
    )
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-trees-init", type=int, default=80)
    parser.add_argument("--n-trees-max", type=int, default=400)
    parser.add_argument("--size", type=float, default=30.0)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    cfg = replace(
        cfg,
        forest=replace(
            cfg.forest,
            n_trees_init=args.n_trees_init,
            n_trees_max=args.n_trees_max,
            size=args.size,
        ),
    )

    history = []  # list of ForestStats

    def on_step(_gen, _forest, stats):
        history.append(stats)

    Forest(cfg, seed=args.seed).run(args.iterations, on_step=on_step)

    print(f"{'gen':>6}  {'N (trees)':>10}  {'M (biomass)':>12}")
    for stats in history:
        if stats.generation % max(1, args.iterations // 10) == 0:
            print(f"{stats.generation:>6}  {stats.n_trees:>10}  {stats.biomass_total:>12.3f}")

    if not args.no_show:
        fig = plot_self_thinning(history)
        fig.show()


if __name__ == "__main__":
    main()
