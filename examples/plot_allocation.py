#!/usr/bin/env python3
"""Drive a single tree and plot per-step allocation diagnostics.

CLI demo for :func:`mechatree.plotting.plot_allocation`. The reusable
plot helper lives in the library; this script is a thin wrapper that
collects :class:`mechatree.simulate.TreeStats` history via the
``on_step`` callback and feeds it to the plot.

Python port of ``../Eloy2017_NatComm_archive/plot_allocation_vs_t.m``.

Run with::

    uv run python examples/plot_allocation.py --iterations 200 --seed 42
"""

from __future__ import annotations

import argparse
from pathlib import Path

import mechatree as mt


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
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    cfg = mt.load_config(args.config)

    history = []  # list of TreeStats

    def on_step(_gen, _tree, stats):
        history.append(stats)

    mt.grow_tree(cfg, n_generations=args.iterations, seed=args.seed, on_step=on_step)

    print(f"{'gen':>4}  {'branches':>8}  {'wind':>6}  {'seeds':>5}  {'pruned':>6}  {'reserve':>8}")
    for stats in history:
        if stats.generation % max(1, args.iterations // 10) == 0:
            print(
                f"{stats.generation:>4}  {stats.n_branches:>8}  "
                f"{stats.wind_amplitude:>6.2f}  {stats.n_seeds:>5}  "
                f"{stats.n_pruned:>6}  "
                f"{stats.reserve / cfg.tree.volume_twig:>8.2f}"
            )

    if not args.no_show:
        fig = mt.plot_allocation(history, volume_twig=cfg.tree.volume_twig)
        fig.show()


if __name__ == "__main__":
    main()
