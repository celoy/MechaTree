#!/usr/bin/env python3
"""Drive a tree to maturity and plot Strahler-order diagnostics.

CLI demo for :func:`mechatree.plotting.plot_strahler_diagnostics`. The
reusable plot helper lives in the library; this script is a thin wrapper
that grows a tree and feeds it to the plot.

Python port of the MATLAB scripts in ``legacy/matlab/``:
``plot_stat_single_tree.m`` + ``Fractal_dim.m`` +
``plot_area_preservation_1tree.m``.

Run with::

    uv run python examples/plot_strahler.py --iterations 300 --seed 42
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
    parser.add_argument(
        "--iterations",
        type=int,
        default=200,
        help="Number of generations (more = bigger tree, better statistics).",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    cfg = mt.load_config(args.config)
    tree = mt.grow_tree(cfg, n_generations=args.iterations, seed=args.seed)

    summary = mt.strahler_summary(tree)
    ratios = mt.leonardo_ratios(tree)
    T = mt.tokunaga_matrix(tree)

    print(
        f"Tree: {tree.get_number_of_branches()} branches, Strahler max order = {summary.max_order}"
    )
    print()
    print(f"{'order':>6}  {'n':>6}  {'<L>':>8}  {'<d>':>8}  {'<A>':>10}")
    for k in range(summary.max_order):
        print(
            f"{k + 1:>6}  {summary.n_branches[k]:>6}  "
            f"{summary.mean_length[k]:>8.3f}  "
            f"{summary.mean_diameter[k]:>8.4f}  "
            f"{summary.mean_area[k]:>10.5f}"
        )
    print()
    if ratios.size:
        print(f"Leonardo ratio at {ratios.size} junctions:")
        print(
            f"  mean = {ratios.mean():.3f}    "
            f"std = {ratios.std():.3f}    "
            f"median = {float(sorted(ratios)[ratios.size // 2]):.3f}"
        )
    print()
    print(f"Tokunaga matrix ({T.shape[0]}x{T.shape[1]}):")
    print(T)

    if not args.no_show:
        fig = mt.plot_strahler_diagnostics(tree)
        fig.show()


if __name__ == "__main__":
    main()
