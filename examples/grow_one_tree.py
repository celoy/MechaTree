#!/usr/bin/env python3
"""Grow one tree from a YAML config.

Loads ``examples/forest.yaml``, runs :func:`mechatree.simulate.grow_tree`
for N generations, prints a summary, and renders the result in 3D.

Run with::

    uv run python examples/grow_one_tree.py
    uv run python examples/grow_one_tree.py --iterations 200 --seed 42

This is the user-facing equivalent of ``legacy_fortran/tree.f90`` — a
single tree growing under wind and light.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from mechatree.config import load_config
from mechatree.plotting import plot_tree_3d
from mechatree.simulate import grow_tree


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "forest.yaml",
        help="YAML config file (defaults to examples/forest.yaml).",
    )
    parser.add_argument("--iterations", type=int, default=100, help="Number of generations.")
    parser.add_argument("--seed", type=int, default=42, help="Master RNG seed.")
    parser.add_argument("--no-show", action="store_true", help="Skip opening the plotly figure.")
    args = parser.parse_args()

    cfg = load_config(args.config)

    history = []

    def on_step(gen, tree):
        if gen % max(1, args.iterations // 10) == 0 or gen == args.iterations - 1:
            history.append(
                (gen, tree.get_number_of_branches(), tree.get_total_leaves(), tree.get_reserve())
            )

    tree = grow_tree(cfg, n_generations=args.iterations, seed=args.seed, on_step=on_step)

    print(f"Ran {args.iterations} generations, seed={args.seed}")
    print(f"  twig dimensions:  L={cfg.tree.twig_length}, d={cfg.tree.twig_diameter}")
    print(f"  Cauchy stiffness: {cfg.tree.cauchy}")
    print()
    print(f"{'gen':>6}  {'branches':>10}  {'leaves':>8}  {'reserve':>10}")
    for gen, nb, nl, r in history:
        print(f"{gen:>6}  {nb:>10}  {nl:>8}  {r:>10.4f}")
    print()
    print(
        f"Final tree: {tree.get_number_of_branches()} branches, "
        f"{tree.get_total_leaves()} leaves, "
        f"trunk diameter = {tree.get_diameter(0):.4f}"
    )

    if not args.no_show:
        fig = plot_tree_3d(tree)
        fig.show()


if __name__ == "__main__":
    main()
