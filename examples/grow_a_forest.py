#!/usr/bin/env python3
"""Build a forest and watch it self-prune under wind.

Loads ``examples/forest.yaml``, drives a :class:`mechatree.forest.Forest`
for N generations, plots biomass / tree count over time, and renders a
top-down view of the final stand.

Run with::

    uv run python examples/grow_a_forest.py
    uv run python examples/grow_a_forest.py --iterations 200 --seed 42 \\
        --n-trees-init 30 --size 50

This is the user-facing equivalent of ``legacy_fortran/Forest.f90`` — a
crowded plot of trees competing for light, with the storm wind culling
the over-extended.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

from mechatree.config import load_config
from mechatree.forest import Forest
from mechatree.plotting import plot_forest_topdown


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
    parser.add_argument("--n-trees-init", type=int, default=20, help="Override initial tree count.")
    parser.add_argument("--size", type=float, default=30.0, help="Override forest radius.")
    parser.add_argument("--n-trees-max", type=int, default=300, help="Override population cap.")
    parser.add_argument("--no-show", action="store_true", help="Skip opening the plotly figures.")
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

    history = []  # (gen, n_trees, biomass, n_branches)
    forest = Forest(cfg, seed=args.seed)
    print(f"Initial population: {len(forest.trees)} trees on a disk of radius {args.size}")

    def on_step(gen, _forest, stats):
        history.append((gen, stats.n_trees, stats.biomass_total, stats.n_branches_total))
        if gen % max(1, args.iterations // 10) == 0:
            print(
                f"  gen={gen:>4}: trees={stats.n_trees:>4}  "
                f"biomass={stats.biomass_total:>8.2f}  "
                f"branches={stats.n_branches_total:>6}  "
                f"+born={stats.n_born:<3} -died={stats.n_died:<3}"
            )

    forest.run(args.iterations, on_step=on_step)

    final = history[-1] if history else (0, 0, 0.0, 0)
    print()
    print(f"Final: trees={final[1]}, biomass={final[2]:.2f}, total branches={final[3]}")

    if not args.no_show:
        import plotly.graph_objects as go

        # Population & biomass over time (twin y-axis via plotly's yaxis2).
        # The top-down view is opened separately because its equal-scale
        # geometry doesn't compose cleanly into a subplot grid.
        gens = [h[0] for h in history]
        trees_over_time = [h[1] for h in history]
        biomass_over_time = [h[2] for h in history]

        fig_history = go.Figure()
        fig_history.add_trace(
            go.Scatter(
                x=gens,
                y=trees_over_time,
                name="trees alive",
                line=dict(color="forestgreen"),  # noqa: C408
                yaxis="y",
            )
        )
        fig_history.add_trace(
            go.Scatter(
                x=gens,
                y=biomass_over_time,
                name="biomass",
                line=dict(color="saddlebrown"),  # noqa: C408
                yaxis="y2",
            )
        )
        fig_history.update_layout(
            title="Population & biomass over time",
            xaxis_title="generation",
            yaxis=dict(  # noqa: C408
                title=dict(text="trees alive", font=dict(color="forestgreen")),  # noqa: C408
                tickfont=dict(color="forestgreen"),  # noqa: C408
            ),
            yaxis2=dict(  # noqa: C408
                title=dict(text="total biomass", font=dict(color="saddlebrown")),  # noqa: C408
                tickfont=dict(color="saddlebrown"),  # noqa: C408
                overlaying="y",
                side="right",
            ),
            paper_bgcolor="white",
            plot_bgcolor="white",
            width=900,
            height=500,
        )
        fig_history.show()

        fig_top = plot_forest_topdown(forest)
        fig_top.update_layout(title=f"Final stand (n={len(forest.trees)})")
        fig_top.show()


if __name__ == "__main__":
    main()
