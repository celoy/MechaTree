#!/usr/bin/env python3
"""Run a Darwinian-island tournament on a Forest (Step 21).

Every tree carries a heritable :class:`~mechatree.evolution.Genome`; seeds
inherit a mutated copy of the parent's genome. Lineages disappear when
no descendant survives long enough to reproduce. The forest *is* the
population — there is no central selection step.

This recipe:

1. Runs a 50-generation tournament from random founders.
2. Plots the per-generation lineage count.
3. Writes ``champions.json`` and prints the surviving species.

Run with::

    uv run python examples/evolve_a_forest.py
    uv run python examples/evolve_a_forest.py --gens 100 --n-init 30 --seed 7

The output ``champions.json`` is consumed unchanged by
``mt.load_champion(path)`` and by ``scripts/strategies_single_tree.py``.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import plotly.graph_objects as go

import mechatree as mt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "forest.yaml",
        help="YAML config (defaults to examples/forest.yaml).",
    )
    parser.add_argument("--gens", type=int, default=50, help="Number of generations.")
    parser.add_argument("--n-init", type=int, default=20, help="Initial number of founders.")
    parser.add_argument("--n-max", type=int, default=100, help="Carrying capacity.")
    parser.add_argument("--seed", type=int, default=42, help="Master RNG seed.")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("out"),
        help="Directory for champions.json + archive snapshots.",
    )
    parser.add_argument(
        "--archive-every",
        type=int,
        default=None,
        help="Write archive_{gen:08d}.json every K generations (default: only the end).",
    )
    parser.add_argument("--no-show", action="store_true", help="Skip the plotly figures.")
    args = parser.parse_args()

    cfg = mt.load_config(args.config)
    cfg = replace(
        cfg,
        forest=replace(cfg.forest, n_trees_init=args.n_init, n_trees_max=args.n_max),
    )

    args.out.mkdir(parents=True, exist_ok=True)
    champions_path = args.out / "champions.json"

    print(
        f"=== evolve_a_forest: {args.gens} gens, "
        f"n_init={args.n_init}, n_max={args.n_max}, seed={args.seed} ==="
    )
    result = mt.run_tournament(
        cfg,
        n_generations=args.gens,
        seed=args.seed,
        archive_every=args.archive_every,
        archive_dir=args.out if args.archive_every is not None else None,
        champions_path=champions_path,
    )

    print()
    print(f"final: {len(result.forest.trees)} trees")
    print(f"       {result.history[-1].n_lineages_alive} lineages alive of {args.n_init} founders")
    print(
        f"       {sum(h.n_born for h in result.history)} births, "
        f"{sum(h.n_died for h in result.history)} deaths over the run"
    )
    print(f"wrote champions to {champions_path}")
    print()
    for sp in mt.load_all_champions(champions_path):
        _safety, _alloc, angles, non_coding = sp
        print(
            f"  species {non_coding['species_id']}: "
            f"{non_coding['n_members']} members, "
            f"champion n_seeds={non_coding['champion_n_seeds']}"
        )

    if args.no_show:
        return

    gens = [s.generation for s in result.history]
    n_lineages = [s.n_lineages_alive for s in result.history]
    n_trees = [s.n_trees for s in result.history]

    fig = mt.figstyle.figure(size="full", aspect=9 / 5)
    fig.add_trace(
        go.Scatter(
            x=gens,
            y=n_lineages,
            name="lineages alive",
            line=dict(color=mt.figstyle.COLORS["red"], width=2),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=gens,
            y=n_trees,
            name="trees alive",
            line=dict(color=mt.figstyle.COLORS["green"], width=2, dash="dash"),
            yaxis="y2",
        )
    )
    fig.update_layout(
        title=f"Darwinian island: {args.gens} gens, seed={args.seed}",
        xaxis_title="generation",
        yaxis=dict(
            title=dict(text="lineages alive", font=dict(color=mt.figstyle.COLORS["red"])),
            tickfont=dict(color=mt.figstyle.COLORS["red"]),
        ),
        yaxis2=dict(
            title=dict(text="trees alive", font=dict(color=mt.figstyle.COLORS["green"])),
            tickfont=dict(color=mt.figstyle.COLORS["green"]),
            overlaying="y",
            side="right",
        ),
    )
    fig.show()


if __name__ == "__main__":
    main()
