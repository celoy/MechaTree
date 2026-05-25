#!/usr/bin/env python3
"""Side-by-side benchmark of the figstyle palette / font choices.

Renders three things so you can A/B them in your browser:

1. **Strahler palettes** — a real species-0 S3 champion tree at 400 yr,
   rendered four times under ``jet`` / ``cool`` / ``parula`` / ``rainbow``.
2. **Font families** — a small line plot duplicated under a few candidate
   font stacks so you can pick one that reads cleanly at 11 pt.
3. **Tick direction and frame** — the same line plot under
   ``ticks="inside"`` (current default) and ``ticks="outside"``
   (SoftMobility default), with and without a 4-sided ``mirror`` frame.

Run with::

    uv run python examples/figstyle_compare.py
    uv run python examples/figstyle_compare.py --no-show     # just print

Pick a winner per group, then flip the corresponding default in
``src/mechatree/plotting/mt.figstyle.py``:

- Strahler palette: ``DEFAULT_STRAHLER_CMAP`` (currently ``"jet"``).
- Font: edit :data:`FONT`.
- Ticks / frame: edit ``_axis_style()``.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import plotly.graph_objects as go

import mechatree as mt

REPO = Path(__file__).resolve().parents[1]


def _sample_curve() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A simple log-decaying sine — three curves to fill the legend slot."""
    x = np.linspace(0, 4 * math.pi, 200)
    y1 = np.sin(x) * np.exp(-x / 8)
    y2 = np.sin(x + math.pi / 3) * np.exp(-x / 8)
    y3 = np.sin(x + 2 * math.pi / 3) * np.exp(-x / 8)
    return x, np.stack([y1, y2, y3])


def compare_strahler_palettes(*, generations: int = 400, species_id: int = 0) -> go.Figure:
    """Render one species-N champion tree under each Strahler palette."""
    cfg = mt.load_config(REPO / "examples" / "forest.yaml")
    from dataclasses import replace

    safety, allocation, angles, _ = mt.load_champion(
        REPO / "data" / "S3_champions.json", species_id=species_id
    )
    cfg = replace(cfg, tree=replace(cfg.tree, **angles))
    tree = mt.grow_tree(
        cfg, n_generations=generations, safety=safety, allocation=allocation, seed=0
    )
    cmaps = ("jet", "cool", "parula", "rainbow")

    fig = mt.figstyle.subplots(
        size="full",
        aspect=4.0,
        rows=1,
        cols=4,
        specs=[[{"type": "scene"}] * 4],
        subplot_titles=cmaps,
    )

    # plot_tree_3d builds its own go.Figure under each palette; we copy its
    # traces and scene config into the corresponding subplot.
    original_cmap = mt.figstyle.get_strahler_cmap()
    try:
        for col, name in enumerate(cmaps, start=1):
            mt.figstyle.set_strahler_cmap(name)
            sub_fig = mt.plot_tree_3d(tree, show_leaves=True)
            for trace in sub_fig.data:
                fig.add_trace(trace, row=1, col=col)
            scene_key = "scene" if col == 1 else f"scene{col}"
            fig.layout[scene_key].update(sub_fig.layout.scene)
    finally:
        mt.figstyle.set_strahler_cmap(original_cmap)

    fig.update_layout(
        title=f"Strahler palettes (species {species_id}, {generations} yr)",
        showlegend=False,
    )
    return fig


def compare_fonts() -> go.Figure:
    """Same curves under several font stacks."""
    candidates = (
        ("Helvetica", "Helvetica, Arial, sans-serif"),
        ("Arial", "Arial, sans-serif"),
        ("Computer Modern", "Computer Modern, Latin Modern Roman, serif"),
        ("system-ui", "system-ui, -apple-system, BlinkMacSystemFont, sans-serif"),
    )

    fig = mt.figstyle.subplots(
        size="full",
        aspect=4 / 3,
        rows=2,
        cols=2,
        subplot_titles=[label for label, _ in candidates],
    )
    x, ys = _sample_curve()
    palette = [mt.figstyle.COLORS["red"], mt.figstyle.COLORS["blue"], mt.figstyle.COLORS["grey"]]

    for idx, (_label, stack) in enumerate(candidates):
        row, col = divmod(idx, 2)
        row += 1
        col += 1
        for k in range(3):
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=ys[k],
                    mode="lines",
                    line=dict(color=palette[k], width=1.5),  # noqa: C408
                    name=f"phase {k}",
                    showlegend=False,
                ),
                row=row,
                col=col,
            )
        fig.update_xaxes(title="x", row=row, col=col, tickfont=dict(family=stack))  # noqa: C408
        fig.update_yaxes(title="y", row=row, col=col, tickfont=dict(family=stack))  # noqa: C408
        # Force the annotation (subplot title) to use this font too.
        ann_idx = idx
        fig.layout.annotations[ann_idx].update(font=dict(family=stack, size=12))  # noqa: C408
    fig.update_layout(title="Font comparison (Helvetica / Arial / CM / system-ui)")
    return fig


def compare_frames() -> go.Figure:
    """ticks inside/outside × frame yes/no."""
    x, ys = _sample_curve()
    cases = (
        ("inside ticks, 4-sided frame", "inside", True),
        ("inside ticks, no top/right", "inside", False),
        ("outside ticks, 4-sided frame", "outside", True),
        ("outside ticks, no top/right", "outside", False),
    )
    fig = mt.figstyle.subplots(
        size="full",
        aspect=4 / 3,
        rows=2,
        cols=2,
        subplot_titles=[label for label, *_ in cases],
    )
    palette = [mt.figstyle.COLORS["red"], mt.figstyle.COLORS["blue"], mt.figstyle.COLORS["grey"]]
    for idx, (_, tick_dir, mirror) in enumerate(cases):
        row, col = divmod(idx, 2)
        row += 1
        col += 1
        for k in range(3):
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=ys[k],
                    mode="lines",
                    line=dict(color=palette[k], width=1.5),  # noqa: C408
                    showlegend=False,
                ),
                row=row,
                col=col,
            )
        fig.update_xaxes(
            ticks=tick_dir, mirror=mirror, showline=True, linecolor="black", row=row, col=col
        )
        fig.update_yaxes(
            ticks=tick_dir, mirror=mirror, showline=True, linecolor="black", row=row, col=col
        )
    fig.update_layout(title="Tick direction × frame")
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-show", action="store_true", help="Don't open browser tabs.")
    parser.add_argument(
        "--generations", type=int, default=400, help="Years to grow the benchmark tree."
    )
    parser.add_argument(
        "--species-id", type=int, default=0, help="Champion species (0 or 1 in S3_champions.json)."
    )
    args = parser.parse_args()

    mt.figstyle.apply()

    figs = {
        "strahler": compare_strahler_palettes(
            generations=args.generations, species_id=args.species_id
        ),
        "fonts": compare_fonts(),
        "frames": compare_frames(),
    }

    print("Built 3 comparison figures:")
    for name in figs:
        print(f"  - {name}")

    if not args.no_show:
        for fig in figs.values():
            fig.show()


if __name__ == "__main__":
    main()
