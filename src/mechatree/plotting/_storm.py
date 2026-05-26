"""Storm-replay plot: multi-panel 3D view of the pre-storm canopy with
pruned branches in black on each fixed-point iteration.

Pairs with :func:`mechatree.wind.replay.run_storm_replay`. One panel
per snapshot; surviving branches in their normal colour, branches that
have already been cut by this iteration overlaid in black on top of
the pre-storm canopy so the eye can see "what just fell".
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from mechatree.wind.replay import StormPreSnapshot, StormSnapshot


def _segments_xyz(
    starts: np.ndarray, axes: np.ndarray, lengths: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build line-segment x/y/z arrays with NaN separators between segments.

    Each branch contributes three points: start, tip, NaN. Same trick
    used by :func:`plot_tree_3d` for the fast ``style="lines"`` path.
    """
    n = starts.shape[0]
    if n == 0:
        empty = np.empty(0, dtype=float)
        return empty, empty, empty
    tips = starts + axes * lengths[:, None]
    xs = np.empty(3 * n, dtype=float)
    ys = np.empty(3 * n, dtype=float)
    zs = np.empty(3 * n, dtype=float)
    xs[0::3] = starts[:, 0]
    xs[1::3] = tips[:, 0]
    xs[2::3] = np.nan
    ys[0::3] = starts[:, 1]
    ys[1::3] = tips[:, 1]
    ys[2::3] = np.nan
    zs[0::3] = starts[:, 2]
    zs[1::3] = tips[:, 2]
    zs[2::3] = np.nan
    return xs, ys, zs


def plot_storm_replay(
    pre: StormPreSnapshot,
    snapshots: list[StormSnapshot],
    *,
    alive_color: str | None = None,
    cut_color: str = "black",
    line_width: float = 2.0,
    cut_line_width: float = 3.0,
    n_cols: int = 3,
):
    """Render a multi-panel storm replay.

    One panel per snapshot. For panel ``k``: surviving branches at the
    end of iteration ``k`` drawn in ``alive_color``; branches that were
    alive at iteration 0 but are gone by iteration ``k`` drawn in
    ``cut_color``. Iteration-0 panel has no cut branches (all alive).

    Parameters
    ----------
    pre
        Pre-storm geometry from :func:`mechatree.wind.replay.run_storm_replay`.
    snapshots
        Per-iteration snapshots from the same call.
    alive_color, cut_color
        Plotly colour strings. ``alive_color`` defaults to
        ``figstyle.COLORS['green']`` if ``None``.
    line_width, cut_line_width
        Plotly line widths for surviving and cut branches respectively.
        Cut branches are drawn slightly thicker so they read clearly
        against the alive ones.
    n_cols
        Subplot grid width. Rows are computed automatically.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    from mechatree.plotting import figstyle

    alive_color = alive_color if alive_color is not None else figstyle.COLORS["green"]

    n_snaps = len(snapshots)
    n_rows = (n_snaps + n_cols - 1) // n_cols
    titles = []
    for s in snapshots:
        if s.iteration == 0:
            titles.append("pre-storm")
        else:
            titles.append(
                f"iter {s.iteration}: −{s.n_pruned_this_iter} branches "
                f"on {s.n_trees_affected_this_iter} trees"
            )

    specs = [[{"type": "scene"} for _ in range(n_cols)] for _ in range(n_rows)]
    fig = make_subplots(
        rows=n_rows,
        cols=n_cols,
        specs=specs,
        subplot_titles=titles,
        horizontal_spacing=0.01,
        vertical_spacing=0.05,
    )
    # Carry the mechatree template through subplots — make_subplots
    # drops the inherited layout, so apply it explicitly.
    figstyle.apply()
    fig.update_layout(template="mechatree", showlegend=False)

    for k, snap in enumerate(snapshots):
        row = k // n_cols + 1
        col = k % n_cols + 1

        # Collect per-tree alive / cut segment arrays.
        alive_x, alive_y, alive_z = [], [], []
        cut_x, cut_y, cut_z = [], [], []
        for tree_geom, mask in zip(pre.per_tree_geometry, snap.alive_mask_per_tree, strict=True):
            s = tree_geom["start"]
            a = tree_geom["axis"]
            L = tree_geom["L"]
            if mask.any():
                xs, ys, zs = _segments_xyz(s[mask], a[mask], L[mask])
                alive_x.append(xs)
                alive_y.append(ys)
                alive_z.append(zs)
            cut = ~mask
            if cut.any():
                xs, ys, zs = _segments_xyz(s[cut], a[cut], L[cut])
                cut_x.append(xs)
                cut_y.append(ys)
                cut_z.append(zs)

        if alive_x:
            fig.add_trace(
                go.Scatter3d(
                    x=np.concatenate(alive_x),
                    y=np.concatenate(alive_y),
                    z=np.concatenate(alive_z),
                    mode="lines",
                    line={"color": alive_color, "width": line_width},
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=row,
                col=col,
            )
        if cut_x:
            fig.add_trace(
                go.Scatter3d(
                    x=np.concatenate(cut_x),
                    y=np.concatenate(cut_y),
                    z=np.concatenate(cut_z),
                    mode="lines",
                    line={"color": cut_color, "width": cut_line_width},
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=row,
                col=col,
            )

        # Lock equal aspect + hide axis cruft per scene; uses the
        # mechatree template's 3D scene defaults otherwise.
        scene_name = f"scene{k + 1}" if k > 0 else "scene"
        fig.update_layout(
            **{
                scene_name: {
                    "aspectmode": "data",
                    "xaxis": {"visible": False},
                    "yaxis": {"visible": False},
                    "zaxis": {"visible": False},
                }
            }
        )

    fig.update_layout(
        height=350 * n_rows,
        width=400 * n_cols,
        margin={"l": 0, "r": 0, "t": 30, "b": 0},
    )
    return fig


__all__ = ["plot_storm_replay"]
