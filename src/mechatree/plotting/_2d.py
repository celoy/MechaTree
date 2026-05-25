"""2D side-view plotting of a property-map ``PyTree`` (self-avoiding-growth demo).

Port of the 2017 intern's ``mod_plot.pyx`` to plotly.
"""

from __future__ import annotations

import math
from pathlib import Path


def plot_2d(
    tree,
    iteration: int | None = None,
    out_dir: str | Path | None = None,
):
    """Render a 2D side-view of ``tree`` as a plotly figure.

    Parameters
    ----------
    tree : PyTree
        The tree to render. Each branch must expose ``x``, ``y``, ``L``
        and ``theta`` properties.
    iteration : int, optional
        Iteration number, used to name the output PNG. Required if
        ``out_dir`` is provided.
    out_dir : path-like, optional
        If provided, the figure is saved to ``out_dir/figNNNN.png``
        (4-digit zero-padded) via ``fig.write_image`` (needs ``kaleido``).
        The directory is created if it does not exist.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go

    from mechatree.plotting import figstyle

    xs: list[float | None] = []
    ys: list[float | None] = []
    for n in range(tree.get_number_of_branches()):
        x1 = tree.get_property(n, "x")
        y1 = tree.get_property(n, "y")
        length = tree.get_property(n, "L")
        theta = tree.get_property(n, "theta")
        x2 = x1 + length * math.sin(theta)
        y2 = y1 + length * math.cos(theta)
        xs.extend([x1, x2, None])
        ys.extend([y1, y2, None])

    fig = figstyle.figure(size="half", aspect=1.0)
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            mode="lines",
            line=dict(color=figstyle.COLORS["red"], width=1),  # noqa: C408
            hoverinfo="skip",
            showlegend=False,
        )
    )
    fig.update_layout(
        xaxis=dict(visible=False, scaleanchor="y", scaleratio=1),  # noqa: C408
        yaxis=dict(visible=False),  # noqa: C408
        margin=dict(l=0, r=0, t=0, b=0),  # noqa: C408
    )

    if out_dir is not None:
        if iteration is None:
            raise ValueError("plot_2d requires `iteration` when `out_dir` is set.")
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        fig.write_image(str(out_path / f"fig{iteration:04d}.png"))

    return fig


__all__ = ["plot_2d"]
