"""2D side-view plotting of a PyTree, used by the self-avoiding-growth demo.

Port of the 2017 intern's ``mod_plot.pyx``. Pure Python.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt


def plot_2d(
    tree,
    iteration: int | None = None,
    ax=None,
    out_dir: str | Path | None = None,
):
    """Render a 2D side-view of ``tree``.

    Parameters
    ----------
    tree : PyTree
        The tree to render. Each branch must expose ``x``, ``y``, ``L`` and
        ``theta`` properties.
    iteration : int, optional
        Iteration number, used to name the output PNG. Required if
        ``out_dir`` is provided.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw onto. A new figure is created if omitted.
    out_dir : path-like, optional
        If provided, the figure is saved to ``out_dir/figNNNN.png`` (4-digit
        zero-padded). The directory is created if it does not exist.

    Returns
    -------
    matplotlib.figure.Figure
        The figure that was drawn on.
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    ax.set_axis_off()

    for n in range(tree.get_number_of_branches()):
        x1 = tree.get_property(n, "x")
        y1 = tree.get_property(n, "y")
        length = tree.get_property(n, "L")
        theta = tree.get_property(n, "theta")

        x2 = x1 + length * math.sin(theta)
        y2 = y1 + length * math.cos(theta)

        ax.plot([x1, x2], [y1, y2], color="red", linewidth=1)

    if out_dir is not None:
        if iteration is None:
            raise ValueError("plot_2d requires `iteration` when `out_dir` is set.")
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path / f"fig{iteration:04d}.png")

    return fig
