"""3D wireframe plotting of a property-map ``PyTree`` (random-growth demo).

Port of the 2017 intern's ``mod_3Dplot.pyx`` to plotly.
"""

from __future__ import annotations

import math


def plot_3d(tree, leaf_length: float = 0.2):
    """Render a 3D wireframe of ``tree``.

    Each branch must expose ``x``, ``y``, ``z``, ``theta`` and ``phi``
    properties. Branches with no children are decorated with a short
    green leaf segment in the tip direction.

    Parameters
    ----------
    tree : PyTree
        The tree to render.
    leaf_length : float
        Length of the green leaf segment drawn past a leaf-branch's tip.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go

    from mechatree.plotting import figstyle

    branch_x: list[float | None] = []
    branch_y: list[float | None] = []
    branch_z: list[float | None] = []
    leaf_x: list[float | None] = []
    leaf_y: list[float | None] = []
    leaf_z: list[float | None] = []

    for n in range(tree.get_number_of_branches()):
        x1 = tree.get_property(n, "x")
        y1 = tree.get_property(n, "y")
        z1 = tree.get_property(n, "z")
        theta = tree.get_property(n, "theta")
        phi = tree.get_property(n, "phi")

        dx = math.sin(theta) * math.cos(phi)
        dy = math.sin(theta) * math.sin(phi)
        dz = math.cos(theta)
        x2, y2, z2 = x1 + dx, y1 + dy, z1 + dz

        branch_x.extend([x1, x2, None])
        branch_y.extend([y1, y2, None])
        branch_z.extend([z1, z2, None])

        if tree.get_number_of_children(n) < 1:
            leaf_x.extend([x2, x2 + leaf_length * dx, None])
            leaf_y.extend([y2, y2 + leaf_length * dy, None])
            leaf_z.extend([z2, z2 + leaf_length * dz, None])

    fig = figstyle.figure_3d(size="half", show_axes=False)
    fig.add_trace(
        go.Scatter3d(
            x=branch_x,
            y=branch_y,
            z=branch_z,
            mode="lines",
            line=dict(color=figstyle.COLORS["black"], width=3),  # noqa: C408
            hoverinfo="skip",
            showlegend=False,
        )
    )
    if leaf_x:
        fig.add_trace(
            go.Scatter3d(
                x=leaf_x,
                y=leaf_y,
                z=leaf_z,
                mode="lines",
                line=dict(color=figstyle.COLORS["green"], width=6),  # noqa: C408
                hoverinfo="skip",
                showlegend=False,
            )
        )
    fig.update_scenes(aspectmode="data")
    return fig


__all__ = ["plot_3d"]
