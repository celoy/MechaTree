"""3D wireframe plotting of a PyTree, used by the random-growth demo.

Port of the 2017 intern's ``mod_3Dplot.pyx``. Pure Python.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt


def plot_3d(tree, ax=None, leaf_length: float = 0.2):
    """Render a 3D wireframe of ``tree``.

    Each branch must expose ``x``, ``y``, ``z``, ``theta`` and ``phi``
    properties. Branches with no children are decorated with a short green
    leaf segment in the tip direction.

    Parameters
    ----------
    tree : PyTree
        The tree to render.
    ax : mpl_toolkits.mplot3d.axes3d.Axes3D, optional
        Existing 3D axes to draw onto. A new figure is created if omitted.
    leaf_length : float
        Length of the green leaf segment drawn past a leaf-branch's tip.

    Returns
    -------
    matplotlib.figure.Figure
        The figure that was drawn on.
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    else:
        fig = ax.figure

    ax.set_axis_off()

    for n in range(tree.get_number_of_branches()):
        x1 = tree.get_property(n, "x")
        y1 = tree.get_property(n, "y")
        z1 = tree.get_property(n, "z")
        theta = tree.get_property(n, "theta")
        phi = tree.get_property(n, "phi")

        dx = math.sin(theta) * math.cos(phi)
        dy = math.sin(theta) * math.sin(phi)
        dz = math.cos(theta)

        x2 = x1 + dx
        y2 = y1 + dy
        z2 = z1 + dz

        ax.plot([x1, x2], [y1, y2], [z1, z2], "k", linewidth=2.5)

        if tree.get_number_of_children(n) < 1:
            xleaf = x2 + leaf_length * dx
            yleaf = y2 + leaf_length * dy
            zleaf = z2 + leaf_length * dz
            ax.plot([x2, xleaf], [y2, yleaf], [z2, zleaf], "g", linewidth=5)

    return fig
