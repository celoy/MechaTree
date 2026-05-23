"""3D / top-down plots that read the typed mechanics fields on a tree.

The earlier :func:`mechatree.plotting.plot_3d` reads property-map keys
(``x``, ``y``, ``z``, ``theta``, ``phi``) used by the legacy ``random_growth``
demos. Trees produced by :func:`mechatree.simulate.grow_tree` and
:class:`mechatree.forest.Forest` store geometry in typed fields
(``location``, ``unit_t``, ``length``, ``diameter``) — these helpers
plot from those directly.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt


def plot_tree_3d(tree, ax=None, leaf_color: str = "green", leaf_size: float = 20.0):
    """Render a 3D view of a mechanics tree.

    Each branch is drawn as a line from ``location`` to ``location + length * unit_t``
    with linewidth scaled by diameter. Leaves (childless branches) are
    accented with a coloured marker.

    Parameters
    ----------
    tree
        A :class:`PyTree` whose branches have typed mechanics fields set
        (as produced by :func:`mechatree.simulate.grow_tree`).
    ax
        Existing 3D axes; a new figure is made if omitted.
    leaf_color
        Marker colour for leaf tips.
    leaf_size
        Marker size for leaf tips.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    else:
        fig = ax.figure

    ax.set_axis_off()

    # Pre-compute diameter range for line-width scaling.
    n = tree.get_number_of_branches()
    if n == 0:
        return fig

    diameters = [tree.get_diameter(i) for i in range(n)]
    d_max = max(diameters) if diameters else 1.0
    d_max = max(d_max, 1e-6)

    leaf_xs, leaf_ys, leaf_zs = [], [], []

    for i in range(n):
        x, y, z = tree.get_location(i)
        tx, ty, tz = tree.get_unit_t(i)
        L = tree.get_length(i)
        d = tree.get_diameter(i)

        x2 = x + L * tx
        y2 = y + L * ty
        z2 = z + L * tz

        lw = max(0.5, 4.0 * d / d_max)
        ax.plot([x, x2], [y, y2], [z, z2], color="saddlebrown", linewidth=lw)

        if tree.get_number_of_children(i) == 0:
            leaf_xs.append(x2)
            leaf_ys.append(y2)
            leaf_zs.append(z2)

    if leaf_xs:
        ax.scatter(leaf_xs, leaf_ys, leaf_zs, c=leaf_color, s=leaf_size, depthshade=False)

    return fig


def plot_forest_topdown(forest, ax=None, biomass_scale: float = 1.0):
    """Top-down view of a forest: each tree is a disk at its trunk position,
    radius proportional to ``sqrt(biomass)``.

    Parameters
    ----------
    forest
        A :class:`mechatree.forest.Forest`.
    ax
        Existing 2D axes; a new figure is made if omitted.
    biomass_scale
        Multiplier for the disk radius.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.figure

    size = forest.config.forest.size
    # Forest boundary circle.
    theta = [i / 100.0 * 2.0 * math.pi for i in range(101)]
    ax.plot(
        [size * math.cos(t) for t in theta],
        [size * math.sin(t) for t in theta],
        color="black",
        linewidth=0.5,
        linestyle="--",
    )
    ax.set_xlim(-1.1 * size, 1.1 * size)
    ax.set_ylim(-1.1 * size, 1.1 * size)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    xs, ys, rs = [], [], []
    for tree in forest.trees:
        x, y, _ = tree.get_location(0)
        # Trunk-volume-based radius; a stand-in for "tree size".
        biomass = 0.0
        for i in range(tree.get_number_of_branches()):
            d = tree.get_diameter(i)
            biomass += tree.get_length(i) * math.pi * d * d / 4.0
        xs.append(x)
        ys.append(y)
        rs.append(biomass_scale * max(0.1, math.sqrt(biomass)))

    if xs:
        ax.scatter(xs, ys, s=[r * 50.0 for r in rs], alpha=0.5, color="forestgreen")

    return fig


__all__ = ["plot_tree_3d", "plot_forest_topdown"]
