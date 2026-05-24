"""Plotly renderers for the typed mechanics fields on a tree.

The earlier matplotlib-based ``plot_tree_3d`` used
:class:`mpl_toolkits.mplot3d.art3d.Line3DCollection`; matplotlib's 3D
backend has no real depth buffer and paints leaves on top of branches
regardless of true depth. Plotly's WebGL renderer fixes that. Returned
figures are :class:`plotly.graph_objects.Figure` — call ``fig.show()``
to open a browser tab, or display inline in a Jupyter notebook.
"""

from __future__ import annotations

import math

import numpy as np

from mechatree.plotting._palette import MAX_STRAHLER, strahler_css


def plot_tree_3d(
    tree,
    *,
    max_strahler: int = MAX_STRAHLER,
    width_min: float = 1.0,
    width_max: float = 6.0,
    show_leaves: bool = False,
    leaf_size_scale: float = 8.0,
    leaf_color: str = "green",
    sun=None,
):
    """Render a 3D view of a mechanics tree as a plotly figure.

    Branches are grouped by Strahler order — one :class:`Scatter3d` line
    trace per order, coloured from the fixed rainbow palette. The
    projection is orthographic and the aspect ratio matches the data, so
    the tree isn't visually distorted.

    Parameters
    ----------
    tree
        A :class:`PyTree` whose branches have typed mechanics fields set
        (as produced by :func:`mechatree.simulate.grow_tree`).
    max_strahler
        Highest Strahler order assumed to occur in any tree. A branch
        with Strahler ``max_strahler`` always renders as red; smaller
        trees just don't use the top buckets.
    width_min, width_max
        Plotly line widths for the thinnest and thickest branch.
    show_leaves
        If ``True``, scatter leaves at branch tips with size and alpha
        linear in absorbed light (recomputed once via the light module).
    leaf_size_scale
        Marker size at peak light when ``show_leaves=True``.
    leaf_color
        Marker colour for leaves (any CSS / plotly colour spec).
    sun
        Optional :class:`mechatree.light.Sun`. Used only when
        ``show_leaves=True``; defaults to ``Sun()``.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go

    n = tree.get_number_of_branches()
    fig = go.Figure()
    if n == 0:
        return fig

    tree.set_strahler()

    starts = np.empty((n, 3), dtype=np.float64)
    ends = np.empty((n, 3), dtype=np.float64)
    diameters = np.empty(n, dtype=np.float64)
    strahler = np.empty(n, dtype=np.int32)
    for i in range(n):
        x, y, z = tree.get_location(i)
        tx, ty, tz = tree.get_unit_t(i)
        L = tree.get_length(i)
        starts[i] = (x, y, z)
        ends[i] = (x + L * tx, y + L * ty, z + L * tz)
        diameters[i] = tree.get_diameter(i)
        strahler[i] = tree.get_strahler(i)

    d_max = max(float(diameters.max()), 1e-6)

    # One trace per Strahler order — plotly's Scatter3d only takes a single
    # line.color per trace, so this is the natural grouping.
    for order in np.unique(strahler):
        mask = strahler == order
        xs: list[float | None] = []
        ys: list[float | None] = []
        zs: list[float | None] = []
        for s, e in zip(starts[mask], ends[mask], strict=True):
            xs.extend([float(s[0]), float(e[0]), None])
            ys.extend([float(s[1]), float(e[1]), None])
            zs.extend([float(s[2]), float(e[2]), None])
        # One representative width per order, taken from the median diameter
        # in that bucket (so the trunk reads thick, twigs thin).
        d_bucket = float(np.median(diameters[mask])) if diameters[mask].size else 0.0
        width = width_min + (width_max - width_min) * (d_bucket / d_max)
        fig.add_trace(
            go.Scatter3d(
                x=xs,
                y=ys,
                z=zs,
                mode="lines",
                line=dict(color=strahler_css(int(order)), width=width),  # noqa: C408
                name=f"Strahler {int(order)}",
                hoverinfo="skip",
                showlegend=False,
            )
        )

    if show_leaves:
        from mechatree.light import Sun, extract_leaves, intercept

        sun = sun if sun is not None else Sun()
        leaves = extract_leaves([tree], n_directions=sun.n_directions)
        if leaves.n_leaves > 0:
            intercept(leaves, sun)
            light = leaves.light_per_direction.mean(axis=1)
            peak = float(light.max())
            t = light / peak if peak > 0 else np.zeros_like(light)
            # Drop fully-dark leaves so the canopy stays crisp.
            keep = t > 0.05
            if keep.any():
                fig.add_trace(
                    go.Scatter3d(
                        x=leaves.location[keep, 0],
                        y=leaves.location[keep, 1],
                        z=leaves.location[keep, 2],
                        mode="markers",
                        marker=dict(  # noqa: C408
                            size=leaf_size_scale * t[keep],
                            color=[_rgba(leaf_color, float(a)) for a in t[keep]],
                        ),
                        name="leaves",
                        hoverinfo="skip",
                        showlegend=False,
                    )
                )

    _apply_iso_layout(fig, np.vstack([starts, ends]))
    return fig


def plot_forest_topdown(forest, *, biomass_scale: float = 1.0):
    """Top-down view of a forest as a plotly figure.

    Each tree is a green disk at its trunk position with area
    proportional to trunk-volume-based biomass. The forest boundary is
    drawn as a dashed circle.

    Parameters
    ----------
    forest
        A :class:`mechatree.forest.Forest`.
    biomass_scale
        Multiplier for the disk radius.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go

    size = forest.config.forest.size
    fig = go.Figure()

    # Boundary circle.
    theta = np.linspace(0, 2 * math.pi, 101)
    fig.add_trace(
        go.Scatter(
            x=size * np.cos(theta),
            y=size * np.sin(theta),
            mode="lines",
            line=dict(color="black", width=1, dash="dash"),  # noqa: C408
            name="boundary",
            hoverinfo="skip",
            showlegend=False,
        )
    )

    xs, ys, sizes = [], [], []
    for tree in forest.trees:
        x, y, _ = tree.get_location(0)
        biomass = 0.0
        for i in range(tree.get_number_of_branches()):
            d = tree.get_diameter(i)
            biomass += tree.get_length(i) * math.pi * d * d / 4.0
        xs.append(x)
        ys.append(y)
        sizes.append(biomass_scale * max(0.1, math.sqrt(biomass)))

    if xs:
        # Plotly marker.size is the diameter in pixels — scale up the
        # geometric radii to something visible.
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="markers",
                marker=dict(  # noqa: C408
                    size=[8.0 * s for s in sizes],
                    color="forestgreen",
                    opacity=0.5,
                ),
                name="trees",
                hoverinfo="skip",
                showlegend=False,
            )
        )

    fig.update_layout(
        xaxis=dict(  # noqa: C408
            range=[-1.1 * size, 1.1 * size],
            scaleanchor="y",
            scaleratio=1,
            title="x",
        ),
        yaxis=dict(range=[-1.1 * size, 1.1 * size], title="y"),  # noqa: C408
        margin=dict(l=40, r=20, t=20, b=40),  # noqa: C408
        paper_bgcolor="white",
        plot_bgcolor="white",
        width=600,
        height=600,
    )
    return fig


# ---- helpers ---------------------------------------------------------------


def _apply_iso_layout(fig, verts: np.ndarray) -> None:
    """Set isometric/orthographic projection with data-true aspect ratio."""
    extents = verts.max(axis=0) - verts.min(axis=0)
    extents = np.where(extents > 0, extents, 1.0)
    aspect = extents / extents.max()
    fig.update_layout(
        scene=dict(  # noqa: C408
            xaxis=dict(visible=False),  # noqa: C408
            yaxis=dict(visible=False),  # noqa: C408
            zaxis=dict(visible=False),  # noqa: C408
            aspectmode="manual",
            aspectratio=dict(x=float(aspect[0]), y=float(aspect[1]), z=float(aspect[2])),  # noqa: C408
            camera=dict(  # noqa: C408
                projection=dict(type="orthographic"),  # noqa: C408
                eye=dict(x=0.9, y=0.9, z=0.45),  # noqa: C408
                center=dict(x=0.0, y=0.0, z=0.0),  # noqa: C408
                up=dict(x=0.0, y=0.0, z=1.0),  # noqa: C408
            ),
        ),
        margin=dict(l=0, r=0, t=0, b=0),  # noqa: C408
        paper_bgcolor="white",
    )


def _rgba(color: str, alpha: float) -> str:
    """Return ``rgba(r,g,b,a)`` for a CSS named colour or ``rgb()`` string."""
    # Hand-rolled lookup for the few names plotly callers might pass; for
    # anything else, assume it's already an rgb()/hex string and append alpha
    # by best effort.
    named = {
        "green": (0, 128, 0),
        "forestgreen": (34, 139, 34),
        "red": (255, 0, 0),
        "black": (0, 0, 0),
    }
    if color in named:
        r, g, b = named[color]
        return f"rgba({r},{g},{b},{alpha:.3f})"
    if color.startswith("#") and len(color) == 7:
        r = int(color[1:3], 16)
        g = int(color[3:5], 16)
        b = int(color[5:7], 16)
        return f"rgba({r},{g},{b},{alpha:.3f})"
    # Fall back to a translucent green if we can't parse.
    return f"rgba(0,128,0,{alpha:.3f})"


__all__ = ["plot_forest_topdown", "plot_tree_3d"]
