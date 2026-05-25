"""Unified plotly styling for MechaTree figures — MATLAB / Nat Comm 2017 look.

Use it from a notebook or script::

    from mechatree.plotting import figstyle

    figstyle.apply()  # register + activate

    fig = figstyle.figure(size="half", aspect=4 / 3)
    fig.add_scatter(x=x, y=y, line=dict(color=figstyle.COLORS["red"]))
    figstyle.save(fig, "fig_demo")  # -> figures/fig_demo.pdf

    fig3 = figstyle.figure_3d(size="full")
    fig = figstyle.subplots(size="full", rows=1, cols=4)

API mirrors :mod:`softmobility.classes.figstyle` (matplotlib counterpart in the
SoftMobility repo) but is plotly-native. Differs from SoftMobility in one
detail: ``ticks="inside"`` (true MATLAB default) rather than ``"outside"``.

The module also ships the four candidate Strahler colormaps from
``../Eloy2017/plot_stat_single_tree.m`` (jet, cool, parula) plus the legacy
MechaTree rainbow, so notebooks can A/B them via :func:`set_strahler_cmap`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

from mechatree.plotting._palette import RAINBOW_STRAHLER

# ---------------------------------------------------------------------------
# Module-level constants — kept name-compatible with softmobility.figstyle.
# ---------------------------------------------------------------------------

COLORS: dict[str, str] = {
    "red": "#DA3B26",
    "red_25": "#F8CCC7",
    "blue": "#0076BA",
    "blue_light": "#7CB9E8",
    "blue_25": "#A0D3FA",
    "black": "#000000",
    "grey": "#7F7F7F",
    "green": "#2CA02C",  # leaves / live tree biomass
    "brown": "#8B4513",  # trunks / dead biomass
}

# (width_px, height_px) at 96 dpi. Widths match SoftMobility's matplotlib SIZES
# (full = 17.4 cm, half = 8.6 cm, third = 5.7 cm); heights default to 4:3.
SIZES: dict[str, tuple[int, int]] = {
    "full": (658, 494),
    "half": (326, 245),
    "third": (218, 164),
}

FONT: dict[str, Any] = dict(family="Helvetica, Arial, sans-serif", size=11, color="black")

# Plotly's renderer applies a hidden 1.25× scale to axis ticklabels (only),
# so ``tickfont.size = 11`` actually paints at 13.75 px while titles at the
# same value paint at 11 px. We compensate so every rendered text element —
# ticklabels, axis titles, subplot titles, plot title — lands at the same
# pixel size. Changing :data:`FONT` automatically updates this.
_TICK_SCALE: float = 1.25


def _tickfont() -> dict[str, Any]:
    """:data:`FONT` with the size scaled to compensate for plotly's 1.25× tick boost."""
    f = dict(FONT)
    f["size"] = FONT["size"] / _TICK_SCALE
    return f


TEMPLATE_NAME: str = "mechatree"


# ---------------------------------------------------------------------------
# Strahler colormaps — 10-stop tables sampled from MATLAB's `colormap(name)`.
# Index = Strahler order k (k = 0..9). The MATLAB script `plot_stat_single_tree.m`
# uses `colormap(cool)` for the Tokunaga/Strahler-order plot (line 37) and
# `colormap(jet)` for the 3D tree keyed by generation (line 58).
# ---------------------------------------------------------------------------

_JET: tuple[tuple[float, float, float], ...] = (
    (0.0000, 0.0000, 0.5000),
    (0.0000, 0.0000, 0.9991),
    (0.0000, 0.3784, 1.0000),
    (0.0000, 0.8333, 1.0000),
    (0.3004, 1.0000, 0.6673),
    (0.6673, 1.0000, 0.3004),
    (1.0000, 0.9012, 0.0000),
    (1.0000, 0.4800, 0.0000),
    (0.9991, 0.0733, 0.0000),
    (0.5000, 0.0000, 0.0000),
)

_COOL: tuple[tuple[float, float, float], ...] = (
    (0.0000, 1.0000, 1.0000),
    (0.1098, 0.8902, 1.0000),
    (0.2196, 0.7804, 1.0000),
    (0.3333, 0.6667, 1.0000),
    (0.4431, 0.5569, 1.0000),
    (0.5569, 0.4431, 1.0000),
    (0.6667, 0.3333, 1.0000),
    (0.7804, 0.2196, 1.0000),
    (0.8902, 0.1098, 1.0000),
    (1.0000, 0.0000, 1.0000),
)

# MATLAB R2014b+ default colormap. matplotlib doesn't ship parula, so these
# stops are pre-computed (k=0 anchored at MATLAB's documented parula start).
_PARULA: tuple[tuple[float, float, float], ...] = (
    (0.2081, 0.1663, 0.5292),
    (0.2422, 0.1504, 0.6603),
    (0.2760, 0.3470, 0.8110),
    (0.2070, 0.5300, 0.8710),
    (0.1220, 0.6580, 0.8360),
    (0.1080, 0.7600, 0.6900),
    (0.4120, 0.8110, 0.4430),
    (0.7370, 0.8110, 0.2130),
    (0.9650, 0.7620, 0.1970),
    (0.9760, 0.9840, 0.0540),
)

STRAHLER_CMAPS: dict[str, tuple[tuple[float, float, float], ...]] = {
    "jet": _JET,
    "cool": _COOL,
    "parula": _PARULA,
    "rainbow": RAINBOW_STRAHLER[1:11] if len(RAINBOW_STRAHLER) >= 11 else RAINBOW_STRAHLER,
}

DEFAULT_STRAHLER_CMAP: str = "jet"

_active_strahler_cmap: str = DEFAULT_STRAHLER_CMAP


def set_strahler_cmap(name: str) -> None:
    """Activate one of :data:`STRAHLER_CMAPS` as the default Strahler palette."""
    if name not in STRAHLER_CMAPS:
        raise ValueError(f"Unknown Strahler cmap {name!r}; choose from {sorted(STRAHLER_CMAPS)}")
    global _active_strahler_cmap
    _active_strahler_cmap = name


def get_strahler_cmap() -> str:
    """Return the name of the active Strahler colormap."""
    return _active_strahler_cmap


def strahler_color(order: int, *, cmap: str | None = None, alpha: float = 1.0) -> str:
    """Return a CSS ``rgb(...)`` / ``rgba(...)`` string for Strahler ``order``.

    Defaults to the active palette set by :func:`set_strahler_cmap`. Order is
    clipped to ``[0, 9]`` — MechaTree's Strahler ladders never exceed that in
    practice (see ``_palette.MAX_STRAHLER``).
    """
    table = STRAHLER_CMAPS[cmap if cmap is not None else _active_strahler_cmap]
    k = max(0, min(len(table) - 1, int(order)))
    r, g, b = table[k]
    if alpha >= 1.0:
        return f"rgb({int(255 * r)},{int(255 * g)},{int(255 * b)})"
    return f"rgba({int(255 * r)},{int(255 * g)},{int(255 * b)},{alpha:.3f})"


# ---------------------------------------------------------------------------
# Template construction.
# ---------------------------------------------------------------------------


def _axis_style() -> dict[str, Any]:
    return dict(
        showgrid=False,
        zeroline=False,
        showline=True,
        linecolor="black",
        linewidth=1,
        mirror=True,
        ticks="inside",
        ticklen=4,
        tickwidth=1,
        tickcolor="black",
        tickfont=_tickfont(),
        # ``standoff`` is the gap (in px) between the tick labels and the axis
        # title. Plotly's default is ~20 px which sits visually far. Tighten
        # to roughly half a line-height (~6 px) so the title hugs the
        # ticklabels — same compact look as MATLAB / SoftMobility figures.
        title=dict(font=FONT, standoff=6),
    )


def _scene_axis_style() -> dict[str, Any]:
    return dict(
        showgrid=False,
        zeroline=False,
        showline=True,
        linecolor="black",
        linewidth=1,
        mirror=True,
        ticks="inside",
        backgroundcolor="white",
        showbackground=False,
        tickfont=_tickfont(),
        title=dict(font=FONT),
    )


def _build_template() -> go.layout.Template:
    axis = _axis_style()
    scene_axis = _scene_axis_style()
    layout = go.Layout(
        paper_bgcolor="white",
        plot_bgcolor="white",
        font=FONT,
        title=dict(font=FONT),
        annotationdefaults=dict(font=FONT),
        margin=dict(l=55, r=15, t=15, b=45),
        colorway=[
            COLORS["red"],
            COLORS["blue"],
            COLORS["grey"],
            COLORS["green"],
            COLORS["brown"],
            COLORS["black"],
        ],
        xaxis=axis,
        yaxis=axis,
        legend=dict(font=FONT, bgcolor="rgba(255,255,255,0)"),
        scene=dict(
            xaxis=scene_axis,
            yaxis=scene_axis,
            zaxis=scene_axis,
            aspectmode="cube",
            camera=dict(
                projection=dict(type="orthographic"),
                eye=dict(x=1.8, y=1.8, z=1.4),
                center=dict(x=0.0, y=0.0, z=0.0),
                up=dict(x=0.0, y=0.0, z=1.0),
            ),
        ),
    )
    return go.layout.Template(layout=layout)


def _shrink_subplot_titles(fig: go.Figure) -> None:
    """Force subplot-title annotations to the FONT size.

    ``plotly.subplots.make_subplots`` hard-codes the subplot-title font
    size at creation time (16 pt), so the template's ``annotationdefaults``
    is ignored for those specific annotations. We walk the annotation
    list and overwrite ``font.size`` / ``font.family`` explicitly.
    """
    for ann in fig.layout.annotations:
        ann.font.family = FONT["family"]
        ann.font.size = FONT["size"]
        ann.font.color = FONT["color"]


def apply() -> None:
    """Register the MechaTree template under :data:`TEMPLATE_NAME` and activate it.

    Re-call after editing module-level constants (:data:`COLORS`, :data:`FONT`)
    to propagate the change.
    """
    pio.templates[TEMPLATE_NAME] = _build_template()
    pio.templates.default = TEMPLATE_NAME


def _ensure_registered() -> None:
    if TEMPLATE_NAME not in pio.templates:
        pio.templates[TEMPLATE_NAME] = _build_template()


# ---------------------------------------------------------------------------
# Canvas factories.
# ---------------------------------------------------------------------------


def _resolve_size(size: str, aspect: float) -> tuple[int, int]:
    if size not in SIZES:
        raise ValueError(f"Unknown size {size!r}; choose from {sorted(SIZES)}")
    width = SIZES[size][0]
    height = int(round(width / aspect))
    return width, height


def figure(size: str = "full", aspect: float = 4 / 3, **layout_kwargs: Any) -> go.Figure:
    """Return a styled :class:`go.Figure` at one of the named paper widths.

    Forces the MechaTree template regardless of whether :func:`apply` has been
    called, so library helpers stay consistent even when the user hasn't.
    """
    _ensure_registered()
    width, height = _resolve_size(size, aspect)
    fig = go.Figure()
    fig.update_layout(template=TEMPLATE_NAME, width=width, height=height, **layout_kwargs)
    return fig


def subplots(
    size: str = "full",
    aspect: float = 4 / 3,
    rows: int = 1,
    cols: int = 1,
    **make_subplots_kwargs: Any,
) -> go.Figure:
    """Multi-panel styled figure. Thin wrapper around :func:`make_subplots`."""
    _ensure_registered()
    width, height = _resolve_size(size, aspect)
    fig = make_subplots(rows=rows, cols=cols, **make_subplots_kwargs)
    fig.update_layout(template=TEMPLATE_NAME, width=width, height=height)
    _shrink_subplot_titles(fig)
    return fig


def figure_3d(
    size: str = "full",
    aspect: float = 1.0,
    show_axes: bool = False,
    **layout_kwargs: Any,
) -> go.Figure:
    """Styled :class:`go.Figure` with an orthographic 3D scene.

    With ``show_axes=False`` (the default) the cube edges, ticks, and labels
    are hidden — appropriate for tree / canopy renders. ``show_axes=True``
    restores the template's axis lines so the scene reads as a labelled plot.
    """
    _ensure_registered()
    width, height = _resolve_size(size, aspect)
    fig = go.Figure()
    layout_extra: dict[str, Any] = dict(layout_kwargs)
    if not show_axes:
        layout_extra.setdefault(
            "scene",
            dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                aspectmode="cube",
                camera=dict(
                    projection=dict(type="orthographic"),
                    eye=dict(x=1.8, y=1.8, z=1.4),
                    up=dict(x=0.0, y=0.0, z=1.0),
                ),
            ),
        )
    fig.update_layout(
        template=TEMPLATE_NAME,
        width=width,
        height=height,
        margin=dict(l=0, r=0, t=0, b=0),
        **layout_extra,
    )
    return fig


# ---------------------------------------------------------------------------
# PDF export.
# ---------------------------------------------------------------------------


def save(fig: go.Figure, name: str, figdir: str | Path = "figures") -> Path:
    """Write ``fig`` as a vector PDF to ``figdir/<name>.pdf`` via kaleido."""
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    path = figdir / f"{name}.pdf"
    fig.write_image(str(path), format="pdf")
    return path


__all__ = [
    "COLORS",
    "DEFAULT_STRAHLER_CMAP",
    "FONT",
    "SIZES",
    "STRAHLER_CMAPS",
    "TEMPLATE_NAME",
    "apply",
    "figure",
    "figure_3d",
    "get_strahler_cmap",
    "save",
    "set_strahler_cmap",
    "strahler_color",
    "subplots",
]
