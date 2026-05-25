"""Tests for :mod:`mechatree.plotting.figstyle`."""

from __future__ import annotations

import re

import plotly.graph_objects as go
import plotly.io as pio
import pytest

from mechatree.plotting import figstyle


def test_apply_registers_template() -> None:
    figstyle.apply()
    assert figstyle.TEMPLATE_NAME in pio.templates
    assert pio.templates.default == figstyle.TEMPLATE_NAME


def test_figure_returns_styled_figure() -> None:
    fig = figstyle.figure(size="half", aspect=4 / 3)
    assert isinstance(fig, go.Figure)
    assert fig.layout.template.layout.paper_bgcolor == "white"
    assert fig.layout.template.layout.plot_bgcolor == "white"
    assert fig.layout.template.layout.xaxis.mirror is True
    assert fig.layout.template.layout.xaxis.ticks == "inside"
    assert fig.layout.template.layout.xaxis.showgrid is False


def test_sizes_known_values() -> None:
    assert figstyle.figure(size="half").layout.width == 326
    assert figstyle.figure(size="full").layout.width == 658
    assert figstyle.figure(size="third").layout.width == 218


def test_unknown_size_raises() -> None:
    with pytest.raises(ValueError, match="Unknown size"):
        figstyle.figure(size="quarter")
    with pytest.raises(ValueError, match="Unknown size"):
        figstyle.subplots(size="quarter")
    with pytest.raises(ValueError, match="Unknown size"):
        figstyle.figure_3d(size="quarter")


def test_colors_are_hex() -> None:
    pattern = re.compile(r"^#[0-9A-Fa-f]{6}$")
    for name, value in figstyle.COLORS.items():
        assert pattern.fullmatch(value), f"COLORS[{name!r}] = {value!r} is not #RRGGBB"


def test_strahler_cmap_switch() -> None:
    assert figstyle.get_strahler_cmap() in figstyle.STRAHLER_CMAPS
    for name in ("jet", "cool", "parula", "rainbow"):
        figstyle.set_strahler_cmap(name)
        assert figstyle.get_strahler_cmap() == name
        # 10 distinct stops; orders 0 and 9 must differ as strings.
        assert figstyle.strahler_color(0) != figstyle.strahler_color(9)
    figstyle.set_strahler_cmap("jet")


def test_strahler_color_unknown_cmap_raises() -> None:
    with pytest.raises(ValueError, match="Unknown Strahler cmap"):
        figstyle.set_strahler_cmap("viridis")
