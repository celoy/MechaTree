"""Plotly ports of the MATLAB analysis scripts in
``legacy/matlab/``.

Each helper takes pre-computed simulation data (a tree, a forest, or a
history of stats objects) and returns the
:class:`plotly.graph_objects.Figure` it drew. The example scripts in
``examples/`` are thin CLI wrappers that drive a simulation and call
these — call ``fig.show()`` to open a browser tab, or display inline in
a Jupyter notebook.

- :func:`plot_self_thinning` — population, biomass, and N-vs-M (log-log).
  Port of ``self_thinning.m``.
- :func:`plot_allocation` — per-generation bookkeeping on a single tree.
  Port of ``plot_allocation_vs_t.m``.
- :func:`plot_strahler_diagnostics` — 2×2 panel of Strahler-order
  statistics + Leonardo's-rule histogram. Port of
  ``plot_stat_single_tree.m`` + ``Fractal_dim.m`` +
  ``plot_area_preservation_1tree.m``.
"""

from __future__ import annotations

from collections.abc import Iterable

from mechatree._core import PyTree


def plot_self_thinning(history):
    """Population, biomass, and N-vs-M log-log on a 1×3 grid.

    Parameters
    ----------
    history:
        Iterable of ``(generation, n_trees, biomass)`` tuples, or
        :class:`mechatree.forest.ForestStats` objects.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go

    from mechatree.plotting import figstyle

    gens, ns, ms = _unpack_thinning(history)

    fig = figstyle.subplots(
        size="full",
        aspect=3.0,
        rows=1,
        cols=3,
        subplot_titles=(
            "Population over time",
            "Biomass over time",
            "Self-thinning (log-log)",
        ),
    )

    fig.add_trace(
        go.Scatter(
            x=gens,
            y=ns,
            mode="lines",
            line=dict(color=figstyle.COLORS["green"]),  # noqa: C408
            showlegend=False,
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=gens,
            y=ms,
            mode="lines",
            line=dict(color=figstyle.COLORS["brown"]),  # noqa: C408
            showlegend=False,
        ),
        row=1,
        col=2,
    )

    pairs = [(n, m) for n, m in zip(ns, ms, strict=True) if n > 0 and m > 0]
    if pairs:
        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines+markers",
                line=dict(color=figstyle.COLORS["black"]),  # noqa: C408
                marker=dict(size=4),  # noqa: C408
                name="data",
                showlegend=False,
            ),
            row=1,
            col=3,
        )
        # Reference line: M ∝ N^{-3/2} (Yoda 1963).
        n_anchor = min(xs)
        m_anchor = max(ys)
        ref_n = [min(xs), max(xs)]
        ref_m = [m_anchor * (n_anchor / n) ** 1.5 for n in ref_n]
        fig.add_trace(
            go.Scatter(
                x=ref_n,
                y=ref_m,
                mode="lines",
                line=dict(color=figstyle.COLORS["red"], dash="dash"),  # noqa: C408
                name="M ∝ N^{-3/2}",
            ),
            row=1,
            col=3,
        )
        fig.update_xaxes(type="log", row=1, col=3)
        fig.update_yaxes(type="log", row=1, col=3)

    fig.update_xaxes(title_text="generation", row=1, col=1)
    fig.update_yaxes(title_text="N (trees alive)", row=1, col=1)
    fig.update_xaxes(title_text="generation", row=1, col=2)
    fig.update_yaxes(title_text="M (total biomass)", row=1, col=2)
    fig.update_xaxes(title_text="N (trees alive)", row=1, col=3)
    fig.update_yaxes(title_text="M (total biomass)", row=1, col=3)

    fig.update_layout(margin=dict(l=60, r=20, t=40, b=50))  # noqa: C408
    return fig


def plot_allocation(stats_history, volume_twig: float = 1.0):
    """Per-generation allocation diagnostics for a single tree.

    Reproduces the columns of the Fortran ``ZAllocation.dat`` file
    (``mod_tools.f90`` ``save_allocation``): branches alive, wind
    amplitude, new seeds dropped, branches pruned, reserve normalised by
    ``volume_twig``. All five overlaid on a semilog y-axis.

    Parameters
    ----------
    stats_history:
        Iterable of :class:`mechatree.simulate.TreeStats` objects.
    volume_twig:
        Base twig volume — ``reserve`` is divided by this to make the
        units match the Fortran plot. Defaults to 1 (raw reserve).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go

    from mechatree.plotting import figstyle

    stats = list(stats_history)
    gens = [s.generation for s in stats]

    # Semilog y can't take zero — clamp at 0.5 so quiet steps stay visible.
    # (Fortran does the same trick.)
    def safe(series):
        return [max(0.5, v) for v in series]

    fig = figstyle.figure(size="full", aspect=9 / 5)
    fig.add_trace(
        go.Scatter(
            x=gens,
            y=safe([s.n_branches for s in stats]),
            mode="lines",
            line=dict(color=figstyle.COLORS["black"]),  # noqa: C408
            name="N branches",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=gens,
            y=safe([s.wind_amplitude for s in stats]),
            mode="lines",
            line=dict(color=figstyle.COLORS["red"]),  # noqa: C408
            name="wind amplitude",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=gens,
            y=safe([s.n_seeds for s in stats]),
            mode="lines",
            line=dict(color=figstyle.COLORS["blue"]),  # noqa: C408
            name="new seeds",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=gens,
            y=safe([s.n_pruned for s in stats]),
            mode="lines+markers",
            line=dict(color=figstyle.COLORS["grey"]),  # noqa: C408
            marker=dict(size=4),  # noqa: C408
            name="N pruned",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=gens,
            y=safe([s.reserve / volume_twig for s in stats]),
            mode="lines",
            line=dict(color=figstyle.COLORS["black"], dash="dash"),  # noqa: C408
            name="reserve / V_twig",
        )
    )

    fig.update_layout(
        title="Per-step allocation",
        xaxis_title="generation",
        yaxis_title="count / amplitude / multiples of V_twig",
        yaxis_type="log",
        margin=dict(l=60, r=20, t=60, b=50),  # noqa: C408
    )
    return fig


def plot_strahler_diagnostics(tree: PyTree):
    """Strahler-order self-similarity panel for one tree.

    Four panels:

    1. Branch count per Strahler order (semilog).
    2. Mean branch length per order (semilog).
    3. Mean cross-sectional area per order (semilog).
    4. Histogram of ``(A_L + A_R) / A_P`` at each binary junction
       (Leonardo's-rule check), with a red dashed line at ratio 1.

    For a perfectly self-similar Horton-Strahler tree the first three
    are straight lines.

    Parameters
    ----------
    tree:
        A :class:`mechatree._core.PyTree`. Strahler order is set as a
        side-effect of running the diagnostics.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import plotly.graph_objects as go

    from mechatree.plotting import figstyle
    from mechatree.stats import leonardo_ratios, strahler_summary

    summary = strahler_summary(tree)
    ratios = leonardo_ratios(tree)

    title4 = f"Area-preservation at {ratios.size} junctions" if ratios.size else "Area-preservation"
    fig = figstyle.subplots(
        size="full",
        aspect=11 / 8,
        rows=2,
        cols=2,
        subplot_titles=(
            "Branch count per order",
            "Mean length per order",
            "Mean area per order",
            title4,
        ),
    )

    orders = list(range(1, summary.max_order + 1))
    marker = dict(color=figstyle.COLORS["black"], size=6)  # noqa: C408

    fig.add_trace(
        go.Scatter(x=orders, y=summary.n_branches, mode="markers", marker=marker, showlegend=False),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=orders, y=summary.mean_length, mode="markers", marker=marker, showlegend=False
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(x=orders, y=summary.mean_area, mode="markers", marker=marker, showlegend=False),
        row=2,
        col=1,
    )

    if ratios.size:
        fig.add_trace(
            go.Histogram(
                x=ratios,
                nbinsx=20,
                marker=dict(  # noqa: C408
                    color=figstyle.COLORS["green"],
                    line=dict(color=figstyle.COLORS["black"], width=1),  # noqa: C408
                ),
                showlegend=False,
            ),
            row=2,
            col=2,
        )
        fig.add_vline(
            x=1.0,
            line=dict(color=figstyle.COLORS["red"], dash="dash"),  # noqa: C408
            annotation_text="Leonardo (=1)",
            row=2,
            col=2,
        )
    else:
        fig.add_annotation(
            text="no binary junctions",
            xref="x4",
            yref="y4",
            x=0.5,
            y=0.5,
            showarrow=False,
        )

    fig.update_xaxes(title_text="Strahler order", row=1, col=1)
    fig.update_yaxes(title_text="N branches", type="log", row=1, col=1)
    fig.update_xaxes(title_text="Strahler order", row=1, col=2)
    fig.update_yaxes(title_text="mean length", type="log", row=1, col=2)
    fig.update_xaxes(title_text="Strahler order", row=2, col=1)
    fig.update_yaxes(title_text="mean cross-section area", type="log", row=2, col=1)
    fig.update_xaxes(title_text="(A_L + A_R) / A_P", row=2, col=2)
    fig.update_yaxes(title_text="count", row=2, col=2)

    fig.update_layout(margin=dict(l=60, r=20, t=60, b=50))  # noqa: C408
    return fig


def _unpack_thinning(history: Iterable) -> tuple[list[int], list[int], list[float]]:
    """Accept either ForestStats or (gen, n, m) tuples."""
    gens: list[int] = []
    ns: list[int] = []
    ms: list[float] = []
    for entry in history:
        if isinstance(entry, tuple):
            g, n, m = entry
        else:
            g = entry.generation
            n = entry.n_trees
            m = entry.biomass_total
        gens.append(int(g))
        ns.append(int(n))
        ms.append(float(m))
    return gens, ns, ms


__all__ = ["plot_allocation", "plot_self_thinning", "plot_strahler_diagnostics"]
