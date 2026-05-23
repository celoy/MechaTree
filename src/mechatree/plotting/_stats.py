"""Reusable plotting helpers — Python ports of the MATLAB scripts in
``../Eloy2017_NatComm_archive/``.

Each helper takes pre-computed simulation data (a tree, a forest, or a
history of stats objects) and returns the matplotlib ``Figure`` it drew
on. The example scripts in ``examples/`` are thin CLI wrappers that
drive a simulation and call these.

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

import matplotlib.pyplot as plt

from mechatree._core import PyTree


def plot_self_thinning(history, ax=None):
    """Population, biomass, and N-vs-M log-log on a 1×3 grid.

    Parameters
    ----------
    history:
        Iterable of ``(generation, n_trees, biomass)`` tuples, or
        :class:`mechatree.forest.ForestStats` objects.
    ax:
        Optional sequence of three matplotlib axes. A new figure is made
        if omitted.

    Returns
    -------
    matplotlib.figure.Figure
    """
    gens, ns, ms = _unpack_thinning(history)

    if ax is None:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    else:
        axes = list(ax)
        if len(axes) != 3:
            raise ValueError("plot_self_thinning needs exactly 3 axes")
        fig = axes[0].figure

    axes[0].plot(gens, ns, color="forestgreen")
    axes[0].set_xlabel("generation")
    axes[0].set_ylabel("N (trees alive)")
    axes[0].set_title("Population over time")

    axes[1].plot(gens, ms, color="saddlebrown")
    axes[1].set_xlabel("generation")
    axes[1].set_ylabel("M (total biomass)")
    axes[1].set_title("Biomass over time")

    # Pairs where N > 0 — log axes can't take zero.
    pairs = [(n, m) for n, m in zip(ns, ms, strict=True) if n > 0 and m > 0]
    if pairs:
        ax_log = axes[2]
        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]
        ax_log.loglog(xs, ys, "o-", color="black", markersize=3)
        # Reference line: M ~ N^{-3/2} — Yoda 1963's self-thinning law.
        n_anchor = min(xs)
        m_anchor = max(ys)
        ref_n = [min(xs), max(xs)]
        ref_m = [m_anchor * (n_anchor / n) ** 1.5 for n in ref_n]
        ax_log.loglog(ref_n, ref_m, "--r", label=r"$M \propto N^{-3/2}$")
        ax_log.legend()
        ax_log.set_xlabel("N (trees alive)")
        ax_log.set_ylabel("M (total biomass)")
        ax_log.set_title("Self-thinning (log-log)")
    else:
        axes[2].text(0.5, 0.5, "no data", ha="center", va="center")
        axes[2].set_axis_off()

    fig.tight_layout()
    return fig


def plot_allocation(stats_history, volume_twig: float = 1.0, ax=None):
    """Per-generation allocation diagnostics for a single tree.

    Reproduces the columns of the Fortran ``ZAllocation.dat`` file
    (``mod_tools.f90`` ``save_allocation``): branches alive, wind
    amplitude, new seeds dropped, branches pruned, reserve normalised
    by ``volume_twig``. All five on one ``semilogy`` plot.

    Parameters
    ----------
    stats_history:
        Iterable of :class:`mechatree.simulate.TreeStats` objects (one
        per generation, in order).
    volume_twig:
        Base twig volume — ``reserve`` is divided by this to make the
        units match the Fortran plot. Defaults to 1 (raw reserve).
    ax:
        Optional axes; a new figure is made if omitted.

    Returns
    -------
    matplotlib.figure.Figure
    """
    stats = list(stats_history)
    gens = [s.generation for s in stats]
    n_branches = [s.n_branches for s in stats]
    wind_amp = [s.wind_amplitude for s in stats]
    n_seeds = [s.n_seeds for s in stats]
    n_pruned = [s.n_pruned for s in stats]
    reserve_units = [s.reserve / volume_twig for s in stats]

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig = ax.figure

    # Semilog can't show zero — clamp at 0.5 so quiet steps stay visible.
    # (Fortran does the same trick.)
    def safe(series):
        return [max(0.5, v) for v in series]

    ax.semilogy(gens, safe(n_branches), "k-", label="N branches")
    ax.semilogy(gens, safe(wind_amp), "r-", label="wind amplitude")
    ax.semilogy(gens, safe(n_seeds), "m-", label="new seeds")
    ax.semilogy(gens, safe(n_pruned), "-co", markersize=3, label="N pruned")
    ax.semilogy(gens, safe(reserve_units), "--k", label="reserve / V_twig")
    ax.set_xlabel("generation")
    ax.set_ylabel("count / amplitude / multiples of V_twig")
    ax.set_title("Per-step allocation")
    ax.legend(loc="best")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    return fig


def plot_strahler_diagnostics(tree: PyTree, ax=None):
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
    ax:
        Optional 2×2 array (or 4-element sequence) of axes.

    Returns
    -------
    matplotlib.figure.Figure
    """
    # Lazy import to avoid pulling numpy at module load.
    from mechatree.stats import leonardo_ratios, strahler_summary

    summary = strahler_summary(tree)
    ratios = leonardo_ratios(tree)

    if ax is None:
        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        axes = axes.ravel()
    else:
        axes = list(ax)
        if len(axes) != 4:
            raise ValueError("plot_strahler_diagnostics needs exactly 4 axes")
        fig = axes[0].figure

    orders = list(range(1, summary.max_order + 1))

    axes[0].semilogy(orders, summary.n_branches, "ok")
    axes[0].set_xlabel("Strahler order")
    axes[0].set_ylabel("N branches")
    axes[0].set_title("Branch count per order")
    axes[0].grid(True, alpha=0.3)

    axes[1].semilogy(orders, summary.mean_length, "ok")
    axes[1].set_xlabel("Strahler order")
    axes[1].set_ylabel("mean length")
    axes[1].set_title("Mean length per order")
    axes[1].grid(True, alpha=0.3)

    axes[2].semilogy(orders, summary.mean_area, "ok")
    axes[2].set_xlabel("Strahler order")
    axes[2].set_ylabel("mean cross-section area")
    axes[2].set_title("Mean area per order")
    axes[2].grid(True, alpha=0.3)

    if ratios.size:
        axes[3].hist(ratios, bins=20, color="forestgreen", edgecolor="black")
        axes[3].axvline(1.0, color="red", linestyle="--", label="Leonardo (=1)")
        axes[3].set_xlabel(r"$(A_\mathrm{L} + A_\mathrm{R}) / A_\mathrm{P}$")
        axes[3].set_ylabel("count")
        axes[3].set_title(f"Area-preservation at {ratios.size} junctions")
        axes[3].legend()
    else:
        axes[3].text(0.5, 0.5, "no binary junctions", ha="center", va="center")
        axes[3].set_axis_off()

    fig.tight_layout()
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
            # Duck-type ForestStats — must have `.generation`, `.n_trees`,
            # `.biomass_total`.
            g = entry.generation
            n = entry.n_trees
            m = entry.biomass_total
        gens.append(int(g))
        ns.append(int(n))
        ms.append(float(m))
    return gens, ns, ms


__all__ = ["plot_allocation", "plot_self_thinning", "plot_strahler_diagnostics"]
