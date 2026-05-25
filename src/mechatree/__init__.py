"""MechaTree — wind- and light-driven tree growth simulator.

The recommended entry point is::

    import mechatree as mt

    cfg = mt.load_config("forest.yaml")
    tree = mt.grow_tree(cfg, n_generations=400, seed=0)
    fig = mt.plot_tree_3d(tree, style="cylinders")

Every symbol exposed here also lives in a focused submodule
(``mechatree.simulate.grow_tree``, ``mechatree.plotting.plot_tree_3d``, …);
both spellings stay working so existing scripts don't break.
"""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from mechatree import evolution
from mechatree._core import PyTree
from mechatree.config import (
    Config,
    ForestConfig,
    GenomeConfig,
    LightConfig,
    TreeConfig,
    WindConfig,
    load_config,
)
from mechatree.evolution import ForestEvolutionResult, Genome, run_tournament
from mechatree.forest import Forest, ForestStats
from mechatree.genome import (
    AllocationModel,
    CallbackAllocation,
    CallbackSafety,
    ConstantAllocation,
    ConstantSafety,
    NeuralAllocation,
    NeuralSafety,
    SafetyModel,
    champion_angles,
    load_all_champions,
    load_champion,
    models_from_config,
)
from mechatree.light import Leaves, Sun
from mechatree.plotting import (
    figstyle,
    plot_2d,
    plot_3d,
    plot_allocation,
    plot_forest_topdown,
    plot_self_thinning,
    plot_strahler_diagnostics,
    plot_tree_3d,
)
from mechatree.simulate import TreeStats, grow_tree
from mechatree.stats import (
    HortonRatios,
    HortonSummary,
    StrahlerSummary,
    distance_to_leaves,
    horton_ratios,
    horton_strahler_counts,
    horton_summary,
    leonardo_ratios,
    mean_distance_to_leaves,
    mean_stream_length,
    strahler_summary,
    tokunaga_matrix,
)

__version__ = "0.0.0.dev0"

# Optional DendroFlow bridge surface — lazy so a bare install never tries to
# import DendroFlow. ``mt.BranchWindBridge`` works iff
# ``pip install 'mechatree[dendroflow]'`` (or a sibling editable install) made
# the ``dendroflow`` package importable.
_DENDROFLOW_NAMES = frozenset(
    {
        "BranchWindBridge",
        "DendroFlowWindParams",
        "make_dendroflow_wind_fn",
    }
)

if TYPE_CHECKING:  # type-checker view only; runtime uses __getattr__
    from mechatree.wind.dendroflow import (  # noqa: F401
        BranchWindBridge,
        DendroFlowWindParams,
        make_dendroflow_wind_fn,
    )


def __getattr__(name: str):
    if name in _DENDROFLOW_NAMES:
        if find_spec("dendroflow") is None:
            raise AttributeError(
                f"mechatree.{name} requires the 'dendroflow' extra; "
                "install with `pip install 'mechatree[dendroflow]'`."
            )
        from mechatree.wind import dendroflow as _df

        return getattr(_df, name)
    raise AttributeError(f"module 'mechatree' has no attribute {name!r}")


__all__ = [
    "AllocationModel",
    "BranchWindBridge",
    "CallbackAllocation",
    "CallbackSafety",
    "Config",
    "ConstantAllocation",
    "ConstantSafety",
    "DendroFlowWindParams",
    "Forest",
    "ForestConfig",
    "ForestEvolutionResult",
    "ForestStats",
    "Genome",
    "GenomeConfig",
    "HortonRatios",
    "HortonSummary",
    "Leaves",
    "LightConfig",
    "NeuralAllocation",
    "NeuralSafety",
    "PyTree",
    "SafetyModel",
    "StrahlerSummary",
    "Sun",
    "TreeConfig",
    "TreeStats",
    "WindConfig",
    "__version__",
    "champion_angles",
    "distance_to_leaves",
    "evolution",
    "figstyle",
    "grow_tree",
    "horton_ratios",
    "horton_strahler_counts",
    "horton_summary",
    "leonardo_ratios",
    "load_all_champions",
    "load_champion",
    "load_config",
    "make_dendroflow_wind_fn",
    "mean_distance_to_leaves",
    "mean_stream_length",
    "models_from_config",
    "plot_2d",
    "plot_3d",
    "plot_allocation",
    "plot_forest_topdown",
    "plot_self_thinning",
    "plot_strahler_diagnostics",
    "plot_tree_3d",
    "run_tournament",
    "strahler_summary",
    "tokunaga_matrix",
]
