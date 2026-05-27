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
    plot_storm_replay,
    plot_strahler_diagnostics,
    plot_tree_3d,
)
from mechatree.simulate import TreeStats, default_wind_fn, grow_tree, make_default_wind_fn
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
from mechatree.wind.distributions import (
    Distribution,
    default_amplitude_sampler,
    uniform_angle_sampler,
)
from mechatree.wind.replay import StormPreSnapshot, StormSnapshot, run_storm_replay

__version__ = "0.0.0.dev0"


__all__ = [
    "AllocationModel",
    "CallbackAllocation",
    "CallbackSafety",
    "Config",
    "ConstantAllocation",
    "ConstantSafety",
    "Distribution",
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
    "StormPreSnapshot",
    "StormSnapshot",
    "StrahlerSummary",
    "Sun",
    "TreeConfig",
    "TreeStats",
    "WindConfig",
    "__version__",
    "champion_angles",
    "default_amplitude_sampler",
    "default_wind_fn",
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
    "make_default_wind_fn",
    "mean_distance_to_leaves",
    "mean_stream_length",
    "models_from_config",
    "plot_2d",
    "plot_3d",
    "plot_allocation",
    "plot_forest_topdown",
    "plot_self_thinning",
    "plot_storm_replay",
    "plot_strahler_diagnostics",
    "plot_tree_3d",
    "run_storm_replay",
    "run_tournament",
    "strahler_summary",
    "tokunaga_matrix",
    "uniform_angle_sampler",
]
