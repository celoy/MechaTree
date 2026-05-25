"""Pin the flat top-level public API exposed by `import mechatree as mt`.

Step 23 — the goal is that notebooks and recipe scripts can call every
genuinely public symbol through the top-level namespace, without
remembering which submodule each lives in. Private helpers
(`_callback_arity`, `_resolve_wind_fn`, `_decode_angles`, `_core` types, …)
must stay buried.

Deep imports (`from mechatree.simulate import grow_tree`) remain working
for back-compat.
"""

from __future__ import annotations

from importlib.util import find_spec

import mechatree as mt
from mechatree.config import load_config as deep_load_config
from mechatree.simulate import grow_tree as deep_grow_tree

# Names that are part of ``mt.__all__`` but only resolve when the
# ``dendroflow`` extra is installed (gated via the lazy ``__getattr__`` in
# :mod:`mechatree`). On a bare install these raise ``AttributeError`` —
# that *is* the documented behaviour — so the public-API test below skips
# them when the extra is unavailable.
_DENDROFLOW_OPTIONAL = frozenset(
    {"BranchWindBridge", "DendroFlowWindParams", "make_dendroflow_wind_fn"}
)
_HAS_DENDROFLOW = find_spec("dendroflow") is not None


def test_all_public_names_resolve():
    """Every name in `mt.__all__` must be importable as `mt.<name>`,
    except the DendroFlow extras when the optional package is missing."""
    missing = [
        name
        for name in mt.__all__
        if not hasattr(mt, name) and not (name in _DENDROFLOW_OPTIONAL and not _HAS_DENDROFLOW)
    ]
    assert not missing, f"missing flat-API names: {missing}"


def test_public_api_contains_expected_symbols():
    expected = {
        # config
        "Config",
        "TreeConfig",
        "LightConfig",
        "ForestConfig",
        "GenomeConfig",
        "load_config",
        # simulate / forest
        "grow_tree",
        "Forest",
        "ForestStats",
        "TreeStats",
        # light
        "Sun",
        "Leaves",
        # genome
        "load_champion",
        "load_all_champions",
        "champion_angles",
        "models_from_config",
        "ConstantSafety",
        "ConstantAllocation",
        "NeuralSafety",
        "NeuralAllocation",
        "CallbackSafety",
        "CallbackAllocation",
        # plotting
        "plot_tree_3d",
        "plot_forest_topdown",
        "plot_2d",
        "plot_3d",
        "plot_self_thinning",
        "plot_allocation",
        "plot_strahler_diagnostics",
        "figstyle",
        # stats
        "horton_summary",
        "strahler_summary",
        "horton_ratios",
        "distance_to_leaves",
        "mean_distance_to_leaves",
        "leonardo_ratios",
        "tokunaga_matrix",
        # core type
        "PyTree",
    }
    missing = expected - set(mt.__all__)
    assert not missing, f"flat API missing required symbols: {sorted(missing)}"


def test_private_helpers_not_exposed():
    """Implementation details from the CLAUDE.md spec must not leak."""
    leaked = [
        name
        for name in (
            "_callback_arity",
            "_resolve_wind_fn",
            "_decode_angles",
            "_champion",
        )
        if hasattr(mt, name)
    ]
    assert not leaked, f"private helpers leaked through flat API: {leaked}"


def test_flat_api_end_to_end():
    """A minimal grow_tree run via flat API only — exercises load_config,
    Config defaults, grow_tree, and plot_tree_3d."""
    cfg = mt.Config()
    tree = mt.grow_tree(cfg, n_generations=5, seed=0)
    assert isinstance(tree, mt.PyTree)
    assert tree.get_number_of_branches() > 0
    fig = mt.plot_tree_3d(tree)
    # plotly is a hard dep so we know the import; we only care that the
    # figure object came out of the flat path.
    assert fig.data is not None


def test_figstyle_module_handle():
    """`mt.figstyle` should be the live submodule, not a copy."""
    from mechatree.plotting import figstyle as deep

    assert mt.figstyle is deep
    assert mt.figstyle.TEMPLATE_NAME == "mechatree"


def test_deep_imports_still_work():
    """Existing `from mechatree.simulate import grow_tree` callers must
    keep working — and resolve to the same function as `mt.grow_tree`."""
    assert mt.grow_tree is deep_grow_tree
    assert mt.load_config is deep_load_config


def test_version_string_exposed():
    assert isinstance(mt.__version__, str)
    assert mt.__version__  # non-empty
