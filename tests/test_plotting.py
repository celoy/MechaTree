import matplotlib.pyplot as plt
import pytest
from matplotlib.figure import Figure

from mechatree import PyTree
from mechatree.config import Config, ForestConfig
from mechatree.forest import Forest
from mechatree.plotting import (
    plot_2d,
    plot_3d,
    plot_allocation,
    plot_forest_topdown,
    plot_self_thinning,
    plot_strahler_diagnostics,
    plot_tree_3d,
)
from mechatree.simulate import TreeStats, grow_tree


@pytest.fixture
def small_2d_tree():
    tree = PyTree({"x": 0.0, "y": 0.0, "theta": 1e-6, "L": 1.0, "grow": 1})
    tree.add_branch(0, {"x": 0.0, "y": 1.0, "theta": 0.3, "L": 0.5, "grow": 1})
    tree.add_branch(0, {"x": 0.0, "y": 1.0, "theta": -0.3, "L": 0.5, "grow": 1})
    return tree


@pytest.fixture
def small_3d_tree():
    tree = PyTree({"x": 0.0, "y": 0.0, "z": 0.0, "theta": 0.0, "phi": 0.0})
    tree.add_branch(0, {"x": 0.0, "y": 0.0, "z": 1.0, "theta": 0.3, "phi": 0.0})
    tree.add_branch(0, {"x": 0.0, "y": 0.0, "z": 1.0, "theta": 0.3, "phi": 3.14})
    return tree


def test_plot_2d_returns_figure(small_2d_tree):
    fig = plot_2d(small_2d_tree)
    assert isinstance(fig, Figure)
    plt.close(fig)


def test_plot_2d_saves_png(small_2d_tree, tmp_path):
    fig = plot_2d(small_2d_tree, iteration=1, out_dir=tmp_path)
    out_file = tmp_path / "fig0001.png"
    assert out_file.exists()
    assert out_file.stat().st_size > 0
    plt.close(fig)


def test_plot_2d_requires_iteration_when_saving(small_2d_tree, tmp_path):
    with pytest.raises(ValueError):
        plot_2d(small_2d_tree, out_dir=tmp_path)


def test_plot_2d_accepts_user_axes(small_2d_tree):
    fig, ax = plt.subplots()
    returned = plot_2d(small_2d_tree, ax=ax)
    assert returned is fig
    plt.close(fig)


def test_plot_3d_returns_figure(small_3d_tree):
    fig = plot_3d(small_3d_tree)
    assert isinstance(fig, Figure)
    # 3D axes were created.
    assert len(fig.axes) == 1
    plt.close(fig)


def test_plot_3d_accepts_user_axes(small_3d_tree):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    returned = plot_3d(small_3d_tree, ax=ax)
    assert returned is fig
    plt.close(fig)


# --- Step 13: typed-field mechanics plots --------------------------------


@pytest.fixture
def grown_tree():
    """A short-grown mechanics tree (uses typed location/unit_t/length fields)."""
    return grow_tree(Config(), n_generations=8, seed=0)


def test_plot_tree_3d_returns_figure(grown_tree):
    fig = plot_tree_3d(grown_tree)
    assert isinstance(fig, Figure)
    assert len(fig.axes) == 1
    plt.close(fig)


def test_plot_tree_3d_accepts_user_axes(grown_tree):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    returned = plot_tree_3d(grown_tree, ax=ax)
    assert returned is fig
    plt.close(fig)


def test_plot_forest_topdown_returns_figure():
    cfg = Config(forest=ForestConfig(size=10.0, n_trees_init=5, n_trees_max=50))
    forest = Forest(cfg, seed=0)
    fig = plot_forest_topdown(forest)
    assert isinstance(fig, Figure)
    plt.close(fig)


def test_plot_forest_topdown_handles_empty_forest():
    """A forest with no trees alive still renders a boundary disk."""
    cfg = Config(forest=ForestConfig(size=10.0, n_trees_init=2, n_trees_max=50))
    forest = Forest(cfg, seed=0)
    forest.trees.clear()
    forest.ages.clear()
    fig = plot_forest_topdown(forest)
    assert isinstance(fig, Figure)
    plt.close(fig)


# --- MATLAB plot ports: stats helpers --------------------------------------


def test_plot_self_thinning_with_tuples():
    history = [(g, 10 + g, 0.5 * (1 + g)) for g in range(20)]
    fig = plot_self_thinning(history)
    assert isinstance(fig, Figure)
    # 1x3 grid.
    assert len(fig.axes) == 3
    plt.close(fig)


def test_plot_self_thinning_with_forest_stats():
    """The helper accepts ForestStats objects via duck-typing."""
    cfg = Config(forest=ForestConfig(size=10.0, n_trees_init=4, n_trees_max=50))
    forest = Forest(cfg, seed=0)
    history = []

    def cb(_gen, _f, stats):
        history.append(stats)

    forest.run(5, on_step=cb)
    fig = plot_self_thinning(history)
    assert isinstance(fig, Figure)
    plt.close(fig)


def test_plot_self_thinning_handles_empty_history():
    """An empty history still produces a (degenerate but valid) figure."""
    fig = plot_self_thinning([])
    assert isinstance(fig, Figure)
    plt.close(fig)


def test_plot_allocation_with_tree_stats():
    history = []

    def cb(_gen, _tree, stats):
        history.append(stats)

    grow_tree(Config(), n_generations=10, seed=0, on_step=cb)
    assert len(history) == 10
    # All entries are TreeStats.
    assert isinstance(history[0], TreeStats)

    fig = plot_allocation(history, volume_twig=0.01)
    assert isinstance(fig, Figure)
    plt.close(fig)


def test_plot_strahler_diagnostics_grown_tree():
    tree = grow_tree(Config(), n_generations=20, seed=0)
    fig = plot_strahler_diagnostics(tree)
    assert isinstance(fig, Figure)
    # 2x2 grid.
    assert len(fig.axes) == 4
    plt.close(fig)


def test_plot_strahler_diagnostics_with_user_axes():
    tree = grow_tree(Config(), n_generations=20, seed=0)
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    returned = plot_strahler_diagnostics(tree, ax=axes.ravel())
    assert returned is fig
    plt.close(fig)
