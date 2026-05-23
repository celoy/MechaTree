import matplotlib.pyplot as plt
import pytest
from matplotlib.figure import Figure

from mechatree import PyTree
from mechatree.plotting import plot_2d, plot_3d


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
