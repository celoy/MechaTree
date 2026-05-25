"""Tests for the DendroFlow wind bridge (Step 17 / DendroFlow M6).

Skipped wholesale when the optional ``dendroflow`` extra isn't installed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

pytest.importorskip("dendroflow")

# Imports below are gated on DendroFlow being present.
from mechatree._core import PyTree  # noqa: E402
from mechatree.config import Config, TreeConfig, WindConfig  # noqa: E402
from mechatree.forest import Forest  # noqa: E402
from mechatree.simulate import grow_tree, make_seed_tree  # noqa: E402
from mechatree.wind.dendroflow import (  # noqa: E402
    DendroFlowWindParams,
    forest_to_cylinders,
    make_dendroflow_wind_fn,
    pytree_to_cylinders,
)


def _small_wind_kwargs():
    # 10 z-layers, H=0.5, domain 0..5 m, uniform U_infty=5 m/s.
    z_centers = tuple(0.25 + 0.5 * i for i in range(10))
    return {
        "U_infty": (5.0,) * 10,
        "z_centers": z_centers,
        "H": 0.5,
        "C_D": 1.0,
    }


# ---------------------------------------------------------------------------
# 1. pytree_to_cylinders shape on a seed tree.
# ---------------------------------------------------------------------------


def test_pytree_to_cylinders_seed_tree():
    cfg = TreeConfig()
    tree = make_seed_tree(cfg)
    tree.reorder()

    cylinders = pytree_to_cylinders(tree, tree_id=1.0)

    assert cylinders.n == 1
    row = cylinders.frame.iloc[0]
    assert row["x_c"] == pytest.approx(0.0)
    assert row["y_c"] == pytest.approx(0.0)
    assert row["z"] == pytest.approx(0.0)
    assert row["nx"] == pytest.approx(0.0)
    assert row["ny"] == pytest.approx(0.0)
    assert row["nz"] == pytest.approx(1.0)
    assert row["D"] == pytest.approx(cfg.twig_diameter)
    assert row["L"] == pytest.approx(cfg.twig_length)
    assert float(row["tree_id"]) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# 2. Bridge contract: (U, 0, 0) and degenerate seed-tree behaviour.
# ---------------------------------------------------------------------------


def test_bridge_returns_three_tuple_with_zero_lateral():
    tree = make_seed_tree(TreeConfig())
    tree.reorder()

    bridge = make_dendroflow_wind_fn(**_small_wind_kwargs())
    rng = np.random.default_rng(0)
    wind = bridge(0, rng, tree)

    assert isinstance(wind, tuple) and len(wind) == 3
    assert wind[1] == 0.0
    assert wind[2] == 0.0
    assert math.isfinite(wind[0])
    # The seed tree contributes a sliver of area; with H=0.5 the thinning
    # factor is ~1 so the canopy-mean ≈ U_infty.
    assert wind[0] == pytest.approx(5.0, rel=1e-3)
    # And the bridge stashed the underlying result for inspection.
    assert bridge.last_result is not None
    assert bridge.last_result.n_branches == 1


# ---------------------------------------------------------------------------
# 3. YAML config wires DendroFlow wind into grow_tree.
# ---------------------------------------------------------------------------


def test_yaml_config_wires_dendroflow_into_grow_tree(tmp_path):
    yaml_text = """
tree:
  twig_length: 1.0
  twig_diameter: 0.1
wind:
  model: dendroflow
  U_infty: [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
  z_centers: [0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75]
  H: 0.5
  C_D: 1.0
n_generations: 3
"""
    path = tmp_path / "tiny.yaml"
    path.write_text(yaml_text)

    cfg = Config.from_yaml(path)
    assert cfg.wind.model == "dendroflow"

    tree = grow_tree(cfg, n_generations=3, seed=42)

    # Smoke check: the tree at least exists and didn't crash. After three
    # generations of primary growth the seedling should have produced
    # additional twigs (deterministic given the seed).
    assert tree.get_number_of_branches() >= 1


# ---------------------------------------------------------------------------
# 4. Forest pools all trees into a single Cylinders.compute() per generation.
# ---------------------------------------------------------------------------


def test_forest_pools_all_trees_into_single_cylinders_call(monkeypatch):
    from dendroflow.wind import BulkThinningBranchWindModel

    calls: list[dict] = []
    real_compute = BulkThinningBranchWindModel.compute

    def recording_compute(self, cylinders, *, wind_params):
        calls.append(
            {
                "n_rows": cylinders.n,
                "tree_ids": sorted(cylinders.frame["tree_id"].unique().tolist()),
            }
        )
        return real_compute(self, cylinders, wind_params=wind_params)

    monkeypatch.setattr(BulkThinningBranchWindModel, "compute", recording_compute)

    cfg = Config(
        wind=WindConfig(
            model="dendroflow",
            **_small_wind_kwargs(),
        ),
    )
    # Override forest count to keep the test fast and deterministic.
    cfg = Config(
        tree=cfg.tree,
        light=cfg.light,
        forest=type(cfg.forest)(
            size=cfg.forest.size,
            n_trees_init=5,
            n_trees_max=cfg.forest.n_trees_max,
            min_branches=cfg.forest.min_branches,
            min_age_for_undersize=cfg.forest.min_age_for_undersize,
            max_age=cfg.forest.max_age,
        ),
        genome=cfg.genome,
        wind=cfg.wind,
        n_generations=cfg.n_generations,
    )

    forest = Forest(cfg, seed=7)
    forest.step(0)

    assert len(calls) == 1, f"expected one compute() per generation, got {len(calls)}"
    assert calls[0]["tree_ids"] == [0.0, 1.0, 2.0, 3.0, 4.0]


# ---------------------------------------------------------------------------
# 5. Arity detection preserves 2-arg wind callables (backwards compatible).
# ---------------------------------------------------------------------------


def test_arity_detection_preserves_two_arg_wind_fn():
    calls_two = []

    def two_arg(gen, rng):
        calls_two.append(gen)
        return (1.0, 0.0, 0.0)

    grow_tree(TreeConfig(), n_generations=2, seed=1, wind_fn=two_arg)
    assert calls_two == [0, 1]

    received_contexts: list[object] = []

    def three_arg(gen, rng, context):
        received_contexts.append(context)
        return (1.0, 0.0, 0.0)

    grow_tree(TreeConfig(), n_generations=2, seed=1, wind_fn=three_arg)
    assert len(received_contexts) == 2
    assert all(isinstance(c, PyTree) for c in received_contexts)


def test_forest_to_cylinders_assigns_tree_ids_by_index():
    cfg = TreeConfig()
    trees = []
    for _ in range(3):
        t = make_seed_tree(cfg)
        t.reorder()
        trees.append(t)

    cylinders = forest_to_cylinders(trees)
    assert cylinders.n == 3
    assert sorted(cylinders.frame["tree_id"].unique().tolist()) == [0.0, 1.0, 2.0]


def test_wind_config_validation():
    # Defaults: model=='default', everything else ignored.
    wc = WindConfig()
    assert wc.model == "default"

    # Missing required fields for dendroflow.
    with pytest.raises(ValueError, match="required"):
        WindConfig(model="dendroflow")

    # Length mismatch.
    with pytest.raises(ValueError, match="same length"):
        WindConfig(
            model="dendroflow",
            U_infty=(5.0, 5.0, 5.0),
            z_centers=(0.25, 0.75),
        )

    # z_centers must be monotone.
    with pytest.raises(ValueError, match="monotone"):
        WindConfig(
            model="dendroflow",
            U_infty=(5.0, 5.0, 5.0),
            z_centers=(0.75, 0.25, 1.25),
        )

    # z_centers must cover z=0.
    with pytest.raises(ValueError, match="covered"):
        WindConfig(
            model="dendroflow",
            U_infty=(5.0, 5.0),
            z_centers=(2.0, 3.0),
            H=0.5,
        )


def test_dendroflow_wind_params_validation():
    # Happy path.
    DendroFlowWindParams(
        U_infty=np.array([5.0, 5.0]),
        z_centers=np.array([0.25, 0.75]),
    )

    # Shape mismatch.
    with pytest.raises(ValueError, match="share shape"):
        DendroFlowWindParams(
            U_infty=np.array([5.0, 5.0, 5.0]),
            z_centers=np.array([0.25, 0.75]),
        )


def test_branch_wind_bridge_handles_empty_forest():
    bridge = make_dendroflow_wind_fn(**_small_wind_kwargs())
    cfg = Config(
        wind=WindConfig(model="dendroflow", **_small_wind_kwargs()),
    )
    forest = Forest(cfg, seed=0)
    forest.trees = []  # simulate mass die-off
    rng = np.random.default_rng(0)
    wind = bridge(0, rng, forest)
    # Falls back to free-stream at the ground.
    assert wind == (5.0, 0.0, 0.0)
