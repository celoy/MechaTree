"""C++ NeuralSafety / NeuralAllocation vs the pure-Python reference.

The Python reference lives in scripts/strategies_single_tree.py — a faithful
port of the matlab `neural_branch` / `neural_reserve` that produces the
contour figures in the Eloy et al. (Nat Commun 2017) lineage. The C++ classes
must agree with it to the last bit at every input point we check; otherwise
running grow_tree with an evolved genome would diverge silently from the
matlab/Fortran reference.

The champion JSON written by the same script is the test fixture.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from mechatree.genome import NeuralAllocation, NeuralSafety, load_all_champions, load_champion

REPO_ROOT = Path(__file__).resolve().parents[1]
CHAMPION_JSON = REPO_ROOT / "data" / "S3_champions.json"


def _load_python_reference():
    """Import scripts/strategies_single_tree.py without poking sys.path."""
    spec = importlib.util.spec_from_file_location(
        "_strategies_ref", REPO_ROOT / "scripts" / "strategies_single_tree.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("_strategies_ref", mod)
    spec.loader.exec_module(mod)
    return mod


pytestmark = pytest.mark.skipif(
    not CHAMPION_JSON.exists(),
    reason=f"champion JSON not present at {CHAMPION_JSON}; run "
    "`uv run python scripts/strategies_single_tree.py --build` first.",
)


@pytest.fixture(scope="module")
def champions():
    return json.loads(CHAMPION_JSON.read_text())["species"]


@pytest.fixture(scope="module")
def reference():
    return _load_python_reference()


# A handful of (nb_leaves, second_input) grid points covering the corners +
# interior of the panel ranges in scripts/strategies_single_tree.py.
SAFETY_POINTS = [
    (0, 0.0),
    (0, 0.5),
    (50, 0.25),
    (100, 0.0),
    (100, 0.5),
    (37, 0.13),
]
ALLOC_POINTS = [
    (0, 0.0),
    (0, 2.0),
    (50, 1.0),
    (100, 0.0),
    (100, 2.0),
    (37, 0.5),
]


@pytest.mark.parametrize("species_id", [0, 1])
def test_neural_safety_matches_python_reference(champions, reference, species_id):
    species = next(s for s in champions if s["species_id"] == species_id)
    weights = np.array(species["nn_branch"])
    cpp = NeuralSafety(weights)
    for nb_leaves, max_stress in SAFETY_POINTS:
        ref = float(
            reference.neural_branch(np.asarray(float(nb_leaves)), np.asarray(max_stress), weights)
        )
        got = cpp.compute(nb_leaves, max_stress)
        assert got == pytest.approx(ref, abs=1e-12), (
            f"species {species_id} at ({nb_leaves}, {max_stress}): C++ {got} != Python {ref}"
        )


@pytest.mark.parametrize("species_id", [0, 1])
def test_neural_allocation_matches_python_reference(champions, reference, species_id):
    species = next(s for s in champions if s["species_id"] == species_id)
    weights = np.array(species["nn_reserve"])
    cpp = NeuralAllocation(weights)
    for nb_leaves, vol_relative in ALLOC_POINTS:
        ps_ref, pl_ref, ph_ref = reference.neural_reserve(
            np.asarray(float(nb_leaves)), np.asarray(vol_relative), weights
        )
        got_ps, got_pl, got_ph = cpp.compute(nb_leaves, vol_relative)
        assert got_ps == pytest.approx(float(ps_ref), abs=1e-12)
        assert got_pl == pytest.approx(float(pl_ref), abs=1e-12)
        assert got_ph == pytest.approx(float(ph_ref), abs=1e-12)


def test_neural_safety_rejects_wrong_length():
    with pytest.raises(ValueError, match="10 weights"):
        NeuralSafety([0.5] * 9)
    with pytest.raises(ValueError, match="10 weights"):
        NeuralSafety([0.5] * 11)


def test_neural_allocation_rejects_wrong_length():
    with pytest.raises(ValueError, match="18 weights"):
        NeuralAllocation([0.5] * 17)


def test_neural_safety_weights_roundtrip():
    w = np.linspace(0.1, 0.9, 10)
    ns = NeuralSafety(w)
    np.testing.assert_array_equal(ns.weights, w)


def test_neural_allocation_weights_roundtrip():
    w = np.linspace(0.05, 0.95, 18)
    na = NeuralAllocation(w)
    np.testing.assert_array_equal(na.weights, w)


def test_load_champion_returns_models_angles_and_non_coding():
    safety, allocation, angles, non_coding = load_champion(CHAMPION_JSON, species_id=0)
    assert isinstance(safety, NeuralSafety)
    assert isinstance(allocation, NeuralAllocation)
    assert non_coding["species_id"] == 0
    assert non_coding["dataset"] == "S3.dat"
    assert "champion_index" in non_coding
    # Angles dict ready for ``replace(cfg.tree, **angles)``.
    assert set(angles) == {"theta1", "theta2", "gamma1", "gamma2"}
    for v in angles.values():
        assert isinstance(v, float)
    # Fortran formula: theta1 = genome(1) * pi/2, theta2 = -genome(2) * pi/2,
    # gamma1 = gamma2 = genome(3) * 2*pi. Both gammas equal.
    assert angles["gamma1"] == angles["gamma2"]
    assert angles["theta1"] > 0
    assert angles["theta2"] < 0


def test_load_champion_unknown_species_raises():
    with pytest.raises(ValueError, match="species_id=99"):
        load_champion(CHAMPION_JSON, species_id=99)


def test_load_all_champions_returns_one_per_species():
    items = load_all_champions(CHAMPION_JSON)
    assert len(items) == 2
    ids = {non_coding["species_id"] for _, _, _, non_coding in items}
    assert ids == {0, 1}
    # Each entry exposes the champion's angles as the third element.
    for _, _, angles, _ in items:
        assert set(angles) == {"theta1", "theta2", "gamma1", "gamma2"}
