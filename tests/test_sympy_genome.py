"""Tests for the SymPy-callable genome (Step 15)."""

from __future__ import annotations

import math

import pytest

from mechatree.config import Config, GenomeConfig, TreeConfig
from mechatree.genome import (
    AllocationModel,
    CallbackAllocation,
    CallbackSafety,
    ConstantAllocation,
    ConstantSafety,
    SafetyModel,
    models_from_config,
)
from mechatree.simulate import grow_tree

pytest.importorskip("sympy")

# Imports below need SymPy; safe to do after the skip.
from mechatree.sympy_genome import (  # noqa: E402
    ALLOCATION_SYMBOLS,
    SAFETY_SYMBOLS,
    _compile,
    sympy_allocation,
    sympy_safety,
)

# ---------------------------------------------------------------------------
# 1. Callback-class smoke (no SymPy involved, exercises the C++ vtable shim)
# ---------------------------------------------------------------------------


def test_pycallbacksafety_dispatches_to_python():
    calls = []

    def fn(nb, ms):
        calls.append((nb, ms))
        return 2.5

    s = CallbackSafety(fn)
    assert s.compute(7, 0.4) == pytest.approx(2.5)
    assert calls == [(7, 0.4)]


def test_pycallbackallocation_dispatches_to_python():
    s = CallbackAllocation(lambda nb, vr: (0.1, 0.5, 0.5))
    assert s.compute(0, 0.0) == (0.1, 0.5, 0.5)


def test_pycallbacksafety_callable_kept_alive_inside_grow_tree():
    # If the C++ side held only a borrowed reference and we didn't keep the
    # Python callable alive, this would crash. The wrapper stashes _py_callable
    # explicitly to guard against that.
    def fn(nb, ms):
        return 3.0 - 0.5 * ms

    s = CallbackSafety(fn)
    a = ConstantAllocation(p_seeds=0.1, p_leaves=0.5, phototropism=0.5)
    tree = grow_tree(TreeConfig(), n_generations=5, seed=0, safety=s, allocation=a)
    assert tree.get_number_of_branches() >= 1


def test_pycallbacksafety_rejects_non_callable():
    with pytest.raises(TypeError, match="callable"):
        CallbackSafety(42)


def test_pycallbackallocation_returning_wrong_length_falls_back_to_zeros(capsys):
    a = CallbackAllocation(lambda nb, vr: (0.1, 0.5))  # length 2!
    p_seeds, p_leaves, phototropism = a.compute(0, 0.0)
    assert (p_seeds, p_leaves, phototropism) == (0.0, 0.0, 0.0)
    captured = capsys.readouterr()
    assert "expected 3" in captured.err


def test_pycallbacksafety_raising_falls_back_to_zero(capsys):
    def boom(nb, ms):
        raise RuntimeError("nope")

    s = CallbackSafety(boom)
    assert s.compute(0, 0.0) == 0.0
    captured = capsys.readouterr()
    assert "nope" in captured.err


# ---------------------------------------------------------------------------
# 2. sympy_safety / sympy_allocation
# ---------------------------------------------------------------------------


def test_sympy_safety_basic():
    s = sympy_safety("3 * tanh(max_stress)")
    assert s.compute(5, 0.0) == pytest.approx(0.0)
    assert s.compute(5, 100.0) == pytest.approx(3.0, abs=1e-9)
    assert s.compute(5, 1.0) == pytest.approx(3 * math.tanh(1.0))


def test_sympy_safety_uses_both_inputs():
    s = sympy_safety("nb_leaves * 0.01 + max_stress")
    assert s.compute(100, 0.5) == pytest.approx(1.5)


def test_sympy_safety_rejects_unknown_symbols():
    with pytest.raises(ValueError, match="unknown symbol"):
        sympy_safety("3 * tanh(x)")


def test_sympy_safety_rejects_empty_string():
    with pytest.raises(ValueError):
        sympy_safety("")


def test_sympy_allocation_all_strings():
    a = sympy_allocation(
        p_seeds="0.1 * tanh(vol_relative)",
        p_leaves="0.5",
        phototropism="0.5 + 0.1 * tanh(nb_leaves / 100)",
    )
    seeds, leaves, photo = a.compute(0, 0.0)
    assert seeds == pytest.approx(0.0)
    assert leaves == pytest.approx(0.5)
    assert photo == pytest.approx(0.5)

    seeds, leaves, photo = a.compute(100, 1.0)
    assert seeds == pytest.approx(0.1 * math.tanh(1.0))
    assert leaves == pytest.approx(0.5)
    assert photo == pytest.approx(0.5 + 0.1 * math.tanh(1.0))


def test_sympy_allocation_mixes_scalars_and_strings():
    a = sympy_allocation(p_seeds=0.1, p_leaves=0.5, phototropism="0.5 + 0 * vol_relative")
    assert a.compute(0, 0.0) == (pytest.approx(0.1), pytest.approx(0.5), pytest.approx(0.5))


def test_sympy_allocation_rejects_bad_symbol():
    with pytest.raises(ValueError, match="unknown symbol"):
        sympy_allocation(p_seeds="0.1 * max_stress", p_leaves=0.5, phototropism=0.5)


def test_sympy_allocation_rejects_non_string_non_numeric():
    with pytest.raises(TypeError):
        sympy_allocation(p_seeds=None, p_leaves=0.5, phototropism=0.5)


# ---------------------------------------------------------------------------
# 3. End-to-end via grow_tree
# ---------------------------------------------------------------------------


def test_sympy_genome_runs_grow_tree():
    s = sympy_safety("3 * tanh(max_stress + 0.1)")
    a = sympy_allocation(
        p_seeds="0.1 * tanh(vol_relative)",
        p_leaves=0.5,
        phototropism="0.5",
    )
    tree = grow_tree(TreeConfig(), n_generations=15, seed=42, safety=s, allocation=a)
    assert tree.get_number_of_branches() >= 1


def test_sympy_genome_matches_constant_when_expressions_are_constant():
    """A constant SymPy expression should produce the same trajectory as
    ConstantSafety/ConstantAllocation. Guards against the callback path
    silently perturbing the simulation."""
    s_const = ConstantSafety(3.0)
    a_const = ConstantAllocation(p_seeds=0.1, p_leaves=0.5, phototropism=0.5)
    t_const = grow_tree(
        TreeConfig(),
        n_generations=20,
        seed=7,
        safety=s_const,
        allocation=a_const,
    )

    s_sym = sympy_safety("3.0")
    a_sym = sympy_allocation(p_seeds="0.1", p_leaves="0.5", phototropism="0.5")
    t_sym = grow_tree(
        TreeConfig(),
        n_generations=20,
        seed=7,
        safety=s_sym,
        allocation=a_sym,
    )

    assert t_const.get_number_of_branches() == t_sym.get_number_of_branches()


# ---------------------------------------------------------------------------
# 4. GenomeConfig / YAML wiring
# ---------------------------------------------------------------------------


def test_genome_config_accepts_string_fields():
    gc = GenomeConfig(safety="3 * tanh(max_stress)", p_seeds=0.1, p_leaves=0.5, phototropism=0.5)
    assert gc.safety == "3 * tanh(max_stress)"


def test_genome_config_rejects_empty_string_field():
    with pytest.raises(ValueError, match="empty string"):
        GenomeConfig(safety="   ")


def test_models_from_config_dispatches_to_sympy_when_any_expression():
    gc = GenomeConfig(safety="3.0", p_seeds=0.1, p_leaves=0.5, phototropism=0.5)
    safety, allocation = models_from_config(gc)
    assert isinstance(safety, SafetyModel)
    assert isinstance(allocation, AllocationModel)
    # Concretely, the SymPy path produces Callback* (not Constant*).
    assert isinstance(safety, CallbackSafety)
    assert isinstance(allocation, CallbackAllocation)


def test_models_from_config_uses_constants_when_all_numeric():
    gc = GenomeConfig(safety=3.0, p_seeds=0.1, p_leaves=0.5, phototropism=0.5)
    safety, allocation = models_from_config(gc)
    assert isinstance(safety, ConstantSafety)
    assert isinstance(allocation, ConstantAllocation)


def test_yaml_string_expression_grows_a_tree(tmp_path):
    yaml_text = """
tree:
  twig_length: 1.0
  twig_diameter: 0.1
genome:
  safety: "3 * tanh(max_stress + 0.1)"
  p_seeds: 0.1
  p_leaves: "0.5"
  phototropism: 0.5
n_generations: 10
"""
    path = tmp_path / "sympy.yaml"
    path.write_text(yaml_text)
    cfg = Config.from_yaml(path)
    assert isinstance(cfg.genome.safety, str)

    tree = grow_tree(cfg, n_generations=5, seed=42)
    assert tree.get_number_of_branches() >= 1


# ---------------------------------------------------------------------------
# 5. _compile internals
# ---------------------------------------------------------------------------


def test_compile_returns_a_callable_with_correct_arity():
    fn_safety = _compile("3 * tanh(max_stress)", SAFETY_SYMBOLS)
    assert callable(fn_safety)
    assert fn_safety(0, 0.0) == pytest.approx(0.0)

    fn_alloc = _compile("0.1 * vol_relative + nb_leaves * 0.0", ALLOCATION_SYMBOLS)
    assert fn_alloc(0, 1.0) == pytest.approx(0.1)


def test_compile_rejects_non_string():
    with pytest.raises(TypeError):
        _compile(3.0, SAFETY_SYMBOLS)
