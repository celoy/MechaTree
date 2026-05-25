"""Tests for the heritable :class:`Genome` (Step 21)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from mechatree.evolution.genome import ALLOCATION_LEN, SAFETY_LEN, Genome
from mechatree.genome import NeuralAllocation, NeuralSafety


def test_random_genome_shape_and_bounds():
    rng = np.random.default_rng(0)
    g = Genome.random(rng, lineage_id=7)
    assert len(g.safety_weights) == SAFETY_LEN
    assert len(g.allocation_weights) == ALLOCATION_LEN
    assert len(g.angle_genes) == 3
    assert g.lineage_id == 7
    for w in g.safety_weights:
        assert 0.0 <= w <= 1.0
    for w in g.allocation_weights:
        assert 0.0 <= w <= 1.0
    for a in g.angle_genes:
        assert 0.0 <= a <= 1.0


def test_genome_validates_lengths():
    with pytest.raises(ValueError):
        Genome(
            safety_weights=tuple([0.5] * (SAFETY_LEN - 1)),
            allocation_weights=tuple([0.5] * ALLOCATION_LEN),
            angle_genes=(0.5, 0.5, 0.5),
        )
    with pytest.raises(ValueError):
        Genome(
            safety_weights=tuple([0.5] * SAFETY_LEN),
            allocation_weights=tuple([0.5] * (ALLOCATION_LEN - 1)),
            angle_genes=(0.5, 0.5, 0.5),
        )
    with pytest.raises(ValueError):
        Genome(
            safety_weights=tuple([0.5] * SAFETY_LEN),
            allocation_weights=tuple([0.5] * ALLOCATION_LEN),
            angle_genes=(0.5, 0.5),  # type: ignore[arg-type]
        )


def test_mutate_is_deterministic_and_returns_new_object():
    parent = Genome.random(np.random.default_rng(0))
    child1 = parent.mutate(np.random.default_rng(123))
    child2 = parent.mutate(np.random.default_rng(123))
    assert child1.safety_weights == child2.safety_weights
    assert child1.allocation_weights == child2.allocation_weights
    assert child1.angle_genes == child2.angle_genes
    # Parent unchanged (frozen).
    assert (
        parent.safety_weights != child1.safety_weights or parent.angle_genes != child1.angle_genes
    )
    # Lineage inherited.
    assert child1.lineage_id == parent.lineage_id


def test_mutate_clips_weights_but_not_angles():
    """A locus at the boundary must stay in [0, 1] for weights and may go
    outside for angle genes (the decoder scales them anyway)."""
    g = Genome(
        safety_weights=tuple([1.0] * SAFETY_LEN),
        allocation_weights=tuple([0.0] * ALLOCATION_LEN),
        angle_genes=(0.0, 0.5, 1.0),
    )
    # Force every locus to mutate with large sigma so we definitely cross
    # a bound.
    rng = np.random.default_rng(0)
    child = g.mutate(rng, sigma=1.0, p_locus=1.0)
    for w in child.safety_weights:
        assert 0.0 <= w <= 1.0
    for w in child.allocation_weights:
        assert 0.0 <= w <= 1.0
    # At least one angle should have escaped the unit interval (huge sigma).
    any_outside = any(a < 0.0 or a > 1.0 for a in child.angle_genes)
    assert any_outside, child.angle_genes


def test_mutation_at_zero_sigma_or_rate_is_identity_modulo_forced_locus():
    """With sigma=0 the forced-mutation guard still fires once but the
    Gaussian noise is 0, so the result equals the parent."""
    g = Genome.random(np.random.default_rng(0))
    child = g.mutate(np.random.default_rng(0), sigma=0.0, p_locus=0.0)
    assert child.safety_weights == g.safety_weights
    assert child.allocation_weights == g.allocation_weights
    assert child.angle_genes == g.angle_genes


def test_to_models_materializes_neural_pair_and_caches():
    g = Genome.random(np.random.default_rng(0))
    s1, a1 = g.to_models()
    s2, a2 = g.to_models()
    assert isinstance(s1, NeuralSafety)
    assert isinstance(a1, NeuralAllocation)
    # Cached — same object both times.
    assert s1 is s2
    assert a1 is a2


def test_tree_angles_decoding_matches_genome_module_formula():
    """``Genome.tree_angles()`` must use the same scaling as the existing
    ``mechatree.genome._decode_angles``."""
    from mechatree.genome import _decode_angles

    g = Genome.random(np.random.default_rng(0))
    full_row = [0.0] * 9
    full_row[6], full_row[7], full_row[8] = g.angle_genes
    expected = _decode_angles(full_row)
    actual = g.tree_angles()
    assert math.isclose(actual["theta1"], expected["theta1"])
    assert math.isclose(actual["theta2"], expected["theta2"])
    assert math.isclose(actual["gamma1"], expected["gamma1"])
    assert math.isclose(actual["gamma2"], expected["gamma2"])


def test_to_from_dict_round_trip():
    g = Genome.random(np.random.default_rng(0), lineage_id=42)
    d = g.to_dict()
    g2 = Genome.from_dict(d)
    assert g2.safety_weights == g.safety_weights
    assert g2.allocation_weights == g.allocation_weights
    assert g2.angle_genes == g.angle_genes
    assert g2.lineage_id == g.lineage_id


def test_from_champion_constructor():
    g = Genome.from_champion(
        safety_weights=[0.1] * SAFETY_LEN,
        allocation_weights=[0.2] * ALLOCATION_LEN,
        angle_genes=[0.3, 0.4, 0.5],
        lineage_id=99,
    )
    assert g.lineage_id == 99
    assert all(w == 0.1 for w in g.safety_weights)
    assert all(w == 0.2 for w in g.allocation_weights)
    assert g.angle_genes == (0.3, 0.4, 0.5)
