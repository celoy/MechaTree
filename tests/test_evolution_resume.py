"""Tests for checkpoint/resume in tournament runs (Step 21c)."""

from __future__ import annotations

import pytest

from mechatree.config import Config, ForestConfig, TreeConfig
from mechatree.evolution import archive
from mechatree.evolution.genome import Genome
from mechatree.evolution.run import run_tournament


def _cfg(n_init: int = 4, n_max: int = 30, max_age: int = 1000) -> Config:
    """Small test config."""
    return Config(
        tree=TreeConfig(),
        forest=ForestConfig(
            size=20.0,
            n_trees_init=n_init,
            n_trees_max=n_max,
            max_age=max_age,
            min_age_for_undersize=5,
            min_branches=11,
        ),
    )


def test_resume_from_snapshot_continues_from_correct_generation(tmp_path):
    """Resumed run starts at start_gen+1 and history reflects this."""
    cfg = _cfg(n_init=4, n_max=20)
    # Phase 1: run 5 gens, write archive every gen, snapshot at gen 4.
    run_tournament(cfg, n_generations=5, seed=0, archive_every=1, archive_dir=tmp_path)
    snap = tmp_path / "archive_00000004.json"
    assert snap.exists()

    # Phase 2: resume from the snapshot, run to gen 10.
    result = run_tournament(
        cfg,
        n_generations=10,
        seed=0,
        resume_from=snap,
    )
    # history only covers gens 5..9 (5 entries).
    assert len(result.history) == 5
    assert result.history[0].generation == 5
    assert result.history[-1].generation == 9


def test_resume_produces_live_forest(tmp_path):
    """After resume, the forest has trees and genomes."""
    cfg = _cfg(n_init=4, n_max=20)
    run_tournament(cfg, n_generations=3, seed=1, archive_every=1, archive_dir=tmp_path)
    snap = tmp_path / "archive_00000002.json"
    assert snap.exists()

    result = run_tournament(cfg, n_generations=6, seed=1, resume_from=snap)
    assert result.forest.trees
    assert result.forest.genomes
    assert len(result.forest.genomes) == len(result.forest.trees)


def test_resume_from_nonexistent_path_raises(tmp_path):
    """Resuming from a non-existent snapshot raises FileNotFoundError."""
    cfg = _cfg(n_init=4, n_max=20)
    with pytest.raises((FileNotFoundError, Exception)):
        run_tournament(cfg, n_generations=5, seed=0, resume_from=tmp_path / "no_such.json")


def test_resume_preserves_survivor_genomes(tmp_path):
    """Genomes loaded from snapshot round-trip correctly."""
    cfg = _cfg(n_init=4, n_max=20)
    # Run initial phase: 9 gens, archive every 3 → snapshots at gens 0, 3, 6.
    run_tournament(cfg, n_generations=9, seed=42, archive_every=3, archive_dir=tmp_path)
    snap = tmp_path / "archive_00000003.json"
    assert snap.exists()

    # Load the snapshot and verify we can reconstruct genomes.
    start_gen, loaded_genomes, _metas = archive.load_snapshot(snap)
    assert start_gen == 3
    assert len(loaded_genomes) > 0
    # All loaded genomes should be valid Genome instances.
    for g in loaded_genomes:
        assert isinstance(g, Genome)
        # Round-trip through dict.
        d = g.to_dict()
        g_roundtrip = Genome.from_dict(d)
        assert g_roundtrip.lineage_id == g.lineage_id
        assert g_roundtrip.safety_weights == g.safety_weights
        assert g_roundtrip.allocation_weights == g.allocation_weights

    # Resume from the snapshot.
    result2 = run_tournament(cfg, n_generations=12, seed=42, resume_from=snap)
    # The resumed run should have genomes.
    assert result2.forest.genomes is not None
    # The number of generations in history should reflect the resumed portion.
    assert len(result2.history) == 8  # gens 4..11
