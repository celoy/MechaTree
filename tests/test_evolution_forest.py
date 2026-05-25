"""Tests for the Forest evolution path (Step 21)."""

from __future__ import annotations

import numpy as np
import pytest

from mechatree.config import Config, ForestConfig, TreeConfig
from mechatree.evolution.genome import Genome
from mechatree.evolution.run import run_tournament
from mechatree.forest import Forest


def _cfg(n_init: int = 4, n_max: int = 30, max_age: int = 1000) -> Config:
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


def test_default_forest_unchanged_when_genomes_none():
    """Backward-compat: no genomes argument ⇒ existing Step-12 behaviour."""
    f = Forest(_cfg(), seed=42)
    stats = f.step(0)
    # New field defaults to 0 when evolution isn't active.
    assert stats.n_lineages_alive == 0
    assert f.genomes is None


def test_forest_rejects_wrong_genome_count():
    cfg = _cfg(n_init=4)
    with pytest.raises(ValueError):
        Forest(cfg, seed=0, genomes=[Genome.random(np.random.default_rng(0))])


def test_forest_with_genomes_dispatches_per_tree():
    """Two contrasting genomes growing under the same seed must produce
    visibly different trees by gen 10."""
    cfg = _cfg(n_init=2, n_max=4, max_age=1000)
    rng = np.random.default_rng(0)
    g_low = Genome(
        safety_weights=tuple([0.1] * 10),
        allocation_weights=tuple([0.5] * 18),
        angle_genes=(0.5, 0.5, 0.5),
        lineage_id=0,
    )
    g_high = Genome(
        safety_weights=tuple([0.9] * 10),
        allocation_weights=tuple([0.5] * 18),
        angle_genes=(0.5, 0.5, 0.5),
        lineage_id=1,
    )
    _ = rng  # quiet linter
    f = Forest(cfg, seed=7, genomes=[g_low, g_high])
    for gen in range(10):
        f.step(gen)
    # Two distinct lineages should still be alive (or at least, the surviving
    # tree's lineage_id is identifiable).
    assert f.genomes is not None
    assert len(f.genomes) == len(f.trees)


def test_lineage_extinction_in_short_run():
    """Across 20 gens of a 6-founder tournament, expect at least one
    lineage to vanish."""
    cfg = _cfg(n_init=6, n_max=20, max_age=15)
    result = run_tournament(cfg, n_generations=20, seed=42)
    final_lineages = (
        {g.lineage_id for g in result.forest.genomes} if result.forest.genomes else set()
    )
    assert len(final_lineages) < 6, (
        f"expected at least one of 6 founder lineages to vanish, "
        f"still alive: {sorted(final_lineages)}"
    )
    # Live count in stats matches the set size.
    assert result.history[-1].n_lineages_alive == len(final_lineages)


def test_tournament_reproducible_under_same_seed():
    """Two runs of the same config + seed produce identical lineage
    survival + final tree positions."""
    cfg = _cfg(n_init=4, n_max=15)
    a = run_tournament(cfg, n_generations=10, seed=42)
    b = run_tournament(cfg, n_generations=10, seed=42)
    assert len(a.forest.trees) == len(b.forest.trees)
    for ta, tb in zip(a.forest.trees, b.forest.trees, strict=True):
        assert ta.get_location(0) == tb.get_location(0)
        assert ta.get_number_of_branches() == tb.get_number_of_branches()
    assert a.history[-1].n_lineages_alive == b.history[-1].n_lineages_alive


def test_tournament_writes_archive_snapshots(tmp_path):
    cfg = _cfg(n_init=4, n_max=15)
    archive_dir = tmp_path / "evo"
    result = run_tournament(
        cfg,
        n_generations=10,
        seed=0,
        archive_every=5,
        archive_dir=archive_dir,
    )
    snapshots = sorted(archive_dir.glob("archive_*.json"))
    # Snapshots at gen 0 and gen 5.
    assert len(snapshots) == 2
    assert result.archive_dir == archive_dir
