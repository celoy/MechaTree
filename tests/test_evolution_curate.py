"""Tests for the champion-curation port (Step 21)."""

from __future__ import annotations

import json

import numpy as np
import pytest

from mechatree.evolution import archive, curate
from mechatree.evolution.genome import ALLOCATION_LEN, SAFETY_LEN, Genome
from mechatree.evolution.run import run_tournament
from mechatree.genome import load_champion
from tests.test_evolution_forest import _cfg


def _genome(angle_x: float, angle_y: float, lineage_id: int = 0) -> Genome:
    """A genome with deterministic weights + the first two angle genes set
    to ``(angle_x, angle_y)`` so the cluster lives where we want."""
    return Genome(
        safety_weights=tuple([0.5] * SAFETY_LEN),
        allocation_weights=tuple([0.5] * ALLOCATION_LEN),
        angle_genes=(angle_x, angle_y, 0.5),
        lineage_id=lineage_id,
    )


def test_kmeans2_separates_well_separated_clusters():
    points = np.array([[0.0, 0.0], [0.05, 0.05], [1.0, 1.0], [0.95, 0.95]])
    labels, _centroids = curate.kmeans2(points)
    assert labels[0] == labels[1]
    assert labels[2] == labels[3]
    assert labels[0] != labels[2]


def test_detect_species_collapses_near_clusters():
    """Points all clumped under the gap threshold ⇒ one species."""
    rng = np.random.default_rng(0)
    points = rng.normal(0.5, 0.02, size=(20, 2))
    labels, centroids = curate.detect_species(points, gap_threshold=0.15)
    assert set(labels.tolist()) == {0}
    assert centroids.shape == (1, 2)


def test_pick_champions_picks_best_per_cluster():
    """Two well-separated clusters of 3 genomes each. In each cluster the
    champion should be the genome with the highest fitness."""
    genomes = [
        _genome(0.1, 0.1, lineage_id=0),  # cluster A
        _genome(0.15, 0.15, lineage_id=1),
        _genome(0.05, 0.05, lineage_id=2),
        _genome(0.9, 0.9, lineage_id=3),  # cluster B
        _genome(0.85, 0.85, lineage_id=4),
        _genome(0.95, 0.95, lineage_id=5),
    ]
    fitness = np.array([1.0, 5.0, 2.0, 3.0, 7.0, 4.0])
    champions = curate.pick_champions(genomes, fitness)
    assert len(champions) == 2
    by_species = {c["species_id"]: c for c in champions}
    # The cluster-A champion is index 1 (fitness 5); cluster-B is index 4 (7).
    assert by_species[0]["champion_index"] in (0, 1, 2)
    assert by_species[1]["champion_index"] in (3, 4, 5)
    # The chosen index per cluster must be the argmax-fitness one.
    cluster_a_indices = [i for i in (0, 1, 2)]
    cluster_b_indices = [i for i in (3, 4, 5)]
    expected_a = max(cluster_a_indices, key=lambda i: fitness[i])
    expected_b = max(cluster_b_indices, key=lambda i: fitness[i])
    assert by_species[0]["champion_index"] == expected_a
    assert by_species[1]["champion_index"] == expected_b


def test_curated_payload_round_trips_through_load_champion(tmp_path):
    cfg = _cfg(n_init=4, n_max=20)
    result = run_tournament(
        cfg,
        n_generations=12,
        seed=42,
        champions_path=tmp_path / "champions.json",
    )
    assert result.champions_path is not None and result.champions_path.exists()
    payload = json.loads(result.champions_path.read_text())
    assert "species" in payload
    assert payload["species"], "no species in curated output"
    # `mt.load_champion` must consume the output as-is.
    for sp in payload["species"]:
        safety, alloc, angles, non_coding = load_champion(
            result.champions_path, species_id=sp["species_id"]
        )
        assert safety is not None
        assert alloc is not None
        assert {"theta1", "theta2", "gamma1", "gamma2"} <= set(angles)
        assert non_coding["species_id"] == sp["species_id"]


def test_archive_snapshot_round_trip(tmp_path):
    cfg = _cfg(n_init=3, n_max=10)
    result = run_tournament(
        cfg,
        n_generations=5,
        seed=0,
        archive_every=1,
        archive_dir=tmp_path,
    )
    snapshots = sorted(tmp_path.glob("archive_*.json"))
    assert snapshots
    gen, genomes, metas = archive.load_snapshot(snapshots[-1])
    assert gen == len(snapshots) - 1
    assert len(genomes) == len(metas) == len(result.forest.trees)
    # The genomes round-trip cleanly.
    for g in genomes:
        assert len(g.safety_weights) == SAFETY_LEN
        assert len(g.allocation_weights) == ALLOCATION_LEN


def test_from_archive_curates(tmp_path):
    cfg = _cfg(n_init=4, n_max=20)
    run_tournament(cfg, n_generations=8, seed=1, archive_every=2, archive_dir=tmp_path)
    snapshots = sorted(tmp_path.glob("archive_*.json"))
    payload = curate.from_archive(snapshots[-1])
    assert payload["species"]
    out_path = tmp_path / "curated.json"
    curate.write(payload, out_path)
    # And load_champion can read it.
    safety, _, _, _ = load_champion(out_path, species_id=payload["species"][0]["species_id"])
    assert safety is not None


def test_pick_champions_rejects_mismatched_fitness():
    genomes = [_genome(0.1, 0.1)]
    with pytest.raises(ValueError):
        curate.pick_champions(genomes, np.array([1.0, 2.0]))
