"""Curation: cluster surviving lineages into "species", pick a champion
per cluster, write a JSON matching the existing :func:`mechatree.genome.load_champion`
schema so the downstream tooling keeps working unchanged.

The 2-cluster k-means + gap-threshold logic is lifted from
``scripts/strategies_single_tree.py`` (which originally built
``data/S3_champions.json`` from the Fortran ``S3.dat`` archive). The
script now imports from this module so there is exactly one
implementation.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from mechatree.evolution.genome import ALLOCATION_LEN, SAFETY_LEN, Genome

# ----- 2-cluster k-means (lifted verbatim from strategies_single_tree.py) --


def kmeans2(points: np.ndarray, n_iter: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """Tiny 2-cluster k-means with farthest-point init. Returns
    ``(labels, centroids)``."""
    points = np.asarray(points, dtype=np.float64)
    c0 = points[0]
    c1 = points[np.argmax(np.linalg.norm(points - c0, axis=1))]
    centroids = np.stack([c0, c1])
    labels = np.zeros(len(points), dtype=int)
    for _ in range(n_iter):
        d = np.linalg.norm(points[:, None, :] - centroids[None, :, :], axis=2)
        new_labels = d.argmin(axis=1)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels
        for k in (0, 1):
            members = points[labels == k]
            if len(members):
                centroids[k] = members.mean(axis=0)
    return labels, centroids


def detect_species(
    tag_genes: np.ndarray, gap_threshold: float = 0.15
) -> tuple[np.ndarray, np.ndarray]:
    """Cluster into 1 or 2 species. Collapses to a single cluster when the
    centroid gap is below ``gap_threshold``."""
    tag_genes = np.asarray(tag_genes, dtype=np.float64)
    if len(tag_genes) < 2:
        return np.zeros(len(tag_genes), dtype=int), tag_genes.mean(axis=0, keepdims=True)
    labels, centroids = kmeans2(tag_genes)
    if np.linalg.norm(centroids[0] - centroids[1]) < gap_threshold:
        return np.zeros(len(tag_genes), dtype=int), tag_genes.mean(axis=0, keepdims=True)
    return labels, centroids


# ----- tag-gene projection for evolving genomes ----------------------------


def tag_genes_from_genomes(genomes: list[Genome]) -> np.ndarray:
    """Project each Genome to a 2-D tag-gene point.

    For the original Fortran tournament these were two extra genes
    explicitly used as cluster coordinates (``mod_tree.f90`` last two
    indices of the genome). Our 31-gene Python port doesn't carry tag
    genes; we project instead onto the first two angle genes, which
    encode the branching geometry and are the strongest visual
    differentiator between species (matches the
    ``branching angles: theta1=..., theta2=...`` lines printed by
    ``notebooks/03_neural_genome.ipynb``).
    """
    if not genomes:
        return np.zeros((0, 2))
    return np.array(
        [(g.angle_genes[0], g.angle_genes[1]) for g in genomes],
        dtype=np.float64,
    )


# ----- champion picking ----------------------------------------------------


def pick_champions(
    genomes: list[Genome],
    fitness: np.ndarray,
    *,
    gap_threshold: float = 0.15,
) -> list[dict[str, Any]]:
    """Cluster genomes by tag-gene proximity, then per-cluster pick the
    individual with the highest ``fitness`` scalar.

    Returns a list of dicts in the same schema as the entries under
    ``species`` in ``data/S3_champions.json`` — so :func:`mechatree.genome.load_champion`
    consumes the output unchanged.
    """
    if len(genomes) != len(fitness):
        raise ValueError(f"genomes ({len(genomes)}) and fitness ({len(fitness)}) length mismatch")
    if not genomes:
        return []

    tag_genes = tag_genes_from_genomes(genomes)
    labels, centroids = detect_species(tag_genes, gap_threshold=gap_threshold)
    species_ids = sorted(np.unique(labels).tolist())

    champions: list[dict[str, Any]] = []
    fitness = np.asarray(fitness, dtype=np.float64)
    for k in species_ids:
        m = labels == k
        masked = np.where(m, fitness, -np.inf)
        rep = int(np.argmax(masked))
        champ = genomes[rep]
        full_row = _full_row_from_genome(champ, fitness[rep], tag_genes[rep])
        champions.append(
            {
                "species_id": int(k),
                "n_members": int(m.sum()),
                "centroid_tag": centroids[k].tolist() if k < len(centroids) else None,
                "champion_index": rep,
                "champion_tag_genes": tag_genes[rep].tolist(),
                "champion_moment_leaves": 0.0,  # not modelled in the Python port
                "champion_n_seeds": int(fitness[rep]),
                "nn_branch": list(champ.safety_weights),
                "nn_reserve": list(champ.allocation_weights),
                "full_row": full_row,
            }
        )
    return champions


def _full_row_from_genome(genome: Genome, fitness: float, tag_genes: np.ndarray) -> list[float]:
    """Build a 50-element ``full_row`` matching the Fortran S3.dat layout
    so ``_decode_angles`` (in :mod:`mechatree.genome`) keeps working.

    Layout (cols are 0-based):

    * 0..5  — physical prefix (zeroed; we don't track those)
    * 6..8  — 3 angle genes
    * 9..18 — 10 nn_branch weights
    * 19..36 — 18 nn_reserve weights
    * 37..43 — 7 reserved zeros
    * 44..45 — 2 tag genes
    * 46     — moment on leaves (zeroed)
    * 47     — N seeds (the fitness scalar)
    * 48..49 — 2 trailing zeros
    """
    row = [0.0] * 50
    row[6], row[7], row[8] = genome.angle_genes
    row[9 : 9 + SAFETY_LEN] = list(genome.safety_weights)
    row[19 : 19 + ALLOCATION_LEN] = list(genome.allocation_weights)
    row[44], row[45] = float(tag_genes[0]), float(tag_genes[1])
    row[46] = 0.0
    row[47] = float(fitness)
    return row


# ----- whole-forest curation ------------------------------------------------


def from_forest(forest: Any, fitness: np.ndarray) -> dict[str, Any]:
    """Curate a champions JSON payload directly from an evolving Forest.

    ``fitness`` is a per-tree scalar (parallel to ``forest.genomes``);
    typically cumulative seeds dropped per lineage, but the curator is
    indifferent — it just argmax-picks per cluster.
    """
    if forest.genomes is None:
        raise ValueError("from_forest requires an evolving Forest (forest.genomes is None)")
    return _payload(forest.genomes, fitness, source="forest")


def from_archive(path: Path) -> dict[str, Any]:
    """Curate from a snapshot file written by
    :func:`mechatree.evolution.archive.write_snapshot`. Fitness is taken
    from the per-tree ``n_branches`` (a stand-in for "lifetime
    reproductive success" inside one snapshot) — pass a snapshot from
    near the end of a run for the most informative comparison.
    """
    from mechatree.evolution import archive  # local: avoid circular import

    _gen, genomes, metas = archive.load_snapshot(Path(path))
    fitness = np.array([m["n_branches"] for m in metas], dtype=np.float64)
    return _payload(genomes, fitness, source=str(path))


def write(payload: dict[str, Any], out_path: Path) -> Path:
    """Write a curated payload to ``out_path`` as JSON. Creates parent
    dirs as needed."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2))
    return out_path


def _payload(genomes: list[Genome], fitness: np.ndarray, *, source: str) -> dict[str, Any]:
    tag_genes = tag_genes_from_genomes(genomes)
    labels, _centroids = detect_species(tag_genes)
    species = pick_champions(genomes, fitness)
    return {
        "dataset": "mechatree.evolution",
        "source_path": source,
        "n_individuals": len(genomes),
        "schema": {
            "col_nnbranch": [9, 9 + SAFETY_LEN],
            "col_nnreserve": [9 + SAFETY_LEN, 9 + SAFETY_LEN + ALLOCATION_LEN],
            "col_tag": [44, 46],
            "col_moment_leaves": 46,
            "col_n_seeds": 47,
            "note": "0-based python col indices; layout mirrors S3.dat so "
            "mechatree.genome.load_champion (which reads full_row[6:9] for "
            "angles and indexes nn_branch / nn_reserve) keeps working "
            "against this file unchanged.",
        },
        "tag_genes": tag_genes.tolist(),
        "species_labels": labels.tolist(),
        "species": species,
    }


__all__ = [
    "detect_species",
    "from_archive",
    "from_forest",
    "kmeans2",
    "pick_champions",
    "tag_genes_from_genomes",
    "write",
]
