"""Archive snapshots for an evolving Forest.

Two writers:

* :func:`write_snapshot` — dumps every ``archive_every`` generations to
  ``out_dir/archive_{gen:08d}.json``. Captures the live population:
  per-tree ``(lineage_id, age, n_branches, biomass, location, genome)``.
* :func:`write_champions` — runs :func:`mechatree.evolution.curate.curate`
  against the final population and writes a ``champions.json`` matching
  the existing :func:`mechatree.genome.load_champion` schema.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from mechatree.evolution.genome import Genome

# Forward reference only — avoid a hard circular import.
# (mechatree.forest -> mechatree.evolution.archive would be a cycle.)


def snapshot_payload(forest: Any, generation: int) -> dict[str, Any]:
    """Build the JSON-ready dict for one generation's archive snapshot.

    Useful when you'd rather hold the snapshot in memory than write it.
    """
    if forest.genomes is None:
        raise ValueError(
            "snapshot_payload requires an evolving Forest "
            "(forest.genomes is None — not in evolution mode)"
        )
    trees_payload = []
    for tree, age, genome in zip(forest.trees, forest.ages, forest.genomes, strict=True):
        x, y, _ = tree.get_location(0)
        n_b = tree.get_number_of_branches()
        biomass = 0.0
        for i in range(n_b):
            d = tree.get_diameter(i)
            biomass += tree.get_length(i) * 0.25 * math.pi * d * d
        trees_payload.append(
            {
                "lineage_id": int(genome.lineage_id),
                "age": int(age),
                "n_branches": int(n_b),
                "n_leaves": int(tree.get_total_leaves()),
                "biomass": float(biomass),
                "location": [float(x), float(y)],
                "genome": genome.to_dict(),
            }
        )
    return {
        "generation": int(generation),
        "n_trees": len(forest.trees),
        "n_lineages_alive": len({g.lineage_id for g in forest.genomes}),
        "trees": trees_payload,
    }


def write_snapshot(forest: Any, generation: int, out_dir: Path) -> Path:
    """Write ``archive_{generation:08d}.json`` under ``out_dir``.

    Returns the written path. Creates ``out_dir`` if missing.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"archive_{generation:08d}.json"
    payload = snapshot_payload(forest, generation)
    path.write_text(json.dumps(payload))
    return path


def load_snapshot(path: Path) -> tuple[int, list[Genome], list[dict[str, Any]]]:
    """Inverse of :func:`write_snapshot`. Returns
    ``(generation, [Genome, ...], [tree_meta, ...])`` where ``tree_meta``
    drops the genome and keeps the spatial/age bookkeeping.
    """
    payload = json.loads(Path(path).read_text())
    genomes = [Genome.from_dict(t["genome"]) for t in payload["trees"]]
    metas = [{k: v for k, v in t.items() if k != "genome"} for t in payload["trees"]]
    return int(payload["generation"]), genomes, metas


__all__ = [
    "load_snapshot",
    "snapshot_payload",
    "write_snapshot",
]
