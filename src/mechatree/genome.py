"""Genome models — the per-branch decision functions consulted by growth.

Two flavours ship today:

* ``ConstantSafety`` / ``ConstantAllocation`` — return fixed scalars regardless
  of input. Useful for testing and minimum-viable simulations.
* ``NeuralSafety`` / ``NeuralAllocation`` — the 3-layer tanh networks evolved
  in Eloy et al. (Nat Commun 2017), driven by 10- and 18-element weight
  vectors. Load champion weights from JSON via :func:`load_champion` or
  :func:`load_all_champions` (e.g. ``data/S3_champions.json``).

The C++ side defines an abstract base class with a virtual ``compute`` method;
both flavours subclass it so the growth code that calls them never changes.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from mechatree._core._core import (
    PyAllocationModel,
    PyConstantAllocation,
    PyConstantSafety,
    PyNeuralAllocation,
    PyNeuralSafety,
    PySafetyModel,
)

# Public aliases — the wrapper-prefix is an implementation detail.
SafetyModel = PySafetyModel
AllocationModel = PyAllocationModel
ConstantSafety = PyConstantSafety
ConstantAllocation = PyConstantAllocation
NeuralSafety = PyNeuralSafety
NeuralAllocation = PyNeuralAllocation


def _champion(payload: dict[str, Any], species_id: int) -> dict[str, Any]:
    for s in payload.get("species", []):
        if int(s["species_id"]) == species_id:
            return s
    available = [int(s["species_id"]) for s in payload.get("species", [])]
    raise ValueError(f"species_id={species_id} not found in champion file; available: {available}")


def load_champion(
    path: str | Path, species_id: int = 0
) -> tuple[NeuralSafety, NeuralAllocation, dict[str, Any]]:
    """Build ``(NeuralSafety, NeuralAllocation, metadata)`` from a champion JSON.

    ``path`` points at a file written by ``scripts/strategies_single_tree.py``
    (or anything with the same schema). ``species_id`` selects which champion
    from the ``species`` list. Metadata holds the bookkeeping fields
    (``champion_index``, ``champion_moment_leaves``, ``champion_n_seeds``,
    ``centroid_tag``, …) so callers can log which genome they loaded.
    """
    payload = json.loads(Path(path).read_text())
    species = _champion(payload, species_id)
    meta = {k: v for k, v in species.items() if k not in ("nn_branch", "nn_reserve", "full_row")}
    meta["dataset"] = payload.get("dataset")
    meta["source_path"] = str(path)
    return NeuralSafety(species["nn_branch"]), NeuralAllocation(species["nn_reserve"]), meta


def load_all_champions(
    path: str | Path,
) -> list[tuple[NeuralSafety, NeuralAllocation, dict[str, Any]]]:
    """Same as :func:`load_champion` but returns every species in the file."""
    payload = json.loads(Path(path).read_text())
    out = []
    for s in payload.get("species", []):
        meta = {k: v for k, v in s.items() if k not in ("nn_branch", "nn_reserve", "full_row")}
        meta["dataset"] = payload.get("dataset")
        meta["source_path"] = str(path)
        out.append((NeuralSafety(s["nn_branch"]), NeuralAllocation(s["nn_reserve"]), meta))
    return out


def models_from_config(gc, base_dir: Path | None = None) -> tuple[SafetyModel, AllocationModel]:
    """Build ``(safety, allocation)`` from a :class:`GenomeConfig`.

    If ``gc.neural_from`` is set, load ``NeuralSafety`` + ``NeuralAllocation``
    from the referenced JSON champion. ``neural_from["path"]`` may be absolute
    or relative; if relative and ``base_dir`` is given, it is resolved against
    ``base_dir`` (typically the YAML file's parent directory), otherwise CWD.

    Otherwise build ``ConstantSafety`` / ``ConstantAllocation`` from the scalar
    fields.
    """
    if gc is None:
        return (
            ConstantSafety(3.0),
            ConstantAllocation(p_seeds=0.1, p_leaves=0.5, phototropism=0.5),
        )
    if gc.neural_from is not None:
        spec = gc.neural_from
        raw_path = Path(spec["path"])
        if not raw_path.is_absolute() and base_dir is not None:
            raw_path = (Path(base_dir) / raw_path).resolve()
        species_id = int(spec.get("species_id", 0))
        safety, allocation, _ = load_champion(raw_path, species_id)
        return safety, allocation
    return (
        ConstantSafety(gc.safety),
        ConstantAllocation(p_seeds=gc.p_seeds, p_leaves=gc.p_leaves, phototropism=gc.phototropism),
    )


__all__ = [
    "AllocationModel",
    "ConstantAllocation",
    "ConstantSafety",
    "NeuralAllocation",
    "NeuralSafety",
    "SafetyModel",
    "load_all_champions",
    "load_champion",
    "models_from_config",
]
