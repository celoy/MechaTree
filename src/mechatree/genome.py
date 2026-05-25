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
import math
from pathlib import Path
from typing import Any

from mechatree._core._core import (
    PyAllocationModel,
    PyCallbackAllocation,
    PyCallbackSafety,
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
CallbackSafety = PyCallbackSafety
CallbackAllocation = PyCallbackAllocation


def _champion(payload: dict[str, Any], species_id: int) -> dict[str, Any]:
    for s in payload.get("species", []):
        if int(s["species_id"]) == species_id:
            return s
    available = [int(s["species_id"]) for s in payload.get("species", [])]
    raise ValueError(f"species_id={species_id} not found in champion file; available: {available}")


def _decode_angles(full_row: list[float] | None) -> dict[str, float] | None:
    """Decode the four branching-angle radians from a full_row.

    Returns ``None`` if ``full_row`` is missing or too short. See
    ``mod_tree.f90:108-111`` for the Fortran formula:
    genome[0..2] (0-based; full_row[6..8] after the dump's 6-col position
    prefix) → (theta1, -theta2, gamma, gamma) scaled by pi/2 and 2*pi.
    """
    if full_row is None or len(full_row) < 9:
        return None
    g1, g2, g3 = float(full_row[6]), float(full_row[7]), float(full_row[8])
    return {
        "theta1": g1 * math.pi / 2.0,
        "theta2": -g2 * math.pi / 2.0,
        "gamma1": g3 * 2.0 * math.pi,
        "gamma2": g3 * 2.0 * math.pi,
    }


def load_champion(
    path: str | Path, species_id: int = 0
) -> tuple[NeuralSafety, NeuralAllocation, dict[str, float], dict[str, Any]]:
    """Build ``(NeuralSafety, NeuralAllocation, angles, non_coding)`` from a champion JSON.

    Returns a 4-tuple so a caller has everything the Fortran genome encodes,
    split into **coding** (drives tree shape) and **non-coding** (bookkeeping
    that identifies which evolved individual you loaded):

    1. ``NeuralSafety`` — coding. 10-weight branch network for safety decisions.
    2. ``NeuralAllocation`` — coding. 18-weight reserve network for
       (p_seeds, p_leaves, phototropism).
    3. ``angles`` — coding. Dict ``{theta1, theta2, gamma1, gamma2}`` (radians)
       decoded from the first three slots of the Fortran genome
       (``mod_tree.f90:108-111``). Drop-in for ``replace(cfg.tree, **angles)``.
    4. ``non_coding`` — bookkeeping fields (``champion_index``,
       ``champion_moment_leaves``, ``champion_n_seeds``, ``centroid_tag``,
       ``species_id``, ``dataset``, ``source_path``, …) so callers can log
       which genome they loaded.

    ``path`` points at a file written by ``scripts/strategies_single_tree.py``
    (or anything with the same schema). ``species_id`` selects which champion
    from the ``species`` list.
    """
    payload = json.loads(Path(path).read_text())
    species = _champion(payload, species_id)
    non_coding = {
        k: v for k, v in species.items() if k not in ("nn_branch", "nn_reserve", "full_row")
    }
    non_coding["dataset"] = payload.get("dataset")
    non_coding["source_path"] = str(path)
    angles = _decode_angles(species.get("full_row"))
    if angles is None:
        raise ValueError(
            f"champion file at {path} has no full_row (or fewer than 9 cols); "
            "cannot decode branching angles."
        )
    return (
        NeuralSafety(species["nn_branch"]),
        NeuralAllocation(species["nn_reserve"]),
        angles,
        non_coding,
    )


def champion_angles(path: str | Path, species_id: int = 0) -> dict[str, float]:
    """Decode the four branching angles for ``species_id`` from a champion JSON.

    Equivalent to ``load_champion(path, species_id)[2]`` but cheaper when the
    NN weights aren't needed (e.g. quick angle-table dumps).
    """
    payload = json.loads(Path(path).read_text())
    species = _champion(payload, species_id)
    angles = _decode_angles(species.get("full_row"))
    if angles is None:
        raise ValueError(
            f"champion file at {path} has no full_row (or fewer than 9 cols); "
            "cannot decode branching angles."
        )
    return angles


def load_all_champions(
    path: str | Path,
) -> list[tuple[NeuralSafety, NeuralAllocation, dict[str, float], dict[str, Any]]]:
    """Same as :func:`load_champion` but returns every species in the file."""
    payload = json.loads(Path(path).read_text())
    out = []
    for s in payload.get("species", []):
        non_coding = {
            k: v for k, v in s.items() if k not in ("nn_branch", "nn_reserve", "full_row")
        }
        non_coding["dataset"] = payload.get("dataset")
        non_coding["source_path"] = str(path)
        angles = _decode_angles(s.get("full_row"))
        if angles is None:
            raise ValueError(
                f"champion file at {path} has no full_row for species "
                f"{s.get('species_id')}; cannot decode branching angles."
            )
        out.append(
            (NeuralSafety(s["nn_branch"]), NeuralAllocation(s["nn_reserve"]), angles, non_coding)
        )
    return out


def models_from_config(
    gc, base_dir: Path | None = None
) -> tuple[SafetyModel, AllocationModel, dict[str, float] | None]:
    """Build ``(safety, allocation, angles)`` from a :class:`GenomeConfig`.

    ``angles`` is ``None`` unless ``gc.neural_from`` is set; in that case
    it carries the champion's four branching angles
    (``{theta1, theta2, gamma1, gamma2}``, radians) so the simulator can
    override ``cfg.tree.theta*`` / ``cfg.tree.gamma*`` automatically. Without
    this, the YAML ``tree:`` block's defaults would silently shadow the
    champion's geometry and the tree wouldn't reproduce the paper's shapes.

    Resolution order:

    1. ``gc.neural_from`` set → ``load_champion(...)`` returns
       ``NeuralSafety`` + ``NeuralAllocation`` + the decoded angles.
       ``neural_from["path"]`` may be absolute or relative; if relative
       and ``base_dir`` is given, it is resolved against ``base_dir``
       (typically the YAML file's parent directory), otherwise CWD.
    2. Any scalar field is a string expression → compile all four fields
       through :mod:`mechatree.sympy_genome` into ``CallbackSafety`` /
       ``CallbackAllocation``. Numeric fields are passed through as
       constants. Requires the ``sympy`` optional extra
       (``pip install 'mechatree[sympy]'``). ``angles`` returned as ``None``.
    3. Otherwise (all scalars) → ``ConstantSafety`` / ``ConstantAllocation``.
       ``angles`` returned as ``None``.
    """
    if gc is None:
        return (
            ConstantSafety(3.0),
            ConstantAllocation(p_seeds=0.1, p_leaves=0.5, phototropism=0.5),
            None,
        )
    if gc.neural_from is not None:
        spec = gc.neural_from
        raw_path = Path(spec["path"])
        if not raw_path.is_absolute() and base_dir is not None:
            raw_path = (Path(base_dir) / raw_path).resolve()
        species_id = int(spec.get("species_id", 0))
        safety, allocation, angles, _ = load_champion(raw_path, species_id)
        return safety, allocation, angles

    has_expr = any(
        isinstance(v, str) for v in (gc.safety, gc.p_seeds, gc.p_leaves, gc.phototropism)
    )
    if has_expr:
        try:
            from mechatree.sympy_genome import sympy_allocation, sympy_safety
        except ImportError as exc:
            raise ImportError(
                "GenomeConfig contains a string expression but SymPy isn't installed. "
                "Install with: pip install 'mechatree[sympy]'"
            ) from exc

        safety_value: SafetyModel
        if isinstance(gc.safety, str):
            safety_value = sympy_safety(gc.safety)
        else:
            safety_value = ConstantSafety(gc.safety)
        allocation_value = sympy_allocation(
            p_seeds=gc.p_seeds,
            p_leaves=gc.p_leaves,
            phototropism=gc.phototropism,
        )
        return safety_value, allocation_value, None

    return (
        ConstantSafety(gc.safety),
        ConstantAllocation(p_seeds=gc.p_seeds, p_leaves=gc.p_leaves, phototropism=gc.phototropism),
        None,
    )


__all__ = [
    "AllocationModel",
    "CallbackAllocation",
    "CallbackSafety",
    "ConstantAllocation",
    "ConstantSafety",
    "NeuralAllocation",
    "NeuralSafety",
    "SafetyModel",
    "load_all_champions",
    "load_champion",
    "models_from_config",
]
