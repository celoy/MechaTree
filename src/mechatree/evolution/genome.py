"""The heritable Genome carried by every tree in an evolving Forest.

A Genome is 31 floats: 3 raw "angle genes" (mapped to ``theta1``,
``theta2``, ``gamma1``, ``gamma2`` via the existing
:func:`mechatree.genome._decode_angles` formula at use time), 10 weights
for the per-branch safety NN, and 18 weights for the reserve allocation
NN. Layout matches the coding part of the Fortran genome in
``legacy/fortran/mod_tree.f90:47-52``.

Mutation operates per-locus with Gaussian noise (Fortran defaults
``sigma=0.005``, ``p_locus=0.05`` from ``legacy/fortran/Evolution.ini``).
Weights are clipped to ``[0, 1]``; angle genes are left unclamped — the
existing decoder scales by π / 2π so values outside ``[0, 1]`` just give
slightly larger angles, no harm done.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from mechatree.genome import NeuralAllocation, NeuralSafety

SAFETY_LEN = 10
ALLOCATION_LEN = 18


@dataclass(frozen=True)
class Genome:
    """One tree's heritable parameters.

    ``lineage_id`` is the founder index; mutation preserves it so a
    sub-population can be tracked back to the initial sapling it descended
    from. ``ForestStats.n_lineages_alive`` is the unique-``lineage_id``
    count over surviving trees.
    """

    safety_weights: tuple[float, ...]  # length SAFETY_LEN
    allocation_weights: tuple[float, ...]  # length ALLOCATION_LEN
    angle_genes: tuple[float, float, float]  # raw [0, 1]-ish g1, g2, g3
    lineage_id: int = 0
    # Cache slot for the materialised C++ models. Populated lazily by
    # ``to_models()``; not a constructor argument.
    _models_cache: tuple[NeuralSafety, NeuralAllocation] | None = field(
        default=None, init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if len(self.safety_weights) != SAFETY_LEN:
            raise ValueError(
                f"safety_weights must have {SAFETY_LEN} elements, got {len(self.safety_weights)}"
            )
        if len(self.allocation_weights) != ALLOCATION_LEN:
            raise ValueError(
                f"allocation_weights must have {ALLOCATION_LEN} elements, "
                f"got {len(self.allocation_weights)}"
            )
        if len(self.angle_genes) != 3:
            raise ValueError(f"angle_genes must have 3 elements, got {len(self.angle_genes)}")

    # ----- construction -------------------------------------------------------

    @classmethod
    def random(cls, rng: np.random.Generator, *, lineage_id: int = 0) -> Genome:
        """Uniform ``[0, 1]`` per-locus — same shape as the Fortran initial
        population (``mod_evolu.f90`` ``init_random`` call)."""
        return cls(
            safety_weights=tuple(rng.random(SAFETY_LEN).tolist()),
            allocation_weights=tuple(rng.random(ALLOCATION_LEN).tolist()),
            angle_genes=tuple(rng.random(3).tolist()),  # type: ignore[arg-type]
            lineage_id=lineage_id,
        )

    @classmethod
    def from_champion(
        cls,
        safety_weights: tuple[float, ...] | list[float],
        allocation_weights: tuple[float, ...] | list[float],
        angle_genes: tuple[float, float, float] | list[float],
        *,
        lineage_id: int = 0,
    ) -> Genome:
        """Build from a champion record (the loader path in
        :mod:`mechatree.genome`)."""
        ag = tuple(float(x) for x in angle_genes)
        if len(ag) != 3:
            raise ValueError(f"angle_genes must have 3 elements, got {len(ag)}")
        return cls(
            safety_weights=tuple(float(x) for x in safety_weights),
            allocation_weights=tuple(float(x) for x in allocation_weights),
            angle_genes=ag,  # type: ignore[arg-type]
            lineage_id=lineage_id,
        )

    # ----- variation ----------------------------------------------------------

    def mutate(
        self,
        rng: np.random.Generator,
        *,
        sigma: float = 0.005,
        p_locus: float = 0.05,
    ) -> Genome:
        """Per-locus Gaussian mutation. Weights clipped to ``[0, 1]``;
        ``angle_genes`` left unclamped. Returns a new Genome with the
        same ``lineage_id``."""
        sw = _mutate_clipped(np.asarray(self.safety_weights), rng, sigma, p_locus)
        aw = _mutate_clipped(np.asarray(self.allocation_weights), rng, sigma, p_locus)
        ag = _mutate_unclamped(np.asarray(self.angle_genes), rng, sigma, p_locus)
        return Genome(
            safety_weights=tuple(sw.tolist()),
            allocation_weights=tuple(aw.tolist()),
            angle_genes=(float(ag[0]), float(ag[1]), float(ag[2])),
            lineage_id=self.lineage_id,
        )

    # ----- accessors ----------------------------------------------------------

    def to_models(self) -> tuple[NeuralSafety, NeuralAllocation]:
        """Materialise the two C++ NN models, cached per Genome instance."""
        if self._models_cache is None:
            models = (
                NeuralSafety(list(self.safety_weights)),
                NeuralAllocation(list(self.allocation_weights)),
            )
            # Cache through frozen-dataclass bypass.
            object.__setattr__(self, "_models_cache", models)
        return self._models_cache  # type: ignore[return-value]

    def tree_angles(self) -> dict[str, float]:
        """Decode the 3 angle genes into the 4 keyword arguments expected
        by :class:`mechatree.config.TreeConfig`. Formula matches
        :func:`mechatree.genome._decode_angles`."""
        g1, g2, g3 = self.angle_genes
        return {
            "theta1": g1 * math.pi / 2.0,
            "theta2": -g2 * math.pi / 2.0,
            "gamma1": g3 * 2.0 * math.pi,
            "gamma2": g3 * 2.0 * math.pi,
        }

    # ----- I/O ----------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "safety_weights": list(self.safety_weights),
            "allocation_weights": list(self.allocation_weights),
            "angle_genes": list(self.angle_genes),
            "lineage_id": self.lineage_id,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Genome:
        return cls(
            safety_weights=tuple(float(x) for x in d["safety_weights"]),
            allocation_weights=tuple(float(x) for x in d["allocation_weights"]),
            angle_genes=tuple(float(x) for x in d["angle_genes"]),  # type: ignore[arg-type]
            lineage_id=int(d.get("lineage_id", 0)),
        )


# ----- internal helpers ------------------------------------------------------


def _mutate_clipped(
    locus: np.ndarray,
    rng: np.random.Generator,
    sigma: float,
    p_locus: float,
) -> np.ndarray:
    mask = rng.random(locus.shape) < p_locus
    if not mask.any():
        # Fortran's "mandatory mutation" guard: if no locus was selected,
        # force exactly one at random so the child differs from the parent.
        idx = int(rng.integers(0, locus.size))
        mask = np.zeros_like(mask)
        mask[idx] = True
    noise = rng.normal(0.0, sigma, size=locus.shape)
    out = np.where(mask, locus + noise, locus)
    return np.clip(out, 0.0, 1.0)


def _mutate_unclamped(
    locus: np.ndarray,
    rng: np.random.Generator,
    sigma: float,
    p_locus: float,
) -> np.ndarray:
    mask = rng.random(locus.shape) < p_locus
    if not mask.any():
        idx = int(rng.integers(0, locus.size))
        mask = np.zeros_like(mask)
        mask[idx] = True
    noise = rng.normal(0.0, sigma, size=locus.shape)
    return np.where(mask, locus + noise, locus)


__all__ = ["ALLOCATION_LEN", "SAFETY_LEN", "Genome"]
