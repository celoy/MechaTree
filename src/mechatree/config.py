"""YAML configuration for the single-tree simulator (Step 11).

Direct translation of the parameters in ``legacy_fortran/Forest.ini``, regrouped
semantically and renamed to ``snake_case``. The Fortran names are noted in
docstrings so cross-referencing the paper or the legacy code stays mechanical.

Only the sections actually consumed by Step 11 (single tree) are validated
strictly; ``forest``, ``evolution``, ``init``, ``io`` are tolerated for
forward-compat with later steps but not parsed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class TreeConfig:
    """Scalar parameters for a single tree's mechanics + growth.

    Mirrors the ``tree:`` section of the YAML and the branching-angle scalars
    that the Fortran ``new_tree`` derives from the genome. Defaults match
    ``legacy_fortran/Forest.ini``.
    """

    # Direct YAML inputs --------------------------------------------------
    twig_length: float = 1.0  # was TwigLength
    twig_diameter: float = 0.1  # was TwigDiameter
    leaf_surface: float = 0.25  # was LeafSurface (S0)
    cauchy: float = 2.0e-5  # was Cauchy (Cy) — runtime value from Forest.ini
    volume_ratio_leaf: float = 8.0  # was VolumeRatioLeaf — runtime Forest.ini value
    maintenance_h: float = 0.02  # was MaintenanceH
    max_branches: int = 10000  # was Nmax (advisory only in Step 11)

    # Branching angles. The Fortran derives these from genes 1..3; we expose
    # them directly. Defaults give a plausible binary-tree shape.
    theta1: float = math.pi / 4  # angle of left daughter (radians)
    theta2: float = -math.pi / 4  # angle of right daughter
    gamma1: float = 0.0  # twist of left daughter
    gamma2: float = math.pi  # twist of right daughter

    def __post_init__(self) -> None:
        for name in (
            "twig_length",
            "twig_diameter",
            "leaf_surface",
            "cauchy",
            "volume_ratio_leaf",
        ):
            v = getattr(self, name)
            if v <= 0.0:
                raise ValueError(f"TreeConfig.{name} must be positive, got {v}")
        if self.maintenance_h < 0.0:
            raise ValueError("TreeConfig.maintenance_h must be non-negative")
        if self.max_branches <= 0:
            raise ValueError("TreeConfig.max_branches must be positive")

    @property
    def volume_twig(self) -> float:
        """``0.25 * pi * twig_length * twig_diameter**2`` — base volume unit."""
        return 0.25 * math.pi * self.twig_length * self.twig_diameter**2

    @property
    def volume_per_leaf(self) -> float:
        """Photosynthate produced per leaf at light=1 — drives secondary_growth."""
        return self.volume_ratio_leaf * self.volume_twig


@dataclass(frozen=True)
class LightConfig:
    """Parameters for the hemispherical light integration."""

    size_leaf: float = 1.0  # was SizeLeaf — 2D shadow-grid cell size
    n_elevations: int = 4  # was Nelev
    n_azimuths: int = 8  # was Nazim

    def __post_init__(self) -> None:
        if self.n_elevations <= 0:
            raise ValueError("LightConfig.n_elevations must be positive")
        if self.n_azimuths <= 0:
            raise ValueError("LightConfig.n_azimuths must be positive")
        if self.size_leaf <= 0.0:
            raise ValueError("LightConfig.size_leaf must be positive")


@dataclass(frozen=True)
class ForestConfig:
    """Parameters for a forest of trees (Step 12).

    ``n_trees_init`` trees are seeded uniformly across a disk of radius
    ``size``. The death rule mirrors ``legacy_fortran/Forest.f90:283`` and
    is configurable so studies of self-thinning can tweak it without
    forking the orchestrator.
    """

    size: float = 100.0  # was SizeForest — disk radius
    n_trees_init: int = 100  # was Ntrees_ini
    n_trees_max: int = 10000  # was Ntrees_max

    # Death rule (Fortran defaults): tree dies if
    #   (n_branches < min_branches AND age > min_age_for_undersize) OR age > max_age
    min_branches: int = 11
    min_age_for_undersize: int = 5
    max_age: int = 1000

    def __post_init__(self) -> None:
        if self.size <= 0.0:
            raise ValueError("ForestConfig.size must be positive")
        if self.n_trees_init <= 0:
            raise ValueError("ForestConfig.n_trees_init must be positive")
        if self.n_trees_max < self.n_trees_init:
            raise ValueError(
                "ForestConfig.n_trees_max must be >= n_trees_init "
                f"({self.n_trees_max} < {self.n_trees_init})"
            )
        if self.min_branches < 1:
            raise ValueError("ForestConfig.min_branches must be >= 1")
        if self.min_age_for_undersize < 0:
            raise ValueError("ForestConfig.min_age_for_undersize must be non-negative")
        if self.max_age <= 0:
            raise ValueError("ForestConfig.max_age must be positive")


@dataclass(frozen=True)
class Config:
    """Top-level config — what ``grow_tree`` and ``Forest`` consume."""

    tree: TreeConfig = field(default_factory=TreeConfig)
    light: LightConfig = field(default_factory=LightConfig)
    forest: ForestConfig = field(default_factory=ForestConfig)
    n_generations: int = 100  # was Ngeneration / Nsteps

    @classmethod
    def from_yaml(cls, path: str | Path) -> Config:
        """Load and validate a config from a YAML file."""
        with open(path) as f:
            data: dict[str, Any] = yaml.safe_load(f) or {}
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Config:
        tree_data = data.get("tree", {}) or {}
        light_data = data.get("light", {}) or {}
        forest_data = data.get("forest", {}) or {}
        # The ``forest`` block carries n_generations in the Fortran .ini; we
        # also accept a top-level ``n_generations`` for users running only the
        # single-tree simulator without a forest block.
        n_gen = data.get("n_generations")
        if n_gen is None:
            n_gen = forest_data.get("n_generations", 100)

        # Filter to known fields so unknown YAML keys don't crash; this keeps
        # the schema forward-compatible with future fields.
        tree_known = {k: v for k, v in tree_data.items() if k in TreeConfig.__dataclass_fields__}
        light_known = {k: v for k, v in light_data.items() if k in LightConfig.__dataclass_fields__}
        forest_known = {
            k: v for k, v in forest_data.items() if k in ForestConfig.__dataclass_fields__
        }
        return cls(
            tree=TreeConfig(**tree_known),
            light=LightConfig(**light_known),
            forest=ForestConfig(**forest_known),
            n_generations=int(n_gen),
        )


def load_config(path: str | Path) -> Config:
    """Convenience entry-point — ``Config.from_yaml(path)``."""
    return Config.from_yaml(path)


__all__ = ["Config", "ForestConfig", "LightConfig", "TreeConfig", "load_config"]
