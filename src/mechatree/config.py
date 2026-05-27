"""YAML configuration for the single-tree simulator (Step 11).

Direct translation of the parameters in ``legacy/fortran/Forest.ini``, regrouped
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
    ``legacy/fortran/Forest.ini``.
    """

    # Direct YAML inputs --------------------------------------------------
    twig_length: float = 1.0  # was TwigLength
    twig_diameter: float = 0.1  # was TwigDiameter
    leaf_surface: float = 0.25  # was LeafSurface (S0)
    # was Cauchy (Cy). Doubled from the Forest.ini runtime value 2.0e-5: the
    # ½ρU² dynamic-pressure factor that the Fortran omitted from the drag was
    # restored in the C++ mechanics (mechanics.cpp wind_force + the leaf term),
    # halving the whole wind load, so 2× Cauchy reproduces the Nat Commun
    # reference stresses bit-for-bit.
    cauchy: float = 4.0e-5
    # was VolumeRatioLeaf; matches the Nat Commun reference V_prod = 4 V_0 l
    volume_ratio_leaf: float = 4.0
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
    """Parameters for the hemispherical light integration.

    ``leaf_transparency`` (``tau``) sets the per-leaf optical opacity used by
    :func:`mechatree.light.interception.intercept`: the i-th leaf from the
    top of a shadow cell receives ``tau**i`` of the incident light.
    ``tau = 0`` recovers the Fortran binary topmost-wins regime; ``tau = 1``
    makes leaves fully transparent; default ``tau = 0.5`` matches the value
    used in Eloy et al. (Nat Commun 2017).
    """

    size_leaf: float = 1.0  # was SizeLeaf — 2D shadow-grid cell size
    n_elevations: int = 4  # was Nelev
    n_azimuths: int = 8  # was Nazim
    leaf_transparency: float = 0.5  # tau in [0, 1]; 0 = binary, 1 = transparent

    def __post_init__(self) -> None:
        if self.n_elevations <= 0:
            raise ValueError("LightConfig.n_elevations must be positive")
        if self.n_azimuths <= 0:
            raise ValueError("LightConfig.n_azimuths must be positive")
        if self.size_leaf <= 0.0:
            raise ValueError("LightConfig.size_leaf must be positive")
        if not (0.0 <= self.leaf_transparency <= 1.0):
            raise ValueError(
                f"LightConfig.leaf_transparency must be in [0, 1], got {self.leaf_transparency}"
            )


@dataclass(frozen=True)
class GenomeConfig:
    """Genome inputs for the simulator.

    The scalar fields (``safety``, ``p_seeds``, ``p_leaves``, ``phototropism``)
    populate :class:`ConstantSafety` / :class:`ConstantAllocation` when no
    neural genome is supplied. Their defaults reproduce the values the Fortran
    ``neural_branch`` / ``neural_reserve`` networks evolved to in Eloy et al.
    (Nat Commun 2017): ``safety = 3`` puts branches at ~(1/3)^(3/2) ≈ 0.19 of
    breaking stress, the value the evolved network settled on.

    ``neural_from``, when set, points at a champion JSON written by
    ``scripts/strategies_single_tree.py`` and selects which species to load.
    When it is present, :class:`NeuralSafety` / :class:`NeuralAllocation` are
    used instead of the scalar fields.

    Step 15: each scalar field may also be a **string** holding a SymPy
    expression in the canonical input names (``nb_leaves`` and ``max_stress``
    for safety; ``nb_leaves`` and ``vol_relative`` for the allocation
    fields). When at least one field is a string, :mod:`mechatree.sympy_genome`
    compiles all four into :class:`CallbackSafety` / :class:`CallbackAllocation`.
    """

    safety: float | str = 3.0  # was neural_branch output (Safety)
    p_seeds: float | str = 0.1  # was neural_reserve output [0]
    p_leaves: float | str = 0.5  # was neural_reserve output [1]
    phototropism: float | str = 0.5  # was neural_reserve output [2]
    neural_from: dict[str, Any] | None = None  # {"path": ..., "species_id": 0}

    def __post_init__(self) -> None:
        # Scalar checks only apply when the field is a number; string
        # expressions are validated at model-build time in
        # :func:`mechatree.sympy_genome._compile`.
        if isinstance(self.safety, (int, float)) and self.safety <= 0.0:
            raise ValueError(f"GenomeConfig.safety must be positive, got {self.safety}")
        if isinstance(self.p_seeds, (int, float)) and self.p_seeds < 0.0:
            raise ValueError(f"GenomeConfig.p_seeds must be non-negative, got {self.p_seeds}")
        if isinstance(self.p_leaves, (int, float)) and self.p_leaves < 0.0:
            raise ValueError(f"GenomeConfig.p_leaves must be non-negative, got {self.p_leaves}")
        if (
            isinstance(self.p_seeds, (int, float))
            and isinstance(self.p_leaves, (int, float))
            and self.p_seeds + self.p_leaves > 1.0
        ):
            raise ValueError(
                "GenomeConfig.p_seeds + p_leaves must be <= 1, "
                f"got {self.p_seeds} + {self.p_leaves}"
            )
        if isinstance(self.phototropism, (int, float)) and not (0.0 <= self.phototropism <= 1.0):
            raise ValueError(
                f"GenomeConfig.phototropism must be in [0, 1], got {self.phototropism}"
            )
        # String fields must be non-empty so they fail loudly rather than
        # silently turning into ``"0"`` at lambdify time.
        for field_name in ("safety", "p_seeds", "p_leaves", "phototropism"):
            v = getattr(self, field_name)
            if isinstance(v, str) and not v.strip():
                raise ValueError(f"GenomeConfig.{field_name}: empty string expression")
        if self.neural_from is not None:
            if not isinstance(self.neural_from, dict):
                raise ValueError(
                    "GenomeConfig.neural_from must be a dict like {'path': ..., 'species_id': 0}"
                )
            if "path" not in self.neural_from:
                raise ValueError("GenomeConfig.neural_from requires a 'path' key")
            if not isinstance(self.neural_from["path"], str):
                raise ValueError(
                    f"GenomeConfig.neural_from['path'] must be a string, "
                    f"got {type(self.neural_from['path']).__name__}"
                )
            species_id = self.neural_from.get("species_id", 0)
            if not isinstance(species_id, int) or isinstance(species_id, bool):
                raise ValueError(
                    f"GenomeConfig.neural_from['species_id'] must be an int, "
                    f"got {type(species_id).__name__}"
                )


@dataclass(frozen=True)
class WindConfig:
    """Wind model selection + parameters.

    Two models, both valid for the top-level ``wind:`` YAML block:

    - ``model: default`` keeps the Fortran-faithful storm wind
      (:func:`mechatree.simulate.default_wind_fn`) — a single, spatially
      uniform direction per generation, no canopy feedback (no screening).
    - ``model: momentum`` selects the native 3-D momentum-wind CFD solver
      in :mod:`mechatree.wind.momentum_wind` — the canonical canopy-aware
      (screening) path. Builds a structured grid from the forest bounding
      box (cell size ``grid_size``, padded by ``momentum_pad_x/y/z``), runs
      the simplified Eloy-style per-branch drag kernel + per-cell momentum
      balance + cross-stream diffusion (single dimensionless knob
      ``momentum_nu_diff``). Inflow is a log-law parameterised by
      ``momentum_ua / z0 / kappa`` by default; set ``momentum_U_uniform``
      for a height-independent constant inflow ``U_in = K`` instead
      (``ua / z0 / kappa`` then ignored). Each branch is scored against its
      own local screened wind in pruning (Step 26a — there is no
      canopy-mean wind mode). See ``docs/momentum_wind_derivation.md``.
      (Step 26: the legacy ``native`` bulk-thinning and ``dendroflow``
      bridges were removed — both collapsed the field to a scalar, which is
      just a constant wind; ``momentum`` supersedes them.)

    Step 24 knobs (canopy-aware wind only):

    - ``max_pruning_iterations`` / ``wind_convergence_eps_rel`` control
      the fixed-point loop that re-evaluates wind after each pruning
      sweep. Exits on the first of zero new cuts, sub-``eps_rel``
      canopy-mean change, or hitting the cap. ``cap=1`` recovers the
      old single-pass behaviour; ``eps_rel=0`` disables the early-exit.

    Step 25 knobs (storm statistics, all wind models):

    - ``amplitude_cdf`` (SymPy CDF in ``a``) drives the per-generation
      storm magnitude in the default wind, and scales the inflow per
      generation in the ``momentum`` model. ``None`` keeps the legacy
      Fortran formula (``a = 0.835 - log(U)/6`` for default; ``a = 1``
      constant scaling for ``momentum``).
    - ``angle_cdf`` (SymPy CDF in ``theta``) drives the storm direction.
      ``None`` keeps the legacy behaviour (``angle = generation rad``
      for default; storm ``+x`` for ``momentum`` unless overridden).
    - ``n_sensing_angles`` (Step 26c) is the number of directions the
      ``momentum`` sensing sweep evaluates each generation: for each, the
      momentum field is solved at a uniform inlet ``U = 1`` (still screened
      through the canopy) and the per-branch max stress over the directions
      drives growth/safety. The directions are drawn from ``angle_cdf``
      (uniform ``[0, 2π)`` when unset). Only used by ``model: momentum``;
      the ``default`` model keeps the legacy hardcoded 4-angle unit-wind
      sweep with no screening.
    """

    model: str = "default"
    # Grid cell size for the momentum model (Step 26b: renamed from the old
    # ``H``). The 3-D cell edge. Default 2 is the O(branches)≈O(cells) cost
    # knee that still resolves the vertical + intra-crown wind gradient.
    grid_size: float = 2.0
    C_D: float = 1.0
    max_pruning_iterations: int = 8
    wind_convergence_eps_rel: float = 0.01
    # Step 25: tunable storm distributions. SymPy CDFs in the documented
    # variable name ('a' for amplitude, 'theta' for angle). ``None``
    # selects the legacy hard-coded sampler so existing YAMLs without
    # these fields keep their behaviour.
    amplitude_cdf: str | None = None
    angle_cdf: str | None = None
    # Step 26c/26d: number of screened directions the momentum sensing sweep
    # evaluates per generation, each a fresh random draw from angle_cdf. 2 is
    # the default (each ~halves the sensing wall-clock vs 4); over generations
    # the random redraw still integrates the directional load.
    n_sensing_angles: int = 2
    # Momentum-wind CFD parameters (only used when model='momentum').
    # Defaults: log-law inflow with ua=0.4, z0=0.1, kappa=0.41; cell
    # padding around the canopy bounding box; single dimensionless
    # cross-stream diffusion knob nu_diff = 0.03 (matches the
    # historical k-ε defaults — see docs/momentum_wind_derivation.md).
    momentum_pad_x: float = 12.0
    momentum_pad_y: float = 2.0
    momentum_pad_z: float = 3.0
    momentum_ua: float = 0.4
    momentum_z0: float = 0.1
    momentum_kappa: float = 0.41
    momentum_nu_diff: float = 0.03
    # Uniform inflow override: when set, the momentum-wind model uses a
    # height-independent constant inflow U_in = momentum_U_uniform
    # and the log-law (ua/z0/kappa) is ignored. ``None`` keeps the log-law.
    momentum_U_uniform: float | None = None

    def __post_init__(self) -> None:
        if self.max_pruning_iterations < 1:
            raise ValueError(
                f"WindConfig.max_pruning_iterations must be >= 1, got {self.max_pruning_iterations}"
            )
        if self.wind_convergence_eps_rel < 0.0:
            raise ValueError(
                "WindConfig.wind_convergence_eps_rel must be non-negative, "
                f"got {self.wind_convergence_eps_rel}"
            )
        if self.n_sensing_angles < 1:
            raise ValueError(
                f"WindConfig.n_sensing_angles must be >= 1, got {self.n_sensing_angles}"
            )
        if self.amplitude_cdf is not None and not (
            isinstance(self.amplitude_cdf, str) and self.amplitude_cdf.strip()
        ):
            raise ValueError(
                "WindConfig.amplitude_cdf must be a non-empty SymPy expression in 'a', or None"
            )
        if self.angle_cdf is not None and not (
            isinstance(self.angle_cdf, str) and self.angle_cdf.strip()
        ):
            raise ValueError(
                "WindConfig.angle_cdf must be a non-empty SymPy expression in 'theta', or None"
            )
        if self.model not in ("default", "momentum"):
            raise ValueError(
                f"WindConfig.model must be one of 'default'/'momentum', got {self.model!r}"
            )
        if self.model == "momentum":
            if self.grid_size <= 0.0:
                raise ValueError(f"WindConfig.grid_size must be positive, got {self.grid_size}")
            if self.C_D <= 0.0:
                raise ValueError(f"WindConfig.C_D must be positive, got {self.C_D}")
            for name, value in (
                ("momentum_pad_x", self.momentum_pad_x),
                ("momentum_pad_y", self.momentum_pad_y),
                ("momentum_pad_z", self.momentum_pad_z),
            ):
                if value < 0.0:
                    raise ValueError(f"WindConfig.{name} must be non-negative, got {value}")
            for name, value in (
                ("momentum_ua", self.momentum_ua),
                ("momentum_z0", self.momentum_z0),
                ("momentum_kappa", self.momentum_kappa),
                ("momentum_nu_diff", self.momentum_nu_diff),
            ):
                if value <= 0.0:
                    raise ValueError(f"WindConfig.{name} must be positive, got {value}")
            if self.momentum_U_uniform is not None and self.momentum_U_uniform <= 0.0:
                raise ValueError(
                    "WindConfig.momentum_U_uniform must be positive when set, "
                    f"got {self.momentum_U_uniform}"
                )


@dataclass(frozen=True)
class ForestConfig:
    """Parameters for a forest of trees (Step 12).

    ``n_trees_init`` trees are seeded uniformly across a disk of radius
    ``size``. The death rule mirrors ``legacy/fortran/Forest.f90:283`` and
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
    genome: GenomeConfig = field(default_factory=GenomeConfig)
    wind: WindConfig = field(default_factory=WindConfig)
    n_generations: int = 100  # was Ngeneration / Nsteps
    # Directory used to resolve relative paths inside the config (e.g.
    # ``genome.neural_from.path``). Populated by ``from_yaml`` with the YAML
    # file's parent dir; ``None`` for programmatically-built configs.
    # ``compare=False, hash=False`` keeps it out of value equality so two
    # configs from different paths but identical content still compare equal.
    base_dir: Path | None = field(default=None, compare=False, hash=False)

    @classmethod
    def from_yaml(cls, path: str | Path) -> Config:
        """Load and validate a config from a YAML file."""
        path = Path(path)
        with open(path) as f:
            data: dict[str, Any] = yaml.safe_load(f) or {}
        return cls.from_dict(data, base_dir=path.parent)

    @classmethod
    def from_dict(cls, data: dict[str, Any], *, base_dir: Path | None = None) -> Config:
        tree_data = data.get("tree", {}) or {}
        light_data = data.get("light", {}) or {}
        forest_data = data.get("forest", {}) or {}
        genome_data = data.get("genome", {}) or {}
        wind_data = data.get("wind", {}) or {}
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
        genome_known = {
            k: v for k, v in genome_data.items() if k in GenomeConfig.__dataclass_fields__
        }
        wind_known = {k: v for k, v in wind_data.items() if k in WindConfig.__dataclass_fields__}
        # Tuple-ify list fields so the frozen dataclass stays hashable.
        for tuple_field in ("U_infty", "z_centers"):
            if tuple_field in wind_known and wind_known[tuple_field] is not None:
                wind_known[tuple_field] = tuple(float(x) for x in wind_known[tuple_field])
        return cls(
            tree=TreeConfig(**tree_known),
            light=LightConfig(**light_known),
            forest=ForestConfig(**forest_known),
            genome=GenomeConfig(**genome_known),
            wind=WindConfig(**wind_known),
            n_generations=int(n_gen),
            base_dir=base_dir,
        )


def load_config(path: str | Path) -> Config:
    """Convenience entry-point — ``Config.from_yaml(path)``."""
    return Config.from_yaml(path)


__all__ = [
    "Config",
    "ForestConfig",
    "GenomeConfig",
    "LightConfig",
    "TreeConfig",
    "WindConfig",
    "load_config",
]
