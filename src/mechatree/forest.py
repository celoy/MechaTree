"""Forest container — many trees on a disk with light competition (Step 12).

Mirrors ``legacy_fortran/Forest.f90``'s main loop. Cross-tree light
competition falls out of the Step-10 light module operating on the union
of every tree's leaves — no special-case code needed.

Per the design principles in CLAUDE.md ("evolution is external"), every
tree in a single Forest shares the same Safety / Allocation models.
Per-tree genome variation comes later when a real Genome class lands.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field, replace

import numpy as np

from mechatree._core import PyTree
from mechatree.config import Config, ForestConfig, TreeConfig
from mechatree.genome import AllocationModel, SafetyModel, models_from_config
from mechatree.growth import primary_growth, requested_growth, secondary_growth
from mechatree.light import Sun, aggregate_onto_trees, extract_leaves, intercept
from mechatree.mechanics import calculate_stresses
from mechatree.pruning import prune
from mechatree.simulate import _callback_arity, _resolve_wind_fn

# See ``mechatree.simulate.WindFn`` for the two accepted arities.
WindFn = Callable[..., tuple[float, float, float]]
OnStep = Callable[[int, "Forest", "ForestStats"], None]


@dataclass
class ForestStats:
    """Per-step snapshot returned from ``Forest.step``.

    Cheap to compute; useful for tracking self-thinning and biomass.
    """

    generation: int
    n_trees: int
    n_trees_big: int  # trees with > 10 branches
    n_branches_total: int
    n_leaves_total: int
    biomass_total: float  # sum of branch volumes
    n_born: int
    n_died: int


def _make_forest_tree(
    tree_cfg: TreeConfig,
    location: tuple[float, float, float],
    angle: float,
    seed: int,
) -> PyTree:
    """Construct a seed tree at ``location`` with trunk's ``unit_b`` rotated
    by ``angle`` radians around the vertical. Per ``Forest.f90:116``."""
    tree = PyTree({})
    tree.set_length(0, tree_cfg.twig_length)
    tree.set_diameter(0, tree_cfg.twig_diameter)
    tree.set_unit_t(0, (0.0, 0.0, 1.0))
    tree.set_unit_b(0, (math.cos(angle), math.sin(angle), 0.0))
    tree.set_location(0, location)
    tree.set_reserve(2.0 * tree_cfg.volume_twig)
    tree.set_seed(seed & 0xFFFFFFFF)
    tree.reorder()
    return tree


@dataclass
class Forest:
    """A growing forest of trees on a circular plot.

    Construct with ``Forest(config, seed=...)`` and either call ``run(n)``
    or drive the loop manually with ``step(generation)``.
    """

    config: Config
    seed: int = 0
    safety: SafetyModel | None = None
    allocation: AllocationModel | None = None
    sun: Sun | None = None
    wind_fn: WindFn | None = None

    # Populated by __post_init__ / step.
    rng: np.random.Generator = field(init=False, repr=False)
    trees: list[PyTree] = field(init=False, default_factory=list, repr=False)
    ages: list[int] = field(init=False, default_factory=list, repr=False)
    _next_tree_seed: int = field(init=False, default=0, repr=False)
    _wind_arity: int = field(init=False, default=2, repr=False)

    def __post_init__(self) -> None:
        self.rng = np.random.default_rng(self.seed)
        if self.safety is None or self.allocation is None:
            default_safety, default_allocation, default_angles = models_from_config(
                self.config.genome, base_dir=self.config.base_dir
            )
            if self.safety is None:
                self.safety = default_safety
            if self.allocation is None:
                self.allocation = default_allocation
            # YAML ``genome.neural_from`` → use the champion's branching angles
            # (Fortran genome[0..2]) instead of the YAML ``tree:`` defaults.
            # See simulate.grow_tree for the matching code path.
            if default_angles is not None:
                self.config = replace(self.config, tree=replace(self.config.tree, **default_angles))
        if self.wind_fn is None:
            self.wind_fn = _resolve_wind_fn(self.config)
        self._wind_arity = _callback_arity(self.wind_fn)
        if self.sun is None:
            lc = self.config.light
            self.sun = Sun(
                n_elevations=lc.n_elevations,
                n_azimuths=lc.n_azimuths,
                size_leaf=lc.size_leaf,
            )
        self._init_population()

    # ---------- construction --------------------------------------------------

    def _next_seed(self) -> int:
        """Derive a fresh, deterministic C++ RNG seed per tree.

        Mixing in ``self.seed`` keeps two Forests with different master seeds
        independent. The counter is bumped on every new tree (initial seeding
        + every later seedling).
        """
        self._next_tree_seed += 1
        # Cheap LCG-style mix; not security-grade, just want decorrelation.
        return (self.seed * 1_000_003 + self._next_tree_seed * 2_654_435_761) & 0xFFFFFFFF

    def _init_population(self) -> None:
        """Place ``n_trees_init`` trees uniformly across a disk of radius
        ``config.forest.size`` (Fortran ``Forest.f90:178``).

        Uses ``sqrt(rx)`` so the density is uniform across the disk's area,
        not just the circumference.
        """
        fc = self.config.forest
        size = fc.size
        for _ in range(fc.n_trees_init):
            rx = float(self.rng.random())
            ry = float(self.rng.random())
            radius = size * math.sqrt(rx)
            theta = ry * 2.0 * math.pi
            angle = float(self.rng.random()) * 2.0 * math.pi
            location = (radius * math.cos(theta), radius * math.sin(theta), 0.0)
            tree = _make_forest_tree(self.config.tree, location, angle, self._next_seed())
            self.trees.append(tree)
            self.ages.append(0)

    # ---------- run loop ------------------------------------------------------

    def step(self, generation: int) -> ForestStats:
        """One generation of the full pipeline.

        Order: light (union of leaves) -> per-tree mechanics + growth ->
        pruning under a common wind -> death -> seedling birth.
        """
        tree_cfg = self.config.tree
        forest_cfg = self.config.forest
        volume_twig = tree_cfg.volume_twig
        volume_per_leaf = tree_cfg.volume_per_leaf

        # 1. Light across the union of leaves — cross-tree competition.
        if self.trees:
            leaves = extract_leaves(self.trees, n_directions=self.sun.n_directions)
            intercept(leaves, self.sun, leaf_transparency=self.config.light.leaf_transparency)
            aggregate_onto_trees(leaves, self.trees)

        # 2. Per-tree mechanics + growth.
        for tree in self.trees:
            calculate_stresses(tree, leaf_drag_S0=tree_cfg.leaf_surface, cauchy=tree_cfg.cauchy)
            requested_growth(tree, self.safety, maintenance_h=tree_cfg.maintenance_h)
            secondary_growth(tree, volume_per_leaf=volume_per_leaf)

        # 3. Pruning under a common wind.
        if self._wind_arity >= 3:
            wind = self.wind_fn(generation, self.rng, self)
        else:
            wind = self.wind_fn(generation, self.rng)
        for tree in self.trees:
            prune(tree, wind=wind, leaf_drag_S0=tree_cfg.leaf_surface, cauchy=tree_cfg.cauchy)
            tree.reorder()
            # Optional: fuse single-child parent->child chains left by pruning
            # into one straight segment (bottom/top + total volume kept).
            # ``collapse_chains_after_prune`` only walks the chains seeded by
            # this generation's cuts (recorded by ``prune`` itself), so the
            # cost is proportional to the number of cuts rather than to the
            # tree's size. Worth trying on long forest runs to keep per-step
            # cost down.
            #
            #   n_pruned = prune(...)
            #   tree.reorder()
            #   if n_pruned > 0:
            #       tree.collapse_chains_after_prune()  # length_max=10.0 by default
            #       tree.reorder()

        # 4. Age + death.
        n_died = 0
        survivors_trees: list[PyTree] = []
        survivors_ages: list[int] = []
        for tree, age in zip(self.trees, self.ages, strict=True):
            new_age = age + 1
            if self._is_dead(tree, new_age, forest_cfg):
                n_died += 1
                continue
            survivors_trees.append(tree)
            survivors_ages.append(new_age)
        self.trees = survivors_trees
        self.ages = survivors_ages

        # 5. Primary growth + seed dispersal.
        n_born = self._grow_and_disperse(generation, volume_twig)

        # 6. Stats.
        n_branches_total = sum(t.get_number_of_branches() for t in self.trees)
        n_leaves_total = sum(t.get_total_leaves() for t in self.trees)
        n_big = sum(1 for t in self.trees if t.get_number_of_branches() > 10)
        biomass = self._biomass()
        return ForestStats(
            generation=generation,
            n_trees=len(self.trees),
            n_trees_big=n_big,
            n_branches_total=n_branches_total,
            n_leaves_total=n_leaves_total,
            biomass_total=biomass,
            n_born=n_born,
            n_died=n_died,
        )

    def run(self, n_generations: int, on_step: OnStep | None = None) -> None:
        """Drive the simulation for ``n_generations`` steps."""
        for gen in range(n_generations):
            stats = self.step(gen)
            if on_step is not None:
                on_step(gen, self, stats)

    # ---------- internal helpers ---------------------------------------------

    def _is_dead(self, tree: PyTree, age: int, forest_cfg: ForestConfig) -> bool:
        """Fortran ``Forest.f90:283``::

        n_branches < min_branches AND age > min_age_for_undersize
          OR age > max_age
        """
        n_branches = tree.get_number_of_branches()
        undersize_and_old = (
            n_branches < forest_cfg.min_branches and age > forest_cfg.min_age_for_undersize
        )
        return undersize_and_old or age > forest_cfg.max_age

    def _grow_and_disperse(self, generation: int, volume_twig: float) -> int:
        """Run primary_growth on every tree, then disperse N_seeds seedlings.

        Mirrors ``Forest.f90:304-329``. The seed count comes from the same
        allocation-model formula used in Step 11.
        """
        tree_cfg = self.config.tree
        forest_cfg = self.config.forest
        size_sq = forest_cfg.size**2

        # Snapshot reserves BEFORE primary_growth so the seed count uses the
        # pre-call R0, matching the Fortran's energy book-keeping.
        seedlings: list[tuple[float, float, float]] = []
        for tree in self.trees:
            reserve_before = tree.get_reserve()
            n_leaves_before = tree.get_total_leaves()

            primary_growth(
                tree,
                self.allocation,
                twig_length=tree_cfg.twig_length,
                twig_diameter=tree_cfg.twig_diameter,
                theta1=tree_cfg.theta1,
                theta2=tree_cfg.theta2,
                gamma1=tree_cfg.gamma1,
                gamma2=tree_cfg.gamma2,
                generation=generation,
            )

            if n_leaves_before <= 0 or reserve_before <= 0.0:
                tree.reorder()
                continue

            vol_relative = reserve_before / n_leaves_before / volume_twig
            p_seeds, _, _ = self.allocation.compute(n_leaves_before, vol_relative)
            n_seeds = int(math.floor(p_seeds * reserve_before / (5.0 * volume_twig)))
            tree.set_reserve(max(0.0, tree.get_reserve() - 5.0 * volume_twig * n_seeds))
            tree.reorder()

            # Cap by remaining global headroom.
            remaining = forest_cfg.n_trees_max - len(self.trees) - len(seedlings)
            n_seeds = max(0, min(n_seeds, remaining))
            if n_seeds == 0:
                continue

            # Pick random leaves; compute landing positions.
            leaf_idxs = tree.leaf_indices()
            if not leaf_idxs:
                continue
            picks = self.rng.integers(0, len(leaf_idxs), size=n_seeds)
            flight_angles = self.rng.random(size=(n_seeds, 2)) * 2.0 * math.pi
            for j in range(n_seeds):
                leaf_idx = leaf_idxs[int(picks[j])]
                base = tree.get_location(leaf_idx)
                # Fortran: dispersal radius = leaf.location.z (the branch base
                # height). Taller leaves throw seeds further.
                dx = base[2] * math.cos(flight_angles[j, 0])
                dy = base[2] * math.sin(flight_angles[j, 0])
                lx = base[0] + dx
                ly = base[1] + dy
                if lx * lx + ly * ly < size_sq:
                    seedlings.append((lx, ly, float(flight_angles[j, 1])))

        # Spawn the surviving seedlings.
        n_born = 0
        for lx, ly, ang in seedlings:
            if len(self.trees) >= forest_cfg.n_trees_max:
                break
            tree = _make_forest_tree(
                tree_cfg, location=(lx, ly, 0.0), angle=ang, seed=self._next_seed()
            )
            self.trees.append(tree)
            self.ages.append(0)
            n_born += 1
        return n_born

    def _biomass(self) -> float:
        """Sum of branch volumes (``length * pi * diameter^2 / 4``) across all
        trees — mirrors ``mod_tree.f90:163`` (``biomass_calculation``)."""
        total = 0.0
        for tree in self.trees:
            n = tree.get_number_of_branches()
            for i in range(n):
                d = tree.get_diameter(i)
                total += tree.get_length(i) * 0.25 * math.pi * d * d
        return total


__all__ = ["Forest", "ForestStats", "WindFn", "OnStep"]
