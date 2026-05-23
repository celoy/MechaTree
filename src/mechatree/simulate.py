"""Single-tree simulator — orchestrates the per-generation pipeline.

Mirrors ``legacy_fortran/tree.f90``'s ``time_step`` loop. Per generation:

1. Light: extract leaves -> intercept -> aggregate.
2. Mechanics: 4-angle stress sweep.
3. Growth requests + secondary growth (diameter allocation).
4. Pruning under a per-generation wind direction.
5. Primary growth (new twig pairs at lit leaves).
6. Reserve depletion for the seeds that would drop to ground.
7. Reorder (refresh ``nb_leaves`` for the next iteration).

The orchestrator is intentionally Python — at most ten C++/NumPy calls per
generation, so the dispatch overhead is microseconds against the millisecond
work inside each phase.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np

from mechatree._core import PyTree
from mechatree.config import Config, TreeConfig
from mechatree.genome import AllocationModel, ConstantAllocation, ConstantSafety, SafetyModel
from mechatree.growth import primary_growth, requested_growth, secondary_growth
from mechatree.light import Sun, aggregate_onto_trees, extract_leaves, intercept
from mechatree.mechanics import calculate_stresses
from mechatree.pruning import prune

# A "wind function" maps (generation, rng) -> a 3-tuple wind vector. The
# rng is the same numpy.random.Generator used for everything else in the
# Python orchestrator; the tree's own C++ RNG (used inside primary_growth
# and prune) is seeded separately.
WindFn = Callable[[int, np.random.Generator], tuple[float, float, float]]

OnStep = Callable[[int, PyTree], None]


def default_wind_fn(generation: int, rng: np.random.Generator) -> tuple[float, float, float]:
    """Fortran-faithful wind direction.

    Per ``tree.f90:193-195``::

        angle = 1.0 * generation
        amplitude = 0.835 - log(rand_uniform) / 6.0
        U = (amplitude * cos(angle), amplitude * sin(angle), 0)

    The amplitude has a long tail — most generations see ~1 (the ``0.835``
    bias), but rare gusts can be 2-3x larger.
    """
    angle = float(generation)
    u = float(rng.random())
    # Guard against `u == 0` producing `-inf` from `log(0)`.
    if u <= 0.0:
        u = np.finfo(np.float64).tiny
    amplitude = 0.835 - math.log(u) / 6.0
    return (amplitude * math.cos(angle), amplitude * math.sin(angle), 0.0)


def make_seed_tree(config: TreeConfig) -> PyTree:
    """Construct the initial trunk per ``legacy_fortran/mod_tree.f90`` ``new_tree``.

    location = (0, 0, 0), unit_t = (0, 0, 1), unit_b = (1, 0, 0).
    length = twig_length, diameter = twig_diameter, reserve = 2 * volume_twig.
    """
    tree = PyTree({})
    tree.set_length(0, config.twig_length)
    tree.set_diameter(0, config.twig_diameter)
    tree.set_unit_t(0, (0.0, 0.0, 1.0))
    tree.set_unit_b(0, (1.0, 0.0, 0.0))
    tree.set_location(0, (0.0, 0.0, 0.0))
    tree.set_reserve(2.0 * config.volume_twig)
    return tree


def grow_tree(
    config: Config | TreeConfig,
    *,
    n_generations: int | None = None,
    seed: int = 0,
    safety: SafetyModel | None = None,
    allocation: AllocationModel | None = None,
    wind_fn: WindFn | None = None,
    sun: Sun | None = None,
    on_step: OnStep | None = None,
) -> PyTree:
    """Grow one tree end-to-end.

    Parameters
    ----------
    config:
        Either a full ``Config`` (TreeConfig + LightConfig + n_generations) or
        just a ``TreeConfig``. The TreeConfig form uses default light + the
        ``n_generations`` keyword argument.
    n_generations:
        Number of generations. Required if ``config`` is a ``TreeConfig``;
        overrides ``Config.n_generations`` if given.
    seed:
        Single integer seed. Drives both the tree's C++ RNG (used by
        primary_growth and prune) and the Python ``numpy.random.Generator``
        (used by ``wind_fn``). Reproducible end-to-end.
    safety, allocation:
        Genome models (mechatree.genome). Default = ConstantSafety(1.0) and
        ConstantAllocation(p_seeds=0.1, p_leaves=0.5, phototropism=0.5).
    wind_fn:
        Optional ``(generation, rng) -> (x, y, z)``. Default = ``default_wind_fn``
        (Fortran-faithful rotating wind with long-tailed amplitude).
    sun:
        Optional ``Sun``. Default = ``Sun()`` (4 elevations x 8 azimuths).
    on_step:
        Optional callback ``(generation, tree)`` called after each iteration —
        for plotting, snapshotting, or aggregating statistics. ``None`` (default)
        skips all callback overhead.
    """
    if isinstance(config, Config):
        tree_cfg = config.tree
        if n_generations is None:
            n_generations = config.n_generations
        if sun is None:
            light_cfg = config.light
            sun = Sun(
                n_elevations=light_cfg.n_elevations,
                n_azimuths=light_cfg.n_azimuths,
                size_leaf=light_cfg.size_leaf,
            )
    else:
        tree_cfg = config
        if n_generations is None:
            raise TypeError("grow_tree: n_generations is required when config is a TreeConfig")

    if sun is None:
        sun = Sun()
    if safety is None:
        safety = ConstantSafety(1.0)
    if allocation is None:
        allocation = ConstantAllocation(p_seeds=0.1, p_leaves=0.5, phototropism=0.5)
    if wind_fn is None:
        wind_fn = default_wind_fn

    tree = make_seed_tree(tree_cfg)
    tree.set_seed(seed & 0xFFFFFFFF)
    tree.reorder()
    rng = np.random.default_rng(seed)

    volume_twig = tree_cfg.volume_twig
    volume_per_leaf = tree_cfg.volume_per_leaf

    for gen in range(n_generations):
        # 1. Light.
        leaves = extract_leaves([tree], n_directions=sun.n_directions)
        intercept(leaves, sun)
        aggregate_onto_trees(leaves, [tree])

        # 2. Mechanics.
        calculate_stresses(tree, leaf_drag_S0=tree_cfg.leaf_surface, cauchy=tree_cfg.cauchy)

        # 3. Growth (reads max_stress + nb_leaves; writes vol_growth, diameter).
        requested_growth(tree, safety, maintenance_h=tree_cfg.maintenance_h)
        secondary_growth(tree, volume_per_leaf=volume_per_leaf)

        # 4. Pruning under per-generation wind.
        wind = wind_fn(gen, rng)
        prune(tree, wind=wind, leaf_drag_S0=tree_cfg.leaf_surface, cauchy=tree_cfg.cauchy)
        tree.reorder()

        # 5. Primary growth (new twig pairs). Snapshot reserve BEFORE the call
        # so the seed-cost mirrors the Fortran's pre-call R0.
        reserve_before = tree.get_reserve()
        n_leaves_before = tree.get_total_leaves()
        primary_growth(
            tree,
            allocation,
            twig_length=tree_cfg.twig_length,
            twig_diameter=tree_cfg.twig_diameter,
            theta1=tree_cfg.theta1,
            theta2=tree_cfg.theta2,
            gamma1=tree_cfg.gamma1,
            gamma2=tree_cfg.gamma2,
            generation=gen,
        )

        # 6. Reserve depletion for the seeds primary_growth "drops" (no forest
        # in Step 11, but the energy cost matches Fortran's tree.f90:205).
        if n_leaves_before > 0 and reserve_before > 0.0:
            vol_relative = reserve_before / n_leaves_before / volume_twig
            p_seeds, _, _ = allocation.compute(n_leaves_before, vol_relative)
            n_seeds = int(math.floor(p_seeds * reserve_before / (5.0 * volume_twig)))
            tree.set_reserve(max(0.0, tree.get_reserve() - 5.0 * volume_twig * n_seeds))

        # 7. Reorder for the next iteration's nb_leaves.
        tree.reorder()

        if on_step is not None:
            on_step(gen, tree)

    return tree


__all__ = ["WindFn", "OnStep", "default_wind_fn", "grow_tree", "make_seed_tree"]
