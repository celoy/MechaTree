"""Single-tree simulator — orchestrates the per-generation pipeline.

Mirrors ``legacy/fortran/tree.f90``'s ``time_step`` loop. Per generation:

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

import inspect
import math
import warnings
from collections.abc import Callable
from dataclasses import dataclass, replace

import numpy as np

from mechatree._core import PyTree
from mechatree.config import Config, LightConfig, TreeConfig, WindConfig
from mechatree.genome import AllocationModel, SafetyModel, models_from_config
from mechatree.growth import primary_growth, requested_growth, secondary_growth
from mechatree.light import Sun, aggregate_onto_trees, extract_leaves, intercept
from mechatree.mechanics import calculate_stresses, calculate_stresses_from_stored_forces
from mechatree.pruning import prune, prune_with_stored_forces

# A "wind function" maps (generation, rng[, context]) -> a 3-tuple wind vector.
# Two arities are accepted (arity detected at call time, same pattern as
# ``OnStep`` below):
#   - ``cb(generation, rng)`` — classic Step-11 shape; ignores the tree state.
#   - ``cb(generation, rng, context)`` — receives the live ``PyTree`` (in
#     ``grow_tree``) or the ``Forest`` (in ``Forest.step``) as the third arg.
#     Used by canopy-aware wind models (e.g. the momentum-wind bridge in
#     ``mechatree.wind.momentum_wind``).
# The ``rng`` is the same numpy.random.Generator used elsewhere in the Python
# orchestrator; the tree's own C++ RNG (used inside primary_growth and prune)
# is seeded separately.
WindFn = Callable[..., tuple[float, float, float]]

# ``on_step`` callback. Two forms are accepted for backward compatibility:
#   - ``cb(generation, tree)`` — original Step-11 signature
#   - ``cb(generation, tree, stats)`` — receives per-step bookkeeping
# Arity is detected automatically.
OnStep = Callable[..., None]


@dataclass
class TreeStats:
    """Per-step bookkeeping for ``grow_tree`` — mirrors the columns of the
    Fortran ``ZAllocation.dat`` file (``mod_tools.f90`` ``save_allocation``).

    Exposed via the 3-arg form of the ``on_step`` callback.
    """

    generation: int
    n_branches: int
    n_leaves: int
    wind: tuple[float, float, float]
    wind_amplitude: float
    n_twigs_created: int  # branches added by primary_growth (always even)
    n_seeds: int  # seeds "dropped to ground" — pays into the reserve
    n_pruned: int  # branches removed by prune (including subtrees)
    reserve: float
    # Step 24: passes through the wind → prune fixed-point loop. Always 1
    # for the default 2-arg wind (the canopy doesn't feed back) and for
    # canopy-aware winds on calm generations. Climbs into 2-4 on storm
    # generations when the surviving canopy materially changes the wind.
    n_wind_iterations: int = 1


def _sense_canopy(
    context,
    trees: list,
    wind_fn: WindFn,
    wind_uses_stored_forces: bool,
    n_sensing_angles: int,
    rng: np.random.Generator,
    *,
    leaf_drag_S0: float,
    cauchy: float,
) -> None:
    """Phase-2 sensing — write per-branch ``max_stress`` for the growth law.

    Step 26c: for the momentum model (``wind_uses_stored_forces``), sweep
    ``n_sensing_angles`` directions sampled from the storm angle distribution.
    Each is a momentum solve at a **uniform inlet U = 1** (the wind scale,
    still screened through the canopy), after which the per-branch screened
    forces drive a stress pass; the per-branch max over directions is kept.
    Sensing and pruning then see the *same* screened field — a sheltered
    branch reinforces against the weaker wind it actually feels.

    For every other wind model (``default`` / user 2-arg callables) this is
    the legacy per-tree uniform 4-angle :func:`calculate_stresses`, leaving
    those paths byte-identical (no extra RNG draws)."""
    if wind_uses_stored_forces and trees:
        angles = wind_fn.sensing_angles(rng, n_sensing_angles)
        if hasattr(wind_fn, "solve_directions"):
            # Step 26e: solve the independent directions in parallel (the
            # GIL-free C++ kernel lets a thread pool overlap them), then write
            # each angle's screened forces + run the stress pass sequentially so
            # the per-branch max-stress reduction stays deterministic.
            _trees, counts, per_angle = wind_fn.solve_directions(context, angles)
            split_at = np.cumsum(counts)[:-1] if counts else np.array([], dtype=int)
            for k, (f_world, w_world) in enumerate(per_angle):
                f_per = np.split(f_world, split_at)
                w_per = np.split(w_world, split_at)
                for tree, ft, wt in zip(trees, f_per, w_per, strict=True):
                    tree.set_segment_forces_batch(ft)
                    tree.set_segment_winds_batch(wt)
                for tree in trees:
                    calculate_stresses_from_stored_forces(
                        tree, leaf_drag_S0=leaf_drag_S0, cauchy=cauchy, reset_max=(k == 0)
                    )
        else:
            for k, theta in enumerate(angles):
                wind_fn.sense(context, theta)
                for tree in trees:
                    calculate_stresses_from_stored_forces(
                        tree, leaf_drag_S0=leaf_drag_S0, cauchy=cauchy, reset_max=(k == 0)
                    )
    else:
        for tree in trees:
            calculate_stresses(tree, leaf_drag_S0=leaf_drag_S0, cauchy=cauchy)


def _prune_to_fixed_point(
    tree: PyTree,
    wind_fn: WindFn,
    wind_arity: int,
    gen: int,
    rng: np.random.Generator,
    *,
    leaf_drag_S0: float,
    cauchy: float,
    max_iterations: int,
    eps_rel: float,
    use_stored_forces: bool = False,
) -> tuple[tuple[float, float, float], int, int, bool]:
    """One generation's wind → prune fixed-point loop. Step 24.

    Returns ``(last_wind, n_pruned_total, n_iterations, cap_hit)``.

    For the default 2-arg ``wind_fn`` (wind independent of canopy), this
    is exactly the Step-11 single-pass behaviour — bit-identical. For a
    3-arg ``wind_fn`` (canopy-aware, e.g. the momentum-wind bridge), the
    loop iterates ``wind_fn → prune`` until one of:

    * the most recent ``prune`` cut nothing (the true fixed point), or
    * the recomputed wind differs from the previous one by less than
      ``eps_rel`` (the canopy barely moved, so iterating again would
      change nothing meaningful — the cheap early exit), or
    * the iteration count hits ``max_iterations`` (safety cap; the
      caller is expected to emit a warning once per run).

    When ``use_stored_forces`` (Step 25c, option B — the momentum-wind
    bridge wrote per-branch forces during the ``wind_fn`` call), each
    sweep uses ``prune_with_stored_forces`` so branches are scored against
    their own local CFD wind. ``wind`` is still tracked for the ε exit.
    """
    canopy_aware = wind_arity >= 3

    def _prune(t: PyTree) -> int:
        if use_stored_forces:
            return prune_with_stored_forces(t, leaf_drag_S0=leaf_drag_S0, cauchy=cauchy)
        return prune(t, wind=wind, leaf_drag_S0=leaf_drag_S0, cauchy=cauchy)

    wind = wind_fn(gen, rng, tree) if canopy_aware else wind_fn(gen, rng)
    n_cut = _prune(tree)
    tree.reorder()
    n_pruned_total = n_cut
    n_iter = 1
    if not canopy_aware or n_cut == 0:
        return wind, n_pruned_total, n_iter, False

    while n_iter < max_iterations:
        prev_wind = wind
        wind = wind_fn(gen, rng, tree)
        n_iter += 1
        if eps_rel > 0.0:
            delta = math.hypot(wind[0] - prev_wind[0], wind[1] - prev_wind[1])
            ref = max(math.hypot(prev_wind[0], prev_wind[1]), 1e-6)
            if delta / ref < eps_rel:
                return wind, n_pruned_total, n_iter, False
        n_cut = _prune(tree)
        tree.reorder()
        n_pruned_total += n_cut
        if n_cut == 0:
            return wind, n_pruned_total, n_iter, False

    # We exhausted ``max_iterations`` and the last sweep still cut at
    # least one branch (the early-exit and zero-cut returns above would
    # have triggered otherwise).
    return wind, n_pruned_total, n_iter, True


def default_wind_fn(generation: int, rng: np.random.Generator) -> tuple[float, float, float]:
    """Fortran-faithful wind direction (legacy, kept byte-identical).

    Per ``tree.f90:193-195``::

        angle = 1.0 * generation
        amplitude = 0.835 - log(rand_uniform) / 6.0
        U = (amplitude * cos(angle), amplitude * sin(angle), 0)

    The amplitude has a long tail — most generations see ~1 (the ``0.835``
    bias), but rare gusts can be 2-3x larger.

    For tunable amplitude / angle distributions (Step 25), build a
    storm-wind closure with :func:`make_default_wind_fn` instead.
    """
    angle = float(generation)
    u = float(rng.random())
    # Guard against `u == 0` producing `-inf` from `log(0)`.
    if u <= 0.0:
        u = np.finfo(np.float64).tiny
    amplitude = 0.835 - math.log(u) / 6.0
    return (amplitude * math.cos(angle), amplitude * math.sin(angle), 0.0)


def make_default_wind_fn(
    *,
    amplitude_sampler: Callable[[np.random.Generator, int], np.ndarray] | None = None,
    angle_sampler: Callable[[np.random.Generator, int], np.ndarray] | None = None,
) -> WindFn:
    """Step 25: build a configurable storm-wind closure.

    Returns a 2-arg ``WindFn(generation, rng) -> (x, y, z)`` that samples
    one amplitude and one angle per call. The resulting wind vector is
    ``(a * cos θ, a * sin θ, 0)`` — same shape as the legacy
    :func:`default_wind_fn`.

    ``amplitude_sampler`` and ``angle_sampler`` are callables
    ``(rng, n) -> ndarray`` of length ``n``. Build them from
    :mod:`mechatree.wind.distributions` (e.g.
    :func:`~mechatree.wind.distributions.default_amplitude_sampler`,
    :func:`~mechatree.wind.distributions.uniform_angle_sampler`) or from
    any user-supplied :class:`~mechatree.wind.distributions.Distribution`
    via ``.sampler()``. When either is ``None`` the closure reproduces
    the legacy Fortran behaviour for that axis (amplitude:
    ``0.835 - log(U)/6``; angle: ``generation`` rad).
    """

    def _wind(generation: int, rng: np.random.Generator) -> tuple[float, float, float]:
        if amplitude_sampler is not None:
            amplitude = float(amplitude_sampler(rng, 1)[0])
        else:
            u = float(rng.random())
            if u <= 0.0:
                u = np.finfo(np.float64).tiny
            amplitude = 0.835 - math.log(u) / 6.0
        angle = float(angle_sampler(rng, 1)[0]) if angle_sampler is not None else float(generation)
        return (amplitude * math.cos(angle), amplitude * math.sin(angle), 0.0)

    return _wind


def make_seed_tree(config: TreeConfig) -> PyTree:
    """Construct the initial trunk per ``legacy/fortran/mod_tree.f90`` ``new_tree``.

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
        Genome models (mechatree.genome). When ``None``, defaults are built
        from ``config.genome`` (a ``GenomeConfig``) — i.e. ``ConstantSafety``
        and ``ConstantAllocation`` populated from YAML. If ``config`` is a
        bare ``TreeConfig`` (no genome block), built-in defaults are used:
        ``ConstantSafety(3.0)`` and
        ``ConstantAllocation(p_seeds=0.1, p_leaves=0.5, phototropism=0.5)``.
    wind_fn:
        Optional wind callable. Two arities are accepted: ``(generation, rng)``
        — classic shape, ignored tree state — or ``(generation, rng, tree)``,
        which receives the live ``PyTree`` as the third argument (used by the
        momentum-wind bridge). When ``None``, the wind callable is resolved from
        ``config.wind``: ``model: default`` (Fortran-faithful rotating wind) or
        ``model: momentum`` (canopy-aware screening).
    sun:
        Optional ``Sun``. Default = ``Sun()`` (4 elevations x 8 azimuths).
    on_step:
        Optional callback ``(generation, tree)`` called after each iteration —
        for plotting, snapshotting, or aggregating statistics. ``None`` (default)
        skips all callback overhead.
    """
    genome_cfg = None
    base_dir = None
    leaf_transparency = LightConfig().leaf_transparency
    if isinstance(config, Config):
        tree_cfg = config.tree
        genome_cfg = config.genome
        base_dir = config.base_dir
        leaf_transparency = config.light.leaf_transparency
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
    if safety is None or allocation is None:
        default_safety, default_allocation, default_angles = models_from_config(
            genome_cfg, base_dir=base_dir
        )
        if safety is None:
            safety = default_safety
        if allocation is None:
            allocation = default_allocation
        # When the YAML config selects a Neural champion via ``genome.neural_from``,
        # ``models_from_config`` also returns the champion's branching angles.
        # Apply them to ``tree_cfg`` so the simulator uses the champion's actual
        # geometry (Fortran genome[0..2]) rather than the YAML's tree.theta*/gamma*
        # defaults. Without this override, the YAML defaults would silently
        # shadow the champion's encoded geometry.
        if default_angles is not None:
            tree_cfg = replace(tree_cfg, **default_angles)
    if wind_fn is None:
        wind_fn = _resolve_wind_fn(config)

    # Pull Step-24 fixed-point knobs from the wind config when available;
    # otherwise use ``WindConfig`` defaults. The 2-arg default wind makes
    # both fields inert anyway, so this is mainly about preserving the
    # bare-``TreeConfig`` entry-point.
    wind_cfg = config.wind if isinstance(config, Config) else WindConfig()
    max_pruning_iterations = wind_cfg.max_pruning_iterations
    wind_eps_rel = wind_cfg.wind_convergence_eps_rel

    cb_arity = _callback_arity(on_step) if on_step is not None else 0
    wind_arity = _callback_arity(wind_fn)
    # Step 25c (option B): the momentum-wind bridge advertises that it writes
    # per-branch forces during its call; when so, prune from those instead of
    # the canopy-mean. Absent on every other wind fn → False → unchanged.
    wind_uses_stored_forces = getattr(wind_fn, "writes_segment_forces", False)

    tree = make_seed_tree(tree_cfg)
    tree.set_seed(seed & 0xFFFFFFFF)
    tree.reorder()
    rng = np.random.default_rng(seed)

    volume_twig = tree_cfg.volume_twig
    volume_per_leaf = tree_cfg.volume_per_leaf
    cap_hit_warned = False

    for gen in range(n_generations):
        # 1. Light.
        leaves = extract_leaves([tree], n_directions=sun.n_directions)
        intercept(leaves, sun, leaf_transparency=leaf_transparency)
        aggregate_onto_trees(leaves, [tree])

        # 2. Mechanics — sensing. For the momentum model this sweeps
        # ``n_sensing_angles`` screened directions (Step 26c); otherwise the
        # legacy uniform 4-angle stress sweep.
        _sense_canopy(
            tree,
            [tree],
            wind_fn,
            wind_uses_stored_forces,
            wind_cfg.n_sensing_angles,
            rng,
            leaf_drag_S0=tree_cfg.leaf_surface,
            cauchy=tree_cfg.cauchy,
        )

        # 3. Growth (reads max_stress + nb_leaves; writes vol_growth, diameter).
        requested_growth(tree, safety, maintenance_h=tree_cfg.maintenance_h)
        secondary_growth(tree, volume_per_leaf=volume_per_leaf)

        # 4. Pruning under per-generation wind. Step 24: when ``wind_fn``
        # is canopy-aware (3-arg), iterate ``wind_fn → prune`` to a
        # fixed point so the surviving tree is scored against the wind
        # it would actually feel post-pruning. For the 2-arg default
        # wind this is exactly a single pass — byte-identical to before.
        wind, n_pruned, n_wind_iter, cap_hit = _prune_to_fixed_point(
            tree,
            wind_fn,
            wind_arity,
            gen,
            rng,
            leaf_drag_S0=tree_cfg.leaf_surface,
            cauchy=tree_cfg.cauchy,
            max_iterations=max_pruning_iterations,
            eps_rel=wind_eps_rel,
            use_stored_forces=wind_uses_stored_forces,
        )
        # cap=1 is the explicit "single-pass mode" knob (recovers the
        # pre-Step-24 behaviour for A/B comparisons), not a cap violation
        # — never warn for it.
        if cap_hit and not cap_hit_warned and max_pruning_iterations > 1:
            warnings.warn(
                f"Wind/pruning fixed-point loop hit cap "
                f"({max_pruning_iterations}) at generation {gen}; "
                "consider raising WindConfig.max_pruning_iterations.",
                RuntimeWarning,
                stacklevel=2,
            )
            cap_hit_warned = True
        # Optional: fuse single-child parent->child chains produced by pruning
        # into one straight segment (bottom/top points + total volume kept;
        # the merged diameter is back-solved from the new length). Pruning is
        # the only event that creates such chains, so guard on n_pruned > 0
        # to skip the no-op pass. The "after_prune" variant only walks the
        # handful of chains seeded by this generation's cuts — much cheaper
        # than the whole-tree ``collapse_single_child_chains``. Useful for
        # long forest simulations where leftover chains accumulate; perturbs
        # the mechanics slightly because merged segments have larger moment
        # arms than the bends they replace.
        # if n_pruned > 0:
        #     tree.collapse_chains_after_prune()  # length_max=10.0 by default
        #     tree.reorder()

        # 5. Primary growth (new twig pairs). Snapshot reserve BEFORE the call
        # so the seed-cost mirrors the Fortran's pre-call R0.
        reserve_before = tree.get_reserve()
        n_leaves_before = tree.get_total_leaves()
        n_twigs = primary_growth(
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

        # 6. Reserve depletion for the seeds primary_growth "drops".
        n_seeds = 0
        if n_leaves_before > 0 and reserve_before > 0.0:
            vol_relative = reserve_before / n_leaves_before / volume_twig
            p_seeds, _, _ = allocation.compute(n_leaves_before, vol_relative)
            n_seeds = int(math.floor(p_seeds * reserve_before / (5.0 * volume_twig)))
            tree.set_reserve(max(0.0, tree.get_reserve() - 5.0 * volume_twig * n_seeds))

        # 7. Reorder for the next iteration's nb_leaves.
        tree.reorder()

        if on_step is not None:
            if cb_arity >= 3:
                stats = TreeStats(
                    generation=gen,
                    n_branches=tree.get_number_of_branches(),
                    n_leaves=tree.get_total_leaves(),
                    wind=wind,
                    wind_amplitude=math.hypot(wind[0], wind[1]),
                    n_twigs_created=n_twigs,
                    n_seeds=n_seeds,
                    n_pruned=n_pruned,
                    reserve=tree.get_reserve(),
                    n_wind_iterations=n_wind_iter,
                )
                on_step(gen, tree, stats)
            else:
                on_step(gen, tree)

    return tree


def _build_storm_samplers(wind_cfg) -> tuple:
    """Resolve ``WindConfig.amplitude_cdf`` / ``angle_cdf`` to samplers.

    Returns ``(amplitude_sampler, angle_sampler)``, either of which may
    be ``None`` if the YAML left the field unset (then the legacy
    Fortran defaults are used by ``make_default_wind_fn`` and friends).
    """
    from mechatree.wind.distributions import Distribution

    amp = None
    if wind_cfg.amplitude_cdf is not None:
        amp = Distribution(
            cdf_expr=wind_cfg.amplitude_cdf,
            var_name="a",
            # Open-ended upper bound; the numerical-fallback path
            # auto-extends until the CDF saturates.
            support=(0.0, math.inf),
        ).sampler()
    ang = None
    if wind_cfg.angle_cdf is not None:
        ang = Distribution(
            cdf_expr=wind_cfg.angle_cdf,
            var_name="theta",
            support=(0.0, 2.0 * math.pi),
        ).sampler()
    return amp, ang


def _resolve_wind_fn(config: Config | TreeConfig) -> WindFn:
    """Pick a wind callable based on the YAML ``wind:`` block.

    For bare ``TreeConfig`` (no wind block at all) returns the legacy
    :func:`default_wind_fn` byte-identically. For ``Config``, dispatches
    on ``wind.model``:

    - ``"default"`` (or no model field): builds a storm-wind closure
      via :func:`make_default_wind_fn`, threading in the Step-25
      ``amplitude_cdf`` / ``angle_cdf`` samplers when set. If neither
      is set the closure reproduces :func:`default_wind_fn` exactly.
    - ``"momentum"``: the native 3-D momentum-wind CFD bridge (the only
      canopy-aware / screening model; Step 26 removed the legacy
      ``native`` and ``dendroflow`` bridges).
    """
    if not isinstance(config, Config):
        return default_wind_fn

    wc = config.wind
    amp_sampler, angle_sampler = _build_storm_samplers(wc)

    if wc.model == "momentum":
        from mechatree.wind.momentum_wind import MomentumWindBridge

        return MomentumWindBridge(
            grid_size=wc.grid_size,
            nu_diff=wc.momentum_nu_diff,
            pad_x=wc.momentum_pad_x,
            pad_y=wc.momentum_pad_y,
            pad_z=wc.momentum_pad_z,
            ua=wc.momentum_ua,
            z0=wc.momentum_z0,
            kappa=wc.momentum_kappa,
            U_uniform=wc.momentum_U_uniform,
            C_D=wc.C_D,
            sensing_threads=wc.momentum_sensing_threads,
            angle_sampler=angle_sampler,
            amplitude_sampler=amp_sampler,
        )

    # ``default`` — Fortran-faithful storm wind, plus optional tunable
    # samplers. When both samplers are None we hand back the legacy
    # `default_wind_fn` so existing seeded runs stay byte-identical.
    if amp_sampler is None and angle_sampler is None:
        return default_wind_fn
    return make_default_wind_fn(
        amplitude_sampler=amp_sampler,
        angle_sampler=angle_sampler,
    )


def _callback_arity(fn: Callable) -> int:
    """Number of positional args the callback expects.

    Used to decide whether to pass ``TreeStats`` to a step callback —
    keeps the Step-11 ``cb(gen, tree)`` signature working while letting
    new code opt into ``cb(gen, tree, stats)``.
    """
    try:
        sig = inspect.signature(fn)
    except (TypeError, ValueError):
        return 2  # builtin / C function — assume the original signature
    return sum(
        1
        for p in sig.parameters.values()
        if p.kind
        in (
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
        )
    )


__all__ = [
    "OnStep",
    "TreeStats",
    "WindFn",
    "default_wind_fn",
    "grow_tree",
    "make_seed_tree",
]
