"""Smoke tests for the native momentum-wind wind bridge (Step 25b).

Pure-NumPy kernel — no optional extras required.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import mechatree as mt
from mechatree.config import Config, ForestConfig, TreeConfig, WindConfig
from mechatree.mechanics import calculate_stresses, calculate_stresses_from_stored_forces
from mechatree.simulate import _resolve_wind_fn
from mechatree.wind._momentum_wind_kernel import compute_momentum_wind
from mechatree.wind.momentum_wind import MomentumWindBridge


def _tiny_warmup_forest(seed: int = 0):
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=6.0, n_trees_init=2, n_trees_max=6),
        n_generations=15,
    )
    forest = mt.Forest(cfg, seed=seed)
    for g in range(15):
        forest.step(g)
    return cfg, forest


# ---------------------------------------------------------------------------
# WindConfig validation
# ---------------------------------------------------------------------------


def test_windconfig_momentum_defaults_validate():
    wc = WindConfig(model="momentum")
    assert wc.model == "momentum"
    assert wc.momentum_nu_diff == pytest.approx(0.03)


def test_windconfig_momentum_rejects_nonpositive_nu_diff():
    with pytest.raises(ValueError, match="momentum_nu_diff"):
        WindConfig(model="momentum", momentum_nu_diff=0.0)


def test_windconfig_momentum_accepts_zero_pad():
    wc = WindConfig(model="momentum", momentum_pad_x=0.0)
    assert wc.momentum_pad_x == 0.0


# ---------------------------------------------------------------------------
# Kernel-level tests (no forest, just arrays)
# ---------------------------------------------------------------------------


def test_kernel_empty_canopy_returns_freestream():
    # No branches → U_out everywhere equals U_infty(z).
    grid_size = 1.0
    bx = np.arange(0.0, 5.0 + grid_size, grid_size)
    by = np.arange(0.0, 5.0 + grid_size, grid_size)
    bz = np.arange(0.0, 5.0 + grid_size, grid_size)
    z_centers = 0.5 * (bz[:-1] + bz[1:])
    U_infty = 2.0 + 0.5 * z_centers  # arbitrary increasing
    result = compute_momentum_wind(
        np.empty((0, 3)),
        np.empty((0, 3)),
        np.empty(0),
        np.empty(0),
        cell_bounds_x=bx,
        cell_bounds_y=by,
        cell_bounds_z=bz,
        grid_size=grid_size,
        U_infty=U_infty,
        nu_diff=0.0,  # disable diffusion so the test is exact
    )
    # Every cell sees free stream (no branches anywhere, no diffusion).
    assert np.allclose(result.U_out, U_infty[:, None, None], atol=1e-9)


def test_kernel_single_branch_makes_wake():
    # One vertical branch at the centre. Cells upstream see U_infty;
    # cell containing branch is thinned; downstream cells stay thinned
    # (no diffusion to recover) until the boundary.
    grid_size = 1.0
    bx = np.arange(0.0, 8.0 + grid_size, grid_size)
    by = np.arange(0.0, 4.0 + grid_size, grid_size)
    bz = np.arange(0.0, 4.0 + grid_size, grid_size)
    z_centers = 0.5 * (bz[:-1] + bz[1:])
    U_infty = np.full_like(z_centers, 1.0)

    # Vertical branch at (3.5, 1.5, z=0..1), unit D L.
    start = np.array([[3.5, 1.5, 0.0]])
    axis = np.array([[0.0, 0.0, 1.0]])
    D = np.array([0.5])
    L = np.array([1.0])
    result = compute_momentum_wind(
        start,
        axis,
        D,
        L,
        cell_bounds_x=bx,
        cell_bounds_y=by,
        cell_bounds_z=bz,
        grid_size=grid_size,
        U_infty=U_infty,
        C_D=1.0,
        nu_diff=0.0,
    )
    # Find the cell containing the branch midpoint and confirm
    # downstream wake exists.
    i_branch = int(np.searchsorted(bx, 3.5, side="right") - 1)
    j_branch = int(np.searchsorted(by, 1.5, side="right") - 1)
    k_branch = int(np.searchsorted(bz, 0.5, side="right") - 1)
    # Upstream of branch: free stream.
    assert result.U_out[k_branch, j_branch, i_branch - 1] == pytest.approx(1.0, abs=1e-9)
    # At branch cell: thinned.
    assert result.U_out[k_branch, j_branch, i_branch] < 1.0
    # Downstream: still thinned (no diffusion to recover).
    assert result.U_out[k_branch, j_branch, i_branch + 1] < 1.0


def test_kernel_diffusion_widens_wake():
    """Higher nu_diff ⇒ wake spreads further laterally."""
    grid_size = 1.0
    bx = np.arange(0.0, 12.0 + grid_size, grid_size)
    by = np.arange(0.0, 8.0 + grid_size, grid_size)
    bz = np.arange(0.0, 4.0 + grid_size, grid_size)
    z_centers = 0.5 * (bz[:-1] + bz[1:])
    U_infty = np.full_like(z_centers, 1.0)
    start = np.array([[3.5, 3.5, 0.0]])
    axis = np.array([[0.0, 0.0, 1.0]])
    D = np.array([0.5])
    L = np.array([1.0])

    res_lo = compute_momentum_wind(
        start,
        axis,
        D,
        L,
        cell_bounds_x=bx,
        cell_bounds_y=by,
        cell_bounds_z=bz,
        grid_size=grid_size,
        U_infty=U_infty,
        nu_diff=0.0,
    )
    res_hi = compute_momentum_wind(
        start,
        axis,
        D,
        L,
        cell_bounds_x=bx,
        cell_bounds_y=by,
        cell_bounds_z=bz,
        grid_size=grid_size,
        U_infty=U_infty,
        nu_diff=0.3,
    )
    # Pick a cell far downstream and just off-axis (j=2 vs j_branch=3).
    i_far = 10
    k = 0
    # With more diffusion, the off-axis cell is more affected (lower U).
    assert res_hi.U_out[k, 2, i_far] < res_lo.U_out[k, 2, i_far]


def test_kernel_returns_per_branch_arrays():
    grid_size = 1.0
    bx = np.arange(0.0, 5.0 + grid_size, grid_size)
    by = np.arange(0.0, 5.0 + grid_size, grid_size)
    bz = np.arange(0.0, 5.0 + grid_size, grid_size)
    z_centers = 0.5 * (bz[:-1] + bz[1:])
    U_infty = np.full_like(z_centers, 1.0)
    # Three branches, mixed orientations.
    start = np.array([[2.0, 2.0, 0.5], [2.0, 2.0, 1.5], [3.0, 2.0, 0.5]])
    axis = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    D = np.array([0.2, 0.2, 0.2])
    L = np.array([1.0, 1.0, 1.0])
    result = compute_momentum_wind(
        start,
        axis,
        D,
        L,
        cell_bounds_x=bx,
        cell_bounds_y=by,
        cell_bounds_z=bz,
        grid_size=grid_size,
        U_infty=U_infty,
        nu_diff=0.0,
    )
    # Per-branch arrays have right length.
    assert result.U_branch.shape == (3,)
    assert result.F_N_branch.shape == (3,)
    assert result.F_D_branch.shape == (3,)
    assert result.F_vec_branch.shape == (3, 3)
    # Horizontal branch (axis = +x) has sin(I) = 0 ⇒ F_D = 0 and F_vec = 0.
    assert result.F_D_branch[1] == 0.0
    assert np.allclose(result.F_vec_branch[1], 0.0)
    # Branch 2 is downstream of branch 0 (same y, z; greater x), so it
    # sits in branch 0's wake and sees less wind ⇒ lower F_N.
    assert result.F_N_branch[2] < result.F_N_branch[0]
    # Branch 0 (upstream) sees free stream U_infty = 1.
    assert result.U_branch[0] == pytest.approx(1.0, abs=1e-9)
    # Vector form sanity: for a vertical branch under +x wind, the
    # full force points in +x (along û_perp = (1, 0, 0) - 0·t = +x̂).
    assert result.F_vec_branch[0, 0] > 0.0
    assert abs(result.F_vec_branch[0, 1]) < 1e-12
    assert abs(result.F_vec_branch[0, 2]) < 1e-12
    # And |F| = F_N (the magnitude form matches the scalar magnitude).
    F_vec_mag = np.linalg.norm(result.F_vec_branch[0])
    assert F_vec_mag == pytest.approx(result.F_N_branch[0], rel=1e-9)
    # Streamwise component = F_D.
    assert result.F_vec_branch[0, 0] == pytest.approx(result.F_D_branch[0], rel=1e-9)


# ---------------------------------------------------------------------------
# Bridge-level tests
# ---------------------------------------------------------------------------


def test_bridge_returns_wind_vector_for_default_storm():
    _cfg, forest = _tiny_warmup_forest(seed=1)
    bridge = MomentumWindBridge(grid_size=1.0, pad_x=4.0, pad_y=2.0, pad_z=2.0)
    rng = np.random.default_rng(0)
    wind = bridge(0, rng, forest)
    assert isinstance(wind, tuple) and len(wind) == 3
    assert math.isfinite(wind[0]) and math.isfinite(wind[1])
    assert abs(wind[1]) < 1e-9  # default storm = +x
    assert wind[2] == 0.0
    assert wind[0] > 0.0


def test_bridge_rotation_recovers_storm_direction():
    _cfg, forest = _tiny_warmup_forest(seed=2)
    bridge = MomentumWindBridge(
        grid_size=1.5,
        pad_x=4.0,
        pad_y=2.0,
        pad_z=2.0,
        angle_sampler=lambda rng, n: np.full(n, math.pi / 4),
    )
    rng = np.random.default_rng(3)
    wind = bridge(0, rng, forest)
    mag = math.hypot(wind[0], wind[1])
    assert mag > 0.0
    # 45° storm ⇒ components nearly equal.
    assert abs(wind[0] - wind[1]) / mag < 0.05
    recovered = math.degrees(math.atan2(wind[1], wind[0]))
    assert abs(recovered - 45.0) < 1.0


def test_bridge_empty_forest_short_circuit():
    from types import SimpleNamespace

    fake_forest = SimpleNamespace(trees=[])
    bridge = MomentumWindBridge()
    rng = np.random.default_rng(0)
    wind = bridge(0, rng, fake_forest)  # type: ignore[arg-type]
    assert wind == (0.0, 0.0, 0.0)


def test_bridge_exposes_last_result():
    _cfg, forest = _tiny_warmup_forest(seed=4)
    bridge = MomentumWindBridge(grid_size=1.5, pad_x=4.0, pad_y=2.0, pad_z=2.0)
    rng = np.random.default_rng(0)
    bridge(0, rng, forest)
    assert bridge.last_result is not None
    U_out = bridge.last_result.U_out
    assert U_out.ndim == 3
    # Per-branch arrays present.
    n_branches = sum(t.get_number_of_branches() for t in forest.trees)
    assert bridge.last_result.U_branch.shape == (n_branches,)
    assert bridge.last_result.F_N_branch.shape == (n_branches,)


# ---------------------------------------------------------------------------
# Config + resolver
# ---------------------------------------------------------------------------


def test_set_forces_batch_round_trip():
    """``PyTree.set_forces_batch`` writes the per-branch force vectors
    that the momentum-wind kernel produces, preparing the Option-B
    plumbing back into mechanics/pruning."""
    _cfg, forest = _tiny_warmup_forest(seed=6)
    tree = forest.trees[0]
    n = tree.get_number_of_branches()
    forces = np.column_stack([np.linspace(0, 1, n), np.linspace(0, -1, n), np.zeros(n)])
    tree.set_forces_batch(forces)
    for i in (0, n // 2, n - 1):
        assert tree.get_force(i) == pytest.approx(tuple(forces[i]), rel=1e-12)


def test_set_forces_batch_validates_shape():
    _cfg, forest = _tiny_warmup_forest(seed=7)
    tree = forest.trees[0]
    n = tree.get_number_of_branches()
    with pytest.raises(ValueError, match="length must equal n_branches"):
        tree.set_forces_batch(np.zeros((n + 5, 3)))
    with pytest.raises(ValueError, match=r"shape \(N, 3\)"):
        tree.set_forces_batch(np.zeros((n,)))


def test_resolver_dispatches_momentum():
    _cfg, _ = _tiny_warmup_forest(seed=5)
    cfg = Config(
        tree=_cfg.tree,
        forest=_cfg.forest,
        wind=WindConfig(model="momentum", grid_size=1.0, momentum_pad_x=4.0),
        n_generations=_cfg.n_generations,
    )
    fn = _resolve_wind_fn(cfg)
    assert isinstance(fn, MomentumWindBridge)


# ---------------------------------------------------------------------------
# Uniform inflow override (U_in = K, independent of z)
# ---------------------------------------------------------------------------


def test_bridge_uniform_inflow_is_constant_in_z():
    """With ``U_uniform=K`` the inlet column is flat in z; the log-law
    default produces a z-varying inlet instead."""
    _cfg, forest = _tiny_warmup_forest(seed=8)
    rng = np.random.default_rng(0)

    uniform = MomentumWindBridge(grid_size=1.5, pad_x=4.0, pad_y=2.0, pad_z=2.0, U_uniform=1.60)
    uniform(0, rng, forest)
    inlet_uniform = uniform.last_result.U_in_grid[:, :, 0]
    # Every inlet cell equals K, independent of height.
    assert np.allclose(inlet_uniform, 1.60, atol=1e-12)

    loglaw = MomentumWindBridge(grid_size=1.5, pad_x=4.0, pad_y=2.0, pad_z=2.0)
    loglaw(0, rng, forest)
    inlet_loglaw = loglaw.last_result.U_in_grid[:, :, 0]
    # Log-law inlet varies with height (not flat).
    assert inlet_loglaw[:, 0].std() > 1e-6


def test_make_momentum_wind_fn_threads_u_uniform():
    from mechatree.wind.momentum_wind import make_momentum_wind_fn

    bridge = make_momentum_wind_fn(U_uniform=2.5)
    assert bridge.U_uniform == pytest.approx(2.5)
    # Default keeps the log-law (U_uniform is None).
    assert make_momentum_wind_fn().U_uniform is None


def test_resolver_threads_u_uniform():
    _cfg, _ = _tiny_warmup_forest(seed=9)
    cfg = Config(
        tree=_cfg.tree,
        forest=_cfg.forest,
        wind=WindConfig(
            model="momentum",
            grid_size=1.0,
            momentum_pad_x=4.0,
            momentum_U_uniform=1.60,
        ),
        n_generations=_cfg.n_generations,
    )
    fn = _resolve_wind_fn(cfg)
    assert isinstance(fn, MomentumWindBridge)
    assert fn.U_uniform == pytest.approx(1.60)


# ---------------------------------------------------------------------------
# Step 25c — Option-B per-branch stored forces
# ---------------------------------------------------------------------------


def _segment_forces(tree):
    n = tree.get_number_of_branches()
    return np.array([tree.get_segment_force(i) for i in range(n)])


def _segment_winds(tree):
    n = tree.get_number_of_branches()
    return np.array([tree.get_segment_wind(i) for i in range(n)])


def test_set_segment_forces_batch_round_trip():
    _cfg, forest = _tiny_warmup_forest(seed=6)
    tree = forest.trees[0]
    n = tree.get_number_of_branches()
    forces = np.column_stack([np.linspace(0, 1, n), np.linspace(0, -1, n), np.full(n, 0.3)])
    # The aggregation slot force_ must be left untouched (non-aliasing).
    force_before = tree.get_force(0)
    tree.set_segment_forces_batch(forces)
    got = _segment_forces(tree)
    assert np.allclose(got, forces, atol=1e-12)
    assert tree.get_force(0) == pytest.approx(force_before)


def test_set_segment_winds_batch_round_trip():
    _cfg, forest = _tiny_warmup_forest(seed=6)
    tree = forest.trees[0]
    n = tree.get_number_of_branches()
    winds = np.column_stack([np.full(n, 0.7), np.full(n, -0.2), np.zeros(n)])
    tree.set_segment_winds_batch(winds)
    assert np.allclose(_segment_winds(tree), winds, atol=1e-12)


def test_set_segment_batches_validate_shape():
    _cfg, forest = _tiny_warmup_forest(seed=7)
    tree = forest.trees[0]
    n = tree.get_number_of_branches()
    with pytest.raises(ValueError, match="length must equal n_branches"):
        tree.set_segment_forces_batch(np.zeros((n + 3, 3)))
    with pytest.raises(ValueError, match=r"shape \(N, 3\)"):
        tree.set_segment_forces_batch(np.zeros((n,)))
    with pytest.raises(ValueError, match="length must equal n_branches"):
        tree.set_segment_winds_batch(np.zeros((n + 3, 3)))
    with pytest.raises(ValueError, match=r"shape \(N, 3\)"):
        tree.set_segment_winds_batch(np.zeros((n, 2)))


@pytest.mark.parametrize("mag", [2.0, 4.0, 8.0])
def test_prune_with_stored_forces_matches_prune(mag):
    """Feeding ``prune_with_stored_forces`` the exact segment drag that
    ``wind_force`` computes for a uniform wind (and that uniform wind as the
    per-branch local wind) must reproduce ``prune``'s cuts bit-for-bit."""
    from mechatree.pruning import prune, prune_with_stored_forces

    tc = TreeConfig()
    wind = (mag * 0.8, mag * 0.6, 0.0)
    a = mt.grow_tree(tc, n_generations=25, seed=11)
    b = mt.grow_tree(tc, n_generations=25, seed=11)
    n = a.get_number_of_branches()
    assert n == b.get_number_of_branches()

    forces = np.empty((n, 3))
    winds = np.empty((n, 3))
    for i in range(n):
        f, _m = b.wind_force(i, wind)
        forces[i] = f
        winds[i] = wind
    b.set_segment_forces_batch(forces)
    b.set_segment_winds_batch(winds)

    a.set_seed(999)
    b.set_seed(999)
    cut_a = prune(a, wind, tc.leaf_surface, tc.cauchy)
    cut_b = prune_with_stored_forces(b, tc.leaf_surface, tc.cauchy)
    assert cut_b == cut_a
    assert a.get_number_of_branches() == b.get_number_of_branches()


def test_bridge_always_writes_segment_forces():
    """Step 26a: per-branch is the only momentum behaviour — the bridge
    always advertises the capability and writes non-zero segment forces."""
    _cfg, forest = _tiny_warmup_forest(seed=8)
    bridge = MomentumWindBridge(grid_size=1.5, pad_x=4.0, pad_y=2.0, pad_z=2.0, U_uniform=4.0)
    assert bridge.writes_segment_forces is True
    bridge(0, np.random.default_rng(0), forest)
    total = 0
    nonzero = 0
    for t in forest.trees:
        f = _segment_forces(t)
        total += f.shape[0]
        nonzero += int(np.any(np.abs(f) > 0, axis=1).sum())
    assert total > 0
    assert nonzero > 0  # at least some branches feel drag


def test_bridge_storm_fixed_within_generation():
    """Step 26a: the storm (θ, amplitude) is sampled once per generation and
    held across the fixed-point loop's repeated same-generation calls; a new
    generation re-samples."""
    _cfg, forest = _tiny_warmup_forest(seed=8)
    rng = np.random.default_rng(0)
    bridge = MomentumWindBridge(
        grid_size=2.0,
        pad_x=4.0,
        U_uniform=1.0,
        angle_sampler=lambda r, n: r.uniform(0.0, 2.0 * math.pi, n),
        amplitude_sampler=lambda r, n: r.uniform(1.0, 3.0, n),
    )
    bridge(0, rng, forest)
    theta0, amp0, dir0 = bridge._storm_theta, bridge._storm_amp, bridge.last_wind_direction
    # Same generation (fixed-point inner iterations): storm must not change.
    for _ in range(4):
        bridge(0, rng, forest)
        assert bridge._storm_theta == theta0
        assert bridge._storm_amp == amp0
        assert bridge.last_wind_direction == dir0
    # New generation re-samples (different draw with high probability).
    bridge(1, rng, forest)
    assert bridge._storm_gen == 1
    assert (bridge._storm_theta, bridge._storm_amp) != (theta0, amp0)


def test_bridge_stored_force_rotation_preserves_magnitude():
    """A 45° storm writes forces rotated into the world frame; the magnitude
    matches the kernel's rotated-frame F_vec (rotation is norm-preserving),
    while a 0° storm writes them verbatim."""
    _cfg, forest = _tiny_warmup_forest(seed=8)
    rng = np.random.default_rng(0)

    # 0°: world == rotated frame, stored == raw kernel output.
    zero = MomentumWindBridge(grid_size=1.5, pad_x=4.0, U_uniform=4.0)
    zero(0, rng, forest)
    stored0 = np.concatenate([_segment_forces(t) for t in forest.trees])
    assert np.allclose(stored0, zero.last_result.F_vec_branch, atol=1e-9)

    # 45°: magnitudes preserved, vectors rotated (not identical).
    theta = math.pi / 4
    rot = MomentumWindBridge(
        grid_size=1.5,
        pad_x=4.0,
        U_uniform=4.0,
        angle_sampler=lambda _rng, _n: np.array([theta]),
    )
    rot(0, rng, forest)
    stored45 = np.concatenate([_segment_forces(t) for t in forest.trees])
    raw45 = rot.last_result.F_vec_branch
    assert np.allclose(np.linalg.norm(stored45, axis=1), np.linalg.norm(raw45, axis=1), atol=1e-9)
    # z-component is untouched by the horizontal rotation.
    assert np.allclose(stored45[:, 2], raw45[:, 2], atol=1e-12)


def test_bridge_per_branch_empty_forest_short_circuits():
    class _Empty:
        trees: list = []

    bridge = MomentumWindBridge()
    assert bridge(0, np.random.default_rng(0), _Empty()) == (0.0, 0.0, 0.0)


def test_forest_dispatches_stored_force_prune_end_to_end():
    """A ``model: momentum`` config routes Forest.step through the
    stored-force prune path, runs cleanly, and actually prunes branches via
    that path (U=2.0 is a regime where trees survive but storms still cut)."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=6.0, n_trees_init=3, n_trees_max=8),
        wind=WindConfig(
            model="momentum",
            grid_size=1.5,
            momentum_pad_x=4.0,
            momentum_U_uniform=2.0,
        ),
        n_generations=15,
    )
    forest = mt.Forest(cfg, seed=3)
    assert forest._wind_uses_stored_forces is True
    total_pruned = 0
    for g in range(15):
        total_pruned += forest.step(g).n_pruned_total
    assert len(forest.trees) >= 1
    assert sum(t.get_number_of_branches() for t in forest.trees) > 0
    assert total_pruned > 0  # the stored-force prune path actually cut


def test_resolver_momentum_writes_segment_forces():
    """The resolved momentum-wind bridge always advertises per-branch forces
    (Step 26a — no flag)."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=6.0, n_trees_init=2, n_trees_max=6),
        wind=WindConfig(model="momentum", grid_size=1.0, momentum_pad_x=4.0),
        n_generations=10,
    )
    fn = _resolve_wind_fn(cfg)
    assert isinstance(fn, MomentumWindBridge)
    assert fn.writes_segment_forces is True


# ---------------------------------------------------------------------------
# Step 26c — sensing unified on the momentum per-branch screened field
# ---------------------------------------------------------------------------


def _max_stress(tree):
    return np.array([tree.get_max_stress(i) for i in range(tree.get_number_of_branches())])


def test_stored_force_sensing_matches_legacy_4angle_sweep():
    """Feeding ``calculate_stresses_from_stored_forces`` the exact ``wind_force``
    drag over the legacy 4 angles (at unit wind) reproduces the legacy
    ``calculate_stresses`` per-branch ``max_stress`` bit-for-bit."""
    tc = TreeConfig()
    a = mt.grow_tree(tc, n_generations=25, seed=11)
    b = mt.grow_tree(tc, n_generations=25, seed=11)
    n = a.get_number_of_branches()
    assert n == b.get_number_of_branches()

    calculate_stresses(a, leaf_drag_S0=tc.leaf_surface, cauchy=tc.cauchy)
    ref = _max_stress(a)

    angles = [math.pi * k / 4.0 for k in (1, 2, 3, 4)]  # the legacy hardcoded sweep
    for idx, theta in enumerate(angles):
        u = (math.cos(theta), math.sin(theta), 0.0)
        forces = np.empty((n, 3))
        winds = np.empty((n, 3))
        for i in range(n):
            f, _m = b.wind_force(i, u)
            forces[i] = f
            winds[i] = u
        b.set_segment_forces_batch(forces)
        b.set_segment_winds_batch(winds)
        calculate_stresses_from_stored_forces(
            b, leaf_drag_S0=tc.leaf_surface, cauchy=tc.cauchy, reset_max=(idx == 0)
        )
    got = _max_stress(b)
    assert np.allclose(got, ref, rtol=1e-12, atol=1e-12)


def test_sense_screens_inflow_below_unit():
    """``bridge.sense`` solves at a uniform inlet U=1 but the kernel screens
    it through the canopy — no branch sees more than 1, and interior branches
    see less."""
    _cfg, forest = _tiny_warmup_forest(seed=8)
    bridge = MomentumWindBridge(grid_size=2.0, pad_x=6.0, U_uniform=1.0)
    bridge.sense(forest, 0.0)
    u_branch = bridge.last_result.U_branch
    assert u_branch.max() <= 1.0 + 1e-9  # screening never amplifies above the inlet
    assert u_branch.min() < 1.0 - 1e-6  # some branches are screened (wake/shadow)


def test_sense_writes_max_stress_and_is_deterministic():
    """A sensing sweep sets non-zero ``max_stress`` and is reproducible."""
    tc = TreeConfig()
    s0, cauchy = tc.leaf_surface, tc.cauchy

    def sense_once(seed):
        _cfg, forest = _tiny_warmup_forest(seed=seed)
        bridge = MomentumWindBridge(grid_size=2.0, pad_x=6.0, U_uniform=2.0)
        for k, theta in enumerate(bridge.sensing_angles(np.random.default_rng(0), 4)):
            bridge.sense(forest, theta)
            for t in forest.trees:
                calculate_stresses_from_stored_forces(
                    t, leaf_drag_S0=s0, cauchy=cauchy, reset_max=(k == 0)
                )
        return np.concatenate([_max_stress(t) for t in forest.trees])

    ms = sense_once(8)
    assert ms.size > 0
    assert np.any(ms > 0.0)
    assert np.all(ms >= 0.0)
    assert np.allclose(ms, sense_once(8))  # same seed → identical


def test_sensing_angles_count_and_uniform_default():
    bridge = MomentumWindBridge(grid_size=2.0)
    angles = bridge.sensing_angles(np.random.default_rng(0), 6)
    assert len(angles) == 6
    assert all(0.0 <= a < 2.0 * math.pi + 1e-9 for a in angles)


def test_grow_tree_momentum_sensing_end_to_end():
    """A ``model: momentum`` config drives grow_tree through the screened
    sensing sweep (Step 26c) and is reproducible."""
    cfg = Config(
        tree=TreeConfig(),
        wind=WindConfig(
            model="momentum",
            grid_size=2.0,
            momentum_pad_x=6.0,
            momentum_U_uniform=1.5,
            n_sensing_angles=3,
        ),
        n_generations=20,
    )
    t1 = mt.grow_tree(cfg, seed=4)
    t2 = mt.grow_tree(cfg, seed=4)
    assert t1.get_number_of_branches() > 0
    assert t1.get_number_of_branches() == t2.get_number_of_branches()


def test_forest_momentum_sensing_end_to_end():
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=10.0, n_trees_init=4, n_trees_max=12),
        wind=WindConfig(
            model="momentum",
            grid_size=2.0,
            momentum_pad_x=6.0,
            momentum_U_uniform=1.5,
            n_sensing_angles=3,
        ),
        n_generations=15,
    )
    forest = mt.Forest(cfg, seed=3)
    for g in range(15):
        forest.step(g)
    assert len(forest.trees) >= 1
    assert sum(t.get_number_of_branches() for t in forest.trees) > 0
