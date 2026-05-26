"""Tests for the native bulk-thinning wind module (Step 25)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from mechatree.config import Config, ForestConfig, TreeConfig, WindConfig
from mechatree.forest import Forest
from mechatree.wind.bulk_thinning import (
    BulkThinningParams,
    BulkThinningWindBridge,
    _edges_from_centers,
    _stack_trees,
    compute_bulk_thinning,
    make_bulk_thinning_wind_fn,
)


def test_params_defaults_are_uniform_no_boundary_layer():
    p = BulkThinningParams.uniform()
    assert np.allclose(p.U_infty, 1.0)
    assert p.H == 0.5
    assert p.C_D == 1.0
    assert p.U_infty.size == p.z_centers.size


def test_params_validation_rejects_bad_inputs():
    with pytest.raises(ValueError, match="share shape"):
        BulkThinningParams(U_infty=np.array([1.0, 2.0]), z_centers=np.array([0.25]))
    with pytest.raises(ValueError, match="strictly monotone"):
        BulkThinningParams(U_infty=np.array([1.0, 1.0]), z_centers=np.array([0.75, 0.25]))
    with pytest.raises(ValueError, match="H must be positive"):
        BulkThinningParams(U_infty=np.array([1.0, 1.0]), z_centers=np.array([0.25, 0.75]), H=0.0)


def test_edges_from_centers_uniform_grid():
    centers = np.array([0.25, 0.75, 1.25, 1.75])
    edges = _edges_from_centers(centers)
    np.testing.assert_allclose(edges, [0.0, 0.5, 1.0, 1.5, 2.0])


def test_compute_with_empty_canopy_returns_free_stream():
    p = BulkThinningParams.uniform()
    res = compute_bulk_thinning(
        start=np.empty((0, 3), dtype=float),
        axis=np.empty((0, 3), dtype=float),
        D=np.empty(0, dtype=float),
        L=np.empty(0, dtype=float),
        branch_offsets=np.zeros(1, dtype=np.int64),
        params=p,
        wind_direction=(1.0, 0.0),
    )
    assert res.n_branches == 0
    np.testing.assert_allclose(res.U_canopy, p.U_infty)
    assert res.canopy_mean == (0.0, 0.0, 0.0)


def test_compute_x_wind_matches_thinning_formula():
    """Single vertical branch under +x wind; the thinning formula should
    reduce the canopy wind from U_infty to U_infty * 0.5 in the
    over-loaded layers (clipped disc)."""
    p = BulkThinningParams.uniform()
    # One vertical branch at origin: start z=0, axis +z, big diameter so it
    # over-loads its z-layer.
    start = np.array([[0.0, 0.0, 0.0]])
    axis = np.array([[0.0, 0.0, 1.0]])
    D = np.array([1.0])
    L = np.array([1.0])
    res = compute_bulk_thinning(
        start,
        axis,
        D,
        L,
        branch_offsets=np.array([0, 1]),
        params=p,
        wind_direction=(1.0, 0.0),
    )
    # Branch is vertical → sin(I) = 1 → area = D*L = 1.0. With H = 0.5
    # and U_infty = 1: F_drag = 0.5 * 1 = 0.5; disc = 1 - 4*0.5/(0.25*1) = 1-8 < 0 → 0;
    # U_canopy = 0.5 * 1 * (1 + 0) = 0.5.
    # The branch z_c = 0 + 0.5*1*1 = 0.5; with z_edges starting at 0 and step 0.5,
    # z_c=0.5 falls into layer 1.
    assert res.U_branch[0] == pytest.approx(0.5)


def test_compute_rotation_preserves_magnitude():
    """Rotating the wind direction by 45° should give the same canopy-mean
    magnitude (the canopy is the same, only the projection differs)."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=15.0, n_trees_init=5, n_trees_max=50),
    )
    f = Forest(cfg, seed=1)
    for g in range(15):
        f.step(g)

    start, axis, D, L, off = _stack_trees(f.trees)
    params = BulkThinningParams.uniform(U=2.0)

    res_x = compute_bulk_thinning(
        start, axis, D, L, branch_offsets=off, params=params, wind_direction=(1.0, 0.0)
    )
    th = math.pi / 4
    res_rot = compute_bulk_thinning(
        start,
        axis,
        D,
        L,
        branch_offsets=off,
        params=params,
        wind_direction=(math.cos(th), math.sin(th)),
    )

    mag_x = math.hypot(res_x.canopy_mean[0], res_x.canopy_mean[1])
    mag_rot = math.hypot(res_rot.canopy_mean[0], res_rot.canopy_mean[1])
    # On an isotropic forest the projected canopy area depends on wind
    # direction (vertical branches → invariant; oblique branches → not).
    # But for a randomly-oriented canopy the means should match closely.
    assert mag_x == pytest.approx(mag_rot, rel=0.1)

    # Direction is recovered correctly.
    recovered = math.atan2(res_rot.canopy_mean[1], res_rot.canopy_mean[0])
    assert recovered == pytest.approx(th, abs=0.02)


def test_bridge_default_returns_streamwise_wind():
    """No angle sampler → wind blows from +x by convention."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=15.0, n_trees_init=3, n_trees_max=20),
    )
    f = Forest(cfg, seed=0)
    for g in range(10):
        f.step(g)

    bridge = make_bulk_thinning_wind_fn()
    rng = np.random.default_rng(0)
    w = bridge(0, rng, f)
    assert w[1] == 0.0
    assert w[2] == 0.0
    assert w[0] > 0.0


def test_bridge_uses_angle_sampler():
    """With an angle sampler the wind direction matches the sampled
    angle exactly (no random jitter)."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=15.0, n_trees_init=3, n_trees_max=20),
    )
    f = Forest(cfg, seed=0)
    for g in range(10):
        f.step(g)

    target_angle = math.radians(60.0)

    def fixed_angle(rng, n):
        return np.array([target_angle])

    bridge = make_bulk_thinning_wind_fn(angle_sampler=fixed_angle)
    w = bridge(0, np.random.default_rng(0), f)
    recovered = math.atan2(w[1], w[0])
    assert recovered == pytest.approx(target_angle, abs=1e-6)


def test_bridge_handles_empty_forest():
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=15.0, n_trees_init=2, n_trees_max=10),
        wind=WindConfig(model="default"),
    )
    f = Forest(cfg, seed=0)
    f.trees = []
    bridge = make_bulk_thinning_wind_fn()
    w = bridge(0, np.random.default_rng(0), f)
    # Empty canopy → free stream.
    assert w == (float(bridge.params.U_infty[0]), 0.0, 0.0)


def test_yaml_model_native_resolves_to_bridge():
    """WindConfig(model='native') routes _resolve_wind_fn to the native bridge."""
    from mechatree.simulate import _resolve_wind_fn

    cfg = Config(wind=WindConfig(model="native"))
    fn = _resolve_wind_fn(cfg)
    assert isinstance(fn, BulkThinningWindBridge)


def test_amplitude_sampler_scales_canopy_wind():
    """Amplitude sampler multiplies U_infty per call, so the canopy-mean
    output scales linearly with the sampled amplitude."""
    cfg = Config(
        tree=TreeConfig(),
        forest=ForestConfig(size=15.0, n_trees_init=3, n_trees_max=20),
    )
    f = Forest(cfg, seed=0)
    for g in range(10):
        f.step(g)

    def amp_5(rng, n):
        return np.array([5.0])

    def amp_1(rng, n):
        return np.array([1.0])

    bridge_1 = make_bulk_thinning_wind_fn(amplitude_sampler=amp_1)
    bridge_5 = make_bulk_thinning_wind_fn(amplitude_sampler=amp_5)
    w1 = bridge_1(0, np.random.default_rng(0), f)
    w5 = bridge_5(0, np.random.default_rng(0), f)
    # At equal U_infty(z), thinning is the same; scaled U_infty just
    # scales everything linearly.
    assert w5[0] == pytest.approx(5.0 * w1[0], rel=0.05)
