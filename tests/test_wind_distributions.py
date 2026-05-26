"""Tests for the Step-25 SymPy-driven wind distributions."""

from __future__ import annotations

import math

import numpy as np
import pytest

from mechatree.wind.distributions import (
    Distribution,
    default_amplitude_distribution,
    default_amplitude_sampler,
    legacy_angle_sampler,
    uniform_angle_distribution,
    uniform_angle_sampler,
)


def test_distribution_validation():
    with pytest.raises(ValueError, match="non-empty"):
        Distribution(cdf_expr="", var_name="a", support=(0, 1))
    with pytest.raises(ValueError, match="identifier"):
        Distribution(cdf_expr="x", var_name="not a name", support=(0, 1))
    with pytest.raises(ValueError, match="lo < hi"):
        Distribution(cdf_expr="x", var_name="x", support=(1, 1))
    with pytest.raises(ValueError, match=">= 16"):
        Distribution(cdf_expr="x", var_name="x", support=(0, 1), n_grid=4)


def test_distribution_rejects_extra_free_symbols():
    d = Distribution(cdf_expr="a + b", var_name="a", support=(0.0, math.inf))
    with pytest.raises(ValueError, match="unexpected free symbols"):
        d.sampler()


def test_default_amplitude_sampler_mean_matches_theory():
    """Shifted exponential with rate 6, shift 0.835 → mean 0.835 + 1/6."""
    s = default_amplitude_sampler()
    rng = np.random.default_rng(0)
    samples = s(rng, 50_000)
    assert samples.mean() == pytest.approx(0.835 + 1 / 6, rel=0.02)
    assert samples.min() >= 0.835 - 1e-9


def test_default_amplitude_via_sympy_matches_in_distribution():
    """Pure-NumPy fast path and SymPy CDF-inversion path produce the
    same distribution (sorted samples should match within Monte-Carlo
    noise)."""
    n = 50_000
    fast = default_amplitude_sampler()(np.random.default_rng(0), n)
    sympy = default_amplitude_distribution().sampler()(np.random.default_rng(0), n)
    # Means agree closely.
    assert abs(fast.mean() - sympy.mean()) < 0.01
    # Quartile-by-quartile agreement.
    for q in (0.25, 0.5, 0.75, 0.95):
        assert np.quantile(fast, q) == pytest.approx(np.quantile(sympy, q), rel=0.02)


def test_uniform_angle_sampler_covers_support():
    s = uniform_angle_sampler()
    rng = np.random.default_rng(7)
    samples = s(rng, 10_000)
    assert samples.min() >= 0.0
    assert samples.max() <= 2 * math.pi
    assert samples.mean() == pytest.approx(math.pi, rel=0.02)


def test_uniform_angle_via_sympy_works():
    d = uniform_angle_distribution()
    samples = d.sampler()(np.random.default_rng(7), 10_000)
    assert samples.min() >= 0.0 - 1e-9
    assert samples.max() <= 2 * math.pi + 1e-9
    assert samples.mean() == pytest.approx(math.pi, rel=0.02)


def test_rayleigh_via_numerical_fallback():
    """A Rayleigh CDF: SymPy symbolic inversion may return the negative
    root; the support-check rejects it and the sampler falls back to
    the numerical lookup."""
    d = Distribution(cdf_expr="1 - exp(-x**2 / 2)", var_name="x", support=(0.0, math.inf))
    samples = d.sampler()(np.random.default_rng(0), 50_000)
    assert samples.min() >= 0.0
    assert samples.mean() == pytest.approx(math.sqrt(math.pi / 2), rel=0.02)


def test_numerical_fallback_for_piecewise_cdf():
    """Piecewise CDFs that SymPy can't invert symbolically should go
    through the np.interp lookup path."""
    d = Distribution(
        cdf_expr="Piecewise((x/2, x < 1), ((x + 1)/4, True))",
        var_name="x",
        support=(0.0, 3.0),
    )
    samples = d.sampler()(np.random.default_rng(0), 20_000)
    assert samples.min() >= 0.0 - 1e-9
    assert samples.max() <= 3.0 + 1e-9


def test_legacy_angle_sampler_returns_fixed_four():
    s = legacy_angle_sampler()
    out = s(np.random.default_rng(0), 4)
    assert out.shape == (4,)
    np.testing.assert_allclose(out, [math.pi / 4, math.pi / 2, 3 * math.pi / 4, math.pi])


def test_uniform_angle_sampler_custom_range():
    s = uniform_angle_sampler(lo=0.0, hi=math.pi / 2)
    samples = s(np.random.default_rng(0), 5_000)
    assert samples.min() >= 0.0
    assert samples.max() <= math.pi / 2 + 1e-9
    assert samples.mean() == pytest.approx(math.pi / 4, rel=0.05)
    with pytest.raises(ValueError, match="hi > lo"):
        uniform_angle_sampler(lo=1.0, hi=0.5)
