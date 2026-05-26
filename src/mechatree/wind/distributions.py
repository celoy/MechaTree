"""SymPy-driven sampling from user-supplied CDFs (Step 25).

Two distributions plug into the simulator:

- **Amplitude distribution** drives the per-generation storm magnitude in
  :func:`mechatree.simulate.default_wind_fn`. The default reproduces the
  Fortran formula ``a = 0.835 - log(U)/6`` (a shifted exponential),
  expressed as its CDF ``F(a) = 1 - exp(-6 * (a - 0.835))``.
- **Angle distribution** drives both the sensing sweep (``angle_samples``
  angles per generation fed into the multi-angle stress calculation in
  ``calculate_stresses``) and the storm direction. Default uniform on
  ``[0, 2π)``.

The CDF is the user-facing primitive. SymPy attempts to invert it
symbolically for inverse-transform sampling; if the symbolic solve fails
or returns a non-real branch, we fall back to a precomputed numerical
lookup table over the support range.

Public API
----------
:class:`Distribution` — a frozen description (CDF expression + variable
name + support) that compiles to a ``Callable[[Generator, int], ndarray]``
on demand via :meth:`Distribution.sampler`.

Built-in factories that match the Fortran/uniform defaults without
requiring SymPy at import time:

- :func:`default_amplitude_distribution` — shifted exponential, mirrors the
  Fortran ``0.835 - log(U)/6``.
- :func:`uniform_angle_distribution` — uniform on ``[0, 2π)``.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

# Plain-Python fast paths so the import-time dependency on SymPy stays
# optional. Built-in distributions never touch SymPy.

Sampler = Callable[[np.random.Generator, int], np.ndarray]


@dataclass(frozen=True)
class Distribution:
    """A 1-D continuous distribution defined by its CDF.

    Parameters
    ----------
    cdf_expr
        SymPy-parseable expression in the variable ``var_name``. Must be
        a monotone non-decreasing function mapping the support range
        into ``[0, 1]``. Examples: ``"a / (2*pi)"`` (uniform on
        ``[0, 2*pi)``), ``"1 - exp(-6 * (a - 0.835))"`` (shifted
        exponential).
    var_name
        Name of the free variable in ``cdf_expr`` (``"a"`` for
        amplitudes, ``"theta"`` for angles, etc.).
    support
        ``(lo, hi)`` numeric range. ``hi`` may be ``math.inf``.
    n_grid
        Grid resolution for the numerical inverse-CDF fallback when
        SymPy cannot invert ``cdf_expr`` symbolically. Higher = more
        accurate samples at extreme tails, at the cost of larger lookup
        tables. Default 1024 is plenty for typical wind distributions.
    """

    cdf_expr: str
    var_name: str
    support: tuple[float, float]
    n_grid: int = 1024

    def __post_init__(self) -> None:
        if not isinstance(self.cdf_expr, str) or not self.cdf_expr.strip():
            raise ValueError("Distribution.cdf_expr must be a non-empty string")
        if not isinstance(self.var_name, str) or not self.var_name.isidentifier():
            raise ValueError(
                f"Distribution.var_name must be a Python identifier, got {self.var_name!r}"
            )
        lo, hi = self.support
        if not (lo < hi):
            raise ValueError(f"Distribution.support must satisfy lo < hi, got {self.support}")
        if self.n_grid < 16:
            raise ValueError(f"Distribution.n_grid must be >= 16, got {self.n_grid}")

    def sampler(self) -> Sampler:
        """Compile and return a ``sample(rng, n) -> ndarray`` callable.

        Lazy on SymPy: only imports it the first time it is called.
        """
        return _build_sampler(self.cdf_expr, self.var_name, self.support, self.n_grid)


def _build_sampler(
    cdf_expr: str,
    var_name: str,
    support: tuple[float, float],
    n_grid: int,
) -> Sampler:
    try:
        import sympy as sp
    except ImportError as err:  # pragma: no cover - exercised only without the extra.
        raise ImportError(
            "mechatree.wind.distributions requires SymPy when the user supplies "
            "a custom CDF. Install with `pip install 'mechatree[sympy]'` or use "
            "one of the built-in distributions (default_amplitude_distribution / "
            "uniform_angle_distribution) that don't require it."
        ) from err

    var = sp.Symbol(var_name, real=True)
    # Parse with a local_dict so the user's ``var_name`` resolves to the
    # same Symbol object we just built — without this, the
    # ``real=True``-flagged Symbol fails to equal the un-flagged Symbol
    # that ``sympify`` would otherwise create from the bare identifier.
    try:
        expr = sp.sympify(cdf_expr, locals={var_name: var})
    except (sp.SympifyError, SyntaxError) as err:
        raise ValueError(f"Distribution.cdf_expr could not be parsed: {cdf_expr!r}") from err

    free = expr.free_symbols
    extra = {s for s in free if s.name != var_name}
    if extra:
        raise ValueError(
            f"Distribution.cdf_expr {cdf_expr!r} has unexpected free symbols "
            f"{sorted(str(s) for s in extra)}; only {var_name!r} is allowed."
        )

    # Symbolic-inversion attempt. Each candidate is round-tripped through
    # the original CDF on a dense grid so we reject solutions that look
    # plausible at the endpoints but pick the wrong complex branch in
    # between (which would silently emit NaN-floored / clipped samples).
    u_sym = sp.Symbol("__u__", real=True)
    sampler_fn: Callable | None = None
    try:
        sols = sp.solve(expr - u_sym, var)
    except (NotImplementedError, ValueError, TypeError):
        sols = []

    cdf_fn_validate = sp.lambdify(var, expr, modules="numpy")
    u_probe = np.linspace(0.01, 0.99, 25)
    lo_s, hi_s = support

    for sol in sols:
        if sol.is_real is False:
            continue
        try:
            lam = sp.lambdify(u_sym, sol, modules="numpy")
        except Exception:  # noqa: BLE001 - lambdify can raise opaque errors
            continue
        # Round-trip validation: F(F^{-1}(u)) ≈ u over the support, AND
        # the chosen branch must stay inside the declared support (so we
        # don't pick e.g. the negative root for a positive-support CDF).
        try:
            with np.errstate(all="ignore"):
                a_probe = np.asarray(lam(u_probe), dtype=complex)
                if np.any(np.abs(a_probe.imag) > 1e-9):
                    continue
                a_probe = a_probe.real
                u_back = np.asarray(cdf_fn_validate(a_probe), dtype=float)
        except Exception:  # noqa: BLE001
            continue
        if not np.all(np.isfinite(u_back)):
            continue
        if np.max(np.abs(u_back - u_probe)) > 1e-6:
            continue
        # Branch must lie inside the support — rejects the negative root
        # of e.g. sqrt(-2*log(1-u)) for a Rayleigh on [0, inf).
        hi_check = hi_s if math.isfinite(hi_s) else 1e18
        if np.any(a_probe < lo_s - 1e-9) or np.any(a_probe > hi_check + 1e-9):
            continue
        sampler_fn = lam
        break

    if sampler_fn is not None:
        lo, hi = support

        def _sample(rng: np.random.Generator, n: int) -> np.ndarray:
            u = rng.random(n)
            with np.errstate(all="ignore"):
                out_complex = np.asarray(sampler_fn(u), dtype=complex)
            out = out_complex.real
            # Defensive clip: keep samples inside the declared support
            # even if numerical noise pushes them a hair outside.
            return np.clip(out, lo, hi if math.isfinite(hi) else lo + 1e18)

        return _sample

    # Numerical fallback: precompute (a_grid, cdf_grid), invert via np.interp.
    lo, hi = support
    if not math.isfinite(hi):
        cdf_fn = sp.lambdify(var, expr, modules="numpy")
        hi = lo + 1.0
        # Double the upper bound until the CDF saturates close to 1.
        for _ in range(60):
            if float(cdf_fn(hi)) >= 1.0 - 1e-6:
                break
            hi *= 2.0
        else:
            raise ValueError(
                f"Could not bracket CDF {cdf_expr!r} on {support}; F(hi) "
                "did not approach 1 after 60 doublings."
            )
    cdf_fn = sp.lambdify(var, expr, modules="numpy")
    grid = np.linspace(lo, hi, n_grid)
    cdf_vals = np.asarray(cdf_fn(grid), dtype=float)
    if not np.all(np.diff(cdf_vals) >= -1e-9):
        raise ValueError(
            f"CDF {cdf_expr!r} is not monotone non-decreasing on {support}; "
            "inverse-CDF sampling requires a monotone CDF."
        )

    def _sample_interp(rng: np.random.Generator, n: int) -> np.ndarray:
        u = rng.random(n)
        return np.asarray(np.interp(u, cdf_vals, grid), dtype=float)

    return _sample_interp


# ---------------------------------------------------------------------------
# Built-in defaults — pure-NumPy fast paths so the simulator's defaults
# don't pull in SymPy.
# ---------------------------------------------------------------------------


def default_amplitude_distribution() -> Distribution:
    """The Fortran-faithful storm amplitude distribution.

    Shifted exponential with mode at ``0.835`` and rate ``6``. CDF::

        F(a) = 1 - exp(-6 * (a - 0.835))   for a >= 0.835

    Mean amplitude ≈ 1.0; rare tails reach ~3 (typical 1-in-1000 storm).
    Inverse: ``a = 0.835 - log(1 - u) / 6``, which equals the original
    Fortran formula ``a = 0.835 - log(u) / 6`` in distribution.
    """
    return Distribution(
        cdf_expr="1 - exp(-6 * (a - 0.835))",
        var_name="a",
        support=(0.835, math.inf),
    )


def default_amplitude_sampler() -> Sampler:
    """Pure-NumPy sampler that reproduces the Fortran amplitude formula.

    Used as the simulator's hot-path default so the import-time SymPy
    dependency is optional. Mathematically equivalent to
    ``default_amplitude_distribution().sampler()`` to within rounding.
    """

    def _sample(rng: np.random.Generator, n: int) -> np.ndarray:
        # Guard u == 0 from -log(0) producing +inf.
        u = rng.random(n)
        u = np.where(u <= 0.0, np.finfo(np.float64).tiny, u)
        return 0.835 - np.log(u) / 6.0

    return _sample


def uniform_angle_distribution(lo: float = 0.0, hi: float = 2.0 * math.pi) -> Distribution:
    """Uniform angle distribution on ``[lo, hi)``. Default ``[0, 2π)``."""
    width = hi - lo
    if width <= 0:
        raise ValueError(f"uniform_angle_distribution requires hi > lo, got ({lo}, {hi})")
    return Distribution(
        cdf_expr=f"(theta - {lo}) / {width}",
        var_name="theta",
        support=(lo, hi),
    )


def uniform_angle_sampler(lo: float = 0.0, hi: float = 2.0 * math.pi) -> Sampler:
    """Pure-NumPy uniform-angle sampler; bypasses SymPy for the default path."""
    if hi <= lo:
        raise ValueError(f"uniform_angle_sampler requires hi > lo, got ({lo}, {hi})")
    width = hi - lo

    def _sample(rng: np.random.Generator, n: int) -> np.ndarray:
        return lo + width * rng.random(n)

    return _sample


def legacy_angle_sampler() -> Sampler:
    """Legacy 4-fixed-angle sampling pattern from the Fortran code.

    Returns the same four angles ``(π/4, π/2, 3π/4, π)`` regardless of
    ``rng``, ignoring ``n`` and always returning four entries. Useful
    only when the user opts into byte-identical reproduction of the
    pre-Step-25 sensing sweep — outside of that one case, prefer
    :func:`uniform_angle_sampler` or a user-supplied
    :class:`Distribution`.

    The ``rng`` and ``n`` arguments are accepted for ``Sampler``-protocol
    compatibility but intentionally ignored.
    """
    fixed = np.array(
        [math.pi / 4, math.pi / 2, 3 * math.pi / 4, math.pi],
        dtype=float,
    )

    def _sample(rng: np.random.Generator, n: int) -> np.ndarray:
        # Always returns the four hardcoded angles; ``n`` is ignored on
        # the legacy path because the count is fixed by the Fortran
        # ``calculate_stresses`` loop.
        del rng, n
        return fixed.copy()

    return _sample


__all__ = [
    "Distribution",
    "Sampler",
    "default_amplitude_distribution",
    "default_amplitude_sampler",
    "legacy_angle_sampler",
    "uniform_angle_distribution",
    "uniform_angle_sampler",
]
