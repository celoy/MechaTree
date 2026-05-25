"""SymPy-backed genome models (Step 15).

Lets users specify the safety factor and allocation splits as closed-form
expressions in YAML, e.g.::

    genome:
      safety: "3 * tanh(max_stress)"
      p_seeds: 0.1
      p_leaves: "0.5 * tanh(nb_leaves / 50)"
      phototropism: 0.5

This module's two factory functions parse the expression(s) with
``sympy.sympify``, validate the free-symbol set against an allow-list,
compile to a fast NumPy callable via ``sympy.lambdify``, and wrap the
result in :class:`mechatree.genome.CallbackSafety` or
:class:`mechatree.genome.CallbackAllocation`. The C++ growth loop then
calls back into the compiled lambda once per branch.

Inputs (per the Fortran reference, ``mod_tree.f90:735`` and ``:771``):

* Safety: ``nb_leaves`` (int → float at evaluation time) and ``max_stress``
  (float). Output should be a non-negative float.
* Allocation: three independent expressions for ``p_seeds``, ``p_leaves``,
  ``phototropism``. Each one takes ``nb_leaves`` and ``vol_relative``.

This module sits between the constant genome (``ConstantSafety`` /
``ConstantAllocation``) and the evolved neural genome
(``NeuralSafety`` / ``NeuralAllocation``). It is purely a research lever —
no evolution, no calibration.

SymPy is an optional dependency. Install with::

    pip install 'mechatree[sympy]'
"""

from __future__ import annotations

from typing import TYPE_CHECKING

try:
    import sympy as _sp
except ImportError as _err:  # pragma: no cover
    raise ImportError(
        "mechatree.sympy_genome requires the SymPy package. "
        "Install with: pip install 'mechatree[sympy]'"
    ) from _err

from mechatree.genome import CallbackAllocation, CallbackSafety

if TYPE_CHECKING:
    pass


SAFETY_SYMBOLS = ("nb_leaves", "max_stress")
ALLOCATION_SYMBOLS = ("nb_leaves", "vol_relative")


def _compile(expr_str: str, allowed: tuple[str, ...]):
    """Parse ``expr_str``, validate its free symbols, return a Python callable.

    Raises ``ValueError`` if the expression references symbols outside the
    allow-list. The returned callable takes positional floats in the same
    order as ``allowed``.
    """
    if not isinstance(expr_str, str):
        raise TypeError(f"_compile expects a string expression; got {type(expr_str).__name__}")
    locals_map = {name: _sp.Symbol(name) for name in allowed}
    try:
        expr = _sp.sympify(expr_str, locals=locals_map)
    except (_sp.SympifyError, SyntaxError) as exc:
        raise ValueError(f"Could not parse SymPy expression {expr_str!r}: {exc}") from exc

    free = {s.name for s in expr.free_symbols}
    unknown = free - set(allowed)
    if unknown:
        raise ValueError(
            f"Expression {expr_str!r} references unknown symbol(s) {sorted(unknown)!r}; "
            f"allowed: {list(allowed)!r}"
        )

    symbols = [locals_map[name] for name in allowed]
    return _sp.lambdify(symbols, expr, modules="numpy")


def sympy_safety(expr: str) -> CallbackSafety:
    """Build a :class:`CallbackSafety` from a SymPy expression.

    The expression may reference ``nb_leaves`` (int per branch) and
    ``max_stress`` (float, the worst peak stress over four wind angles).
    Other symbols raise ``ValueError`` at construction time. SymPy
    functions like ``tanh``, ``exp``, ``Max``, ``Min``, ``Piecewise`` are
    available out of the box::

        safety = sympy_safety("3 * tanh(0.5 * max_stress)")
        tree = grow_tree(cfg, n_generations=100, seed=42, safety=safety)
    """
    fn = _compile(expr, SAFETY_SYMBOLS)

    def _call(nb_leaves: int, max_stress: float) -> float:
        return float(fn(nb_leaves, max_stress))

    # Attach the source expression to the closure so callers can introspect
    # which expression a given bridge wraps (the C++ cdef class itself has
    # fixed slots, so we hang it off the Python callable instead).
    _call.expr = expr  # type: ignore[attr-defined]
    return CallbackSafety(_call)


def sympy_allocation(
    p_seeds: str | float,
    p_leaves: str | float,
    phototropism: str | float,
) -> CallbackAllocation:
    """Build a :class:`CallbackAllocation` from three SymPy expressions.

    Each argument is either a scalar (treated as a constant) or a string
    expression in ``nb_leaves`` / ``vol_relative``. Outputs are *not*
    renormalised here — if the caller wants ``p_seeds + p_leaves <= 1``,
    they must clamp the expressions themselves. (The Fortran neural net
    does its own renormalisation; the constant model trusts the YAML
    author to write valid numbers; we follow the latter convention.)

    ::

        alloc = sympy_allocation(
            p_seeds="0.1 * tanh(vol_relative)",
            p_leaves=0.5,
            phototropism="0.5 + 0.1 * tanh(nb_leaves / 100)",
        )
    """
    fns = []
    for label, value in (
        ("p_seeds", p_seeds),
        ("p_leaves", p_leaves),
        ("phototropism", phototropism),
    ):
        if isinstance(value, (int, float)):
            const_value = float(value)
            fns.append(lambda nb, vr, _v=const_value: _v)
        elif isinstance(value, str):
            fns.append(_compile(value, ALLOCATION_SYMBOLS))
        else:
            raise TypeError(
                f"sympy_allocation.{label} must be a number or string expression; "
                f"got {type(value).__name__}"
            )

    f_seeds, f_leaves, f_photo = fns

    def _call(nb_leaves: int, vol_relative: float) -> tuple[float, float, float]:
        return (
            float(f_seeds(nb_leaves, vol_relative)),
            float(f_leaves(nb_leaves, vol_relative)),
            float(f_photo(nb_leaves, vol_relative)),
        )

    _call.exprs = (p_seeds, p_leaves, phototropism)  # type: ignore[attr-defined]
    return CallbackAllocation(_call)


__all__ = ["sympy_allocation", "sympy_safety"]
