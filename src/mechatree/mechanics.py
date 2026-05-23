"""Mechanics — wind force on a branch and stress propagation up the tree.

These functions read and write the typed mechanics fields on each Branch
(force, moment, stress, max_stress). They do not touch the property map.
"""

from mechatree._core._core import PyTree


def wind_force(tree: PyTree, branch_index: int, wind):
    """Return (force, moment) on a single branch under wind velocity ``wind``.

    ``wind`` is a 3-tuple in the same units used to set ``leaf_drag_S0`` and
    ``cauchy`` for ``calculate_stresses``. When the wind is parallel to the
    branch's ``unit_t`` (no projected area), both outputs are zero.
    """
    return tree.wind_force(branch_index, wind)


def calculate_stresses(tree: PyTree, leaf_drag_S0: float, cauchy: float) -> None:
    """Sweep four horizontal wind angles and write per-branch ``max_stress``.

    ``leaf_drag_S0`` is the leaf-surface drag coefficient (Fortran ``S0``);
    ``cauchy`` is the material stiffness constant (Fortran ``Cy``). The
    caller must have called ``tree.reorder()`` since the last structural
    change.
    """
    tree.calculate_stresses(leaf_drag_S0, cauchy)


__all__ = ["wind_force", "calculate_stresses"]
