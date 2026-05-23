"""Growth — primary, secondary and requested growth routines."""

from mechatree._core._core import PyTree
from mechatree.genome import AllocationModel, SafetyModel


def requested_growth(tree: PyTree, safety: SafetyModel, maintenance_h: float) -> None:
    """Compute per-branch ``vol_growth``, ``vol_summed``, ``maintenance_vol``.

    Reads ``max_stress`` (set by ``mechanics.calculate_stresses``) and
    ``nb_leaves`` (set by ``tree.reorder()``). Mutates the typed fields.
    """
    tree.requested_growth(safety, maintenance_h)


def secondary_growth(tree: PyTree, volume_per_leaf: float) -> None:
    """Allocate photosynthate (``light * volume_per_leaf``) along leaf-to-root
    chains. Grows diameters; leftover photosynthate feeds ``tree.reserve``.
    """
    tree.secondary_growth(volume_per_leaf)


def primary_growth(
    tree: PyTree,
    alloc: AllocationModel,
    twig_length: float,
    twig_diameter: float,
    theta1: float,
    theta2: float,
    gamma1: float,
    gamma2: float,
    generation: int,
) -> int:
    """Spawn new twig branches at the most-lit leaves.

    Returns the number of branches actually added (always even — twigs come
    in pairs). The allocation model is consulted once for the tree; per-leaf
    light values drive which leaves get selected.
    """
    return tree.primary_growth(
        alloc,
        twig_length,
        twig_diameter,
        theta1,
        theta2,
        gamma1,
        gamma2,
        generation,
    )


__all__ = ["requested_growth", "secondary_growth", "primary_growth"]
