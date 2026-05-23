"""Light interception — port of ``legacy_fortran/mod_tree.f90`` light routines.

Independent of ``PyTree`` at the algorithm layer: ``intercept`` takes only
a ``Leaves`` (struct of arrays) and a ``Sun`` (direction grid). The only
``PyTree`` coupling is at the I/O edges — ``extract_leaves`` reads tree
geometry, ``aggregate_onto_trees`` writes the per-leaf scalar back.
"""

from mechatree.light.interception import aggregate_onto_trees, intercept
from mechatree.light.leaves import Leaves, extract_leaves
from mechatree.light.sun import Sun

__all__ = [
    "Leaves",
    "Sun",
    "aggregate_onto_trees",
    "extract_leaves",
    "intercept",
]
