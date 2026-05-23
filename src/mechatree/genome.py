"""Genome models — the per-branch decision functions consulted by growth.

Step 9 ships two constant-returning built-ins. The C++ side defines an
abstract base class with a virtual ``compute`` method; a future neural-net
subclass drops in by inheritance without touching the growth code that
calls it.
"""

from mechatree._core._core import (
    PyAllocationModel,
    PyConstantAllocation,
    PyConstantSafety,
    PySafetyModel,
)

# Public aliases — the wrapper-prefix is an implementation detail.
SafetyModel = PySafetyModel
AllocationModel = PyAllocationModel
ConstantSafety = PyConstantSafety
ConstantAllocation = PyConstantAllocation


__all__ = [
    "SafetyModel",
    "AllocationModel",
    "ConstantSafety",
    "ConstantAllocation",
]
