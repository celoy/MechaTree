"""Wind models for MechaTree.

The canonical canopy-aware (screening) wind model is the native 3-D
**momentum-wind** CFD bridge in :mod:`mechatree.wind.momentum_wind`
(:class:`MomentumWindBridge` / :func:`make_momentum_wind_fn`), backed by the
pure-NumPy kernel in :mod:`mechatree.wind._momentum_wind_kernel`. Select it
via the YAML ``wind:`` block::

    wind:
      model: momentum
      grid_size: 2.0
      momentum_U_uniform: 1.6

Storm statistics (:mod:`mechatree.wind.distributions`) and the storm-replay
diagnostic (:mod:`mechatree.wind.replay`) are model-agnostic helpers.

Step 26 removed the legacy ``native`` bulk-thinning and ``dendroflow``
bridges: both collapsed the 3-D field to a single scalar (equivalent to a
constant wind), which the per-branch momentum model supersedes.
"""

from __future__ import annotations
