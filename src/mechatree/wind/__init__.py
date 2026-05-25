"""Wind models for MechaTree (Step 17).

The only entry point today is the DendroFlow bridge — a thin wrapper around
DendroFlow's ``BulkThinningBranchWindModel`` that turns a live ``PyTree`` /
``Forest`` into a streamwise canopy-mean wind vector suitable for MechaTree's
per-generation ``wind_fn`` callback.

The bridge module is import-gated on the optional ``dendroflow`` extra so
``import mechatree`` keeps working without DendroFlow installed. To use it::

    pip install 'mechatree[dendroflow]'
    # Then in YAML:
    #   wind:
    #     model: dendroflow
    #     U_infty: [...]
    #     z_centers: [...]
    #     H: 0.5
    #     C_D: 1.0

or, programmatically::

    from mechatree.wind import make_dendroflow_wind_fn

    wind_fn = make_dendroflow_wind_fn(U_infty=[...], z_centers=[...], H=0.5)
    grow_tree(cfg, n_generations=50, wind_fn=wind_fn)
"""

from __future__ import annotations


def __getattr__(name: str):
    # Lazy re-export: pulling in ``dendroflow`` should only happen when the
    # bridge symbols are actually requested. Lets ``import mechatree.wind``
    # remain side-effect-free for users who don't have the extra installed.
    if name in {
        "BranchWindBridge",
        "DendroFlowWindParams",
        "make_dendroflow_wind_fn",
        "pytree_to_cylinders",
        "forest_to_cylinders",
    }:
        from mechatree.wind import dendroflow as _df

        return getattr(_df, name)
    raise AttributeError(f"module 'mechatree.wind' has no attribute {name!r}")


__all__ = [
    "BranchWindBridge",
    "DendroFlowWindParams",
    "forest_to_cylinders",
    "make_dendroflow_wind_fn",
    "pytree_to_cylinders",
]
