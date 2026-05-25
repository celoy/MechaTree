"""Grow one tree using DendroFlow's lean wind model (Step 17 / DendroFlow M6).

Reads ``examples/dendroflow_wind.yaml`` and drives :func:`grow_tree` with the
DendroFlow bridge wired in via the YAML's ``wind:`` block. Per generation,
prints ``(gen, n_branches, n_leaves, wind_amplitude)`` and confirms the
streamwise canopy mean never exceeds the free-stream maximum (a sanity floor
for the bulk-thinning model).

Usage::

    uv pip install -e .[dev,dendroflow]
    uv pip install -e ../DendroFlow   # until DendroFlow ships to PyPI
    uv run python examples/dendroflow_wind.py
"""

from __future__ import annotations

from pathlib import Path

from mechatree.config import Config
from mechatree.simulate import TreeStats, grow_tree

THIS_DIR = Path(__file__).resolve().parent
YAML_PATH = THIS_DIR / "dendroflow_wind.yaml"


def main() -> None:
    cfg = Config.from_yaml(YAML_PATH)
    free_stream_max = max(cfg.wind.U_infty or ())

    print(f"DendroFlow wind: U_infty in [{min(cfg.wind.U_infty):.1f}, {free_stream_max:.1f}]")
    print(f"  {len(cfg.wind.U_infty)} z-layers, H={cfg.wind.H}, C_D={cfg.wind.C_D}")
    print()
    print(f"{'gen':>4} {'n_br':>6} {'n_lvs':>6} {'wind':>8}")

    def on_step(gen: int, tree, stats: TreeStats) -> None:
        print(
            f"{stats.generation:>4d} "
            f"{stats.n_branches:>6d} "
            f"{stats.n_leaves:>6d} "
            f"{stats.wind_amplitude:>8.3f}"
        )
        # Sanity floor: bulk-thinning never produces canopy wind > free stream.
        assert stats.wind_amplitude <= free_stream_max + 1e-9, (
            f"canopy mean {stats.wind_amplitude} exceeded free stream {free_stream_max}"
        )

    tree = grow_tree(cfg, seed=42, on_step=on_step)
    print()
    print(f"final: {tree.get_number_of_branches()} branches, height ≈ {_height(tree):.2f}")


def _height(tree) -> float:
    n = tree.get_number_of_branches()
    if n == 0:
        return 0.0
    return max(tree.get_location(i)[2] for i in range(n))


if __name__ == "__main__":
    main()
