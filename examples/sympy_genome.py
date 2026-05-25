"""Grow one tree with a SymPy-callable genome (Step 15).

Reads ``examples/sympy_genome.yaml`` and grows a tree with safety + allocation
defined as closed-form SymPy expressions in YAML. The expressions are
compiled once via ``sympy.lambdify`` and called per-branch from the C++
growth loop through the new ``CallbackSafety`` / ``CallbackAllocation``
vtable subclasses.

Usage::

    uv pip install -e .[dev,sympy]
    uv run python examples/sympy_genome.py
"""

from __future__ import annotations

from pathlib import Path

from mechatree.config import Config
from mechatree.simulate import TreeStats, grow_tree

THIS_DIR = Path(__file__).resolve().parent
YAML_PATH = THIS_DIR / "sympy_genome.yaml"


def main() -> None:
    cfg = Config.from_yaml(YAML_PATH)
    print("SymPy genome:")
    print(f"  safety:       {cfg.genome.safety!r}")
    print(f"  p_seeds:      {cfg.genome.p_seeds!r}")
    print(f"  p_leaves:     {cfg.genome.p_leaves!r}")
    print(f"  phototropism: {cfg.genome.phototropism!r}")
    print()
    print(f"{'gen':>4} {'n_br':>6} {'n_lvs':>6} {'reserve':>10}")

    def on_step(gen: int, tree, stats: TreeStats) -> None:
        if gen % 10 == 0 or gen == cfg.n_generations - 1:
            print(
                f"{stats.generation:>4d} "
                f"{stats.n_branches:>6d} "
                f"{stats.n_leaves:>6d} "
                f"{stats.reserve:>10.4f}"
            )

    tree = grow_tree(cfg, seed=42, on_step=on_step)
    print()
    print(
        f"final: {tree.get_number_of_branches()} branches, "
        f"trunk diameter = {tree.get_diameter(0):.4f}"
    )


if __name__ == "__main__":
    main()
