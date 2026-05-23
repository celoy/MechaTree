#!/usr/bin/env python3
"""Generate a Blender Python script that renders a simulated tree.

Workflow:

1. Run this script to grow a tree and write ``render_my_tree.py``.
2. Open Blender, then run::

       /Applications/Blender.app/Contents/MacOS/Blender --background \\
           --python render_my_tree.py

   This loads the geometry, renders to PNG, and (optionally) saves a
   ``.blend`` you can open interactively to tweak materials / lighting.

The exporter does not require Blender to be installed — that's only
needed to run the generated script.

Run with::

    uv run python examples/render_blender.py
    uv run python examples/render_blender.py --iterations 200 --output my_tree
"""

from __future__ import annotations

import argparse
from pathlib import Path

from mechatree.config import load_config
from mechatree.export import to_blender_script
from mechatree.simulate import grow_tree


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "forest.yaml",
        help="YAML config file.",
    )
    parser.add_argument("--iterations", type=int, default=100, help="Number of generations.")
    parser.add_argument("--seed", type=int, default=42, help="Master RNG seed.")
    parser.add_argument(
        "--output",
        type=str,
        default="render_my_tree",
        help=(
            "Output base name (writes <output>.py for Blender and the script "
            "will save <output>.png + <output>.blend when run)."
        ),
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    tree = grow_tree(cfg, n_generations=args.iterations, seed=args.seed)
    print(
        f"Grew a tree with {tree.get_number_of_branches()} branches "
        f"after {args.iterations} generations."
    )

    out_base = Path(args.output).resolve()
    script_path = out_base.with_suffix(".py")
    png_path = out_base.with_suffix(".png")
    blend_path = out_base.with_suffix(".blend")

    to_blender_script(
        tree,
        script_path,
        render_path=png_path,
        save_blend_path=blend_path,
        image_resolution=(1280, 720),
    )

    print(f"Wrote: {script_path}")
    print()
    print("Run this in Blender to render the tree (Blender 4.x):")
    print()
    print(f"  /Applications/Blender.app/Contents/MacOS/Blender --background --python {script_path}")
    print()
    print(f"It will write: {png_path}")
    print(f"         and: {blend_path}")


if __name__ == "__main__":
    main()
