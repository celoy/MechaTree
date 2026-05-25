#!/usr/bin/env python3
"""Random 3D growth demo.

Builds a tree by repeatedly attaching random children to random branches,
occasionally pruning, and renders the final result with :func:`plot_3d`.

Port of the 2017 intern's ``random_growth.py``.
"""

from __future__ import annotations

import argparse
import random

import numpy as np

import mechatree as mt


def run(n_generations: int, p_remove: float = 0.05, max_children: int = 4) -> mt.PyTree:
    trunk = {"x": 0, "y": 0, "z": 0, "theta": 0, "phi": 0}
    tree = mt.PyTree(trunk)

    for _ in range(n_generations):
        n_branches = tree.get_number_of_branches()
        mother = int(random.random() * (n_branches - 1))
        nbchild = int(random.random() * (max_children - tree.get_number_of_children(mother)))

        xmom = tree.get_property(mother, "x")
        ymom = tree.get_property(mother, "y")
        zmom = tree.get_property(mother, "z")
        thetamom = tree.get_property(mother, "theta")
        phimom = tree.get_property(mother, "phi")

        x = xmom + np.sin(thetamom) * np.cos(phimom)
        y = ymom + np.sin(thetamom) * np.sin(phimom)
        z = zmom + np.cos(thetamom)

        if z >= 1:
            for _k in range(nbchild):
                theta = np.pi * 2.0 * random.random()
                phi = np.pi * 0.5 * random.random()
                tree.add_branch(mother, {"x": x, "y": y, "z": z, "theta": theta, "phi": phi})

        if random.random() > (1.0 - p_remove) and tree.get_number_of_branches() > 1:
            n = tree.get_number_of_branches()
            branch2rmv = int(random.random() * (n - 2)) + 1
            tree.remove_branch(branch2rmv)

    return tree


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=300, help="Number of generations.")
    parser.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility.")
    parser.add_argument("--no-show", action="store_true", help="Skip opening the plotly figure.")
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    tree = run(args.iterations)
    print(f"Final tree size: {tree.get_number_of_branches()}")

    fig = mt.plot_3d(tree)
    if not args.no_show:
        fig.show()


if __name__ == "__main__":
    main()
