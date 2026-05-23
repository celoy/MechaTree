#!/usr/bin/env python3
"""Sap-transport resource-allocation demo.

Leaves produce a fixed resource; every branch consumes proportional to its
volume; branches whose reserve goes negative are pruned; surplus reserves
feed creation of new child branches. Volume of ancestors is incremented
recursively whenever a new branch is grown.

Port of the 2017 intern's ``sap_transport.py``.
"""

from __future__ import annotations

import argparse
import math
import random
import sys

from mechatree import PyTree


def increment_volume(tree: PyTree, index: int) -> None:
    if index < tree.get_number_of_branches():
        tree.set_property(index, "volume", tree.get_property(index, "volume") + 1)
        if tree.has_parent(index) == 1:
            increment_volume(tree, tree.get_parent_index(index))


def run(sim_time: int, alpha: float = 2.0, cost: float = 1.0) -> PyTree | None:
    trunk = {
        "x": 1e-6,
        "y": 1e-6,
        "theta": 1e-6,
        "volume": 1,
        "reserve": 0,
        "production": 10,
        "leaf": 1,
    }
    tree = PyTree(trunk)

    for _t in range(sim_time):
        # Leaves -> trunk.
        for i in reversed(range(tree.get_number_of_branches())):
            if tree.get_property(i, "leaf") == 1:
                tree.set_property(i, "reserve", tree.get_property(i, "production"))

            maintenance = alpha * tree.get_property(i, "volume")
            tree.set_property(i, "reserve", tree.get_property(i, "reserve") - maintenance)

            if tree.get_property(i, "reserve") < 0:
                tree.remove_branch(i)
                if tree.get_number_of_branches() == 0:
                    print("Tree died.")
                    return None
            elif tree.has_parent(i) == 1:
                parent_index = tree.get_parent_index(i)
                tree.set_property(parent_index, "reserve", tree.get_property(i, "reserve"))
                tree.set_property(i, "reserve", 0)

        # Trunk -> leaves.
        for i in range(tree.get_number_of_branches()):
            excedent = int(math.floor(tree.get_property(i, "reserve") - cost))
            if excedent <= 0:
                continue

            nc = tree.get_number_of_children(i)
            target = int(round(random.gauss(2.5, 0.5)))
            if nc < target:
                x = tree.get_property(i, "x") + math.sin(tree.get_property(i, "theta"))
                y = tree.get_property(i, "y") + math.cos(tree.get_property(i, "theta"))
                theta = random.uniform(-math.pi * 0.5, 0.5 * math.pi)
                tree.add_branch(
                    i,
                    {
                        "x": x,
                        "y": y,
                        "theta": theta,
                        "volume": 1,
                        "reserve": 0,
                        "production": 1,
                        "leaf": 1,
                    },
                )
                if tree.get_property(i, "leaf") == 1:
                    tree.set_property(i, "leaf", 0)
                increment_volume(tree, i)

    return tree


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=100, help="Simulation time.")
    parser.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility.")
    parser.add_argument("--no-show", action="store_true", help="Unused; accepted for symmetry.")
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    tree = run(args.iterations)
    if tree is None:
        sys.exit(0)
    print(f"Final tree size: {tree.get_number_of_branches()} branches.")


if __name__ == "__main__":
    main()
