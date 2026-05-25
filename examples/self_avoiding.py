#!/usr/bin/env python3
"""Self-avoiding 2D growth demo.

Grows a coral-like ramified structure where branches stop growing when they
get too close to other branches. A snapshot is saved at the final iteration
(use ``--snapshot-every`` to capture intermediate frames for a video).

Port of the 2017 intern's ``self_avoiding_modules.py``.
"""

from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

import mechatree as mt
from mechatree.geometry import distance_test  # not part of the flat surface yet


def run(
    n_steps: int,
    out_dir: Path,
    snapshot_every: int | None,
    dL: float = 0.1,
    p_b: float = 0.2,
    mean_angle: float = math.pi / 3,
    std_angle: float = math.pi / 10,
) -> mt.PyTree:
    influence_radius = 5 * dL
    L_b = 2 * influence_radius

    trunk = {"x": 0, "y": 0, "theta": 1e-6, "L": 2.0, "grow": 1}
    coral = mt.PyTree(trunk)

    for i in range(n_steps):
        if snapshot_every is not None and i % snapshot_every == 0:
            mt.plot_2d(coral, iteration=i, out_dir=out_dir)

        # Growing loop.
        for n in range(coral.get_number_of_branches()):
            if coral.get_property(n, "grow") != 1:
                continue

            y_0 = coral.get_property(n, "y")
            L = coral.get_property(n, "L")
            theta = coral.get_property(n, "theta")
            if y_0 + L * math.cos(theta) < 0.5:
                coral.set_property(n, "grow", 0)
                continue

            if distance_test(coral, n, influence_radius) == 1:
                coral.set_property(n, "L", L + dL)
            else:
                coral.set_property(n, "grow", 0)

        # Branching loop.
        nb = coral.get_number_of_branches()
        n = 0
        while n < nb:
            if coral.get_property(n, "grow") == 1:
                L = coral.get_property(n, "L")
                if random.random() <= p_b and L_b < L:
                    angle = random.gauss(mean_angle, std_angle)
                    theta = coral.get_property(n, "theta")
                    thetaplus = theta + angle / 2
                    thetaminus = theta - angle / 2

                    x = coral.get_property(n, "x")
                    y = coral.get_property(n, "y")
                    x_0 = x + L * math.sin(theta)
                    y_0 = y + L * math.cos(theta)

                    coral.add_branch(
                        n, {"x": x_0, "y": y_0, "theta": thetaplus, "L": 0.0, "grow": 1}
                    )
                    coral.add_branch(
                        n, {"x": x_0, "y": y_0, "theta": thetaminus, "L": 0.0, "grow": 1}
                    )
                    coral.set_property(n, "grow", 0)
                    n += 2
                    nb = coral.get_number_of_branches()
            n += 1

    # Always emit a final snapshot.
    mt.plot_2d(coral, iteration=n_steps, out_dir=out_dir)
    return coral


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=200, help="Number of timesteps.")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("SimResults/selfAvoidImages"),
        help="Directory for snapshot PNGs.",
    )
    parser.add_argument(
        "--snapshot-every",
        type=int,
        default=None,
        help="Save a snapshot every N steps (default: only the final frame).",
    )
    parser.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility.")
    parser.add_argument("--no-show", action="store_true", help="Unused; accepted for symmetry.")
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    coral = run(args.iterations, args.out_dir, args.snapshot_every)
    print(f"Final coral size: {coral.get_number_of_branches()} branches; PNGs in {args.out_dir}")


if __name__ == "__main__":
    main()
