"""Wall time of ``intercept`` over varying leaf counts.

Light interception is the dominant cost of a simulation generation once
trees have non-trivial canopies — so this benchmark exists to catch
regressions in the vectorised path.

The cloud is a random ``[-15, 15] x [-15, 15] x [0, 20]`` box; the sun is
the default ``Sun()`` (4 elevations x 8 azimuths = 32 directions). Best
of 3 wall times after a warmup call.
"""

import time

import numpy as np

from mechatree.light import Leaves, Sun, intercept


def _make_leaves(n: int, rng: np.random.Generator, n_dir: int) -> Leaves:
    locs = np.zeros((n, 3))
    locs[:, 0] = rng.uniform(-15.0, 15.0, n)
    locs[:, 1] = rng.uniform(-15.0, 15.0, n)
    locs[:, 2] = rng.uniform(0.0, 20.0, n)
    return Leaves(
        location=locs,
        branch_index=np.arange(n, dtype=np.int32),
        tree_index=np.zeros(n, dtype=np.int32),
        light_per_direction=np.zeros((n, n_dir), dtype=np.float64),
    )


def main() -> None:
    sun = Sun()
    rng = np.random.default_rng(0)
    sizes = [100, 500, 1000, 2500, 5000, 10000]

    print(f"{'N':>6}  {'best (ms)':>10}  {'ms/dir':>8}  {'us/leaf':>8}")
    for n in sizes:
        leaves = _make_leaves(n, rng, sun.n_directions)
        intercept(leaves, sun)  # warmup
        best = float("inf")
        for _ in range(3):
            t0 = time.perf_counter()
            intercept(leaves, sun)
            best = min(best, (time.perf_counter() - t0) * 1e3)
        ms_per_dir = best / sun.n_directions
        us_per_leaf = best * 1e3 / n
        print(f"{n:>6}  {best:>10.2f}  {ms_per_dir:>8.3f}  {us_per_leaf:>8.2f}")


if __name__ == "__main__":
    main()
