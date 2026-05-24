"""Pre-sampled rainbow palette used to colour branches by Strahler order.

Sampled from matplotlib's ``rainbow`` cmap at positions ``k / 9`` for
``k ∈ [0, 9]`` so the runtime has no matplotlib dependency. The Nat Comm
2017 paper convention treats 0 as purple and 9 as red; order 1 lands at
position 1/9 (pale purple) and order 9 at position 1 (red).
"""

from __future__ import annotations

MAX_STRAHLER = 9

# RGB tuples in [0, 1]; index k = Strahler order k.
RAINBOW_STRAHLER: tuple[tuple[float, float, float], ...] = (
    (0.5000, 0.0000, 1.0000),  # 0 – deep purple (unused; Strahler starts at 1)
    (0.2804, 0.3382, 0.9852),  # 1 – pale purple / blue
    (0.0608, 0.6365, 0.9411),  # 2 – blue
    (0.1667, 0.8660, 0.8660),  # 3 – cyan
    (0.3863, 0.9841, 0.7674),  # 4 – green
    (0.6137, 0.9841, 0.6412),  # 5 – pale green
    (0.8333, 0.8660, 0.5000),  # 6 – yellow-green
    (1.0000, 0.6365, 0.3382),  # 7 – orange
    (1.0000, 0.3382, 0.1716),  # 8 – red-orange
    (1.0000, 0.0000, 0.0000),  # 9 – red
)


def strahler_rgb(order: int) -> tuple[float, float, float]:
    """Return the RGB tuple for Strahler order ``order`` (clipped to [0, 9])."""
    return RAINBOW_STRAHLER[max(0, min(MAX_STRAHLER, int(order)))]


def strahler_css(order: int, alpha: float = 1.0) -> str:
    """Return a CSS ``rgb()`` / ``rgba()`` string for the given Strahler order."""
    r, g, b = strahler_rgb(order)
    if alpha >= 1.0:
        return f"rgb({int(255 * r)},{int(255 * g)},{int(255 * b)})"
    return f"rgba({int(255 * r)},{int(255 * g)},{int(255 * b)},{alpha:.3f})"


__all__ = ["MAX_STRAHLER", "RAINBOW_STRAHLER", "strahler_css", "strahler_rgb"]
