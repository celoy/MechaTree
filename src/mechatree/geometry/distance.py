"""Self-avoidance distance check for 2D branch growth.

Port of the 2017 intern's ``mod_dist.pyx``. Pure Python; uses stdlib ``math``.
"""

from __future__ import annotations

import math


def distance_test(tree, i: int, influence_radius: float) -> int:
    """Decide whether branch ``i`` can grow without colliding with another branch.

    The branch's own parent and brother(s) are excluded from the check (they
    necessarily share an endpoint with branch ``i``).

    Returns
    -------
    int
        ``1`` if the tip of branch ``i`` is far enough from every other branch
        in ``tree``, otherwise ``0``.

    Notes
    -----
    Preserves the exact perpendicular-distance formula from the published
    Nat. Commun. (2017) simulator, including its use of ``1 / tan(theta)`` —
    which is undefined at ``theta = 0``. The archive avoids this by
    initialising the trunk with ``theta = 1e-6``; callers should do the same.
    """

    x_0 = tree.get_property(i, "x")
    y_0 = tree.get_property(i, "y")
    length_i = tree.get_property(i, "L")
    theta_i = tree.get_property(i, "theta")

    x_tip = x_0 + length_i * math.sin(theta_i)
    y_tip = y_0 + length_i * math.cos(theta_i)

    n_branches = tree.get_number_of_branches()

    excluded: set[int] = {i}
    if tree.has_parent(i) == 1:
        excluded.add(tree.get_parent_index(i))
        brothers = tree.get_brothers_index(i)
        # PyTree.get_brothers_index returns: 0 (no brothers), int (one brother),
        # or list[int] (many). Normalise.
        if isinstance(brothers, list):
            excluded.update(brothers)
        elif brothers != 0:
            excluded.add(brothers)

    for j in range(n_branches):
        if j in excluded:
            continue

        xj = tree.get_property(j, "x")
        yj = tree.get_property(j, "y")
        theta_j = tree.get_property(j, "theta")
        length_j = tree.get_property(j, "L")

        cot = 1.0 / math.tan(theta_j)
        d = math.fabs((x_tip - xj) * cot - (y_tip - yj)) / math.sqrt(cot * cot + 1.0)

        xf = xj + length_j * math.sin(theta_j)
        yf = yj + length_j * math.cos(theta_j)

        hi = math.sqrt((x_tip - xj) ** 2 + (y_tip - yj) ** 2)
        hf = math.sqrt((x_tip - xf) ** 2 + (y_tip - yf) ** 2)

        di = math.sqrt(max(hi * hi - d * d, 0.0))
        df = math.sqrt(max(hf * hf - d * d, 0.0))

        if di < length_j and df < length_j:
            if d <= influence_radius:
                return 0
        else:
            if hi <= influence_radius or hf <= influence_radius:
                return 0

    return 1
