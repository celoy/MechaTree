/*
 * Light interception — C++ kernel for the per-direction shadow sort.
 *
 * Ports the body of legacy_fortran/mod_tree.f90:219 (light_interception)
 * across the whole sun direction grid. One call covers every (leaf, direction)
 * pair; the Python wrapper in mechatree.light.interception.intercept just
 * passes the leaf-position buffer and sun arrays through.
 *
 * The previous Python implementation allocated ~15 NumPy arrays per
 * direction (rotated coords, cell keys, lexsort, unique, diff, repeat,
 * tau**depth). At island scale (R = 200 L, ~430k leaves, 32 sun directions)
 * that meant ~150 MB of allocations per generation — the dominant cost.
 * This kernel does the same work with one stack of per-direction sort
 * entries, reused across all directions, and writes `tau**depth` directly
 * into the caller's `light_per_direction` output buffer.
 *
 * Equivalence with the Python reference is exact up to the chosen
 * within-cell tie-break (ascending leaf index), which already matched the
 * Fortran's stable sort.
 */

#ifndef MECHATREE_LIGHT_H_
#define MECHATREE_LIGHT_H_

#include <cstddef>

// Run the shadow sort for every (leaf, direction) pair.
//
//   leaf_locations          — (n_leaves, 3) row-major, float64
//   sun_elev, sun_azim      — (n_directions,) sun-direction angles, radians
//   size_leaf               — cell width in world units (Fortran SizeLeaf)
//   leaf_transparency (tau) — per-leaf transmittance in [0, 1]
//                             (i-th leaf from the top gets tau**i)
//   light_per_direction     — (n_leaves, n_directions) row-major output,
//                             OVERWRITTEN in place
//
// Pre-condition: light_per_direction has length n_leaves * n_directions.
// Pre-condition: tau in [0, 1] (caller validates).
void light_intercept(const double* leaf_locations,
                     std::size_t n_leaves,
                     const double* sun_elev,
                     const double* sun_azim,
                     std::size_t n_directions,
                     double size_leaf,
                     double leaf_transparency,
                     double* light_per_direction);

#endif  // MECHATREE_LIGHT_H_
