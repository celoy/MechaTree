/*
 * Pruning — wind-driven stochastic branch removal.
 *
 * Ported from legacy/fortran/mod_tree.f90:850 (subroutine pruning). Unlike
 * `calculate_stresses` which sweeps four angles, pruning evaluates one
 * specific wind direction U — the "what cuts under this gust" question.
 *
 * Per Fortran reference: P_fail = 1 - exp(-(d/0.1)^2 * stress^10). A branch
 * is cut when a uniform [0,1) draw is below P_fail. The trunk is never cut.
 */

#ifndef PRUNING_H_
#define PRUNING_H_

#include <array>

#include "tree.h"

// Run one pruning pass. Returns the number of branches removed (including
// every descendant of any branch that was directly cut).
int prune(Tree& tree,
          const std::array<double, 3>& wind,
          double leaf_drag_S0,
          double cauchy);

// Step 25c (option B): same pruning pass, but instead of recomputing each
// branch's woody-segment drag from a single canopy-mean wind via
// `wind_force`, it reads the per-branch force the momentum-wind CFD already
// stored on each branch (`Branch::segment_force_`) and the per-branch local
// wind (`Branch::segment_wind_`, used for the leaf-cluster drag term on
// terminals). The bridge must have populated those fields for the current
// branch set before calling this. Otherwise byte-for-byte identical to
// `prune` (same RNG draw structure, Weibull test, removal tail).
int prune_with_stored_forces(Tree& tree,
                             double leaf_drag_S0,
                             double cauchy);

#endif
