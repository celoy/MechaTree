/*
 * Pruning — wind-driven stochastic branch removal.
 *
 * Ported from legacy_fortran/mod_tree.f90:850 (subroutine pruning). Unlike
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

#endif
