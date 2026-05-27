/*
 * Mechanics — wind force on a branch and stress propagation up the tree.
 *
 * Ported from legacy/fortran/mod_tree.f90:613 (wind_force) and :642
 * (calculate_stresses). The 4-angle sweep that calculate_stresses runs
 * internally matches the Fortran reference; max_stress on each branch is the
 * worst of the four.
 */

#ifndef MECHANICS_H_
#define MECHANICS_H_

#include <array>

#include "branch.h"
#include "tree.h"

// Force and moment exerted by wind of velocity V on a single branch.
//
// `force` is the drag force; `moment` is the moment about the branch's base
// (its location, i.e. the parent's tip). When the wind is aligned with the
// branch's unit_t, force and moment are both zero (the Fortran formula would
// divide by zero; we guard against that here).
void wind_force(const Branch& b,
                const std::array<double, 3>& V,
                std::array<double, 3>& force,
                std::array<double, 3>& moment);

// Sweep four horizontal wind angles (pi/4, pi/2, 3pi/4, pi). For each:
// compute per-branch stress by propagating leaf and child forces up to the
// trunk. After all four angles, every branch's `max_stress` holds the worst
// stress it experienced.
//
// `leaf_drag_S0` is the leaf-surface drag coefficient (Fortran S0). `cauchy`
// is the material stiffness constant (Fortran Cy). The caller must have
// called `tree.reorder()` since the last structural change so that the
// children lists are walked in the right order.
void calculate_stresses(Tree& tree, double leaf_drag_S0, double cauchy);

// Step 26c: one sensing angle's stress pass from the per-branch forces the
// momentum-wind bridge pre-stored on each branch (`segment_force_` woody
// drag + `segment_wind_` local wind for the leaf-cluster term), instead of
// recomputing `wind_force` from a single uniform wind. The caller runs one
// momentum solve + one call per sensing angle; pass `reset_max=true` on the
// first angle (seed `max_stress`) and `false` afterwards (accumulate the
// per-branch max over angles). Lets sensing see the same canopy screening
// that pruning does, so a sheltered branch reinforces against the weaker
// wind it actually feels.
void calculate_stresses_from_stored_forces(Tree& tree, double leaf_drag_S0,
                                           double cauchy, bool reset_max);

#endif
