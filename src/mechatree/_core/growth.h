/*
 * Growth — three free functions ported from mod_tree.f90:
 *   - requested_growth      (mod_tree.f90:704)
 *   - secondary_growth      (mod_tree.f90:817)
 *   - primary_growth        (mod_tree.f90:425)
 *
 * `requested_growth` and `primary_growth` consult a polymorphic genome model
 * (genome.h) — Step 9 ships only the Constant subclasses; a Neural subclass
 * is a self-contained later step.
 */

#ifndef GROWTH_H_
#define GROWTH_H_

#include "genome.h"
#include "tree.h"

// Compute per-branch `vol_growth`, `vol_summed`, `maintenance_vol` based on
// each branch's max_stress (set by `calculate_stresses`) and the genome's
// safety factor. Caller must have run `tree.reorder()` since the last
// structural change (this routine reads nb_leaves).
void requested_growth(Tree& tree, const SafetyModel& safety, double maintenance_h);

// Allocate photosynthate (light * volume_per_leaf, summed over leaves) along
// the leaf-to-root chains. Grows branch diameters; leftover photosynthate
// feeds the tree's reserve pool.
void secondary_growth(Tree& tree, double volume_per_leaf);

// Spawn new twig branches at the most-lit leaves. Returns the number of new
// branches actually created (always even — twigs are added in pairs).
//
// The allocation model is consulted once for the tree, yielding
// (p_seeds, p_leaves, phototropism). Branching angles theta1/theta2 and
// gamma1/gamma2 are passed in directly — they are physical tree-shape
// parameters, separable from the per-step allocation decision.
int primary_growth(
    Tree& tree,
    const AllocationModel& alloc,
    double twig_length, double twig_diameter,
    double theta1, double theta2,
    double gamma1, double gamma2,
    int generation);

#endif
