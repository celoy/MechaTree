/*
 * Pruning implementation. See pruning.h for the public contract.
 */

#include "pruning.h"

#include <array>
#include <cmath>
#include <random>
#include <unordered_set>
#include <vector>

#include "branch.h"
#include "mechanics.h"  // wind_force

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double NORM_EPS = 1e-12;
constexpr double D_REF = 0.1;       // diameter scale for vol_relative
constexpr double M_EXP = 10.0;      // stress exponent in P_fail

inline std::array<double, 3> add(const std::array<double, 3>& a,
                                  const std::array<double, 3>& b) {
    return {{a[0] + b[0], a[1] + b[1], a[2] + b[2]}};
}

inline std::array<double, 3> scale(double s, const std::array<double, 3>& a) {
    return {{s * a[0], s * a[1], s * a[2]}};
}

inline std::array<double, 3> cross(const std::array<double, 3>& a,
                                    const std::array<double, 3>& b) {
    return {{a[1] * b[2] - a[2] * b[1],
             a[2] * b[0] - a[0] * b[2],
             a[0] * b[1] - a[1] * b[0]}};
}

inline double norm(const std::array<double, 3>& a) {
    return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

// Shared removal tail for prune / prune_with_stored_forces. Captures the
// parent pointers of every cut branch (they survive index-shifting removal),
// removes the subtrees, records the chain-start candidates, and returns the
// number of branches actually removed (including descendants).
int apply_cuts(Tree& tree, const std::vector<int>& to_cut) {
    if (to_cut.empty()) {
        tree.setLastPruneParents(std::vector<Branch*>{});
        return 0;
    }

    // Capture parent pointers BEFORE removeBranches: indices will shift, but
    // parent pointers (the chain-start candidates) survive the removal.
    // Dedupe in pointer space so a single parent that loses several children
    // is recorded once.
    std::unordered_set<Branch*> parent_set;
    for (int idx : to_cut) {
        if (idx < 0 || idx >= tree.getNumberOfBranches()) continue;
        Branch* b = tree.getBranch(idx);
        Branch* p = b->getParent();
        if (p != nullptr) {
            parent_set.insert(p);
        }
    }

    const int before = tree.getNumberOfBranches();
    tree.removeBranches(to_cut);

    tree.setLastPruneParents(
        std::vector<Branch*>(parent_set.begin(), parent_set.end()));

    return before - tree.getNumberOfBranches();
}

}  // namespace


int prune(Tree& tree,
          const std::array<double, 3>& wind,
          double leaf_drag_S0,
          double cauchy) {
    // Always start with an empty parents list — every prune call refreshes
    // it, including the early-out paths.
    tree.setLastPruneParents(std::vector<Branch*>{});

    const int N = tree.getNumberOfBranches();
    if (N <= 1) return 0;  // never cut the trunk in isolation

    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    auto& rng = tree.rng();
    std::vector<double> rndm(static_cast<std::size_t>(N));
    for (int i = 0; i < N; ++i) {
        rndm[static_cast<std::size_t>(i)] = uniform(rng);
    }

    const double U_norm = norm(wind);
    // Leaf-cluster drag ½ S0 |U| U (½ρU² convention; see mechanics.cpp
    // wind_force — the ½ is uniform across woody + leaf so doubling Cauchy
    // reproduces the no-½ Fortran reference exactly).
    const std::array<double, 3> S0U_norm_U = scale(0.5 * leaf_drag_S0 * U_norm, wind);

    std::vector<int> to_cut;

    std::array<double, 3> force, moment;
    for (int i = N - 1; i >= 0; --i) {
        Branch* b = tree.getBranch(i);
        const bool is_leaf = !b->hasChildren();

        wind_force(*b, wind, force, moment);

        std::array<double, 3> T_total, M_total;
        if (is_leaf) {
            T_total = add(force, S0U_norm_U);
            const auto L_unit_t = scale(b->getLength(), b->getUnitT());
            const auto torqueU = cross(L_unit_t, S0U_norm_U);
            M_total = add(moment, torqueU);
        } else {
            std::array<double, 3> child_T{{0.0, 0.0, 0.0}};
            std::array<double, 3> child_M{{0.0, 0.0, 0.0}};
            for (const Branch* c : b->getChildren()) {
                child_T = add(child_T, c->getForce());
                child_M = add(child_M, c->getMoment());
            }
            T_total = add(force, child_T);
            const auto L_unit_t = scale(b->getLength(), b->getUnitT());
            const auto torqueU = cross(L_unit_t, child_T);
            M_total = add(add(moment, child_M), torqueU);
        }

        b->setForce(T_total[0], T_total[1], T_total[2]);
        b->setMoment(M_total[0], M_total[1], M_total[2]);

        const double d = b->getDiameter();
        if (d <= NORM_EPS) continue;  // diameter is degenerate; skip safely

        const auto bend_moment = cross(b->getUnitT(), M_total);
        const double stress = 16.0 / PI * cauchy * norm(bend_moment) / (d * d * d);
        b->setStress(stress);

        if (i == 0) continue;  // never cut the trunk

        const double vol_relative = (d / D_REF) * (d / D_REF);
        // stress^M_EXP may overflow to +inf; -vol_relative * inf = -inf;
        // exp(-inf) = 0; p_fail = 1. That's the right answer (definite cut).
        const double p_fail = 1.0 - std::exp(-vol_relative * std::pow(stress, M_EXP));
        if (p_fail > rndm[static_cast<std::size_t>(i)]) {
            to_cut.push_back(i);
        }
    }

    return apply_cuts(tree, to_cut);
}


int prune_with_stored_forces(Tree& tree,
                             double leaf_drag_S0,
                             double cauchy) {
    // Always start with an empty parents list — every prune call refreshes
    // it, including the early-out paths.
    tree.setLastPruneParents(std::vector<Branch*>{});

    const int N = tree.getNumberOfBranches();
    if (N <= 1) return 0;  // never cut the trunk in isolation

    // Same RNG draw structure as prune() so an equivalent force/wind setup
    // reproduces prune()'s cuts bit-for-bit.
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    auto& rng = tree.rng();
    std::vector<double> rndm(static_cast<std::size_t>(N));
    for (int i = 0; i < N; ++i) {
        rndm[static_cast<std::size_t>(i)] = uniform(rng);
    }

    std::vector<int> to_cut;

    for (int i = N - 1; i >= 0; --i) {
        Branch* b = tree.getBranch(i);
        const bool is_leaf = !b->hasChildren();

        // Per-branch woody-segment drag from the CFD (option B), replacing
        // the wind_force recompute. The segment moment about the base is the
        // same lever as wind_force uses: cross(0.5 L t, F_seg).
        const auto& F_seg = b->getSegmentForce();
        const auto L_unit_t = scale(b->getLength(), b->getUnitT());
        const auto moment_seg = cross(scale(0.5, L_unit_t), F_seg);

        std::array<double, 3> T_total, M_total;
        if (is_leaf) {
            // Leaf-cluster drag uses the branch's own local CFD wind.
            const auto& w = b->getSegmentWind();
            // Leaf-cluster drag ½ S0 |w| w (½ρU² convention; see wind_force).
            const std::array<double, 3> leaf_drag =
                scale(0.5 * leaf_drag_S0 * norm(w), w);
            T_total = add(F_seg, leaf_drag);
            const auto torqueU = cross(L_unit_t, leaf_drag);
            M_total = add(moment_seg, torqueU);
        } else {
            std::array<double, 3> child_T{{0.0, 0.0, 0.0}};
            std::array<double, 3> child_M{{0.0, 0.0, 0.0}};
            for (const Branch* c : b->getChildren()) {
                child_T = add(child_T, c->getForce());
                child_M = add(child_M, c->getMoment());
            }
            T_total = add(F_seg, child_T);
            const auto torqueU = cross(L_unit_t, child_T);
            M_total = add(add(moment_seg, child_M), torqueU);
        }

        // Reuse force_/moment_ as the leaves-to-trunk aggregation accumulator
        // (read back by this branch's parent). segment_* stays untouched.
        b->setForce(T_total[0], T_total[1], T_total[2]);
        b->setMoment(M_total[0], M_total[1], M_total[2]);

        const double d = b->getDiameter();
        if (d <= NORM_EPS) continue;  // diameter is degenerate; skip safely

        const auto bend_moment = cross(b->getUnitT(), M_total);
        const double stress = 16.0 / PI * cauchy * norm(bend_moment) / (d * d * d);
        b->setStress(stress);

        if (i == 0) continue;  // never cut the trunk

        const double vol_relative = (d / D_REF) * (d / D_REF);
        const double p_fail = 1.0 - std::exp(-vol_relative * std::pow(stress, M_EXP));
        if (p_fail > rndm[static_cast<std::size_t>(i)]) {
            to_cut.push_back(i);
        }
    }

    return apply_cuts(tree, to_cut);
}
