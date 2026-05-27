/*
 * Mechanics implementation. See mechanics.h for the public contract.
 *
 * Geometry note: in the Fortran reference, the variable `costheta` is
 * actually sin(angle between branch and wind) — the projection factor for
 * cylinder drag perpendicular to flow. We follow the same naming so the
 * cross-reference to mod_tree.f90 stays mechanical.
 */

#include "mechanics.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace {

constexpr double PI = 3.14159265358979323846;

inline std::array<double, 3> add(const std::array<double, 3>& a,
                                  const std::array<double, 3>& b) {
    return {{a[0] + b[0], a[1] + b[1], a[2] + b[2]}};
}

inline std::array<double, 3> sub(const std::array<double, 3>& a,
                                  const std::array<double, 3>& b) {
    return {{a[0] - b[0], a[1] - b[1], a[2] - b[2]}};
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

constexpr double NORM_EPS = 1e-12;

}  // namespace


void wind_force(const Branch& b,
                const std::array<double, 3>& V,
                std::array<double, 3>& force,
                std::array<double, 3>& moment) {
    const double L = b.getLength();
    const double d = b.getDiameter();
    const auto& t = b.getUnitT();

    const double V_norm = norm(V);
    if (V_norm < NORM_EPS) {
        force = {{0.0, 0.0, 0.0}};
        moment = {{0.0, 0.0, 0.0}};
        return;
    }
    const std::array<double, 3> u = scale(1.0 / V_norm, V);
    const std::array<double, 3> Nn = cross(t, u);
    const double sin_theta = norm(Nn);
    if (sin_theta < NORM_EPS) {
        // Wind parallel to the branch — no projected area, no drag.
        force = {{0.0, 0.0, 0.0}};
        moment = {{0.0, 0.0, 0.0}};
        return;
    }
    const std::array<double, 3> n = scale(1.0 / sin_theta, Nn);
    const std::array<double, 3> bvec = cross(n, t);

    // |F| = ½ ρ C_D D L U² sin²θ — the standard ½ρU² dynamic-pressure drag
    // (ρ = C_D = 1 here). The Fortran reference (mod_tree.f90:637) and the
    // original port DROPPED the ½, folding it into the Cauchy calibration;
    // restored here so this matches the momentum-wind CFD kernel's convention
    // (which already carries ½ρU²) and is physically correct. The default
    // ``TreeConfig.cauchy`` was doubled to absorb this (see config.py), so the
    // default model still reproduces the Nat Commun reference bit-for-bit.
    const double scale_f = 0.5 * V_norm * V_norm * d * L * sin_theta * sin_theta;
    force = scale(scale_f, bvec);
    moment = cross(scale(0.5 * L, t), force);
}


void calculate_stresses(Tree& tree, double leaf_drag_S0, double cauchy) {
    const int N = tree.getNumberOfBranches();

    // Reset max_stress on every branch.
    for (int i = 0; i < N; ++i) {
        tree.getBranch(i)->setMaxStress(0.0);
    }

    std::array<double, 3> force, moment, torqueU, bend_moment;

    for (int angle = 1; angle <= 4; ++angle) {
        const double theta = PI * static_cast<double>(angle) / 4.0;
        const std::array<double, 3> U{{std::cos(theta), std::sin(theta), 0.0}};
        const double U_norm = norm(U);
        // Leaf-cluster drag ½ S0 |U| U. The ½ matches the ½ρU² dynamic
        // pressure restored in wind_force; on the leaf term it's a convention
        // (the tuned S0 absorbs any constant), applied uniformly so doubling
        // Cauchy reproduces the no-½ Fortran reference exactly.
        const std::array<double, 3> S0U_norm_U = scale(0.5 * leaf_drag_S0 * U_norm, U);

        // Walk branches leaves-to-trunk. In our depth-first ordering, every
        // descendant has a higher index than its parent, so a reverse walk
        // sees all children of branch i before branch i itself.
        for (int i = N - 1; i >= 0; --i) {
            Branch* b = tree.getBranch(i);
            const bool is_leaf = !b->hasChildren();

            wind_force(*b, U, force, moment);

            std::array<double, 3> T_total;
            std::array<double, 3> M_total;
            if (is_leaf) {
                // Leaf adds its own drag (S0 * U * |U|) and its torque on the
                // branch tip.
                T_total = add(force, S0U_norm_U);
                const std::array<double, 3> L_unit_t = scale(b->getLength(), b->getUnitT());
                torqueU = cross(L_unit_t, S0U_norm_U);
                M_total = add(moment, torqueU);
            } else {
                // Internal branch sums forces & moments from immediate
                // children; their force/moment fields were populated this
                // iteration of the angle loop because we walk in reverse.
                std::array<double, 3> child_T{{0.0, 0.0, 0.0}};
                std::array<double, 3> child_M{{0.0, 0.0, 0.0}};
                for (const Branch* c : b->getChildren()) {
                    child_T = add(child_T, c->getForce());
                    child_M = add(child_M, c->getMoment());
                }
                T_total = add(force, child_T);
                const std::array<double, 3> L_unit_t = scale(b->getLength(), b->getUnitT());
                torqueU = cross(L_unit_t, child_T);
                M_total = add(add(moment, child_M), torqueU);
            }

            b->setForce(T_total[0], T_total[1], T_total[2]);
            b->setMoment(M_total[0], M_total[1], M_total[2]);

            bend_moment = cross(b->getUnitT(), M_total);
            const double d = b->getDiameter();
            double stress = 0.0;
            if (d > NORM_EPS) {
                stress = 16.0 / PI * cauchy * norm(bend_moment) / (d * d * d);
            }
            b->setStress(stress);
            b->setMaxStress(std::max(stress, b->getMaxStress()));
        }
    }
}


void calculate_stresses_from_stored_forces(Tree& tree, double leaf_drag_S0,
                                           double cauchy, bool reset_max) {
    // Step 26c: one sensing angle's stress pass from the per-branch screened
    // forces the momentum-wind bridge stored (`segment_force_`) + the local
    // wind it felt (`segment_wind_`, for the leaf-cluster drag). Mirrors the
    // leaves-to-trunk aggregation in `calculate_stresses` but reads the
    // stored forces instead of recomputing `wind_force` from a uniform wind,
    // so a sheltered branch is sensed against the weaker wind it actually
    // feels. The caller drives one solve + one call per sensing angle;
    // `reset_max` is True on the first angle (seed max_stress) and False
    // afterwards (accumulate the per-branch max over angles).
    const int N = tree.getNumberOfBranches();
    for (int i = N - 1; i >= 0; --i) {
        Branch* b = tree.getBranch(i);
        const bool is_leaf = !b->hasChildren();

        const auto& F_seg = b->getSegmentForce();
        const auto L_unit_t = scale(b->getLength(), b->getUnitT());
        const auto moment_seg = cross(scale(0.5, L_unit_t), F_seg);

        std::array<double, 3> T_total, M_total;
        if (is_leaf) {
            const auto& w = b->getSegmentWind();
            // Leaf-cluster drag ½ S0 |w| w (½ρU² convention; see wind_force).
            const std::array<double, 3> leaf_drag = scale(0.5 * leaf_drag_S0 * norm(w), w);
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

        b->setForce(T_total[0], T_total[1], T_total[2]);
        b->setMoment(M_total[0], M_total[1], M_total[2]);

        const auto bend_moment = cross(b->getUnitT(), M_total);
        const double d = b->getDiameter();
        double stress = 0.0;
        if (d > NORM_EPS) {
            stress = 16.0 / PI * cauchy * norm(bend_moment) / (d * d * d);
        }
        b->setStress(stress);
        if (reset_max) {
            b->setMaxStress(stress);
        } else {
            b->setMaxStress(std::max(stress, b->getMaxStress()));
        }
    }
}
