/*
 * Pruning implementation. See pruning.h for the public contract.
 */

#include "pruning.h"

#include <array>
#include <cmath>
#include <random>
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

}  // namespace


int prune(Tree& tree,
          const std::array<double, 3>& wind,
          double leaf_drag_S0,
          double cauchy) {
    const int N = tree.getNumberOfBranches();
    if (N <= 1) return 0;  // never cut the trunk in isolation

    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    auto& rng = tree.rng();
    std::vector<double> rndm(static_cast<std::size_t>(N));
    for (int i = 0; i < N; ++i) {
        rndm[static_cast<std::size_t>(i)] = uniform(rng);
    }

    const double U_norm = norm(wind);
    const std::array<double, 3> S0U_norm_U = scale(leaf_drag_S0 * U_norm, wind);

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

    if (to_cut.empty()) return 0;

    const int before = tree.getNumberOfBranches();
    tree.removeBranches(to_cut);
    return before - tree.getNumberOfBranches();
}
