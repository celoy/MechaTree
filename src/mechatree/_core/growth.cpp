/*
 * Growth implementation. See growth.h for the public contract.
 *
 * All three functions are direct ports of the Fortran reference; line
 * numbers in the comments below refer to legacy/fortran/mod_tree.f90.
 */

#include "growth.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <vector>

#include "branch.h"

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double DEG_TO_RAD_10 = 0.174532925;  // ~10° — Fortran uses this as
                                               // the std-dev of angle noise

inline std::array<double, 3> cross(const std::array<double, 3>& a,
                                    const std::array<double, 3>& b) {
    return {{a[1] * b[2] - a[2] * b[1],
             a[2] * b[0] - a[0] * b[2],
             a[0] * b[1] - a[1] * b[0]}};
}

// Daughter unit vectors — ported from mod_tree.f90:530 (daughter_unit_vectors).
//   t = cos(theta) * mother.unit_t + sin(theta) * (mother.unit_t × mother.unit_b)
//   n = cos(gamma) * mother.unit_b + sin(gamma) * (t × mother.unit_b)
void daughter_unit_vectors(
        const std::array<double, 3>& m_t,
        const std::array<double, 3>& m_b,
        double theta, double gamma,
        std::array<double, 3>& t_out,
        std::array<double, 3>& n_out) {
    const auto b1 = cross(m_t, m_b);
    const double ct = std::cos(theta), st = std::sin(theta);
    t_out = {{ct * m_t[0] + st * b1[0],
              ct * m_t[1] + st * b1[1],
              ct * m_t[2] + st * b1[2]}};
    const auto b2 = cross(t_out, m_b);
    const double cg = std::cos(gamma), sg = std::sin(gamma);
    n_out = {{cg * m_b[0] + sg * b2[0],
              cg * m_b[1] + sg * b2[1],
              cg * m_b[2] + sg * b2[2]}};
}

}  // namespace


void requested_growth(Tree& tree, const SafetyModel& safety, double maintenance_h) {
    const int N = tree.getNumberOfBranches();

    // Walk root-to-leaves: a child's vol_summed reads its parent's vol_summed.
    // In depth-first order, parents always precede their descendants.
    for (int i = 0; i < N; ++i) {
        Branch* b = tree.getBranch(i);
        const double d = b->getDiameter();
        const double vol_actual = 0.25 * PI * d * d * b->getLength();
        const double maintenance_vol = PI * d * maintenance_h;
        b->setMaintenanceVol(maintenance_vol);

        const double s = safety.compute(b->getNbLeaves(), b->getMaxStress());
        const double vol_wished = s * vol_actual *
                                  std::pow(b->getMaxStress(), 2.0 / 3.0);
        const double vol_growth = std::max(0.0, vol_wished - vol_actual) + maintenance_vol;
        b->setVolGrowth(vol_growth);

        const int parent_index = tree.getParentIndex(i);
        if (parent_index < 0) {
            b->setVolSummed(vol_growth);
        } else {
            const Branch* parent = tree.getBranch(parent_index);
            const int parent_nb = parent->getNbLeaves();
            const double parent_summed = parent->getVolSummed();
            // Guard against div-by-zero in degenerate cases (e.g. a tree
            // where reorder() has not been called).
            const double share = (parent_nb > 0)
                ? parent_summed * static_cast<double>(b->getNbLeaves())
                                / static_cast<double>(parent_nb)
                : 0.0;
            b->setVolSummed(share + vol_growth);
        }
    }
}


void secondary_growth(Tree& tree, double volume_per_leaf) {
    const auto leaves = tree.leafIndices();
    for (int leaf_idx : leaves) {
        Branch* leaf = tree.getBranch(leaf_idx);
        const double leaf_vol_summed = leaf->getVolSummed();
        const double photosynth = leaf->getLight() * volume_per_leaf;
        const double vol_growth_branches = std::min(photosynth, leaf_vol_summed);
        tree.addReserve(photosynth - vol_growth_branches);

        if (leaf_vol_summed <= 0.0 || leaf->getNbLeaves() <= 0) {
            // Nothing to grow — degenerate input.
            continue;
        }

        // Leaf itself.
        {
            const double fraction = leaf->getVolGrowth()
                / static_cast<double>(leaf->getNbLeaves())
                / leaf_vol_summed;
            const double growth = vol_growth_branches * fraction
                - leaf->getMaintenanceVol() / static_cast<double>(leaf->getNbLeaves());
            const double d0 = leaf->getDiameter();
            const double new_d2 = growth / leaf->getLength() * 4.0 / PI + d0 * d0;
            // Fortran uses 1e-3 floor on the leaf and 0 on ancestors. Same here.
            leaf->setDiameter(std::sqrt(std::max(1e-3, new_d2)));
        }

        // Walk leaf-to-root.
        int cur = leaf_idx;
        while (true) {
            const int parent_idx = tree.getParentIndex(cur);
            if (parent_idx < 0) break;
            Branch* b = tree.getBranch(parent_idx);
            if (b->getNbLeaves() <= 0 || leaf_vol_summed <= 0.0) {
                cur = parent_idx;
                continue;
            }
            const double fraction = b->getVolGrowth()
                / static_cast<double>(b->getNbLeaves())
                / leaf_vol_summed;
            const double growth = vol_growth_branches * fraction
                - b->getMaintenanceVol() / static_cast<double>(b->getNbLeaves());
            const double d0 = b->getDiameter();
            const double new_d2 = growth / b->getLength() * 4.0 / PI + d0 * d0;
            b->setDiameter(std::sqrt(std::max(0.0, new_d2)));
            cur = parent_idx;
        }
    }
}


int primary_growth(
        Tree& tree,
        const AllocationModel& alloc,
        double twig_length, double twig_diameter,
        double theta1, double theta2,
        double gamma1, double gamma2,
        int generation) {
    (void)generation;  // mirrored from Fortran signature; currently unused
                       // on the C++ side (no per-branch generation field).

    const double volume_twig = 0.25 * PI * twig_length * twig_diameter * twig_diameter;
    if (volume_twig <= 0.0) return 0;

    const auto leaf_indices = tree.leafIndices();
    const int n_leaves = static_cast<int>(leaf_indices.size());
    if (n_leaves == 0) return 0;

    const double vol_relative = tree.getReserve()
                                / static_cast<double>(n_leaves) / volume_twig;
    double p_seeds = 0.0, p_leaves = 0.0, phototropism = 0.0;
    alloc.compute(n_leaves, vol_relative, p_seeds, p_leaves, phototropism);

    const int n_new_leaves = static_cast<int>(std::floor(
        p_leaves * tree.getReserve() / (2.0 * volume_twig)));

    // Sort leaves by light (ascending — dimmest first). `sorted_pos` maps the
    // rank-position (0 = dimmest) to its index in `leaf_indices`.
    std::vector<int> sorted_pos(n_leaves);
    for (int i = 0; i < n_leaves; ++i) sorted_pos[i] = i;
    std::sort(sorted_pos.begin(), sorted_pos.end(), [&](int a, int b) {
        return tree.getBranch(leaf_indices[a])->getLight()
             < tree.getBranch(leaf_indices[b])->getLight();
    });

    const int n = std::min(n_leaves, n_new_leaves);
    const int k = n + static_cast<int>(std::lround(
        (1.0 - phototropism) * static_cast<double>(n_leaves - n)));

    if (k <= 0 || n == 0) return 0;

    // Random permutation perm[0..k-1] over [1..k], matching Fortran rperm2.
    std::vector<int> perm(static_cast<std::size_t>(k));
    for (int i = 0; i < k; ++i) perm[static_cast<std::size_t>(i)] = i + 1;
    std::shuffle(perm.begin(), perm.end(), tree.rng());

    std::normal_distribution<double> noise(0.0, 1.0);
    auto& rng = tree.rng();

    // Stage decisions before mutating the tree — `addBranchWithGeometry`
    // shifts indices > parent_index, so we process spawns in descending
    // parent-index order to keep earlier indices valid.
    struct Spawn {
        int parent_index;
        double t1, t2, g1, g2;
    };
    std::vector<Spawn> spawns;
    spawns.reserve(static_cast<std::size_t>(n));

    double remaining_reserve = tree.getReserve();
    for (int i = 0; i < n; ++i) {
        // Fortran: j = ind(n_leaves - perm[i] + 1) — pick the (perm[i])-th
        // brightest leaf (1 = brightest).
        const int pick = perm[static_cast<std::size_t>(i)];
        const int rank_pos = n_leaves - pick;
        if (rank_pos < 0 || rank_pos >= n_leaves) continue;
        const int leaf_array_idx = sorted_pos[static_cast<std::size_t>(rank_pos)];
        const int leaf_branch_idx = leaf_indices[static_cast<std::size_t>(leaf_array_idx)];

        if (remaining_reserve < 2.0 * volume_twig) break;
        remaining_reserve -= 2.0 * volume_twig;

        Branch* leaf = tree.getBranch(leaf_branch_idx);
        // "Above ground" check — Fortran formula preserved verbatim.
        const auto& loc = leaf->getLocation();
        const auto& ut = leaf->getUnitT();
        const double tip_z = loc[2] + leaf->getLength() * ut[2];
        const bool ground = (tip_z > 0.99 * leaf->getLength());
        if (!ground) continue;

        const double t1 = theta1 + DEG_TO_RAD_10 * noise(rng);
        const double t2 = theta2 + DEG_TO_RAD_10 * noise(rng);
        const double g1 = gamma1 + DEG_TO_RAD_10 * noise(rng);
        const double g2 = gamma2 + DEG_TO_RAD_10 * noise(rng);
        spawns.push_back({leaf_branch_idx, t1, t2, g1, g2});
    }

    tree.setReserve(remaining_reserve);

    // Descending so earlier insertions don't invalidate later parents.
    std::sort(spawns.begin(), spawns.end(),
              [](const Spawn& a, const Spawn& b) {
                  return a.parent_index > b.parent_index;
              });

    int created = 0;
    for (const auto& s : spawns) {
        Branch* mother = tree.getBranch(s.parent_index);
        const auto m_t = mother->getUnitT();
        const auto m_b = mother->getUnitB();

        std::array<double, 3> t_l, n_l;
        daughter_unit_vectors(m_t, m_b, s.t1, s.g1, t_l, n_l);
        tree.addBranchWithGeometry(
            s.parent_index, twig_length, twig_diameter,
            t_l[0], t_l[1], t_l[2], n_l[0], n_l[1], n_l[2]);
        ++created;

        // After the first insertion, mother is still at s.parent_index; its
        // unit_t / unit_b are unchanged. The first child sat momentarily at
        // s.parent_index + 1; the second insertion shifts it to + 2.
        std::array<double, 3> t_r, n_r;
        daughter_unit_vectors(m_t, m_b, s.t2, s.g2, t_r, n_r);
        tree.addBranchWithGeometry(
            s.parent_index, twig_length, twig_diameter,
            t_r[0], t_r[1], t_r[2], n_r[0], n_r[1], n_r[2]);
        ++created;

        // The mother is no longer a leaf — clear its light (mirrors Fortran).
        tree.getBranch(s.parent_index)->setLight(0.0);
    }

    return created;
}
