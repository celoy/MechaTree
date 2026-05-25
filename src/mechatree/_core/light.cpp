#include "light.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

// One entry per leaf used by the per-direction sort.
//
// `cell_key` is the integer-binned (x_cell, y_cell) packed into one int64
// (upper 32 bits = x_cell, lower 32 bits = y_cell as unsigned). Two leaves
// in the same cell share a key; sort by (key ascending) groups them.
// `neg_z` is `-Zprime` so the largest Z' sorts first within a group.
// `original_idx` is the leaf's row in light_per_direction (also the final
// tie-breaker so the sort is deterministic).
struct LeafEntry {
    std::int64_t cell_key;
    double       neg_z;
    int          original_idx;
};

inline bool entry_less(const LeafEntry& a, const LeafEntry& b) {
    if (a.cell_key != b.cell_key) return a.cell_key < b.cell_key;
    if (a.neg_z    != b.neg_z)    return a.neg_z    < b.neg_z;
    return a.original_idx < b.original_idx;
}

}  // namespace

void light_intercept(const double* leaf_locations,
                     std::size_t n_leaves,
                     const double* sun_elev,
                     const double* sun_azim,
                     std::size_t n_directions,
                     double size_leaf,
                     double leaf_transparency,
                     double* light_per_direction)
{
    if (n_leaves == 0 || n_directions == 0) {
        return;
    }
    const double inv_size = 1.0 / size_leaf;
    const double tau = leaf_transparency;

    // Each direction is independent (same leaf input, distinct output column),
    // so we fan the outer loop out to one direction per thread. Each thread
    // owns its own ``entries`` scratch buffer — there is no cross-thread
    // sharing past the read-only input arrays. When the build lacks OpenMP
    // (e.g. wheels built without libomp), the ``parallel for`` pragma is a
    // no-op and the loop runs serially.
    #pragma omp parallel
    {
        std::vector<LeafEntry> entries(n_leaves);

    #pragma omp for schedule(static) nowait
    for (std::size_t k = 0; k < n_directions; ++k) {
        const double elev = sun_elev[k];
        const double azim = sun_azim[k];
        const double cos_e = std::cos(elev);
        const double sin_e = std::sin(elev);
        const double cos_a = std::cos(azim);
        const double sin_a = std::sin(azim);

        // Pass 1: rotate every leaf into the sun frame, bin, build entry.
        // Matches mod_tree.f90:240–243 verbatim.
        for (std::size_t i = 0; i < n_leaves; ++i) {
            const double X0 = leaf_locations[3 * i + 0];
            const double Y0 = leaf_locations[3 * i + 1];
            const double Z0 = leaf_locations[3 * i + 2];

            const double Xp     =  X0 * cos_a + Y0 * sin_a;
            const double Xprime =  Xp * cos_e + Z0 * sin_e;
            const double Yprime = -X0 * sin_a + Y0 * cos_a;
            const double Zprime = -Xp * sin_e + Z0 * cos_e;

            // Fortran nint(): banker-friendly rounding to nearest integer.
            // std::lround is the standard C++ analogue.
            const std::int64_t x_cell = static_cast<std::int64_t>(
                std::lround(Xprime * inv_size));
            const std::int64_t y_cell = static_cast<std::int64_t>(
                std::lround(Yprime * inv_size));

            // Pack (x_cell, y_cell) into a single int64 key. Each (x, y) maps
            // to a unique key as long as |y_cell| < 2^31 (always true at the
            // forest scales we care about); the high 32 bits carry x_cell
            // sign-extended, the low 32 hold y_cell reinterpreted as
            // unsigned, so distinct pairs never collide.
            const std::int64_t cell_key =
                (x_cell << 32) |
                (static_cast<std::int64_t>(static_cast<std::uint32_t>(y_cell)));

            entries[i] = LeafEntry{ cell_key, -Zprime, static_cast<int>(i) };
        }

        // Sort by (cell_key, neg_z, original_idx). std::sort is introsort
        // (O(n log n)); the comparator is a hot inner loop so we keep it
        // small and branch-predictable.
        std::sort(entries.begin(), entries.end(), entry_less);

        // Walk the sorted entries, assign in-cell depth (topmost leaf gets
        // depth 0 → light = 1; the i-th below it gets light = tau^i), and
        // write directly into the output buffer in original-leaf order.
        //
        // Running-multiply trick: since depth advances by exactly 1 per step
        // within a cell, ``light *= tau`` gives ``tau^depth`` without any
        // ``std::pow`` calls. At island scale this eliminates ~14M pow()
        // calls per generation.
        std::int64_t current_key = entries[0].cell_key - 1;  // sentinel mismatch
        double light = 1.0;
        for (std::size_t s = 0; s < n_leaves; ++s) {
            if (entries[s].cell_key != current_key) {
                current_key = entries[s].cell_key;
                light = 1.0;
            } else {
                light *= tau;
            }
            light_per_direction[entries[s].original_idx * n_directions + k] = light;
        }
    }  // end of #pragma omp for
    }  // end of #pragma omp parallel
}
