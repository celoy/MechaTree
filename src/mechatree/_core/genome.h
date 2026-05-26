/*
 * Genome models — polymorphic branch-function objects evaluated inline per
 * branch by growth/primary_growth. Step 9 ships the constant subclasses; a
 * NeuralSafety / NeuralAllocation port of the Fortran 3-layer tanh networks
 * (mod_tree.f90:735 neural_branch, :771 neural_reserve) is a self-contained
 * later step that drops in by subclassing — no changes to the growth code
 * that calls these.
 */

#ifndef GENOME_H_
#define GENOME_H_

#include <algorithm>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class SafetyModel {
public:
    virtual ~SafetyModel() = default;
    virtual double compute(int nb_leaves, double max_stress) const = 0;
};

class ConstantSafety final : public SafetyModel {
public:
    explicit ConstantSafety(double value) : value_(value) {}
    double compute(int, double) const override { return value_; }
    double value() const { return value_; }
private:
    double value_;
};

class AllocationModel {
public:
    virtual ~AllocationModel() = default;
    // Writes p_seeds, p_leaves, phototropism via out-references. The Fortran
    // reference (neural_reserve) returns three values in [0, 1], with the
    // normalisation p_seeds + p_leaves <= 1 enforced by the caller (see the
    // growth port). This contract is mirrored here.
    virtual void compute(int nb_leaves, double vol_relative,
                         double& p_seeds, double& p_leaves,
                         double& phototropism) const = 0;
};

class ConstantAllocation final : public AllocationModel {
public:
    ConstantAllocation(double p_seeds, double p_leaves, double phototropism)
        : p_seeds_(p_seeds), p_leaves_(p_leaves), phototropism_(phototropism) {}
    void compute(int, double,
                 double& p_seeds, double& p_leaves,
                 double& phototropism) const override {
        p_seeds = p_seeds_;
        p_leaves = p_leaves_;
        phototropism = phototropism_;
    }
    double pSeeds() const { return p_seeds_; }
    double pLeaves() const { return p_leaves_; }
    double phototropism() const { return phototropism_; }
private:
    double p_seeds_, p_leaves_, phototropism_;
};

// 3-layer tanh networks ported from legacy/fortran/mod_tree.f90:735 (neural_branch)
// and :771 (neural_reserve). Each gene g in [0,1] is decoded into a network
// weight via tan((g - 0.5) * π * 0.99). The 3×2 input matrix M1 has two entries
// pinned to zero by evolutionary constraint (M1[0,0] and M1[2,1]).
// Inputs to the net are (0.01 * nb_leaves, max_stress) for safety and
// (0.01 * nb_leaves, vol_relative) for allocation — the 0.01 scaling lives
// inside compute() so growth.cpp keeps passing raw nb_leaves.

class NeuralSafety final : public SafetyModel {
public:
    static constexpr int N_WEIGHTS = 10;
    explicit NeuralSafety(const double* weights) {
        for (int i = 0; i < N_WEIGHTS; ++i) {
            weights_[i] = weights[i];
            toto_[i] = std::tan((weights[i] - 0.5) * M_PI * 0.99);
        }
        toto_[0] = 0.0;  // M1(0,0) constraint
        toto_[5] = 0.0;  // M1(2,1) constraint  (toto[5] = M1 col-major last entry of col 1)
    }
    double compute(int nb_leaves, double max_stress) const override {
        const double x0 = 0.01 * static_cast<double>(nb_leaves);
        const double x1 = max_stress;
        double Z[3];
        for (int i = 0; i < 3; ++i) {
            Z[i] = toto_[i] * x0 + toto_[3 + i] * x1;
        }
        double Zp[4];
        for (int i = 0; i < 3; ++i) {
            Zp[i] = std::tanh(5.0 * Z[i]) / 3.0;
        }
        Zp[3] = 1.0 / 3.0;
        double F = 0.0;
        for (int j = 0; j < 4; ++j) {
            F += toto_[6 + j] * Zp[j];
        }
        return std::max(0.0, F + 1.0);
    }
    const double* weights() const { return weights_; }
private:
    double weights_[N_WEIGHTS];
    double toto_[N_WEIGHTS];
};

class NeuralAllocation final : public AllocationModel {
public:
    static constexpr int N_WEIGHTS = 18;
    explicit NeuralAllocation(const double* weights) {
        for (int i = 0; i < N_WEIGHTS; ++i) {
            weights_[i] = weights[i];
            toto_[i] = std::tan((weights[i] - 0.5) * M_PI * 0.99);
        }
        toto_[0] = 0.0;  // M1(0,0)
        toto_[5] = 0.0;  // M1(2,1)
    }
    void compute(int nb_leaves, double vol_relative,
                 double& p_seeds, double& p_leaves,
                 double& phototropism) const override {
        const double x0 = 0.01 * static_cast<double>(nb_leaves);
        const double x1 = vol_relative;
        double Z[3];
        for (int i = 0; i < 3; ++i) {
            Z[i] = toto_[i] * x0 + toto_[3 + i] * x1;
        }
        double Zp[4];
        for (int i = 0; i < 3; ++i) {
            Zp[i] = std::tanh(5.0 * Z[i]) / 3.0;
        }
        Zp[3] = 1.0 / 3.0;
        // M2 is (3 outputs, 4 hidden) column-major: M2[i,j] = toto_[6 + j*3 + i].
        // F[0] -> p_leaves, F[1] -> p_seeds, F[2] -> phototropism (matlab order).
        double F[3] = {0.0, 0.0, 0.0};
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 3; ++i) {
                F[i] += toto_[6 + j * 3 + i] * Zp[j];
            }
        }
        double pl = std::min(std::max(0.0, F[0] + 2.0), 4.0) / 4.0;
        double ps = std::min(std::max(0.0, F[1] + 2.0), 4.0) / 4.0;
        const double ph = std::min(std::max(0.0, F[2] + 2.0), 4.0) / 4.0;
        const double s = ps + pl;
        if (s > 1.0) {
            pl /= s;
            ps /= s;
        }
        p_seeds = ps;
        p_leaves = pl;
        phototropism = ph;
    }
    const double* weights() const { return weights_; }
private:
    double weights_[N_WEIGHTS];
    double toto_[N_WEIGHTS];
};

// Callback-driven subclasses — let Python-side code (e.g. a SymPy-compiled
// lambda) plug into the C++ vtable that growth.cpp calls. The callback takes
// the same inputs the abstract `compute()` does, plus an opaque user-data
// pointer (in practice, a borrowed PyObject* that the caller keeps alive).
//
// The Cython wrapper marks its shim `with gil` so we can call back into
// Python safely. The GIL is held throughout requested_growth / primary_growth
// today, so the only cost is the extra wrapping.

// Function-pointer typedefs at namespace scope so Cython's extern declaration
// can name them directly.
typedef double (*safety_callback_fn)(int nb_leaves, double max_stress,
                                     void* user_data);
typedef void (*allocation_callback_fn)(
    int nb_leaves, double vol_relative,
    double* p_seeds, double* p_leaves, double* phototropism,
    void* user_data);

class CallbackSafety final : public SafetyModel {
public:
    CallbackSafety(safety_callback_fn fn, void* user_data)
        : fn_(fn), user_data_(user_data) {}
    double compute(int nb_leaves, double max_stress) const override {
        return fn_(nb_leaves, max_stress, user_data_);
    }
private:
    safety_callback_fn fn_;
    void* user_data_;
};

class CallbackAllocation final : public AllocationModel {
public:
    CallbackAllocation(allocation_callback_fn fn, void* user_data)
        : fn_(fn), user_data_(user_data) {}
    void compute(int nb_leaves, double vol_relative,
                 double& p_seeds, double& p_leaves,
                 double& phototropism) const override {
        fn_(nb_leaves, vol_relative, &p_seeds, &p_leaves, &phototropism, user_data_);
    }
private:
    allocation_callback_fn fn_;
    void* user_data_;
};

#endif
