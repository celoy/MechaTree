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

#endif
