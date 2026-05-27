/*
 * Declaration of the Branch class.
 *
 * A Branch is the elementary constituent of a Tree: a node carrying named
 * scalar properties, a parent link, and a children list. Branches do NOT own
 * their children — the enclosing Tree owns all Branch* lifetimes (see tree.h).
 */

#ifndef BRANCH_H_
#define BRANCH_H_

#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

class Branch {
private:
    // unordered_map: O(1) average property lookup vs O(log N) for std::map.
    // Iteration order doesn't matter — no caller relies on it.
    std::unordered_map<std::string, double> properties;

    // Single nullable pointer rather than a vector — Branch invariants only
    // ever allow one parent. Saves a heap allocation and indirection per
    // Branch.
    Branch* parent = nullptr;

    std::vector<Branch*> children;
    int strahler = 0;
    int horton = 0;

    // typed mechanics fields (Step 9). Default-initialised so existing
    // dict-built branches keep working without explicit construction.
    double length_ = 1.0;
    double diameter_ = 0.0;
    double light_ = 0.0;
    double stress_ = 0.0;
    double max_stress_ = 0.0;
    double vol_growth_ = 0.0;
    double vol_summed_ = 0.0;
    double maintenance_vol_ = 0.0;
    int    nb_leaves_ = 0;
    std::array<double, 3> location_{{0.0, 0.0, 0.0}};
    std::array<double, 3> unit_t_{{0.0, 0.0, 1.0}};
    std::array<double, 3> unit_b_{{1.0, 0.0, 0.0}};
    std::array<double, 3> force_{{0.0, 0.0, 0.0}};
    std::array<double, 3> moment_{{0.0, 0.0, 0.0}};
    // Step 25c: per-branch wind state pre-stored by the momentum-wind
    // bridge (option B). `segment_force_` is the branch's own woody drag
    // F_seg from the CFD; `segment_wind_` is the local wind it felt
    // (magnitude from the CFD, along the storm direction). Both are kept
    // separate from `force_`/`moment_` because those slots are reused as
    // the leaves-to-trunk aggregation accumulator during prune.
    std::array<double, 3> segment_force_{{0.0, 0.0, 0.0}};
    std::array<double, 3> segment_wind_{{0.0, 0.0, 0.0}};

public:
    Branch();
    ~Branch() = default;

    // properties
    void addProperty(const std::string& name, double value);
        // throws std::invalid_argument if `name` already exists
    void setProperty(const std::string& name, double value);
        // throws std::out_of_range if `name` does not exist
    void setProperties(const std::unordered_map<std::string, double>& props);
    double getProperty(const std::string& name) const;
        // throws std::out_of_range if `name` does not exist

    // parent
    bool hasParent() const;
    bool hasParent(const Branch* test) const;
    void addParent(Branch* p);
        // throws std::invalid_argument if this branch already has a parent
    void removeParent(const Branch* p);
        // throws std::invalid_argument if `p` is not the current parent
    void removeParent();
        // no-op if this branch has no parent (e.g. the trunk)
    Branch* getParent() const;
        // returns nullptr if this branch has no parent

    // children
    bool hasChildren() const;
    bool hasChild(const Branch* test) const;
    void addChild(Branch* ch);
        // throws std::invalid_argument if `ch` is already a child
    void removeChild(const Branch* ch);
        // throws std::invalid_argument if `ch` is not a child
    void removeChildren();
    const std::vector<Branch*>& getChildren() const;
    std::vector<Branch*> getBrothers() const;

    // hierarchy
    void setStrahler(int strahler_index);
    int getStrahler() const;
    void setHorton(int horton_index);
    int getHorton() const;

    // typed mechanics fields (Step 9)
    //
    // The property map above is the user-extension bag; these typed fields
    // are the hot path read/written by mechanics, growth and pruning. The
    // two storages are independent: a branch built with a `{"length": 1.0}`
    // dict has `getProperty("length") == 1.0` and `getLength() == 0.0`.

    double getLength() const { return length_; }
    void   setLength(double v) { length_ = v; }
    double getDiameter() const { return diameter_; }
    void   setDiameter(double v) { diameter_ = v; }
    double getLight() const { return light_; }
    void   setLight(double v) { light_ = v; }
    double getStress() const { return stress_; }
    void   setStress(double v) { stress_ = v; }
    double getMaxStress() const { return max_stress_; }
    void   setMaxStress(double v) { max_stress_ = v; }
    double getVolGrowth() const { return vol_growth_; }
    void   setVolGrowth(double v) { vol_growth_ = v; }
    double getVolSummed() const { return vol_summed_; }
    void   setVolSummed(double v) { vol_summed_ = v; }
    double getMaintenanceVol() const { return maintenance_vol_; }
    void   setMaintenanceVol(double v) { maintenance_vol_ = v; }
    int    getNbLeaves() const { return nb_leaves_; }
    void   setNbLeaves(int v) { nb_leaves_ = v; }

    // 3D-vector accessors.
    //
    // - `getLocation()` returns a const reference for C++ callers (the
    //   mechanics walks in PR2 read locations thousands of times per
    //   generation).
    // - `locationAt(i)` reads a single component — used by the Cython
    //   wrapper to avoid binding std::array.
    // - `setLocation(x, y, z)` writes three doubles atomically.
    const std::array<double, 3>& getLocation() const { return location_; }
    double locationAt(std::size_t i) const { return location_[i]; }
    void   setLocation(double x, double y, double z) { location_ = {{x, y, z}}; }

    const std::array<double, 3>& getUnitT() const { return unit_t_; }
    double unitTAt(std::size_t i) const { return unit_t_[i]; }
    void   setUnitT(double x, double y, double z) { unit_t_ = {{x, y, z}}; }

    const std::array<double, 3>& getUnitB() const { return unit_b_; }
    double unitBAt(std::size_t i) const { return unit_b_[i]; }
    void   setUnitB(double x, double y, double z) { unit_b_ = {{x, y, z}}; }

    const std::array<double, 3>& getForce() const { return force_; }
    double forceAt(std::size_t i) const { return force_[i]; }
    void   setForce(double x, double y, double z) { force_ = {{x, y, z}}; }

    const std::array<double, 3>& getMoment() const { return moment_; }
    double momentAt(std::size_t i) const { return moment_[i]; }
    void   setMoment(double x, double y, double z) { moment_ = {{x, y, z}}; }

    // Step 25c: pre-stored per-branch wind state (see fields above).
    const std::array<double, 3>& getSegmentForce() const { return segment_force_; }
    double segmentForceAt(std::size_t i) const { return segment_force_[i]; }
    void   setSegmentForce(double x, double y, double z) { segment_force_ = {{x, y, z}}; }

    const std::array<double, 3>& getSegmentWind() const { return segment_wind_; }
    double segmentWindAt(std::size_t i) const { return segment_wind_[i]; }
    void   setSegmentWind(double x, double y, double z) { segment_wind_ = {{x, y, z}}; }
};

#endif
