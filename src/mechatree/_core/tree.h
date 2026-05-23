/*
 * Declaration of the Tree class.
 *
 * Tree owns the Branch* lifetimes — the destructor deletes every branch in
 * `tree_branches`. Branches are stored in a depth-first order: a parent is
 * always immediately followed by its subtree, so [parent_index,
 * last_descendant_index] is a contiguous range.
 *
 * `branch_to_index` mirrors `tree_branches` for O(1) reverse lookup
 * (Branch* -> index). It is maintained by addBranch/removeBranch and is the
 * reason topology walks (getParentIndex, getBrothersIndex, getChildrenIndex)
 * run in O(degree) instead of O(N).
 */

#ifndef TREE_H_
#define TREE_H_

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "branch.h"

class Tree {
private:
    std::vector<Branch*> tree_branches;
    std::unordered_map<const Branch*, int> branch_to_index;
    std::map<int, int> Strahler_distribution;
    std::map<int, int> Horton_distribution;

    // shift branch_to_index entries with stored index in [from, INF) by
    // `delta` (positive when a branch was inserted; negative when erased).
    void shift_indices(int from, int delta);

public:
    Tree();
    explicit Tree(const std::unordered_map<std::string, double>& trunk_props);
    ~Tree();

    // Non-copyable: Tree owns raw Branch* and the default copy would
    // double-free at destruction time.
    Tree(const Tree&) = delete;
    Tree& operator=(const Tree&) = delete;

    // general information
    int getNumberOfBranches() const;
    const std::map<int, int>& getStrahlerDistribution() const;
    const std::map<int, int>& getHortonDistribution() const;
    std::map<int, double> meanAggregativePropS(const std::string& name) const;
    std::map<int, double> meanAggregativePropH(const std::string& name) const;

    // access branches
    Branch* getTrunk() const;
    Branch* getSummit() const;
    Branch* getBranch(int branch_index) const;
        // throws std::out_of_range if `branch_index` is invalid
    int getIndex(const Branch* branch) const;
        // returns -1 if `branch` is not in this tree; O(1)

    // hierarchy
    void setStrahler();
    int getStrahler(int index) const;
        // throws std::out_of_range if `index` is invalid
    void setHorton();
    int getHorton(int index) const;
        // throws std::out_of_range if `index` is invalid

    // family information
    int getLastDescendantIndex(int ancestor_index) const;
        // throws std::out_of_range if `ancestor_index` is invalid
    int getParentIndex(int child_index) const;
        // returns -1 if `child_index` is the trunk;
        // throws std::out_of_range if `child_index` is invalid
    std::vector<int> getBrothersIndex(int branch_index) const;
        // throws std::out_of_range if `branch_index` is invalid
    std::vector<int> getChildrenIndex(int parent_index) const;
        // throws std::out_of_range if `parent_index` is invalid
    bool hasParent(int index) const;
        // throws std::out_of_range if `index` is invalid
    int getNumberOfChildren(int parent_index) const;
        // throws std::out_of_range if `parent_index` is invalid

    // modify tree
    void addBranch(int parent_index);
        // throws std::out_of_range if `parent_index` is invalid
    void addBranch(int parent_index, const std::unordered_map<std::string, double>& props);
        // throws std::out_of_range if `parent_index` is invalid
    void removeBranch(int branch_index);
        // throws std::out_of_range if `branch_index` is invalid

    // branch properties
    void addProperty(int index, const std::string& name, double value);
        // throws std::out_of_range if `index` is invalid
        // throws std::invalid_argument if `name` already exists
    void setProperty(int index, const std::string& name, double value);
        // throws std::out_of_range if `index` is invalid or `name` is missing
    double getProperty(int index, const std::string& name) const;
        // throws std::out_of_range if `index` is invalid or `name` is missing
};

#endif
