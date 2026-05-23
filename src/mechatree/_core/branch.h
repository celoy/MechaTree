/*
 * Declaration of the Branch class.
 *
 * A Branch is the elementary constituent of a Tree: a node carrying named
 * scalar properties, a parent link, and a children list. Branches do NOT own
 * their children — the enclosing Tree owns all Branch* lifetimes (see tree.h).
 */

#ifndef BRANCH_H_
#define BRANCH_H_

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
};

#endif
