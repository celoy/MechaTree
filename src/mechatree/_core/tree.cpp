/*
 * Tree class implementation. See tree.h for the public contract.
 */

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "tree.h"

namespace {

// Bounds-check `index` against `n`. Throws std::out_of_range on failure.
void check_index(int index, int n, const char* where) {
    if (index < 0) {
        throw std::out_of_range(std::string(where) + ": index must be non-negative");
    }
    if (index >= n) {
        throw std::out_of_range(std::string(where) + ": index out of range");
    }
}

}  // namespace

// ---------- constructors / destructor ----------------------------------------

Tree::Tree() {
    auto* trunk = new Branch();
    tree_branches.push_back(trunk);
    branch_to_index[trunk] = 0;
}

Tree::Tree(const std::unordered_map<std::string, double>& trunk_props) {
    auto* trunk = new Branch();
    trunk->setProperties(trunk_props);
    tree_branches.push_back(trunk);
    branch_to_index[trunk] = 0;
}

Tree::~Tree() {
    for (Branch* b : tree_branches) {
        delete b;
    }
}

// ---------- index map maintenance --------------------------------------------

void Tree::shift_indices(int from, int delta) {
    // Walk tree_branches from `from` to end; update branch_to_index in lockstep
    // so the map mirrors the vector after a structural change.
    const int n = static_cast<int>(tree_branches.size());
    for (int i = from; i < n; ++i) {
        branch_to_index[tree_branches[static_cast<std::size_t>(i)]] = i;
    }
    (void)delta;  // semantic hint only; the walk above sets the absolute value
}

// ---------- general information ----------------------------------------------

int Tree::getNumberOfBranches() const {
    return static_cast<int>(tree_branches.size());
}

const std::map<int, int>& Tree::getStrahlerDistribution() const {
    return Strahler_distribution;
}

const std::map<int, int>& Tree::getHortonDistribution() const {
    return Horton_distribution;
}

std::map<int, double> Tree::meanAggregativePropS(const std::string& name) const {
    std::map<int, double> mean_prop;
    for (const Branch* b : tree_branches) {
        mean_prop[b->getStrahler()] += b->getProperty(name);
    }
    for (const auto& [order, count] : Strahler_distribution) {
        auto it = mean_prop.find(order);
        if (it != mean_prop.end() && count > 0) {
            it->second /= count;
        }
    }
    return mean_prop;
}

std::map<int, double> Tree::meanAggregativePropH(const std::string& name) const {
    std::map<int, double> mean_prop;
    for (const Branch* b : tree_branches) {
        mean_prop[b->getHorton()] += b->getProperty(name);
    }
    for (const auto& [order, count] : Horton_distribution) {
        auto it = mean_prop.find(order);
        if (it != mean_prop.end() && count > 0) {
            it->second /= count;
        }
    }
    return mean_prop;
}

// ---------- access branches --------------------------------------------------

Branch* Tree::getTrunk() const {
    return tree_branches.front();
}

Branch* Tree::getSummit() const {
    return tree_branches.back();
}

Branch* Tree::getBranch(int branch_index) const {
    check_index(branch_index, static_cast<int>(tree_branches.size()), "Tree::getBranch");
    return tree_branches[static_cast<std::size_t>(branch_index)];
}

int Tree::getIndex(const Branch* branch) const {
    // O(1) average via the unordered_map mirror of tree_branches.
    auto it = branch_to_index.find(branch);
    return it == branch_to_index.end() ? -1 : it->second;
}

// ---------- hierarchy --------------------------------------------------------

void Tree::setStrahler() {
    Strahler_distribution.clear();
    for (auto rit = tree_branches.rbegin(); rit != tree_branches.rend(); ++rit) {
        Branch* b = *rit;
        if (!b->hasChildren()) {
            b->setStrahler(1);
            Strahler_distribution[1] += 1;
            continue;
        }

        const auto& kids = b->getChildren();
        if (kids.size() != 2) {
            throw std::runtime_error(
                "Tree::setStrahler: only binary trees are supported "
                "(branch has " + std::to_string(kids.size()) + " children)");
        }

        const int s1 = kids[0]->getStrahler();
        const int s2 = kids[1]->getStrahler();
        if (s1 == s2) {
            b->setStrahler(s1 + 1);
            Strahler_distribution[s1 + 1] += 1;
        } else {
            b->setStrahler(std::max(s1, s2));
            // contiguous branches of the same order count as one
        }
    }
}

int Tree::getStrahler(int index) const {
    check_index(index, static_cast<int>(tree_branches.size()), "Tree::getStrahler");
    return tree_branches[static_cast<std::size_t>(index)]->getStrahler();
}

void Tree::setHorton() {
    Horton_distribution.clear();
    if (Strahler_distribution.empty()) {
        setStrahler();
    }

    Branch* trunk = getTrunk();
    const int trunk_order = trunk->getStrahler();
    trunk->setHorton(trunk_order);
    Horton_distribution[trunk_order] = 1;

    for (Branch* b : tree_branches) {
        const auto& kids = b->getChildren();
        if (kids.empty()) {
            continue;
        }
        if (kids.size() != 2) {
            throw std::runtime_error(
                "Tree::setHorton: only binary trees are supported "
                "(branch has " + std::to_string(kids.size()) + " children)");
        }

        const int s1 = kids[0]->getStrahler();
        const int s2 = kids[1]->getStrahler();
        const int hortondad = b->getHorton();

        if (s1 >= s2) {
            kids[0]->setHorton(hortondad);
            kids[1]->setHorton(s2);
            Horton_distribution[s2] += 1;
        } else {
            kids[1]->setHorton(hortondad);
            kids[0]->setHorton(s1);
            Horton_distribution[s1] += 1;
        }
        // When the two children have equal order we arbitrarily pass the
        // parent's Horton order down the first child.
    }
}

int Tree::getHorton(int index) const {
    check_index(index, static_cast<int>(tree_branches.size()), "Tree::getHorton");
    return tree_branches[static_cast<std::size_t>(index)]->getHorton();
}

// ---------- family information -----------------------------------------------

int Tree::getLastDescendantIndex(int ancestor_index) const {
    const int N = static_cast<int>(tree_branches.size());
    check_index(ancestor_index, N, "Tree::getLastDescendantIndex");

    Branch* ancestor = tree_branches[static_cast<std::size_t>(ancestor_index)];
    if (ancestor_index == 0) {
        // The trunk's last descendant is always the last branch.
        return N - 1;
    }
    if (!ancestor->hasChildren()) {
        return ancestor_index;
    }

    // Walk siblings of the ancestor: the first sibling whose index is
    // greater than ancestor_index marks (sibling_index - 1) as the last
    // descendant. If no such sibling exists, recurse into the parent.
    const std::vector<int> brothers = getBrothersIndex(ancestor_index);
    for (int brother_index : brothers) {
        if (brother_index > ancestor_index) {
            return brother_index - 1;
        }
    }
    return getLastDescendantIndex(getParentIndex(ancestor_index));
}

int Tree::getParentIndex(int child_index) const {
    const int N = static_cast<int>(tree_branches.size());
    check_index(child_index, N, "Tree::getParentIndex");
    if (child_index == 0) {
        return -1;  // trunk has no parent
    }
    return getIndex(tree_branches[static_cast<std::size_t>(child_index)]->getParent());
}

std::vector<int> Tree::getBrothersIndex(int branch_index) const {
    const int N = static_cast<int>(tree_branches.size());
    check_index(branch_index, N, "Tree::getBrothersIndex");

    std::vector<int> indices;
    for (const Branch* brother : getBranch(branch_index)->getBrothers()) {
        indices.push_back(getIndex(brother));
    }
    return indices;
}

std::vector<int> Tree::getChildrenIndex(int parent_index) const {
    const int N = static_cast<int>(tree_branches.size());
    check_index(parent_index, N, "Tree::getChildrenIndex");

    std::vector<int> indices;
    for (const Branch* child : getBranch(parent_index)->getChildren()) {
        indices.push_back(getIndex(child));
    }
    return indices;
}

bool Tree::hasParent(int index) const {
    check_index(index, static_cast<int>(tree_branches.size()), "Tree::hasParent");
    return tree_branches[static_cast<std::size_t>(index)]->hasParent();
}

int Tree::getNumberOfChildren(int parent_index) const {
    check_index(parent_index, static_cast<int>(tree_branches.size()), "Tree::getNumberOfChildren");
    return static_cast<int>(
        tree_branches[static_cast<std::size_t>(parent_index)]->getChildren().size());
}

// ---------- modify tree ------------------------------------------------------

void Tree::addBranch(int parent_index) {
    const int N = static_cast<int>(tree_branches.size());
    check_index(parent_index, N, "Tree::addBranch");

    Branch* parent = tree_branches[static_cast<std::size_t>(parent_index)];
    auto* new_branch = new Branch();
    parent->addChild(new_branch);
    tree_branches.insert(tree_branches.begin() + parent_index + 1, new_branch);
    shift_indices(parent_index + 1, +1);
}

void Tree::addBranch(int parent_index, const std::unordered_map<std::string, double>& props) {
    const int N = static_cast<int>(tree_branches.size());
    check_index(parent_index, N, "Tree::addBranch");

    Branch* parent = tree_branches[static_cast<std::size_t>(parent_index)];
    auto* new_branch = new Branch();
    new_branch->setProperties(props);
    parent->addChild(new_branch);
    tree_branches.insert(tree_branches.begin() + parent_index + 1, new_branch);
    shift_indices(parent_index + 1, +1);
}

void Tree::removeBranch(int branch_index) {
    const int N = static_cast<int>(tree_branches.size());
    check_index(branch_index, N, "Tree::removeBranch");

    Branch* branch2remove = tree_branches[static_cast<std::size_t>(branch_index)];

    // Unlink from parent's children list while the parent pointer is valid.
    if (branch_index != 0) {
        branch2remove->getParent()->removeChild(branch2remove);
    }

    const int last_descendant_index = getLastDescendantIndex(branch_index);

    // Tree owns the Branch* lifetimes — free the subtree (and drop the index
    // map entries) before erasing from the vector.
    for (int i = branch_index; i <= last_descendant_index; ++i) {
        Branch* b = tree_branches[static_cast<std::size_t>(i)];
        branch_to_index.erase(b);
        delete b;
    }

    if (last_descendant_index == N - 1) {
        tree_branches.erase(tree_branches.begin() + branch_index, tree_branches.end());
    } else {
        tree_branches.erase(
            tree_branches.begin() + branch_index,
            tree_branches.begin() + last_descendant_index + 1);
        // Survivors past the removed range now sit at lower indices.
        shift_indices(branch_index, branch_index - last_descendant_index - 1);
    }
}

// ---------- branch properties ------------------------------------------------

void Tree::addProperty(int index, const std::string& name, double value) {
    check_index(index, static_cast<int>(tree_branches.size()), "Tree::addProperty");
    tree_branches[static_cast<std::size_t>(index)]->addProperty(name, value);
}

double Tree::getProperty(int index, const std::string& name) const {
    check_index(index, static_cast<int>(tree_branches.size()), "Tree::getProperty");
    return tree_branches[static_cast<std::size_t>(index)]->getProperty(name);
}

void Tree::setProperty(int index, const std::string& name, double value) {
    check_index(index, static_cast<int>(tree_branches.size()), "Tree::setProperty");
    tree_branches[static_cast<std::size_t>(index)]->setProperty(name, value);
}
