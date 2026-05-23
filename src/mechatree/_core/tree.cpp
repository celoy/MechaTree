/*
 * Tree class implementation. See tree.h for the public contract.
 */

#include <algorithm>
#include <functional>
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

int Tree::addBranchWithGeometry(
        int parent_index,
        double length, double diameter,
        double unit_t_x, double unit_t_y, double unit_t_z,
        double unit_b_x, double unit_b_y, double unit_b_z) {
    const int N = static_cast<int>(tree_branches.size());
    check_index(parent_index, N, "Tree::addBranchWithGeometry");

    Branch* parent = tree_branches[static_cast<std::size_t>(parent_index)];

    // child.location = parent.location + parent.length * parent.unit_t
    const auto& p_loc = parent->getLocation();
    const auto& p_t   = parent->getUnitT();
    const double p_len = parent->getLength();
    const double cx = p_loc[0] + p_len * p_t[0];
    const double cy = p_loc[1] + p_len * p_t[1];
    const double cz = p_loc[2] + p_len * p_t[2];

    auto* new_branch = new Branch();
    new_branch->setLength(length);
    new_branch->setDiameter(diameter);
    new_branch->setUnitT(unit_t_x, unit_t_y, unit_t_z);
    new_branch->setUnitB(unit_b_x, unit_b_y, unit_b_z);
    new_branch->setLocation(cx, cy, cz);

    parent->addChild(new_branch);
    const int new_index = parent_index + 1;
    tree_branches.insert(tree_branches.begin() + new_index, new_branch);
    shift_indices(new_index, +1);
    return new_index;
}

void Tree::removeBranch(int branch_index) {
    const int N = static_cast<int>(tree_branches.size());
    check_index(branch_index, N, "Tree::removeBranch");

    Branch* branch2remove = tree_branches[static_cast<std::size_t>(branch_index)];

    // CRITICAL: compute the deletion range BEFORE unlinking from the parent.
    // getLastDescendantIndex walks the parent's children to find the next
    // sibling that ends `branch_index`'s subtree. Calling parent.removeChild
    // first can leave the parent with zero children, at which point the
    // recursion in getLastDescendantIndex collapses to the parent's own
    // index (smaller than branch_index) and the deletion loop becomes a
    // no-op — leaving the subtree stranded in tree_branches with a
    // dangling parent pointer.
    const int last_descendant_index = getLastDescendantIndex(branch_index);

    if (branch_index != 0) {
        branch2remove->getParent()->removeChild(branch2remove);
    }

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

void Tree::removeBranches(std::vector<int> branch_indices) {
    std::sort(branch_indices.begin(), branch_indices.end(), std::greater<int>());
    branch_indices.erase(
        std::unique(branch_indices.begin(), branch_indices.end()),
        branch_indices.end());

    // Snapshot the Branch* per requested index BEFORE any removal. Iterating
    // tree_branches[idx] inside the loop is wrong: a previous removeBranch
    // can shrink the vector enough that the SAME idx now points to a
    // different branch. By capturing pointers up front and dispatching on
    // pointer presence in branch_to_index, we tolerate cascading subtree
    // removals (a parent and a descendant both in the input list).
    std::vector<Branch*> targets;
    targets.reserve(branch_indices.size());
    for (int idx : branch_indices) {
        if (idx < 0 || idx >= static_cast<int>(tree_branches.size())) {
            throw std::out_of_range("Tree::removeBranches: index out of range");
        }
        targets.push_back(tree_branches[static_cast<std::size_t>(idx)]);
    }

    for (Branch* b : targets) {
        auto it = branch_to_index.find(b);
        if (it == branch_to_index.end()) {
            // Already cascade-deleted as part of an ancestor's subtree.
            continue;
        }
        removeBranch(it->second);
    }
}

// ---------- nb_leaves / leaf indices -----------------------------------------

int Tree::reorder() {
    // Walk from highest depth-first index back to the trunk. Every branch's
    // descendants sit at higher indices than itself, so by the time we touch
    // index i all of its children's nb_leaves are already finalised.
    const int N = static_cast<int>(tree_branches.size());
    for (int i = N - 1; i >= 0; --i) {
        Branch* b = tree_branches[static_cast<std::size_t>(i)];
        const auto& kids = b->getChildren();
        if (kids.empty()) {
            b->setNbLeaves(1);
        } else {
            int sum = 0;
            for (const Branch* c : kids) {
                sum += c->getNbLeaves();
            }
            b->setNbLeaves(sum);
        }
    }
    return N == 0 ? 0 : tree_branches.front()->getNbLeaves();
}

int Tree::getNbLeaves() const {
    if (tree_branches.empty()) {
        return 0;
    }
    return tree_branches.front()->getNbLeaves();
}

std::vector<int> Tree::leafIndices() const {
    std::vector<int> out;
    const int N = static_cast<int>(tree_branches.size());
    out.reserve(static_cast<std::size_t>(N));
    for (int i = 0; i < N; ++i) {
        if (!tree_branches[static_cast<std::size_t>(i)]->hasChildren()) {
            out.push_back(i);
        }
    }
    return out;
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
