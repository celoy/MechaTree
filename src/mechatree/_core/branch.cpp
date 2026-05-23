/*
 * Branch class implementation. See branch.h for the public contract.
 */

#include <algorithm>
#include <stdexcept>
#include <string>

#include "branch.h"

// ---------- constructor ------------------------------------------------------

Branch::Branch() = default;

// ---------- properties -------------------------------------------------------

void Branch::addProperty(const std::string& name, double value) {
    auto [it, inserted] = properties.emplace(name, value);
    if (!inserted) {
        throw std::invalid_argument(
            "Branch::addProperty: property '" + name + "' already exists; use setProperty");
    }
}

void Branch::setProperty(const std::string& name, double value) {
    auto it = properties.find(name);
    if (it == properties.end()) {
        throw std::out_of_range(
            "Branch::setProperty: no such property '" + name + "'");
    }
    it->second = value;
}

void Branch::setProperties(const std::unordered_map<std::string, double>& props) {
    properties = props;
}

double Branch::getProperty(const std::string& name) const {
    auto it = properties.find(name);
    if (it == properties.end()) {
        throw std::out_of_range(
            "Branch::getProperty: no such property '" + name + "'");
    }
    return it->second;
}

// ---------- parent -----------------------------------------------------------

bool Branch::hasParent() const {
    return parent != nullptr;
}

bool Branch::hasParent(const Branch* test) const {
    return parent == test;
}

void Branch::addParent(Branch* p) {
    if (parent != nullptr) {
        throw std::invalid_argument(
            "Branch::addParent: branch already has a parent");
    }
    parent = p;
}

void Branch::removeParent(const Branch* p) {
    if (parent != p) {
        throw std::invalid_argument(
            "Branch::removeParent: pointer is not the current parent");
    }
    parent = nullptr;
}

void Branch::removeParent() {
    parent = nullptr;
}

Branch* Branch::getParent() const {
    return parent;
}

// ---------- children ---------------------------------------------------------

bool Branch::hasChildren() const {
    return !children.empty();
}

bool Branch::hasChild(const Branch* test) const {
    return std::find(children.begin(), children.end(), test) != children.end();
}

void Branch::addChild(Branch* ch) {
    if (hasChild(ch)) {
        throw std::invalid_argument(
            "Branch::addChild: pointer is already a child");
    }
    children.push_back(ch);
    ch->addParent(this);
}

void Branch::removeChild(const Branch* ch) {
    auto it = std::find(children.begin(), children.end(), ch);
    if (it == children.end()) {
        throw std::invalid_argument(
            "Branch::removeChild: pointer is not a child");
    }
    children.erase(it);
}

void Branch::removeChildren() {
    children.clear();
}

const std::vector<Branch*>& Branch::getChildren() const {
    return children;
}

std::vector<Branch*> Branch::getBrothers() const {
    std::vector<Branch*> brothers;
    if (parent == nullptr) {
        return brothers;
    }
    for (Branch* sib : parent->getChildren()) {
        if (sib != this) {
            brothers.push_back(sib);
        }
    }
    return brothers;
}

// ---------- hierarchy --------------------------------------------------------

void Branch::setStrahler(int strahler_index) {
    strahler = strahler_index;
}

int Branch::getStrahler() const {
    return strahler;
}

void Branch::setHorton(int horton_index) {
    horton = horton_index;
}

int Branch::getHorton() const {
    return horton;
}
