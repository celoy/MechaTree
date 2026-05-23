#cython: language_level=3
#cython: wraparound=False
#cython: boundscheck=False
#cython: cdivision=True
#distutils: language=c++

"""Python wrapping of the C++ Tree class. See tree.h for the C++ contract."""

from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

from mechatree._core.cytree cimport Tree


cdef class PyTree:
    """Python view of a C++ Tree.

    Bounds violations and missing-property lookups raise Python exceptions
    (`IndexError` / `KeyError` analogues surfaced as `RuntimeError` from the
    C++ side via Cython's `except +` machinery).
    """

    cdef Tree* c_tree

    def __cinit__(self, dict trunk_props):
        cdef unordered_map[string, double] properties
        cdef bytes bkey
        for key, value in trunk_props.items():
            bkey = key.encode("utf-8")
            properties[bkey] = value
        self.c_tree = new Tree(properties)

    def __dealloc__(self):
        # Pair the `new Tree(...)` in __cinit__ with a matching delete; without
        # this every PyTree leaks its entire branch graph on garbage collection.
        if self.c_tree is not NULL:
            del self.c_tree
            self.c_tree = NULL

    # ---------- properties ---------------------------------------------------

    def add_property(self, int index, str name, double value):
        """Add a new property to the branch at `index`. Raises if `name` already exists."""
        self.c_tree.addProperty(index, name.encode("utf-8"), value)

    def get_property(self, int index, str name):
        """Return the value of property `name` on the branch at `index`."""
        return self.c_tree.getProperty(index, name.encode("utf-8"))

    def set_property(self, int index, str name, double value):
        """Set the (already-existing) property `name` on the branch at `index`."""
        self.c_tree.setProperty(index, name.encode("utf-8"), value)

    # ---------- topology -----------------------------------------------------

    def get_number_of_branches(self):
        """Total number of branches in the tree."""
        return self.c_tree.getNumberOfBranches()

    def add_branch(self, int index, dict props):
        """Add a new branch as a child of the branch at `index`.

        The new branch is inserted at `index + 1` in the depth-first ordering;
        existing branches with greater indices are shifted up by one.
        """
        cdef unordered_map[string, double] properties
        cdef bytes bkey
        for key, value in props.items():
            bkey = key.encode("utf-8")
            properties[bkey] = value
        self.c_tree.addBranch(index, properties)

    def remove_branch(self, int index):
        """Remove the branch at `index` and the entire subtree rooted at it."""
        self.c_tree.removeBranch(index)

    # ---------- family relations --------------------------------------------

    def get_last_descendant_index(self, int ancestor_index):
        """Index of the last descendant of the branch at `ancestor_index`."""
        return self.c_tree.getLastDescendantIndex(ancestor_index)

    def get_parent_index(self, int child_index):
        """Index of the parent of the branch at `child_index`, or -1 if it is the trunk."""
        return self.c_tree.getParentIndex(child_index)

    def get_brothers_index(self, int index):
        """List of indices of the brothers of the branch at `index`.

        Returns an empty list when the branch has no brothers (e.g. the trunk
        or an only child).
        """
        cdef vector[int] brothers = self.c_tree.getBrothersIndex(index)
        return [brothers[i] for i in range(brothers.size())]

    def get_children_index(self, int index):
        """List of indices of the children of the branch at `index`.

        Returns an empty list when the branch has no children (a leaf).
        """
        cdef vector[int] children = self.c_tree.getChildrenIndex(index)
        return [children[i] for i in range(children.size())]

    def has_parent(self, int index):
        """True if the branch at `index` has a parent (i.e. is not the trunk)."""
        return self.c_tree.hasParent(index)

    def get_number_of_children(self, int parent_index):
        """Number of direct children of the branch at `parent_index`."""
        return self.c_tree.getNumberOfChildren(parent_index)

    # ---------- Strahler / Horton classification -----------------------------

    def set_strahler(self):
        """Compute the Strahler order of every branch."""
        self.c_tree.setStrahler()

    def get_strahler(self, int index):
        """Strahler order of the branch at `index`."""
        return self.c_tree.getStrahler(index)

    def get_strahler_distribution(self):
        """Dict mapping Strahler order -> number of branches of that order."""
        cdef map[int, int] dist = self.c_tree.getStrahlerDistribution()
        return {it.first: it.second for it in dist}

    def set_horton(self):
        """Compute the Horton order of every branch."""
        self.c_tree.setHorton()

    def get_horton(self, int index):
        """Horton order of the branch at `index`."""
        return self.c_tree.getHorton(index)

    def get_horton_distribution(self):
        """Dict mapping Horton order -> number of branches of that order."""
        cdef map[int, int] dist = self.c_tree.getHortonDistribution()
        return {it.first: it.second for it in dist}

    def mean_agg_prop_s(self, str name):
        """Dict mapping Strahler order -> mean of property `name` over that order."""
        cdef map[int, double] means = self.c_tree.meanAggregativePropS(name.encode("utf-8"))
        return {it.first: it.second for it in means}

    def mean_agg_prop_h(self, str name):
        """Dict mapping Horton order -> mean of property `name` over that order."""
        cdef map[int, double] means = self.c_tree.meanAggregativePropH(name.encode("utf-8"))
        return {it.first: it.second for it in means}
