#cython: language_level=3
#cython: wraparound=False
#cython: boundscheck=False
#cython: cdivision=True
#distutils: language=c++

"""Python wrapping of the C++ Tree class. See tree.h for the C++ contract."""

from cython.operator cimport dereference

from libc.stdint cimport uint32_t
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

from mechatree._core.cytree cimport (
    AllocationModel,
    Branch,
    ConstantAllocation,
    ConstantSafety,
    SafetyModel,
    Tree,
    array3d,
    cpp_calculate_stresses,
    cpp_primary_growth,
    cpp_prune,
    cpp_requested_growth,
    cpp_secondary_growth,
    cpp_wind_force,
)


cdef inline tuple _as_xyz(value, str where):
    """Validate that `value` is a 3-element numeric tuple/sequence."""
    if not isinstance(value, tuple):
        try:
            value = tuple(value)
        except TypeError:
            raise TypeError(f"{where}: expected a 3-element tuple of floats")
    if len(value) != 3:
        raise ValueError(f"{where}: expected a 3-element tuple, got {len(value)}")
    return value


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

    # ---------- typed mechanics fields (Step 9) ------------------------------
    #
    # The `properties` map exposed via add_property / get_property / set_property
    # is the user-extension bag. These typed accessors are the hot path read
    # and written by mechanics, growth and pruning. The two storages are
    # independent — `add_branch(0, {"length": 1.0})` does NOT set the typed
    # length field; use `set_length(idx, 1.0)` or `add_branch_with_geometry`
    # for that.

    def get_length(self, int index):
        return self.c_tree.getBranch(index).getLength()

    def set_length(self, int index, double value):
        self.c_tree.getBranch(index).setLength(value)

    def get_diameter(self, int index):
        return self.c_tree.getBranch(index).getDiameter()

    def set_diameter(self, int index, double value):
        self.c_tree.getBranch(index).setDiameter(value)

    def get_light(self, int index):
        return self.c_tree.getBranch(index).getLight()

    def set_light(self, int index, double value):
        self.c_tree.getBranch(index).setLight(value)

    def get_stress(self, int index):
        return self.c_tree.getBranch(index).getStress()

    def set_stress(self, int index, double value):
        self.c_tree.getBranch(index).setStress(value)

    def get_max_stress(self, int index):
        return self.c_tree.getBranch(index).getMaxStress()

    def set_max_stress(self, int index, double value):
        self.c_tree.getBranch(index).setMaxStress(value)

    def get_vol_growth(self, int index):
        return self.c_tree.getBranch(index).getVolGrowth()

    def set_vol_growth(self, int index, double value):
        self.c_tree.getBranch(index).setVolGrowth(value)

    def get_vol_summed(self, int index):
        return self.c_tree.getBranch(index).getVolSummed()

    def set_vol_summed(self, int index, double value):
        self.c_tree.getBranch(index).setVolSummed(value)

    def get_maintenance_vol(self, int index):
        return self.c_tree.getBranch(index).getMaintenanceVol()

    def set_maintenance_vol(self, int index, double value):
        self.c_tree.getBranch(index).setMaintenanceVol(value)

    def get_nb_leaves(self, int index):
        return self.c_tree.getBranch(index).getNbLeaves()

    def set_nb_leaves(self, int index, int value):
        self.c_tree.getBranch(index).setNbLeaves(value)

    # ---------- 3D vector fields --------------------------------------------

    def get_location(self, int index):
        cdef Branch* b = self.c_tree.getBranch(index)
        return (b.locationAt(0), b.locationAt(1), b.locationAt(2))

    def set_location(self, int index, value):
        cdef tuple xyz = _as_xyz(value, "set_location")
        self.c_tree.getBranch(index).setLocation(
            <double>xyz[0], <double>xyz[1], <double>xyz[2])

    def get_unit_t(self, int index):
        cdef Branch* b = self.c_tree.getBranch(index)
        return (b.unitTAt(0), b.unitTAt(1), b.unitTAt(2))

    def set_unit_t(self, int index, value):
        cdef tuple xyz = _as_xyz(value, "set_unit_t")
        self.c_tree.getBranch(index).setUnitT(
            <double>xyz[0], <double>xyz[1], <double>xyz[2])

    def get_unit_b(self, int index):
        cdef Branch* b = self.c_tree.getBranch(index)
        return (b.unitBAt(0), b.unitBAt(1), b.unitBAt(2))

    def set_unit_b(self, int index, value):
        cdef tuple xyz = _as_xyz(value, "set_unit_b")
        self.c_tree.getBranch(index).setUnitB(
            <double>xyz[0], <double>xyz[1], <double>xyz[2])

    def get_force(self, int index):
        cdef Branch* b = self.c_tree.getBranch(index)
        return (b.forceAt(0), b.forceAt(1), b.forceAt(2))

    def set_force(self, int index, value):
        cdef tuple xyz = _as_xyz(value, "set_force")
        self.c_tree.getBranch(index).setForce(
            <double>xyz[0], <double>xyz[1], <double>xyz[2])

    def get_moment(self, int index):
        cdef Branch* b = self.c_tree.getBranch(index)
        return (b.momentAt(0), b.momentAt(1), b.momentAt(2))

    def set_moment(self, int index, value):
        cdef tuple xyz = _as_xyz(value, "set_moment")
        self.c_tree.getBranch(index).setMoment(
            <double>xyz[0], <double>xyz[1], <double>xyz[2])

    # ---------- tree-level reserve pool --------------------------------------

    def get_reserve(self):
        """Primary-growth reserve pool for the tree."""
        return self.c_tree.getReserve()

    def set_reserve(self, double value):
        self.c_tree.setReserve(value)

    # ---------- mechanics + growth + pruning (Step 9, PR2) -------------------

    def set_seed(self, uint32_t seed):
        """Seed the tree's RNG. A single integer reproduces a run end-to-end."""
        self.c_tree.setSeed(seed)

    def reorder(self):
        """Recompute per-branch `nb_leaves` and return the total leaf count.

        Must be called after any structural change (add_branch_with_geometry,
        prune, primary_growth, etc.) before reading `nb_leaves` or running
        `requested_growth` / `secondary_growth`.
        """
        return self.c_tree.reorder()

    def get_total_leaves(self):
        """Total leaf count (childless branches). Set by `reorder()`."""
        return self.c_tree.getNbLeaves()

    def leaf_indices(self):
        """Indices of childless branches in depth-first order."""
        cdef vector[int] v = self.c_tree.leafIndices()
        return [v[i] for i in range(v.size())]

    def wind_force(self, int index, V):
        """Wind force and moment on a single branch.

        Returns ``(force_xyz, moment_xyz)`` — two 3-tuples.
        """
        cdef tuple xyz = _as_xyz(V, "wind_force (V)")
        cdef array3d Va, force, moment
        Va[0] = <double>xyz[0]
        Va[1] = <double>xyz[1]
        Va[2] = <double>xyz[2]
        cpp_wind_force(
            dereference(self.c_tree.getBranch(index)), Va, force, moment)
        return (
            (force[0], force[1], force[2]),
            (moment[0], moment[1], moment[2]),
        )

    def calculate_stresses(self, double leaf_drag_S0, double cauchy):
        """Sweep four horizontal wind angles; set `max_stress` per branch."""
        cpp_calculate_stresses(dereference(self.c_tree), leaf_drag_S0, cauchy)

    def requested_growth(self, PySafetyModel safety not None, double maintenance_h):
        """Compute per-branch growth requests. Reads `max_stress` and `nb_leaves`."""
        if safety._model is NULL:
            raise ValueError("PySafetyModel is not initialised; use a subclass")
        cpp_requested_growth(
            dereference(self.c_tree),
            dereference(safety._model),
            maintenance_h)

    def secondary_growth(self, double volume_per_leaf):
        """Allocate photosynthate along leaf-to-root chains; grow diameters."""
        cpp_secondary_growth(dereference(self.c_tree), volume_per_leaf)

    def primary_growth(
            self,
            PyAllocationModel alloc not None,
            double twig_length,
            double twig_diameter,
            double theta1,
            double theta2,
            double gamma1,
            double gamma2,
            int generation,
    ):
        """Spawn new twig branches at the most-lit leaves.

        Returns the number of branches added (always even — twigs in pairs).
        """
        if alloc._model is NULL:
            raise ValueError("PyAllocationModel is not initialised; use a subclass")
        return cpp_primary_growth(
            dereference(self.c_tree),
            dereference(alloc._model),
            twig_length, twig_diameter,
            theta1, theta2, gamma1, gamma2,
            generation)

    def prune(self, wind, double leaf_drag_S0, double cauchy):
        """Stochastically remove branches under wind direction `wind`.

        Returns the number of branches removed (including all descendants of
        any directly-cut branch).
        """
        cdef tuple xyz = _as_xyz(wind, "prune (wind)")
        cdef array3d Wa
        Wa[0] = <double>xyz[0]
        Wa[1] = <double>xyz[1]
        Wa[2] = <double>xyz[2]
        return cpp_prune(dereference(self.c_tree), Wa, leaf_drag_S0, cauchy)

    # ---------- geometric branch addition ------------------------------------

    def add_branch_with_geometry(
            self,
            int parent_index,
            double length,
            double diameter,
            unit_t,
            unit_b,
    ):
        """Add a child branch initialised with typed mechanics fields.

        The child's `location` is derived from the parent's
        `location + length * unit_t`. The caller must ensure the parent's
        location, length and unit_t are already set; the trunk is set up by
        the orchestrator before any children are added.

        Returns the new branch's index in depth-first order.
        """
        cdef tuple t_xyz = _as_xyz(unit_t, "add_branch_with_geometry (unit_t)")
        cdef tuple b_xyz = _as_xyz(unit_b, "add_branch_with_geometry (unit_b)")
        return self.c_tree.addBranchWithGeometry(
            parent_index, length, diameter,
            <double>t_xyz[0], <double>t_xyz[1], <double>t_xyz[2],
            <double>b_xyz[0], <double>b_xyz[1], <double>b_xyz[2])


# ============================================================================
# Genome models (Step 9, PR2)
#
# C++ polymorphism is preserved at the Cython boundary: the base class owns
# a `SafetyModel*` / `AllocationModel*` (heap-allocated, virtual destructor).
# Subclasses construct the concrete impl in their __cinit__. PyTree growth
# methods accept the base type and dereference the pointer to a reference.
# ============================================================================

cdef class PySafetyModel:
    """Abstract base — instantiate a concrete subclass like `PyConstantSafety`."""

    cdef SafetyModel* _model

    def __cinit__(self):
        self._model = NULL

    def __init__(self, *args, **kwargs):
        if type(self) is PySafetyModel:
            raise TypeError(
                "PySafetyModel is abstract; instantiate a concrete subclass")

    def __dealloc__(self):
        if self._model is not NULL:
            del self._model
            self._model = NULL


cdef class PyConstantSafety(PySafetyModel):
    """Returns a constant safety factor regardless of inputs.

    Useful for testing and minimum-viable demos; a future `PyNeuralSafety`
    subclass will swap in the Fortran 3-layer tanh network.
    """

    def __cinit__(self, double value):
        self._model = new ConstantSafety(value)

    @property
    def value(self):
        return (<ConstantSafety*>self._model).value()


cdef class PyAllocationModel:
    """Abstract base — instantiate a concrete subclass like `PyConstantAllocation`."""

    cdef AllocationModel* _model

    def __cinit__(self):
        self._model = NULL

    def __init__(self, *args, **kwargs):
        if type(self) is PyAllocationModel:
            raise TypeError(
                "PyAllocationModel is abstract; instantiate a concrete subclass")

    def __dealloc__(self):
        if self._model is not NULL:
            del self._model
            self._model = NULL


cdef class PyConstantAllocation(PyAllocationModel):
    """Returns constant (p_seeds, p_leaves, phototropism) regardless of inputs."""

    def __cinit__(self, double p_seeds, double p_leaves, double phototropism):
        self._model = new ConstantAllocation(p_seeds, p_leaves, phototropism)

    @property
    def p_seeds(self):
        return (<ConstantAllocation*>self._model).pSeeds()

    @property
    def p_leaves(self):
        return (<ConstantAllocation*>self._model).pLeaves()

    @property
    def phototropism(self):
        return (<ConstantAllocation*>self._model).phototropism()
