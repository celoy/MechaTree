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
    CallbackAllocation,
    CallbackSafety,
    ConstantAllocation,
    ConstantSafety,
    NeuralAllocation,
    NeuralSafety,
    SafetyModel,
    Tree,
    allocation_callback_fn,
    array3d,
    cpp_calculate_stresses,
    cpp_calculate_stresses_from_stored_forces,
    cpp_light_intercept,
    cpp_momentum_solve,
    cpp_momentum_solve_world,
    cpp_primary_growth,
    cpp_prune,
    cpp_prune_with_stored_forces,
    cpp_requested_growth,
    cpp_secondary_growth,
    cpp_wind_force,
    safety_callback_fn,
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

    def get_segment_force(self, int index):
        """Per-branch CFD woody-segment force (Step 25c, option B)."""
        cdef Branch* b = self.c_tree.getBranch(index)
        return (b.segmentForceAt(0), b.segmentForceAt(1), b.segmentForceAt(2))

    def get_segment_wind(self, int index):
        """Per-branch local CFD wind (Step 25c, option B)."""
        cdef Branch* b = self.c_tree.getBranch(index)
        return (b.segmentWindAt(0), b.segmentWindAt(1), b.segmentWindAt(2))

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

    def get_leaf_tips_batch(self):
        """Return ``(tips, branch_indices)`` for every leaf in depth-first order.

        ``tips`` is an ``(n_leaves, 3)`` ``float64`` ndarray whose rows are
        each leaf branch's tip position (``location + length * unit_t``).
        ``branch_indices`` is an ``(n_leaves,)`` ``int32`` ndarray matching
        :meth:`leaf_indices`.

        Batched replacement for the per-leaf
        ``[get_location, get_unit_t, get_length]`` pattern in
        :func:`mechatree.light.extract_leaves`. Cuts ~85 % off the
        ``extract_leaves`` wall-clock at island scale.
        """
        import numpy as np
        cdef vector[int] leaves = self.c_tree.leafIndices()
        cdef int n = leaves.size()
        cdef int i, bi
        cdef Branch* b
        cdef double L
        tips = np.empty((n, 3), dtype=np.float64)
        branch_idx = np.empty(n, dtype=np.int32)
        cdef double[:, ::1] tips_view = tips
        cdef int[::1] idx_view = branch_idx
        for i in range(n):
            bi = leaves[i]
            b = self.c_tree.getBranch(bi)
            L = b.getLength()
            tips_view[i, 0] = b.locationAt(0) + L * b.unitTAt(0)
            tips_view[i, 1] = b.locationAt(1) + L * b.unitTAt(1)
            tips_view[i, 2] = b.locationAt(2) + L * b.unitTAt(2)
            idx_view[i] = bi
        return tips, branch_idx

    def set_lights_batch(self, branch_indices, values):
        """Vectorised :meth:`set_light` for a batch of branches.

        ``branch_indices`` and ``values`` are 1-D arrays of equal length.
        Replacement for the per-leaf Python loop in
        :func:`mechatree.light.aggregate_onto_trees`.
        """
        import numpy as np
        idx_arr = np.ascontiguousarray(branch_indices, dtype=np.int32)
        val_arr = np.ascontiguousarray(values, dtype=np.float64)
        cdef int[::1] idx_view = idx_arr
        cdef double[::1] val_view = val_arr
        cdef int n = idx_view.shape[0]
        if val_view.shape[0] != n:
            raise ValueError(
                f"set_lights_batch: branch_indices and values length mismatch "
                f"({n} vs {val_view.shape[0]})"
            )
        cdef int i
        for i in range(n):
            self.c_tree.getBranch(idx_view[i]).setLight(val_view[i])

    def set_forces_batch(self, forces):
        """Write an ``(N, 3)`` array of per-branch force vectors to the
        tree's branches in depth-first index order.

        Used by the Step-25b momentum-wind wind bridge to plumb the
        per-branch ``F_seg`` vectors computed in the CFD solve back
        into the tree's stress / pruning chain (option B in
        :doc:`/momentum_wind_derivation`), avoiding the per-branch
        recomputation of the force kernel that the legacy
        :func:`calculate_stresses` does.

        ``forces`` must have shape ``(n_branches, 3)``.
        """
        import numpy as np
        arr = np.ascontiguousarray(forces, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[1] != 3:
            raise ValueError(
                f"set_forces_batch: forces must have shape (N, 3); got {arr.shape}"
            )
        cdef int n_branches = self.c_tree.getNumberOfBranches()
        if arr.shape[0] != n_branches:
            raise ValueError(
                f"set_forces_batch: forces length must equal n_branches "
                f"({arr.shape[0]} vs {n_branches})"
            )
        cdef double[:, ::1] arr_view = arr
        cdef int i
        cdef Branch* b
        for i in range(n_branches):
            b = self.c_tree.getBranch(i)
            b.setForce(arr_view[i, 0], arr_view[i, 1], arr_view[i, 2])

    def set_segment_forces_batch(self, forces):
        """Write an ``(N, 3)`` array of per-branch woody-segment forces to
        the tree's branches in depth-first index order (Step 25c, option B).

        Unlike :meth:`set_forces_batch` (which writes the ``force_`` slot
        that :func:`prune` reuses as its leaves-to-trunk aggregation
        accumulator), this writes the dedicated ``segment_force_`` field
        that :func:`prune_with_stored_forces` reads as each branch's own
        drag contribution. Populated by the momentum-wind bridge with the
        per-branch ``F_vec`` from the CFD solve.

        ``forces`` must have shape ``(n_branches, 3)``.
        """
        import numpy as np
        arr = np.ascontiguousarray(forces, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[1] != 3:
            raise ValueError(
                f"set_segment_forces_batch: forces must have shape (N, 3); got {arr.shape}"
            )
        cdef int n_branches = self.c_tree.getNumberOfBranches()
        if arr.shape[0] != n_branches:
            raise ValueError(
                f"set_segment_forces_batch: forces length must equal n_branches "
                f"({arr.shape[0]} vs {n_branches})"
            )
        cdef double[:, ::1] arr_view = arr
        cdef int i
        cdef Branch* b
        for i in range(n_branches):
            b = self.c_tree.getBranch(i)
            b.setSegmentForce(arr_view[i, 0], arr_view[i, 1], arr_view[i, 2])

    def set_segment_winds_batch(self, winds):
        """Write an ``(N, 3)`` array of per-branch local wind vectors to the
        tree's branches in depth-first index order (Step 25c, option B).

        Stored in the dedicated ``segment_wind_`` field;
        :func:`prune_with_stored_forces` uses it for the leaf-cluster drag
        term on terminal branches (``leaf_drag_S0 * |w| * w``). Populated by
        the momentum-wind bridge with each branch's local CFD wind
        (magnitude ``U_branch`` along the storm direction).

        ``winds`` must have shape ``(n_branches, 3)``.
        """
        import numpy as np
        arr = np.ascontiguousarray(winds, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[1] != 3:
            raise ValueError(
                f"set_segment_winds_batch: winds must have shape (N, 3); got {arr.shape}"
            )
        cdef int n_branches = self.c_tree.getNumberOfBranches()
        if arr.shape[0] != n_branches:
            raise ValueError(
                f"set_segment_winds_batch: winds length must equal n_branches "
                f"({arr.shape[0]} vs {n_branches})"
            )
        cdef double[:, ::1] arr_view = arr
        cdef int i
        cdef Branch* b
        for i in range(n_branches):
            b = self.c_tree.getBranch(i)
            b.setSegmentWind(arr_view[i, 0], arr_view[i, 1], arr_view[i, 2])

    def get_branch_data_batch(self):
        """Return ``(start, axis, diameter, length)`` for every branch in
        depth-first order.

        ``start`` and ``axis`` are ``(n_branches, 3)`` ``float64`` arrays.
        ``diameter`` and ``length`` are ``(n_branches,)`` ``float64``
        arrays. Batched accessor for the per-branch geometry the
        momentum-wind bridge (:class:`mechatree.wind.momentum_wind.MomentumWindBridge`)
        pools across trees each generation, replacing a per-branch Python
        loop over ``[get_location, get_unit_t, get_diameter, get_length]``.

        The Step-24 fixed-point loop calls the bridge (and hence this)
        once per inner iteration; at 12 k branches the Python loop was
        ~6 ms, this Cython version drops it under 0.5 ms. Mirrors the
        Phase-A pattern from Step 21b (``get_leaf_tips_batch``).
        """
        import numpy as np
        cdef int n = self.c_tree.getNumberOfBranches()
        cdef int i
        cdef Branch* b
        start = np.empty((n, 3), dtype=np.float64)
        axis = np.empty((n, 3), dtype=np.float64)
        diameter = np.empty(n, dtype=np.float64)
        length = np.empty(n, dtype=np.float64)
        cdef double[:, ::1] start_view = start
        cdef double[:, ::1] axis_view = axis
        cdef double[::1] d_view = diameter
        cdef double[::1] L_view = length
        for i in range(n):
            b = self.c_tree.getBranch(i)
            start_view[i, 0] = b.locationAt(0)
            start_view[i, 1] = b.locationAt(1)
            start_view[i, 2] = b.locationAt(2)
            axis_view[i, 0] = b.unitTAt(0)
            axis_view[i, 1] = b.unitTAt(1)
            axis_view[i, 2] = b.unitTAt(2)
            d_view[i] = b.getDiameter()
            L_view[i] = b.getLength()
        return start, axis, diameter, length

    def collapse_single_child_chains(self, double length_max=10.0):
        """Fuse every maximal single-child run into one straight segment.

        Walks the tree depth-first. Any branch whose only child has only one
        child (and so on) starts a *chain*: a strip of single-child parents.
        Each strip is replaced by one segment running straight from the
        strip's first branch to the strip's tip. Volume (sum of
        ``pi/4 * d**2 * L`` over the merged branches) is preserved and the
        new diameter is back-solved from the new (Euclidean) length.

        ``length_max`` (default 10.0) caps the resulting merged segment's
        Euclidean length. A chain is only extended while the prospective
        merged length stays at or below ``length_max``; beyond that point
        the chain is truncated and the next-down branch becomes the merged
        segment's child.

        Returns the number of branches absorbed. Does *not* call
        ``reorder()`` — the caller is expected to do so before the next
        phase that reads ``nb_leaves`` (e.g. ``requested_growth``).

        O(N) — walks every branch. Prefer ``collapse_chains_after_prune``
        for the steady-state case where only a handful of chains formed in
        the most recent pruning pass.
        """
        return self.c_tree.collapseSingleChildChains(length_max)

    def collapse_chains_after_prune(self, double length_max=10.0):
        """Targeted variant: collapse only the chains seeded by the
        most recent ``prune`` call.

        ``prune`` records the parents of every cut subtree on the tree.
        This method walks up from each of those parents to the chain start
        and runs the same merge logic as ``collapse_single_child_chains``,
        with the same ``length_max`` semantics. O(P + total chain length)
        where P is the cut-parent set's size — typically a handful per
        generation versus the whole tree.

        Returns the number of branches absorbed.
        """
        return self.c_tree.collapseChainsAfterPrune(length_max)

    def get_last_prune_parents(self):
        """Parents of the subtrees removed by the most recent ``prune`` call.

        Empty when no pruning has occurred. Indices are computed at call
        time, so they remain valid even if later operations
        (``primary_growth``, ``add_branch_with_geometry``, etc.) have
        shifted indices since the prune. Mainly useful for benchmarking /
        introspection — ``collapse_chains_after_prune`` reads the same
        state internally via pointers.
        """
        cdef vector[int] v = self.c_tree.getLastPruneParentIndices()
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

    def calculate_stresses_from_stored_forces(
            self, double leaf_drag_S0, double cauchy, bint reset_max):
        """One sensing angle's stress pass from the pre-stored per-branch
        forces (Step 26c). Reads ``segment_force_`` / ``segment_wind_`` (set
        by the momentum-wind bridge for this angle), aggregates leaves-to-
        trunk and sets per-branch ``stress``. ``reset_max=True`` seeds
        ``max_stress`` (first sensing angle); ``False`` accumulates the
        per-branch max over angles.
        """
        cpp_calculate_stresses_from_stored_forces(
            dereference(self.c_tree), leaf_drag_S0, cauchy, reset_max)

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

    def prune_with_stored_forces(self, double leaf_drag_S0, double cauchy):
        """Prune using the per-branch forces pre-stored on each branch.

        Step 25c (option B). Reads each branch's ``segment_force_`` (the
        woody-segment drag from the momentum-wind CFD) and ``segment_wind_``
        (its local wind, used for the leaf-cluster drag term) instead of
        recomputing the drag from a single canopy-mean wind. The caller
        (the momentum-wind bridge) must have populated those fields for the
        current branch set via :meth:`set_segment_forces_batch` /
        :meth:`set_segment_winds_batch` before calling this.

        Returns the number of branches removed (including all descendants of
        any directly-cut branch).
        """
        return cpp_prune_with_stored_forces(
            dereference(self.c_tree), leaf_drag_S0, cauchy)

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

    def compute(self, int nb_leaves, double max_stress) -> float:
        """Evaluate the safety factor for a single branch. Dispatches through
        the C++ vtable to the concrete subclass."""
        if self._model is NULL:
            raise ValueError("PySafetyModel is not initialised; use a subclass")
        return self._model.compute(nb_leaves, max_stress)


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

    def compute(self, int nb_leaves, double vol_relative):
        """Evaluate the allocation model for a tree.

        Returns ``(p_seeds, p_leaves, phototropism)``. Dispatches through the
        C++ vtable to the concrete subclass.
        """
        cdef double p_seeds = 0.0
        cdef double p_leaves = 0.0
        cdef double phototropism = 0.0
        if self._model is NULL:
            raise ValueError("PyAllocationModel is not initialised; use a subclass")
        self._model.compute(nb_leaves, vol_relative, p_seeds, p_leaves, phototropism)
        return (p_seeds, p_leaves, phototropism)


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


cdef class PyNeuralSafety(PySafetyModel):
    """3-layer tanh `neural_branch` network (10 weights, gene-domain [0, 1]).

    Mirrors ``legacy/fortran/mod_tree.f90:735``. Each gene is decoded with
    ``tan((g - 0.5) * pi * 0.99)`` then folded into M1 (3x2) and M2 (1x4)
    Fortran-order matrices with two entries pinned to zero by evolutionary
    constraint. The 0.01 scaling on ``nb_leaves`` lives inside the C++ forward
    pass; growth.cpp keeps passing raw integer leaf counts.
    """

    def __cinit__(self, weights):
        cdef double[10] buf
        cdef int i
        if len(weights) != 10:
            raise ValueError(
                f"PyNeuralSafety expects 10 weights, got {len(weights)}")
        for i in range(10):
            buf[i] = float(weights[i])
        self._model = new NeuralSafety(buf)

    @property
    def weights(self):
        import numpy as np
        cdef const double* w = (<NeuralSafety*>self._model).weights()
        cdef int i
        out = np.empty(10, dtype=np.float64)
        for i in range(10):
            out[i] = w[i]
        return out


cdef class PyNeuralAllocation(PyAllocationModel):
    """3-layer tanh `neural_reserve` network (18 weights, gene-domain [0, 1]).

    Mirrors ``legacy/fortran/mod_tree.f90:771``. 3-output net feeding back
    ``(p_seeds, p_leaves, phototropism)`` with the per-output clip and the
    ``p_seeds + p_leaves > 1`` renormalisation from the Fortran reference.
    """

    def __cinit__(self, weights):
        cdef double[18] buf
        cdef int i
        if len(weights) != 18:
            raise ValueError(
                f"PyNeuralAllocation expects 18 weights, got {len(weights)}")
        for i in range(18):
            buf[i] = float(weights[i])
        self._model = new NeuralAllocation(buf)

    @property
    def weights(self):
        import numpy as np
        cdef const double* w = (<NeuralAllocation*>self._model).weights()
        cdef int i
        out = np.empty(18, dtype=np.float64)
        for i in range(18):
            out[i] = w[i]
        return out


# ----------------------------------------------------------------------------
# Callback-driven safety / allocation (Step 15)
#
# These thread a Python callable through the C++ vtable that growth.cpp calls.
# The C++ side stores a function pointer + an opaque user-data void*; the
# Cython shim casts that void* back to the Python object and invokes it. The
# `with gil` annotation lets the shim re-acquire the GIL if growth is ever
# called from a `with nogil:` block in the future (today it is not).
# ----------------------------------------------------------------------------

cdef double _safety_callback(int nb_leaves, double max_stress,
                             void* user_data) noexcept with gil:
    cdef object py_callable = <object>user_data
    try:
        return float(py_callable(nb_leaves, max_stress))
    except BaseException as exc:
        # The C++ growth loop has no error-propagation channel; the safest
        # fallback is to return a benign value and let the Python caller pick
        # up the exception on the next callable.compute() they make directly.
        # Print to stderr so the failure isn't silent.
        import sys
        sys.stderr.write(
            f"PyCallbackSafety callback raised: {type(exc).__name__}: {exc}\n"
        )
        return 0.0


cdef void _allocation_callback(int nb_leaves, double vol_relative,
                               double* p_seeds, double* p_leaves,
                               double* phototropism,
                               void* user_data) noexcept with gil:
    cdef object py_callable = <object>user_data
    cdef tuple result
    try:
        result = tuple(py_callable(nb_leaves, vol_relative))
    except BaseException as exc:
        import sys
        sys.stderr.write(
            f"PyCallbackAllocation callback raised: "
            f"{type(exc).__name__}: {exc}\n"
        )
        p_seeds[0] = 0.0
        p_leaves[0] = 0.0
        phototropism[0] = 0.0
        return
    if len(result) != 3:
        import sys
        sys.stderr.write(
            f"PyCallbackAllocation callback returned {len(result)} values; "
            f"expected 3 (p_seeds, p_leaves, phototropism)\n"
        )
        p_seeds[0] = 0.0
        p_leaves[0] = 0.0
        phototropism[0] = 0.0
        return
    p_seeds[0] = float(result[0])
    p_leaves[0] = float(result[1])
    phototropism[0] = float(result[2])


cdef class PyCallbackSafety(PySafetyModel):
    """Safety model backed by an arbitrary Python callable.

    ``fn`` is invoked as ``fn(nb_leaves: int, max_stress: float) -> float``
    every time the C++ growth loop needs a safety factor. Use this to plug
    in a SymPy-compiled expression (:func:`mechatree.sympy_genome.sympy_safety`),
    a hand-written closure, or any other Python decision function.

    The callable is kept alive by the wrapper for the lifetime of the model.
    """

    cdef object _py_callable

    def __cinit__(self, fn):
        if not callable(fn):
            raise TypeError("PyCallbackSafety expects a callable")
        self._py_callable = fn
        # Pass a borrowed PyObject* — _py_callable keeps it alive.
        self._model = new CallbackSafety(
            <safety_callback_fn>_safety_callback,
            <void*>fn,
        )


cdef class PyCallbackAllocation(PyAllocationModel):
    """Allocation model backed by an arbitrary Python callable.

    ``fn`` is invoked as ``fn(nb_leaves: int, vol_relative: float) ->
    (p_seeds, p_leaves, phototropism)`` every time the C++ allocation step
    needs reserve splits. Returning a sequence other than length-3 logs an
    error and falls back to ``(0, 0, 0)``.
    """

    cdef object _py_callable

    def __cinit__(self, fn):
        if not callable(fn):
            raise TypeError("PyCallbackAllocation expects a callable")
        self._py_callable = fn
        self._model = new CallbackAllocation(
            <allocation_callback_fn>_allocation_callback,
            <void*>fn,
        )


# ---- light interception kernel (Step 21b, Phase B) -------------------------


def light_intercept_kernel(
    double[:, ::1] leaf_locations not None,
    double[::1] sun_elev not None,
    double[::1] sun_azim not None,
    double size_leaf,
    double leaf_transparency,
    double[:, ::1] light_per_direction not None,
):
    """C++-backed body of :func:`mechatree.light.intercept`.

    ``leaf_locations`` must be a contiguous ``(n_leaves, 3)`` float64 view;
    ``sun_elev`` / ``sun_azim`` must each be contiguous length-``n_directions``
    float64 views; ``light_per_direction`` is the contiguous
    ``(n_leaves, n_directions)`` float64 output, written in place.

    Replaces the per-direction NumPy lexsort + unique + diff + repeat dance
    in the Python implementation. Cuts ``intercept`` wall-clock by ~5–10×
    at island scale (Eloy et al., Nat Commun 2017 config).
    """
    cdef Py_ssize_t n_leaves = leaf_locations.shape[0]
    cdef Py_ssize_t n_directions = sun_elev.shape[0]
    if leaf_locations.shape[1] != 3:
        raise ValueError(
            f"leaf_locations must have 3 columns, got {leaf_locations.shape[1]}"
        )
    if sun_azim.shape[0] != n_directions:
        raise ValueError(
            f"sun_elev/sun_azim length mismatch: {n_directions} vs {sun_azim.shape[0]}"
        )
    if light_per_direction.shape[0] != n_leaves or light_per_direction.shape[1] != n_directions:
        raise ValueError(
            f"light_per_direction shape ({light_per_direction.shape[0]}, "
            f"{light_per_direction.shape[1]}) != ({n_leaves}, {n_directions})"
        )
    if not (0.0 <= leaf_transparency <= 1.0):
        raise ValueError(
            f"leaf_transparency must be in [0, 1], got {leaf_transparency}"
        )
    if n_leaves == 0 or n_directions == 0:
        return
    cpp_light_intercept(
        &leaf_locations[0, 0],
        <size_t>n_leaves,
        &sun_elev[0],
        &sun_azim[0],
        <size_t>n_directions,
        size_leaf,
        leaf_transparency,
        &light_per_direction[0, 0],
    )


# ---- momentum-wind kernel (Step 26e) ---------------------------------------


def momentum_solve_kernel(
    double[:, ::1] start not None,
    double[:, ::1] axis not None,
    double[::1] D not None,
    double[::1] L not None,
    double[::1] cell_bounds_x not None,
    double[::1] cell_bounds_y not None,
    double[::1] cell_bounds_z not None,
    double grid_size,
    double[::1] U_infty not None,
    double C_D,
    double nu_diff,
    bint diffusion_per_line,
    double[::1] U_out not None,
    double[::1] U_in_grid not None,
    double[::1] F_D_cell_grid not None,
    double[::1] U_branch not None,
    double[::1] F_N_branch not None,
    double[::1] F_D_branch not None,
    double[:, ::1] F_vec_branch not None,
    double cos_theta=1.0,
    double sin_theta=0.0,
    double[::1] canopy_mean_out=None,
):
    """GIL-free C++ body of one momentum-wind solve (Step 26e).

    Numerically equivalent to
    :func:`mechatree.wind._momentum_wind_kernel.compute_momentum_wind` (the
    readable NumPy reference) to atol 1e-10. The C++ kernel is declared
    ``nogil`` so the independent ``n_sensing_angles`` solves can run on a
    thread pool without GIL contention.

    All array views must be C-contiguous float64. ``start`` / ``axis`` are
    ``(n, 3)``; ``D`` / ``L`` are ``(n,)``; the three ``cell_bounds_*`` are the
    grid boundaries; ``U_infty`` is ``(Nz,)``. The grid outputs ``U_out`` /
    ``U_in_grid`` / ``F_D_cell_grid`` are flat ``(Nz*Ny*Nx,)`` buffers (row-major
    ``(Nz, Ny, Nx)``), written in place; the per-branch outputs are ``(n,)`` /
    ``F_vec_branch`` is ``(n, 3)``.
    """
    cdef Py_ssize_t n = start.shape[0]
    cdef Py_ssize_t nbx = cell_bounds_x.shape[0]
    cdef Py_ssize_t nby = cell_bounds_y.shape[0]
    cdef Py_ssize_t nbz = cell_bounds_z.shape[0]
    if start.shape[1] != 3 or axis.shape[1] != 3:
        raise ValueError("start/axis must have shape (n, 3)")
    if axis.shape[0] != n or D.shape[0] != n or L.shape[0] != n:
        raise ValueError("start/axis/D/L must share the same length n")
    if F_vec_branch.shape[0] != n or F_vec_branch.shape[1] != 3:
        raise ValueError("F_vec_branch must have shape (n, 3)")
    if canopy_mean_out is not None and canopy_mean_out.shape[0] < 1:
        raise ValueError("canopy_mean_out must have length >= 1 when provided")
    # Per-branch pointers are only dereferenced inside the kernel's per-branch
    # loops (skipped when n == 0), so pass NULL for an empty canopy rather than
    # indexing a zero-length view. The grid is still solved → free-stream U_out.
    cdef double* start_ptr = NULL
    cdef double* axis_ptr = NULL
    cdef double* D_ptr = NULL
    cdef double* L_ptr = NULL
    cdef double* U_branch_ptr = NULL
    cdef double* F_N_ptr = NULL
    cdef double* F_D_ptr = NULL
    cdef double* F_vec_ptr = NULL
    cdef double* canopy_ptr = NULL
    if n > 0:
        start_ptr = &start[0, 0]
        axis_ptr = &axis[0, 0]
        D_ptr = &D[0]
        L_ptr = &L[0]
        U_branch_ptr = &U_branch[0]
        F_N_ptr = &F_N_branch[0]
        F_D_ptr = &F_D_branch[0]
        F_vec_ptr = &F_vec_branch[0, 0]
    if canopy_mean_out is not None:
        canopy_ptr = &canopy_mean_out[0]
    with nogil:
        cpp_momentum_solve(
            start_ptr,
            axis_ptr,
            D_ptr,
            L_ptr,
            <size_t>n,
            &cell_bounds_x[0],
            <size_t>nbx,
            &cell_bounds_y[0],
            <size_t>nby,
            &cell_bounds_z[0],
            <size_t>nbz,
            grid_size,
            &U_infty[0],
            C_D,
            nu_diff,
            <int>diffusion_per_line,
            &U_out[0],
            &U_in_grid[0],
            &F_D_cell_grid[0],
            U_branch_ptr,
            F_N_ptr,
            F_D_ptr,
            F_vec_ptr,
            cos_theta,
            sin_theta,
            canopy_ptr,
        )


def momentum_solve_world_kernel(
    double[:, ::1] start not None,
    double[:, ::1] axis not None,
    double[::1] D not None,
    double[::1] L not None,
    double theta,
    double grid_size,
    double pad_x,
    double pad_y,
    double pad_z,
    double U_uniform,
    double ua,
    double z0,
    double kappa,
    double amp,
    double C_D,
    double nu_diff,
    bint diffusion_per_line,
    double[:, ::1] F_world not None,
    double[:, ::1] w_world not None,
):
    """GIL-free consolidated sensing solve (Step 26f).

    Rotation, grid build, inflow, the column march, and the world-frame force
    rotation all happen in C++, so the bridge just pools the *unrotated*
    geometry and passes the storm angle ``theta``. ``U_uniform`` >= 0 selects a
    uniform inlet; a negative value selects the log-law (``ua`` / ``z0`` /
    ``kappa``). Outputs ``F_world`` / ``w_world`` are ``(n, 3)`` written in place.
    """
    cdef Py_ssize_t n = start.shape[0]
    if start.shape[1] != 3 or axis.shape[1] != 3:
        raise ValueError("start/axis must have shape (n, 3)")
    if axis.shape[0] != n or D.shape[0] != n or L.shape[0] != n:
        raise ValueError("start/axis/D/L must share the same length n")
    if F_world.shape[0] != n or F_world.shape[1] != 3:
        raise ValueError("F_world must have shape (n, 3)")
    if w_world.shape[0] != n or w_world.shape[1] != 3:
        raise ValueError("w_world must have shape (n, 3)")
    if n == 0:
        return
    cdef double* start_ptr = &start[0, 0]
    cdef double* axis_ptr = &axis[0, 0]
    with nogil:
        cpp_momentum_solve_world(
            start_ptr,
            axis_ptr,
            &D[0],
            &L[0],
            <size_t>n,
            theta,
            grid_size,
            pad_x,
            pad_y,
            pad_z,
            U_uniform,
            ua,
            z0,
            kappa,
            amp,
            C_D,
            nu_diff,
            <int>diffusion_per_line,
            &F_world[0, 0],
            &w_world[0, 0],
        )
