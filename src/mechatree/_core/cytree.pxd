#distutils: language=c++

from libc.stdint cimport uint32_t
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

# Named-typedef workaround for std::array<double, 3> — Cython's templated
# array<T, N> syntax doesn't parse cleanly. This buys us a stack-allocatable
# 3-vector with operator[] access.
cdef extern from "<array>" namespace "std" nogil:
    cdef cppclass array3d "std::array<double, 3>":
        double& operator[](size_t)
        array3d()

cdef extern from "genome.h" nogil:
    cdef cppclass SafetyModel:
        double compute(int nb_leaves, double max_stress)

    cdef cppclass ConstantSafety(SafetyModel):
        ConstantSafety(double value) except +
        double value()

    cdef cppclass AllocationModel:
        void compute(int nb_leaves, double vol_relative,
                     double& p_seeds, double& p_leaves, double& phototropism)

    cdef cppclass ConstantAllocation(AllocationModel):
        ConstantAllocation(double p_seeds, double p_leaves, double phototropism) except +
        double pSeeds()
        double pLeaves()
        double phototropism()

    cdef cppclass NeuralSafety(SafetyModel):
        NeuralSafety(const double* weights) except +
        const double* weights()

    cdef cppclass NeuralAllocation(AllocationModel):
        NeuralAllocation(const double* weights) except +
        const double* weights()

    # Callback-driven subclasses — bridge Python callables (e.g. a SymPy-
    # compiled lambda) into the C++ vtable that growth.cpp calls.
    ctypedef double (*safety_callback_fn)(int nb_leaves, double max_stress,
                                          void* user_data)
    ctypedef void   (*allocation_callback_fn)(
        int nb_leaves, double vol_relative,
        double* p_seeds, double* p_leaves, double* phototropism,
        void* user_data)

    cdef cppclass CallbackSafety(SafetyModel):
        CallbackSafety(safety_callback_fn fn, void* user_data) except +

    cdef cppclass CallbackAllocation(AllocationModel):
        CallbackAllocation(allocation_callback_fn fn, void* user_data) except +

cdef extern from "branch.h" nogil:
    cdef cppclass Branch:
        # typed mechanics fields (Step 9) — getters/setters declared inline
        # in branch.h. Only the Cython-friendly per-component accessors and
        # three-double setters are mirrored here; the const-ref getters are
        # for C++ callers only.
        double getLength()
        void   setLength(double v)
        double getDiameter()
        void   setDiameter(double v)
        double getLight()
        void   setLight(double v)
        double getStress()
        void   setStress(double v)
        double getMaxStress()
        void   setMaxStress(double v)
        double getVolGrowth()
        void   setVolGrowth(double v)
        double getVolSummed()
        void   setVolSummed(double v)
        double getMaintenanceVol()
        void   setMaintenanceVol(double v)
        int    getNbLeaves()
        void   setNbLeaves(int v)

        double locationAt(size_t i)
        void   setLocation(double x, double y, double z)
        double unitTAt(size_t i)
        void   setUnitT(double x, double y, double z)
        double unitBAt(size_t i)
        void   setUnitB(double x, double y, double z)
        double forceAt(size_t i)
        void   setForce(double x, double y, double z)
        double momentAt(size_t i)
        void   setMoment(double x, double y, double z)
        double segmentForceAt(size_t i)
        void   setSegmentForce(double x, double y, double z)
        double segmentWindAt(size_t i)
        void   setSegmentWind(double x, double y, double z)

cdef extern from "tree.h" nogil:
    cdef cppclass Tree:
        Tree(unordered_map[string, double] trunk_props) except +

        int getNumberOfBranches()
        const map[int, int]& getStrahlerDistribution()
        const map[int, int]& getHortonDistribution()
        map[int, double] meanAggregativePropS(string name) except +
        map[int, double] meanAggregativePropH(string name) except +

        Branch* getBranch(int branch_index) except +

        void setStrahler() except +
        int getStrahler(int index) except +
        void setHorton() except +
        int getHorton(int index) except +

        int getLastDescendantIndex(int ancestor_index) except +
        int getParentIndex(int child_index) except +
        vector[int] getBrothersIndex(int index) except +
        vector[int] getChildrenIndex(int index) except +
        bint hasParent(int index) except +
        int getNumberOfChildren(int parent_index) except +

        void addBranch(int index, unordered_map[string, double] branch_props) except +
        int  addBranchWithGeometry(
            int parent_index,
            double length, double diameter,
            double unit_t_x, double unit_t_y, double unit_t_z,
            double unit_b_x, double unit_b_y, double unit_b_z) except +
        void removeBranch(int index) except +

        void addProperty(int index, string name, double value) except +
        double getProperty(int index, string name) except +
        void setProperty(int index, string name, double value) except +

        double getReserve()
        void   setReserve(double v)
        void   addReserve(double v)

        void   setSeed(uint32_t seed)
        int    reorder()
        int    getNbLeaves()
        vector[int] leafIndices()
        void   removeBranches(vector[int] indices) except +
        int    collapseSingleChildChains(double length_max) except +
        int    collapseChainsAfterPrune(double length_max) except +
        vector[int] getLastPruneParentIndices()

cdef extern from "mechanics.h" nogil:
    void cpp_wind_force "wind_force"(
        const Branch& b,
        const array3d& V,
        array3d& force,
        array3d& moment) except +
    void cpp_calculate_stresses "calculate_stresses"(
        Tree& tree, double leaf_drag_S0, double cauchy) except +
    void cpp_calculate_stresses_from_stored_forces "calculate_stresses_from_stored_forces"(
        Tree& tree, double leaf_drag_S0, double cauchy, bint reset_max) except +

cdef extern from "growth.h" nogil:
    void cpp_requested_growth "requested_growth"(
        Tree& tree, const SafetyModel& safety, double maintenance_h) except +
    void cpp_secondary_growth "secondary_growth"(
        Tree& tree, double volume_per_leaf) except +
    int  cpp_primary_growth "primary_growth"(
        Tree& tree,
        const AllocationModel& alloc,
        double twig_length, double twig_diameter,
        double theta1, double theta2,
        double gamma1, double gamma2,
        int generation) except +

cdef extern from "pruning.h" nogil:
    int cpp_prune "prune"(
        Tree& tree,
        const array3d& wind,
        double leaf_drag_S0,
        double cauchy) except +
    int cpp_prune_with_stored_forces "prune_with_stored_forces"(
        Tree& tree,
        double leaf_drag_S0,
        double cauchy) except +

cdef extern from "light.h" nogil:
    void cpp_light_intercept "light_intercept"(
        const double* leaf_locations,
        size_t n_leaves,
        const double* sun_elev,
        const double* sun_azim,
        size_t n_directions,
        double size_leaf,
        double leaf_transparency,
        double* light_per_direction) except +

cdef extern from "momentum.h" nogil:
    void cpp_momentum_solve "momentum_solve"(
        const double* start,
        const double* axis,
        const double* D,
        const double* L,
        size_t n,
        const double* cell_bounds_x,
        size_t nbx,
        const double* cell_bounds_y,
        size_t nby,
        const double* cell_bounds_z,
        size_t nbz,
        double grid_size,
        const double* U_infty,
        double C_D,
        double nu_diff,
        int diffusion_per_line,
        double* U_out,
        double* U_in_grid,
        double* F_D_cell_grid,
        double* U_branch,
        double* F_N_branch,
        double* F_D_branch,
        double* F_vec_branch,
        double cos_theta,
        double sin_theta,
        double* canopy_mean_out) except +
    void cpp_momentum_solve_world "momentum_solve_world"(
        const double* start,
        const double* axis,
        const double* D,
        const double* L,
        size_t n,
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
        int diffusion_per_line,
        double* F_world,
        double* w_world) except +
