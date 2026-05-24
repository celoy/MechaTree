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

cdef extern from "mechanics.h" nogil:
    void cpp_wind_force "wind_force"(
        const Branch& b,
        const array3d& V,
        array3d& force,
        array3d& moment) except +
    void cpp_calculate_stresses "calculate_stresses"(
        Tree& tree, double leaf_drag_S0, double cauchy) except +

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
