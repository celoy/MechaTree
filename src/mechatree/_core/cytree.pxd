#distutils: language=c++

from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

cdef extern from "tree.h" nogil:
    cdef cppclass Tree:
        Tree(unordered_map[string, double] trunk_props) except +

        int getNumberOfBranches()
        const map[int, int]& getStrahlerDistribution()
        const map[int, int]& getHortonDistribution()
        map[int, double] meanAggregativePropS(string name) except +
        map[int, double] meanAggregativePropH(string name) except +

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
        void removeBranch(int index) except +

        void addProperty(int index, string name, double value) except +
        double getProperty(int index, string name) except +
        void setProperty(int index, string name, double value) except +
