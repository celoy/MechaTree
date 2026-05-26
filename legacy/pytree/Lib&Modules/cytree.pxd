#distutils: language=c++

import sys
import cython
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp cimport bool

cdef extern from "tree.h" nogil:
  cdef cppclass Tree:
    Tree(map[string, double] trunk_props) except +

    int getNumberOfBranches()
    map[int,int] getStrahlerDistribution()
    map[int,int] getHortonDistribution()
    map[int,double] meanAggregativePropS(string name)
    map[int,double] meanAggregativePropH(string name)

    void setStrahler()
    int getStrahler(int index)
    void setHorton()
    int getHorton(int index)

    int getLastDescendantIndex(int ancestor_index)
    int getParentIndex(int child_index)
    vector[int] getBrothersIndex(int index)
    vector[int] getChildrenIndex(int index)
    int hasParent(int index)
    int getNumberOfChildren(int parent_index)

    void addBranch(int index,map[string, double] branch_props)
    void removeBranch(int index)

    void addProperty(int index, string name, double value)
    double getProperty(int index, string name)
    void setProperty(int index, string name, double value)
