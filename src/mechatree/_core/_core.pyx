#cython: language_level=3
#cython: wraparound=False
#cython: boundscheck=False
#cython: cdivision=True
#distutils: language=c++

import sys
import cython
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.vector cimport vector

from mechatree._core.cytree cimport Tree

cdef class PyTree:
    """
    PyTree class wraps the C++ Tree class. The only attribute of PyTree is a
    pointer "*ctree" to an instance of the Tree class.
    """
    cdef Tree *c_tree

    def __cinit__(self, dict trunk_props):
        """
        Wrapping of the non-nullary constructor of the Tree class. Trunk's
        properties are given by the dictionnary "trunk_props". The contempt of
        "trunk_props" is copied in a C++ map named "properties". Then the
        constructor of the C++ Tree class is called with "properties" as the
        argument. The attribut "*c_tree" is initialized as a pointer to a Tree
        object with a trunk possessing the properties contained in "trunk_props".
        This function is used to initialize a PyTree object.
        """
        cdef map[string,double] properties
        cdef int i
        cdef tuple tuple_props
        cdef list tuple_props_list=list(trunk_props.items())
        cdef bytes btuple_props
        for i in range(len(tuple_props_list)):
          tuple_props=tuple_props_list[i]
          btuple_props=str2bytes(tuple_props[0])
          properties[btuple_props]=tuple_props[1]
        self.c_tree=new Tree(properties)

    def add_property(self,int index,str name, double value):
        """
        Adds a property named "name" with value "value" to the branch of PyTree
        with index "index". If the branch doesn't exist, nothing is done.
        """
        cdef string bname=str2bytes(name)
        self.c_tree.addProperty(index,bname,value)

    def get_property(self,int index,str name):
        """
        Returns the value of the property "name" of the branch of PyTree with
        index "index". If the branch, or the property don't exist, this method
        returns 0.
        """
        cdef string bname=str2bytes(name)
        return self.c_tree.getProperty(index,bname)

    def set_property(self,int index, str name, double value):
        """
        Sets the value of the already existing property "name" of the branch of
        PyTree with index "index" to "value". If the branch, or the property
        don't exist, this method returns 0.
        """
        cdef string bname=str2bytes(name)
        self.c_tree.setProperty(index,bname,value)

#------------------------------------------------------------------------------

    def get_number_of_branches(self):
        """
        Returns the number of branches of PyTree.
        """
        return self.c_tree.getNumberOfBranches()

    def add_branch(self,int index, dict props):#parent index
        """
        Adds a new branch with the properties contained in the dictionnary
        "props" to PyTree. The branch is added at the index "index+1", where
        "index" is the index of the parent. The creation of family links is
        dealt with in C++ routines. If there is no branch located at "index",
        nothing is done.
        """
        cdef map[string,double] properties
        cdef list tuple_props_list=list(props.items())
        cdef int i
        cdef tuple tuple_props
        cdef bytes btuple_props
        for i in range(len(tuple_props_list)):
          tuple_props=tuple_props_list[i]
          btuple_props=str2bytes(tuple_props[0])
          properties[btuple_props]=tuple_props[1]
        self.c_tree.addBranch(index,properties)

    def remove_branch(self,int index):
        """
        Removes the branch of PyTree with index "index" and all of its
        descendance. The removal of family links is dealt with in C++ routines.
        If there is no branch located at "index", nothing is done.
        """
        self.c_tree.removeBranch(index)

#------------------------------------------------------------------------------

    def get_last_descendant_index(self,int ancestor_index):
        """
        Returns the index of the last descendant of the PyTree branch located at
        "ancestor_index". If there is no branch located at "ancestor_index",
        this method returns 0.
        """
        cdef int ret=self.c_tree.getLastDescendantIndex(ancestor_index)
        return ret

    def get_parent_index(self,int child_index):
        """
        Returns the index of the parent of the PyTree branch located at
        "child_index". If the trunk is located at "child_index", this method
        returns -1. If there is no branch located at "child_index", this method
        returns -2.
        """
        cdef int ret=self.c_tree.getParentIndex(child_index)
        return ret

    def get_brothers_index(self, int index):
        """
        Returns a list containing the indexes of the brothers of the PyTree
        branch located at "index". If there is no branch located at "index",
        this method returns an empty list.
        """
        cdef vector[int] brothers_cpp
        cdef list brothers_py=[]
        cdef int i
        brothers_cpp=self.c_tree.getBrothersIndex(index)
        if len(brothers_cpp)>0:
          for i in range(len(brothers_cpp)):
            brothers_py.append(brothers_cpp[i])
          if len(brothers_py)==1:
            return brothers_py[0]
          else:
            return brothers_py
        else:
           return 0

    def get_children_index(self, int index):
        """
        Returns a list containing the indexes of the children of the PyTree
        branch located at "index". If there is no branch located at "index",
        this method returns an empty list.
        """
        cdef vector[int] children_cpp
        cdef list children_py=[]
        cdef int i
        children_cpp=self.c_tree.getChildrenIndex(index)
        if len(children_cpp)>0:
          for i in range(len(children_cpp)):
            children_py.append(children_cpp[i])

          return children_py
        else:
          return 0

#------------------------------------------------------------------------------

    def has_parent(self,int index):
        """
        Returns 1 if the PyTree branch located at "index" has a parent, and 0 if
        it doesn't. If there is no branch located at "index", this method
        returns -1.
        """
        return self.c_tree.hasParent(index)

    def get_number_of_children(self,int parent_index):
        """
        Returns the number of children of the PyTree branch located at
        "parent_index". If there is no branch located at "parent_index",
        this method returns -1.
        """
        return self.c_tree.getNumberOfChildren(parent_index)

#------------------------------------------------------------------------------

    def set_strahler(self):
        """
        Effectuates the Strahler classification of PyTree.
        """
        self.c_tree.setStrahler()

    def get_strahler(self,int index):
        """
        Returns the Strahler order of the PyTree branch located at "index". If
        there is no branch located at "index", this method returns 0.
        """
        return self.c_tree.getStrahler(index)

    def get_strahler_distribution(self):
        """
        Returns a dictionnary where the keys are the Strahler orders and the
        values the number of branches per order.
        """
        cdef map[int,int] distCpp=self.c_tree.getStrahlerDistribution()
        cdef dict distPy={}
        for it in distCpp:
          distPy[it.first]=it.second
        return distPy

    def set_horton(self):
        """
        Effectuates the Horton classification of PyTree.
        """
        self.c_tree.setHorton()

    def get_horton(self,int index):
        """
        Returns the Horton order of the PyTree branch located at "index". If
        there is no branch located at "index", this method returns 0.
        """
        return self.c_tree.getHorton(index)

    def get_horton_distribution(self):
        """
        Returns a dictionnary where the keys are the Horton orders and the
        values the number of branches per order.
        """
        cdef map[int,int] distCpp=self.c_tree.getHortonDistribution()
        cdef dict distPy={}
        for it in distCpp:
          distPy[it.first]=it.second
        return distPy

    def mean_agg_prop_s(self,str name):
        """
        Returns a dictionnary where the keys are the Strahler orders and the
        values the mean value of a branch property named "name" per order.
        """
        cdef string bname=str2bytes(name)
        cdef map[int,double] meanCpp=self.c_tree.meanAggregativePropS(bname)
        cdef dict meanPy={}
        for it in meanCpp:
          meanPy[it.first]=it.second
        return meanPy

    def mean_agg_prop_h(self,str name):
        """
        Returns a dictionnary where the keys are the Horton orders and the
        values the mean value of a branch property named "name" per order.
        """
        cdef bname=str2bytes(name)
        cdef map[int,double] meanCpp=self.c_tree.meanAggregativePropH(bname)
        cdef dict meanPy={}
        for it in meanCpp:
          meanPy[it.first]=it.second
        return meanPy

#------------------------------------------------------------------------------


cdef str2bytes(str s):
  new_s=s.encode('utf-8')
  return new_s
