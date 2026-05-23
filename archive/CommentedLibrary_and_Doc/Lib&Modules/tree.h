/*
This is the declaration of th Tree class
*/

#ifndef TREE_H_
#define TREE_H_

#include <memory>
#include <vector>
#include <map>
#include"branch.h"

using namespace std;

class Tree{

  /*
  This is the formalization of the concept of a tree. The tree is a container of
  pointers to Branch objects. It is also characterized by its hierarchical
  organization ("Strahler_distribution" & "Horton_distribution")
  */

private:
      /*
      All the attributes of the class are private and hence cannot be accesed
      from external programs.

      1) "tree_branches" is a vector stocking pointers to Branch objects. The
      branches can be dinamically allocated.

      2) "Strahler_distribution" is a map containing the information about the
      number of branches per strahler order. The key is the order (an integer)
      and the value is the number of branches for each order (also an integer).

      3) "Horton_distribution" is a map containing the information about the
      number of branches per horton order. The key is the order (an integer)
      and the value is the number of branches for each order (also an integer).

      */

      vector<Branch*> tree_branches;

      map<int,int> Strahler_distribution;
      map<int,int> Horton_distribution;

public:
      /*
      All the methods of the class are public and are designed to manipulate
      the private attributes. This methods can be accessed from an external
      program.

      The methods can be classified in the following way:

      1) Class constructors
      2) Tree structural information
      3) Access branches
      4) Interact with hierarchical classification
      5) Interact with family relations between branches
      6) Interact with the tree structure (add & remove branches)
      7) Interact with branch properties (add, get & set)

      */

//////CONSTRUCTORS/////////////////////////////////////////////////////////////
      Tree();
        /* Nullary constructor */
      Tree(map<string, double> trunk_props);
        /* Non nullary constructor initializes trunk properties. */

//////TREE GENERAL INFORMATION/////////////////////////////////////////////////
      int getNumberOfBranches();
        /* Returns the number of branches of the tree. */
      map<int,int> getStrahlerDistribution();
        /* Returns the number of branches per strahler order. Two contiguous
        branches of same order are considered as the same branch. */
      map<int,int> getHortonDistribution();
        /* Returns the number of branches per Horton order. Two contiguous
        branches of same order are considered as the same branch. */
      map<int,double> meanAggregativePropS(string name);
        /* Returns the mean of an aggregative property (for example the length)
        in function of the branch strahler order. */
      map<int,double> meanAggregativePropH(string name);
        /* Returns the mean of an aggregative property (for example the length)
        in function of branch horton order. */

//////ACCESS BRANCHES//////////////////////////////////////////////////////////
      Branch* getTrunk();
        /* Returns first branch of "tree_branches" vector. */
      Branch* getSummit();
        /* Returns last branch of "tree_branches" vector. */
      Branch* getBranch(int branch_index);
        /* Returns the element of "tre_branches" vector at the position
        "index". */
      int getIndex(Branch* branch);
        /* Returns the position of "branch" (its index) in the "tree_branches"
        vector. */

//////BRANCH HIERARCHY/////////////////////////////////////////////////////////
      void setStrahler();
        /* Sets each branch's Strahler order. */
      int getStrahler(int index);
        /* Returns Strahler order of the branch at the position "index" of
        "tree_branches" vector. */

      void setHorton();
        /* Sets each branch's Horton order. */
      int getHorton(int index);
        /* Returns Horton order of the branch at the position "index" of
        "tree_branches" vector. */

//////BRANCH FAMILY INFORMATION////////////////////////////////////////////////
      int getLastDescendantIndex(int ancestor_index);
        /* Returns the index of the last descendant of the branch located at
        "ancestor_index". */
      int getParentIndex(int child_index);
        /* Returns the index of the parent of the branch located at
        "child_index". */
      vector<int> getBrothersIndex(int branch_index);
        /* Returns a vector containing the indexes of the brothers of the branch
        located at "branch_index". */
      vector<int> getChildrenIndex(int parent_index);
        /* Returns a vector containing the indexes of the children of the branch
        located at "parent_index". */
      int hasParent(int index);
        /* Checks if the branch located at "index" has a parent. Returns 1 if it
         is true, 0 if it is false. */
      int getNumberOfChildren(int parent_index);
        /* Returns the number of children of the branch located at
        "parent_index". */

//////MODIFY TREE_BRANCHES VECTOR//////////////////////////////////////////////
      void addBranch(int parent_index);
        /* Adds a branch  at "parent_index+1" without initializing its
        properties. The new branch is the child of the branch located at
        "parent_index". */
      void addBranch(int parent_index,map<string,double> props);
        /* Adds a branch and initializes its properties as "props" map. */
      void removeBranch(int branch_index);
        /* Removes the branch located at "branch_index" and all its descendants.
         */

//////MODIFY BRANCH PROPERTIES/////////////////////////////////////////////////
      void addProperty(int index,string name, double value);
        /* Creates and initializes property for the branch located at "index".
        The new property's name is "name" and its value is "value". */
      void setProperty(int index,string name, double value);
        /* Actualizes the value of the already existing property "name" of
        the branch located at "index" to the value "value"  */
      double getProperty(int index,string name);
        /* Returns the value of the property "name" of the branch located at
        "index". */

};

#endif
