/*
This is the implementation of the tree class. Comments in this file aim to
document how the functions work. For a brief description of what the functions
do you can refer to "tree.h".
*/

#include <memory>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

#include "tree.h"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*
1--CONSTRUCTORS
2--TREE GENERAL INFORMATION
3--ACCES BRANCHES
4--BRANCH HIERARCHY
5--BRANCH FAMILY INFORMATION
6--MODIFY TREE_BRANCHES VECTOR
7--MODIFY BRANCH PROPERTIES
*/
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//////1--CONSTRUCTORS//////////////////////////////////////////////////////////

Tree::Tree():tree_branches(vector<Branch*>()){
    Branch* trunk=new Branch();
      /* Creating a pointer to the first branch of the tree */
    tree_branches.insert(tree_branches.begin(),trunk);
      /* Adding "trunk" at the begin of the vector "tree_branches" */
}

/*****************************************************************************/

Tree::Tree(map<string, double> trunk_props):tree_branches(vector<Branch*>()){
    Branch* trunk=new Branch();
      /* Creating a pointer to the first branch of the tree */
    trunk->setProperties(trunk_props);
      /* Initializing "trunk" properties as the map "trunk_props" */
    tree_branches.insert(tree_branches.begin(),trunk);
      /* Adding "trunk" at the begin of the vector "tree_branches" */
}



//////TREE GENERAL INFORMATION/////////////////////////////////////////////////

int Tree::getNumberOfBranches(){
    int number_of_branches=tree_branches.size();
      /* The number of branches is equal to the number of elements of the
      "tree_branches" vector */
    return number_of_branches;
}

/*****************************************************************************/

map <int,int> Tree::getStrahlerDistribution(){
    return Strahler_distribution;
      /* Returns the map "Strahler_distribution" */
}

/*****************************************************************************/

map   <int,int> Tree::getHortonDistribution(){
    return Horton_distribution;
      /* Returns the map "Horton_distribution" */
}

/*****************************************************************************/

map<int,double> Tree::meanAggregativePropS(string name){
    /* "name" is the name of the aggregative property we are interested in */
    map<int,double> mean_prop;
      /* Declaring the return variable */
    map<int,int> dist=getStrahlerDistribution();
      /* "dist" receives the map "Strahler_distribution" */

    for(vector<Branch*>::iterator it=tree_branches.begin();it!=tree_branches.end();it++){
      /* The vector "tree_branches" is traversed element by element with an
      iterator pointing to each element from the beginning to the end */
      double value=(*it)->getProperty(name);
        /* "value" receives the value of the property "name" of the branch of
        "tree_branches" pointed by "it" */
      int order=(*it)->getStrahler();
        /* "order" receives the strahler order of the branch of "tree_branches"
        pointed by "it" */

      mean_prop[order]=mean_prop[order]+value;
        /* The sum of the values of the property "name" is stocked in the map
        "mean_prop". The value is the sum and the key is the order */
    }

    for(map<int,int>::iterator it=dist.begin();it!=dist.end();it++){
      /* The map "dist" is traversed element by element with an iterator
      pointing to each element from the beginning to the end */
      int order=it->first;
        /* The strahler order is the key of the element of the map pointed by
        "it" */
      map<int,double>::iterator ite=mean_prop.find(order);
        /* The iterator "ite" points to the element of the map "mean_prop" with
        key "order" */
      int Nk=it->second;
        /* By definition of the strahler distribution the number of branches of
        order k is given by the value associated with the key k. "it" points to
        an element of the map "dist" which is the strahler distribution. */
      double lk=ite->second;
        /* "ite" points to the element of the map "mean_prop" with key equal to
        the order k ("it->second"). "lk" is the sum of the lengths of all the
        branches of order k */
      ite->second=lk/Nk;
        /* The mean length of the branches of order k is given dividing lk by Nk
        (the number of branches of order k). The value of the element pointed
        by "ite" is changed by the mean value of the property */

    }

    return mean_prop;
      /* At the end of the procedure "mean_props" contains the mean of an
      aggregative property in function of branches orders, the function returns
      this map */
}

/*****************************************************************************/

map<int,double> Tree::meanAggregativePropH(string name){
    /* "name" is the name of the aggregative property we are interested in */
    map<int,double> mean_prop;
      /* Declaring the return variable */
    map<int,int> dist=getHortonDistribution();
      /* "dist" receives the map "Horton_distribution" */

    for(vector<Branch*>::iterator it=tree_branches.begin();it!=tree_branches.end();it++){
      /* The vector "tree_branches" is traversed element by element with an
      iterator pointing to each element from the beginning to the end */
      double value=(*it)->getProperty(name);
        /* "value" receives the value of the property "name" of the branch of
        "tree_branches" pointed by "it" */
      int order=(*it)->getHorton();
        /* "order" receives the horton order of the branch of "tree_branches"
        pointed by "it" */

      mean_prop[order]=mean_prop[order]+value;
        /* The sum of the values of the property "name" is stocked in the map
        "mean_prop". The value is the sum and the key is the order */
    }

    for(map<int,int>::iterator it=dist.begin();it!=dist.end();it++){
      /* The map "dist" is traversed element by element with an iterator
      pointing to each element from the beginning to the end */
      int order=it->first;
        /* The strahler order is the key of the element of the map pointed by
        "it" */
      map<int,double>::iterator ite=mean_prop.find(order);
        /* The iterator "ite" points to the element of the map "mean_prop" with
        key "order" */
      int Nk=it->second;
        /* By definition of the strahler distribution the number of branches of
        order k is given by the value associated with the key k. "it" points to
        an element of the map "dist" which is the strahler distribution. */
      double lk=ite->second;
        /* "ite" points to the element of the map "mean_prop" with key equal to
        the order k ("it->second"). "lk" is the sum of the lengths of all the
        branches of order k */
      ite->second=lk/Nk;
        /* The mean length of the branches of order k is given dividing lk by Nk
        (the number of branches of order k). The value of the element pointed
        by "ite" is changed by the mean value of the property */

    }

    return mean_prop;
      /* At the end of the procedure "mean_props" contains the mean of an
      aggregative property in function of branches orders, the function returns
      this map */
}



//////2--ACCESS BRANCHES///////////////////////////////////////////////////////

Branch* Tree::getTrunk(){
    return tree_branches[0];
      /* Returning the first branch (index 0) of the "tree_branches" vector */
}

/*****************************************************************************/

Branch* Tree::getSummit(){
    int N=tree_branches.size();
      /* N-1 is the index of the last branch of "tree_branches" vector */
    return tree_branches[N-1];
}

/*****************************************************************************/

Branch* Tree::getBranch(int index){
  if(index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches in the tree */
    if(index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      unsigned int u_index=index;
      return tree_branches.at(u_index);
        /* The function returns the element of "tree_branches" located at
        "index" after finding it */
    }
    else{
      /* If the index of the branch to find is bigger than the number of
      branches in the tree, display an error message and return a pointer to a
      null branch */
      cerr << "Tree::getBranch Error: trying to acces a branch at a non-existing index. Returning a default branch.\n";
      Branch* branch2return;
      branch2return=new Branch();
      return branch2return;
    }
  }
  else{
    /* If the index of the branch to find is negative display an error message
    and return a pointer to a null branch */
    cerr << "Tree::getBranch Error: the branch index can't be negative. Returning a default branch.\n";
    Branch* branch2return;
    branch2return=new Branch();
    return branch2return;
  }
}

/*****************************************************************************/

int Tree::getIndex(Branch* branch){
    int index=0;
      /* Declaring the return variable */
    for(vector<Branch*>::iterator it=tree_branches.begin();it!=tree_branches.end();it++){
      /* The vector "tree_branches" is traversed element by element with an
      iterator pointing to each element from the beginning to the end */
      if((*it)==branch){
        /* If the element of "tree_branches" pointed by "it" is the branch we
        are looking for ... */
        return index;
          /* ... the function returns index */
      }
      /* If the "it" wasn't pointing to the desired branch, increment index
      and repeat the procedure */
      index++;
    }
    /* If the function hasn't already returned something it means that "branch"
    isn't at "tree_branches". Then display an error message and return -1 */
    cerr<< "Tree::getIndex Error: trying to get the index of a branch we can't find. Returning 0.\n";
    return -1;
}


//////3--BRANCH HIERARCHY//////////////////////////////////////////////////////

void Tree::setStrahler(){ /*/!\ CASE NBER OF CHILDREN==1 NOT PRESENT. CAN BE GENERALIZED TO MORE CHILDREN THAN 2 */

    for(vector<Branch*>::reverse_iterator rit=tree_branches.rbegin();rit!=tree_branches.rend();rit++){
      /* The vector "tree_branches" is traversed backwards element by element
      with an iterator pointing to each element from the beginning to the end */
      if((*rit)->hasChildren()==false){
        /* If the considered branch hasn't any children...*/
        (*rit)->setStrahler(1);
          /*... then it is a leaf and its order is equal to 1 */
        Strahler_distribution[1]=Strahler_distribution[1]+1;
          /* The number of branches of order equal to 1 is incremented */
      }

      else{
        /* If the considered branch has at least one children... */

        vector<Branch*> children=(*rit)->getChildren();
          /* Retrieving the children of the considered branch */
        if(children.size()==2){
          /* If the considered branch has 2 children */
          int strahler1=children[0]->getStrahler();
            /* Retrieving the strahler order of the first one */
          int strahler2=children[1]->getStrahler();
            /* Retrieving the strahler order of the second one */
          if(strahler1==strahler2){
            /*If both children have the same order... */
            (*rit)->setStrahler(strahler1+1);
              /* ... then the strahler order of the considered branch is the
              strahler order of the children plus one */

            Strahler_distribution[strahler1+1]=Strahler_distribution[strahler1+1]+1;
              /* The number of branches of order "strahler+1" is incremented */
          }

          else{
            /* If children have different orders... */
            int strahler=max(strahler1,strahler2);
              /* Retrieving the maximum of the children orders */
            (*rit)->setStrahler(strahler);
              /* Setting the considered branch strahler's order to the children
              greatest order */

            /* There is no need of incrementing number of branches of order
            strahler since two contiguous branches of the same order are
            considered as the same branch. */
          }
        }
        else{
          /* If the branch has more than 2 children display an error message
          and do nothing */
          cerr<< "Tree::setStrahler Error: trying to set strahler order but tree has more than 2 children per branching.\n ";
        }
      }
    }

}

/*****************************************************************************/

int Tree::getStrahler(int index){
  if(index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches in the tree */
    if(index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      unsigned int i=index;
      return tree_branches.at(i)->getStrahler();
        /* The function returns the strahler order of the  element of
        "tree_branches" located at "index" after finding it */
    }
    else{
      /* If the index of the branch to find is bigger than the number of
      branches in the tree, display an error message and return 0 */
      cerr<<"Tree::getStrahler Error: trying to acces a branch at a non-existing index. Returning 0.\n ";
      return 0;
    }
  }
  else{
    /* If the index of the branch to find is negative display an error message
    and return 0 */
    cerr << "Tree::getStrahler Error: the branch index can't be negative. Returning 0.\n";
    return 0;
  }
}


/*****************************************************************************/

void Tree::setHorton(){

    if (Strahler_distribution.empty()==true){
      /* If the Strahler classification wasn't done previously ... */
      setStrahler();
        /* Effectuate the Strahler classification since it is the first step of
        the Horton classification */
    }

    Branch* trunk=getTrunk();
      /* "trunk" is the first branch of "tree_branches" */
    int horton_trunk=trunk->getStrahler();
      /* Retrieving the Strahler order of "trunk" */
    trunk->setHorton(horton_trunk);
      /* For the trunk, Strahler and Horton orders are the same */

    Horton_distribution[horton_trunk]=1;
      /* Iinitializing the Horton distribution for the trunk's order */

    for(vector<Branch*>::iterator it=tree_branches.begin();it!=tree_branches.end();it++){
      /* The vector "tree_branches" is traversed element by element with an
      iterator pointing to each element from the beginning to the end */
      vector<Branch*> children=(*it)->getChildren();
        /* Retrieving the children of the considered branch */

      if(children.size()>0){
        /* If the considered branch has any children... */

        if(children.size()==2){
          /* If the considered branch has 2 children */

          int strahler1=children[0]->getStrahler();
            /* Retrieving the strahler order of the first child of the
            "tree_branches" element pointed by "it" */
          int strahler2=children[1]->getStrahler();
            /* Retrieving the strahler order of the first child of the
            "tree_branches" element pointed by "it" */
          int hortondad=(*it)->getHorton();
            /* Retrieving the horton order of the considered branch (the
            "tree_branches" element pointed by "it") */

          if(strahler1>=strahler2){
            /* If the strahler order of the first child is bigger or equal
            than the one of the second */

            children[0]->setHorton(hortondad);
              /* The first child gets its dad's horton order */
            children[1]->setHorton(strahler2);
              /* The second child keeps its strahler order as its horton order */

            Horton_distribution[strahler2]=Horton_distribution[strahler2]+1;
              /* Incrementing the number of branches of horton order strahler2 */

            /* There is no need of actualizing Horton_distribution for order
            hortondad since contiguous branches of same order are considered
            as the same branch. */
          }

          else if(strahler1<strahler2){
            /* If the strahler order of the second child is bigger than the one
            of the first */
            children[1]->setHorton(hortondad);
              /* The second child gets its dad's horton order */
            children[0]->setHorton(strahler1);
              /* The first child keeps its strahler as its horton order */

            Horton_distribution[strahler1]=Horton_distribution[strahler1]+1;
              /* Incrementing the number of branches of horton-order strahler1 */
          }
        }

        /* When the two children have the same order we choose arbitrarily
        whom receives the same horton order as its dad. In this case it is
        the first one but it could be the second or it could be chosen randomly
        */

        else{
          /* If the considered branch hasn't exactly two children display an
          error message and do nothing */
          cerr<< "Tree::setHorton Error: trying to set horton order but tree has more than 2 children per branching.\n ";
        }
      }
    }
}

/*****************************************************************************/

int Tree::getHorton(int index){
  if(index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      unsigned int i=index;
      return tree_branches.at(i)->getHorton();
        /* The function returns the strahler order of the  element of
        "tree_branches" located at "index" after finding it */
    }
    else{
      /* If the index of the branch to find is bigger than the number of
      branches in the tree, display an error message and return 0 */
      cerr<< "Tree::getHorton Error: trying to acces a branch at a non-existing index. Returning 0.\n ";
      return 0;
    }
  }
  else{
    /* If the index of the branch to find is negative display an error message
    and return 0 */
    cerr << "Tree::getHorton Error: the branch index can't be negative. Returning 0.\n";
    return 0;
  }
}



//////4--BRANCH FAMILY INFORMATION/////////////////////////////////////////////

int Tree::getLastDescendantIndex(int ancestor_index){
  if(ancestor_index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(ancestor_index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      Branch* bancestor=tree_branches.at(ancestor_index);
        /* "bancestor" is the branch located at "ancestor_index" */
      if(ancestor_index==0){
        /* If "bancestor" is the trunk... */
        int last_descendant_index=N-1;
        return last_descendant_index;
          /* ... then its last descendant is the last element of "tree_branches"
          , hence its index is "N-1" */
      }

      else if(bancestor->hasChildren()==false){
        /* If "bancestor" is a leaf... */
        return ancestor_index;
          /* ... then its last ancestor is itself, its index is returned */
      }

      else{
        /* If the branch isn't the trunk or a leaf... */
        vector<int> brothers_indexes=Tree::getBrothersIndex(ancestor_index);
          /* Retrieving the brothers of "bancestor" */
        for(vector<int>::iterator it=brothers_indexes.begin();it!=brothers_indexes.end();it++){
          /* Traversing the brothers of the "bancestor" */
          int brother_index=(*it);
            /* "it" points to the index of "bancestor" brother's */
          if(brother_index>ancestor_index){
            /* If the brother index is bigger than "bancestor's" index... */
            int last_descendant_index=brother_index-1;
            return last_descendant_index;
              /* ... then by construction of the indexation the
              "last_descendant_index" is "brother_index-1" */
          }
        }
        /* If the function hasn't returned anything yet it is because the
        index of the brother was smaller than the one of "bancestor", making it
        impossible to determine the "last_descendant_index" with this aproach.
        Nevertheless, the index of the last descendant of "bancestor" is the
        as the one of its father! Hence we can do the same test for "bancestor"
        fathers. Eventually, if we get until the trunk it means that the index
        of the last descendant of bancestor is "N-1" */
        int parent_index=Tree::getParentIndex(ancestor_index);
          /* Retrieving the index of "bancestor's" parent */
        return Tree::getLastDescendantIndex(parent_index);
          /* Returning the last descendant index of "bancestor" parent */

      }
    }
    else{
      /* If "bancestor" index is bigger than the number of branches in
      "tree_branches" display an error message and return -1 */
      cerr<< "Tree::getLastDescendantIndex Error: trying to find the last descendant of a branch at a non-existing index. Returning -1.\n";
      return -1;
    }
  }
  else{
    /* If "bancestor" index is negative then display an error message and return
    -1 */
    cerr << "Tree::getLastDescendantIndex Error: the branch index can't be negative. Returning -1.\n";
    return -1;
  }
}

/*****************************************************************************/

int Tree::getParentIndex(int child_index){
  if(child_index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(child_index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      if(child_index==0){
        /* If the considered branch is the trunk then it hasn't a parent. the
        index -1 is returned */
        return -1;
      }
      unsigned int u_child_index=child_index;
      Branch* parent=tree_branches.at(u_child_index)->getParent();
        /* Retrieving the parent of the considered branch */
      int parent_index=Tree::getIndex(parent);
        /* Retrieving and returning the index of the parent */
      return parent_index;
    }
    else{
      /* If the child index was bigger than the number of branches in the tree
      display an error message and return -2 */
      cerr<< "Tree::getParentIndex Error: trying to get the index of a branch at a non-existing index. Returning -2.\n ";
      return -2;
    }
  }
  else{
    /* If the child index was negative display an error message and return -2 */
    cerr << "Tree::getParentIndex Error: the branch index can't be negative. Returning -2.\n";
    return -2;
  }
}

/*****************************************************************************/

vector<int> Tree::getBrothersIndex(int branch_index){
  vector<int> brothers_index;
    /* Declaring the return variable */

  if(branch_index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(branch_index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      Branch* this_branch=Tree::getBranch(branch_index);
        /* Retrieving the branch at "branch_index" */
      vector<Branch*> brothers=this_branch->getBrothers();
        /* Retrieving branch's brothers */

      for(vector<Branch*>::iterator it=brothers.begin();it!=brothers.end();it++){
          /* Traversing a vector containing the brothers of the considered
          branch */
          int this_index=Tree::getIndex(*it);
            /* Retrieving the index of a brother */
          brothers_index.push_back(this_index);
            /* Stocking the index in the return variable */
        }
      return brothers_index;
    }

    else{
      /* If the considered branch's index was bigger than the number of branches
      in "tree_branches" display an error message and return an empty vector */
      cerr<< "Tree::getBrothersIndex Error: trying to get the brothers of a branch at a non-existing index. Returning an empty vector.\n";
      return brothers_index;
    }
  }

  else{
    /* If the considered branch's index was negative display an error message
    and return an empty vector */
    cerr << "Tree::getBrothersIndex Error: the branch index can't be negative. Returning an empty vector.\n";
    return brothers_index;
  }
}

/*****************************************************************************/

vector<int> Tree::getChildrenIndex(int parent_index){
  vector<int> children_indexes;
  /* Declaring the return variable */
  if(parent_index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if (parent_index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      Branch* this_branch=Tree::getBranch(parent_index);
        /* Retrieving the branch located at "parent_index" */
      vector<Branch*> children=this_branch->getChildren();
        /* Retrieving the children of the branch located at "parent_index" */

      for(vector<Branch*>::iterator it=children.begin(); it!=children.end();it++){
        /* Traversing element by element the vector containing the children
        of the considered branch */
        int this_index=Tree::getIndex(*it);
          /* Retrieving the indexes of the children */
        children_indexes.push_back(this_index);
          /* Stocking the children indexes in "children_indexes" vector */
      }
      return children_indexes;
    }

    else{
      /* If the considered branch's index was bigger than the number of branches
      in "tree_branches" display an error message and return an empty vector */
      cerr<< "Tree::getChildrenIndex Error: trying to get the children of a branch at a non-existing index. Returning an empty vector.\n";
      return children_indexes; //returning empty vector

    }
  }
  else{
    /* If the considered branch's index was negative display an error message
    and return an empty vector */
    cerr << "Tree::getChildrenIndex Error: the branch index can't be negative. Returning an empty vector.\n";
    return children_indexes;
  }
}

/*****************************************************************************/

int Tree::hasParent(int index){
  if(index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      unsigned int i=index;
      Branch* branch2test=tree_branches.at(i);
        /* Retrieving the branch located at "index" */
      if (branch2test->hasParent()==true){
        /* If the considered branch has parents then return 1 */
        return 1;
      }
      else{
        /* If the considered branch hasn't parents then return 0 */
        return 0;
      }
    }
    else{
      /* If the considered branch's index was bigger than the number of branches
      in "tree_branches" display an error message and return an empty vector */
      cerr<< "Tree::hasParent Error: trying to check if a branch at a non-existing index has a parent. Returning -1.\n";
      return -1;
    }
  }
  else{
    /* If the considered branch's index was negative display an error message
    and return an empty vector */
    cerr << "Tree:hasParent Error: branch index can't be negative. Returning -1.\n";
    return -1;
  }
}

/*****************************************************************************/

int Tree::getNumberOfChildren(int parent_index){
  if(parent_index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(parent_index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      unsigned int u_parent_index=parent_index;
      Branch* parent=tree_branches.at(u_parent_index);
        /* Retrieving the branch located at "parent_index" */
      vector<Branch*> children=parent->getChildren();
        /* Retrieving the children of the considered branch */
      int nb_child=children.size();
      return nb_child;
        /* Retrieving and returning the number of children of the considered
        branch */
    }
    else{
      /* If the considered branch's index was bigger than the number of branches
      in "tree_branches" display an error message and return an empty vector */
      cerr<< "Tree::getNumberOfChildren Error: trying to know how many children has a branch at a non-existing index. Returning -1.\n";
      return -1;
    }
  }
  else{
    /* If the considered branch's index was negative display an error message
    and return an empty vector */
    cerr << "Tree::getNumberOfChildren Error: the branch index can't be negative. Returning -1.\n";
    return -1;

  }
}


//////5--MODIFY TREE///////////////////////////////////////////////////////////

void Tree::addBranch(int parent_index){
  if(parent_index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(parent_index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      unsigned int u_parent_index=parent_index;
      Branch* mom=tree_branches.at(u_parent_index);
        /* Retrieving the futur parent "mom" */
      Branch* branch2insert=new Branch();
        /* Creation of the branch to insert */
      mom->addChild(branch2insert);
        /* Setting the family relationship between "mom" and the new branch */
      vector<Branch*>::iterator begin=tree_branches.begin();
        /* Creating iterators to manipulate vector elements */
      unsigned int inserting_position=u_parent_index+1;
        /* Initializing "inserting_position" */
      tree_branches.insert(begin+inserting_position,branch2insert);
        /* Inserting the new branch */
    }
    else{
      /* If the parent's index was bigger than the number of branches in
      "tree_branches" display an error message and do nothing */
      cerr<< "Tree::addBranch Error: trying to add a branch but the futur parent index doesn't exist.\n";
    }
  }
  else{
    /* If the parent's index was negative display an error message and return
    an empty vector */
    cerr << "Tree::addBranch Error: the branch index can't be negative.\n";
  }
}

/*****************************************************************************/

void Tree::addBranch(int parent_index,map<string,double> props){

  if(parent_index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(parent_index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      unsigned int u_parent_index=parent_index;
      Branch* mom=tree_branches.at(u_parent_index);
        /* Retrieving the futur parent "mom" */
      Branch* branch2insert=new Branch();
        /* Creation of the branch to insert */
      branch2insert->setProperties(props);
        /* Setting new branch properties as the map "props" */
      mom->addChild(branch2insert);
        /* Setting the family relationship between "mom" and the new branch */
      vector<Branch*>::iterator begin=tree_branches.begin();
        /* Creating iterators to manipulate vector elements */
      unsigned int inserting_position=u_parent_index+1;
        /* Initializing "inserting_position" */
      tree_branches.insert(begin+inserting_position,branch2insert);
        /* Inserting the new branch */
    }
    else{
      /* If the parent's index was bigger than the number of branches in
      "tree_branches" display an error message and do nothing */
      cerr<< "Tree::addBranch Error: trying to add a branch but the futur parent index doesn't exist.\n";
    }
  }
  else{
    /* If the parent's index was negative display an error message and do
    nothing */
    cerr << "Tree::addBranch Error: the branch index can't be negative.\n";
  }
}

/*****************************************************************************/

void Tree::removeBranch(int branch_index){
  if(branch_index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(branch_index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      unsigned int u_branch_index=branch_index;
      Branch* branch2remove=tree_branches.at(u_branch_index);
        /* Retrieving the branch to remove */
      branch2remove->removeParent();
        /* Removing parenthood link, if "branch2remove" is the trunk nothing
        will be done */
      if(branch_index!=0){
        /* Removing childhood link between parent and branch located at
        "branch_index" */
        Branch* branch2removeParent=branch2remove->getParent();
          /* Retrieving the parent of the branch to remove! */
        branch2removeParent->removeChild(branch2remove);
          /* Removing childhood link*/
      }

      vector<Branch*>::iterator begin=tree_branches.begin();
        /* Creating iterators to manipulate vector elements */
      int last_descendant_index=Tree::getLastDescendantIndex(branch_index);
        /* Retrieving the index of the considered branch's last descendant to
        remove the whole subtree */
      branch2remove->removeDescendants();
        /* First all the family links of the subtree are removed, then all the
        subtree is removed from "tree_branches" */

      if(last_descendant_index==N-1){
        tree_branches.erase(tree_branches.begin()+u_branch_index,tree_branches.end());
      }

      else{
        unsigned int u_last_descendant_index=last_descendant_index;
        tree_branches.erase(tree_branches.begin()+u_branch_index,tree_branches.begin()+u_last_descendant_index+1);
        /* Removing the considered branch and its descendants from the
        "tree_branches" vector*/
      }
    }
    else{
      /* If the considered branch's index was bigger than the number of branches
       in "tree_branches" display an error message and do nothing */
      cerr << "Tree::removeBranch Error: trying to remove a branch at a non-existing index.\n";
    }
  }
  else{
    /* If the considered branch's index was negative display an error message
    and do nothing */
    cerr << "Tree::removeBranch Error: the branch index can't be negative.\n";
  }
}



//////6--BRANCH PROPERTIES/////////////////////////////////////////////////////

void Tree::addProperty(int index,string name,double value){
  if(index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      unsigned int i=index;
      tree_branches.at(i)->addProperty(name,value);
        /* Adding the property named "name" with value "value" to the branch
        located at "i". The function "getProperty" of the branch class has its
        own error message if the property named "name" doesn't exist  */
    }
    else{
      /* If the considered branch's index was bigger than the number of branches
       in "tree_branches" display an error message and do nothing */
      cerr << "Tree::addProperty Error: trying to add a property to a branch at a non-existing index.\n";
    }
  }
  else{
    /* If the considered branch's index was negative display an error message
    and do nothing */
    cerr << "Tree::addProperty Error: the branch index can't be negative.\n";
  }
}

/*****************************************************************************/

double Tree::getProperty(int index,string name){
  if(index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      unsigned int i=index;
      return tree_branches.at(i)->getProperty(name);
        /* Returning the value of the property named "name" of the branch
        located at "i". The function "getProperty" of the branch class has its
        own error message if the property named "name" doesn't exist */
    }
    else{
      /* If the considered branch's index was bigger than the number of branches
       in "tree_branches" display an error message and return 0 */
      cerr << "Tree::getProperty Error: trying to get property "<< name <<" from a branch at the non-existing index "<< index<<". Returning 0.\n";
      return 0;
    }
  }
  else{
    /* If the considered branch's index was negative display an error message
    and return 0 */
    cerr << "Tree::getProperty Error: the branch index can't be negative. Returning 0.\n";
    return 0;
  }
}

/*****************************************************************************/

void Tree::setProperty(int index,string name,double value){
  if(index>=0){
    /* Negative indexes aren't possible, since the index is the position of the
    branch in "tree_branches" vector */
    int N=tree_branches.size();
      /* N is the number of branches */
    if(index<N){
      /* A branch with an index bigger than the number of branches in the tree
      isn't a part of the tree */
      unsigned int i=index;
      tree_branches.at(i)->setProperty(name,value);
        /* Setting the value of the property named "name" of the branch located
        at "index" to "value"*/
    }
    else{
      /* If the considered branch's index was bigger than the number of branches
       in "tree_branches" display an error message and do nothing */
      cerr << "Tree::setProperty Error: trying to set a property from a branch at a non-existing index. Returning 0.\n";
    }
  }
  else{
    /* If the considered branch's index was negative display an error message
    and do nothing */
    cerr << "Tree::setProperty Error: the branch index can't be negative. \n";
  }
}
