/*
This is the implementation of the branch class. Comments in this file aim to
document how the functions work. For a brief description of what the functions
do you can refer to "branch.h".
*/

#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <iostream>
#include <sstream>
#include <map>
#include <stdio.h>

using namespace std;

#include "branch.h"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*

Here is the order of implementation of the functions in this file:

1--CONSTRUCTOR
    -Branch()
2--PROPERTIES
    -void addProperty(string name, double value )
    -void setProperty(string name, double value)
    -void setProperties(map<string,double> props)
    -double getProperty(string name)
3--PARENT
    -bool hasParent()
    -bool hasParent(Branch*test)
    -void addParent(Branch*p)
    -void removeParent(Branch*p)
    -void removeParent()
    -Branch* getParent()
4--CHILDREN
    -bool hasChildren()
    -bool hasChild(Branch*test)
    -void addChild(Branch*ch)
    -void removeChild(Branch*ch)
    -void removeChildren()
    -void removeDescendants()
    -vector<Branch*> getChildren()
    -vector<Branch*> getBrothers()
5--BRANCH HIERARCHY
    -void setStrahler(int strahler_index)
    -int getStrahler()
    -void setHorton(int horton_index)
    -int getHorton()

*/
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//////1--CONSTRUCTOR///////////////////////////////////////////////////////////

Branch::Branch():properties(map<string,double>()),parent(vector<Branch*>()),
                children(vector<Branch*>())
        {
          /*Nullary constructor all the container-attributes are initialized as
          empty containers. Strahler and Horton??? */
        }



//////2--PROPERTIES////////////////////////////////////////////////////////////

void Branch::addProperty(string name, double value){
  map<string,double>::iterator it=properties.find(name);
    /* "it" points to the element of the map "properties" with "name" as key, if
    there is no element in the map verifying this condition, "it" points to the
    "end" of the map: properties.end() */
  if(it==properties.end()){
    /* If there is no property named "name"... */
    properties[name]=value;
      /* ...add a new property. New element in the map "properties": the key
      "name" is associated with the value "value" */
  }
  else{
    /* If the property already exists do nothing and display a message error */
    cerr << "Branch::addProperty Error: trying to add an already existing property. Use set instead.\n";
  }
}


/*****************************************************************************/

void Branch::setProperty(string name, double value){
  map<string, double>::iterator it=properties.find(name);
    /* "it" points to the element of the map "properties" with "name" as key, if
    there is no element in the map verifying this condition, "it" points to the
    "end" of the map: properties.end() */
  if(it!=properties.end()){
    /* If the property exists... */
    it->second=value;
      /* The value of the map element pointed by it is changed to "value". Here
      "it->second" is the value of the map element pointed by "it", "it->first"
      would be its key. */
  }
  else{
    /* If the property doesn't exist display an error message and do nothing. */
    cout << "Branch::setProperty Error: trying to set the value of a property that I'm unable to find. Please check the spelling of the property's name. \n";
  }
}

/*****************************************************************************/

void Branch::setProperties(map<string,double> props){
  properties=props;
  /* The map "properties" receives the map "props" */
}

/*****************************************************************************/

double Branch::getProperty(string name){
  map<string,double>::iterator it=properties.find(name);
    /* "it" points to the element of the map "properties" with "name" as key, if
    there is no element in the map verifying this condition, "it" points to the
    "end" of the map: properties.end() */
  if(it!=properties.end()){
    /* If the property exists... */
    double value=it->second;
      /* "value" gets the value associated with the key "name". Here "it->second"
      is the value of the map element pointed by "it", "it->first" would be its
      key. */
    return value;
      /* The value of the property is returned */
  }
  else{
    /* If the property doesn't exist display an error message and return 0 */
    cerr << "Branch::getProperty Error: trying to get the value of a property that I'm unable to find. Please check the spelling of the property's name. Returning 0.\n";
    return 0;
  }
}


//////3--PARENT////////////////////////////////////////////////////////////////

bool Branch::hasParent(){
  if(parent.size()==0){
    /* If there are no elements in the vector "parent", then the branch has no
    parent, the function returns false */
    return false;
  }
	else{
    /* If the vector "parent" is not empty, then the branch has a parent, the
    function returns true */
    return true;
  }
}

/*****************************************************************************/

bool Branch::hasParent(Branch* test){
  bool ret = false;
    /* The return value is set to false by default */
  for(vector<Branch*>::iterator it=parent.begin();it!=parent.end();it++){
    /* The vector parent is traversed element by element with an iterator
    pointing to each element from the beginning to the end */
    if((*it)==test){
      /* If it points to "test", then "test" is an element of the vector
      "parent" and hence the parent of the branch is indeed "test" */
      ret = true;
        /* The return value is changed to true */
      break;
        /* Getting out of the for loop */
    }
  }
  return ret;
    /* The function returns the return value */
}

/*****************************************************************************/

void Branch::addParent(Branch* mom){
  if(hasParent(mom)==false){
    /* If "mom" isn't already the parent of the branch... */
    if(parent.size()==0){
      /* If there are no elements in the vector "parent", then the branch has no
      parent and we can proceed with the adding. Multiple parents aren't allowed
      in a tree-like network */
		  parent.push_back(mom);
        /* "mom" is added at the end of the "parent" vector */
	  }
    else{
      /* If the branch has already a parent display an error message and do
      nothing */
      cout << "Branch::addParent Error: trying to add a parent to a branch that already has one.\n";
    }
  }
  else{
    /* If the "mom" is already the parent of the branch display an error message
    and do nothing */
    cout << "Branch::addParent Error: trying to add a parent twice.\n";
  }
}

/*****************************************************************************/

void Branch::removeParent(Branch* mom){
  if(hasParent(mom)==true){
    /* If "mom" is the parent of the branch we can proceed with the removal */
		for(vector<Branch*>::iterator it=parent.begin();it!=parent.end();it++){
      /* The vector parent is traversed element by element with an iterator
      pointing to each element from the beginning to the end */
			if((*it)==mom){
        /* If it points to "test", then "test" is an element of the vector
        "parent" and hence the parent of the branch is indeed "mom" */
				parent.erase(it);
          /* Erase the map element pointed by "it" from "parent" vector */
				break;
          /* Getting out of the for loop */
			}
		}
	}
  else{
    /* If "mom" isn't the parent of the branch display an error message and do
    nothing */
		cout << "Branch::removeParent Error: trying to remove a parent that doesn't exist.\n";
  }
}

/*****************************************************************************/

void Branch::removeParent(){
  if(this->hasParent()==true){
    /* If this branch has a parent we can proceed with the removal */
    vector<Branch*>::iterator it=parent.begin();
      /* "it" points to the first element of the "parent" vector. Since there
      is only one element in "parent", "it" points to the parent of the branch
      */
    parent.erase(it);
      /* Erase the vector "parent" element pointed by it, hence the parent of
      the branch */
  }
  else{
    /* If the branch hasn't any parent do nothing */
    cout << "Removing the trunk.\n"; // WHY THIS MESSAGE????
  }
}

/*****************************************************************************/

Branch* Branch::getParent(){
  return parent[0];
    /* Return the first element of the "parent" vector: the parent of this
    branch */
}



//////4--CHILDREN//////////////////////////////////////////////////////////////

bool Branch::hasChildren(){
  if(children.size()==0){
    /* If there are no elements in the "children" vector (the branch has no
    children) return false */
  	return false;
  }
	else{
    /* If "children" vector isn't empty (the branch has children) return true */
    return true;
  }
}

/*****************************************************************************/

bool Branch::hasChild(Branch* test){
  bool ret = false;
  /* The return variable is set to false by default */
  for(vector<Branch*>::iterator it=children.begin();it!=children.end();it++){
    /* The vector "children" is traversed element by element with an iterator
    pointing to each element from the beginning to the end */
    if((*it)==test){
      /* If the "parent" vector element pointed by "it" is "test", then "test"
      is a child of the branch */
      ret = true;
        /* The return variable is changed to true */
      break;
        /* Getting out of the for loop */
    }
  }
  return ret;
    /* The function returns the return variable */
}

/*****************************************************************************/

void Branch::addChild(Branch* ch){
  if(hasChild(ch)==false){
    /* If "ch" isn't already a child of this branch we can proceed with the
    adding */
		children.push_back(ch);
      /* "ch" is added at the end of "children" vector */
		ch->addParent(this);
      /* The parenthood link is established with the function "addParent" */
	}
  else{
    /* If "ch" was already a child of this branch display an error message and
    do nothing */
		cout << "Branch::addChild Error: trying to add a child but the child has been already added.\n";
  }
}

/*****************************************************************************/

void Branch::removeChild(Branch* ch){
  if(hasChild(ch)==true){
    /* If "ch" is a child oh this branch we can proceed with the removal */
		for(vector<Branch*>::iterator it=children.begin();it!=children.end();it++){
      /* The vector "children" is traversed element by element with an iterator
      pointing to each element from the beginning to the end */
			if((*it)==ch){
        /* If the "parent" vector element pointed by "it" is "ch", then "ch"
        is a child of the branch */
				children.erase(it);
          /* Removing the child: the "children" vector element pointed by "it"
          is erased from the vector */
				break;
          /* Getting out of the for loop */
			}
		}
	}
  else{
    /* If "ch" isn't a child of this branch display an error message and do
    nothing */
	  cout << "Branch::removeChild Error: trying to remove a child that doesn't exist.\n";
  }
}

/*****************************************************************************/

void Branch::removeChildren(){
  if(this->hasChildren()==true){
    /* If this branch has children we can proceed with the removal */
    children.erase(children.begin(),children.end());
      /* All the elements of the "children" vector are erased */
  }
  else{
    /* If this branch hasn't any children display an error message and do
    nothing */
    cout << "Branch::removeChildren Error: trying to remove all the children of a branch but tha branch hasn't any.\n";
  }
}

/*****************************************************************************/

void Branch::removeDescendants(){
  /* This is a recursive function traversing all the descendance of the branch
  for which the function is called. At each step the children of the branches
  are removed and then the parenting links are removed. */
  if (this!=0){
    /* If this branch is not a NULL pointer (it exists!) we can proceed */
    if(this->hasChildren()==true){
      /* If this branch has children we can proceed */
      for(vector<Branch*>::iterator it=children.begin();it!=children.end();it++){
        /* The vector "children" is traversed element by element with an
        iterator pointing to each element from the beginning to the end */
        (*it)->removeDescendants();
          /* The function "removeDescendants" is also called for the children
          pointed by "it" */
      }
    removeChildren();
      /* The children of this branch are removed */
    }
    for(vector<Branch*>::iterator it=children.begin();it!=children.end();it++){
      /* The vector "children" is traversed element by element with an
      iterator pointing to each element from the beginning to the end */
      (*it)->removeParent();
        /* The parenting links are removed */
    }
  }
}

/*****************************************************************************/

vector<Branch*> Branch::getChildren(){
  return children;
    /* Returns the "children" vector */
}

/*****************************************************************************/

vector<Branch*> Branch::getBrothers(){
  vector<Branch*> brothers;
    /* "brothers" is the return varaible */
  vector<Branch*> potential_brothers;
    /* "potential_brothers" is the same vector as "brothers", except it will
    contain this branch also */
  if (parent.size()!=0){
    /* If this branch has parents... (If not it is the trunk, hence it has no
    brothers) */
    potential_brothers=parent[0]->children;
      /* "potential_brothers" receives the vector "children" of the parent of
      this branch  */
    for (vector<Branch*>::iterator it=potential_brothers.begin();it!=potential_brothers.end();it++){
      /* The vector "potential_brothers" is traversed element by element with an
      iterator pointing to each element from the beginning to the end */
      if ((*it)!=this){
        /* When the "potential_brothers" element pointed by "it" isn't this
        branch, we add it to "brothers" vector */
        brothers.push_back((*it));
      }
    }
  }
  return brothers;
    /* The function returns "brothers" */
}



//////5--BRANCH HIERARCHY//////////////////////////////////////////////////////
void Branch::setStrahler(int strahler_index){
  strahler=strahler_index;
    /* The attribute "strahler" receives the value "strahler_index" */
}

/*****************************************************************************/

int Branch::getStrahler(){
    return strahler;
      /* Returns the attribute "strahler" */
}

/*****************************************************************************/

void Branch::setHorton(int horton_index){
    horton=horton_index;
      /* The attribute "horton" receives the value "horton_index" */
}

/*****************************************************************************/

int Branch::getHorton(){
    return horton;
      /* Returns the attribute "horton" */
}
