/*
This is the declaration of the branch class.
*/

#ifndef BRANCH_H_
#define BRANCH_H_

#include <memory>
#include <string>
#include <vector>
#include <map>
using namespace std;

class Branch{

  /*
  This is a formalization of the concept of a branch. Branches are seen as the
  elementary constituents of a tree. Branches are characterized by a set of
  properties, a set of family links (with its parents and children), and a set
  of orders characterizing their place in the tree hierarchy. It can be seen as
  the generalization of the branch class in a binary tree.
  */

private:
      /*
      All the attributes of the class are private and hence cannot be accesed
      from external programs.

      1) "properties" is a map containing the set of properties of a branch. A
      property is given by a couple name-value. The name is of type string and
      the value of type double.

      2) "parent" is a vector containing only one pointer to a branch which
      is the parent of the current branch. The pointer is stocked in a vector
      in order to use the memory management of that built-in class.

      3) "children" is a vector containing pointers to the branches which are
      the children of the current branch. The children can be dinamically
      allocated.

      4) "strahler" is an integer representing the strahler order of the current
      branch

      5) "horton" is an integer representing the horton order of the current
      branch

      */

      map<string, double> properties;

      vector<Branch*> parent;
      vector<Branch*> children;

      int strahler;
      int horton;

public:

      /*
      All the methods of the class are public and are designed to manipulate
      the private attributes. This methods can be accessed from an external
      program.

      The methods can be classified in the following way:

      1) Class constructor
      2) Interact with branch properties (add, get & set)
      3) Interact with branch parent (has, add, remove & get)
      4) Interact with branch children (has, add, remove & get)
      5) Interact with branch orders (set & get)

      */

//////CONSTRUCTOR//////////////////////////////////////////////////////////////
      Branch();
        /* nullary constructor */
//////PROPERTIES///////////////////////////////////////////////////////////////
      void addProperty(string name, double value );
        /* adds the property "name" with value "value" */
      void setProperty(string name, double value);
        /* changes the value of property "name" to "value"*/
      void setProperties(map<string,double> props);
        /* attribute properties receives the map "props" */
      double getProperty(string name);
        /* returns the value of the property "name" */

//////PARENT///////////////////////////////////////////////////////////////////
      bool hasParent();
        /* returns True if branch has parent, False if it doesn't */
      bool hasParent(Branch*test);
        /* returns True if "*test" is the parent, False if it isn't */
      void addParent(Branch*p);
        /* adds "*p" to the vector "parent" */
      void removeParent(Branch*p);
        /* removes "*p" from the vector "parent" */
      void removeParent();
        /* removes the parent */
      Branch* getParent();
        /* returns a pointer to the parent */

//////CHILDREN/////////////////////////////////////////////////////////////////
      bool hasChildren();
        /* returns True if branch has children, False if it doesn't */
      bool hasChild(Branch*test);
        /* returns True if "*test" is a child, False if it isn't */
      void addChild(Branch*ch);
        /* adds "*ch" to the vector "children" */
      void removeChild(Branch*ch);
        /* removes "*ch" from the vector "children" */
      void removeChildren();
        /* removes all the elements of the vector "children" */
      void removeDescendants();
        /* removes all the elements of the vector "children" and repeats the
         operation for all the children in the descendance */
      vector<Branch*> getChildren();
        /* returns the vector "children" */
      vector<Branch*> getBrothers();
        /* returns a vector containing pointers to all the brothers of the
         branch */

//////BRANCH HIERARCHY/////////////////////////////////////////////////////////
      void setStrahler(int strahler_index);
        /* sets the attribute "strahler" to "strahler_index" */
      int getStrahler();
        /* returns the value of the attribute "strahler" */
      void setHorton(int horton_index);
        /* sets the attribute "horton" to "horton_index" */
      int getHorton();
        /* returns the value of the attribute "horton" */

};

#endif
