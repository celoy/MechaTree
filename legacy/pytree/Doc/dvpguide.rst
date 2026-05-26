=====================
Developper guide
=====================

Branch class
============

A **Branch** object is described by a set of properties of type *double* and by a set of pointers to its children and parent. The methods of the class allow us to *add*, *set* or *get* properties of a **Branch** object, and to *add*, *get* or *remove* family links beetween objects.  

.. cpp:function:: void addProperty(string name, double value)
   
   Adds a property named "name" and with value "value" to the set of properties of the **Branch** object.

.. cpp:function:: void setProperty(string name, double value)
   
   If a property named "name" exists, this function updates its value to "value". 

.. cpp:function:: void setProperties(map<string,double> props)	

   Sets the set of properties to the map "props".

.. cpp:function:: double getProperty(string name)

   If the property "name" exists, it returns its value.

.. cpp:function:: bool hasParent()

   Returns *true* if the **vector** "parent" isn't empty and *false* if it is.

.. cpp:function:: bool hasParent(Branch *test)

   Returns *true* if the **vector** "parent" contains "test" and *false* if it doesn't.

.. cpp:function:: void addParent(Branch *p)

   Checks if the **vector** "parent" contains "p" and adds it if it doesn't.

.. cpp:function:: void removeParent(Branch *p)

   Checks if the **vector** "parent" contains "p" and removes it if it does.	   

.. cpp:function:: void removeParent()  

   Checks if the **vector** "parent" isn't empty and if it isn't the function removes the first element of the **vector**.

.. cpp:function:: Branch* getParent()

   Returns the first element of the **vector** "parent".

.. cpp:function:: bool hasChildren()

   Returns *true* if the **vector** "children" isn't empty and *false* if it does.

.. cpp:function:: bool hasChild(Branch* test)

   Returns *true* if the **vector** "children" contains "test" and *false* if it doesn't.

.. cpp:function:: void addChild(Branch* ch)

   Checks if the **vector** "children" contains "ch" and adds it if it doesn't.

.. cpp:function:: void removeChild(Branch* ch)

   Checks if the **vector** "children" contains "ch" and removes it if it does.

.. cpp:function:: void removeChildren()

   Checks if the **vector** "children" isn't empty and if it isn't the function removes all the elements from the **vector**.

.. cpp:function:: void removeDescendants()

   Removes all the descendants of the **Branch** object.

.. cpp:function:: vector<Branch*> getChildren()

   Returns the **vector** "children" of the **Branch** object.

.. cpp:function:: vector<Branch*> getBrothers()

   Returns a **vector** containing all the brothers of the **Branch** object.

.. cpp:function:: void setStrahler(int strahler_index)

   Sets the attribute "strahler" value to "strahler_index".

.. cpp:function:: int getStrahler()

   Returns the **int** "strahler".

.. cpp:function:: void setHorton(int horton_index)

   Sets the attribute "horton" value to "horton_index".

.. cpp:function:: int getHorton()

   Returns the **int** "horton".



Tree class
==========

A **Tree** object is a **vector** containing pointers to **Branch** objects. Each branch of the **Tree** can be identified with an index, which is its position on the vector.

.. cpp:function:: Tree(map<string,double> trunk_props)

   Non nullary constructor. The **Tree** object is initialized as a single branch with "trunk_props" as properties. 

.. cpp:function:: int getNumberOfBranches()

   Returns the number of elements on the **vector** "tree_branches".

.. cpp:function:: map<int,int> getStrahlerDistribution()

   Returns a **map** containing the number of branches per order when the branches are ordered using the Strahler classification. The order is given by the *key* and the number of branches by the *value*. Two contiguous branches of same order become te same branch.

.. cpp:function:: map<int,int> getHortonDistribution()

   Returns a **map** containing the number of branches per order when the branches are ordered using the Horton classification. The order is given by the *key* and the number of branches by the *value*. Two contiguous branches of same order become te same branch.

.. cpp:function:: map<int,double> meanAggregativePropS(string name)

   Returns a **map** containing the mean value of an aggregative property "name" (for example the branch length) in function of the order of the branches, when the branches are ordered using the Strahler classification. The order is given by the *key* and the mean value of the aggregative property by the *value*.

.. function:: map<int,double> meanAggregativePropH(string name)  

   Returns a **map** containing the mean value of an aggregative property "name" (for example the branch length) in function of the order of the branches, when the branches are ordered using the Horton classification. The order is given by the *key* and the mean value of the aggregative property by the *value*.

.. cpp:function:: Branch* getTrunk()

   Returns the first element of the **vector** "tree_branches".

.. cpp:function:: Branch* getSummit()

   Returns the last element of the **vector** "tree_branches".	

.. cpp:function:: Branch* getBranch(int branch_index)

   Returns the element at the position "branch_index" of the **vector** "tree_branches". If "branch_index" doesn't exist, the function returns a pointer to a nullary-constructed **Branch** object.

.. function:: int getIndex(Branch* branch)

   If the **vector** "tree_branches" contains the element "branch", its position on the **vector** is  returned. Else, the function returns 0.

.. cpp:function:: void setStrahler()

   Sets the Strahler order for every element of the **vector** "tree_branches".

.. cpp:function:: int getStrahler(int index)

   Returns the Strahler order of the element of "tree_branches" at the position "index". If "index" doesn't exist, the function returns 0.

.. cpp:function:: void setHorton()

   Sets the Horton order for every element of the **vector** "tree_branches".

.. cpp:function:: int getHorton(int index)

   Returns the Horton order of the element of "tree_branches" at the position "index". If "index" doesn't exist, the function returns 0.

.. cpp:function:: int getLastDescendantIndex(int ancestor_index)

   Returns the maximum of the indexes of the descendants of the "tree_branches" element located at "ancestor_index". If "ancestor_index" doesn't exist, the function returns 0.

.. cpp:function:: int getParentIndex(int child_index)

   Returns the index of the first element of the **vector** parent of the branch located at "child_index". If "child_index" is the index of the first branch, the function returns -1. If it doesn't exists the function returns -2.

.. cpp:function:: vector<int> getBrothersIndex(int branch_index)

   Returns a **vector** containing the indexes of the brothers of the branch located at "branch_index". If "branch_index" doesn't exist, the function returns an empty **vector**.

.. cpp:function:: vector<int> getChildrenIndex(int parent_index)

   Returns a **vector** containing the indexes of the children of the branch located at "parent_index". If "parent_index" doesn't exist, the function returns an empty **vector**.	

.. cpp:function:: int hasParent(int index)

   Returns 1 if the branch at "index" has a parent and 0 if it doesn't. If "index" doesn't exist, the function returns -1.

.. cpp:function:: int getNumberOfChildren(int parent_index)

   Returns the number of elements of the **vector** "children" of the branch located at "parent_index". If "parent_index" doesn't exist, the function returns -1.

.. cpp:function:: void addBranch(int parent_index)

   Adds a branch without any properties to the **vector** "tree_branches". The new element is added at the position "parent_index+1". The parenting links between the branches are also generated here. If "parent_index" doesn't exist, nothing is done.

.. cpp:function:: void addBranch(int parent_index, map<string,double> props)

   Adds a branch with the properties "props" to the **vector** "tree_branches". The new element is added at the position "parent_index+1". The parenting links between the branches are also generated here. If "parent_index" doesn't exist, nothing is done.

.. cpp:function:: void removeBranch(int branch_index)

   Removes the element at the position "branch_index" of the **vector** "tree_branches" and all of its descendants. Because of the way the branches are added to the tree, this comes down to erase all the elements of "tree_branches" located between "bracnh_index" and its last descendant index. The parenting links between the branches are removed here. If "branch_index" doesn't exist, nothing is done.  

.. cpp:function:: void addProperty(int index, string name, double value)   

   Adds an element of *key*  "name" and *value* "value" to the **map** "properties" of the branch located at "index". If "index" doesn't exist, nothing is done.

.. cpp:function:: double getProperty(int index, string name)

   Returns the value of the property "name" of the branch located at "index". If "index" doesn't exist, the function returns 0.

.. cpp:function:: void setProperty(int index, string name, double value)

   Sets the value of the property "name" of the branch located at "index" at "value". If "index" doesn't exist, nothing is done.


Interfacing with Cython, PyTree class in depth
==============================================
