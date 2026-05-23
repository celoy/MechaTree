=================
Users guide
=================

In this section, we will only expose practical information allowing the reader to make his first steps with PyTreeLib avoiding to enter in the implementation details of the library. This knowledge should be sufficient to use the plateform with purely simulation interests. The curious reader, or the one wanting to contribute with the developpement of the library, will found the C++ implementation details in the section "Developpers guide". 

PyTree class
============

.. py:function:: add_property(self, int i, str prop, double val) 

   Adds a property named ``prop`` with value ``val``, to the branch of PyTree with index ``i``. If the branch doesn't exist, nothing is done.

.. py:function:: get_property(self, int i, str prop)

   Returns the value of the property named ``prop`` of the branch of PyTree with index ``i``. If the branch, or the property, don't exist, this method returns 0.

.. py:function:: set_property(self, int i, str prop, double val)

   Sets the value of the already existing property named ``prop`` of the branch of PyTree with index ``i`` to the value ``val``. If the branch, or the property don't exist, this method returns 0.

.. py:function:: get_number_of_branches(self)

   Returns the number of branches of PyTree.

.. py:function:: add_branch(self, int i, dict props)
   
   Adds a new branch with the properties contained in the ``dictionnary props`` to PyTree. The branch is added at the index ``i+1``, where ``i`` is the index of the parent. The creation of family links is dealt with in C++ routines. If there is no branch located at ``i``, nothing is done.

.. py:function:: remove_branch(self, int i)
   
   Removes the branch of PyTree with index ``i`` and all of its descendance. The removal of family links is dealt with in C++ routines. If there is no branch located at ``i``, nothing is done.

.. py:function:: get_last_descendant_index(self, int i)

   Returns the index of the last descendant of the PyTree branch located at ``i``. If there is no branch located at ``i``, this method returns 0.

.. py:function:: get_parent_index(self, int i)

   Returns the index of the parent of the PyTree branch located at ``i``. If the trunk is located at ``i``, this method returns -1. If there is no branch located at ``i``, this method returns -2.

.. py:function:: get_brothers_index(self, int i)

   Returns a list containing the indexes of the brothers of the PyTree branch located at ``i``. If there is no branch located at ``i``, this method returns an empty list.

.. py:function:: get_children_index(self, int i)

   Returns a list containing the indexes of the children of the PyTree branch located at ``i``. If there is no branch located at ``i``, this method returns an empty list.

.. py:function:: has_parent(self, int i)

   Returns 1 if the PyTree branch located at ``i`` has a parent, and 0 if it doesn't. If there is no branch located at ``i``, this method returns -1.

.. py:function:: get_number_of_children(self, int i)

   Returns the number of children of the PyTree branch located at ``i``. If there is no branch located at ``i``, this method returns -1.

.. py:function:: set_strahler(self)

   Effectuates the Strahler classification of PyTree.

.. py:function:: get_strahler(self, int i)

   Returns the Strahler order of the PyTree branch located at ``i``. If there is no branch located at ``i``, this method returns 0.

.. py:function:: get_strahler_distribution(self)

   Returns a ``dictionnary`` where the keys are the Strahler orders and the values the number of branches per order.

.. py:function:: set_horton(self)

   Effectuates the Horton classification of PyTree.

.. py:function:: get_horton(self, int i)

   Returns the Horton order of the PyTree branch located at ``i``. If there is no branch located at ``i``, this method returns 0.

.. py:function:: get_horton_distribution(self)

   Returns a ``dictionnary`` where the keys are the Horton orders and the values the number of branches per order.

.. py:function:: mean_agg_prop_s(self, str prop)

   Returns a ``dictionnary`` where the keys are the Strahler orders and the values the mean value of a branch property named ``prop`` per order.

.. py:function:: mean_agg_prop_h(self, str prop)

   Returns a ``dictionnary`` where the keys are the Horton orders and the values the mean value of a branch property named ``prop`` per order.

Basic tutorial
==============

Import the library
------------------

To use PyTree objects from Python you must import the ``PyTree`` module:

``from pytree import PyTree``


Create a tree
-------------

To initialize an instance of the PyTree class we need to provide the set of properties of the first branch in the form of a Python ``dictionnary``:

``properties={"prop1":val1, "prop2":val2, ... , "propN":valN}``

where the names of the properties are written between inverted commas and their value are real numbers. The PyTree object named ``tree`` can now be created:

``tree=PyTree(properties)``.

Add and remove branches
-----------------------

As with the first branch, to add new branches to the tree, we must provide a set of properties in the form of a ``dictionnary``. We must also pass the index n of the parent as the first argument of the ``add_branch`` function:

``tree.add_branch(n,properties)``.  

To remove a branch, we just hace to indicate the index of the branch we desire to remove: 

``tree.remove_branch(n)``

this will automatically remove all its descendance.

Traverse the tree
-----------------

The simplest way to traverse the tree is to iterate over its branches. To do that you can acces the number of branches of the tree with the function ``get_number_of_branches``:

``for i in range tree.get_number_of_branches()``

Acces and modify branches properties
------------------------------------

We provide three methods to interact with the branches properties. To access the value of a property we can use the function ``get_property``:

``tree.get_property(n,"prop")``

where ``n`` is the index of the branch and ``"prop"`` the name of the property to access. To modify the value of a property we can use the function ``set_property``:

``tree.set_property(n,"prop",val)``

where ``val`` is the new value of the property and a real number. Finally to add a new property to a branch we can use the function ``add_property``: 

``tree.add_property(n,"prop",val)``.

Code examples
=============



