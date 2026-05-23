===============
Developer guide
===============

This section documents the C++ classes underlying MechaTree's Cython
bindings. The headers live in ``src/mechatree/_core/``.


Branch class
============

A ``Branch`` object is described by a set of properties of type ``double``
and by a set of pointers to its children and parent. The methods of the
class allow us to add, set, or get properties of a ``Branch``, and to add,
get, or remove family links between objects.

.. cpp:function:: void addProperty(string name, double value)

   Adds a property named ``name`` with value ``value`` to the set of
   properties of the ``Branch`` object.

.. cpp:function:: void setProperty(string name, double value)

   If a property named ``name`` exists, updates its value to ``value``.

.. cpp:function:: void setProperties(map<string,double> props)

   Sets the entire property map to ``props``.

.. cpp:function:: double getProperty(string name)

   If the property ``name`` exists, returns its value.

.. cpp:function:: bool hasParent()

   Returns ``true`` if the ``parent`` vector isn't empty, ``false`` otherwise.

.. cpp:function:: bool hasParent(Branch *test)

   Returns ``true`` if the ``parent`` vector contains ``test``.

.. cpp:function:: void addParent(Branch *p)

   Adds ``p`` to the ``parent`` vector if not already present.

.. cpp:function:: void removeParent(Branch *p)

   Removes ``p`` from the ``parent`` vector if present.

.. cpp:function:: void removeParent()

   Removes the first element of the ``parent`` vector if non-empty.

.. cpp:function:: Branch* getParent()

   Returns the first element of the ``parent`` vector.

.. cpp:function:: bool hasChildren()

   Returns ``true`` if the ``children`` vector isn't empty.

.. cpp:function:: bool hasChild(Branch* test)

   Returns ``true`` if the ``children`` vector contains ``test``.

.. cpp:function:: void addChild(Branch* ch)

   Adds ``ch`` to the ``children`` vector if not already present.

.. cpp:function:: void removeChild(Branch* ch)

   Removes ``ch`` from the ``children`` vector if present.

.. cpp:function:: void removeChildren()

   Removes all elements from the ``children`` vector.

.. cpp:function:: void removeDescendants()

   Removes all descendants of the ``Branch`` object.

.. cpp:function:: vector<Branch*> getChildren()

   Returns the ``children`` vector.

.. cpp:function:: vector<Branch*> getBrothers()

   Returns a vector containing all the brothers of the ``Branch`` object.

.. cpp:function:: void setStrahler(int strahler_index)

   Sets the ``strahler`` attribute to ``strahler_index``.

.. cpp:function:: int getStrahler()

   Returns the ``strahler`` attribute.

.. cpp:function:: void setHorton(int horton_index)

   Sets the ``horton`` attribute to ``horton_index``.

.. cpp:function:: int getHorton()

   Returns the ``horton`` attribute.


Tree class
==========

A ``Tree`` object is a vector containing pointers to ``Branch`` objects.
Each branch can be identified by an index — its position in the vector.

.. cpp:function:: Tree(map<string,double> trunk_props)

   Non-nullary constructor. The ``Tree`` object is initialized as a single
   branch (the trunk) with ``trunk_props`` as its properties.

.. cpp:function:: int getNumberOfBranches()

   Returns the number of elements in the ``tree_branches`` vector.

.. cpp:function:: map<int,int> getStrahlerDistribution()

   Returns a map of (Strahler order → branch count). Two contiguous
   branches of the same order are merged into one.

.. cpp:function:: map<int,int> getHortonDistribution()

   Returns a map of (Horton order → branch count). Two contiguous branches
   of the same order are merged into one.

.. cpp:function:: map<int,double> meanAggregativePropS(string name)

   Returns a map of (Strahler order → mean value of property ``name``).

.. cpp:function:: map<int,double> meanAggregativePropH(string name)

   Returns a map of (Horton order → mean value of property ``name``).

.. cpp:function:: Branch* getTrunk()

   Returns the first element of ``tree_branches``.

.. cpp:function:: Branch* getSummit()

   Returns the last element of ``tree_branches``.

.. cpp:function:: Branch* getBranch(int branch_index)

   Returns the element at position ``branch_index``. If invalid, returns a
   pointer to a nullary-constructed ``Branch``.

.. cpp:function:: int getIndex(Branch* branch)

   Returns the position of ``branch`` in ``tree_branches``, or 0 if absent.

.. cpp:function:: void setStrahler()

   Sets the Strahler order for every element of ``tree_branches``.

.. cpp:function:: int getStrahler(int index)

   Returns the Strahler order of the branch at ``index``. Returns 0 if
   ``index`` is invalid.

.. cpp:function:: void setHorton()

   Sets the Horton order for every element of ``tree_branches``.

.. cpp:function:: int getHorton(int index)

   Returns the Horton order of the branch at ``index``. Returns 0 if
   ``index`` is invalid.

.. cpp:function:: int getLastDescendantIndex(int ancestor_index)

   Returns the maximum of the indexes of the descendants of the branch at
   ``ancestor_index``. Returns 0 if ``ancestor_index`` is invalid.

.. cpp:function:: int getParentIndex(int child_index)

   Returns the index of the first element of the ``parent`` vector of the
   branch located at ``child_index``. If ``child_index`` is the trunk,
   returns -1. If invalid, returns -2.

.. cpp:function:: vector<int> getBrothersIndex(int branch_index)

   Returns a vector containing the indexes of the brothers of the branch
   at ``branch_index``. Returns an empty vector if invalid.

.. cpp:function:: vector<int> getChildrenIndex(int parent_index)

   Returns a vector containing the indexes of the children of the branch
   at ``parent_index``. Returns an empty vector if invalid.

.. cpp:function:: int hasParent(int index)

   Returns 1 if the branch at ``index`` has a parent, 0 if not. Returns
   -1 if ``index`` is invalid.

.. cpp:function:: int getNumberOfChildren(int parent_index)

   Returns the number of children of the branch at ``parent_index``.
   Returns -1 if invalid.

.. cpp:function:: void addBranch(int parent_index)

   Adds a property-less branch at position ``parent_index + 1`` and wires
   up parent/child links. If ``parent_index`` is invalid, nothing is done.

.. cpp:function:: void addBranch(int parent_index, map<string,double> props)

   Adds a branch with the properties ``props`` at position
   ``parent_index + 1`` and wires up parent/child links. If
   ``parent_index`` is invalid, nothing is done.

.. cpp:function:: void removeBranch(int branch_index)

   Removes the element at ``branch_index`` and all of its descendants
   (i.e. erases all elements between ``branch_index`` and its last
   descendant). The parent/child links are removed too. If
   ``branch_index`` is invalid, nothing is done.

.. cpp:function:: void addProperty(int index, string name, double value)

   Adds an entry ``name`` → ``value`` to the property map of the branch at
   ``index``. If ``index`` is invalid, nothing is done.

.. cpp:function:: double getProperty(int index, string name)

   Returns the value of the property ``name`` of the branch at ``index``.
   Returns 0 if ``index`` is invalid.

.. cpp:function:: void setProperty(int index, string name, double value)

   Sets the value of the property ``name`` of the branch at ``index`` to
   ``value``. If ``index`` is invalid, nothing is done.
