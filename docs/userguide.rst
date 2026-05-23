==========
User guide
==========

This section is a practical guide to your first steps with MechaTree.
For the full method-by-method reference, see :doc:`api/mechatree`.
For the C++ internals, see :doc:`dvpguide`.


Basic tutorial
==============

Import the library
------------------

To use ``PyTree`` objects from Python, import the class from the package:

.. code-block:: python

   from mechatree import PyTree


Create a tree
-------------

To initialize an instance of ``PyTree``, provide the set of properties of
the first branch (the trunk) as a Python dictionary:

.. code-block:: python

   properties = {"length": 1.0, "radius": 0.1}
   tree = PyTree(properties)

The property names are strings; their values are floating-point numbers.


Add and remove branches
-----------------------

As with the trunk, to add new branches we provide a set of properties and
the index ``n`` of the parent:

.. code-block:: python

   tree.add_branch(n, properties)

To remove a branch, we indicate the index of the branch to remove:

.. code-block:: python

   tree.remove_branch(n)

This automatically removes all of its descendants as well.


Traverse the tree
-----------------

The simplest way to traverse the tree is to iterate over its branches.
``get_number_of_branches`` returns the current count:

.. code-block:: python

   for i in range(tree.get_number_of_branches()):
       ...


Access and modify branch properties
-----------------------------------

To access the value of a property:

.. code-block:: python

   tree.get_property(n, "prop")

To modify the value of a property:

.. code-block:: python

   tree.set_property(n, "prop", val)

To add a new property to an existing branch:

.. code-block:: python

   tree.add_property(n, "prop", val)


Strahler and Horton classifications
-----------------------------------

The library provides two stream-order classifications for analyzing tree
topology. After populating the tree, call ``set_strahler()`` or
``set_horton()`` once to compute the orders, then query individual
branches with ``get_strahler(i)`` / ``get_horton(i)``, or get a full
distribution with ``get_strahler_distribution()`` /
``get_horton_distribution()``.

.. code-block:: python

   tree.set_strahler()
   print(tree.get_strahler_distribution())


Examples
========

Runnable demos live under ``examples/`` at the repository root:

- ``examples/random_growth.py`` — 3D random growth, rendered with
  :func:`mechatree.plotting.plot_3d`.
- ``examples/self_avoiding.py`` — 2D coral-like growth with
  self-avoidance, using :func:`mechatree.geometry.distance_test` and
  saving snapshots via :func:`mechatree.plotting.plot_2d`.
- ``examples/sap_transport.py`` — resource-allocation model where
  branches grow and prune based on a sap-reserve balance.

Each script exposes ``--iterations N`` and ``--seed N`` flags:

.. code-block:: bash

   uv run python examples/random_growth.py --iterations 100 --seed 42
