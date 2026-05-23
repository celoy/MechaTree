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


Simulating a single tree
========================

The :func:`mechatree.simulate.grow_tree` orchestrator reproduces
``legacy_fortran/tree.f90``: starting from a one-twig seedling it
iterates light → stresses → growth → pruning → branching → reserve
accounting once per generation. Configuration lives in a YAML file with
defaults from the original Forest.ini.

.. code-block:: python

   from mechatree.config import load_config
   from mechatree.simulate import grow_tree

   cfg = load_config("examples/forest.yaml")
   tree = grow_tree(cfg, n_generations=100, seed=42)

   print(f"{tree.get_number_of_branches()} branches, "
         f"{tree.get_total_leaves()} leaves")

To inspect intermediate state, pass an ``on_step`` callback:

.. code-block:: python

   def cb(generation, tree):
       print(generation, tree.get_number_of_branches())

   grow_tree(cfg, n_generations=100, seed=42, on_step=cb)

A 3-D rendering uses :func:`mechatree.plotting.plot_tree_3d`, which
reads the typed mechanics fields:

.. code-block:: python

   from mechatree.plotting import plot_tree_3d
   plot_tree_3d(tree)


Simulating a forest
===================

For many trees competing for light, use :class:`mechatree.forest.Forest`
— a direct port of ``legacy_fortran/Forest.f90``. Trees are placed
uniformly across a disk; light interception runs over the union of all
their leaves, so cross-tree shading happens automatically.

.. code-block:: python

   from mechatree.config import load_config
   from mechatree.forest import Forest

   cfg = load_config("examples/forest.yaml")
   forest = Forest(cfg, seed=42)

   def cb(generation, forest, stats):
       print(generation, stats.n_trees, stats.biomass_total)

   forest.run(n_generations=100, on_step=cb)

Each call to ``forest.step(gen)`` returns a
:class:`mechatree.forest.ForestStats` with branch / leaf totals, biomass,
and births / deaths for the step.

The death rule (Fortran default: a tree dies if
``n_branches < 11 AND age > 5``, OR ``age > 1000``) is configurable on
``ForestConfig`` — handy for studying self-thinning sensitivity.


Customizing the simulation
==========================

``grow_tree`` (and ``Forest``) accept optional overrides on top of the
YAML config:

- ``safety`` — a :class:`mechatree.genome.SafetyModel`; default is
  ``ConstantSafety(1.0)``. Pass ``ConstantSafety(value)`` to bias the
  tree toward more or less mechanical headroom.
- ``allocation`` — a :class:`mechatree.genome.AllocationModel`; default
  is ``ConstantAllocation(p_seeds=0.1, p_leaves=0.5, phototropism=0.5)``.
- ``wind_fn`` — a Python callable ``(generation, rng) -> (x, y, z)``.
  The default mirrors the Fortran's rotating storm with a long-tailed
  amplitude.
- ``sun`` — a :class:`mechatree.light.Sun`. The default 4×8 grid samples
  the Lambert hemisphere; ``Sun.from_arrays(elev, azim)`` lets you
  specify arbitrary directions.

A worked end-to-end example showing all four levers lives in
``examples/custom_simulation.py``.

.. note::

   Non-constant ``SafetyModel`` / ``AllocationModel`` subclasses are not
   yet pluggable from Python — the C++ side dispatches through a virtual
   ``compute()`` and needs a concrete subclass. A neural-network genome
   port is a self-contained later step; until it lands, your runtime
   knob is the constant value plus the wind / sun functions.


Examples
========

Runnable demos live under ``examples/`` at the repository root.

Simulation tutorials (Step 11+):

- ``examples/grow_one_tree.py`` — load a YAML config, grow a tree, plot
  it in 3-D.
- ``examples/grow_a_forest.py`` — drive a :class:`Forest`, plot biomass
  over time and a top-down map of the final stand.
- ``examples/custom_simulation.py`` — three side-by-side runs with
  different wind / sun / genome settings.
- ``examples/forest.yaml`` — annotated configuration with the Fortran
  parameter names in comments.

Legacy topology demos (Step 3):

- ``examples/random_growth.py`` — 3D random growth, rendered with
  :func:`mechatree.plotting.plot_3d`.
- ``examples/self_avoiding.py`` — 2D coral-like growth with
  self-avoidance, using :func:`mechatree.geometry.distance_test` and
  saving snapshots via :func:`mechatree.plotting.plot_2d`.
- ``examples/sap_transport.py`` — resource-allocation model where
  branches grow and prune based on a sap-reserve balance.

Each script exposes ``--iterations N`` and ``--seed N`` flags:

.. code-block:: bash

   uv run python examples/grow_one_tree.py --iterations 100 --seed 42
