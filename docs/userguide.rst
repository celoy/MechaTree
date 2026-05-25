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
- ``wind_fn`` — a Python callable returning a 3-tuple wind vector. Two
  arities are accepted: ``(generation, rng) -> (x, y, z)`` for stateless
  winds, or ``(generation, rng, context) -> (x, y, z)`` for canopy-aware
  winds that need access to the live ``PyTree`` (in ``grow_tree``) or
  ``Forest`` (in ``Forest.step``). The default mirrors the Fortran's
  rotating storm with a long-tailed amplitude. See the
  :ref:`canopy-aware-wind` section for an out-of-the-box canopy-aware
  option.
- ``sun`` — a :class:`mechatree.light.Sun`. The default 4×8 grid samples
  the Lambert hemisphere; ``Sun.from_arrays(elev, azim)`` lets you
  specify arbitrary directions.

A worked end-to-end example showing all four levers lives in
``examples/custom_simulation.py``.

Beyond the constants, ``mechatree.genome`` also ships
:class:`mechatree.genome.NeuralSafety` and
:class:`mechatree.genome.NeuralAllocation` — direct ports of the
three-layer tanh networks evolved in Eloy et al. 2017. Load a champion
genome with :func:`mechatree.genome.load_champion`:

.. code-block:: python

   from mechatree.genome import load_champion

   safety, allocation, meta = load_champion(
       "data/S3_champions.json", species_id=0
   )
   tree = grow_tree(cfg, n_generations=100, seed=42,
                    safety=safety, allocation=allocation)

Or set ``genome.neural_from`` in the YAML config and ``load_config``
will build the neural models for you.

For a closed-form genome that sits between the constants and the neural
nets, MechaTree ships a SymPy bridge (Step 15). Any ``genome:`` scalar
in the YAML can be replaced with a string expression in the canonical
input symbols — ``nb_leaves`` and ``max_stress`` for ``safety``,
``nb_leaves`` and ``vol_relative`` for the allocation fields:

.. code-block:: yaml

   genome:
     safety:        "3 * tanh(max_stress + 0.1)"
     p_seeds:       "0.1 * tanh(vol_relative)"
     p_leaves:      0.5
     phototropism:  "0.5 + 0.1 * tanh(nb_leaves / 100)"

The bridge is gated behind the ``sympy`` optional extra
(``pip install 'mechatree[sympy]'``); expressions are parsed with
``sympy.sympify``, compiled with ``sympy.lambdify``, and dispatched per
branch via :class:`mechatree.genome.CallbackSafety` /
:class:`mechatree.genome.CallbackAllocation` — Cython-side
``with gil`` shims hand control back to Python from the C++ growth
loop. Programmatic constructors live in :mod:`mechatree.sympy_genome`:

.. code-block:: python

   from mechatree.sympy_genome import sympy_safety, sympy_allocation

   safety = sympy_safety("3 * tanh(max_stress + 0.1)")
   allocation = sympy_allocation(
       p_seeds="0.1 * tanh(vol_relative)",
       p_leaves=0.5,
       phototropism="0.5 + 0.1 * tanh(nb_leaves / 100)",
   )
   tree = grow_tree(cfg, n_generations=100, seed=42,
                    safety=safety, allocation=allocation)

See ``examples/sympy_genome.py`` + ``examples/sympy_genome.yaml`` for a
runnable recipe.


.. _canopy-aware-wind:

Canopy-aware wind via DendroFlow
================================

For a wind model that accounts for the canopy thinning the inflow
profile, MechaTree ships a bridge to `DendroFlow
<https://github.com/celoy/DendroFlow>`_'s lean
``BulkThinningBranchWindModel``. The bridge is gated behind an optional
extra::

   uv pip install -e '.[dendroflow]'
   # DendroFlow isn't on PyPI yet — install from a sibling checkout:
   uv pip install -e ../DendroFlow

YAML drives the wiring — add a ``wind:`` block:

.. code-block:: yaml

   wind:
     model: dendroflow
     U_infty:   [3.0, 4.0, 5.0, 5.8, 6.4, 6.9, 7.3, 7.6, 7.8, 8.0]
     z_centers: [0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75]
     H: 0.5
     C_D: 1.0

Per generation, the bridge snapshots the live tree (or every tree in a
``Forest``) into a DendroFlow ``Cylinders``, runs the 1-D bulk-thinning
solver, and returns the canopy-mean streamwise wind as
``(Ū, 0, 0)``. For ``Forest``, every tree's branches go into a single
``Cylinders`` call so the wake/shadow is captured. Python factory if
you'd rather skip the YAML:

.. code-block:: python

   from mechatree.wind import make_dendroflow_wind_fn

   wind_fn = make_dendroflow_wind_fn(
       U_infty=[3, 4, 5, 6, 7, 8],
       z_centers=[0.25, 0.75, 1.25, 1.75, 2.25, 2.75],
       H=0.5,
   )
   tree = grow_tree(cfg, n_generations=100, seed=42, wind_fn=wind_fn)

See ``examples/dendroflow_wind.py`` for a runnable recipe.


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
- ``examples/plot_strahler.py`` — Strahler-order diagnostics
  (self-similarity, Leonardo's rule, Tokunaga matrix).
- ``examples/dendroflow_wind.py`` + ``dendroflow_wind.yaml`` —
  canopy-aware wind via the DendroFlow bridge (Step 17). Requires the
  ``dendroflow`` optional extra.
- ``examples/sympy_genome.py`` + ``sympy_genome.yaml`` — closed-form
  safety + allocation via SymPy expressions (Step 15). Requires the
  ``sympy`` optional extra.
- ``examples/forest.yaml`` — annotated configuration with the Fortran
  parameter names in comments.

Annotated notebooks for the same tutorials live under ``notebooks/``,
with prose between the cells:

- ``notebooks/01_grow_one_tree.ipynb`` — companion to
  ``examples/grow_one_tree.py``.
- ``notebooks/02_forest_under_wind.ipynb`` — population & biomass
  dynamics, self-thinning.
- ``notebooks/03_neural_genome.ipynb`` — load an evolved S3 champion
  and compare side-by-side with the default constant genome.
- ``notebooks/04_custom_growth_law.ipynb`` — plug in custom ``wind_fn``,
  ``Sun``, ``ConstantSafety``, ``ConstantAllocation``.
- ``notebooks/05_strahler_diagnostics.ipynb`` — self-similar branching
  analysis.
- ``notebooks/06_fractal_dimension.ipynb`` — reproduces SI Fig. S8 of
  Eloy et al. 2017 for one of the S3 champion genomes: Horton ratios
  :math:`R_n, R_l, R_d, R_a, D`, time evolution, branch tapering, and
  area conservation.

Install the notebook extra (``uv pip install -e ".[notebooks]"``) and
launch with ``uv run jupyter lab``. Output cells are stripped by
``nbstripout`` before commit; re-execute to populate the plots.

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
