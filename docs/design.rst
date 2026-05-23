==============================
Simulator port — design & API
==============================

This document is the design contract for the simulator port (modernization
steps 9–13: mechanics, growth, pruning, light, single-tree loop, forest).
It pins down data shapes, the YAML configuration schema, the public API
for composing tree functionalities, the genome callback signatures, and
the interface between the C++ core and the Python orchestrators. No
production code is written at this stage — the design exists so it can
be reviewed before any C++ is touched.

The reference implementation is Fortran 90, preserved verbatim under
``legacy_fortran/``. The physics is documented in Eloy *et al.*,
*Nat. Commun.* 8:1014 (2017).

Out of scope at this design stage:

* **Evolution** (``EvoluAlgo.f90``, ``mod_evolu.f90``) — per the design
  principles, evolution sits outside the core library and is deferred to
  a separate package or notebook.
* **GPU / vectorised light** — premature; the first port is CPU-only.
* **Removal of the existing ``properties`` map on ``Branch``** — kept
  indefinitely for user extensions.


Data shapes
===========

Branch (C++ struct, extended)
-----------------------------

The current ``Branch`` (see ``src/mechatree/_core/branch.h``) stores
topology — parent, children, Strahler / Horton orders — plus a
``std::unordered_map<std::string, double>`` for user-defined properties.
Step 9 adds **typed fields** for the mechanics state that the Fortran
inner loops read on every branch:

.. list-table::
   :header-rows: 1
   :widths: 22 25 28 25

   * - Field
     - C++ type
     - Fortran source
     - Role
   * - ``length``, ``diameter``
     - ``double``
     - ``branch%length``, ``%diameter``
     - allometry, growth
   * - ``light``
     - ``double``
     - ``branch%light``
     - per-step photosynthetic input ∈ [0, 1]
   * - ``stress``, ``max_stress``
     - ``double``
     - ``branch%stress``, ``%max_stress``
     - per-wind-direction + worst-of-4
   * - ``vol_growth``, ``vol_summed``
     - ``double``
     - ``branch%vol_growth``, ``%vol_summed``
     - secondary-growth bookkeeping
   * - ``nb_leaves``
     - ``int``
     - ``branch%nb_leaves``
     - leaf count under this branch
   * - ``location``
     - ``std::array<double, 3>``
     - ``branch%location``
     - tip position
   * - ``unit_t``, ``unit_b``
     - ``std::array<double, 3>``
     - ``branch%unit_t``, ``%unit_b``
     - tangent / binormal frame
   * - ``force``, ``moment``
     - ``std::array<double, 3>``
     - ``branch%T``, ``%M``
     - mechanics intermediates

Two rules govern the choice:

1. **The typed fields are the hot path.** Mechanics, light, growth and
   pruning never read the property map.
2. **The property map stays.** It is the extension point for user
   simulations and keeps the existing topology examples
   (``examples/random_growth.py``, ``self_avoiding.py``, ``sap_transport.py``)
   working unchanged.


Leaves (struct of arrays, new)
------------------------------

Leaves are *not* modelled as a Python ``list[Leaf]``. The light module
loops over thousands of leaves across 32 directions every generation; a
list of Python objects would burn most of its time on attribute lookup.
The new representation is a struct of NumPy arrays:

.. code-block:: python

   from dataclasses import dataclass
   import numpy as np

   @dataclass
   class Leaves:
       location: np.ndarray             # (n, 3),       float64
       diameter: np.ndarray             # (n,),         float64
       transparency: np.ndarray         # (n,),         float64
       branch_index: np.ndarray         # (n,),         int32
       light_per_direction: np.ndarray  # (n, n_dir),   float64

``branch_index`` is the back-pointer that lets ``aggregate_onto_trees``
write light back onto the right ``Branch``. ``light_per_direction`` is
populated in place by ``intercept``.

Leaves are extracted on demand from a tree (or concatenated across trees
in a forest) by ``extract_leaves(tree | trees) -> Leaves``. The light
module takes a ``Leaves`` — it does not depend on ``PyTree``. That is
what makes "light is decoupled from PyTree" mechanical, not aspirational.

A single ``Leaf`` namedtuple is offered as a convenience for unit tests
and small examples; it is not used in the hot path.


YAML configuration schema
=========================

The schema is a direct translation of ``legacy_fortran/Forest.ini``,
regrouped semantically and renamed to ``snake_case``. The Fortran names
appear as comments so cross-referencing the paper or legacy code stays
easy.

.. code-block:: yaml

   tree:
     twig_length: 1.0              # was TwigLength
     twig_diameter: 0.1            # was TwigDiameter
     leaf_surface: 0.25            # was LeafSurface
     cauchy: 2.0e-5                # was Cauchy
     volume_ratio_leaf: 8.0        # was VolumeRatioLeaf
     maintenance_h: 0.02           # was MaintenanceH
     max_branches: 10000           # was Nmax

   light:
     size_leaf: 1.0                # was SizeLeaf; light-grid cell size
     n_elevations: 4               # was hard-coded in Fortran
     n_azimuths: 8                 # was hard-coded in Fortran

   forest:
     n_generations: 3000           # was Ngeneration
     size: 100.0                   # was SizeForest
     n_trees_init: 10000           # was Ntrees_ini
     n_trees_max: 100000           # was Ntrees_max

   evolution:
     p_mutation: 0.05              # was Pmutation
     std_mutation: 0.005           # was STDmutation

   io:
     save_file: 'ZZZ'              # was save_file
     save_rate: 1000               # was save_rate

   init:
     tree_init: 'BestShort'        # was tree_init
     f2_min: 0.0                   # was F2min
     f2_max: 100000.0              # was F2max
     from_forest: true             # was FromForest

Validation is one frozen dataclass per group, plus a top-level
``Config`` aggregator with ``Config.from_yaml(path)``. Range checks live
in each dataclass's ``__post_init__``. **No new runtime dependency** —
no pydantic, no jsonschema. An annotated example ships at
``examples/forest.yaml``.


Genome callbacks
================

The two Fortran decision functions ``neural_branch`` and ``neural_reserve``
(``legacy_fortran/mod_tree.f90:735`` and ``:771``) become a ``Genome``
protocol with two tiers: a documentation-grade per-branch API and a
batched API that is the actual hot path.

.. code-block:: python

   from typing import Protocol
   import numpy as np

   class Genome(Protocol):
       # Per-branch API — documentation surface, used by ConstantGenome
       # and by tests. Default implementations of the batched methods
       # fall back to looping these.
       def safety(self, nb_leaves: int, max_stress: float) -> float: ...
       def allocation(
           self, nb_leaves: int, vol_relative: float
       ) -> tuple[float, float, float]: ...
       def mutate(self) -> "Genome": ...

       # Batched API — called once per generation, not once per branch.
       # The C++ orchestrator hands the genome NumPy arrays gathered
       # across all branches in one phase; the genome returns arrays.
       def safety_batch(
           self,
           nb_leaves: np.ndarray,
           max_stress: np.ndarray,
       ) -> np.ndarray: ...
       def allocation_batch(
           self,
           nb_leaves: np.ndarray,
           vol_relative: np.ndarray,
       ) -> tuple[np.ndarray, np.ndarray, np.ndarray]: ...

Two built-ins ship in ``mechatree.genome``:

* ``NeuralGenome(nn_branch, nn_reserve)`` — direct port of the Fortran
  3-layer tanh networks (10 and 18 weights). **The weights are passed
  into C++ at construction and the per-branch NN evaluation runs inline
  in C++**, with no Python callback in the loop. This is the default
  path and the one that matches Fortran performance.
* ``ConstantGenome(safety_val, p_seeds, p_leaves, phototropism)`` — for
  tests and minimum-viable demos. Also resolves to a C++ struct, also
  bypasses Python.

User-supplied Python genomes opt into the slow path. The C++ side
detects a non-built-in ``Genome`` and falls back to the batched
callback: one Python call per phase per generation, never per branch.
At three thousand generations and a tens-of-microseconds-per-call
NumPy forward pass, that adds tens of milliseconds — negligible.
Naive per-branch Python dispatch is explicitly **not supported**
(documented as such).

``mutate`` is part of the protocol from day one so the deferred
evolution work (out of scope here) is a free composition rather than a
refactor.


Public API for building a tree
==============================

The "functionalities" of a tree (mechanics rule, light response, growth
law, pruning rule, genome) are bundled into a frozen ``TreeConfig``
dataclass. The single-tree entry point reads from it:

.. code-block:: python

   from dataclasses import dataclass
   from typing import Callable, Optional
   from mechatree._core import PyTree

   @dataclass(frozen=True)
   class TreeConfig:
       # scalar parameters sourced from YAML
       twig_length: float
       twig_diameter: float
       leaf_surface: float
       cauchy: float
       volume_ratio_leaf: float
       maintenance_h: float
       max_branches: int

       # functionality plug-ins; each defaults to a Fortran-faithful
       # built-in from mechatree.mechanics / .light / .growth / .pruning
       genome: Genome
       mechanics: Callable[..., None]
       light: Callable[..., None]
       growth_law: Callable[..., None]
       pruning_rule: Callable[..., None]

   def grow_tree(
       config: TreeConfig,
       *,
       n_generations: int,
       seed: Optional[int] = None,
   ) -> PyTree:
       """Single-tree main loop. Mirrors legacy_fortran/tree.f90."""

A minimal run:

.. code-block:: python

   from mechatree.config import load_config
   from mechatree.simulate import grow_tree
   from mechatree.genome import NeuralGenome

   cfg = load_config("examples/forest.yaml").to_tree_config(
       genome=NeuralGenome.from_archive("BestShort"),
   )
   tree = grow_tree(cfg, n_generations=300, seed=42)

Swapping any functionality is a single keyword argument — the
composition mechanism. There is no class hierarchy and no plug-in
registry; each functionality is just a callable (or a ``Genome``).


Free functions: module map
==========================

All mechanics, growth, pruning, light and tree-shape routines are **free
functions in submodules**, not methods on ``PyTree``. Each takes a
``PyTree`` (or a ``Leaves``) as its first argument and mutates it in
place. None of them touches the property map; they read and write the
typed mechanics fields.

The Python module structure is a Pythonic facade — the implementations
live in ``src/mechatree/_core/`` and are reached through Cython entry
points.

.. list-table::
   :header-rows: 1
   :widths: 25 45 30

   * - Module
     - Function
     - Fortran reference
   * - ``mechatree.mechanics``
     - ``wind_force(branch, wind) -> (force, moment)``
     - ``mod_tree.f90:613``
   * - ``mechatree.mechanics``
     - ``calculate_stresses(tree, leaf_surface, cauchy) -> None``
     - ``mod_tree.f90:642``
   * - ``mechatree.growth``
     - ``requested_growth(tree, genome, maintenance_h) -> None``
     - ``mod_tree.f90:704``
   * - ``mechatree.growth``
     - ``secondary_growth(tree, volume_per_leaf) -> None``
     - ``mod_tree.f90:817``
   * - ``mechatree.growth``
     - ``primary_growth(tree, genome, twig_length, twig_diameter, generation) -> int``
     - ``mod_tree.f90:425``
   * - ``mechatree.pruning``
     - ``prune(tree, wind, leaf_surface, cauchy) -> int``
     - ``mod_tree.f90:850``
   * - ``mechatree.light``
     - ``intercept(leaves, sun) -> None``
     - ``mod_tree.f90:219``
   * - ``mechatree.light``
     - ``aggregate_onto_trees(leaves, trees) -> None``
     - ``mod_tree.f90:273``
   * - ``mechatree.tree_ops``
     - ``extract_leaves(tree | trees) -> Leaves``
     - ``leaves_extract``
   * - ``mechatree.tree_ops``
     - ``reorder(tree) -> None``
     - ``order_tree``


C++ / Python boundary
=====================

What lives where, and why:

.. list-table::
   :header-rows: 1
   :widths: 18 52 30

   * - Layer
     - Responsibilities
     - Rationale
   * - C++
     - Topology, mechanics struct fields, load-propagation walks,
       per-direction stress accumulation, shadow rasterisation for light,
       and the default genome NNs (weights stored as C++ arrays,
       evaluated inline). All per-branch inner loops.
     - Matches the Fortran inner-loop budget.
   * - Cython
     - Typed accessors on ``PyTree`` for each mechanics field, with
       NumPy-style batch read/write. Thin wrappers that release the GIL
       across each C++ phase call. The batched genome callback
       trampoline (gather NumPy arrays from C++, call Python once,
       scatter results back).
     - One marshalling layer per *phase*, never per branch.
   * - Python
     - YAML loading and validation. Per-generation step sequencing.
       Leaf extraction (which calls into C++). Saving, plotting,
       examples.
     - Glue code; correctness over speed.

The load-propagation walks live in C++. The per-generation step
ordering — the ``for gen in range(...)`` loop — lives in Python. The
Python orchestrator calls into C++ once per phase per generation; with
at most ten phases that is microseconds of dispatch overhead per
generation, negligible compared with the millisecond-scale work inside
each phase.

A short note on memory: ``tests/test_pytree_memory.py`` is currently
marked ``xfail strict`` (the C++ ``~Tree()`` does not yet delete its
branches). The mechanics-fields work in Step 9 touches that destructor;
the leak fix happens alongside.


Performance budget
==================

The design commits to an explicit target so Step 9 has a pass/fail
criterion:

   A single tree, 3000 generations, default configuration (built-in
   ``NeuralGenome``, built-in mechanics, light and pruning) runs
   **within 2× of the Fortran reference** on the same machine.
   User-supplied Python genomes through the batched API: **within 3×**.
   Naive per-branch Python genomes are not supported.

Estimated breakdown on a steady-state tree of ~5000 branches over 3000
generations:

.. list-table::
   :header-rows: 1
   :widths: 36 32 16 16

   * - Phase
     - Per-generation work
     - Cost / gen
     - Total
   * - Stress (4 wind dirs × tree walk)
     - 4 × N C++ float ops
     - ~50 μs
     - 150 ms
   * - Light interception (32 dirs × rasterise)
     - 32 × L grid ops, C++
     - ~1 ms
     - 3 s
   * - Default genome NN evaluation
     - 2 × N forward passes, C++ inline
     - ~30 μs
     - 90 ms
   * - Secondary / primary growth + pruning
     - N tree walks, C++
     - ~100 μs
     - 300 ms
   * - Python orchestrator overhead
     - ~10 C++ calls @ 1 μs
     - ~10 μs
     - 30 ms
   * - **Total**
     -
     -
     - **~4 s**

The Fortran reference is roughly 2–3 s at this scale, so a 2× ceiling
is plausible but not free. Step 9's ``benchmarks/simulation_bench.md``
will gate this with side-by-side Fortran and Python numbers on a fixed
seed; ``benchmarks/baseline.md`` (which already exists for the topology
benchmarks) sets the precedent.

Items deliberately outside the budget:

* **Property map on Branch** — ~56 bytes per branch, off the hot path.
  At 10 000 branches it does not blow the L3 cache and is never read by
  inner loops.
* **3D vectors as ``std::array<double, 3>``** — stack-allocated and
  vectorisable; no worse than Fortran's ``real(3)``.
* **Struct-of-arrays leaves** — NumPy ops where Python touches them;
  raw pointers via Cython memoryviews where C++ does.
* **YAML loading** — once per run, microseconds.

Risks that are inside the budget and need attention at Step 9:

* **Per-branch Python callbacks** — forbidden by design; only the
  batched genome API is supported. A user who insists pays for it
  knowingly.
* **``extract_leaves`` allocating a fresh ``Leaves`` each generation** —
  ~80 KB allocation. If this shows up in profiling, a ``Leaves.reset(n)``
  buffer-reuse method is the fix.
* **Tree-mutation during pruning** — ``removeBranch`` is O(N) per call
  because of the ``branch_to_index`` shift. 100 prunings per generation
  on a 10 000-branch tree is borderline; a bulk-remove API may be
  needed.


Single-tree orchestration loop
==============================

Python pseudo-code mirroring ``legacy_fortran/tree.f90:175-238``:

.. code-block:: python

   def grow_tree(config, *, n_generations, seed):
       tree = _make_seed_tree(config)
       rng = np.random.default_rng(seed)
       for gen in range(n_generations):
           leaves = extract_leaves(tree)
           light.intercept(leaves, sun=DEFAULT_SUN)
           light.aggregate_onto_trees(leaves, [tree])
           mechanics.calculate_stresses(
               tree, config.leaf_surface, config.cauchy,
           )
           growth.requested_growth(
               tree, config.genome, config.maintenance_h,
           )
           growth.secondary_growth(
               tree, config.volume_ratio_leaf * twig_vol,
           )
           pruning.prune(tree, wind=DEFAULT_WIND)
           tree_ops.reorder(tree)
           growth.primary_growth(
               tree, config.genome,
               config.twig_length, config.twig_diameter,
               generation=gen,
           )
           tree_ops.reorder(tree)
       return tree

One deliberate departure from Fortran: the Fortran loop runs a second
``light`` pass at the end of each iteration on the new leaves. The
Python loop omits it because the **next** iteration re-extracts leaves
and runs light again at the top. This is a behavioural change, not a
bug; Step 11's regression test pins it down.


Forest container (sketch, for Step 12)
======================================

Just enough to show the design is internally consistent. The forest is
a lightweight container around a list of trees and a leaf concatenation
step; cross-tree light competition falls out of running ``intercept``
on the union of leaves.

.. code-block:: python

   @dataclass
   class Forest:
       trees: list[PyTree]
       positions: np.ndarray  # (n, 2), ground positions
       ages: np.ndarray       # (n,)
       config: ForestConfig

   def step(forest):
       leaves = extract_leaves(forest.trees)  # concatenated
       light.intercept(leaves, sun=DEFAULT_SUN)
       light.aggregate_onto_trees(leaves, forest.trees)
       for t in forest.trees:
           mechanics.calculate_stresses(t, ...)
           growth.requested_growth(t, t.genome, ...)
           growth.secondary_growth(t, ...)
           pruning.prune(t, ...)
           tree_ops.reorder(t)
       _kill_old_or_small(forest)
       _birth_seedlings(forest)


Review checklist
================

When this document is reviewed, every item in the CLAUDE.md Step 8
list should be one click away:

1. ``Leaf`` and ``Branch`` data shapes → `Data shapes`_.
2. YAML configuration schema → `YAML configuration schema`_.
3. Public API for constructing trees with chosen functionalities →
   `Public API for building a tree`_.
4. Callback signatures for genome decisions → `Genome callbacks`_.
5. Interface between the C++ core and the Python orchestrators →
   `C++ / Python boundary`_ and `Performance budget`_.
