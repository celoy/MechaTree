===========================================
Customize the library: add external modules
===========================================

Because of its generic design, the library is meant to support external
modules containing the specificities of different kinds of growth. For
example, to simulate the growth of a tree under mechanical constraints,
you may want branches to carry already-defined properties such as their
deformation or stress, and an arsenal of adapted functions for computing
those properties. Accessing these features from external modules keeps the
core simulation programs lean.

Two such modules ship today: a plotting subpackage and a geometry
subpackage that implements a 2D self-avoidance test. They were ported
from Diego Bengochea Paz's 2017 Cython modules to pure Python.


``mechatree.plotting``
======================

.. automodule:: mechatree.plotting._2d
   :members:

.. automodule:: mechatree.plotting._3d
   :members:


``mechatree.geometry``
======================

.. automodule:: mechatree.geometry.distance
   :members:


Adding your own subpackage
==========================

The recommended approach is to add a new subdirectory under
``src/mechatree/`` containing your Python (or Cython) sources, plus an
``__init__.py`` that re-exports the public API. For Cython sources,
declare an additional ``Extension`` in ``setup.py`` and rebuild with
``uv pip install -e ".[dev]"``.
