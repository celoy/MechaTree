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

.. note::

   The intern's 2017 port shipped a small set of optional Cython modules
   (``mod_plot``, ``mod_dist``, ``mod_3Dplot``) and example scripts for
   self-avoiding growth and sap transport. These are not yet ported into
   the modern ``mechatree`` package — they remain available under
   ``archive/`` in the repository for reference.

   The mechanism for declaring and shipping user-supplied Cython modules
   will be revisited as part of a later release. The recommended approach
   today is to add your own Cython source to ``src/mechatree/_core/`` (or
   a sibling subpackage), declare it as an additional ``Extension`` in
   ``setup.py``, and rebuild with ``uv pip install -e ".[dev]"``.
