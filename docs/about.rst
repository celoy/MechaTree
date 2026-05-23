==================
About this library
==================

MechaTree aims to be a useful tool for the simulation of theoretical models
describing the growth of branching structures. The core of the library is
written in C++ but the user interacts with it through Python. The interface
between C++ and Python is provided by Cython. Our intention is to provide a
platform where users less familiar with numerical tools can simulate the
growth of ramified structures in a simple and straightforward way, without
losing the characteristic performance of low-level programming languages.

One of the purposes of MechaTree is to be generic enough to allow the
simulation of any kind of growth process involving ramified structures —
such as trees, corals, or neurons. This is achieved by letting the user
decide which properties of the structure are of interest and how they
evolve. The platform makes it straightforward to interact with a structure
(the tree) consisting of a set of elements (the branches) linked by family
relations and carrying user-defined attributes.

If you are interested in contributing to the development of this library,
we advise you to read the whole documentation. If you are just interested
in using the library and do not need to understand its internals in depth,
you can stop after the :doc:`userguide` section.

Acknowledgments
===============

The original Python + Cython + C++ port of the simulator was written in 2017
by **Diego Bengochea Paz** (ORCID `0000-0002-0835-3981
<https://orcid.org/0000-0002-0835-3981>`_) during an internship. That code is
preserved under ``archive/`` for provenance, and the modernization in
``src/mechatree/`` is a direct descendant of his work.
