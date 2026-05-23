===============
Getting started
===============


Overview of the library
=======================

MechaTree's core consists of two C++ classes: the **Branch** class and the
**Tree** class. Each instance of the Branch class is characterized by a set
of family relations with other instances of the class (pointers to parent
and children) and a set of user-defined properties of type **double**. An
instance of the Tree class is a container (a C++ ``std::vector``) of
pointers to Branch objects.

The Tree class is wrapped by the **PyTree** class using Cython. This means
that ``PyTree`` is a Python version of the C++ Tree class. Many of the
methods of the Tree class have an equivalent in the PyTree class, which
calls the C++ implementation under the hood.

A user writes a Python script where ``PyTree`` objects are initialized and
manipulated via the wrapped methods. These methods act by calling the
equivalent methods of the C++ Tree class. Hence, from a Python script you
interact with the C++ core of the library.


Installing MechaTree
====================

MechaTree requires Python 3.10+ and a C++17-capable compiler (recent gcc,
clang, or MSVC). On macOS, the system clang shipped with Xcode is sufficient.
On Linux, ``build-essential`` provides everything needed.

We recommend `uv <https://docs.astral.sh/uv/>`_ for environment management:

.. code-block:: bash

   # one-time
   brew install uv

   # in the repo root
   uv venv --python 3.12
   uv pip install -e ".[dev,docs]"

This compiles the Cython/C++ extension (``mechatree._core._core``) into the
project's ``.venv``. To verify the install:

.. code-block:: python

   from mechatree import PyTree
   t = PyTree({"length": 1.0, "radius": 0.1})
   print(t.get_number_of_branches())   # -> 1
