#!/usr/bin/env python3
# -*- coding: utf-8 -*

from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

"""
Compilation of the library core
"""

"""
The name of the module to import is 'pytree'. The core of the library is written
in the files "pytree.pyx", "branch.cpp" and "tree.cpp". You might want to add
other C++ modules you coded yourself. We don't advice you to do this unless you
are comfortable with C++ language and Cython wrapping. If you want to add new
modules in a simpler way you can code them in cython ".pyx" files and compile
them with the 'setup_mods.py' file.
"""

ext_modules = [
    Extension(
        "pytree",
        sources=["pytree.pyx", "branch.cpp", "tree.cpp"],
        #extra_compile_args=["-fopenmp","-fPIC","-g"],
        #extra_link_args=["-fopenmp"],
        language="c++",
    )
]

for e in ext_modules:
    e.pyrex_directives = {"boundscheck": False}

setup(
    name ='Tree_library',
    ext_modules = cythonize(ext_modules),
    include_dirs = [numpy.get_include()]
)
