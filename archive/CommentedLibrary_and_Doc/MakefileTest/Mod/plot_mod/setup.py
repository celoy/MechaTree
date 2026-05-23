#!/usr/bin/env python3
# -*- coding: utf-8 -*

from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

"""
Compilation of the library modules
"""

"""
The name of the module to import is 'pytreemods'. To compile the files
containing the external modules you have to add its names after the 'source'
command:
    'source=[filename1, filename2, ..., filenameN]'
External modules are ".pyx" Cython files. You can create your own if you want.
"""

ext_modules = [
    Extension(
        "mod_plot",
        sources=["mod_plot.pyx"],
    )
]

for e in ext_modules:
    e.pyrex_directives = {"boundscheck": False}

setup(
    name ='PyTree_mods',
    ext_modules = cythonize(ext_modules),
    include_dirs = [numpy.get_include()]
)
