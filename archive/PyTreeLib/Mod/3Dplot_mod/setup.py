#!/usr/bin/env python3
# -*- coding: utf-8 -*

from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

"""
Compilation of a library module

The name of the module to import is 'mod_3Dplot'. To compile the file
containing the external module you have to add its name after the 'source'
command.
"""

ext_modules = [
    Extension(
        "mod_3Dplot",
        sources=["mod_3Dplot.pyx"],
    )
]

for e in ext_modules:
    e.pyrex_directives = {"boundscheck": False}

setup(
    name ='PyTree_3Dplot',
    ext_modules = cythonize(ext_modules),
    include_dirs = [numpy.get_include()]
)
