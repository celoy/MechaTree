import sys
from pathlib import Path

from Cython.Build import cythonize
from setuptools import setup
from setuptools.extension import Extension

CORE = Path("src/mechatree/_core")

if sys.platform == "win32":
    cxx_flags = ["/std:c++17", "/O2"]
else:
    cxx_flags = ["-std=c++17", "-O3", "-fvisibility=hidden"]

extensions = [
    Extension(
        name="mechatree._core._core",
        sources=[
            str(CORE / "_core.pyx"),
            str(CORE / "branch.cpp"),
            str(CORE / "tree.cpp"),
            str(CORE / "mechanics.cpp"),
            str(CORE / "growth.cpp"),
            str(CORE / "pruning.cpp"),
        ],
        include_dirs=[str(CORE)],
        language="c++",
        extra_compile_args=cxx_flags,
    ),
]

setup(
    ext_modules=cythonize(
        extensions,
        include_path=["src"],
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
        },
    ),
)
