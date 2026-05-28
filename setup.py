import os
import subprocess
import sys
import tempfile
from pathlib import Path

from Cython.Build import cythonize
from setuptools import setup
from setuptools.extension import Extension

CORE = Path("src/mechatree/_core")

if sys.platform == "win32":
    cxx_flags = ["/std:c++17", "/O2"]
else:
    cxx_flags = ["-std=c++17", "-O3", "-fvisibility=hidden"]


def _detect_openmp() -> tuple[list[str], list[str]]:
    """Return (extra_compile_args, extra_link_args) to enable OpenMP, or
    empty lists when not available.

    The C++ kernel's ``#pragma omp`` directives are a no-op without these
    flags, so silent fallback gives a correct serial build.

    Per platform:

    * **Windows / MSVC** — uses ``/openmp:llvm`` for the LLVM runtime, which
      supports the ``parallel for`` directive shape used in ``light.cpp``.
      Available in Visual Studio 2019 16.10+ (i.e. every supported CI
      image).
    * **Linux / GCC** — ``-fopenmp`` works out of the box; manylinux wheels
      built via cibuildwheel will pick this up automatically.
    * **macOS** — Apple Clang doesn't ship OpenMP by default; libomp from
      Homebrew is the standard route. We probe the brew prefix and fall
      back to serial when libomp is not installed.

    Override with ``MECHATREE_NO_OPENMP=1`` to force a serial build (useful
    for reproducible-bytes-across-CPU debugging).
    """
    if os.environ.get("MECHATREE_NO_OPENMP", "").strip() not in ("", "0", "false"):
        return [], []

    if sys.platform == "win32":
        return ["/openmp:llvm"], []

    if sys.platform == "darwin":
        # Apple Clang needs libomp explicitly. brew --prefix returns the
        # cellar root for ARM (/opt/homebrew) or Intel (/usr/local).
        try:
            prefix = (
                subprocess.check_output(["brew", "--prefix", "libomp"], stderr=subprocess.DEVNULL)
                .decode()
                .strip()
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            return [], []
        inc = Path(prefix) / "include"
        lib = Path(prefix) / "lib"
        if not (inc / "omp.h").exists():
            return [], []
        return (
            ["-Xpreprocessor", "-fopenmp", f"-I{inc}"],
            [f"-L{lib}", "-lomp", f"-Wl,-rpath,{lib}"],
        )

    # Linux / other Unix-ish. Probe by trying to compile a tiny program.
    return _try_gcc_openmp()


def _try_gcc_openmp() -> tuple[list[str], list[str]]:
    src = "#include <omp.h>\nint main(){return omp_get_max_threads();}"
    with tempfile.TemporaryDirectory() as d:
        src_path = Path(d) / "omp_test.cpp"
        src_path.write_text(src)
        out_path = Path(d) / "omp_test"
        try:
            subprocess.check_call(
                [os.environ.get("CXX", "c++"), "-fopenmp", str(src_path), "-o", str(out_path)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            return [], []
    return ["-fopenmp"], ["-fopenmp"]


_omp_compile, _omp_link = _detect_openmp()
if _omp_compile:
    print(f"[mechatree.setup] OpenMP enabled: compile={_omp_compile} link={_omp_link}")
else:
    print("[mechatree.setup] OpenMP not detected — `light_intercept` will run serial.")

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
            str(CORE / "light.cpp"),
            str(CORE / "momentum.cpp"),
        ],
        include_dirs=[str(CORE)],
        language="c++",
        extra_compile_args=cxx_flags + _omp_compile,
        extra_link_args=_omp_link,
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
