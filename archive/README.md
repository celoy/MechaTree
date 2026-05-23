# archive/

**Reference-only**. These two folders are the intern's 2017 Python+Cython+C++ port of the original Fortran code (see `legacy_fortran/`). They are kept here for provenance — to preserve the lineage of the Python port — and are **not** built, run, or imported by the modern package under `src/mechatree/`.

The two folders are near-duplicates (single typo-fix difference between their core sources). They were both committed because each retains material the other lacks:

- **`PyTreeLib/`** — tidier `Lib/` (core C++/Cython) + `Mod/` (auxiliary Cython modules) split, with example scripts (`random_growth.py`, `sap_transport.py`, `self_avoiding_modules.py`).
- **`CommentedLibrary_and_Doc/`** — same code under `Lib&Modules/`, plus a full Sphinx documentation site under `Doc/` (rst sources + built HTML), plus a `MakefileTest/` test-build layout. The `pytree.pyx` here received a typo fix in Feb 2025.

When the modernization gets to porting actual code into `src/mechatree/`, the `Doc/` rst sources here are the seed for `docs/` and the source files here are the canonical reference.

## Rules

- **Do not edit** anything under `archive/`. Treat it as frozen history.
- **Do not import** from `archive/` in `src/mechatree/` or `tests/`.
- If a build artifact (`.so`, `build/`, `__pycache__/`) appears under `archive/`, it is already gitignored.
