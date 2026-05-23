# legacy_fortran/

Reference-only Fortran90 sources from the original tree-growth simulator:

> Eloy, C., Fournier, M., Lacointe, A. et al. *Wind loads and competition for light sculpt trees into self-similar structures.* **Nat Commun 8, 1014 (2017).** https://doi.org/10.1038/s41467-017-00995-6

These files are kept for **provenance and cross-reference** with the Python port under `src/mechatree/`. They are **not** built, run, or imported by the modern package.

## Contents

- `Forest.f90`, `EvoluAlgo.f90`, `tree.f90`, `stat.f90` — top-level programs
- `mod_tree.f90`, `mod_evolu.f90`, `mod_tools.f90`, `precis_mod.f90` — modules
- `Makefile` — original gfortran build (untested in current toolchains)
- `Forest.ini`, `Evolution.ini` — simulation parameter files

## Where the rest of the original archive lives

Large simulation outputs, MATLAB analysis scripts, compiled binaries, datasets, and figures from the paper were moved out of this repository to:

`/Users/Ch/Documents/Python/Eloy2017_NatComm_archive/`

(That archive is local-only and not committed anywhere.)
