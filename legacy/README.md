# legacy/

**Reference-only.** Everything in this directory is intellectual history from the 2017 Nat Commun paper:

> Eloy, C., Fournier, M., Lacointe, A. et al. *Wind loads and competition for light sculpt trees into self-similar structures.* **Nat Commun 8, 1014 (2017).** https://doi.org/10.1038/s41467-017-00995-6

Kept here for provenance and cross-reference. The modern Python package under [`src/mechatree/`](../src/mechatree/) is a direct descendant of this material but does **not** build, run, or import from anything below.

## Layout

| Path | What it is |
| --- | --- |
| `fortran/` | The original Fortran90 simulator that produced the paper's figures. |
| `matlab/` | The MATLAB analysis + plotting scripts used to make the paper's figures. |
| `pdf/` | The paper itself and its Supplementary Information. |
| `pytree/` | Diego Bengochea Paz's 2017 Python + Cython + C++ port of the Fortran code, with Sphinx documentation. |

## `fortran/`

Fortran90 sources of the original simulator:

- `Forest.f90`, `EvoluAlgo.f90`, `tree.f90`, `stat.f90` — top-level programs
- `mod_tree.f90`, `mod_evolu.f90`, `mod_tools.f90`, `precis_mod.f90` — modules
- `Makefile` — original gfortran build (untested in current toolchains)
- `Forest.ini`, `Evolution.ini` — simulation parameter files

These are the line-numbers cross-referenced from comments + docstrings throughout `src/mechatree/` (e.g. ``legacy/fortran/mod_tree.f90:735`` for `neural_branch`).

## `matlab/`

MATLAB analysis + plotting scripts from the paper workflow. The modern Python equivalents live under [`examples/`](../examples/) / [`src/mechatree/stats.py`](../src/mechatree/stats.py) / [`src/mechatree/plotting/`](../src/mechatree/plotting/); the table below maps the originals to their ports.

| MATLAB | Python equivalent |
| --- | --- |
| `self_thinning.m` | [`examples/plot_self_thinning.py`](../examples/plot_self_thinning.py) |
| `plot_allocation_vs_t.m` | [`examples/plot_allocation.py`](../examples/plot_allocation.py) |
| `plot_stat_single_tree.m`, `Fractal_dim.m`, `plot_area_preservation_1tree.m` | [`examples/plot_strahler.py`](../examples/plot_strahler.py), [`notebooks/06_fractal_dimension.ipynb`](../notebooks/06_fractal_dimension.ipynb), [`src/mechatree/stats.py`](../src/mechatree/stats.py) |
| `plot_single_tree_3D.m`, `plot_single_tree_3D_flat.m`, `plotForest2D.m` | [`src/mechatree/plotting/`](../src/mechatree/plotting/) |
| `strategies_single_tree.m`, `strategies_Forest.m` | [`scripts/strategies_single_tree.py`](../scripts/strategies_single_tree.py), [`src/mechatree/evolution/curate.py`](../src/mechatree/evolution/curate.py) |
| `neural_branch.m`, `neural_reserve.m` | [`src/mechatree/_core/genome.h`](../src/mechatree/_core/genome.h) (C++ port; Python wrappers in [`src/mechatree/genome.py`](../src/mechatree/genome.py)) |
| `ART_Allometry_simplified.m` | analytical reference, no direct port |

## `pdf/`

- `s41467-017-00995-6.pdf` — the paper.
- `SI11.pdf` — the Supplementary Information.

## `pytree/`

Diego Bengochea Paz's 2017 internship port (ORCID [0000-0002-0835-3981](https://orcid.org/0000-0002-0835-3981)). Three subdirectories:

- `Lib&Modules/` — the core C++/Cython sources of the original port.
- `Doc/` — the full Sphinx documentation site for the original port (rst sources + built HTML). Seed material for the modern [`docs/`](../docs/).
- `MakefileTest/` — a test-build layout from the 2017 work.

A second near-duplicate port (`PyTreeLib/`) used to sit alongside this one. It was dropped during the 2026-05-26 legacy consolidation; its three example scripts (`random_growth.py`, `sap_transport.py`, `self_avoiding_modules.py`) had already been modernized as [`examples/random_growth.py`](../examples/random_growth.py), [`examples/sap_transport.py`](../examples/sap_transport.py), [`examples/self_avoiding.py`](../examples/self_avoiding.py).

## Rules

- **Do not edit** anything under `legacy/`. Treat it as frozen history.
- **Do not import** from `legacy/` in `src/mechatree/` or `tests/`.
- Build artifacts (`.so`, `build/`, `__pycache__/`) under `legacy/` are gitignored.
