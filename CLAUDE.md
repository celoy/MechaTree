# MechaTree — guide for Claude sessions

## What this project is

MechaTree is a Python (with Cython + C++) port of a Fortran90 simulator of **tree growth in response to wind and light**. The model captures self-pruning, allometric scaling, and self-similar tree architectures that emerge from the mechanical and photosynthetic constraints. It accompanies:

> Eloy, C., Fournier, M., Lacointe, A. et al. *Wind loads and competition for light sculpt trees into self-similar structures.* **Nat Commun 8, 1014 (2017).** https://doi.org/10.1038/s41467-017-00995-6

## Lineage

1. **Fortran90 reference code** — written by the user, simulations published in Nat Commun (2017). Sources preserved in `legacy_fortran/`.
2. **Intern's Python port (2017)** — Python + Cython + C++ library that handles tree architectures. Lightly touched up Feb 2025. Preserved verbatim in `archive/`.
3. **Modernization (in progress)** — the work happening now. Canonical code lives in `src/mechatree/`. Step 1 (this commit) is scaffold-only; the actual port from `archive/` into `src/mechatree/` happens in later steps.

## Where code lives

| Path | Status | Rule |
| --- | --- | --- |
| `src/mechatree/` | **Canonical**, currently a stub | Write new code here |
| `tests/` | Active | Add tests here |
| `archive/PyTreeLib/` | Reference only | Do not edit. Do not import from. |
| `archive/CommentedLibrary_and_Doc/` | Reference only | Do not edit. Do not import from. |
| `legacy_fortran/` | Reference only | Do not edit. Not built or run by Python. |

## Build & test commands

The project uses **uv** for environment management (replaces pyenv + pip + venv).

```bash
# One-time setup
brew install uv                          # if not already installed
uv venv --python 3.12                    # creates .venv/
uv pip install -e ".[dev]"               # installs mechatree + dev deps
uv run pre-commit install                # install git hooks

# Day-to-day
uv run pytest                            # run tests
uv run ruff check .                      # lint
uv run ruff format .                     # format
```

## Conventions

- **Python**: 3.10+. Type hints encouraged. Ruff-formatted (line length 100).
- **Naming**: scientific naming may mirror the Fortran originals where it aids cross-reference with the paper.
- **Provenance**: keep `archive/` and `legacy_fortran/` intact — they are intellectual history.

## What NOT to do

- Don't edit anything under `archive/` or `legacy_fortran/`.
- Don't import from `archive/` in `src/mechatree/` or `tests/`.
- Don't re-introduce the 15 GB of simulation outputs / zips / MATLAB scripts that lived in the old `FORESTArticlePNAS/` folder — they now live at `/Users/Ch/Documents/Python/Eloy2017_NatComm_archive/` (outside the repo, local-only).
- Don't be misled by the legacy folder name `FORESTArticlePNAS` — the paper landed at **Nature Communications**, not PNAS. The "PNAS" name reflects an earlier submission target.

## Step-by-step modernization plan

The user prefers staged steps. **Step 1 (done): repo scaffold.** Future steps (in rough order):

- **Step 2**: Choose & wire the C++/Cython build (probably stays setuptools+Cython, or migrate to scikit-build-core). Port the core `branch.h/cpp`, `tree.h/cpp`, `pytree.pyx`, `cytree.pxd` into `src/mechatree/`. Get `uv pip install -e .` building the extension.
- **Step 3**: Port auxiliary Cython modules (`mod_3Dplot`, `mod_dist`, `mod_plot`) and the example scripts (`random_growth.py`, etc.) into `src/mechatree/` as proper subpackages.
- **Step 4**: Real tests against the C++ extension; possibly numerical regression tests against Fortran outputs.
- **Step 5**: Migrate Sphinx docs from `archive/CommentedLibrary_and_Doc/Doc/` into `docs/` and update.
- **Step 6**: Real CI matrix (multi-OS, cibuildwheel for wheels).
- **Step 7**: First PyPI release.
