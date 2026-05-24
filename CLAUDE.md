# MechaTree — guide for Claude sessions

## What this project is

MechaTree is a Python (with Cython + C++) port of a Fortran90 simulator of **tree growth in response to wind and light**. The model captures self-pruning, allometric scaling, and self-similar tree architectures that emerge from the mechanical and photosynthetic constraints. It accompanies:

> Eloy, C., Fournier, M., Lacointe, A. et al. *Wind loads and competition for light sculpt trees into self-similar structures.* **Nat Commun 8, 1014 (2017).** https://doi.org/10.1038/s41467-017-00995-6

## Lineage

1. **Fortran90 reference code** — written by the user, simulations published in Nat Commun (2017). Sources preserved in `legacy_fortran/`.
2. **Diego Bengochea Paz's Python port (2017)** — Python + Cython + C++ library that handles tree architectures, written by Diego Bengochea Paz (ORCID [0000-0002-0835-3981](https://orcid.org/0000-0002-0835-3981)) during a 2017 internship. Lightly touched up Feb 2025. Preserved verbatim in `archive/`.
3. **Modernization (in progress)** — canonical code lives in `src/mechatree/`. Steps 1–6 complete: scaffold, C++/Cython port, auxiliary modules, tests, Sphinx docs, CI matrix + cibuildwheel. The simulator port (Steps 8–12) is also done: mechanics + growth in the C++ core, a light-interception subpackage, single-tree orchestrator, and a forest container. Genome configurability is in: Step 14 (constants in YAML) and Step 16 (`NeuralSafety` / `NeuralAllocation`, loadable from a champion JSON written by the verification script) are landed. The plotting layer has migrated from matplotlib to plotly (no matplotlib runtime dependency). Remaining work is the PyPI release (Step 7), forester-facing notebook tutorials (Step 13), and the optional SymPy-callable genome (Step 15).

## Where code lives

| Path | Status | Rule |
| --- | --- | --- |
| `src/mechatree/_core/` | **Canonical** — compiled tree-architecture core (`PyTree`, branches, Strahler/Horton) + genome models (`SafetyModel` / `AllocationModel` and the `Constant*` / `Neural*` subclasses in `genome.h`) | Hot path; C++/Cython |
| `src/mechatree/geometry/` | **Canonical** — geometric helpers (distance, etc.) | Pure Python |
| `src/mechatree/plotting/` | **Canonical** — 2D / 3D plotly renderers | Pure Python |
| `examples/` | Active | Recipe-style scripts + commented YAMLs for users |
| `scripts/` | Active | Off-recipe tooling (e.g. `strategies_single_tree.py` — verifies an evolved genome against the 2017 paper and writes a champion JSON) |
| `data/` | Active | Small reference datasets shared by scripts + tests + examples (e.g. `S3_champions.json`, the per-species champion weights extracted from the Fortran tournament dump) |
| `tests/` | Active | Add tests here |
| `benchmarks/` | Active | Micro-benchmarks + `baseline.md` |
| `docs/` | Active | Sphinx sources |
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

## Who this is for

The intended end user is a **scientist working in forestry**, or someone adopting MechaTree for a **course or research project**. They are not assumed to know Fortran or C++.

API implications:
- Pythonic identifiers in the public surface; cross-reference Fortran names in docstrings/comments only.
- Errors over silent sentinels (raise `IndexError`/`ValueError`, don't return `-1`).
- Examples in `examples/` should read like recipes — runnable end-to-end with one command, minimal boilerplate.
- Config in YAML/JSON, not hand-edited Python.

## Design principles for the simulator port

These principles came out of the user's roadmap notes. They guide all Step-8+ work; revisit before re-litigating.

- **Composition over inheritance.** A tree is configured with a list of "functionalities" (mechanics, light response, growth law, pruning rule) chosen at construction. Avoid a single giant `PyTree` subclass; prefer pluggable hooks.
- **Mechanics is central.** Load propagation (wind force → bending moments → stresses up the tree graph) lives in the compiled core. Hot path; do not push to Python.
- **Light interception is detached.** Separate module. Input: a list of `Leaf` records (location, diameter, transparency). Output: light absorbed per leaf. Decoupled from `PyTree` so it can apply to a single tree, a forest, or any leaf cloud.
- **Genome as callback.** Fortran's `neural_branch` / `neural_reserve` decision functions become user-supplied Python callables (or a small set of named built-ins). Trees do not own a genome class; they accept decision callables.
- **Evolution is external.** The evolutionary algorithm (`EvoluAlgo.f90`, `mod_evolu.f90`) is **not** part of the core library. Possibly a separate package or notebook later.
- **Forest is a container.** Many trees, spatial layout, lifecycle (birth → growth → death). Cross-tree light competition handled by the detached light module operating on the union of leaves.
- **Config in YAML/JSON.** Replace Fortran `.ini` files. Schema validated; an example file ships with the package.

## What NOT to do

- Don't edit anything under `archive/` or `legacy_fortran/`.
- Don't import from `archive/` in `src/mechatree/` or `tests/`.
- Don't re-introduce the 15 GB of simulation outputs / zips / MATLAB scripts that lived in `/Users/Ch/Documents/Python/Eloy2017_NatComm_archive/` (outside the repo, local-only).
- Don't create git commits on `main` or push to `main`. The user manages all commits on the main branch themselves — stage changes and surface diffs, but leave `git commit`/`git push` to them. (This overrides Claude Code's default willingness to commit when explicitly asked, unless the user reconfirms in the current session.)

## Step-by-step modernization plan

The user prefers staged steps.

### Done

- **Step 1** — Repo scaffold (`pyproject.toml`, `src/mechatree/` package, pre-commit, initial CI).
- **Step 2** — Wired the setuptools + Cython build; ported `branch.{h,cpp}`, `tree.{h,cpp}`, `_core.pyx`, `cytree.pxd` into `src/mechatree/_core/`. `uv pip install -e .` builds the extension.
- **Step 3** — Auxiliary modules: `src/mechatree/geometry/` (distance), `src/mechatree/plotting/` (2D / 3D), and recipe-style scripts under `examples/`.
- **Step 4** — Real tests against the C++ extension under `tests/`, including a memory regression test for the C++ leak fix.
- **Step 5** — Sphinx docs migrated into `docs/`.
- **Step 6** — Multi-OS CI matrix (Linux / macOS arm64 / macOS x86_64 / Windows × Python 3.10–3.13) + cibuildwheel pipeline producing wheels and sdist as artifacts.
- **Step 8 — Design doc & API sketch.** `docs/design.rst` captures the `Leaf` / `Branch` shapes, YAML schema, public API and callback signatures, and the C++/Python boundary.
- **Step 9 — Mechanics + growth in the core.** `legacy_fortran/mod_tree.f90` ported into `src/mechatree/_core/`: `wind_force`, `calculate_stresses` ([mechanics.cpp](src/mechatree/_core/mechanics.cpp)), `requested_growth`, `primary_growth`, `secondary_growth` ([growth.cpp](src/mechatree/_core/growth.cpp)), `pruning`, `cut_branch` ([pruning.cpp](src/mechatree/_core/pruning.cpp)). Genome callbacks (`SafetyModel`, `AllocationModel`) plug in via the C++ vtable from the Cython boundary.
- **Step 10 — Light interception module.** `src/mechatree/light/` ports `light_interception` / `light_on_trees`; pure Python operating on a `Leaves` collection + `Sun` model. Decoupled from `PyTree`.
- **Step 11 — Single-tree simulator.** `mechatree.simulate.grow_tree(config)` runs the `{light → stresses → growth → pruning → seeds/leaves → ordering}` loop. Example in `examples/single_tree.py`.
- **Step 12 — Forest container.** `mechatree.forest.Forest` owns N trees with spatial layout, birth/death lifecycle, cross-tree light competition through the Step-10 module on the concatenated leaf list.
- **Step 14 — Genome constants in YAML.** Scalars `safety`, `p_seeds`, `p_leaves`, `phototropism` live under a `genome:` block in YAML, populated through `GenomeConfig` (`src/mechatree/config.py`). Default `safety = 3.0` per the Fortran neural-net finding. `grow_tree` / `Forest` build `ConstantSafety` / `ConstantAllocation` from those values when no explicit model is passed.
- **Step 16 — `NeuralSafety` / `NeuralAllocation`.** The Fortran 3-layer tanh networks (`mod_tree.f90:735 neural_branch`, `:771 neural_reserve`) are ported as header-only C++ classes in [src/mechatree/_core/genome.h](src/mechatree/_core/genome.h), wrapped at the Cython boundary as `PyNeuralSafety` / `PyNeuralAllocation`, and re-exported from `mechatree.genome`. Champion weight vectors are loaded with `mechatree.genome.load_champion(path, species_id)` from a JSON file written by [scripts/strategies_single_tree.py](scripts/strategies_single_tree.py); the reference dataset is [data/S3_champions.json](data/S3_champions.json). YAML configs can select a neural genome via `genome.neural_from: {path, species_id}` (path resolved relative to the YAML file). The C++ port is verified against the Python reference to `atol=1e-12` at every grid point.

### Up next

- **Step 7 — First PyPI release.** Tag `v0.1.0`; release the current "tree architecture + simulator + neural genome" library. CHANGELOG should note scope: ships `PyTree` topology + geometry + plotting + mechanics + light + single-tree + forest + `Constant*` / `Neural*` genome models. Evolution is deferred to a separate package.
- **Step 13 — Forester-facing examples & tutorials.** Annotated Jupyter notebooks under `notebooks/`, all rendering through **plotly inline** (the library has no matplotlib dependency). Each notebook walks through the science, the YAML config, and the code, with prose between cells. Concrete plan:
  - `notebooks/01_grow_one_tree.ipynb` — load `examples/forest.yaml`, call `grow_tree`, display the Strahler-coloured 3D tree with `plot_tree_3d` (with leaves on/off toggles). Companion to [examples/grow_one_tree.py](examples/grow_one_tree.py).
  - `notebooks/02_forest_under_wind.ipynb` — drive a `Forest`, animate population/biomass curves alongside the top-down view (`plot_forest_topdown`), show self-thinning emerging (`plot_self_thinning`). Companion to [examples/grow_a_forest.py](examples/grow_a_forest.py).
  - `notebooks/03_neural_genome.ipynb` — load an S3 champion via `mechatree.genome.load_champion` (or YAML `genome.neural_from`), compare a tree grown with the constant default genome against the species-0 and species-1 champions side-by-side. Uses `data/S3_champions.json`.
  - `notebooks/04_custom_growth_law.ipynb` — plug in user-supplied `wind_fn` and `sun` callables, vary `ConstantSafety` / `ConstantAllocation`, observe how trunk diameter, canopy depth and pruning rate respond. Companion to [examples/custom_simulation.py](examples/custom_simulation.py).
  - `notebooks/05_strahler_diagnostics.ipynb` — grow a mature tree, render `plot_strahler_diagnostics` (4-panel: branch count, length, area per Strahler order + Leonardo's-rule histogram), discuss the self-similarity result from the paper. Companion to [examples/plot_strahler.py](examples/plot_strahler.py).

  Ship a commented example YAML alongside the notebooks; refresh `docs/userguide.rst` to point at them.
- **Step 15 — Callable genome from YAML via SymPy.** Allow `genome.safety` (and allocation fields) to be either a scalar *or* a string expression like `"3 * tanh(max_stress)"`. Parse with `sympy.sympify`, compile with `sympy.lambdify`, wrap as `SymPySafety` / `SymPyAllocation` subclasses of `SafetyModel` / `AllocationModel`. Inputs: the same `(nb_leaves, max_stress)` / `(nb_leaves, vol_relative)` pairs the Fortran NNs took. Validate against an allow-list of free symbols. Useful as a research tool sitting between `Constant*` and `Neural*`.

### Out of scope (do not re-litigate unless explicitly asked)

- **Evolution** (`legacy_fortran/EvoluAlgo.f90`, `legacy_fortran/mod_evolu.f90`) — deferred, possibly a separate package later. Per the design principle above, evolution is external to the core library.
