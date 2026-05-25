# MechaTree — guide for Claude sessions

## What this project is

MechaTree is a Python (with Cython + C++) port of a Fortran90 simulator of **tree growth in response to wind and light**. The model captures self-pruning, allometric scaling, and self-similar tree architectures that emerge from the mechanical and photosynthetic constraints. It accompanies:

> Eloy, C., Fournier, M., Lacointe, A. et al. *Wind loads and competition for light sculpt trees into self-similar structures.* **Nat Commun 8, 1014 (2017).** https://doi.org/10.1038/s41467-017-00995-6

## Lineage

1. **Fortran90 reference code** — written by the user, simulations published in Nat Commun (2017). Sources preserved in `legacy_fortran/`.
2. **Diego Bengochea Paz's Python port (2017)** — Python + Cython + C++ library that handles tree architectures, written by Diego Bengochea Paz (ORCID [0000-0002-0835-3981](https://orcid.org/0000-0002-0835-3981)) during a 2017 internship. Lightly touched up Feb 2025. Preserved verbatim in `archive/`.
3. **Modernization (in progress)** — canonical code lives in `src/mechatree/`. See the [Modernization roadmap](#modernization-roadmap) table below for live milestone status.

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
- **Evolution lives behind a clean boundary.** The evolutionary algorithm (`EvoluAlgo.f90`, `mod_evolu.f90`) gets ported as part of Step 21. Whether it ships as `mechatree.evolution` or as a sister package is one of Step 21's open design questions; either way it must import the simulator, never the other way round.
- **Forest is a container.** Many trees, spatial layout, lifecycle (birth → growth → death). Cross-tree light competition handled by the detached light module operating on the union of leaves.
- **Config in YAML/JSON.** Replace Fortran `.ini` files. Schema validated; an example file ships with the package.

## What NOT to do

- Don't edit anything under `archive/` or `legacy_fortran/`.
- Don't import from `archive/` in `src/mechatree/` or `tests/`.
- Don't re-introduce the 15 GB of simulation outputs / zips / MATLAB scripts that lived in `/Users/Ch/Documents/Python/Eloy2017_NatComm_archive/` (outside the repo, local-only).
- Don't create git commits on `main` or push to `main`. The user manages all commits on the main branch themselves — stage changes and surface diffs, but leave `git commit`/`git push` to them. (This overrides Claude Code's default willingness to commit when explicitly asked, unless the user reconfirms in the current session.)

## Modernization roadmap

The user prefers staged steps. The table is the source of truth — flip a row's status as work lands and add a `### Step N notes (closed YYYY-MM-DD)` subsection below.

|    | Step                                | Status     | Summary |
|----|-------------------------------------|------------|---------|
| 1  | Repo scaffold                       | ✅ done    | `pyproject.toml`, `src/mechatree/` package skeleton, pre-commit, initial CI. |
| 2  | C++/Cython port                     | ✅ done    | setuptools + Cython build; `branch.{h,cpp}`, `tree.{h,cpp}`, `_core.pyx`, `cytree.pxd` under [src/mechatree/_core/](src/mechatree/_core/). |
| 3  | Auxiliary modules                   | ✅ done    | [src/mechatree/geometry/](src/mechatree/geometry/) (distance), [src/mechatree/plotting/](src/mechatree/plotting/) (2D/3D), recipe scripts under [examples/](examples/). |
| 4  | Real C++ tests                      | ✅ done    | tests/ exercises the compiled extension, incl. a memory-leak regression. |
| 5  | Sphinx docs                         | ✅ done    | docs/ scaffolded with Furo. |
| 6  | CI matrix + cibuildwheel            | ✅ done    | Linux/macOS arm64/macOS x86_64/Windows × py3.10–3.13 + wheel/sdist artifacts. |
| 7  | Pre-release shake-out (TestPyPI)    | ⬜ pending | `v0.1.0a1` candidate: `uv build` wheel + sdist, smoke-install in a clean venv, dry-run publish to TestPyPI via `workflow_dispatch`. Validates packaging before burning the real PyPI name. |
| 8  | Design doc & API sketch             | ✅ done    | `docs/design.rst` — Leaf/Branch shapes, YAML schema, public API, callback signatures, C++/Python boundary. |
| 9  | Mechanics + growth in core          | ✅ done    | `wind_force`, `calculate_stresses` ([mechanics.cpp](src/mechatree/_core/mechanics.cpp)); growth ([growth.cpp](src/mechatree/_core/growth.cpp)); `pruning` / `cut_branch` ([pruning.cpp](src/mechatree/_core/pruning.cpp)). Genome callbacks plug in via a C++ vtable. |
| 10 | Light interception module           | ✅ done    | [src/mechatree/light/](src/mechatree/light/) — pure-Python `Leaves`/`Sun`, decoupled from `PyTree`. |
| 11 | Single-tree simulator               | ✅ done    | `mechatree.simulate.grow_tree(config)` runs the `{light → stresses → growth → pruning → seeds/leaves → ordering}` loop. |
| 12 | Forest container                    | ✅ done    | `mechatree.forest.Forest` — N trees, spatial layout, birth/death lifecycle, cross-tree light competition. |
| 13 | Forester notebooks                  | ✅ done    | Five plotly-inline tutorials under [notebooks/](notebooks/), each a literate companion to an `examples/` script. `docs/userguide.rst` cross-references them. |
| 14 | Genome constants in YAML            | ✅ done    | `safety`, `p_seeds`, `p_leaves`, `phototropism` under `genome:` block; default `safety = 3.0` per the Fortran NN result. |
| 15 | SymPy callable genome               | ✅ done    | YAML `genome:` scalars accept string expressions; SymPy parses + lambdifies, the C++ growth loop dispatches via new `CallbackSafety` / `CallbackAllocation` vtable subclasses. |
| 16 | Neural genome (`Neural*`)           | ✅ done    | header-only C++ port of the Fortran 3-layer tanh nets ([genome.h](src/mechatree/_core/genome.h)), verified to atol=1e-12. Champion JSON loader. |
| 17 | DendroFlow wind bridge              | ✅ done    | [src/mechatree/wind/dendroflow.py](src/mechatree/wind/dendroflow.py) wraps DendroFlow's `BulkThinningBranchWindModel` as a 3-arg `WindFn`. YAML `wind:` block selects it. M6 of DendroFlow. |
| 18 | Stable PyPI 0.1.0 release           | ⬜ pending | After Step 17 lands. Tag `v0.1.0`, flip the TestPyPI dry-run from Step 7 to a real-PyPI publish job. Mirrors DendroFlow's M7. |
| 19 | Leaf transparency in light model    | ⬜ pending | Wire a single global scalar `light.leaf_transparency` (default reproduces today's binary topmost-wins behaviour) into [`mechatree.light.interception`](src/mechatree/light/interception.py): replace `max(1 - shadow, 0)` with `tau ** shadow` so light at depth `i` in a shadow cell is `tau^i` (1 → 0.5 → 0.25 → … for `tau=0.5`). New field on [`LightConfig`](src/mechatree/config.py), surfaced via YAML `light:`. ~15 LOC + a test that `tau=0` recovers the current binary behaviour and `tau=1` makes leaves fully transparent. |
| 20 | Mesh3d cylindrical tree renderer    | ⬜ pending | Today's [`plot_tree_3d`](src/mechatree/plotting/_mechanics.py) uses pixel-width plotly `Scatter3d` lines bucketed by Strahler order, so individual branch diameters are flattened and the canopy looks pencil-thin compared to the Blender renders in the paper. Add a `style="cylinders"` path that extrudes each branch as a tapered `Mesh3d` cylinder (parent diameter → child diameter), with the existing `style="lines"` kept as the fast default. ~60 LOC + a small cylinder-primitive helper. |
| 21 | Evolution port (`EvoluAlgo` / `mod_evolu`) | ⬜ pending | Port the Fortran evolutionary algorithm so the Nat Comms tournament simulations are reproducible end-to-end in Python: per-generation fitness evaluation via `grow_tree` on a population of `NeuralSafety` / `NeuralAllocation` champions, tournament selection, neural-net weight mutation + crossover, and a CSV/JSON dump matching `S3.dat` so [`scripts/strategies_single_tree.py`](scripts/strategies_single_tree.py) can keep loading champions the same way. Open design questions to settle before coding: separate package vs in-tree `mechatree.evolution` module; whether to parallelise fitness eval via `multiprocessing` or keep it sequential; whether to mirror the Fortran tournament structure verbatim or restructure around a `numpy`-friendly `Population` dataclass. Reverses the "Evolution is external" call in *Design principles for the simulator port*. |

### Step 1 notes (closed 2026-05-23)

Repo scaffold landed in commit `chore: scaffold modern repo structure (Step 1)`. Establishes `pyproject.toml`, `src/mechatree/` layout, pre-commit + ruff config, and an initial GitHub Actions CI workflow.

### Step 2 notes (closed 2026-05-23)

`feat: port Cython/C++ library + migrate Sphinx docs (Step 2)` wired the setuptools + Cython build and ported `branch.{h,cpp}` / `tree.{h,cpp}` / `_core.pyx` / `cytree.pxd` into [src/mechatree/_core/](src/mechatree/_core/). `uv pip install -e .` builds the extension on macOS/Linux/Windows.

### Step 3 notes (closed 2026-05-23)

Pure-Python helpers around the C++ core:
- [src/mechatree/geometry/](src/mechatree/geometry/) — distance/orientation utilities used by Step-10 light + Step-11 growth.
- [src/mechatree/plotting/](src/mechatree/plotting/) — 2D/3D renderers, now plotly-only after the May-2026 migration.
- [examples/](examples/) — recipe-style scripts.

### Step 4 notes (closed 2026-05-23)

Tests exercise the compiled extension directly; a dedicated `test_pytree_memory.py` regression locks in the C++ destructor fix that stopped each `PyTree` from leaking its branch graph on GC.

### Step 5 notes (closed 2026-05-23)

Sphinx scaffolding under [docs/](docs/) using the Furo theme. Real narrative content is still scheduled with Step 13.

### Step 6 notes (closed 2026-05-24)

`.github/workflows/ci.yml` runs lint + pytest across Linux / macOS arm64 / macOS x86_64 / Windows × Python 3.10–3.13. `cibuildwheel` builds platform wheels + an sdist as CI artifacts; Windows-specific bugs were ironed out in `ci bug` / `CI windows bug` commits the same day.

### Step 8 notes (closed 2026-05-23)

`docs/design.rst` captures the `Leaf` / `Branch` shapes, the YAML schema, the public API, callback signatures, and the C++/Python boundary. Drafted before code so Steps 9–12 had a target.

### Step 9 notes (closed 2026-05-23)

Mechanics + growth + pruning lowered into the compiled core. Genome callbacks (`SafetyModel`, `AllocationModel`) plug in via a C++ vtable from the Cython boundary — see `genome.h`. Files:
- `wind_force`, `calculate_stresses` → [src/mechatree/_core/mechanics.cpp](src/mechatree/_core/mechanics.cpp)
- `requested_growth`, `primary_growth`, `secondary_growth` → [src/mechatree/_core/growth.cpp](src/mechatree/_core/growth.cpp)
- `pruning`, `cut_branch` → [src/mechatree/_core/pruning.cpp](src/mechatree/_core/pruning.cpp)

### Step 10 notes (closed 2026-05-23)

`light interception ported` — [src/mechatree/light/](src/mechatree/light/) hosts `light_interception` / `light_on_trees`, pure-Python on a `Leaves` collection + `Sun` model. Decoupling lets the same module score a single tree, a `Forest`, or any leaf cloud (e.g. LiDAR).

### Step 11 notes (closed 2026-05-23)

`single-tree simulator` — `mechatree.simulate.grow_tree(config)` walks the Fortran `{light → stresses → growth → pruning → seeds/leaves → ordering}` loop. The orchestrator is intentionally Python — at most ten C++/NumPy calls per generation. Example at [examples/grow_one_tree.py](examples/grow_one_tree.py).

### Step 12 notes (closed 2026-05-23)

`forest port` — `mechatree.forest.Forest` owns N trees with disk layout, birth/death lifecycle, and cross-tree light competition through the Step-10 module on the concatenated leaf list. Wind is currently shared across all trees in a generation.

### Step 14 notes (closed 2026-05-24)

`default safety updated` — scalars `safety`, `p_seeds`, `p_leaves`, `phototropism` live under a `genome:` block in YAML, populated via [`GenomeConfig`](src/mechatree/config.py). Default `safety = 3.0` reproduces the value the Fortran neural net evolved to. `grow_tree` / `Forest` build `ConstantSafety` / `ConstantAllocation` from those values when no explicit model is passed.

### Step 16 notes (closed 2026-05-24)

`import genome from Nat Comms` — Fortran 3-layer tanh networks (`mod_tree.f90:735 neural_branch`, `:771 neural_reserve`) ported as header-only C++ in [src/mechatree/_core/genome.h](src/mechatree/_core/genome.h), wrapped at the Cython boundary as `PyNeuralSafety` / `PyNeuralAllocation`, and re-exported from `mechatree.genome`. Champion weight vectors load via `mechatree.genome.load_champion(path, species_id)` from a JSON file written by [scripts/strategies_single_tree.py](scripts/strategies_single_tree.py); the reference dataset is [data/S3_champions.json](data/S3_champions.json). YAML configs can select a neural genome via `genome.neural_from: {path, species_id}` (path resolved relative to the YAML file). The C++ port is verified against the Python reference to `atol=1e-12` at every grid point.

### Step 13 notes (closed 2026-05-25)

Five plotly-inline Jupyter tutorials under [notebooks/](notebooks/), each a literate-programming companion to an `examples/` script. Markdown prose between cells walks through the science; code cells run end-to-end against the canonical YAML config. Verified by executing all five via `jupyter nbconvert --execute`; outputs stripped back out via `nbstripout` per the [notebooks/README.md](notebooks/README.md) convention.

- [01_grow_one_tree.ipynb](notebooks/01_grow_one_tree.ipynb) — companion to [examples/grow_one_tree.py](examples/grow_one_tree.py). Walks the per-generation loop and renders the Strahler-coloured 3D canopy + leaves.
- [02_forest_under_wind.ipynb](notebooks/02_forest_under_wind.ipynb) — companion to [examples/grow_a_forest.py](examples/grow_a_forest.py). Population & biomass over time + final top-down stand + the Yoda −3/2 self-thinning curve.
- [03_neural_genome.ipynb](notebooks/03_neural_genome.ipynb) — loads the two S3 champions from [data/S3_champions.json](data/S3_champions.json) and renders them side-by-side against the default constant genome.
- [04_custom_growth_law.ipynb](notebooks/04_custom_growth_law.ipynb) — companion to [examples/custom_simulation.py](examples/custom_simulation.py). Plug-in `wind_fn` / `Sun` / `ConstantSafety` / `ConstantAllocation`.
- [05_strahler_diagnostics.ipynb](notebooks/05_strahler_diagnostics.ipynb) — companion to [examples/plot_strahler.py](examples/plot_strahler.py). Horton–Strahler scaling, Leonardo's rule, Tokunaga matrix.
- [06_fractal_dimension.ipynb](notebooks/06_fractal_dimension.ipynb) — reproduces SI Fig. S8 of Eloy et al. 2017 from the species-0 S3 champion. Five panels: 3D tree at 400 yr, log-log Horton ratios `R_n / R_l / R_d / R_a` + fractal dim `D`, time evolution from 25 yr to 500 yr, branch tapering scatter, area conservation. Added 2026-05-25 alongside two new public helpers in [src/mechatree/stats.py](src/mechatree/stats.py): `horton_summary(tree)` (per-Horton-stream means — needed because every MechaTree segment is a unit twig, so `strahler_summary.mean_length` is constant and only the chain length varies) and `horton_ratios(summary)` (log-linear fit returning the four ratios + `D = log R_n / log R_l`). Found and worked around a C++ staleness bug along the way: `Tree::setHorton` only calls `setStrahler` when `Strahler_distribution` is empty, so once Strahler is computed it never refreshes; `horton_summary` always forces `set_strahler()` first. Regression test in [tests/test_stats.py](tests/test_stats.py).

[docs/userguide.rst](docs/userguide.rst) gained a `.. _canopy-aware-wind:` section pointing at Step 17's DendroFlow bridge and corrected the `wind_fn` type signature (2-arg or 3-arg). A dedicated canopy-aware-wind notebook is a sensible follow-up but wasn't required for Step 13's planned scope.

### Step 15 notes (closed 2026-05-25)

Closed-form genome expressions via SymPy. Two paths landed together:

- **C++ side**: new vtable subclasses `CallbackSafety` / `CallbackAllocation` in [src/mechatree/_core/genome.h](src/mechatree/_core/genome.h). Each holds a function pointer + opaque `void* user_data`. The Cython shims (`_safety_callback` / `_allocation_callback`) cast `user_data` back to the Python callable and invoke it `with gil`. Errors in the Python callback are swallowed (zero-fallback + stderr log) because the C++ growth loop has no exception channel.
- **Python side**: [src/mechatree/sympy_genome.py](src/mechatree/sympy_genome.py) exposes `sympy_safety(expr)` and `sympy_allocation(p_seeds, p_leaves, phototropism)` factories. They `sympify` the expression, validate the free-symbol set against the allow-list (`nb_leaves`, `max_stress` for safety; `nb_leaves`, `vol_relative` for allocation), and `lambdify(modules="numpy")` into a Python callable that the C++ side calls per branch.
- **YAML**: `GenomeConfig.safety` / `p_seeds` / `p_leaves` / `phototropism` are now `float | str`. `models_from_config` dispatches: any field is a string → SymPy path; otherwise unchanged.
- **Optional extra**: `mechatree[sympy]` brings in `sympy>=1.12`. Without it, the SymPy path raises `ImportError` with the install hint.
- Tests at [tests/test_sympy_genome.py](tests/test_sympy_genome.py) (23 cases) — exercises the C++ vtable shim independently of SymPy, then the SymPy parsing/validation, end-to-end through `grow_tree`, and the YAML → `Config` → `models_from_config` round-trip. Recipe: [examples/sympy_genome.py](examples/sympy_genome.py) + [examples/sympy_genome.yaml](examples/sympy_genome.yaml).

### Step 17 notes (closed 2026-05-25)

[src/mechatree/wind/dendroflow.py](src/mechatree/wind/dendroflow.py) wraps DendroFlow's `BulkThinningBranchWindModel` (M4 of DendroFlow) as a streamwise canopy-mean `WindFn`. Wiring:
- New optional extra `mechatree[dendroflow]` ([pyproject.toml](pyproject.toml)) — DendroFlow isn't on PyPI yet, so the README hint is `uv pip install -e ../DendroFlow`.
- `WindFn` widened to `Callable[..., (x,y,z)]` and arity-detected at the call site, mirroring the `on_step` pattern. 2-arg callables still work.
- YAML `wind:` block with `model: dendroflow` (or `default`) selects the bridge automatically in `grow_tree` and `Forest` ([src/mechatree/config.py](src/mechatree/config.py)).
- Forest pools every tree's branches into a **single** `Cylinders` per generation; the resulting `canopy_mean` becomes the shared wind vector applied to all trees in pruning.
- Tests at [tests/test_wind_dendroflow.py](tests/test_wind_dendroflow.py); example recipe at [examples/dendroflow_wind.py](examples/dendroflow_wind.py) + [examples/dendroflow_wind.yaml](examples/dendroflow_wind.yaml).

### Notes on pending steps

- **Step 7 / Step 18** — release engineering (TestPyPI dry-run + real PyPI). No design work outstanding.
- **Step 19 — leaf transparency** — small, no design work outstanding; spec is the `tau ** shadow` formula.
- **Step 20 — Mesh3d cylindrical renderer** — small, no design work outstanding; `style="cylinders"` toggle on `plot_tree_3d`.
- **Step 21 — Evolution port** — non-trivial; settle the open design questions in the row before coding (package boundary, parallelism strategy, tournament shape).

## Out of scope (do not re-litigate unless explicitly asked)

- *(empty for now — the previous "Evolution is external" entry was reinstated as roadmap Step 21 on 2026-05-25 at the user's request.)*

## How to update this file

**Update CLAUDE.md after substantive changes.** When a roadmap step lands:

1. Flip its row in the table from `⬜ pending` / `🚧 wip` to `✅ done`.
2. Add a `### Step N notes (closed YYYY-MM-DD)` subsection below the table — one short paragraph + bullets, in the shape used for closed steps above.
3. Add a corresponding entry to [CHANGELOG.md](CHANGELOG.md) under `## [Unreleased]`, using Keep-a-Changelog `### Added` / `### Changed` / `### Fixed` headings.

For new work introduced mid-stream, add a fresh table row before doing the implementation so the roadmap stays the source of truth. Trivial changes (typos, single-line fixes, formatting-only) don't need an update; anything that adds/changes a module, public function, YAML key, or test file does.
