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
| 19 | Leaf transparency in light model    | ✅ done    | Global scalar `light.leaf_transparency` (`tau`) wired into [`mechatree.light.interception`](src/mechatree/light/interception.py): light at depth `i` in a shadow cell is `tau^i`. Default `tau = 0.5` matches Eloy et al. (2017); `tau = 0` recovers the Fortran binary topmost-wins regime; `tau = 1` makes leaves fully transparent. New field on [`LightConfig`](src/mechatree/config.py), surfaced via YAML `light:`. |
| 20 | Mesh3d cylindrical tree renderer    | ⬜ pending | Today's [`plot_tree_3d`](src/mechatree/plotting/_mechanics.py) uses pixel-width plotly `Scatter3d` lines bucketed by Strahler order, so individual branch diameters are flattened and the canopy looks pencil-thin compared to the Blender renders in the paper. Add a `style="cylinders"` path that extrudes each branch as a tapered `Mesh3d` cylinder (parent diameter → child diameter), with the existing `style="lines"` kept as the fast default. ~60 LOC + a small cylinder-primitive helper. |
| 21 | Evolution port (`EvoluAlgo` / `mod_evolu`) | ⬜ pending | Port the Fortran evolutionary algorithm so the Nat Comms tournament simulations are reproducible end-to-end in Python: per-generation fitness evaluation via `grow_tree` on a population of `NeuralSafety` / `NeuralAllocation` champions, tournament selection, neural-net weight mutation + crossover, and a CSV/JSON dump matching `S3.dat` so [`scripts/strategies_single_tree.py`](scripts/strategies_single_tree.py) can keep loading champions the same way. Open design questions to settle before coding: separate package vs in-tree `mechatree.evolution` module; whether to parallelise fitness eval via `multiprocessing` or keep it sequential; whether to mirror the Fortran tournament structure verbatim or restructure around a `numpy`-friendly `Population` dataclass. Reverses the "Evolution is external" call in *Design principles for the simulator port*. |
| 22 | Unified figure style                | ✅ done    | [`mechatree.plotting.figstyle`](src/mechatree/plotting/figstyle.py) registers a MATLAB-look plotly template (white bg, 4-sided box, inside ticks, Helvetica 11 pt, palette in `COLORS`, sized canvases in `SIZES`) under `pio.templates["mechatree"]`. All seven plotting helpers + the 6 forester notebooks + [examples/grow_a_forest.py](examples/grow_a_forest.py) route through it; `figstyle.figure` / `figstyle.subplots` / `figstyle.figure_3d` / `figstyle.save` mirror the SoftMobility API. Ships four Strahler colormaps (`jet`, `cool`, `parula`, `rainbow`) selectable via `set_strahler_cmap`; notebook 06 renders a 4-panel benchmark so the eventual default can be picked visually. CI lint also fixed (gated `[tool.uv.sources].dendroflow` on `extra = "dendroflow"`). |
| 23 | Flat top-level API (`import mechatree as mt`) | ⬜ pending | Curate a flat public surface so notebooks/users can do `import mechatree as mt` and then call `mt.load_config(...)`, `mt.load_champion(...)`, `mt.grow_tree(...)`, `mt.Forest(...)`, `mt.plot_tree_3d(...)`, `mt.figstyle.apply()`, etc., instead of remembering the subpackage path (`from mechatree.simulate import grow_tree`, …). User-facing motivation: "I can choose which functions are really useful" — only the genuinely public API is re-exported at the top level; private helpers (`_callback_arity`, `_resolve_wind_fn`, `_decode_angles`, `_core` types, etc.) stay buried. Concretely: edit [`src/mechatree/__init__.py`](src/mechatree/__init__.py) to re-export `load_config`, `Config`, `TreeConfig`, `LightConfig`, `ForestConfig`, `GenomeConfig`, `grow_tree`, `Forest`, `ForestStats`, `Sun`, `Leaves`, `load_champion`, `load_all_champions`, `champion_angles`, `models_from_config`, `ConstantSafety`, `ConstantAllocation`, `NeuralSafety`, `NeuralAllocation`, `CallbackSafety`, `CallbackAllocation`, `plot_tree_3d`, `plot_forest_topdown`, `plot_2d`, `plot_3d`, `plot_self_thinning`, `plot_allocation`, `plot_strahler_diagnostics`, `figstyle`, `horton_summary`, `strahler_summary`, `horton_ratios`, `distance_to_leaves`, `mean_distance_to_leaves`, `leonardo_ratios`, `tokunaga_matrix`. Sweep notebooks 01–06 + examples to use `mt.*` form; keep deep imports working for backward compat. Open question: include or omit the `dendroflow`-extra surface (e.g. `mt.BranchWindBridge`)? Likely include behind an `if find_spec(...) is not None` guard so a bare install doesn't error on import. |

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

**Notebook 06 `<l>` follow-up (2026-05-25):** the original SI Fig. S8(b) `<l>` trace was using `HortonSummary.mean_length`, which is the per-Horton-stream chain length (sum of unit segments per stream), not the per-branch recursive distance-to-leaves the paper actually plots. Added two helpers in [src/mechatree/stats.py](src/mechatree/stats.py): `distance_to_leaves(tree)` — length-aware per-branch recursion: terminals → `length / 2`, internals → `length + nb_leaves`-weighted mean of children's distance. Mirrors Fortran `b%distance_leaves` in `legacy_fortran/mod_tree.f90:1174-1203`, with the Fortran's `+1.0` increment replaced by the branch's actual length so the metric stays correct when post-pruning chain merger ([tree.cpp:502](src/mechatree/_core/tree.cpp#L502)) grows branches past unit length. `mean_distance_to_leaves(tree)` aggregates per Horton rank. `horton_ratios` gained a `mean_length_override` kwarg and a `max_rank` cap; notebook 06 threads `mean_distance_to_leaves` through the former and sets the latter to 7 (matches SI Fig. S12, where the top ranks sit on too few branches to be stable). `HortonSummary.mean_length` keeps its original meaning (per-stream chain length) — useful for histograms of stream geometry, not for the paper's `R_l`.

**`volume_ratio_leaf` default 8.0 → 4.0 (2026-05-25):** the original MechaTree default `volume_ratio_leaf = 8.0` mirrored `legacy_fortran/Forest.ini` (an older PNAS-submission tarball). The Nat Commun revision used `VolumeRatioLeaf = 4.0d0` per `~/Documents/Arbres/FORTRAN/ART_Revision2/Forest_reference.ini` and `…/ART_Revision2b/Forest_reference.ini` (the figure caption in SI Fig. S12 confirms `V_prod = 4 V_0 l`). At 200 yr with the species-0 S3 champion the over-production halves: 77k → 8k branches, matching the paper's ~10⁴. [config.py:36](src/mechatree/config.py#L36) + [examples/forest.yaml](examples/forest.yaml) + [tests/test_config.py](tests/test_config.py) updated; `legacy_fortran/Forest.ini` left untouched per the provenance rule.

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

### Step 19 notes (closed 2026-05-25)

Per-cell light interception is now `light = tau ** depth`, where `depth` is each leaf's 0-indexed rank in Z'-descending order within its shadow cell (topmost leaf has depth 0).
- **Algorithm**: [src/mechatree/light/interception.py](src/mechatree/light/interception.py) computes per-leaf depth in one vectorised pass — `np.unique(return_index=True)` on the lexsort-sorted cell keys gives each group's start, `np.repeat` broadcasts that back, and `depth = arange(n) - group_start_sorted`. The final write is `leaves.light_per_direction[:, k] = leaf_transparency ** depth`. NumPy follows IEEE `0**0 = 1, 0**k = 0` for `k >= 1`, so `tau = 0` recovers the previous binary topmost-wins behaviour exactly with no special-case branch.
- **Default `tau = 0.5`** (Eloy et al., Nat Commun 2017). This is a behavioural change vs. the prior binary default: understorey leaves now contribute photosynthate they didn't before, so quantitative tree size / leaf count under the default config shifts. Structural tests (reproducibility under seed, callback firing, etc.) are unaffected and still pass.
- **YAML / API**: [`LightConfig.leaf_transparency`](src/mechatree/config.py) is the canonical knob, surfaced via the YAML `light:` block. Validated to `[0, 1]` in `LightConfig.__post_init__` and again in `intercept()`. Threaded through [`mechatree.simulate.grow_tree`](src/mechatree/simulate.py) and [`mechatree.forest.Forest.step`](src/mechatree/forest.py); [`mechatree.plotting`](src/mechatree/plotting/_mechanics.py) uses the function default.
- **Tests**: [tests/test_light.py](tests/test_light.py) adds `test_intercept_default_tau_attenuates_shadowed_leaf` (default 0.5 ⇒ shadowed leaf sees 0.5), `test_intercept_tau_zero_recovers_binary`, `test_intercept_tau_one_makes_leaves_transparent`, `test_intercept_tau_half_geometric_attenuation` (4-leaf stack ⇒ `[1, 1/2, 1/4, 1/8]`), and validation. [tests/test_config.py](tests/test_config.py) locks in `LightConfig().leaf_transparency == 0.5` and the `[0, 1]` validator.

### Step 22 notes (closed 2026-05-25)

Unified figure styling, MATLAB / Nat Commun 2017 look, plotly-native.

- **Module**: [src/mechatree/plotting/figstyle.py](src/mechatree/plotting/figstyle.py) registers a `go.layout.Template` under `pio.templates["mechatree"]` and exposes `apply()` / `figure(size, aspect)` / `subplots(...)` / `figure_3d(...)` / `save(fig, name)` plus the `COLORS`, `SIZES`, `FONT` dicts. API names match [SoftMobility's `figstyle.py`](../SoftMobility/softmobility/classes/figstyle.py); template shape mirrors [SoftMobility's `figstyle_plotly_legacy.py`](../SoftMobility/softmobility/drafts/figstyle_plotly_legacy.py) but with `ticks="inside"` (true MATLAB default) instead of `"outside"`. The 4-sided `mirror=True` frame, white background, Helvetica 11 pt, and `colorway` defaults are baked into the template, so any call to `figstyle.figure(...)` returns a styled canvas without the caller having to invoke `apply()`.
- **Strahler colormaps**: four 10-stop tables sampled from MATLAB's `colormap(name)` — `jet` (default), `cool` (literal match for the `colormap(cool)` line at [`../Eloy2017_NatComm_archive/plot_stat_single_tree.m:37`](../Eloy2017_NatComm_archive/plot_stat_single_tree.m)), `parula` (post-R2014b default), and `rainbow` (legacy MechaTree palette from `_palette.py`). Switch via `set_strahler_cmap("...")`. The benchmark cell in [notebooks/06_fractal_dimension.ipynb](notebooks/06_fractal_dimension.ipynb) renders the species-0 champion under all four side-by-side so the visual winner can be picked before flipping `DEFAULT_STRAHLER_CMAP`.
- **Wire-up**: every plotly call in the library now routes through `figstyle.figure*` ([_2d.py](src/mechatree/plotting/_2d.py), [_3d.py](src/mechatree/plotting/_3d.py), [_mechanics.py](src/mechatree/plotting/_mechanics.py), [_stats.py](src/mechatree/plotting/_stats.py)). [examples/grow_a_forest.py](examples/grow_a_forest.py)'s inline twin-axis figure was migrated as well. All 6 notebooks call `figstyle.apply()` in their first code cell; notebooks 02 / 03 / 04 / 06 also swap their direct `make_subplots` for `figstyle.subplots`. Hardcoded `"forestgreen"` / `"saddlebrown"` / `"red"` / `"blue"` / `"black"` / `"magenta"` / `"cyan"` literals from `_stats.py` and `grow_a_forest.py` were replaced with `figstyle.COLORS[...]` entries.
- **Comparison script**: [examples/figstyle_compare.py](examples/figstyle_compare.py) — three side-by-side benchmarks (Strahler palettes ×4, font families ×4, ticks-inside vs ticks-outside × 4-sided frame vs no-top-right). `uv run python examples/figstyle_compare.py` opens three browser tabs; pick winners then edit `DEFAULT_STRAHLER_CMAP` / `FONT` / `_axis_style()` in `figstyle.py`.
- **Tests**: [tests/test_figstyle.py](tests/test_figstyle.py) — 7 cases covering template registration, sized canvas dimensions, axis attributes (`mirror=True`, `ticks="inside"`, `showgrid=False`), `COLORS` hex validation, and the Strahler cmap switch.
- **CI fix**: [pyproject.toml](pyproject.toml)'s `[tool.uv.sources].dendroflow` block was unconditional, so `uv pip install -e ".[dev]"` on the GH Actions runner tried (and failed) to resolve `../DendroFlow`. Gated on `extra = "dendroflow"` so the sibling-checkout source is only consulted when the dendroflow extra is requested.

### Notes on pending steps

- **Step 7 / Step 18** — release engineering (TestPyPI dry-run + real PyPI). No design work outstanding.
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
