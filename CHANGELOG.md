# Changelog

All notable changes to this project will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added — Step 15: SymPy-callable genome

- [`mechatree._core.PyCallbackSafety`](src/mechatree/_core/_core.pyx) / `PyCallbackAllocation` — new Cython wrappers around C++ `CallbackSafety` / `CallbackAllocation` vtable subclasses. Each holds a Python callable; the Cython shims (`_safety_callback` / `_allocation_callback`) marshal `(nb_leaves, max_stress)` / `(nb_leaves, vol_relative)` back into Python under `with gil` so the C++ growth loop can call arbitrary Python decision functions. Errors raised by the callback fall back to `0` and log to stderr, since the C++ growth loop has no exception channel.
- [`mechatree.genome.CallbackSafety`](src/mechatree/genome.py) / `CallbackAllocation` — public aliases of the above; sit alongside `Constant*` / `Neural*`.
- [`mechatree.sympy_genome`](src/mechatree/sympy_genome.py) — `sympy_safety(expr)` and `sympy_allocation(p_seeds, p_leaves, phototropism)` factories. Expressions are `sympify`'d, validated against an allow-list of free symbols (`nb_leaves`, `max_stress` for safety; `nb_leaves`, `vol_relative` for allocation), then `lambdify(modules="numpy")`-compiled and wrapped in `CallbackSafety` / `CallbackAllocation`.
- [`mechatree.config.GenomeConfig`](src/mechatree/config.py) fields `safety` / `p_seeds` / `p_leaves` / `phototropism` accept `float | str`. When any is a string, [`models_from_config`](src/mechatree/genome.py) dispatches to the SymPy path; otherwise unchanged.
- [`examples/sympy_genome.py`](examples/sympy_genome.py) + [`examples/sympy_genome.yaml`](examples/sympy_genome.yaml) — runnable recipe.
- [`tests/test_sympy_genome.py`](tests/test_sympy_genome.py) — 23 tests (`pytest.importorskip("sympy")`-gated) covering the C++ vtable shim, SymPy parsing/validation, end-to-end `grow_tree`, and YAML wiring.
- `pyproject.toml`: new optional extra `sympy = ["sympy>=1.12"]`. Without it, `mechatree.sympy_genome` raises `ImportError` with the install hint.

### Changed — Step 15

- `GenomeConfig.__post_init__` skips the numeric-range checks (positivity, sum-≤-1, phototropism-in-[0,1]) for string fields — those expressions are validated by SymPy at model-build time instead.
- [`docs/userguide.rst`](docs/userguide.rst) gained a SymPy-genome subsection under *Customizing the simulation* and a new entry in *Examples*.

### Added — Step 13: Forester-facing notebooks

- Five plotly-inline Jupyter tutorials under [`notebooks/`](notebooks/), each a literate-programming companion to an `examples/` script. End-to-end execution verified via `jupyter nbconvert --execute`; checked-in `.ipynb` files are output-stripped per the [`notebooks/README.md`](notebooks/README.md) convention.
  - [`01_grow_one_tree.ipynb`](notebooks/01_grow_one_tree.ipynb) — single-tree pipeline + Strahler-coloured 3D + leaf overlay.
  - [`02_forest_under_wind.ipynb`](notebooks/02_forest_under_wind.ipynb) — population/biomass dynamics, top-down stand, Yoda −3/2 self-thinning.
  - [`03_neural_genome.ipynb`](notebooks/03_neural_genome.ipynb) — S3 champions vs constant genome, side-by-side.
  - [`04_custom_growth_law.ipynb`](notebooks/04_custom_growth_law.ipynb) — plug-in `wind_fn` / `Sun` / `ConstantSafety` / `ConstantAllocation`.
  - [`05_strahler_diagnostics.ipynb`](notebooks/05_strahler_diagnostics.ipynb) — Horton–Strahler scaling + Leonardo's rule + Tokunaga matrix.

### Changed — Step 13

- [`docs/userguide.rst`](docs/userguide.rst): corrected the `wind_fn` type signature in *Customizing the simulation* to cover both the 2-arg and the 3-arg shapes introduced in Step 17. Added a new *Canopy-aware wind via DendroFlow* section with a YAML + Python recipe. Listed `examples/dendroflow_wind.py` / `.yaml` under *Examples*.

### Added — Step 17: DendroFlow wind bridge (M6 of DendroFlow)

- [`mechatree.wind.dendroflow`](src/mechatree/wind/dendroflow.py) — `BranchWindBridge` wraps DendroFlow's `BulkThinningBranchWindModel` as a 3-arg `WindFn(generation, rng, context)`. `pytree_to_cylinders` / `forest_to_cylinders` snapshot the current geometry into a DendroFlow `Cylinders`; the per-generation `canopy_mean` becomes the streamwise wind vector applied during pruning. `BranchWindBridge.last_result` exposes the underlying `BranchWindResult` for diagnostics.
- [`mechatree.config.WindConfig`](src/mechatree/config.py) — new `wind:` block in YAML. `model: dendroflow` selects the bridge in `grow_tree` and `Forest` automatically; `U_infty`, `z_centers`, `H`, `C_D`, `z_representative` configure DendroFlow's `BulkThinningBranchWindModel`. Validation: `z_centers` must be monotone and cover `z = 0` so the trunk base isn't dropped.
- [`examples/dendroflow_wind.py`](examples/dendroflow_wind.py) + [`examples/dendroflow_wind.yaml`](examples/dendroflow_wind.yaml) — minimal recipe.
- [`tests/test_wind_dendroflow.py`](tests/test_wind_dendroflow.py) — DendroFlow-gated test module (`pytest.importorskip("dendroflow")`) covering the cylinders adapter, the bridge contract, YAML wiring, forest pooling, arity backwards-compat, and `WindConfig` validation.
- `pyproject.toml`: new optional extra `dendroflow = ["dendroflow>=0.1"]`. DendroFlow isn't on PyPI yet — for local dev use `uv pip install -e ../DendroFlow` alongside.

### Changed — Step 17

- `mechatree.simulate.WindFn` type widened from `Callable[[int, np.random.Generator], (x,y,z)]` to `Callable[..., (x,y,z)]`. Two arities are now accepted (detected at call time via the same `_callback_arity` helper used for `on_step`):
  - `wind_fn(generation, rng)` — classic shape, still the default.
  - `wind_fn(generation, rng, context)` — receives the live `PyTree` (in `grow_tree`) or the `Forest` (in `Forest.step`).
- `grow_tree` and `Forest.step` now pass the third `context` argument when the wind callable accepts it. Existing 2-arg callables work unchanged.
- `Forest` lost its hard import of `default_wind_fn`; both `grow_tree` and `Forest.__post_init__` route through `mechatree.simulate._resolve_wind_fn(config)` to pick between the Fortran-faithful default and the DendroFlow bridge based on the YAML `wind:` block.

### Documentation — Step 17

- [`CLAUDE.md`](CLAUDE.md) restructured to a DendroFlow-style milestone table at the top of the roadmap section, with per-step `### Step N notes (closed YYYY-MM-DD)` subsections below. Added Step 17 row + a new Step 7 (TestPyPI shake-out) and Step 18 (stable PyPI 0.1.0) split. Added a "How to update this file" section codifying the table-flip / CHANGELOG cadence.

### Added — Step 1 (foundation, retained for first-release notes)

- Modern repo scaffold: `pyproject.toml` (setuptools + Cython build backend), `src/mechatree/` package skeleton, `tests/` with an import smoke test, ruff + pre-commit + pytest configuration, GitHub Actions CI, `CLAUDE.md` orientation guide.
- `legacy_fortran/` — Fortran90 sources, `Makefile`, and `.ini` parameter files from Eloy et al., *Nat Commun* 8:1014 (2017), kept as a read-only reference.
- `archive/` — preserved the 2017 Python + Cython + C++ port by **Diego Bengochea Paz** (ORCID [0000-0002-0835-3981](https://orcid.org/0000-0002-0835-3981)) — both `PyTreeLib/` and the renamed `CommentedLibrary_and_Doc/` — for provenance.
- `benchmarks/` — standalone micro-benchmarks (creation, removal, topology walk, RSS) with a `baseline.md` capturing before/after numbers.
- New test coverage for `_core`: `tests/test_pytree_remove.py`, `tests/test_pytree_classification.py`, `tests/test_pytree_topology.py`, and a `tests/test_pytree_memory.py` regression for the C++ leak.
- `slow` pytest marker registered in `pyproject.toml` for large-tree tests.

### Fixed
- `Tree` and `Branch` had no destructors; `PyTree` had no `__dealloc__`. Every PyTree leaked its entire branch graph on garbage collection. Added `~Tree()` (deletes every branch), an explicit defaulted `~Branch()`, and `PyTree.__dealloc__` so RSS now plateaus instead of growing linearly with the number of trees created.
- `Tree::removeBranch` erased pointers from `tree_branches` without `delete`-ing them. Now deletes the entire subtree first, then erases.
- `Branch::removeDescendants` had a second loop that ran after the first cleared `children`, so grandchild parent-links were never reset. Removed the method entirely; `Tree::removeBranch` now does the teardown directly.
- `setStrahler()` and `setHorton()` accumulated into their distribution maps across repeated calls. Both now clear the map at the start so calls are idempotent.
- `Tree::getBranch` allocated and returned a fresh `new Branch()` on out-of-range — both a leak and a silent corruption. Now throws `std::out_of_range`.
- `Branch` constructor left `strahler` and `horton` uninitialized (UB on read before `setStrahler`/`setHorton`). Default-initialized to 0 in the header.

### Changed
- **Breaking**: out-of-range indices and missing properties now raise Python exceptions (surfaced from C++ `std::out_of_range` / `std::invalid_argument` via Cython `except +`) instead of printing to stderr and returning sentinel values (`-1`, `-2`, `0`, default Branch). Callers that relied on those sentinels need to switch to try/except.
- **Breaking**: `PyTree.get_brothers_index` and `PyTree.get_children_index` now always return a `list` (empty if no brothers / no children). Previously `get_brothers_index` returned a bare `int` when there was exactly one brother and `0` when none; `get_children_index` returned `0` when none.
- Modernized the C++/Cython core to C++17 / Cython 3 idioms: dropped `using namespace std;` from headers, qualified everything with `std::`, marked read-only methods `const`, replaced hand-rolled iterator loops with range-for + `auto`, switched container getters (`getChildren`, `getStrahlerDistribution`, `getHortonDistribution`) to return `const&`, made `Tree` non-copyable to prevent accidental double-frees of the owned `Branch*`s. No behavior change beyond the breaking notes above.

### Removed (from the repo tree, not from disk)
- The 15 GB `FORESTArticlePNAS/` folder (simulation outputs, zips, MATLAB scripts, compiled 2014 binaries) was moved out to `/Users/Ch/Documents/Python/Eloy2017_NatComm_archive/` — local-only, not committed.
