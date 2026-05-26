# Changelog

All notable changes to this project will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed — Legacy folder consolidation (2026-05-26)

- Consolidated legacy material under a single top-level `legacy/` directory.
  - `legacy_fortran/` → `legacy/fortran/` (history preserved via `git mv`).
  - `archive/CommentedLibrary_and_Doc/` → `legacy/pytree/` (history preserved).
  - New `legacy/matlab/` holds the MATLAB analysis scripts from the paper workflow (previously local-only at `~/Documents/Python/Eloy2017_NatComm_archive/*.m`).
  - New `legacy/pdf/` holds the paper (`s41467-017-00995-6.pdf`) and the Supplementary Information (`SI11.pdf`).
  - The near-duplicate `archive/PyTreeLib/` was dropped; its three example scripts already have modern equivalents under `examples/` (`random_growth.py`, `sap_transport.py`, `self_avoiding.py`).
- New consolidated [`legacy/README.md`](legacy/README.md) replaces the two old per-folder READMEs and maps each MATLAB script to its modern Python port.
- Path references swept across source comments, docstrings, examples, benchmarks, [README.md](README.md), [CONTRIBUTING.md](CONTRIBUTING.md), [docs/about.rst](docs/about.rst), [docs/design.rst](docs/design.rst), [docs/userguide.rst](docs/userguide.rst), [benchmarks/baseline.md](benchmarks/baseline.md), and [CLAUDE.md](CLAUDE.md).
- Config: [pyproject.toml](pyproject.toml) `[tool.ruff].extend-exclude` and [.pre-commit-config.yaml](.pre-commit-config.yaml) `exclude:` both retargeted to `^legacy/`.
- No behavioural change in the modern package — comment/docstring-level rewrites only.

### Added — Step 25: Tunable wind statistics + native bulk-thinning + storm replay + explainer notebook

- [`mechatree.wind.distributions`](src/mechatree/wind/distributions.py) — SymPy-driven inverse-CDF sampling. `Distribution(cdf_expr, var_name, support).sampler()` compiles to a `(rng, n) -> ndarray` callable. Symbolic inversion via `sympy.solve` validated by round-tripping through the original CDF on a 25-point probe (rejects branches that pass an endpoint check but pick the wrong complex root). Falls back to a 1024-point `np.interp` lookup when symbolic inversion fails. Built-in pure-NumPy fast paths skip SymPy entirely: `default_amplitude_sampler` reproduces the Fortran `0.835 - log(U)/6`, `uniform_angle_sampler` is uniform on `[0, 2π)`, `legacy_angle_sampler` returns the original 4 fixed sensing angles.
- [`mechatree.wind.bulk_thinning`](src/mechatree/wind/bulk_thinning.py) — self-contained native port of DendroFlow's `BulkThinningBranchWindModel`. Operates directly on `(start, axis, D, L)` arrays from [`PyTree.get_branch_data_batch`](src/mechatree/_core/_core.pyx) — no DataFrame round-trip, no DendroFlow dependency for the common case. The actuator-disk per-z-layer thinning formula is line-for-line equivalent to [DendroFlow's branch.py](../DendroFlow/src/dendroflow/wind/branch.py). New types: `BulkThinningParams` (default `uniform()` is `U = 1` over `z ∈ [0, 50]`, `H = 0.5` — no boundary layer, no opt-in), `BulkThinningResult` (per-branch `U_branch`, per-layer `U_canopy`, plus `canopy_mean` rotated into the world frame), `BulkThinningWindBridge` (stateful `WindFn` accepting optional angle + amplitude samplers).
- **Rotation to face the storm**: the bridge samples a per-iteration storm direction `θ` from the configured `angle_sampler`; the canopy is conceptually rotated so wind = `+x`, which collapses to `cos I = axis · wind_direction` so no rotated-coordinate allocations are needed. Mirrors the light-interception pattern. The legacy [DendroFlow bridge](src/mechatree/wind/dendroflow.py) still solves in world `+x` only (deferred follow-up); the native module is the recommended path when a tunable storm direction is wanted.
- [`mechatree.wind.replay`](src/mechatree/wind/replay.py) — `run_storm_replay(forest, wind_fn, …)` replays the Step-24 fixed-point loop pass-by-pass, snapshotting per-iteration cuts, per-tree affected counts, and per-branch alive masks (matched back to the immutable pre-storm geometry by branch base position, which `Tree::removeBranches` preserves bit-identically for survivors). Pairs with [`mechatree.plotting.plot_storm_replay`](src/mechatree/plotting/_storm.py) — multi-panel 3D figure where surviving branches render in green and pruned branches in **black** on top of the pre-storm canopy.
- New `WindConfig` fields: `amplitude_cdf: str | None`, `angle_cdf: str | None`, `angle_samples: int = 4`. Leaving them unset preserves today's behaviour byte-identically. Validation rejects empty strings.
- New `wind.model = "native"` option in `_resolve_wind_fn` ([src/mechatree/simulate.py](src/mechatree/simulate.py)). `make_default_wind_fn(amplitude_sampler, angle_sampler)` builds a tunable storm-wind closure that drops back to the Fortran defaults on either axis when not configured.
- [`notebooks/07_wind_models.ipynb`](notebooks/07_wind_models.ipynb) — explainer notebook. Five sections: mechanics intro + sensing sweep, storm distributions (amplitude + angle), native bulk-thinning with rotation demo, Step-24-in-pictures (grow 10 trees with zero wind for 80 generations, hit with a single big storm, replay iter-by-iter with pruned branches in black), and a prevailing-wind lab section.
- Top-level API: `mt.Distribution`, `mt.default_amplitude_sampler`, `mt.uniform_angle_sampler`, `mt.make_default_wind_fn`, `mt.BulkThinningParams`, `mt.BulkThinningResult`, `mt.BulkThinningWindBridge`, `mt.make_bulk_thinning_wind_fn`, `mt.StormPreSnapshot`, `mt.StormSnapshot`, `mt.run_storm_replay`, `mt.plot_storm_replay`, `mt.default_wind_fn`. Flat surface now 67 names.
- 27 new tests across [tests/test_wind_distributions.py](tests/test_wind_distributions.py) (10), [tests/test_wind_bulk_thinning.py](tests/test_wind_bulk_thinning.py) (11), and [tests/test_wind_replay.py](tests/test_wind_replay.py) (6). 368/368 tests green; ruff clean.

### Changed — Step 25

- `_resolve_wind_fn` ([src/mechatree/simulate.py](src/mechatree/simulate.py)) gained the `model: native` branch and now builds amplitude/angle samplers from the new YAML CDF fields before dispatching to whichever wind model is selected. The Fortran `default_wind_fn` is kept unchanged; `make_default_wind_fn(...)` is the new factory for a tunable storm-wind closure.
- `WindConfig.model` choices widened to `default | native | dendroflow` (was `default | dendroflow`). Validation message updated.

### Deferred to Step 25b / DendroFlow

- Extend the C++ `calculate_stresses` to accept a runtime angle list so sensing also honours `angle_cdf` / `angle_samples` (the YAML knob is in place but currently ignored sensing-side; sensing still uses the legacy 4 hardcoded angles).
- Rotation in the legacy DendroFlow bridge (the native module is the recommended path when storm direction matters).
- Per-branch wind colour overlay in `plot_tree_3d` from `BulkThinningResult.U_branch`.
- DendroFlow cell-size sensitivity study + single-cylinder k-ε flow visualisation (lives in the DendroFlow repo).

### Added — Step 24: Coupled wind ↔ pruning fixed-point loop

- [`mechatree.config.WindConfig`](src/mechatree/config.py) gained two knobs: `max_pruning_iterations: int = 8` caps the wind ↔ prune inner loop; `wind_convergence_eps_rel: float = 0.01` controls the relative-tolerance early-exit on the canopy-mean wind delta. Validation in `__post_init__` (`max_pruning_iterations >= 1`, `wind_convergence_eps_rel >= 0`).
- [`mechatree.simulate._prune_to_fixed_point`](src/mechatree/simulate.py) — helper that wraps `wind_fn → prune` in a fixed-point iteration when the wind callable is 3-arg (canopy-aware). Exits on the first of `n_cut == 0` (true fixed point), `|wind_new − wind_old| / |wind_old| < ε_rel` (cheap early-exit), or the cap. For 2-arg default wind it short-circuits to a single pass — byte-identical to before. Wired into `grow_tree`; new `TreeStats.n_wind_iterations` field reports the per-step iteration count.
- [`mechatree.forest.Forest.step`](src/mechatree/forest.py) gained the forest-wide variant of the same loop (pool all trees → wind → prune every tree → check the sum). New `ForestStats.n_wind_iterations` field. Cap-hit emits one `RuntimeWarning` per `Forest.run` call (gated by `_cap_hit_warned`, reset at the top of `run`). Per-tree per-iter cost stays bounded thanks to the Step 21b-style batched accessor below.
- [`mechatree._core.PyTree.get_branch_data_batch`](src/mechatree/_core/_core.pyx) — new Cython accessor that returns `(start, axis, diameter, length)` arrays for every branch in one call. Used by [`pytree_to_cylinders`](src/mechatree/wind/dendroflow.py) and [`forest_to_cylinders`](src/mechatree/wind/dendroflow.py) instead of the per-branch Python `[get_location, get_unit_t, get_diameter, get_length]` round-trip. **5–7× faster** at 12 k–129 k branches; also speeds up today's single-pass DendroFlow path by 24–49 % as a free win.
- [`benchmarks/bench_wind_fixed_point.py`](benchmarks/bench_wind_fixed_point.py) — Phase-0 sizing harness used to gate the design before implementation. Measures `forest_to_cylinders` + `BulkThinningBranchWindModel.compute` cost in isolation at 10 / 100 / 1000 (and optionally 20 k) trees, projects per-step overhead at 1/2/4/8 inner iterations, and probes the ε-tolerance triggering on a synthetic "sparse storm" (drop 0.5 % of branches → measure pooled-wind delta).
- [`tests/test_wind_fixed_point.py`](tests/test_wind_fixed_point.py) — 13 cases covering: default 2-arg wind stays single-pass; canopy-aware calm exits in one iter; storms iterate within the cap; cap=1 recovers single-pass; higher cap cuts at least as many branches as cap=1 under intensifying wind (the correct physical direction for bulk-thinning models — see note below); cap-hit warns exactly once per run; cap=1 never warns; ε early-exit stops iter 2 before its prune sweep when wind delta is small; ε=0 never iterates fewer than ε=0.01; direct unit test of `_prune_to_fixed_point`; soft wall-clock guard at 10-tree scale (`< 3×` cap=1 baseline).
- [`tests/test_pytree_api.py`](tests/test_pytree_api.py) — 2 new cases pin `get_branch_data_batch` against the per-branch path on a grown tree + a seed-only tree.
- [`tests/test_config.py`](tests/test_config.py) — 4 new cases pin the WindConfig defaults / validation / YAML round-trip for the Step-24 knobs.
- [`examples/dendroflow_wind.py`](examples/dendroflow_wind.py) + [`examples/dendroflow_wind.yaml`](examples/dendroflow_wind.yaml) — per-step printout now shows `iters` + `pruned`; YAML documents the two new knobs.
- [`docs/userguide.rst`](docs/userguide.rst) — new "Coupled wind ↔ pruning fixed-point loop (Step 24)" subsection under canopy-aware wind, including the directional note on bulk-thinning models.

### Changed — Step 24

- Default `forest_to_cylinders` / `pytree_to_cylinders` code path: replaced the per-branch Python loop with one batched Cython call per tree + a single NumPy `concatenate`. Cuts the pool cost from 8 ms → 1.2 ms at 10-tree / 12 k-branch scale, 71 ms → 13.5 ms at 1000-tree / 129 k-branch scale.
- Per-generation order in `grow_tree` and `Forest.step`: the `wind → prune` block is now an iterating helper rather than a single pair of calls. For 2-arg default wind the iteration runs once (byte-identical to before); for 3-arg canopy-aware winds it iterates to the cap, the ε threshold, or `n_cut == 0`.

### Note — Direction of the fix for bulk-thinning canopy models

For DendroFlow's `BulkThinningBranchWindModel` (the only canopy-aware wind currently shipped via the bridge) and any model where thinning the canopy *raises* the mean wind on the survivors, Step 24's fixed-point loop converges to a canopy that is **smaller** than the single-pass would produce, not larger. The single-pass under-prunes because it scores the surviving canopy against the weaker pre-pruning wind. The original roadmap description ("at storm tails it over-prunes") was directionally inverted for bulk-thinning models; it would still hold for a wind model where thinning *lowers* the mean wind on survivors (e.g. a shelter-style model). The fix is correct in either case: re-evaluating the wind after each cut sweep settles on the self-consistent canopy ↔ wind state.

### Added — Step 21b: Island-scale light pipeline (Phase A + B + C)

- [benchmarks/bench_forest_scale.py](benchmarks/bench_forest_scale.py) — per-phase wall-clock probe for `Forest.step` at any scale. Supports `--full` (paper config: R = 200 L, N_init = 20 000, 1000 gens), `--profile` (cProfile dump), and extrapolation to paper scale from a small-N probe.
- [src/mechatree/_core/_core.pyx](src/mechatree/_core/_core.pyx): new Cython accessors `PyTree.get_leaf_tips_batch()` (returns `(tips: (n, 3) float64, branch_indices: (n,) int32)` for every leaf in one call) and `PyTree.set_lights_batch(branch_indices, values)`. Replace the per-leaf Python loops in [`mechatree.light.extract_leaves`](src/mechatree/light/leaves.py) and [`mechatree.light.aggregate_onto_trees`](src/mechatree/light/interception.py). Measured **12.5×** faster `extract_leaves` at 1000 trees / **8.5×** at 20 k trees; **3×** faster `aggregate_onto_trees`.
- [src/mechatree/_core/light.h](src/mechatree/_core/light.h) + [light.cpp](src/mechatree/_core/light.cpp): new C++ kernel `light_intercept(leaf_locations, sun_elev, sun_azim, size_leaf, tau, light_per_direction)` ported from the Python `intercept`'s per-direction NumPy lexsort path. One scratch `std::vector<LeafEntry>` reused across all directions; cell key packed as `(x_cell << 32) | uint32(y_cell)` so no hash table is needed. Wired through Cython as `mechatree._core._core.light_intercept_kernel`. [`mechatree.light.intercept`](src/mechatree/light/interception.py) calls it; the public signature is unchanged. Measured `intercept` **2.2×** faster at island scale.
- [setup.py](setup.py) builds `light.cpp` alongside the existing core sources.
- 4 new tests in [tests/test_pytree_api.py](tests/test_pytree_api.py) pin the new accessors (equivalence with the per-leaf path, empty-tree handling, round-trip set/get, length-mismatch raises). The existing 26 cases in `test_light.py` continue to pass against the C++ kernel — including `tau = 0` (binary topmost-wins), `tau = 1` (transparent), within-cell ties, single leaf, empty leaves.
- **Phase C — parallel-over-directions in the C++ kernel** (`light.cpp`): wrapped the outer per-direction loop in `#pragma omp parallel` with a thread-local `std::vector<LeafEntry>` scratch. Each direction writes a disjoint column of the output, so no critical sections or atomics are needed. Also replaced the per-leaf `std::pow(tau, depth)` with a running multiply since depth advances by exactly 1 along the sorted walk — eliminates ~14M `pow` calls per generation at island scale.
- [setup.py](setup.py) gained an OpenMP probe that picks per-platform flags (`/openmp:llvm` on MSVC, brew-prefix libomp on macOS, `-fopenmp` on Linux probed via a test-compile) and falls back silently to serial when OpenMP isn't available. Set `MECHATREE_NO_OPENMP=1` to force serial for reproducibility debugging.
- Measured thread scaling on Apple M-series, 8 cores, n = 200k leaves × 32 dirs: 1 thread = 445 ms/call, 2 = 233 ms (1.9×), 4 = 120 ms (3.7×), 8 = 70 ms (**6.4×**).
- **Net result**: at the paper config (R = 200 L, N_init = 20 000) mean wall-clock dropped from **1524 ms/gen → 336 ms/gen** (−77.9 %). Projected 1000-gen total: **~5.6 min on this Mac**, down from ~25–40 min. Steady-state phase mix shifted from 86 % light + 14 % rest to 35 % light + 65 % rest; `intercept` alone is now 19 % of step time.

### Added — Step 21: Per-tree Darwinian dynamics on the Forest

- New package [`mechatree.evolution`](src/mechatree/evolution/) implementing the in-silico Darwinian island the Fortran `EvoluAlgo` was approximating: every tree carries a heritable [`Genome`](src/mechatree/evolution/genome.py), seeds inherit a mutated copy at birth, and lineages disappear naturally when no descendant survives to reproduce. The forest *is* the population — no central selection step.
- [`mechatree.evolution.Genome`](src/mechatree/evolution/genome.py) — frozen dataclass holding 10 safety-NN weights + 18 allocation-NN weights + 3 raw angle genes (decoded to `theta1/theta2/gamma1/gamma2` via the existing `_decode_angles` formula). `Genome.random(rng)`, `Genome.mutate(rng, sigma=0.005, p_locus=0.05)` (per-locus Gaussian, Fortran defaults), `Genome.to_models()` materialises a cached `NeuralSafety` / `NeuralAllocation` pair, `Genome.tree_angles()` returns the radian dict for `replace(cfg.tree, **angles)`.
- [`mechatree.evolution.run_tournament`](src/mechatree/evolution/run.py) — high-level entry point. Returns a `ForestEvolutionResult` (forest, fitness, per-gen history, optional champions path + archive dir). Default fitness = per-tree `n_branches` as a stand-in for lifetime reproductive success.
- [`mechatree.evolution.archive`](src/mechatree/evolution/archive.py) — per-generation JSON snapshot writer/loader (`archive_{gen:08d}.json` with `{generation, n_trees, n_lineages_alive, trees: [{lineage_id, age, n_branches, biomass, location, genome}, ...]}`).
- [`mechatree.evolution.curate`](src/mechatree/evolution/curate.py) — `kmeans2`, `detect_species` lifted verbatim from `scripts/strategies_single_tree.py`; plus `pick_champions(genomes, fitness)`, `from_forest`, `from_archive`, `write`. Output schema matches `data/S3_champions.json` so [`mt.load_champion`](src/mechatree/genome.py) consumes it unchanged.
- [`mechatree.forest.Forest`](src/mechatree/forest.py) gained optional `genomes=[Genome, ...]` + `mutation_sigma` / `mutation_p_locus` kwargs. When `genomes` is provided, `Forest.step` dispatches `requested_growth` / `primary_growth` / the seed-count `compute(...)` call through each tree's own NN pair, threads per-tree angles into `primary_growth`, filters `self.genomes` lock-step with `self.trees` / `self.ages` through the death pass, and mutates the parent's genome at every seedling birth. New `ForestStats.n_lineages_alive` field reports surviving lineages each generation.
- Flat API ([src/mechatree/__init__.py](src/mechatree/__init__.py)) re-exports `mt.Genome`, `mt.run_tournament`, `mt.ForestEvolutionResult`, plus the `mt.evolution` submodule handle so `mt.evolution.curate.from_forest(...)` and `mt.evolution.archive.write_snapshot(...)` work without an extra import. Total flat-API surface now 54 names.
- [`tests/test_evolution_genome.py`](tests/test_evolution_genome.py) (9 cases), [`tests/test_evolution_forest.py`](tests/test_evolution_forest.py) (6), [`tests/test_evolution_curate.py`](tests/test_evolution_curate.py) (7) — pin: random + from_dict round-trip + length validation + bounds clipping; mutation is deterministic per seed and never modifies parent; `Forest(genomes=None)` is byte-identical to Step-12 behaviour; per-tree dispatch produces divergent trajectories; lineage extinction inside 20 gens of a 6-founder run; same-seed reproducibility; archive snapshot round-trip; `pick_champions` argmax per cluster; `from_archive` → `load_champion` round-trip.
- [`benchmarks/bench_evolution.py`](benchmarks/bench_evolution.py) — 20-founder, 100-cap, 50-gen Darwinian island runs in ~470 ms total (~9 ms/gen) on the dev mac; mutation overhead is **2 %** of total run time, `to_models` (NN materialisation, cached) is 0.2 %.
- [`examples/evolve_a_forest.py`](examples/evolve_a_forest.py) — CLI recipe: 50-gen tournament from random founders, writes `champions.json`, prints surviving species, plots `n_lineages_alive` + `n_trees` over time.

### Changed — `scripts/strategies_single_tree.py` shares the curate primitives

- The local `kmeans2` + `detect_species` implementations in [`scripts/strategies_single_tree.py`](scripts/strategies_single_tree.py) were replaced with `from mechatree.evolution.curate import ...`. CLI surface (`--build`, `--dat`, `--json`, `--save-dir`, `--no-show`) is unchanged; the script still builds `data/S3_champions.json` from `S3.dat` the same way it did before.

### Added — Step 23: Flat top-level API (`import mechatree as mt`)

- [`src/mechatree/__init__.py`](src/mechatree/__init__.py) re-exports the curated public surface so users can call `mt.load_config(...)`, `mt.grow_tree(...)`, `mt.Forest(...)`, `mt.plot_tree_3d(...)`, `mt.figstyle.apply()`, `mt.horton_ratios(...)`, etc. without remembering each subpackage path. 49 names total across `config` / `simulate` / `forest` / `light` / `genome` / `plotting` / `stats`, plus the `figstyle` submodule handle.
- DendroFlow surface (`mt.BranchWindBridge`, `mt.make_dendroflow_wind_fn`, `mt.DendroFlowWindParams`) resolved lazily via module-level `__getattr__` — a bare install never tries to import DendroFlow, and asking for the symbol when the extra is missing raises `AttributeError` with the `pip install 'mechatree[dendroflow]'` hint.
- [`tests/test_public_api.py`](tests/test_public_api.py): 7 cases pinning the contract — all `mt.__all__` names resolve, the expected set is present, private helpers (`_callback_arity`, `_resolve_wind_fn`, `_decode_angles`, `_champion`) stay buried, end-to-end `grow_tree` through the flat path works, `mt.figstyle` is the live submodule, deep imports still resolve to the same objects, and `__version__` is exposed.
- [`src/mechatree/genome.py`](src/mechatree/genome.py): `champion_angles` added to `__all__` (the helper was already importable but missing from the list).

### Changed — Notebooks 01–06 + recipe examples use the flat API

- All six notebooks ([01_grow_one_tree.ipynb](notebooks/01_grow_one_tree.ipynb), [02_forest_under_wind.ipynb](notebooks/02_forest_under_wind.ipynb), [03_neural_genome.ipynb](notebooks/03_neural_genome.ipynb), [04_custom_growth_law.ipynb](notebooks/04_custom_growth_law.ipynb), [05_strahler_diagnostics.ipynb](notebooks/05_strahler_diagnostics.ipynb), [06_fractal_dimension.ipynb](notebooks/06_fractal_dimension.ipynb)) now lead with `import mechatree as mt` and call every public function through `mt.*` (`mt.load_config`, `mt.grow_tree`, `mt.Forest`, `mt.plot_tree_3d`, `mt.figstyle.apply()`, `mt.horton_summary`, …).
- Nine recipe-style examples converted to the same form: [grow_one_tree.py](examples/grow_one_tree.py), [grow_a_forest.py](examples/grow_a_forest.py), [custom_simulation.py](examples/custom_simulation.py), [plot_strahler.py](examples/plot_strahler.py), [plot_allocation.py](examples/plot_allocation.py), [plot_self_thinning.py](examples/plot_self_thinning.py), [sympy_genome.py](examples/sympy_genome.py), [dendroflow_wind.py](examples/dendroflow_wind.py), [figstyle_compare.py](examples/figstyle_compare.py).
- Historical/coral-growth examples ([random_growth.py](examples/random_growth.py), [self_avoiding.py](examples/self_avoiding.py), [sap_transport.py](examples/sap_transport.py), [render_blender.py](examples/render_blender.py)) flipped to `mt.PyTree` / `mt.plot_2d` / `mt.plot_3d` / `mt.grow_tree`, retaining deep imports for symbols outside the curated set (`mechatree.geometry.distance_test`, `mechatree.export.to_blender_script`).
- Existing deep imports continue to work — `from mechatree.simulate import grow_tree` resolves to the same callable as `mt.grow_tree` (asserted in tests).

### Added — Step 20: Mesh3d cylindrical tree renderer

- [`mechatree.plotting.plot_tree_3d`](src/mechatree/plotting/_mechanics.py) gained a `style="cylinders"` opt-in path that extrudes each branch as a tapered `Mesh3d` cylinder — bottom radius from the parent's diameter, top radius from the branch's own diameter — so the trunk-to-canopy taper reads volumetrically instead of the previous pixel-width `Scatter3d` line bucketed by Strahler order. Default `style="lines"` is unchanged. Sides per cylinder configurable via `n_sides` (default 8). One `Mesh3d` trace per Strahler order, coloured via the active `figstyle.strahler_color` palette.
- New private helper `_cylinder_mesh(starts, ends, r_base, r_top, *, n_sides=8)` in [`src/mechatree/plotting/_mechanics.py`](src/mechatree/plotting/_mechanics.py) — vectorised, returns plotly-ready `(vertices, i, j, k)` arrays. Orthonormal-frame construction picks the world axis least parallel to the tangent so the cross product stays well-conditioned for any branch orientation.
- [`notebooks/06_fractal_dimension.ipynb`](notebooks/06_fractal_dimension.ipynb) panel (a) now renders the mature champion with `style="cylinders"`; the previous `width_max=18.0` workaround and the comment explaining cylinders weren't available were dropped.
- [`tests/test_plotting.py`](tests/test_plotting.py): five new cases — `style="cylinders"` returns a figure containing at least one `Mesh3d` trace, empty-tree path works, `style="bogus"` raises `ValueError`, and direct `_cylinder_mesh` unit tests pin vertex count, ring radii / z-coordinates, face-index ranges, and the empty-batch return shape.

### Changed — Champion genome loader now returns coding angles too (**breaking**)

- [`mechatree.genome.load_champion`](src/mechatree/genome.py) now returns a 4-tuple `(safety, allocation, angles, non_coding)` instead of the previous 3-tuple `(safety, allocation, meta)`. The third element `angles` is a dict `{theta1, theta2, gamma1, gamma2}` (radians) decoded from `full_row[6:9]` via the Fortran formula in `mod_tree.f90:108-111`. Previously these angle genes were silently discarded and the simulator fell back to the YAML's `tree.theta*` / `tree.gamma*` defaults, so champion trees grew with the wrong branching geometry. With the fix, the species-0 S3 champion at gen 400 reaches ~9 600 branches (was ~7 300), in line with the paper's ~10⁴.
- The fourth element renamed `meta` → `non_coding` (biology metaphor: the bookkeeping fields `species_id` / `centroid_tag` / `champion_index` etc. are the non-coding part of the champion record; angles + NN weights are the coding part).
- [`mechatree.genome.load_all_champions`](src/mechatree/genome.py) returns 4-tuples per species accordingly.
- [`mechatree.genome.models_from_config`](src/mechatree/genome.py) now returns 3-tuple `(safety, allocation, angles)`. For non-Neural paths (`Constant*`, SymPy callbacks) `angles` is `None`; for `neural_from:` it carries the champion's decoded angles.
- [`mechatree.simulate.grow_tree`](src/mechatree/simulate.py) and [`mechatree.forest.Forest`](src/mechatree/forest.py): when `models_from_config` returns non-None angles (YAML neural_from path), the simulator auto-applies them to `cfg.tree` before running, so champion runs use the champion's geometry by default. Explicit `safety=`/`allocation=` callers still need `replace(cfg.tree, **angles)` themselves.
- [`mechatree.genome.champion_angles`](src/mechatree/genome.py): new lightweight helper that decodes the four angles without instantiating the NN models. Equivalent to `load_champion(path, species_id)[2]`.
- [`examples/forest.yaml`](examples/forest.yaml): `light.leaf_transparency: 0.5` made explicit (was previously implicit via `LightConfig` default). Comment on the `tree.theta*` / `gamma*` block pointing at the champion-override pattern.
- [`tests/test_genome_neural.py`](tests/test_genome_neural.py): `test_load_champion_returns_models_angles_and_non_coding` asserts the 4-tuple contract, the Fortran-formula sign convention (`theta1 > 0`, `theta2 < 0`, `gamma1 == gamma2`), and that `load_all_champions` yields a 4-tuple per species.
- All callers updated: [`tests/test_simulate.py`](tests/test_simulate.py), [`tests/test_stats.py`](tests/test_stats.py), [`tests/test_sympy_genome.py`](tests/test_sympy_genome.py), [`examples/figstyle_compare.py`](examples/figstyle_compare.py), [`notebooks/03_neural_genome.ipynb`](notebooks/03_neural_genome.ipynb), [`notebooks/06_fractal_dimension.ipynb`](notebooks/06_fractal_dimension.ipynb).

### Added — `horton_strahler_counts` (paper-faithful # branches per rank)

- [`mechatree.stats.horton_strahler_counts`](src/mechatree/stats.py) — port of the Fortran ``NsegmentsS`` array writer in [`legacy_fortran/mod_tree.f90:1052-1070`](legacy_fortran/mod_tree.f90), the data behind the SI Fig. S8(b) "number of branches per Strahler order" trace. Implements the classical Horton-Strahler stream count: ``N(1) = number of leaves``; ``N(w+1) = number of internal branches whose two children share Strahler order w``. Unlike `strahler_summary.n_branches` (per-unit-segment count, which double-counts long order-W chains along the trunk) and the C++ `horton_summary.n_branches` (a buggy "inheriting child" rule that mis-allocates leaves above bifurcations), this matches the paper's ratios cleanly: ``N(w) / N(w+1) ≈ R_n ≈ 3.5`` across consecutive ranks for the species-0 S3 champion at gen 400.
- [`mechatree.stats.horton_ratios`](src/mechatree/stats.py) gained an `n_branches_override` kwarg (parallel to the existing `mean_length_override`). Used by [`notebooks/06_fractal_dimension.ipynb`](notebooks/06_fractal_dimension.ipynb) so the panel-(b) `R_n` fit and markers come from the same series (`horton_strahler_counts`) — previously the fit anchored on `summary.n_branches` (Horton stream count) while the markers plotted `strahler_summary.n_branches`, so the dashed reference line floated above the data points.
- [`tests/test_stats.py`](tests/test_stats.py) — three new cases: `horton_strahler_counts_perfect_binary_tree` (depth-3 tree → `[8, 4, 2, 1]`), `_n1_equals_leaves_on_grown_tree` (the invariant the paper relies on, and the simulator now passes), `_empty_tree_returns_zero_length_array` (smoke).

### Changed — Notebook 06 panel (b) wired to `horton_strahler_counts`

- [`notebooks/06_fractal_dimension.ipynb`](notebooks/06_fractal_dimension.ipynb): the "# branches" trace and the `R_n` fit both consume `horton_strahler_counts(tree)` (= the Fortran `NsegmentsS`). Cell 3's `snapshot_cb` now stores the array under `"n_w"` per snapshot so the 25-yr open markers + panel-(c) time evolution use the same metric. The `tree.collapse_single_child_chains()` workaround that I added in the previous round to make `strahler_summary.n_branches[0]` match leaves is removed — the new count is correct regardless of chain artifacts.

### Changed — Figstyle axis title closer to ticklabels

- [`src/mechatree/plotting/figstyle.py`](src/mechatree/plotting/figstyle.py) `_axis_style()`: set `title.standoff = 6` (px). Plotly's default standoff is ~20 px which leaves the axis title visually detached from the ticklabels; 6 px matches the SoftMobility / MATLAB compact layout.

### Added — Step 22: Unified figure style (plotly, MATLAB look)

- [`mechatree.plotting.figstyle`](src/mechatree/plotting/figstyle.py) — single-file plotly counterpart of [`SoftMobility/softmobility/classes/figstyle.py`](../SoftMobility/softmobility/classes/figstyle.py). Registers a `go.layout.Template` under `pio.templates["mechatree"]` with white background, 4-sided `mirror=True` black frame, `ticks="inside"` (true MATLAB default), Helvetica 11 pt, `colorway` from `COLORS`, and an orthographic 3D scene. Public API: `apply()`, `figure(size, aspect)`, `subplots(...)`, `figure_3d(...)`, `save(fig, name)`, `set_strahler_cmap(name)`, `strahler_color(order)`, plus the `COLORS` / `SIZES` / `FONT` dicts. Same names as the SoftMobility module so the muscle-memory transfers.
- Four 10-stop Strahler colormaps sampled from MATLAB's `colormap(name)`: `"jet"` (default, matching the `colormap('jet')` line in `../Eloy2017_NatComm_archive/plot_stat_single_tree.m:58`), `"cool"` (literal match for the `colormap('cool')` line at `:37`), `"parula"` (post-R2014b default), `"rainbow"` (legacy `_palette.RAINBOW_STRAHLER`). Pick via `set_strahler_cmap("...")`.
- [`examples/figstyle_compare.py`](examples/figstyle_compare.py) — A/B benchmark script. Renders the species-0 S3 champion under all four Strahler palettes, then the same line plot under four font stacks (Helvetica / Arial / Computer Modern / system-ui), then the same plot under all four (ticks-inside vs ticks-outside) × (4-sided frame vs no top/right) combinations. `uv run python examples/figstyle_compare.py` opens three browser tabs.
- [`tests/test_figstyle.py`](tests/test_figstyle.py) — 7 cases covering template registration, sized canvas dimensions, axis attributes (`mirror=True`, `ticks="inside"`, `showgrid=False`), `COLORS` hex validation, the Strahler cmap switch, and unknown-cmap error.
- New "Strahler palette benchmark" cell in [`notebooks/06_fractal_dimension.ipynb`](notebooks/06_fractal_dimension.ipynb): renders the mature champion under all four palettes in a 1×4 plotly subplot grid.

### Changed — Step 22

- All seven plotly helpers in [`src/mechatree/plotting/`](src/mechatree/plotting/) (`plot_2d`, `plot_3d`, `plot_tree_3d`, `plot_forest_topdown`, `plot_self_thinning`, `plot_allocation`, `plot_strahler_diagnostics`) now build their figures through `figstyle.figure*` and pull colors from `figstyle.COLORS` instead of named CSS strings (`"forestgreen"`, `"saddlebrown"`, `"magenta"`, `"cyan"`, etc.) — the per-figure `update_layout(paper_bgcolor="white", ...)` boilerplate is gone. Public signatures unchanged except `plot_tree_3d`'s `leaf_color` argument defaults to `None` (resolved to `figstyle.COLORS["green"]` at call time).
- [`examples/grow_a_forest.py`](examples/grow_a_forest.py) inline twin-axis figure rewritten through `figstyle.figure`.
- All 6 forester notebooks call `figstyle.apply()` in their first code cell; notebooks 02, 03, 04, and 06 also swap their direct `plotly.subplots.make_subplots` calls for `figstyle.subplots(...)` to inherit the template + sizing.

### Fixed — CI lint job

- [`pyproject.toml`](pyproject.toml): `[tool.uv.sources].dendroflow` was unconditional, so `uv pip install -e ".[dev]"` on the GH Actions runner tried (and failed) to resolve the editable path `../DendroFlow` (which does not exist on CI). Gated on `extra = "dendroflow"` so the sibling-checkout source is only consulted when the optional `dendroflow` extra is requested. Local-developer workflow `uv pip install -e ".[dendroflow]"` is unchanged.

### Changed — `volume_ratio_leaf` default 8.0 → 4.0

- [`mechatree.config.TreeConfig.volume_ratio_leaf`](src/mechatree/config.py) default flipped from `8.0` to `4.0`, matching the Nat Commun reference `Forest_reference.ini` (kept locally under `~/Documents/Arbres/FORTRAN/ART_Revision2/` and `…/ART_Revision2b/`) and the SI Fig. S12 caption (`V_prod = 4 V_0 l`). The legacy `legacy_fortran/Forest.ini` (the older PNAS-submission tarball) still has `VolumeRatioLeaf = 8.0d0`, which is what prior MechaTree code mirrored; that file is unedited per the `legacy_fortran/` provenance rule. Likely fixes a long-standing branch-count over-production: at 200 yr with the species-0 S3 champion, branch count drops from ~77k (VRL=8) to ~8k (VRL=4), in line with the paper's ~10⁴.
- [`examples/forest.yaml`](examples/forest.yaml) explicit value updated to 4.0.
- [`tests/test_config.py`](tests/test_config.py) default assertions updated to 4.0.

### Fixed — Notebook 06 `<l>` mean length

- [`notebooks/06_fractal_dimension.ipynb`](notebooks/06_fractal_dimension.ipynb) was plotting `summary.mean_length` as the panel-(b) `<l>` curve and feeding it to the `R_l` fit. `HortonSummary.mean_length` is the **per-Horton-stream chain length** (sum of unit segments in the stream), not the per-branch recursive distance-to-leaves the paper plots in SI Fig. S8(b). The notebook now uses two new helpers in [`mechatree.stats`](src/mechatree/stats.py) instead:
  - `distance_to_leaves(tree)` — length-aware per-branch recursive distance: terminal branches get `length / 2`; internal branches get `length + nb_leaves`-weighted mean of children's distance. Mirrors `b%distance_leaves` in `legacy_fortran/mod_tree.f90:1174-1203` (`save_area`), generalized from the Fortran's `+1.0` constant to use the actual branch length so the metric stays correct when post-pruning chain-merger (`Tree::removeBranches` at [tree.cpp:502](src/mechatree/_core/tree.cpp#L502)) grows branches past unit length.
  - `mean_distance_to_leaves(tree)` — per-Horton-rank mean of the above. Drop-in for `summary.mean_length` when reproducing the paper figure.
- [`mechatree.stats.horton_ratios`](src/mechatree/stats.py) — new `mean_length_override: np.ndarray | None` kwarg swaps in an alternative length series for the `R_l` fit while leaving the count / diameter / area fits untouched. Notebook 06 passes `mean_distance_to_leaves(tree)` through it. Also gained a `max_rank: int | None` kwarg that caps the fit support; notebook 06 sets it to `7` to mirror Eloy et al. 2017 SI Fig. S12, which fits only the lower ranks because the highest ranks sit on too few branches to be stable.
- [`tests/test_stats.py`](tests/test_stats.py) — six new cases: trunk-only, single binary split (`[1.5, 0.5, 0.5]`), depth-3 perfect binary, unbalanced subtree (verifies the `nb_leaves` weighting on uneven subtrees), per-Horton-rank aggregation, and the `mean_length_override` override path through `horton_ratios`.

### Added — Step 19: leaf transparency in the light model

- [`mechatree.config.LightConfig`](src/mechatree/config.py) — new `leaf_transparency` field (alias `tau`, in `[0, 1]`). Surfaced via YAML `light:` block; defaults to `0.5` per Eloy et al. (Nat Commun 2017).
- [`mechatree.light.intercept`](src/mechatree/light/interception.py) — replaces the binary topmost-wins write with a depth-rank-based `tau ** depth` write. The i-th leaf from the top of a shadow cell receives `tau**i` of the incident light. `tau = 0` recovers the Fortran binary regime (since NumPy follows IEEE `0**0 = 1`); `tau = 1` makes leaves fully transparent. Optional third argument with the same default; threaded through from `LightConfig` in `grow_tree` and `Forest.step`.
- [`tests/test_light.py`](tests/test_light.py) — four new tests cover `tau = 0` / `tau = 0.5` (default) / `tau = 1` and the geometric attenuation on a 4-leaf stack. Two pre-existing cell-binning tests updated to the new default.

### Changed — Step 19

- Default light interception switches from binary topmost-wins to graded `tau = 0.5`. Quantitative growth output under the default config shifts because understorey leaves now contribute photosynthate.

### Added — Notebook 06 (FigS8) + Horton stream stats

- [`notebooks/06_fractal_dimension.ipynb`](notebooks/06_fractal_dimension.ipynb) — reproduces SI Fig. S8 of Eloy et al. 2017 from the species-0 S3 champion. Five panels: 3D mature tree (400 yr), Horton ratios `R_n, R_l, R_d, R_a` and fractal dimension `D = log R_n / log R_l` (open markers = 25 yr juvenile, filled = 400 yr mature), time evolution of those ratios from 0–500 yr, per-segment branch-tapering scatter with 2/3 / 1 / 3/2 reference slopes, and per-bifurcation area-conservation scatter (Leonardo's rule, mean ≈ 0.95).
- [`mechatree.stats.horton_summary`](src/mechatree/stats.py) + [`HortonSummary`](src/mechatree/stats.py) — per-Horton-stream length / count / diameter / area, mirroring the existing per-segment `strahler_summary`. Needed because every MechaTree segment is a unit twig, so per-segment mean length is identically 1; only the per-Horton-chain length grows with rank.
- [`mechatree.stats.horton_ratios`](src/mechatree/stats.py) + [`HortonRatios`](src/mechatree/stats.py) — log-linear fit across consecutive Strahler / Horton ranks recovering the four geometric ratios and the fractal dimension. Mirrors `Fractal_dim.m` from the Eloy et al. 2017 MATLAB archive. `drop_top=True` by default to exclude the trunk (where `N_W=1` flattens the slope).

### Fixed — Step 13 follow-up

- [`mechatree.stats.horton_summary`](src/mechatree/stats.py) forces `tree.set_strahler()` before `tree.set_horton()`. The C++ `Tree::setHorton` only calls `setStrahler` when `Strahler_distribution` is empty, so once Strahler is computed it never refreshes — every `on_step` snapshot during `grow_tree` saw stale rank-1 / rank-2 labels and collapsed to `max_order == 2`. Regression test at [`tests/test_stats.py::test_horton_summary_refreshes_strahler_on_grown_tree`](tests/test_stats.py).

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
