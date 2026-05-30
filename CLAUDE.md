# MechaTree — guide for Claude sessions

## What this project is

MechaTree is a Python (with Cython + C++) port of a Fortran90 simulator of **tree growth in response to wind and light**. The model captures self-pruning, allometric scaling, and self-similar tree architectures that emerge from the mechanical and photosynthetic constraints. It accompanies:

> Eloy, C., Fournier, M., Lacointe, A. et al. *Wind loads and competition for light sculpt trees into self-similar structures.* **Nat Commun 8, 1014 (2017).** https://doi.org/10.1038/s41467-017-00995-6

## Lineage

1. **Fortran90 reference code** — written by the user, simulations published in Nat Commun (2017). Sources preserved in `legacy/fortran/`.
2. **Diego Bengochea Paz's Python port (2017)** — Python + Cython + C++ library that handles tree architectures, written by Diego Bengochea Paz (ORCID [0000-0002-0835-3981](https://orcid.org/0000-0002-0835-3981)) during a 2017 internship. Lightly touched up Feb 2025. Preserved verbatim in `legacy/pytree/`.
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
| `legacy/` | Reference only | All 2017-paper provenance in one place: `legacy/fortran/` (Fortran90 sources), `legacy/matlab/` (analysis scripts), `legacy/pdf/` (paper + SI), `legacy/pytree/` (Diego's Python+Cython+C++ port, was `archive/CommentedLibrary_and_Doc/`). Do not edit. Do not import from. |

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
- **Provenance**: keep `legacy/` intact — it is intellectual history.

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

- Don't edit anything under `legacy/`.
- Don't import from `legacy/` in `src/mechatree/` or `tests/`.
- Don't re-introduce the rest of the 15 GB of simulation outputs / zips that live in `/Users/Ch/Documents/Python/Eloy2017_NatComm_archive/` (outside the repo, local-only). Only the MATLAB analysis scripts and the paper + SI PDFs were brought in — under `legacy/matlab/` and `legacy/pdf/`.
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
| 7  | Pre-release shake-out (TestPyPI)    | ✅ done    | Version bumped to `0.1.0a1` (pyproject.toml + __init__.py). CHANGELOG updated. [`.github/workflows/wheels.yml`](.github/workflows/wheels.yml) extended with `publish-testpypi` job (OIDC trusted publisher, `workflow_dispatch` only). Ready for manual tag + dispatch. |
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
| 20 | Mesh3d cylindrical tree renderer    | ✅ done    | [`plot_tree_3d`](src/mechatree/plotting/_mechanics.py) gained a `style="cylinders"` path that extrudes each branch as a tapered `Mesh3d` cylinder (parent diameter → own diameter). `style="lines"` remains the fast default. Notebook 06 panel (a) switched over. |
| 21 | Per-tree Darwinian dynamics on the Forest | ✅ done    | [`mechatree.evolution`](src/mechatree/evolution/) — `Genome` dataclass (31 floats: 3 angle genes + 10 NN-safety + 18 NN-allocation, per-locus Gaussian mutation), per-tree dispatch in [`Forest`](src/mechatree/forest.py) when `genomes=[...]` provided, seeds inherit a mutated copy of the parent's genome, and curation lifted from [`scripts/strategies_single_tree.py`](scripts/strategies_single_tree.py) into `mechatree.evolution.curate` writes a `champions.json` consumed unchanged by [`mt.load_champion`](src/mechatree/genome.py). Flat API: `mt.Genome`, `mt.run_tournament`, `mt.evolution`. Per the user's intent, this is the in-silico Darwinian island the Fortran was approximating — *not* a port of NSGA-style external tournament + MPI. Reverses the "Evolution is external" call. |
| 21b | Island-scale prep: profile + batch + parallelise the light pipeline | ✅ done    | Profiled a 1000-yr forest at the Nat Commun config (R = 200 L, N_init = 20 000). Found `light.intercept` + `extract_leaves` dominated wall-clock (86 %). **Phase A** added batched Cython accessors `PyTree.get_leaf_tips_batch` and `set_lights_batch`, cutting per-leaf Python loops in [`extract_leaves`](src/mechatree/light/leaves.py) and [`aggregate_onto_trees`](src/mechatree/light/interception.py). **Phase B** ported the per-direction `intercept` loop to a new C++ kernel ([light.h](src/mechatree/_core/light.h) / [light.cpp](src/mechatree/_core/light.cpp)) wrapped via `_core.light_intercept_kernel`. **Phase C** parallelised that kernel over sun directions with OpenMP (each direction is independent so the fan-out is embarrassingly parallel) and replaced `pow(tau, depth)` with a running multiply. Net: 1000-gen × 20 000-tree forest projects to **~5.6 min** (was ~25–40 min). |
| 21c | Full-scale Nat Commun tournament reproduction | ✅ done | Checkpoint/resume infrastructure + parallel replicates runner. `mechatree.evolution.run.run_tournament` now accepts `resume_from` parameter; [`scripts/run_tournament.py`](scripts/run_tournament.py) drives N independent replicates via `multiprocessing.Pool`, combines survivors, and curates champions. [`scripts/validate_tournament.py`](scripts/validate_tournament.py) validates output against `data/S3_champions.json` (required: ≥1 species, centroid gap > 0.15; optional: match Fortran species count, centroid proximity). Tournament-scale YAML: [`examples/tournament_natcomm.yaml`](examples/tournament_natcomm.yaml) (R=200 L, N_init=20000, n_gens=100000, matches paper parameters). Ready for full 100k-gen run. |
| 22 | Unified figure style                | ✅ done    | [`mechatree.plotting.figstyle`](src/mechatree/plotting/figstyle.py) registers a MATLAB-look plotly template (white bg, 4-sided box, inside ticks, Helvetica 11 pt, palette in `COLORS`, sized canvases in `SIZES`) under `pio.templates["mechatree"]`. All seven plotting helpers + the 6 forester notebooks + [examples/grow_a_forest.py](examples/grow_a_forest.py) route through it; `figstyle.figure` / `figstyle.subplots` / `figstyle.figure_3d` / `figstyle.save` mirror the SoftMobility API. Ships four Strahler colormaps (`jet`, `cool`, `parula`, `rainbow`) selectable via `set_strahler_cmap`; notebook 06 renders a 4-panel benchmark so the eventual default can be picked visually. CI lint also fixed (gated `[tool.uv.sources].dendroflow` on `extra = "dendroflow"`). |
| 23 | Flat top-level API (`import mechatree as mt`) | ✅ done    | [`src/mechatree/__init__.py`](src/mechatree/__init__.py) re-exports the curated public surface so notebooks/users do `import mechatree as mt; mt.load_config(...)`, `mt.grow_tree(...)`, `mt.plot_tree_3d(...)`, `mt.figstyle.apply()`, etc. Notebooks 01–06 + all recipe examples swept to use `mt.*`. DendroFlow surface (`mt.BranchWindBridge`, `mt.make_dendroflow_wind_fn`, `mt.DendroFlowWindParams`) is gated via lazy `__getattr__`. Deep imports keep working for back-compat. |
| 25 | Tunable wind statistics + native bulk-thinning + storm replay + explainer notebook | ✅ done    | Unified the scattered wind story: a SymPy-driven [`Distribution`](src/mechatree/wind/distributions.py) inverse-CDF sampler (`amplitude_cdf` / `angle_cdf` YAML knobs, Fortran defaults preserved), a self-contained native bulk-thinning model that rotates the canopy to face the storm, a pass-by-pass [storm-replay diagnostic](src/mechatree/wind/replay.py) + plot, and the [notebook 07](notebooks/07_wind_models.ipynb) explainer. (Native/dendroflow models later removed in Step 26.) |
| 24 | Coupled wind ↔ pruning fixed-point loop | ✅ done    | `wind → prune` iterates to a fixed point when the wind is canopy-aware (3-arg), exiting on zero new cuts, a sub-`ε_rel` canopy-mean change, or the `max_pruning_iterations` cap (knobs on [`WindConfig`](src/mechatree/config.py); cap=1 recovers single-pass byte-identically, 2-arg wind is a no-op). Forest-wide loop; new `n_wind_iterations` stat. Gated/sped by a batched Cython accessor [`PyTree.get_branch_data_batch`](src/mechatree/_core/_core.pyx) (5–7× faster canopy pool). |
| 25c | Option-B C++ plumbing: per-branch F_seg → prune | ✅ done | Per-branch CFD force now reaches pruning instead of a canopy-mean re-broadcast. Two new `Branch` fields (`segment_force_`, `segment_wind_`) written by the actuator-disk bridge, read by a new `prune_with_stored_forces(tree, leaf_drag_S0, cauchy)` C++ entry that skips the per-branch `wind_force` recompute and uses each branch's own local wind for the leaf-cluster drag too. Opt-in via `wind.actuator_disk_use_per_branch_forces`. Verified equivalent to legacy `prune` when fed the canopy-mean force, and the A/B bench shows the real per-branch divergence (exposed crowns prune more). See Step 25c notes below. |
| 25d | Forest-scale perf benchmark + per-phase profile | ✅ done    | Realistic forest-step profile with canopy-aware wind. Script: [`benchmarks/bench_actuator_disk_at_scale.py`](benchmarks/bench_actuator_disk_at_scale.py). Result at the requested paper-scale probe (R = 100, n_trees_init = 1000, 100 gens): **total wall-clock 2.8 s** (28 ms/gen mean). **Wind kernel dominates 93 %** of wall-clock; light / mechanics / growth / pruning combined are < 5 %. Fixed-point iter count is well-behaved: **mean 2.23 iters/gen, median 2, max 7** out of the cap-8. Histogram: `{1 iter: 2 gens, 2 iters: 77 gens, 3 iters: 20 gens, 7 iters: 1 gen}` — one big-storm outlier. Net of the actuator-disk kernel: 223 calls across 100 gens, average 11.7 ms/call. Most trees die early in the run (1000 init → 61 alive at gen 100, 1620 final branches) because the initial canopy hasn't adapted yet — open question whether this matches the Nat Commun reference; defer to Step 21c validation. **Optimisation pointer**: since the wind kernel is 93 % of wall-clock, the next perf step is at the kernel level (Numba/OpenMP/SIMD on the column march), not on the rest of `Forest.step`. |
| 25e | (nu_diff, H) sensitivity sweep | ✅ done    | Originally requested as "(l_mix, H) sweep to see the impact on profiles". After the Occam's razor pass collapsed k-ε's `(I_turb, l_mix, C_μ)` into the single `nu_diff` knob, this is the (nu_diff, H) sensitivity sweep. Script: included in [`bench_actuator_disk_at_scale.py`](benchmarks/bench_actuator_disk_at_scale.py) sweep section (5 × 3 grid on a R = 50, 218-tree, 50k-branch warmup forest). Findings: **H is the cost knob** (single call: H = 0.5 → ~57 ms, H = 1.0 → ~13 ms, H = 2.0 → ~7 ms — a 60 % drop per H doubling; the 4× cell count → 4× compute scaling). **`nu_diff` is the wake-shape knob**: at H = 1 the min `U_branch/U_∞` rises from 0.018 (`nu_diff = 0`, sharp wake) to 0.197 (`nu_diff = 0.3`, spread wake) — diffusion fills in the deficit. ⟨U_branch/U_∞⟩ varies less (0.51–0.62 across the sweep). `nu_diff` doesn't change wall-clock noticeably (same arithmetic, just a different coefficient). Recommended defaults: `H = 1.0` for canopy at the unit-twig scale; `nu_diff = 0.03` as a moderate-diffusion baseline. The wake-FWHM column in the sweep table runs into the lateral boundary at small H — should re-bench on a larger pad_y if a clean FWHM measurement is needed. |
| 25b | Native 3-D actuator-disk wind + Nat Comms force kernel + speed pass | ✅ done    | Replaces the 1-D bulk-thinning bridge as the canonical canopy-aware wind path. **(a)** [`mechatree.wind._actuator_disk_kernel`](src/mechatree/wind/_actuator_disk_kernel.py) — pure-NumPy column-march CFD: per-cell drag from the Nat Commun vector kernel `F_seg = ½ρU²DL ‖t×u‖² (n×t)` (simplifies to `F_D = ½ C_D D L sin³i · U²` streamwise with `ρ = 1`), per-cell actuator-disk update `U_out = ½U_in[1 + √(1 - 4F_D/(H²U_in²))]` (the +-root of the §1.2 quadratic in Loukas's PDF), implicit 4-neighbour cross-stream diffusion with a single `nu_diff` knob that subsumes the legacy k-ε closure `(I_turb, l_mix, C_μ)` after the algebra collapses to `ν_t = α U_∞ I_turb l_mix`. Skips DendroFlow's cylinder subdivision (each MechaTree twig fits in one cell with `H ≥ 1`). Per-branch `U_loc / F_N / F_D / F_vec` arrays exposed on `ActuatorDiskResult` for Option-B propagation. **(b)** [`mechatree.wind.actuator_disk.ActuatorDiskWindBridge`](src/mechatree/wind/actuator_disk.py) — 3-arg `WindFn`: rotates canopy to face storm, builds grid from bounding box, log-law inflow `U(z) = (u_a/κ) log(z/z_0)`, returns canopy-mean wind in world frame. **(c)** New `WindConfig.model = "actuator_disk"` + fields (`actuator_disk_nu_diff / pad_x / pad_y / pad_z / ua / z0 / kappa`) wired through `_resolve_wind_fn`. `examples/actuator_disk_wind.yaml` is the YAML template. **(d)** Derivation: [`docs/actuator_disk_derivation.md`](docs/actuator_disk_derivation.md) — Nat Comms force kernel + Loukas's §1.1-§1.3 conservation-law derivation of the actuator-disk formula, with the eq-19 typo (missing `U_in²` in the denominator) flagged + corrected. **(e)** Bench: [`benchmarks/bench_actuator_disk.py`](benchmarks/bench_actuator_disk.py) — wall-clock vs grid size + forest size. After the `np.pad` rewrite, scratch pre-allocation, and skip-empty-columns: **2.7× faster** than the first cut. Notebook-scale (5k branches, 67k cells) runs in **2.7 ms**; 533k cells in 9.6 ms. **(f)** Step-25 follow-up (b) closed: per-branch wind colourmap in notebook 07 §3.5 now uses the native kernel directly (no DendroFlow dependency). **(g)** Option-B Cython half: [`PyTree.set_forces_batch`](src/mechatree/_core/_core.pyx) writes the per-branch F_vec from the kernel back into the tree's branch.force_ storage, prepping the C++ side. Deferred to Step 25c: the C++ `prune_with_stored_forces` variant + A-vs-B perf compare. |
| 26 | **Momentum-wind overhaul** (rename + consolidate + sensing unification + scale) | ✅ done | Umbrella for the 2026-05-27 rework (sub-steps 26a–26f, all landed). The screening-aware CFD wind is now physically consistent — sensing + pruning share one per-branch screened field, ½ρU² restored — and fast enough for a Darwinian island tournament: legacy scalar-mean models removed, `actuator_disk → momentum` rename, storm fixed per generation, GIL-free C++ kernel + parallel sensing. Validated at island scale (R = 200 / 20 000 founders ≈ 0.24 s/gen `default`, 0.89 s/gen `momentum`). |
| 26a | Storm-fixed-per-gen + drop canopy-mean wind mode | ✅ done | Sample the storm `(θ, amplitude)` **once per generation** and hold it across the fixed-point pruning iterations (today the bridge re-samples every `__call__`, so a "big-storm" gen is really N *different* storms — inflating the iteration count and biasing big-storm pruning). Make per-branch the **only** momentum behaviour; delete the canopy-mean-as-wind path + `use_per_branch_forces` flag (a 3-D solve averaged to one scalar == constant wind; pay CFD cost for nothing). Keep the mean only as the ε convergence thermometer. **Iteration-criterion study deferred to after 26c** — sensing must be screening-consistent before the pruning iteration count is physically meaningful (per user 2026-05-27). **Test deliverable (this step):** storm identical across fixed-point iters; canopy-mean mode gone; suite green; a *preliminary* iteration count (caveated as pre-sensing). |
| 26b | Rename `actuator_disk → momentum`, `H → grid_size` | ✅ done | Mechanical rename across kernel/bridge/config/`model:` string/docs/tests/bench/notebook 07; shared `WindConfig.H → grid_size` (default **2**, the `O(branches)≈O(cells)` knee). Scaling bench confirmed: cost knee at gs=2 (≥5 saves nothing past the branch-work floor), wake preserved at gs=2 vs washed-out by gs≥5; `max(branch length)=1.0` so gs=2 is binning-safe. 401 tests green. Legacy `native`/`dendroflow` bridge *params* kept as `H` (their fate is the open question). See Step 26b notes. |
| 26c | Unify sensing on the momentum-wind per-branch field | ✅ done | Sensing now uses the same screened per-branch field as pruning. New C++ `calculate_stresses_from_stored_forces` (leaves-to-trunk aggregation + `max_stress` from `segment_force_`/`segment_wind_`, no Weibull/cut; `reset_max` flag). The `momentum` sensing sweep solves `n_sensing_angles` directions (from `angle_cdf`) at a **uniform inlet U = 1** — the kernel still screens it down the canopy — and keeps the per-branch max stress; `default` keeps the legacy unit-wind 4-angle sweep. Also surfaced + fixed the legacy `wind_force` missing-½ρU² bug (separate note below). 389 tests green. See Step 26c notes. |
| 26d | Profile + island-scale validation (parallelism: measured, deferred) | ✅ done | Measured before building: the `n_sensing_angles` sensing solves are the bottleneck (68 % of wind), but in-process thread fan-out is 0.67× on the GIL-bound NumPy march, so a GIL-free kernel rewrite was deferred to 26e. Validation gate passed — momentum+sensing is ~linear in branches (≈0.5 s/gen at island scale) and a cold-start screening-aware tournament survives + selects. |
| 26e | GIL-free Cython/OpenMP momentum kernel | ✅ done | Ported the pure-NumPy `compute_momentum_wind` column march to a `nogil` C++ kernel ([momentum.{h,cpp}](src/mechatree/_core/momentum.cpp)) + Cython binding `momentum_solve_kernel`. Equivalence-gate passes (native ≡ NumPy to atol 1e-10; actual drift ~1e-15). **5.8× faster per isolated solve**; the parallelizable sensing sweep (68 % of wind) goes **numpy 43.6 ms → native serial 21.8 ms (2.0×) → native + 4 threads 10.8 ms (4.0×)** at 54 k branches — the 26d thread result (**0.67×**) is flipped. Parallelism is a stdlib `ThreadPoolExecutor` over the `nogil` kernel (no OpenMP dependency — works on any runner), driven by [`MomentumWindBridge.solve_directions`](src/mechatree/wind/momentum_wind.py) + `WindConfig.momentum_sensing_threads`. Test-driven: port → equivalence-test → parallelise → re-measure. See Step 26e notes. |
| 26f | Momentum per-solve overhead reduction (consolidate glue into the kernel) | ✅ done | Pushed rotation / grid build / inflow / canopy-mean / F_vec world-rotation into the C++ kernel — 26e profiling had shown the kernel was only ~29 % of a real solve (5.5 ms at 53 k branches), the rest Python glue. Two entries share one march core: `momentum_solve` extended additively with `cos/sin theta` (F_vec → world frame) + nullable `canopy_mean_out` (pruning path, keeps full grid outputs for `last_result`); new `momentum_solve_world` (sensing — pooled *unrotated* geometry + `theta`, builds the grid internally via an `np.arange` replica, returns only `F_world`/`w_world`, fully `nogil`). **Per-solve: sensing 5.5 → 1.94 ms (2.8×), pruning 5.5 → 2.94 ms (1.9×).** Thread scaling on the sensing sweep is now real (1.4× @2, **1.85× @4** — sublinear because the column march is memory-bandwidth-bound, not GIL-bound). Equivalence-gated (native ≡ NumPy storm-frame oracle; `momentum_solve_world` ≡ the old Python pipeline incl. the `np.arange` replica, atol 1e-10). See Step 26f notes. |

> **Detailed per-step implementation notes** live in [docs/roadmap_history.md](docs/roadmap_history.md); the granular Keep-a-Changelog log is in [CHANGELOG.md](CHANGELOG.md). The table above is the milestone source of truth.

## What's left & future directions

The modernization is functionally **complete** — Steps 1–26 and 21c are all ✅ done. The infrastructure for the Nat Commun tournament is in place. Step 7 (TestPyPI pre-release) is ready. What remains: Step 18 (PyPI stable release) and running the full tournament.

### Pending roadmap steps

- **Step 18 — Stable PyPI 0.1.0 release.** Depends on Step 7 (complete). After a successful TestPyPI dry-run, tag `v0.1.0` (non-alpha), modify [`.github/workflows/wheels.yml`](.github/workflows/wheels.yml) to publish to real PyPI (add a `publish-pypi` job with `repository-url: https://upload.pypi.org/legacy/`), then trigger manually.
- **Step 21c full tournament run** (infrastructure complete, execution pending). Run the 100k-gen tournament at paper scale (R = 200 L, N_init = 20 000, multiple independent replicates) using `python scripts/run_tournament.py --config examples/tournament_natcomm.yaml --n-gens 100000 --n-replicates 3 --out-dir out/natcomm_tournament`. Validation against `data/S3_champions.json` via `python scripts/validate_tournament.py`. At island scale, a single forest-generation takes ~0.24 s (`default` wind) to ~0.89 s (`momentum` wind), so a 100k-gen run is 6.7–24.7 h single-threaded; parallelism across replicates can reduce wall-clock. Deferred: per-Strahler ratio comparison to SI Fig. S8 for deeper validation.

### Small deferred follow-ups

Each is recorded in its closed-step note in [docs/roadmap_history.md](docs/roadmap_history.md); surfaced here so they aren't lost:

- **Coarse-grid diffusion guard (Step 26b):** the momentum kernel's cross-stream diffusion stencil raises `IndexError` when the canopy is shorter than one cell (`Nz < 2`); only bites at very coarse `grid_size`. Clamp to `Nz ≥ 2` or skip z-diffusion when `Nz < 2`.
- **OpenMP-over-trees in pruning (Step 24):** trees are independent given a fixed wind, and the Step-21b build infra (`setup.py` OpenMP probe + serial fallback) is already in place — but per-tree Python overhead may dominate, so profile before building.
- **`default`-model directional sensing (Step 25):** the `momentum` model's sensing sweep already honours `angle_cdf` / `n_sensing_angles` (Step 26c); the `default` model still uses 4 hardcoded angles. Extend the C++ `calculate_stresses` to accept a runtime angle list if directional sensing is wanted under `default` too.
- **Per-branch "stress floor" cache (Step 24):** the fixed-point loop re-rolls the Weibull test on survivors each iteration, which can cut a few branches that would physically survive; a `Branch::p_fail_survived_` cache would remove the artifact. Deferred until a concrete use case demands it.

### Recommended direction (for the user to steer)

1. **Ship it (Steps 7 → 18).** The port is feature-complete; releasing to PyPI is the natural v1 milestone that finally serves the forestry / course audience this project is built for. Low risk.
2. **Scientific validation (Step 21c).** Run the full tournament and validate against Eloy et al. (2017) — this closes the loop with the paper and is the project's reason for existing.
3. **Close the docs gap.** An evolution notebook (Step 21's was deferred) plus a momentum-wind island tutorial would cover the two newest capabilities for end users.

## Out of scope (do not re-litigate unless explicitly asked)

- *(empty for now — the previous "Evolution is external" entry was reinstated as roadmap Step 21 on 2026-05-25 at the user's request.)*

## How to update this file

**Update CLAUDE.md after substantive changes.** When a roadmap step lands:

1. Flip its row in the table from `⬜ pending` / `🚧 wip` to `✅ done`, keeping the Summary cell short (1–2 sentences — the long-form detail goes in step 2).
2. Append a `### Step N notes (closed YYYY-MM-DD)` subsection to [docs/roadmap_history.md](docs/roadmap_history.md) — one short paragraph + bullets, in the shape used for the closed steps there. (Detailed notes no longer live inline in this file.)
3. Add a corresponding entry to [CHANGELOG.md](CHANGELOG.md) under `## [Unreleased]`, using Keep-a-Changelog `### Added` / `### Changed` / `### Fixed` headings.
4. If the step closes or surfaces forward work, update the **What's left & future directions** section above.

For new work introduced mid-stream, add a fresh table row before doing the implementation so the roadmap stays the source of truth. Trivial changes (typos, single-line fixes, formatting-only) don't need an update; anything that adds/changes a module, public function, YAML key, or test file does.
