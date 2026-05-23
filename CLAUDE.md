# MechaTree — guide for Claude sessions

## What this project is

MechaTree is a Python (with Cython + C++) port of a Fortran90 simulator of **tree growth in response to wind and light**. The model captures self-pruning, allometric scaling, and self-similar tree architectures that emerge from the mechanical and photosynthetic constraints. It accompanies:

> Eloy, C., Fournier, M., Lacointe, A. et al. *Wind loads and competition for light sculpt trees into self-similar structures.* **Nat Commun 8, 1014 (2017).** https://doi.org/10.1038/s41467-017-00995-6

## Lineage

1. **Fortran90 reference code** — written by the user, simulations published in Nat Commun (2017). Sources preserved in `legacy_fortran/`.
2. **Diego Bengochea Paz's Python port (2017)** — Python + Cython + C++ library that handles tree architectures, written by Diego Bengochea Paz (ORCID [0000-0002-0835-3981](https://orcid.org/0000-0002-0835-3981)) during a 2017 internship. Lightly touched up Feb 2025. Preserved verbatim in `archive/`.
3. **Modernization (in progress)** — canonical code lives in `src/mechatree/`. Steps 1–6 complete: scaffold, C++/Cython port, auxiliary modules, tests, Sphinx docs, CI matrix + cibuildwheel. Now entering the simulator-port phase (Steps 8+ below).

## Where code lives

| Path | Status | Rule |
| --- | --- | --- |
| `src/mechatree/_core/` | **Canonical** — compiled tree-architecture core (`PyTree`, branches, Strahler/Horton) | Hot path; C++/Cython |
| `src/mechatree/geometry/` | **Canonical** — geometric helpers (distance, etc.) | Pure Python |
| `src/mechatree/plotting/` | **Canonical** — 2D / 3D matplotlib renderers | Pure Python |
| `examples/` | Active | Recipe-style scripts for users |
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
- Don't re-introduce the 15 GB of simulation outputs / zips / MATLAB scripts that lived in the old `FORESTArticlePNAS/` folder — they now live at `/Users/Ch/Documents/Python/Eloy2017_NatComm_archive/` (outside the repo, local-only).
- Don't be misled by the legacy folder name `FORESTArticlePNAS` — the paper landed at **Nature Communications**, not PNAS. The "PNAS" name reflects an earlier submission target.
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

### Up next

- **Step 7 — First PyPI release.** Tag `v0.1.0`; release the current "tree architecture" library. CHANGELOG note explains scope: this release ships the `PyTree` topology + geometry + plotting; the full mechanical/light simulator follows in subsequent releases.

### Simulator port (replaces what was previously a single TODO)

- **Step 8 — Design doc & API sketch.** Write `docs/design.rst` capturing: the `Leaf` and `Branch` data shapes, the YAML config schema (translated from `legacy_fortran/Forest.ini`), the public API for constructing trees with chosen functionalities, the callback signatures for genome decisions, and the interface between the C++ core and the Python orchestrators. No production code yet — this step exists so the design is reviewed before any C++ is touched.
- **Step 9 — Mechanics + growth in the core.** Port load propagation, stress calculation, primary/secondary growth, and pruning from `legacy_fortran/mod_tree.f90` (`wind_force`, `calculate_stresses`, `requested_growth`, `primary_growth`, `secondary_growth`, `pruning`, `cut_branch`, `new_branches`) into `src/mechatree/_core/`. Expose as methods on `Tree`/`PyTree` *or* as free functions taking a tree — decision deferred to Step 8. Add the genome-callback hook at the Cython boundary. Numerical regression tests against Fortran reference outputs.
- **Step 10 — Light interception module.** New subpackage `src/mechatree/light/` (Python first; promote to Cython if profiling demands). Port `light_interception` and `light_on_trees` from `legacy_fortran/mod_tree.f90`. Public API takes a `Leaves` collection + a sun model; returns light per leaf. Independent of `PyTree`.
- **Step 11 — Single-tree simulator.** Python orchestrator replicating `legacy_fortran/tree.f90`'s main loop: load YAML config, build a tree, iterate {light → stresses → growth → pruning → seeds/leaves → ordering → save}. Ship as `mechatree.simulate.grow_tree(config)`. Regression test vs Fortran `tree.f90` on a fixed seed.
- **Step 12 — Forest container.** Replicate `legacy_fortran/Forest.f90`'s main loop. New `Forest` class — owns N trees, spatial positions, birth/death lifecycle. Cross-tree light competition uses the Step-10 module on the concatenated leaf list. Regression test vs Fortran `Forest.f90`.
- **Step 13 — Forester-facing examples & tutorials.** Jupyter notebooks: "Grow one tree from a YAML config", "Build a forest and watch it self-prune under wind", "Plug in your own growth law". Ship a commented example YAML. Update `docs/userguide.rst`.

### Out of scope (do not re-litigate unless explicitly asked)

- **Evolution** (`legacy_fortran/EvoluAlgo.f90`, `legacy_fortran/mod_evolu.f90`) — deferred, possibly a separate package later. Per the design principle above, evolution is external to the core library.
