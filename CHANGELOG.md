# Changelog

All notable changes to this project will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
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
