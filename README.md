# MechaTree

**Mechanical simulation of tree growth in response to wind and light.**

MechaTree is a Python (with Cython + C++) port of the Fortran90 simulator behind the paper:

> Eloy, C., Fournier, M., Lacointe, A. et al. *Wind loads and competition for light sculpt trees into self-similar structures.* **Nat Commun 8, 1014 (2017).** https://doi.org/10.1038/s41467-017-00995-6

The model evolves tree architectures under biomechanical (wind loads, self-supporting limits) and photosynthetic (competition for light) constraints, producing self-pruning behavior, allometric scaling, and self-similar branching.

## Status

Pre-release / port in progress. The repository currently contains:

- `src/mechatree/` — the modern Python package (skeleton; real code lands in Step 2).
- `legacy_fortran/` — the Fortran90 sources from the 2017 paper, for reference.
- `archive/` — Diego Bengochea Paz's 2017 Python+Cython+C++ port, kept for provenance. **Not** imported by the modern package.

## Install

```bash
# Coming in Step 2 (once the Cython extension is wired up).
# Until then, only the import smoke test passes.
```

## Develop

```bash
brew install uv                          # one-time, if needed
uv venv --python 3.12
uv pip install -e ".[dev]"
uv run pre-commit install
uv run pytest
```

See [CLAUDE.md](./CLAUDE.md) for the project's full layout, conventions, and step-by-step modernization plan.

## Acknowledgments

The original Python + Cython + C++ port of the simulator was written in 2017 by **Diego Bengochea Paz** (ORCID [0000-0002-0835-3981](https://orcid.org/0000-0002-0835-3981)) during an internship. That code is preserved under `archive/` for provenance, and the ongoing modernization in `src/mechatree/` is a direct descendant of his work.

## License

[GNU GPL v3 or later](./LICENSE).

## Citation

If you use MechaTree or the underlying model in academic work, please cite:

```bibtex
@article{eloy2017wind,
  title   = {Wind loads and competition for light sculpt trees into self-similar structures},
  author  = {Eloy, Christophe and Fournier, Meriem and Lacointe, Andr{\'e} and Moulia, Bruno},
  journal = {Nature Communications},
  volume  = {8},
  number  = {1},
  pages   = {1014},
  year    = {2017},
  doi     = {10.1038/s41467-017-00995-6}
}
```
