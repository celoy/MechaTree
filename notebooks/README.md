# Notebooks

Forester-facing tutorials that walk through the science alongside runnable
code. Each notebook is a literate-programming companion to a script in
[../examples/](../examples/); the science prose lives in the markdown
cells, the simulation calls live in the code cells.

All figures render inline through plotly. The library has no matplotlib
runtime dependency.

## Install

```bash
uv pip install -e ".[notebooks]"
```

That extra brings in `jupyter`, `ipykernel`, and `nbstripout`. The
`nbstripout` pre-commit hook strips output cells on commit, so the
checked-in `.ipynb` files are clean — re-execute to see plots.

## Launch

```bash
uv run jupyter lab
```

Or run a single notebook end-to-end without opening the UI:

```bash
uv run jupyter nbconvert --to notebook --execute --inplace \
    notebooks/01_grow_one_tree.ipynb
```

## Index

| Notebook | Companion example | What it shows |
| --- | --- | --- |
| [01_grow_one_tree.ipynb](01_grow_one_tree.ipynb) | [grow_one_tree.py](../examples/grow_one_tree.py) | Load a YAML config, grow one tree, render in 3D with Strahler colouring and leaf overlay. |
| [02_forest_under_wind.ipynb](02_forest_under_wind.ipynb) | [grow_a_forest.py](../examples/grow_a_forest.py) | Drive a forest, watch population/biomass dynamics, see Yoda's −3/2 self-thinning emerge. |
| [03_neural_genome.ipynb](03_neural_genome.ipynb) | [grow_one_tree.py](../examples/grow_one_tree.py) | Compare the default constant genome against the two S3 champions from Eloy et al. 2017. |
| [04_custom_growth_law.ipynb](04_custom_growth_law.ipynb) | [custom_simulation.py](../examples/custom_simulation.py) | Plug in user-supplied `wind_fn`, `Sun`, `ConstantSafety`, `ConstantAllocation`. |
| [05_strahler_diagnostics.ipynb](05_strahler_diagnostics.ipynb) | [plot_strahler.py](../examples/plot_strahler.py) | Grow a mature tree, read off Horton-Strahler scaling, Leonardo's rule, Tokunaga matrix. |
