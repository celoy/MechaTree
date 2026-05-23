# Sphinx documentation Makefile for MechaTree
#
# IMPORTANT: autodoc imports the `mechatree` package, which loads a compiled
# Cython/C++ extension. Run `uv pip install -e ".[dev,docs]"` first.

SPHINXOPTS    ?= -W
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help html clean

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

html:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS)

clean:
	rm -rf "$(BUILDDIR)"
