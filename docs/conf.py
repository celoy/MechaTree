from importlib.metadata import version as _v

project = "MechaTree"
author = "Christophe Eloy"
copyright = "2017-2026, Christophe Eloy"
release = _v("mechatree")
version = ".".join(release.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = ".rst"
master_doc = "index"

html_theme = "furo"
html_static_path = ["_static"]
html_theme_options = {
    "source_repository": "https://github.com/celoy/MechaTree/",
    "source_branch": "main",
    "source_directory": "docs/",
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
