"""Exporters that convert ``PyTree`` / ``Forest`` to external 3D formats.

Currently supported:

- :func:`mechatree.export.blender.to_blender_script` writes a standalone
  Blender Python script that, when run via ``blender --python script.py``,
  constructs the tree geometry, applies materials, and renders.
"""

from mechatree.export.blender import to_blender_script

__all__ = ["to_blender_script"]
