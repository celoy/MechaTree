"""Tests for the Blender script generator.

We can't execute the generated script without Blender installed, so the
tests focus on: the script is well-formed Python, contains the expected
geometry payload, and references the right output paths.
"""

import ast
import json

from mechatree.config import Config
from mechatree.export import to_blender_script
from mechatree.simulate import grow_tree


def test_script_is_valid_python(tmp_path):
    """The generated script must parse as Python — otherwise Blender
    can't run it."""
    tree = grow_tree(Config(), n_generations=10, seed=0)
    script = to_blender_script(tree, tmp_path / "out.py")
    code = script.read_text()
    ast.parse(code)


def test_script_carries_all_branches(tmp_path):
    """Every branch in the source tree appears in the script's payload."""
    tree = grow_tree(Config(), n_generations=15, seed=0)
    n_expected = tree.get_number_of_branches()

    script = to_blender_script(tree, tmp_path / "out.py")
    code = script.read_text()

    # Find the embedded PAYLOAD JSON literal.
    start = code.index('PAYLOAD = json.loads(r"""') + len('PAYLOAD = json.loads(r"""')
    end = code.index('""")', start)
    payload = json.loads(code[start:end])
    assert len(payload["branches"]) == n_expected
    # Trunk's location field must be a 3-element list.
    assert len(payload["branches"][0]["location"]) == 3
    # Strahler orders populated.
    assert all(b["strahler"] >= 1 for b in payload["branches"])


def test_script_records_render_and_save_paths(tmp_path):
    tree = grow_tree(Config(), n_generations=5, seed=0)
    render = tmp_path / "render.png"
    blend = tmp_path / "scene.blend"

    script = to_blender_script(
        tree,
        tmp_path / "out.py",
        render_path=render,
        save_blend_path=blend,
        image_resolution=(640, 480),
    )
    code = script.read_text()

    # The paths show up resolved (posix form) in the embedded JSON.
    assert render.resolve().as_posix() in code
    assert blend.resolve().as_posix() in code
    assert '"image_resolution": [\n    640,\n    480\n  ]' in code


def test_script_works_for_a_forest(tmp_path):
    """A Forest's trees are concatenated into the payload."""
    from mechatree.config import ForestConfig
    from mechatree.forest import Forest

    cfg = Config(forest=ForestConfig(size=10.0, n_trees_init=4, n_trees_max=50))
    forest = Forest(cfg, seed=0)
    forest.step(0)  # one step so trees have some structure

    script = to_blender_script(forest, tmp_path / "out.py")
    code = script.read_text()

    start = code.index('PAYLOAD = json.loads(r"""') + len('PAYLOAD = json.loads(r"""')
    end = code.index('""")', start)
    payload = json.loads(code[start:end])
    expected = sum(t.get_number_of_branches() for t in forest.trees)
    assert len(payload["branches"]) == expected
    assert payload["tree_count"] == len(forest.trees)
