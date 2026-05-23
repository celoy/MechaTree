def test_import_package():
    import mechatree

    assert mechatree.__version__


def test_pytree_importable():
    from mechatree import PyTree

    assert PyTree is not None


def test_pytree_instantiable():
    from mechatree import PyTree

    t = PyTree({"length": 1.0, "radius": 0.1})
    assert t.get_number_of_branches() == 1
