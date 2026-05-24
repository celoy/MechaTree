"""Tests for the YAML config loader (Step 11)."""

import math
from pathlib import Path

import pytest

from mechatree.config import Config, GenomeConfig, LightConfig, TreeConfig, load_config


def test_treeconfig_defaults():
    cfg = TreeConfig()
    assert cfg.twig_length == 1.0
    assert cfg.twig_diameter == 0.1
    assert cfg.leaf_surface == 0.25
    assert cfg.cauchy == pytest.approx(2.0e-5)
    assert cfg.volume_ratio_leaf == 8.0
    assert cfg.maintenance_h == 0.02
    assert cfg.theta1 == pytest.approx(math.pi / 4)
    assert cfg.theta2 == pytest.approx(-math.pi / 4)


def test_treeconfig_derived_volumes():
    cfg = TreeConfig(twig_length=2.0, twig_diameter=0.1, volume_ratio_leaf=4.0)
    expected_vol_twig = 0.25 * math.pi * 2.0 * 0.01
    assert cfg.volume_twig == pytest.approx(expected_vol_twig)
    assert cfg.volume_per_leaf == pytest.approx(4.0 * expected_vol_twig)


def test_treeconfig_rejects_nonpositive_scalars():
    with pytest.raises(ValueError):
        TreeConfig(twig_length=0.0)
    with pytest.raises(ValueError):
        TreeConfig(twig_diameter=-0.1)
    with pytest.raises(ValueError):
        TreeConfig(cauchy=0.0)
    with pytest.raises(ValueError):
        TreeConfig(max_branches=0)


def test_treeconfig_accepts_zero_maintenance():
    """maintenance_h can be zero (just not negative)."""
    cfg = TreeConfig(maintenance_h=0.0)
    assert cfg.maintenance_h == 0.0
    with pytest.raises(ValueError):
        TreeConfig(maintenance_h=-0.01)


def test_lightconfig_defaults():
    lc = LightConfig()
    assert lc.size_leaf == 1.0
    assert lc.n_elevations == 4
    assert lc.n_azimuths == 8


def test_lightconfig_validation():
    with pytest.raises(ValueError):
        LightConfig(n_elevations=0)
    with pytest.raises(ValueError):
        LightConfig(n_azimuths=-1)
    with pytest.raises(ValueError):
        LightConfig(size_leaf=0.0)


def test_config_from_yaml_example_file():
    """The example config in examples/forest.yaml loads cleanly."""
    cfg = load_config(Path(__file__).parent.parent / "examples" / "forest.yaml")
    assert cfg.tree.cauchy == pytest.approx(2.0e-5)
    assert cfg.tree.volume_ratio_leaf == 8.0
    assert cfg.light.n_elevations == 4
    assert cfg.n_generations > 0


def test_config_from_dict_ignores_unknown_keys():
    """The YAML may contain forest/evolution blocks for later steps — those
    must not crash the loader."""
    data = {
        "tree": {"twig_length": 0.5, "unknown_key": 42},
        "light": {"n_elevations": 2},
        "forest": {"n_generations": 50, "size": 100.0},
        "evolution": {"p_mutation": 0.1},
    }
    cfg = Config.from_dict(data)
    assert cfg.tree.twig_length == 0.5
    assert cfg.light.n_elevations == 2
    # forest.n_generations is honoured as a fallback.
    assert cfg.n_generations == 50


def test_config_top_level_n_generations_overrides_forest():
    data = {"n_generations": 200, "forest": {"n_generations": 50}}
    cfg = Config.from_dict(data)
    assert cfg.n_generations == 200


def test_config_empty_yaml_uses_all_defaults():
    """An empty YAML (or a YAML where every section is omitted) yields
    the dataclass defaults."""
    cfg = Config.from_dict({})
    assert cfg.tree == TreeConfig()
    assert cfg.light == LightConfig()
    assert cfg.genome == GenomeConfig()


def test_genomeconfig_defaults():
    gc = GenomeConfig()
    assert gc.safety == 3.0
    assert gc.p_seeds == 0.1
    assert gc.p_leaves == 0.5
    assert gc.phototropism == 0.5


def test_genomeconfig_validation():
    with pytest.raises(ValueError):
        GenomeConfig(safety=0.0)
    with pytest.raises(ValueError):
        GenomeConfig(p_seeds=-0.01)
    with pytest.raises(ValueError):
        GenomeConfig(p_leaves=-0.01)
    with pytest.raises(ValueError):
        GenomeConfig(p_seeds=0.6, p_leaves=0.6)
    with pytest.raises(ValueError):
        GenomeConfig(phototropism=1.5)


def test_config_from_dict_reads_genome_block():
    data = {"genome": {"safety": 4.2, "p_seeds": 0.2, "phototropism": 0.0}}
    cfg = Config.from_dict(data)
    assert cfg.genome.safety == pytest.approx(4.2)
    assert cfg.genome.p_seeds == pytest.approx(0.2)
    # Unset fields fall back to defaults.
    assert cfg.genome.p_leaves == pytest.approx(0.5)
    assert cfg.genome.phototropism == pytest.approx(0.0)


def test_config_yaml_example_carries_genome_block():
    cfg = load_config(Path(__file__).parent.parent / "examples" / "forest.yaml")
    assert cfg.genome.safety == pytest.approx(3.0)
    assert cfg.genome.p_seeds == pytest.approx(0.1)
    assert cfg.genome.p_leaves == pytest.approx(0.5)
    assert cfg.genome.phototropism == pytest.approx(0.5)
