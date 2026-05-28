"""Tests for the YAML config loader (Step 11)."""

import math
from pathlib import Path

import pytest

from mechatree.config import Config, GenomeConfig, LightConfig, TreeConfig, WindConfig, load_config


def test_treeconfig_defaults():
    cfg = TreeConfig()
    assert cfg.twig_length == 1.0
    assert cfg.twig_diameter == 0.1
    assert cfg.leaf_surface == 0.25
    assert cfg.cauchy == pytest.approx(4.0e-5)
    assert cfg.volume_ratio_leaf == 4.0
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
    # Step 19: default tau matches Eloy et al. (Nat Commun 2017).
    assert lc.leaf_transparency == 0.5


def test_lightconfig_validation():
    with pytest.raises(ValueError):
        LightConfig(n_elevations=0)
    with pytest.raises(ValueError):
        LightConfig(n_azimuths=-1)
    with pytest.raises(ValueError):
        LightConfig(size_leaf=0.0)
    with pytest.raises(ValueError, match="leaf_transparency"):
        LightConfig(leaf_transparency=-0.1)
    with pytest.raises(ValueError, match="leaf_transparency"):
        LightConfig(leaf_transparency=1.1)


def test_config_from_yaml_example_file():
    """The example config in examples/forest.yaml loads cleanly."""
    cfg = load_config(Path(__file__).parent.parent / "examples" / "forest.yaml")
    assert cfg.tree.cauchy == pytest.approx(4.0e-5)
    assert cfg.tree.volume_ratio_leaf == 4.0
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


def test_genomeconfig_neural_from_default_is_none():
    assert GenomeConfig().neural_from is None


def test_genomeconfig_neural_from_accepts_well_formed_dict():
    gc = GenomeConfig(neural_from={"path": "champions.json", "species_id": 1})
    assert gc.neural_from == {"path": "champions.json", "species_id": 1}


def test_genomeconfig_neural_from_validation():
    with pytest.raises(ValueError, match="must be a dict"):
        GenomeConfig(neural_from="not a dict")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="requires a 'path' key"):
        GenomeConfig(neural_from={"species_id": 0})
    with pytest.raises(ValueError, match="must be a string"):
        GenomeConfig(neural_from={"path": 42})
    with pytest.raises(ValueError, match="species_id.*must be an int"):
        GenomeConfig(neural_from={"path": "x.json", "species_id": "0"})


def test_config_from_dict_reads_neural_block():
    data = {"genome": {"neural_from": {"path": "../data/S3_champions.json", "species_id": 1}}}
    cfg = Config.from_dict(data, base_dir=Path("/tmp/somewhere"))
    assert cfg.genome.neural_from == {
        "path": "../data/S3_champions.json",
        "species_id": 1,
    }
    assert cfg.base_dir == Path("/tmp/somewhere")


def test_config_base_dir_not_in_equality():
    """base_dir is metadata about provenance, not part of the config's value."""
    a = Config.from_dict({}, base_dir=Path("/a"))
    b = Config.from_dict({}, base_dir=Path("/b"))
    assert a == b
    assert a.base_dir != b.base_dir


def test_config_from_yaml_stashes_base_dir(tmp_path):
    yml = tmp_path / "cfg.yaml"
    yml.write_text("genome:\n  safety: 2.5\n")
    cfg = Config.from_yaml(yml)
    assert cfg.base_dir == tmp_path


# ---------------------------------------------------------------------------
# Step 24: WindConfig fields for the wind ↔ pruning fixed-point loop.
# ---------------------------------------------------------------------------


def test_windconfig_step24_defaults():
    wc = WindConfig()
    assert wc.max_pruning_iterations == 8
    assert wc.wind_convergence_eps_rel == pytest.approx(0.01)


def test_windconfig_max_pruning_iterations_validation():
    # 1 is the minimum (recovers the old single-pass behaviour).
    WindConfig(max_pruning_iterations=1)
    with pytest.raises(ValueError, match="max_pruning_iterations"):
        WindConfig(max_pruning_iterations=0)
    with pytest.raises(ValueError, match="max_pruning_iterations"):
        WindConfig(max_pruning_iterations=-1)


def test_windconfig_eps_rel_validation():
    # 0 is allowed (disables the early-exit).
    WindConfig(wind_convergence_eps_rel=0.0)
    with pytest.raises(ValueError, match="wind_convergence_eps_rel"):
        WindConfig(wind_convergence_eps_rel=-0.01)


def test_windconfig_step24_yaml_roundtrip(tmp_path):
    yml = tmp_path / "cfg.yaml"
    yml.write_text(
        "wind:\n  model: default\n  max_pruning_iterations: 3\n  wind_convergence_eps_rel: 0.05\n"
    )
    cfg = Config.from_yaml(yml)
    assert cfg.wind.max_pruning_iterations == 3
    assert cfg.wind.wind_convergence_eps_rel == pytest.approx(0.05)


# Momentum-wind uniform-inflow override (U_in = K, independent of z).


def test_windconfig_momentum_u_uniform_default_is_none():
    wc = WindConfig(model="momentum")
    assert wc.momentum_U_uniform is None


def test_windconfig_momentum_u_uniform_accepts_positive():
    wc = WindConfig(model="momentum", momentum_U_uniform=1.6)
    assert wc.momentum_U_uniform == pytest.approx(1.6)


def test_windconfig_momentum_u_uniform_rejects_nonpositive():
    with pytest.raises(ValueError, match="momentum_U_uniform"):
        WindConfig(model="momentum", momentum_U_uniform=0.0)
    with pytest.raises(ValueError, match="momentum_U_uniform"):
        WindConfig(model="momentum", momentum_U_uniform=-2.0)


def test_windconfig_momentum_u_uniform_yaml_roundtrip(tmp_path):
    yml = tmp_path / "cfg.yaml"
    yml.write_text("wind:\n  model: momentum\n  grid_size: 1.0\n  momentum_U_uniform: 1.60\n")
    cfg = Config.from_yaml(yml)
    assert cfg.wind.model == "momentum"
    assert cfg.wind.momentum_U_uniform == pytest.approx(1.60)


def test_windconfig_n_sensing_angles_default_and_validation():
    assert WindConfig().n_sensing_angles == 4
    with pytest.raises(ValueError, match="n_sensing_angles"):
        WindConfig(n_sensing_angles=0)


def test_windconfig_n_sensing_angles_yaml_roundtrip(tmp_path):
    yml = tmp_path / "cfg.yaml"
    yml.write_text("wind:\n  model: momentum\n  n_sensing_angles: 8\n")
    cfg = Config.from_yaml(yml)
    assert cfg.wind.n_sensing_angles == 8
