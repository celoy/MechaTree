"""Smoke tests for the example scripts under ``examples/``.

Each script is run as a subprocess with ``--iterations 5 --seed 0`` to keep
runtime tiny and deterministic.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = REPO_ROOT / "examples"


def _run(script: str, *extra: str, cwd: Path) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    return subprocess.run(
        [sys.executable, str(EXAMPLES / script), "--iterations", "5", "--seed", "0", *extra],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )


def test_random_growth_runs(tmp_path):
    result = _run("random_growth.py", "--no-show", cwd=tmp_path)
    assert "Final tree size" in result.stdout


def test_self_avoiding_runs(tmp_path):
    out_dir = tmp_path / "snaps"
    try:
        result = _run("self_avoiding.py", "--out-dir", str(out_dir), cwd=tmp_path)
    except subprocess.CalledProcessError as e:
        if "chrome" in e.stderr.lower() or "kaleido" in e.stderr.lower():
            pytest.skip("Chrome not available for kaleido")
        raise
    pngs = list(out_dir.glob("fig*.png"))
    assert pngs, "self_avoiding should produce at least one snapshot"
    assert "Final coral size" in result.stdout


def test_sap_transport_runs(tmp_path):
    # sap_transport prints either "Tree died." or "Final tree size: N branches.";
    # either is a successful run.
    result = _run("sap_transport.py", cwd=tmp_path)
    assert "Tree died" in result.stdout or "Final tree size" in result.stdout


# --- Step 13: simulation tutorials ----------------------------------------


def test_grow_one_tree_runs(tmp_path):
    result = _run("grow_one_tree.py", "--no-show", cwd=tmp_path)
    assert "Final tree:" in result.stdout


def test_grow_a_forest_runs(tmp_path):
    result = _run(
        "grow_a_forest.py",
        "--no-show",
        "--n-trees-init",
        "5",
        "--n-trees-max",
        "30",
        "--size",
        "15",
        cwd=tmp_path,
    )
    assert "Final:" in result.stdout


def test_custom_simulation_runs(tmp_path):
    result = _run("custom_simulation.py", "--no-show", cwd=tmp_path)
    # All three labels should appear in the summary table.
    assert "default" in result.stdout
    assert "steady-west" in result.stdout
    assert "calm-then-storm" in result.stdout


# --- MATLAB plot ports + Blender exporter ---------------------------------


def test_plot_self_thinning_runs(tmp_path):
    result = _run(
        "plot_self_thinning.py",
        "--no-show",
        "--n-trees-init",
        "5",
        "--n-trees-max",
        "30",
        "--size",
        "15",
        cwd=tmp_path,
    )
    assert "M (biomass)" in result.stdout


def test_plot_allocation_runs(tmp_path):
    result = _run("plot_allocation.py", "--no-show", cwd=tmp_path)
    assert "reserve" in result.stdout


def test_plot_strahler_runs(tmp_path):
    result = _run("plot_strahler.py", "--no-show", cwd=tmp_path)
    assert "Strahler max order" in result.stdout


def test_render_blender_writes_script(tmp_path):
    out_base = tmp_path / "my_render"
    result = _run("render_blender.py", "--output", str(out_base), cwd=tmp_path)
    assert "Wrote:" in result.stdout
    assert out_base.with_suffix(".py").exists()
