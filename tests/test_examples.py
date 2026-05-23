"""Smoke tests for the example scripts under ``examples/``.

Each script is run as a subprocess with ``--iterations 5 --seed 0`` to keep
runtime tiny and deterministic.
"""

import os
import subprocess
import sys
from pathlib import Path

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
    result = _run("self_avoiding.py", "--out-dir", str(out_dir), cwd=tmp_path)
    pngs = list(out_dir.glob("fig*.png"))
    assert pngs, "self_avoiding should produce at least one snapshot"
    assert "Final coral size" in result.stdout


def test_sap_transport_runs(tmp_path):
    # sap_transport prints either "Tree died." or "Final tree size: N branches.";
    # either is a successful run.
    result = _run("sap_transport.py", cwd=tmp_path)
    assert "Tree died" in result.stdout or "Final tree size" in result.stdout
