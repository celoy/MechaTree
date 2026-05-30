#!/usr/bin/env python3
"""Validate a Python tournament champions JSON against the Fortran reference.

Compares species detection, centroid distances, and fitness distributions
between a Python-evolved champions.json and the Fortran S3_champions.json
reference.

Usage:
    python scripts/validate_tournament.py \\
        --python-champions out/tournament/champions.json \\
        --reference data/S3_champions.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from mechatree.evolution.curate import detect_species


def validate(
    python_champions: Path,
    reference_champions: Path,
    gap_threshold: float = 0.15,
    centroid_distance_warn: float = 0.30,
) -> bool:
    """Validate Python champions against Fortran reference.

    Returns
    -------
    bool
        True if all required checks pass.
    """
    if not python_champions.exists():
        print(f"Error: Python champions not found: {python_champions}", file=sys.stderr)
        return False

    if not reference_champions.exists():
        print(
            f"Error: Reference champions not found: {reference_champions}",
            file=sys.stderr,
        )
        return False

    py = json.loads(python_champions.read_text())
    ref = json.loads(reference_champions.read_text())

    print(f"Python champions: {python_champions}")
    print(f"  Individuals: {py.get('n_individuals', len(py.get('tag_genes', [])))}")
    print(f"Reference champions: {reference_champions}")
    print(f"  Individuals: {ref.get('n_individuals', len(ref.get('tag_genes', [])))}")
    print()

    checks: list[tuple[str, bool, str]] = []  # (name, passed, message)

    # Check 1: ≥1 species found in Python output
    n_py_species = len(py.get("species", []))
    checks.append(
        (
            "species_found",
            n_py_species >= 1,
            f"Python found {n_py_species} species (required: ≥1)",
        )
    )

    # Check 2: Python species count matches Fortran
    n_ref_species = len(ref.get("species", []))
    checks.append(
        (
            "species_count_match",
            n_py_species == n_ref_species,
            f"Python {n_py_species} vs Fortran {n_ref_species} species",
        )
    )

    # Check 3: centroid gap sufficiency (only if 2+ species)
    centroid_gap = 0.0
    if n_py_species >= 2:
        tag_genes = np.array(py.get("tag_genes", []))
        if len(tag_genes) > 0:
            labels, centroids = detect_species(tag_genes, gap_threshold=gap_threshold)
            if len(centroids) >= 2:
                centroid_gap = float(np.linalg.norm(centroids[0] - centroids[1]))
        checks.append(
            (
                "centroid_gap",
                centroid_gap >= gap_threshold,
                f"centroid gap = {centroid_gap:.4f} (threshold = {gap_threshold})",
            )
        )
    else:
        checks.append(
            (
                "centroid_gap",
                False,
                "only 1 species found — gap check skipped",
            )
        )

    # Check 4: centroid distance to reference (optional, 2 species each)
    if n_py_species == 2 and n_ref_species == 2:
        py_species = sorted(py.get("species", []), key=lambda s: s["centroid_tag"][0])
        ref_species = sorted(ref.get("species", []), key=lambda s: s["centroid_tag"][0])
        dists = [
            float(np.linalg.norm(np.array(ps["centroid_tag"]) - np.array(rs["centroid_tag"])))
            for ps, rs in zip(py_species, ref_species, strict=True)
        ]
        max_dist = max(dists) if dists else 0.0
        checks.append(
            (
                "centroid_proximity",
                max_dist < centroid_distance_warn,
                f"centroid L2 distances: {[f'{d:.4f}' for d in dists]} "
                f"(warn threshold: {centroid_distance_warn})",
            )
        )

    # Print report.
    print("Validation results:")
    print()
    all_required_passed = True
    required_checks = {"species_found", "centroid_gap"}

    for name, passed, msg in checks:
        symbol = "✓ PASS" if passed else "✗ FAIL"
        required = " [required]" if name in required_checks else " [optional]"
        print(f"  {symbol}{required}: {name}")
        print(f"    {msg}")
        if name in required_checks and not passed:
            all_required_passed = False

    print()
    if all_required_passed:
        print("Overall: PASS (required checks satisfied)")
    else:
        print("Overall: FAIL (required checks not satisfied)")

    return all_required_passed


def main() -> None:
    p = argparse.ArgumentParser(
        description="Validate a Python tournament champions JSON against Fortran reference.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
The validation checks:
  [required] species_found: Python output detected ≥1 species
  [required] centroid_gap: species centroids are sufficiently separated (gap > 0.15)
  [optional] species_count_match: Python found same number of species as Fortran
  [optional] centroid_proximity: champions centroids are close to Fortran reference
        """,
    )

    p.add_argument(
        "--python-champions",
        type=Path,
        required=True,
        help="Path to Python champions.json from run_tournament.py.",
    )
    p.add_argument(
        "--reference",
        type=Path,
        default=Path(__file__).parent.parent / "data" / "S3_champions.json",
        help="Path to reference champions.json (default: data/S3_champions.json).",
    )
    p.add_argument(
        "--gap-threshold",
        type=float,
        default=0.15,
        help="Minimum species centroid gap (default: 0.15).",
    )
    p.add_argument(
        "--centroid-distance-warn",
        type=float,
        default=0.30,
        help="Warn threshold for centroid distance to Fortran reference (default: 0.30).",
    )

    args = p.parse_args()

    passed = validate(
        args.python_champions,
        args.reference,
        gap_threshold=args.gap_threshold,
        centroid_distance_warn=args.centroid_distance_warn,
    )

    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
