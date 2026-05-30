#!/usr/bin/env python3
"""Run an evolutionary tournament across multiple independent replicates.

Each replicate is a separate Forest with its own RNG seed. The independent
runs are parallelized using multiprocessing.Pool (spawn context on macOS).
After all replicates complete, surviving genomes are combined and curated
for champions.

Usage:
    python scripts/run_tournament.py \\
        --config examples/tournament_natcomm.yaml \\
        --n-gens 1000 \\
        --n-replicates 3 \\
        --archive-every 100 \\
        --out-dir out/tournament
"""

from __future__ import annotations

import argparse
import multiprocessing
import os
import sys
from pathlib import Path

import numpy as np

from mechatree.config import load_config
from mechatree.evolution import curate
from mechatree.evolution.genome import Genome
from mechatree.evolution.run import run_tournament


def _run_one_replicate(args: tuple) -> dict:
    """Worker function for multiprocessing.Pool.

    Must be a top-level function (not a closure) for pickling on macOS spawn.
    Returns a dict with genome dicts and fitness arrays.
    """
    (
        replicate_id,
        config_path,
        n_gens,
        seed,
        archive_every,
        archive_dir,
        resume_from,
    ) = args

    # Reload config in worker process.
    config = load_config(Path(config_path))

    # Run tournament for this replicate.
    result = run_tournament(
        config,
        n_gens,
        seed=seed,
        archive_every=archive_every,
        archive_dir=archive_dir,
        resume_from=resume_from,
    )

    # Extract surviving genomes and fitness for merging later.
    genomes = result.forest.genomes or []
    fitness = result.fitness.tolist() if len(genomes) > 0 else []

    return {
        "replicate_id": replicate_id,
        "n_genomes": len(genomes),
        "genomes_dicts": [g.to_dict() for g in genomes],
        "fitness": fitness,
    }


def main() -> None:
    p = argparse.ArgumentParser(
        description="Run an evolutionary tournament across multiple independent replicates.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Smoke test: 2 replicates, 10 gens each
  python scripts/run_tournament.py --n-gens 10 --n-replicates 2 --n-workers 2

  # Full Nat Commun tournament: 3 replicates, 100k gens each, 12 workers
  python scripts/run_tournament.py \\
    --config examples/tournament_natcomm.yaml \\
    --n-gens 100000 \\
    --n-replicates 3 \\
    --n-workers 12 \\
    --archive-every 1000 \\
    --out-dir out/tournament_natcomm
        """,
    )

    p.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent.parent / "examples" / "tournament_natcomm.yaml",
        help="Config YAML file (default: examples/tournament_natcomm.yaml).",
    )
    p.add_argument(
        "--n-gens",
        type=int,
        default=None,
        help="Number of generations (overrides YAML).",
    )
    p.add_argument(
        "--n-replicates",
        type=int,
        default=1,
        help="Number of independent replicates (default: 1).",
    )
    p.add_argument(
        "--archive-every",
        type=int,
        default=1000,
        help="Write archive snapshot every K generations (default: 1000).",
    )
    p.add_argument(
        "--resume-from",
        type=Path,
        default=None,
        help="Resume all replicates from an archive_XXXXXXXX.json snapshot.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("out/tournament"),
        help="Output directory for archives and champions.json (default: out/tournament).",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Base RNG seed (replicate i uses seed+i, default: 0).",
    )
    p.add_argument(
        "--n-workers",
        type=int,
        default=None,
        help="Process pool size (default: min(n_replicates, cpu_count)).",
    )
    p.add_argument(
        "--champions-out",
        type=Path,
        default=None,
        help="Where to write combined champions.json (default: OUT_DIR/champions.json).",
    )
    p.add_argument(
        "--no-progress",
        action="store_true",
        help="Suppress progress messages.",
    )

    args = p.parse_args()

    # Load and validate config.
    if not args.config.exists():
        print(f"Error: config file not found: {args.config}", file=sys.stderr)
        sys.exit(1)

    config = load_config(args.config)
    n_gens = args.n_gens if args.n_gens is not None else config.n_generations

    if not args.no_progress:
        print(f"Tournament config: {args.config}")
        print(f"  Generations: {n_gens}")
        print(f"  Replicates: {args.n_replicates}")
        print(f"  Base seed: {args.seed}")
        print(f"  Output dir: {args.out_dir}")
        if args.resume_from:
            print(f"  Resuming from: {args.resume_from}")

    # Create replicate dirs.
    args.out_dir.mkdir(parents=True, exist_ok=True)
    replicate_dirs = [args.out_dir / f"replicate_{i:02d}" for i in range(args.n_replicates)]
    for d in replicate_dirs:
        d.mkdir(parents=True, exist_ok=True)

    # Build worker arguments: (replicate_id, config_path, n_gens, seed, ...).
    worker_args = [
        (
            i,
            str(args.config.absolute()),
            n_gens,
            args.seed + i,
            args.archive_every,
            str(replicate_dirs[i].absolute()),
            str(args.resume_from.absolute()) if args.resume_from else None,
        )
        for i in range(args.n_replicates)
    ]

    # Run replicates in parallel.
    n_workers = args.n_workers or min(args.n_replicates, os.cpu_count() or 1)
    if not args.no_progress:
        print(f"  Workers: {n_workers}")
        print()

    ctx = multiprocessing.get_context("spawn")
    try:
        with ctx.Pool(processes=n_workers) as pool:
            if args.no_progress:
                results = pool.map(_run_one_replicate, worker_args)
            else:
                results = pool.imap_unordered(_run_one_replicate, worker_args)
                results = list(results)  # Consume iterator to show progress.
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        sys.exit(1)

    if not args.no_progress:
        print(f"\nAll {args.n_replicates} replicate(s) completed.")

    # Combine genomes and fitness across all replicates.
    all_genomes: list[Genome] = []
    all_fitness: list[float] = []
    total_survivors = 0

    for result in results:
        replicate_id = result["replicate_id"]
        n_genomes = result["n_genomes"]
        total_survivors += n_genomes

        if n_genomes > 0:
            all_genomes.extend(Genome.from_dict(d) for d in result["genomes_dicts"])
            all_fitness.extend(result["fitness"])

        if not args.no_progress:
            print(f"  Replicate {replicate_id}: {n_genomes} survivors")

    if not all_genomes:
        print("Warning: no surviving genomes across all replicates.", file=sys.stderr)
        sys.exit(1)

    if not args.no_progress:
        print(f"\nTotal survivors: {total_survivors}")
        print(f"Combined population: {len(all_genomes)} genomes")

    # Curate champions from combined population.
    combined_fitness = np.array(all_fitness, dtype=np.float64)
    payload = curate._payload(
        all_genomes,
        combined_fitness,
        source=str(args.out_dir),
    )

    # Write champions.
    champions_path = args.champions_out or args.out_dir / "champions.json"
    champions_path = curate.write(payload, champions_path)

    if not args.no_progress:
        n_species = len(payload.get("species", []))
        print(f"\nWrote champions JSON ({len(all_genomes)} individuals, {n_species} species)")
        print(f"  → {champions_path}")


if __name__ == "__main__":
    main()
