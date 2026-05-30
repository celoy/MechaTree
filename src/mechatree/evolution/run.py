"""High-level entry point: ``run_tournament(cfg, n_generations, seed)``.

Drives the Step-12 :class:`mechatree.forest.Forest` in evolution mode
(``genomes=[...]``), tallies per-lineage lifetime seeds as fitness, and
optionally writes per-step archive snapshots + a final champion JSON
(consumed by :func:`mechatree.genome.load_champion`).

The fitness accumulator is intentionally tiny: each tree's ``n_branches``
at the *moment of death* tallies into ``self._fitness[lineage_id]``, plus
the final-snapshot ``n_branches`` for the lineages still alive. This is
the natural Darwinian fitness in an island model — "how big a tree this
genome managed to grow in the time it had".
"""

from __future__ import annotations

import dataclasses
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from mechatree.config import Config
from mechatree.evolution import archive, curate
from mechatree.evolution.genome import Genome
from mechatree.forest import Forest, ForestStats

OnStep = Callable[[int, "ForestEvolutionResult", ForestStats], None]


@dataclass
class ForestEvolutionResult:
    """Return value of :func:`run_tournament`.

    * ``forest``     — the live Forest at the end of the run.
    * ``fitness``    — per-tree-alive scalar (cumulative ``n_branches``)
                       parallel to ``forest.genomes``.
    * ``history``    — list of ``ForestStats`` per generation.
    * ``champions_path`` — written ``champions.json`` (None if
                       ``write_champions=False``).
    * ``archive_dir`` — directory holding per-step snapshots (None if
                       ``archive_every`` was None).
    """

    forest: Forest
    fitness: np.ndarray
    history: list[ForestStats]
    champions_path: Path | None = None
    archive_dir: Path | None = None


def run_tournament(
    config: Config,
    n_generations: int,
    *,
    seed: int = 0,
    initial_genomes: list[Genome] | None = None,
    mutation_sigma: float = 0.005,
    mutation_p_locus: float = 0.05,
    archive_every: int | None = None,
    archive_dir: Path | str | None = None,
    champions_path: Path | str | None = None,
    resume_from: Path | str | None = None,
    on_step: OnStep | None = None,
) -> ForestEvolutionResult:
    """Run a Darwinian-island tournament on a Forest.

    Parameters
    ----------
    config
        The standard :class:`mechatree.config.Config`. ``config.forest``
        controls island size + carrying capacity. The shared
        ``config.genome.safety``/``allocation`` are ignored (each tree
        carries its own).
    n_generations
        Number of forest steps.
    seed
        Master RNG seed for the Forest + initial genome sampling.
    initial_genomes
        Founder genomes (length ``config.forest.n_trees_init``). When
        ``None``, samples uniformly random genomes via
        :meth:`Genome.random`, each assigned a unique ``lineage_id`` in
        ``[0, n_trees_init)``.
    mutation_sigma, mutation_p_locus
        Per-locus Gaussian mutation parameters; Fortran defaults from
        ``legacy/fortran/Evolution.ini``.
    archive_every
        If set, write ``archive_{gen:08d}.json`` every K generations
        under ``archive_dir``.
    archive_dir
        Where to write archive snapshots (required if ``archive_every``
        is set).
    champions_path
        If set, after the loop curate the final population and write a
        champion JSON to this path.
    resume_from
        Path to an ``archive_XXXXXXXX.json`` snapshot to resume from.
        When set, the run starts at ``generation + 1`` with the surviving
        population from that snapshot. Overrides ``initial_genomes``.
    on_step
        Optional ``(generation, result, stats)`` callback.

    Returns
    -------
    ForestEvolutionResult
        On a resumed run, ``history`` only covers the generations from
        ``start_gen`` to ``n_generations - 1`` (i.e., the resumed portion).
        To reconstruct a full run history, concatenate the history from
        prior archive snapshots with ``result.history``.
    """
    rng = np.random.default_rng(seed)

    start_gen = 0
    n_init = config.forest.n_trees_init

    if resume_from is not None:
        # Load the snapshot and resume from the next generation.
        resume_path = Path(resume_from)
        start_gen, founders, _metas = archive.load_snapshot(resume_path)
        start_gen += 1  # snapshot contains generation N; resume from N+1
        # Patch config to accept the (typically smaller) surviving population.
        config = dataclasses.replace(
            config,
            forest=dataclasses.replace(config.forest, n_trees_init=len(founders)),
        )
    elif initial_genomes is None:
        founders = [Genome.random(rng, lineage_id=i) for i in range(n_init)]
    else:
        if len(initial_genomes) != n_init:
            raise ValueError(
                f"initial_genomes length {len(initial_genomes)} != n_trees_init {n_init}"
            )
        founders = list(initial_genomes)

    forest = Forest(
        config=config,
        seed=seed,
        genomes=founders,
        mutation_sigma=mutation_sigma,
        mutation_p_locus=mutation_p_locus,
    )

    history: list[ForestStats] = []
    archive_path = Path(archive_dir) if archive_dir is not None else None
    if archive_every is not None and archive_path is None:
        raise ValueError("archive_every is set but archive_dir is None")

    # Place-holder so we can return a populated result from the callback.
    result = ForestEvolutionResult(
        forest=forest,
        fitness=np.zeros(len(forest.genomes), dtype=np.float64),
        history=history,
        archive_dir=archive_path,
    )

    for gen in range(start_gen, n_generations):
        stats = forest.step(gen)
        history.append(stats)
        if archive_every is not None and gen % archive_every == 0:
            archive.write_snapshot(forest, gen, archive_path)  # type: ignore[arg-type]
        if on_step is not None:
            result.fitness = _alive_fitness(forest)
            on_step(gen, result, stats)

    result.fitness = _alive_fitness(forest)

    if champions_path is not None and forest.genomes:
        payload = curate.from_forest(forest, result.fitness)
        result.champions_path = curate.write(payload, Path(champions_path))

    return result


def _alive_fitness(forest: Forest) -> np.ndarray:
    """Per-alive-tree fitness scalar. Cumulative ``n_branches`` proxies
    lifetime reproductive success (more branches → more leaves → more
    seeds dropped) without requiring a parallel ledger of dead lineages."""
    n = len(forest.trees)
    if n == 0:
        return np.zeros(0)
    return np.array(
        [tree.get_number_of_branches() for tree in forest.trees],
        dtype=np.float64,
    )


__all__ = ["ForestEvolutionResult", "run_tournament"]
