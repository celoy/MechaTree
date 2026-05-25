"""Per-tree evolutionary dynamics on a Forest (Step 21).

The Python port replaces the Fortran's external NSGA-style tournament with
the simpler model the simulator was actually approximating: an *in silico*
Darwinian island where every tree carries a heritable :class:`Genome`,
seeds inherit a mutated copy of the parent's genome, and species disappear
when no tree carrying that lineage survives long enough to reproduce.

The forest **is** the population — there is no central selection step.

Entry points::

    from mechatree.evolution import Genome, run_tournament, curate

    result = run_tournament(cfg, n_generations=50, seed=0)
    champions = curate.from_forest(result.forest, result.fitness)
"""

from __future__ import annotations

from mechatree.evolution import archive, curate
from mechatree.evolution.genome import Genome
from mechatree.evolution.run import ForestEvolutionResult, run_tournament

__all__ = [
    "ForestEvolutionResult",
    "Genome",
    "archive",
    "curate",
    "run_tournament",
]
