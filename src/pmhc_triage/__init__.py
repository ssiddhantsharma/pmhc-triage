"""pmhc-triage: HLA-coverage-adjusted addressable-population estimates for pMHC targets.

Public API is intentionally small and grows as modules land. Currently:

- :mod:`pmhc_triage.provenance` -- ``Provenance`` / ``Sourced`` (every number is sourced)
- :mod:`pmhc_triage.hla` -- diploid population coverage (the load-bearing core)
"""

from __future__ import annotations

from .afnd import FrequencyTable, load_afnd_frequencies
from .antigen import resolve_entrez, variant_frequency
from .burden import load_burden_table, manual_incidence
from .score import PopulationScore, TargetScore, score_target
from .hla import coverage_by_population, normalize_allele, parse_locus, population_coverage
from .opentargets import (
    associated_targets,
    resolve_disease,
    resolve_target,
    tractability,
)
from .peptides import mutant_peptides, parse_substitution
from .pipeline import TargetSpec, preflight, run_target
from .provenance import Provenance, Sourced, today_iso
from .report import write_provenance_json, write_results_csv
from .sequences import fetch_uniprot_sequence

__version__ = "0.0.1"

__all__ = [
    "Provenance",
    "Sourced",
    "today_iso",
    "population_coverage",
    "coverage_by_population",
    "parse_locus",
    "normalize_allele",
    "fetch_uniprot_sequence",
    "mutant_peptides",
    "parse_substitution",
    "load_afnd_frequencies",
    "FrequencyTable",
    "resolve_target",
    "resolve_disease",
    "tractability",
    "associated_targets",
    "variant_frequency",
    "resolve_entrez",
    "manual_incidence",
    "load_burden_table",
    "score_target",
    "TargetScore",
    "PopulationScore",
    "TargetSpec",
    "run_target",
    "preflight",
    "write_results_csv",
    "write_provenance_json",
    "__version__",
]
