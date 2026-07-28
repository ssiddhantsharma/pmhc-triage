"""pmhc-triage: HLA-coverage-adjusted addressable-population estimates for pMHC targets.

Public API is intentionally small and grows as modules land. Currently:

- :mod:`pmhc_triage.provenance` -- ``Provenance`` / ``Sourced`` (every number is sourced)
- :mod:`pmhc_triage.hla` -- diploid population coverage (the load-bearing core)
"""

from __future__ import annotations

from .afnd import FrequencyTable, load_afnd_frequencies
from .antigen import (
    expression_positive_fraction,
    resolve_entrez,
    variant_frequency,
    variant_frequency_multi,
)
from .burden import bundled_diseases, load_bundled, load_burden_table, manual_incidence
from .caching import cached_client
from .fetch import (
    COMMON_CLASS_I,
    DEFAULT_POPULATIONS,
    FetchResult,
    ReferencePopulation,
    fetch_frequencies,
    parse_afnd_table,
)
from .hla import (
    coverage_by_population,
    hla_class,
    normalize_allele,
    parse_locus,
    population_coverage,
)
from .identity import (
    align_populations,
    canonical_disease,
    canonical_population,
    disease_matches,
    suggest_match,
)
from .montecarlo import (
    MCResult,
    effective_n_interval,
    sample_antigen_fraction,
    sample_coverage,
)
from .opentargets import (
    associated_targets,
    resolve_disease,
    resolve_target,
    tractability,
)
from .peptides import mutant_peptides, parse_substitution
from .pipeline import TargetSpec, preflight, run_target
from .presentation import (
    manual_presenting_alleles,
    predict_presenting_alleles,
    select_presenting,
    sweep_presenting_alleles,
)
from .provenance import Provenance, Sourced, today_iso
from .report import write_provenance_json, write_results_csv
from .score import PopulationScore, TargetScore, score_target
from .sensitivity import ThresholdRow, threshold_sensitivity
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
    "hla_class",
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
    "variant_frequency_multi",
    "expression_positive_fraction",
    "resolve_entrez",
    "manual_incidence",
    "load_burden_table",
    "load_bundled",
    "bundled_diseases",
    "score_target",
    "TargetScore",
    "PopulationScore",
    "TargetSpec",
    "run_target",
    "preflight",
    "write_results_csv",
    "write_provenance_json",
    "manual_presenting_alleles",
    "predict_presenting_alleles",
    "select_presenting",
    "sweep_presenting_alleles",
    "MCResult",
    "effective_n_interval",
    "sample_antigen_fraction",
    "sample_coverage",
    "threshold_sensitivity",
    "ThresholdRow",
    "fetch_frequencies",
    "parse_afnd_table",
    "FetchResult",
    "ReferencePopulation",
    "DEFAULT_POPULATIONS",
    "COMMON_CLASS_I",
    "canonical_population",
    "canonical_disease",
    "disease_matches",
    "suggest_match",
    "align_populations",
    "cached_client",
    "__version__",
]
