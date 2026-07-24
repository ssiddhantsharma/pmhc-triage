"""pmhc-triage: HLA-coverage-adjusted addressable-population estimates for pMHC targets.

Public API is intentionally small and grows as modules land. Currently:

- :mod:`pmhc_triage.provenance` -- ``Provenance`` / ``Sourced`` (every number is sourced)
- :mod:`pmhc_triage.hla` -- diploid population coverage (the load-bearing core)
"""

from __future__ import annotations

from .hla import coverage_by_population, normalize_allele, parse_locus, population_coverage
from .provenance import Provenance, Sourced, today_iso

__version__ = "0.0.1"

__all__ = [
    "Provenance",
    "Sourced",
    "today_iso",
    "population_coverage",
    "coverage_by_population",
    "parse_locus",
    "normalize_allele",
    "__version__",
]
