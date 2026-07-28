"""Threshold-sensitivity report: how coverage / effective-N move with the cutoff.

The presentation threshold (which alleles count as "presenting") is the single most
arbitrary knob in the pipeline -- a looser cutoff admits more alleles and mechanically
inflates HLA coverage and therefore effective-N. The tool's contract is to *surface*
levers, not bury them in a single default. This module turns a threshold sweep into a
small table so the reader sees the whole curve and picks a defensible point themselves.

Given a scored (peptide, allele) table and a set of thresholds, it reports, per
threshold: the presenting-allele set and the resulting per-population HLA coverage
(via the same Bui-2006 diploid method as everywhere else). Optionally multiplies
through a fixed ``incidence x antigen_fraction`` to show effective-N sensitivity.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Sequence

from .hla import population_coverage
from .presentation import sweep_presenting_alleles


@dataclass
class ThresholdRow:
    """One threshold's presenting alleles + resulting coverage (and optional effective-N)."""

    threshold: float
    alleles: list[str]
    coverage_by_population: dict[str, float | None]
    effective_n_by_population: dict[str, float | None] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "threshold": self.threshold,
            "n_alleles": len(self.alleles),
            "alleles": list(self.alleles),
            "coverage_by_population": dict(self.coverage_by_population),
            "effective_n_by_population": dict(self.effective_n_by_population),
        }


def threshold_sensitivity(
    df,
    thresholds: Sequence[float],
    freqs_by_population: Mapping[str, Mapping[str, float]],
    *,
    allele_col: str = "sample_name",
    score_col: str = "presentation_percentile",
    lower_is_better: bool = True,
    incidence_by_population: Mapping[str, float] | None = None,
    antigen_fraction: float | None = None,
) -> list[ThresholdRow]:
    """Coverage (and optionally effective-N) across presentation thresholds.

    ``df`` is a scored (peptide, allele) table (e.g. MHCflurry output). ``thresholds``
    is the set to sweep (e.g. ``[0.5, 2.0, 5.0]``). ``freqs_by_population`` maps
    population -> {allele: freq}. If both ``incidence_by_population`` and
    ``antigen_fraction`` are given, the point effective-N is included per threshold so
    the reader sees how the final number swings with the cutoff.

    Coverage is ``None`` for a population where no presenting allele has a known
    frequency (surfaced, never 0), so effective-N is ``None`` there too.
    """
    alleles_by_threshold = sweep_presenting_alleles(
        df,
        thresholds,
        allele_col=allele_col,
        score_col=score_col,
        lower_is_better=lower_is_better,
    )
    rows: list[ThresholdRow] = []
    for t in sorted(alleles_by_threshold):
        alleles = alleles_by_threshold[t]
        cov_by_pop: dict[str, float | None] = {}
        eff_by_pop: dict[str, float | None] = {}
        for pop, freqs in freqs_by_population.items():
            cov = population_coverage(alleles, freqs, pop)
            cov_by_pop[pop] = cov.value
            if (
                incidence_by_population is not None
                and antigen_fraction is not None
                and cov.value is not None
                and incidence_by_population.get(pop) is not None
            ):
                eff_by_pop[pop] = round(
                    incidence_by_population[pop] * antigen_fraction * cov.value, 2
                )
            else:
                eff_by_pop[pop] = None
        rows.append(ThresholdRow(float(t), alleles, cov_by_pop, eff_by_pop))
    return rows
