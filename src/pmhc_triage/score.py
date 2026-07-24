"""Combine the factors into an effective addressable-N -- the exit gate.

This is the single place where the pipeline decides whether a number can be emitted
at all. Three rules, all "fail loud, never fake":

1. **Join guard.** ``effective_N(population) = incidence x antigen_fraction x
   hla_coverage`` is computed ONLY per matching population label. Incidence and HLA
   coverage must refer to the *same* population, or they are simply not combined --
   there is no code path that multiplies incidence["India"] by coverage["Europe"].
2. **Propagate missing.** If any required factor for a population is missing
   (``None``), that population's ``effective_n`` is ``None`` *with the reasons* --
   never silently ``0``.
3. **One provenance log.** Every factor's provenance is collected into a single
   auditable log for the whole result.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

from .provenance import Sourced


@dataclass
class PopulationScore:
    """Effective-N for one population, with the factors and reasons it (didn't) compute."""

    population: str
    effective_n: float | None
    incidence: Sourced[float] | None
    antigen_fraction: Sourced[float]
    hla_coverage: Sourced[float] | None
    reasons: list[str] = field(default_factory=list)

    @property
    def computable(self) -> bool:
        return self.effective_n is not None

    def to_dict(self) -> dict[str, Any]:
        def val(s):
            return None if s is None else s.value

        af_extra = self.antigen_fraction.extra if self.antigen_fraction else {}
        ci = None
        if "ci95_low" in af_extra and "ci95_high" in af_extra:
            ci = [af_extra["ci95_low"], af_extra["ci95_high"]]
        return {
            "population": self.population,
            "effective_n": self.effective_n,
            "computable": self.computable,
            "incidence": val(self.incidence),
            "antigen_fraction": val(self.antigen_fraction),
            "antigen_n": af_extra.get("denominator"),
            "antigen_fraction_ci95": ci,
            "hla_coverage": val(self.hla_coverage),
            "reasons": list(self.reasons),
        }


@dataclass
class TargetScore:
    gene: str
    variant: str
    disease: str
    per_population: dict[str, PopulationScore]
    provenance_log: list[dict]
    warnings: list[str] = field(default_factory=list)

    def rows(self) -> list[dict]:
        """Flat rows (one per population) for CSV / a report -- factors never hidden."""
        return [
            {"gene": self.gene, "variant": self.variant, "disease": self.disease, **ps.to_dict()}
            for ps in self.per_population.values()
        ]

    def to_dict(self) -> dict[str, Any]:
        return {
            "target": {"gene": self.gene, "variant": self.variant, "disease": self.disease},
            "per_population": {p: ps.to_dict() for p, ps in self.per_population.items()},
            "warnings": list(self.warnings),
            "provenance_log": self.provenance_log,
        }


def score_target(
    *,
    gene: str,
    variant: str,
    disease: str,
    antigen_fraction: Sourced[float],
    incidence_by_population: Mapping[str, Sourced[float]],
    coverage_by_population: Mapping[str, Sourced[float]],
    tractability_context: Sourced[list] | None = None,
) -> TargetScore:
    """Combine factors into per-population effective-N, enforcing the join guard.

    ``antigen_fraction`` is treated as population-agnostic (its source cohort is in
    its provenance). ``incidence_by_population`` and ``coverage_by_population`` are
    keyed by population label and combined only where the *same* label appears in
    both with a present value.
    """
    warnings: list[str] = []
    inc_pops, cov_pops = set(incidence_by_population), set(coverage_by_population)
    only_inc, only_cov = inc_pops - cov_pops, cov_pops - inc_pops
    if only_inc:
        warnings.append(
            f"populations with incidence but no HLA coverage (not scored -- label mismatch?): {sorted(only_inc)}"
        )
    if only_cov:
        warnings.append(
            f"populations with HLA coverage but no incidence (not scored -- label mismatch?): {sorted(only_cov)}"
        )

    provenance_log: list[dict] = [{"factor": "antigen_fraction", **antigen_fraction.to_dict()}]

    per_population: dict[str, PopulationScore] = {}
    for pop in sorted(inc_pops | cov_pops):
        inc = incidence_by_population.get(pop)
        cov = coverage_by_population.get(pop)
        reasons: list[str] = []

        if inc is None:
            reasons.append(f"no incidence for population {pop!r}")
        elif inc.is_missing:
            reasons.append(f"incidence missing for {pop!r}: {'; '.join(inc.warnings) or 'no value'}")
        if cov is None:
            reasons.append(f"no HLA coverage for population {pop!r}")
        elif cov.is_missing:
            reasons.append(f"HLA coverage missing for {pop!r}: {'; '.join(cov.warnings) or 'no value'}")
        if antigen_fraction.is_missing:
            reasons.append(f"antigen fraction missing: {'; '.join(antigen_fraction.warnings) or 'no value'}")

        if inc is not None:
            provenance_log.append({"factor": "incidence", "population": pop, **inc.to_dict()})
        if cov is not None:
            provenance_log.append({"factor": "hla_coverage", "population": pop, **cov.to_dict()})

        if reasons:
            effective_n = None
        else:
            effective_n = round(inc.value * antigen_fraction.value * cov.value, 2)

        per_population[pop] = PopulationScore(
            population=pop,
            effective_n=effective_n,
            incidence=inc,
            antigen_fraction=antigen_fraction,
            hla_coverage=cov,
            reasons=reasons,
        )

    if tractability_context is not None:
        provenance_log.append({"factor": "tractability_context", **tractability_context.to_dict()})

    return TargetScore(gene, variant, disease, per_population, provenance_log, warnings)
