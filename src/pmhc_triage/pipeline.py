"""Orchestrate the factors into a scored target, plus a fail-fast preflight.

``run_target`` wires the modules together: antigen (cBioPortal) + tractability
context (Open Targets) + HLA coverage (your AFND-format frequencies) + burden
(your incidence table/values) -> the score gate. It never fabricates: any factor
that can't be fetched arrives at the score step as a surfaced ``Sourced(None)``.

``preflight`` is the clean early exit -- it checks the joins *before* running so a
misspelled study, a population absent from the frequency file, or a variant whose
wild-type residue doesn't match UniProt is caught up front rather than producing a
confident-looking but meaningless number.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import httpx

from .afnd import load_afnd_frequencies
from .antigen import check_study, resolve_entrez, variant_frequency
from .burden import load_burden_table, manual_incidence
from .hla import coverage_by_population, normalize_allele
from .opentargets import resolve_target, tractability
from .peptides import parse_substitution
from .presentation import manual_presenting_alleles
from .score import TargetScore, score_target
from .sequences import fetch_uniprot_sequence


@dataclass
class TargetSpec:
    gene: str
    variant: str
    disease: str
    study: str
    alleles: list[str]
    populations: list[str]
    freqs_path: str | None = None
    burden_path: str | None = None
    burden_manual: dict[str, float] | None = None
    burden_source: str | None = None
    uniprot: str | None = None

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "TargetSpec":
        return cls(
            gene=d["gene"],
            variant=d["variant"],
            disease=d["disease"],
            study=d["study"],
            alleles=list(d["alleles"]),
            populations=list(d["populations"]),
            freqs_path=d.get("freqs"),
            burden_path=d.get("burden") if isinstance(d.get("burden"), str) else None,
            burden_manual=d.get("burden") if isinstance(d.get("burden"), dict) else None,
            burden_source=d.get("burden_source"),
            uniprot=d.get("uniprot"),
        )


def run_target(spec: TargetSpec, *, client: httpx.Client | None = None) -> TargetScore:
    """Fetch every factor and combine via the score gate. Missing factors surface, never fake."""
    antigen = variant_frequency(spec.study, spec.variant, gene=spec.gene, client=client)

    tract = None
    tid = resolve_target(spec.gene, client=client)
    if not tid.is_missing:
        tract = tractability(tid.value, client=client)

    alleles = manual_presenting_alleles(spec.alleles)
    allele_list = alleles.value or []

    if spec.freqs_path:
        ft = load_afnd_frequencies(spec.freqs_path)
        freqs_by_pop = {p: ft.get(p) for p in spec.populations}
    else:
        freqs_by_pop = {p: {} for p in spec.populations}
    coverage = coverage_by_population(allele_list, freqs_by_pop)

    incidence: dict = {}
    if spec.burden_path:
        loaded = load_burden_table(spec.burden_path, spec.disease)
        incidence = {p: loaded[p] for p in spec.populations if p in loaded}
    elif spec.burden_manual:
        for pop, cases in spec.burden_manual.items():
            incidence[pop] = manual_incidence(
                spec.disease, pop, cases, source=spec.burden_source or "inline config value"
            )

    return score_target(
        gene=spec.gene,
        variant=spec.variant,
        disease=spec.disease,
        antigen_fraction=antigen,
        incidence_by_population=incidence,
        coverage_by_population=coverage,
        tractability_context=tract,
        allele_source=alleles,
    )


def preflight(spec: TargetSpec, *, client: httpx.Client | None = None) -> list[str]:
    """Return a list of issues that would make the run meaningless. Empty list = OK."""
    issues: list[str] = []

    if resolve_entrez(spec.gene, client=client).is_missing:
        issues.append(f"gene {spec.gene!r} not found in cBioPortal")

    if not check_study(spec.study, client=client).value:
        issues.append(f"study {spec.study!r} not found or has no sequenced samples")

    if spec.freqs_path:
        ft = load_afnd_frequencies(spec.freqs_path)
        norm_alleles = [normalize_allele(a) for a in spec.alleles]
        for pop in spec.populations:
            if pop not in ft.by_population:
                issues.append(f"population {pop!r} absent from frequency file (join would fail)")
                continue
            present = {normalize_allele(k) for k in ft.get(pop)}
            missing = [a for a in norm_alleles if a not in present]
            if missing:
                issues.append(f"alleles with no frequency in {pop!r}: {missing}")
    else:
        issues.append("no frequency file supplied -> HLA coverage will be missing")

    if spec.burden_path:
        loaded = load_burden_table(spec.burden_path, spec.disease)
        for pop in spec.populations:
            if pop not in loaded:
                issues.append(f"no incidence for {pop!r} (disease {spec.disease!r}) in burden file")
    elif not spec.burden_manual:
        issues.append("no burden source supplied -> incidence will be missing")

    if spec.uniprot:
        seq = fetch_uniprot_sequence(spec.uniprot, client=client)
        if seq.is_missing:
            issues.append(f"could not fetch UniProt {spec.uniprot!r} to validate variant WT")
        else:
            try:
                wt, pos, _ = parse_substitution(spec.variant)
                if pos > len(seq.value):
                    issues.append(f"variant {spec.variant} position beyond UniProt {spec.uniprot} length")
                elif seq.value[pos - 1] != wt:
                    issues.append(
                        f"variant {spec.variant} WT mismatch: UniProt {spec.uniprot} has "
                        f"{seq.value[pos - 1]} at position {pos} (wrong isoform/accession?)"
                    )
            except ValueError as exc:
                issues.append(str(exc))

    return issues
