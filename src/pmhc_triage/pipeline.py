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

from dataclasses import dataclass
from typing import Any

import httpx

from .afnd import load_afnd_frequencies
from .antigen import (
    check_study,
    expression_positive_fraction,
    resolve_entrez,
    study_cancer_type,
    variant_frequency,
    variant_frequency_multi,
)
from .burden import load_bundled, load_burden_table, manual_incidence
from .hla import coverage_by_population, normalize_allele
from .identity import disease_matches, suggest_match
from .montecarlo import MCResult, effective_n_interval
from .opentargets import resolve_target, tractability
from .peptides import mutant_peptides, parse_substitution
from .presentation import manual_presenting_alleles, predict_presenting_alleles
from .provenance import Provenance, Sourced
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
    predict_alleles: bool = False
    presentation_threshold: float = 2.0
    studies: list[str] | None = None   # pool antigen fraction across these (default [study])
    variants: list[str] | None = None  # pool across these variants (default [variant])
    antigen_mode: str = "mutation"     # "mutation" (cBioPortal variant) or "expression" (RNA-seq)
    expression_threshold: float = 1.0  # z-score cutoff for expression mode
    mc: bool = True                    # Monte-Carlo interval on effective_N (antigen+coverage)
    mc_draws: int = 20000              # MC draws
    mc_seed: int = 0                   # MC seed (reproducible intervals)
    incidence_rel_sd: float | None = None  # optional USER-owned incidence relative SD

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "TargetSpec":
        return cls(
            gene=d["gene"],
            variant=d.get("variant", ""),
            disease=d["disease"],
            study=d["study"],
            alleles=list(d.get("alleles", [])),
            populations=list(d["populations"]),
            freqs_path=d.get("freqs"),
            burden_path=d.get("burden") if isinstance(d.get("burden"), str) else None,
            burden_manual=d.get("burden") if isinstance(d.get("burden"), dict) else None,
            burden_source=d.get("burden_source"),
            uniprot=d.get("uniprot"),
            predict_alleles=bool(d.get("predict_alleles", False)),
            presentation_threshold=float(d.get("presentation_threshold", 2.0)),
            studies=list(d["studies"]) if d.get("studies") else None,
            variants=list(d["variants"]) if d.get("variants") else None,
            antigen_mode=d.get("antigen_mode", "mutation"),
            expression_threshold=float(d.get("expression_threshold", 1.0)),
            mc=bool(d.get("mc", True)),
            mc_draws=int(d.get("mc_draws", 20000)),
            mc_seed=int(d.get("mc_seed", 0)),
            incidence_rel_sd=(
                float(d["incidence_rel_sd"]) if d.get("incidence_rel_sd") is not None else None
            ),
        )


def _predicted_alleles(spec: "TargetSpec", freqs_by_pop: dict, client) -> Sourced:
    """Predict presenting alleles via MHCflurry: UniProt seq -> peptides -> predict over
    the frequency file's allele panel. Any failure surfaces (missing), never fabricates."""
    if not spec.uniprot:
        return Sourced(None, Provenance(source="allele prediction")).warn(
            "predict_alleles requires a uniprot accession"
        )
    seq = fetch_uniprot_sequence(spec.uniprot, client=client)
    if seq.is_missing:
        return Sourced(None, seq.provenance).warn("UniProt fetch failed; cannot predict alleles")
    peps = mutant_peptides(seq.value, spec.variant)
    if peps.is_missing:
        return Sourced(None, peps.provenance, warnings=list(peps.warnings))
    panel = sorted({a for pop in freqs_by_pop.values() for a in pop})
    if not panel:
        return Sourced(None, Provenance(source="allele prediction")).warn(
            "no alleles in frequency file to predict over"
        )
    return predict_presenting_alleles(peps.value, panel, threshold_percentile=spec.presentation_threshold)


def run_target(spec: TargetSpec, *, client: httpx.Client | None = None) -> TargetScore:
    """Fetch every factor and combine via the score gate. Missing factors surface, never fake."""
    if spec.antigen_mode == "expression":
        antigen = expression_positive_fraction(
            spec.study, spec.gene, threshold=spec.expression_threshold, client=client
        )
    else:
        studies = spec.studies or [spec.study]
        variants = spec.variants or [spec.variant]
        if len(studies) > 1 or len(variants) > 1:
            antigen = variant_frequency_multi(studies, variants, gene=spec.gene, client=client)
        else:
            antigen = variant_frequency(spec.study, spec.variant, gene=spec.gene, client=client)

    tract = None
    tid = resolve_target(spec.gene, client=client)
    if not tid.is_missing:
        tract = tractability(tid.value, client=client)

    if spec.freqs_path:
        ft = load_afnd_frequencies(spec.freqs_path)
        freqs_by_pop = {p: ft.get(p) for p in spec.populations}
        sizes_by_pop = {p: ft.get_sample_sizes(p) for p in spec.populations}
    else:
        freqs_by_pop = {p: {} for p in spec.populations}
        sizes_by_pop = {p: {} for p in spec.populations}

    if spec.predict_alleles:
        alleles = _predicted_alleles(spec, freqs_by_pop, client)
    else:
        alleles = manual_presenting_alleles(spec.alleles)
    allele_list = alleles.value or []
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
    else:
        # No user burden supplied -> fall back to the shipped cited starter bundle
        # (World + per-region country-proxy figures). Missing disease/population simply
        # doesn't populate that key, and the score gate surfaces it as missing (never 0).
        for pop in spec.populations:
            s = load_bundled(spec.disease, population=pop)
            if not s.is_missing:
                incidence[pop] = s

    mc_by_pop = _monte_carlo(spec, antigen, incidence, coverage, allele_list, freqs_by_pop, sizes_by_pop)

    return score_target(
        gene=spec.gene,
        variant=spec.variant,
        disease=spec.disease,
        antigen_fraction=antigen,
        incidence_by_population=incidence,
        coverage_by_population=coverage,
        tractability_context=tract,
        allele_source=alleles,
        mc_by_population=mc_by_pop,
    )


def _monte_carlo(spec, antigen, incidence, coverage, allele_list, freqs_by_pop, sizes_by_pop):
    """Per-population Monte-Carlo interval, only where all factors are present.

    Runs only for populations whose incidence, coverage, and antigen fraction are all
    non-missing (the same gate the score step uses to emit a point estimate) and where
    the antigen fraction carries numerator/denominator. Never fabricates an interval.
    """
    if not spec.mc or antigen.is_missing:
        return {}
    numer = antigen.extra.get("numerator")
    denom = antigen.extra.get("denominator")
    if numer is None or not denom:
        return {}
    out: dict[str, MCResult] = {}
    for pop in spec.populations:
        inc = incidence.get(pop)
        cov = coverage.get(pop)
        if inc is None or inc.is_missing or cov is None or cov.is_missing:
            continue
        out[pop] = effective_n_interval(
            incidence=inc.value,
            antigen_numerator=int(numer),
            antigen_denominator=int(denom),
            covering_alleles=allele_list,
            allele_freqs=freqs_by_pop.get(pop, {}),
            sample_sizes=sizes_by_pop.get(pop, {}),
            n_draws=spec.mc_draws,
            seed=spec.mc_seed,
            incidence_rel_sd=spec.incidence_rel_sd,
        )
    return out


def preflight(spec: TargetSpec, *, client: httpx.Client | None = None) -> list[str]:
    """Return a list of issues that would make the run meaningless. Empty list = OK."""
    issues: list[str] = []

    if resolve_entrez(spec.gene, client=client).is_missing:
        issues.append(f"gene {spec.gene!r} not found in cBioPortal")

    if not check_study(spec.study, client=client).value:
        issues.append(f"study {spec.study!r} not found or has no sequenced samples")
    else:
        ct = study_cancer_type(spec.study, client=client)
        if not ct.is_missing and not disease_matches(spec.disease, ct.value):
            issues.append(
                f"study {spec.study!r} is {ct.value!r} but --disease is {spec.disease!r} "
                "-- study and disease may not correspond (wrong study for this neoantigen?)"
            )

    if spec.freqs_path:
        ft = load_afnd_frequencies(spec.freqs_path)
        norm_alleles = [normalize_allele(a) for a in spec.alleles]
        for pop in spec.populations:
            if pop not in ft.by_population:
                hint = suggest_match(pop, ft.by_population.keys())
                suffix = f" -- did you mean {hint!r}? (label alias)" if hint else ""
                issues.append(f"population {pop!r} absent from frequency file (join would fail){suffix}")
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
                hint = suggest_match(pop, loaded.keys())
                suffix = f" -- did you mean {hint!r}? (label alias)" if hint else ""
                issues.append(f"no incidence for {pop!r} (disease {spec.disease!r}) in burden file{suffix}")
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
