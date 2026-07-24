"""Which HLA alleles present the epitope -> feeds HLA coverage.

Two paths:

- **manual** (default, no heavy deps): you supply the presenting alleles (e.g. from
  the literature or IEDB). This is what lets the KRAS example run without
  TensorFlow.
- **MHCflurry** (optional): predict presenting alleles from candidate peptides over
  a reference allele panel. Requires the ``presentation`` extra (pulls TensorFlow).
  If MHCflurry isn't installed, this returns a *surfaced missing* result with an
  install hint -- it never crashes and never silently drops to the manual path.

The peptide->allele thresholding logic (:func:`select_presenting`) is a pure,
tested function independent of MHCflurry, so the selection rule is verifiable even
without the optional dependency installed.
"""

from __future__ import annotations

from typing import Iterable

from .hla import normalize_allele, parse_locus
from .provenance import Provenance, Sourced, today_iso


def manual_presenting_alleles(
    alleles: Iterable[str], *, source: str = "user-specified presenting alleles"
) -> Sourced[list[str]]:
    """Normalize, de-duplicate, and validate a hand-supplied allele set.

    Malformed alleles (no locus) are surfaced as warnings and excluded, rather than
    silently kept and breaking coverage downstream.
    """
    prov = Provenance(source=source, query_date=today_iso(), method="manually supplied allele set")
    valid: list[str] = []
    seen: set[str] = set()
    warnings: list[str] = []
    for raw in alleles:
        a = normalize_allele(raw)
        try:
            parse_locus(a)
        except ValueError:
            warnings.append(f"excluded malformed allele {raw!r}")
            continue
        if a not in seen:
            seen.add(a)
            valid.append(a)
    result = Sourced(valid, prov)
    result.warnings.extend(warnings)
    if not valid:
        result.value = None
        result.warn("no valid alleles supplied")
    return result


def select_presenting(
    df,
    *,
    allele_col: str = "allele",
    score_col: str = "percentile",
    threshold: float = 2.0,
    lower_is_better: bool = True,
) -> list[str]:
    """Pure selection: alleles whose best score over the peptides passes ``threshold``.

    ``df`` is any table with an allele column and a per-(peptide, allele) score
    column. With ``lower_is_better`` (percentile rank), an allele presents if its
    *minimum* score across peptides is ``<= threshold``. Returned alleles are
    normalized and de-duplicated.
    """
    keep: set[str] = set()
    for allele, group in df.groupby(allele_col):
        best = group[score_col].min() if lower_is_better else group[score_col].max()
        passes = best <= threshold if lower_is_better else best >= threshold
        if passes:
            keep.add(normalize_allele(str(allele)))
    return sorted(keep)


def predict_presenting_alleles(
    peptides: Iterable[str],
    allele_panel: Iterable[str],
    *,
    threshold_percentile: float = 2.0,
    _predictor=None,
) -> Sourced[list[str]]:
    """Optional MHCflurry path: predict which panel alleles present any peptide.

    Returns a surfaced *missing* result (never a crash) if MHCflurry is not
    installed or prediction fails. ``_predictor`` may be injected for testing.

    NOTE: the MHCflurry call is written against its documented
    ``Class1PresentationPredictor`` API but has not been live-verified in this
    build (kept dependency-free by default); verify column names against your
    installed MHCflurry version. The selection rule itself is tested via
    :func:`select_presenting`.
    """
    prov = Provenance(
        source="MHCflurry Class1PresentationPredictor",
        query_date=today_iso(),
        method=f"presentation percentile <= {threshold_percentile} over the allele panel",
    )
    peptides = list(peptides)
    allele_panel = list(allele_panel)

    if _predictor is None:
        try:
            from mhcflurry import Class1PresentationPredictor
        except ImportError:
            return Sourced(None, prov).warn(
                "mhcflurry not installed; run `pip install 'pmhc-triage[presentation]'` "
                "(pulls TensorFlow) or supply alleles manually"
            )
        try:
            _predictor = Class1PresentationPredictor.load()
        except Exception as exc:  # pragma: no cover - environment dependent
            return Sourced(None, prov).warn(f"could not load MHCflurry predictor: {exc}")

    try:
        df = _predictor.predict(peptides=peptides, alleles=allele_panel)
    except Exception as exc:  # pragma: no cover - environment dependent
        return Sourced(None, prov).warn(f"MHCflurry prediction failed: {exc}")

    # Normalize to (allele, percentile); MHCflurry uses 'affinity_percentile' (lower=better).
    score_col = "affinity_percentile" if "affinity_percentile" in getattr(df, "columns", []) else "percentile"
    alleles = select_presenting(df, score_col=score_col, threshold=threshold_percentile, lower_is_better=True)
    return Sourced(alleles, prov)
