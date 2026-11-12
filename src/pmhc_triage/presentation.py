"""Which HLA alleles present the epitope -> feeds HLA coverage.

manual (default): you supply presenting alleles (literature/IEDB). MHCflurry (optional
[presentation] extra): predict from candidate peptides over an allele panel; if not
installed it returns a surfaced-missing result with an install hint, never crashes. The
threshold selection rule (select_presenting) is a pure, tested function.
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


def sweep_presenting_alleles(
    df,
    thresholds,
    *,
    allele_col: str = "sample_name",
    score_col: str = "presentation_percentile",
    lower_is_better: bool = True,
) -> dict[float, list[str]]:
    """Presenting-allele set at each threshold, making the threshold dependence visible.

    The presentation threshold is a lever, not a fact: a looser cutoff admits more
    alleles and inflates HLA coverage. Rather than commit to one arbitrary point, sweep
    it. ``df`` is a scored (peptide, allele) table (e.g. MHCflurry output); returns
    ``{threshold: [alleles]}``. Pure; no MHCflurry needed to test the selection rule.
    """
    return {
        float(t): select_presenting(
            df,
            allele_col=allele_col,
            score_col=score_col,
            threshold=t,
            lower_is_better=lower_is_better,
        )
        for t in thresholds
    }


def _load_predictor(_predictor=None):
    """Return ``(predictor, error)``. Shared by the allele- and peptide-oriented queries so
    the import/load failure modes are handled in exactly one place."""
    if _predictor is not None:
        return _predictor, None
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError:
        return None, ("mhcflurry not installed; run `pip install 'pmhc-triage[presentation]'` "
                      "(pulls TensorFlow) or supply alleles manually")
    try:
        return Class1PresentationPredictor.load(), None
    except Exception as exc:  # pragma: no cover - environment dependent
        return None, f"could not load MHCflurry predictor: {exc}"


def presented_peptides(
    peptides: Iterable[str],
    allele: str,
    *,
    threshold_percentile: float = 2.0,
    _predictor=None,
) -> Sourced[set[str]]:
    """Which of ``peptides`` the ONE ``allele`` is predicted to present.

    The transpose of ``predict_presenting_alleles``: that asks which alleles present a
    peptide, this asks which peptides an allele presents. Used to tell a self-scan hit that
    could actually be displayed from one that merely matches on sequence.

    Returns a surfaced *missing* result (value None) when MHCflurry is unavailable, never a
    crash and never an empty set that would read as "nothing is presented".
    """
    prov = Provenance(
        source="MHCflurry Class1PresentationPredictor",
        query_date=today_iso(),
        method=f"presentation percentile <= {threshold_percentile} for {allele}",
    )
    peptides = list(dict.fromkeys(peptides))
    if not peptides:
        return Sourced(set(), prov)

    predictor, err = _load_predictor(_predictor)
    if predictor is None:
        return Sourced(None, prov).warn(err)
    try:
        df = predictor.predict(peptides=peptides, alleles={allele: [allele]})
    except Exception as exc:  # pragma: no cover - environment dependent
        return Sourced(None, prov).warn(f"MHCflurry prediction failed: {exc}")

    col = "presentation_percentile"
    if col not in df or "peptide" not in df:
        return Sourced(None, prov).warn(f"MHCflurry output missing '{col}'/'peptide' column")
    keep = df[df[col] <= threshold_percentile]["peptide"]
    return Sourced(set(keep.tolist()), prov)


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

    Live-verified against MHCflurry 2.2.1: each panel allele is passed as its own
    single-allele "sample" (``{allele: [allele]}``) so we get a per-allele
    ``presentation_percentile`` (lower = stronger); an allele presents if its best
    (minimum) percentile over the peptides is ``<= threshold_percentile``.
    """
    prov = Provenance(
        source="MHCflurry Class1PresentationPredictor",
        query_date=today_iso(),
        method=f"presentation percentile <= {threshold_percentile} over the allele panel",
    )
    peptides = list(peptides)
    allele_panel = list(allele_panel)

    _predictor, err = _load_predictor(_predictor)
    if _predictor is None:
        return Sourced(None, prov).warn(err)

    try:
        # one allele per "sample" -> per-allele presentation_percentile
        df = _predictor.predict(peptides=peptides, alleles={a: [a] for a in allele_panel})
    except Exception as exc:  # pragma: no cover - environment dependent
        return Sourced(None, prov).warn(f"MHCflurry prediction failed: {exc}")

    alleles = select_presenting(
        df,
        allele_col="sample_name",
        score_col="presentation_percentile",
        threshold=threshold_percentile,
        lower_is_better=True,
    )
    return Sourced(alleles, prov)
