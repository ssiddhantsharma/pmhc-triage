"""HLA population coverage: fraction of a population carrying >=1 presenting allele.

Diploid Hardy-Weinberg method of Bui et al. 2006 (the IEDB Population Coverage tool).
For locus L with covering-allele frequency sum p_L, an individual is uncovered at L
iff both alleles are non-covering; assuming independent loci:

    coverage = 1 - prod_L (1 - p_L) ** 2

Frequencies are summed per locus, not across (that would ignore diploidy). Missing
frequencies are excluded and surfaced, never treated as 0.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Iterable, Mapping

from .provenance import Provenance, Sourced, today_iso

_AFND = "Allele Frequency Net Database (AFND)"
_METHOD = "Bui et al. 2006 Hardy-Weinberg diploid population coverage"


def normalize_allele(allele: str) -> str:
    """Canonicalize an allele string for matching.

    Strips an optional ``HLA-`` prefix and surrounding whitespace, leaving e.g.
    ``A*02:01``. Matching between the covering set and the frequency table must go
    through this so ``HLA-A*02:01`` and ``A*02:01`` are treated as the same allele.
    """
    a = allele.strip()
    if a.upper().startswith("HLA-"):
        a = a[4:]
    return a


_CLASS_I_LOCI = {"A", "B", "C"}


def hla_class(allele: str) -> str:
    """Return 'I', 'II', or 'unknown' for an allele (by locus)."""
    try:
        locus = parse_locus(allele)
    except ValueError:
        return "unknown"
    if locus in _CLASS_I_LOCI:
        return "I"
    if locus.startswith(("DR", "DQ", "DP")):
        return "II"
    return "unknown"


def parse_locus(allele: str) -> str:
    """Return the locus (gene) of an allele, e.g. A*02:01 -> A. Raises ValueError if malformed."""
    a = normalize_allele(allele)
    if "*" not in a:
        raise ValueError(f"malformed HLA allele (no '*'): {allele!r}")
    return a.split("*", 1)[0]


def population_coverage(
    covering_alleles: Iterable[str],
    allele_freqs: Mapping[str, float],
    population: str,
    *,
    url: str | None = None,
) -> Sourced[float]:
    """Fraction of ``population`` carrying >=1 of ``covering_alleles``.

    Parameters
    ----------
    covering_alleles
        HLA alleles predicted to present the epitope (any format accepted by
        :func:`normalize_allele`).
    allele_freqs
        Mapping of allele -> frequency for this one population (0..1), e.g. from
        AFND. Frequencies are per-chromosome allele frequencies, not phenotypic.
    population
        Population/cohort name, recorded in the provenance.
    url
        Optional AFND URL for the frequency data, recorded in the provenance.

    Returns
    -------
    Sourced[float]
        Coverage in ``[0, 1]``. ``value is None`` (missing) iff *none* of the
        covering alleles have a known frequency in this population. Alleles absent
        from ``allele_freqs`` are excluded and reported in ``warnings``.
    """
    prov = Provenance(source=_AFND, url=url, query_date=today_iso(), method=_METHOD)

    # drop non-numeric freqs; the allele then reads as missing (surfaced, not a crash)
    norm_freqs: dict[str, float] = {}
    for k, v in allele_freqs.items():
        try:
            norm_freqs[normalize_allele(k)] = float(v)
        except (TypeError, ValueError):
            pass

    # de-dup so "HLA-A*02:01" and "A*02:01" don't both add into the per-locus sum
    covering = list(dict.fromkeys(normalize_allele(a) for a in covering_alleles))

    if not covering:
        return Sourced(None, prov).warn("no covering alleles provided")

    per_locus: dict[str, float] = defaultdict(float)
    loci_with_data: set[str] = set()
    missing: list[str] = []
    malformed: list[str] = []
    invalid: list[str] = []
    for allele in covering:
        try:
            locus = parse_locus(allele)
        except ValueError:
            malformed.append(allele)
            continue
        if allele not in norm_freqs:
            missing.append(allele)
            continue
        f = norm_freqs[allele]
        if not math.isfinite(f) or f < 0.0 or f > 1.0:
            invalid.append(allele)  # exclude NaN/inf/out-of-range, surface it
            continue
        per_locus[locus] += f
        loci_with_data.add(locus)

    result: Sourced[float]
    if not loci_with_data:
        result = Sourced(None, prov).warn(
            f"no frequency data in {population!r} for any covering allele"
        )
    else:
        not_covered = 1.0
        for locus in loci_with_data:
            p = per_locus[locus]
            if p > 1.0:
                p = 1.0  # clamp noisy AFND sums >1 so (1-p) can't go negative
            not_covered *= (1.0 - p) ** 2
        coverage = 1.0 - not_covered
        result = Sourced(round(coverage, 6), prov)
        if any(per_locus[loc] > 1.0 for loc in loci_with_data):
            result.warn("covering-allele frequency sum exceeded 1.0 at a locus; clamped")

    if missing:
        result.warn(
            "excluded (no frequency in "
            f"{population!r}, NOT treated as 0): {', '.join(sorted(missing))}"
        )
    if invalid:
        result.warn(
            f"excluded (non-finite or out-of-range frequency, NOT treated as 0): {', '.join(sorted(invalid))}"
        )
    if malformed:
        result.warn(
            f"excluded (malformed allele, no locus): {', '.join(sorted(malformed))}"
        )
    if any(hla_class(a) == "II" for a in covering):
        result.warn(
            "class II alleles present: DR (via DRB1) is a valid proxy (invariant alpha chain), "
            "but DQ/DP are alpha-beta heterodimers -- coverage over the beta locus alone ignores "
            "alpha-chain pairing and is approximate"
        )
    return result


def coverage_by_population(
    covering_alleles: Iterable[str],
    freqs_by_population: Mapping[str, Mapping[str, float]],
    *,
    urls: Mapping[str, str] | None = None,
) -> dict[str, Sourced[float]]:
    """Run :func:`population_coverage` across several populations.

    ``freqs_by_population`` maps population name -> {allele: freq}. Returns a dict
    keyed by the same population names.
    """
    covering = list(covering_alleles)
    urls = urls or {}
    return {
        pop: population_coverage(covering, freqs, pop, url=urls.get(pop))
        for pop, freqs in freqs_by_population.items()
    }
