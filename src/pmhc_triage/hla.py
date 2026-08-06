"""HLA population coverage -- the load-bearing module.

Given the set of HLA alleles that present an epitope (from ``presentation``) and
per-population allele frequencies (from AFND), compute the fraction of a population
that carries at least one presenting allele, using the diploid Hardy-Weinberg
method of Bui et al. 2006 (Immunogenetics; the method behind the IEDB Population
Coverage tool).

Method
------
For one locus ``L`` whose covering alleles have frequency sum ``p_L`` (the chance a
random chromosome carries a covering allele at ``L``), an individual is *not*
covered at ``L`` iff *both* of their alleles are non-covering::

    P(not covered at L) = (1 - p_L) ** 2        # Hardy-Weinberg, diploid

Assuming loci are independent (the standard simplification), coverage combines as::

    P(not covered anywhere) = prod_L (1 - p_L) ** 2
    coverage                = 1 - prod_L (1 - p_L) ** 2

Two failure modes this module refuses to commit
------------------------------------------------
1. Do **not** sum allele frequencies as a proxy for coverage -- that ignores
   diploidy and multi-locus combination and overstates reach.
2. Do **not** treat a missing allele frequency as ``0`` -- that silently
   understates coverage and can make a real target look unreachable. Missing
   frequencies are excluded from the maths *and surfaced as warnings*.

The result is returned as a :class:`~pmhc_triage.provenance.Sourced` float so the
provenance (AFND, population, method, query date) and any warnings travel with it.
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
    """Return the locus (gene) of an allele, e.g. ``A*02:01`` -> ``A``.

    Raises ``ValueError`` if the allele has no ``*`` separating gene from field --
    that indicates a malformed input, which we surface rather than guess at.
    """
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

    # Parse frequencies defensively: a non-numeric value is dropped (the allele then
    # reads as "no frequency" and is surfaced), never crashes the computation.
    norm_freqs: dict[str, float] = {}
    for k, v in allele_freqs.items():
        try:
            norm_freqs[normalize_allele(k)] = float(v)
        except (TypeError, ValueError):
            pass

    # De-duplicate covering alleles (order-preserving). An allele listed twice -- e.g.
    # "HLA-A*02:01" and "A*02:01", which normalize to the same -- must NOT double-count
    # its frequency into the per-locus sum (that silently inflates coverage).
    covering = list(dict.fromkeys(normalize_allele(a) for a in covering_alleles))

    if not covering:
        return Sourced(None, prov).warn("no covering alleles provided")

    # Sum covering-allele frequencies per locus; track missing / malformed / invalid.
    per_locus: dict[str, float] = defaultdict(float)
    loci_with_data: set[str] = set()
    missing: list[str] = []
    malformed: list[str] = []
    invalid: list[str] = []
    for allele in covering:
        try:
            locus = parse_locus(allele)
        except ValueError:
            malformed.append(allele)  # surface, don't crash the whole run for one bad allele
            continue
        if allele not in norm_freqs:
            missing.append(allele)
            continue
        f = norm_freqs[allele]
        if not math.isfinite(f) or f < 0.0 or f > 1.0:
            invalid.append(allele)  # NaN/inf/out-of-range freq -> excluded + surfaced, never laundered
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
                # Noisy/overlapping AFND data can push a sum slightly over 1; clamp
                # and surface rather than silently produce a negative (1-p).
                p = 1.0
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
