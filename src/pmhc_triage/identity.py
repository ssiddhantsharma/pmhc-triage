"""Identity normalization for the join hazard: population + disease label matching.

The single biggest silent-error risk in this tool is a **label mismatch** across
sources -- burden's "European", the frequency file's "Europe", and a requested
"europe" are the same thing, but string-compared they don't join, and the product
silently drops. This module canonicalizes only *true synonyms* (spelling / case /
spacing / a short curated alias list). It deliberately does **not** map a country
to a region or an ethnicity to a population -- those are modeling choices the user
must make explicitly, not silent merges.

The philosophy stays "surface, don't silently fix": these helpers are used to
*warn and suggest* (e.g. "'European' in freqs matches requested 'Europe'"), not to
rewrite the scoring path behind the user's back.
"""

from __future__ import annotations

from typing import Iterable

# Conservative: only obvious synonyms of the SAME label. No country->region merges.
POPULATION_ALIASES = {
    "europe": "Europe",
    "european": "Europe",
    "east asia": "EastAsia",
    "east asian": "EastAsia",
    "eastern asia": "EastAsia",
    "south asia": "SouthAsia",
    "south asian": "SouthAsia",
    "africa": "Africa",
    "african": "Africa",
    "north america": "NorthAmerica",
    "north american": "NorthAmerica",
    "world": "World",
    "global": "World",
}

DISEASE_ALIASES = {
    "pdac": "pancreatic ductal adenocarcinoma",
    "crc": "colorectal carcinoma",
    "nsclc": "non-small cell lung carcinoma",
    "hcc": "hepatocellular carcinoma",
    "aml": "acute myeloid leukemia",
}


def _key(name: str) -> str:
    return " ".join(str(name).strip().lower().split())


def canonical_population(name: str) -> str:
    """Canonical population label; known synonyms mapped, otherwise the trimmed input."""
    return POPULATION_ALIASES.get(_key(name), str(name).strip())


def canonical_disease(name: str) -> str:
    """Canonical disease label; known abbreviations expanded, otherwise the trimmed input."""
    return DISEASE_ALIASES.get(_key(name), str(name).strip())


def suggest_match(target: str, candidates: Iterable[str]) -> str | None:
    """Return a candidate that canonically matches ``target`` (via aliasing), else None.

    Used to turn an unhelpful "population not found" into "did you mean 'European'?".
    """
    ct = canonical_population(target)
    for cand in candidates:
        if canonical_population(cand) == ct and str(cand).strip() != str(target).strip():
            return cand
    return None


def align_populations(a: Iterable[str], b: Iterable[str]) -> dict[str, list[str]]:
    """Compare two population key sets by canonical form.

    Returns ``{"matched": [...], "only_a": [...], "only_b": [...]}`` where matched
    lists the canonical labels present in both.
    """
    ca = {canonical_population(x): x for x in a}
    cb = {canonical_population(x): x for x in b}
    return {
        "matched": sorted(set(ca) & set(cb)),
        "only_a": sorted(set(ca) - set(cb)),
        "only_b": sorted(set(cb) - set(ca)),
    }
