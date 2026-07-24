"""Antigen-positive fraction from cBioPortal (public mutation data; ODbL).

For a missense neoantigen (e.g. KRAS G12D), the "antigen-positive fraction" is the
share of profiled tumors of a given type carrying that exact protein change:

    numerator   = distinct sequenced samples with proteinChange == variant
    denominator = size of the study's ``*_sequenced`` sample list
    fraction    = numerator / denominator

Computed live from the cBioPortal REST API (verified endpoints). Data is ODbL
(attribution + share-alike); some individual studies additionally restrict
commercial use, noted per study in cBioPortal and left to the user to check. We
query at runtime and never bundle cBioPortal data.

Note the scope boundary from the spec: this clean path is for **mutation**
antigens. Expression antigens (e.g. PRAME) need an RNA-seq threshold and a
different source -- out of this function; supply that fraction manually.
"""

from __future__ import annotations

from dataclasses import replace

import httpx

from .provenance import Provenance, Sourced, fetched_at_or_today, today_iso
from .uncertainty import wilson_ci

API = "https://www.cbioportal.org/api"
_SOURCE = "cBioPortal REST"


def _request(method, url, client, timeout, *, json=None):
    owns = client is None
    client = client or httpx.Client(timeout=timeout)
    try:
        resp = client.request(method, url, json=json, headers={"Accept": "application/json"})
    except httpx.HTTPError as exc:
        return None, f"cBioPortal request failed: {exc}", None
    finally:
        if owns:
            client.close()
    if resp.status_code != 200:
        return None, f"HTTP {resp.status_code} for {url}", None
    return resp.json(), None, fetched_at_or_today(resp)


def resolve_entrez(gene: str, *, client: httpx.Client | None = None, timeout: float = 30.0) -> Sourced[int]:
    """Resolve a HUGO gene symbol to an NCBI Entrez id (e.g. ``KRAS`` -> ``3845``)."""
    url = f"{API}/genes/{gene}"
    prov = Provenance(source=_SOURCE, url=url, query_date=today_iso(),
                      method="cBioPortal /genes lookup -> entrezGeneId")
    data, err, fa = _request("GET", url, client, timeout)
    if err:
        return Sourced(None, prov).warn(err)
    eid = (data or {}).get("entrezGeneId")
    if eid is None:
        return Sourced(None, prov).warn(f"no Entrez id for gene {gene!r}")
    return Sourced(int(eid), replace(prov, query_date=fa))


def check_study(study_id: str, *, client: httpx.Client | None = None, timeout: float = 30.0) -> Sourced[bool]:
    """Preflight: does the cBioPortal study exist (and does its _sequenced list resolve)?"""
    url = f"{API}/sample-lists/{study_id}_sequenced"
    prov = Provenance(source=_SOURCE, url=url, query_date=today_iso(),
                      method="cBioPortal study preflight (_sequenced sample list)")
    data, err, fa = _request("GET", url, client, timeout)
    if err:
        return Sourced(False, prov).warn(err)
    n = len((data or {}).get("sampleIds", []))
    if n == 0:
        return Sourced(False, prov).warn(f"{study_id!r} has no sequenced samples")
    return Sourced(True, replace(prov, query_date=fa))


def study_cancer_type(study_id: str, *, client: httpx.Client | None = None, timeout: float = 30.0) -> Sourced[str]:
    """The cancer type a cBioPortal study represents (for study<->disease preflight)."""
    url = f"{API}/studies/{study_id}"
    prov = Provenance(source=_SOURCE, url=url, query_date=today_iso(),
                      method="cBioPortal study cancerType")
    data, err, fa = _request("GET", url, client, timeout)
    if err:
        return Sourced(None, prov).warn(err)
    ct = (data or {}).get("cancerType") or {}
    name = ct.get("name") or (data or {}).get("cancerTypeId")
    if not name:
        return Sourced(None, prov).warn(f"no cancerType for study {study_id!r}")
    return Sourced(str(name), replace(prov, query_date=fa))


def _counts(study_id, changes, entrez, sample_list, client, timeout):
    """(numerator, denominator, fetched_at, error) for samples carrying ANY of ``changes``.

    numerator = distinct sequenced samples with proteinChange in the set; denominator
    = size of the sequenced sample list. Shared by single- and multi-study paths.
    """
    profile = f"{study_id}_mutations"
    sl_data, err, _ = _request("GET", f"{API}/sample-lists/{sample_list}", client, timeout)
    if err:
        return 0, 0, None, err
    denom = len((sl_data or {}).get("sampleIds", []))
    if denom == 0:
        return 0, 0, None, f"sample list {sample_list!r} is empty or missing (denominator 0)"
    body = {"sampleListId": sample_list, "entrezGeneIds": [entrez]}
    muts, err, fa = _request(
        "POST", f"{API}/molecular-profiles/{profile}/mutations/fetch", client, timeout, json=body
    )
    if err:
        return 0, denom, None, err
    numer = len({m["sampleId"] for m in (muts or []) if m.get("proteinChange") in changes})
    return numer, denom, fa, None


def variant_frequency(
    study_id: str,
    protein_change: str,
    *,
    gene: str = "KRAS",
    entrez: int | None = None,
    sample_list: str | None = None,
    client: httpx.Client | None = None,
    timeout: float = 30.0,
) -> Sourced[float]:
    """Fraction of a study's sequenced tumors carrying ``protein_change`` in ``gene``.

    Returns a :class:`~pmhc_triage.provenance.Sourced` float. A genuine *zero*
    (gene sequenced but no samples carry the change) returns ``0.0`` with a note --
    distinct from *missing* (``None``), which means we could not compute it (bad
    study id, empty sample list, API error). The numerator/denominator are recorded
    in the provenance method so the fraction is auditable, never a bare number.
    """
    profile = f"{study_id}_mutations"
    sample_list = sample_list or f"{study_id}_sequenced"
    context = f"{protein_change} in {gene} across {sample_list}"

    def prov(method: str) -> Provenance:
        return Provenance(source=_SOURCE, url=API, query_date=today_iso(), method=method)

    if entrez is None:
        resolved = resolve_entrez(gene, client=client, timeout=timeout)
        if resolved.is_missing:
            return Sourced(None, prov(context)).warn(
                f"could not resolve Entrez id for {gene!r}: {resolved.warnings}"
            )
        entrez = resolved.value

    numer, denom, fa, err = _counts(study_id, {protein_change}, entrez, sample_list, client, timeout)
    if err:
        return Sourced(None, prov(context)).warn(err)
    fraction = round(numer / denom, 6)

    ci_low, ci_high = wilson_ci(numer, denom)
    result = Sourced(
        fraction,
        replace(
            prov(
                f"{protein_change}: {numer}/{denom} sequenced samples in {study_id} "
                f"(cBioPortal ODbL; profile {profile})"
            ),
            query_date=fa,
        ),
    )
    result.extra = {
        "numerator": numer,
        "denominator": denom,
        "ci95_low": round(ci_low, 6),
        "ci95_high": round(ci_high, 6),
    }
    if numer == 0:
        result.warn(
            f"zero samples with {protein_change} in {gene} for {study_id} "
            "(a real 0.0, not missing -- check the study/variant if unexpected)"
        )
    return result


def variant_frequency_multi(
    studies,
    protein_changes,
    *,
    gene: str = "KRAS",
    entrez: int | None = None,
    client: httpx.Client | None = None,
    timeout: float = 30.0,
) -> Sourced[float]:
    """Pooled antigen-positive fraction across multiple studies and/or variants.

    numerator = distinct samples carrying ANY of ``protein_changes``, pooled across
    ``studies``; denominator = pooled sequenced samples. Per-study breakdown is in
    ``extra["per_study"]``. Studies that error are **surfaced** (listed, skipped),
    never silently dropped. Pooling assumes comparable cohorts -- flagged as a warning.
    """
    changes = set(protein_changes)
    label = "/".join(sorted(changes))
    studies = list(studies)

    def prov(method):
        return Provenance(source=_SOURCE, url=API, query_date=today_iso(), method=method)

    context = f"{label} in {gene} pooled across {len(studies)} studies"
    if entrez is None:
        resolved = resolve_entrez(gene, client=client, timeout=timeout)
        if resolved.is_missing:
            return Sourced(None, prov(context)).warn(f"could not resolve Entrez id for {gene!r}")
        entrez = resolved.value

    per_study, skipped = [], []
    total_n = total_d = 0
    last_fa = None
    for s in studies:
        numer, denom, fa, err = _counts(s, changes, entrez, f"{s}_sequenced", client, timeout)
        if err:
            skipped.append(f"{s}: {err}")
            continue
        per_study.append({"study": s, "numerator": numer, "denominator": denom,
                          "fraction": round(numer / denom, 6)})
        total_n += numer
        total_d += denom
        last_fa = fa or last_fa

    if total_d == 0:
        r = Sourced(None, prov(context)).warn("no usable studies (all errored or empty)")
        for sk in skipped:
            r.warn(f"skipped {sk}")
        return r

    fraction = round(total_n / total_d, 6)
    lo, hi = wilson_ci(total_n, total_d)
    result = Sourced(
        fraction,
        replace(prov(f"{label}: {total_n}/{total_d} pooled across {len(per_study)} studies "
                     "(cBioPortal ODbL)"), query_date=last_fa or today_iso()),
    )
    result.extra = {
        "numerator": total_n,
        "denominator": total_d,
        "ci95_low": round(lo, 6),
        "ci95_high": round(hi, 6),
        "variants": sorted(changes),
        "per_study": per_study,
    }
    result.warn("pooled across studies; assumes comparable cohorts (heterogeneity not modeled)")
    for sk in skipped:
        result.warn(f"skipped {sk}")
    if total_n == 0:
        result.warn(f"zero samples with {label} across studies (a real 0.0)")
    return result
