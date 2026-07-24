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

import httpx

from .provenance import Provenance, Sourced, today_iso

API = "https://www.cbioportal.org/api"
_SOURCE = "cBioPortal REST"


def _request(method, url, client, timeout, *, json=None):
    owns = client is None
    client = client or httpx.Client(timeout=timeout)
    try:
        resp = client.request(method, url, json=json, headers={"Accept": "application/json"})
    except httpx.HTTPError as exc:
        return None, f"cBioPortal request failed: {exc}"
    finally:
        if owns:
            client.close()
    if resp.status_code != 200:
        return None, f"HTTP {resp.status_code} for {url}"
    return resp.json(), None


def resolve_entrez(gene: str, *, client: httpx.Client | None = None, timeout: float = 30.0) -> Sourced[int]:
    """Resolve a HUGO gene symbol to an NCBI Entrez id (e.g. ``KRAS`` -> ``3845``)."""
    url = f"{API}/genes/{gene}"
    prov = Provenance(source=_SOURCE, url=url, query_date=today_iso(),
                      method="cBioPortal /genes lookup -> entrezGeneId")
    data, err = _request("GET", url, client, timeout)
    if err:
        return Sourced(None, prov).warn(err)
    eid = (data or {}).get("entrezGeneId")
    if eid is None:
        return Sourced(None, prov).warn(f"no Entrez id for gene {gene!r}")
    return Sourced(int(eid), prov)


def check_study(study_id: str, *, client: httpx.Client | None = None, timeout: float = 30.0) -> Sourced[bool]:
    """Preflight: does the cBioPortal study exist (and does its _sequenced list resolve)?"""
    url = f"{API}/sample-lists/{study_id}_sequenced"
    prov = Provenance(source=_SOURCE, url=url, query_date=today_iso(),
                      method="cBioPortal study preflight (_sequenced sample list)")
    data, err = _request("GET", url, client, timeout)
    if err:
        return Sourced(False, prov).warn(err)
    n = len((data or {}).get("sampleIds", []))
    if n == 0:
        return Sourced(False, prov).warn(f"{study_id!r} has no sequenced samples")
    return Sourced(True, prov)


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

    # denominator: number of sequenced samples
    sl_data, err = _request("GET", f"{API}/sample-lists/{sample_list}", client, timeout)
    if err:
        return Sourced(None, prov(context)).warn(err)
    denom = len((sl_data or {}).get("sampleIds", []))
    if denom == 0:
        return Sourced(None, prov(context)).warn(
            f"sample list {sample_list!r} is empty or missing (denominator 0)"
        )

    # numerator: distinct samples with the exact protein change
    body = {"sampleListId": sample_list, "entrezGeneIds": [entrez]}
    muts, err = _request(
        "POST", f"{API}/molecular-profiles/{profile}/mutations/fetch", client, timeout, json=body
    )
    if err:
        return Sourced(None, prov(context)).warn(err)
    positive = {m["sampleId"] for m in (muts or []) if m.get("proteinChange") == protein_change}
    numer = len(positive)
    fraction = round(numer / denom, 6)

    result = Sourced(
        fraction,
        prov(
            f"{protein_change}: {numer}/{denom} sequenced samples in {study_id} "
            f"(cBioPortal ODbL; profile {profile})"
        ),
    )
    if numer == 0:
        result.warn(
            f"zero samples with {protein_change} in {gene} for {study_id} "
            "(a real 0.0, not missing -- check the study/variant if unexpected)"
        )
    return result
