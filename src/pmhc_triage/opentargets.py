"""Thin Open Targets GraphQL client -- consumed as CONTEXT, never a fork.

Guardrail: pmhc-triage does **not** reimplement Open Targets' association / L2G
engine. It calls the public GraphQL API (data is CC0) for two things that feed
*around* the pMHC/HLA spine:

- ``tractability`` -- small-molecule / antibody druggability buckets. **Context
  only**: it is the wrong modality for pMHC and must never drive the addressable-N
  score.
- ``associated_targets`` -- top targets for a disease, powering the optional
  ``discovery`` front-end that proposes candidate genes to run through the pipeline.

Endpoint verified live: https://api.platform.opentargets.org/api/v4/graphql
Citing Open Targets is courtesy (CC0).
"""

from __future__ import annotations

from dataclasses import replace

import httpx

from .provenance import Provenance, Sourced, fetched_at_or_today, today_iso

ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
_SOURCE = "Open Targets Platform GraphQL"


def _post(
    query: str,
    variables: dict,
    client: httpx.Client | None,
    timeout: float,
) -> tuple[dict | None, str | None, str | None]:
    """POST a GraphQL query. Returns ``(data, error, fetched_at)`` (fetched_at None on error)."""
    owns = client is None
    client = client or httpx.Client(timeout=timeout)
    try:
        resp = client.post(ENDPOINT, json={"query": query, "variables": variables})
    except httpx.HTTPError as exc:
        return None, f"request to Open Targets failed: {exc}", None
    finally:
        if owns:
            client.close()
    if resp.status_code != 200:
        return None, f"HTTP {resp.status_code} from Open Targets", None
    body = resp.json()
    if body.get("errors"):
        return None, f"GraphQL errors: {body['errors']}", None
    return body.get("data"), None, fetched_at_or_today(resp)


def resolve_target(symbol: str, *, client: httpx.Client | None = None, timeout: float = 30.0) -> Sourced[str]:
    """Resolve an approved gene symbol (e.g. ``KRAS``) to an Ensembl gene id."""
    prov = Provenance(
        source=_SOURCE,
        url=ENDPOINT,
        query_date=today_iso(),
        method=f"search(entityNames=[target]) for {symbol!r} -> Ensembl gene id",
    )
    q = 'query($q:String!){ search(queryString:$q, entityNames:["target"]){ hits { id name entity } } }'
    data, err, fa = _post(q, {"q": symbol}, client, timeout)
    if err:
        return Sourced(None, prov).warn(err)
    hits = (data or {}).get("search", {}).get("hits", [])
    exact = [h for h in hits if h.get("name", "").upper() == symbol.upper()]
    chosen = (exact or hits or [None])[0]
    if not chosen:
        return Sourced(None, prov).warn(f"no target hit for {symbol!r}")
    result = Sourced(chosen["id"], replace(prov, query_date=fa))
    if not exact:
        result.warn(f"no exact symbol match; using top hit {chosen.get('name')} ({chosen['id']})")
    return result


def resolve_disease(name: str, *, client: httpx.Client | None = None, timeout: float = 30.0) -> Sourced[str]:
    """Resolve a disease name to an EFO/MONDO id (top search hit)."""
    prov = Provenance(
        source=_SOURCE,
        url=ENDPOINT,
        query_date=today_iso(),
        method=f"search(entityNames=[disease]) for {name!r} -> EFO/MONDO id",
    )
    q = 'query($q:String!){ search(queryString:$q, entityNames:["disease"]){ hits { id name entity } } }'
    data, err, fa = _post(q, {"q": name}, client, timeout)
    if err:
        return Sourced(None, prov).warn(err)
    hits = (data or {}).get("search", {}).get("hits", [])
    if not hits:
        return Sourced(None, prov).warn(f"no disease hit for {name!r}")
    top = hits[0]
    result = Sourced(top["id"], replace(prov, query_date=fa))
    if top.get("name", "").lower() != name.lower():
        result.warn(f"using top hit {top.get('name')!r} ({top['id']}) for query {name!r}")
    return result


def tractability(ensembl_id: str, *, client: httpx.Client | None = None, timeout: float = 30.0) -> Sourced[list[dict]]:
    """Tractability buckets for a target. CONTEXT ONLY -- never a pMHC-score driver."""
    prov = Provenance(
        source=_SOURCE,
        url=ENDPOINT,
        query_date=today_iso(),
        method="target.tractability (SM/AB buckets; context only, not a pMHC-score driver)",
    )
    q = 'query($id:String!){ target(ensemblId:$id){ approvedSymbol tractability { label modality value } } }'
    data, err, fa = _post(q, {"id": ensembl_id}, client, timeout)
    if err:
        return Sourced(None, prov).warn(err)
    target = (data or {}).get("target")
    if not target:
        return Sourced(None, prov).warn(f"no target for {ensembl_id!r}")
    return Sourced(target.get("tractability", []), replace(prov, query_date=fa))


def associated_targets(
    disease_id: str,
    *,
    top: int = 20,
    client: httpx.Client | None = None,
    timeout: float = 30.0,
) -> Sourced[list[dict]]:
    """Top associated targets for a disease -- the discovery front-end (chained, not cloned).

    Returns rows of ``{ensembl_id, symbol, association_score}``. This is Open
    Targets' association scoring, consumed verbatim; pmhc-triage does not recompute it.
    """
    prov = Provenance(
        source=_SOURCE,
        url=ENDPOINT,
        query_date=today_iso(),
        method=f"disease.associatedTargets top {top} (Open Targets association score, consumed verbatim)",
    )
    q = (
        "query($id:String!,$n:Int!){ disease(efoId:$id){ name "
        "associatedTargets(page:{index:0,size:$n}){ count rows { score target { id approvedSymbol } } } } }"
    )
    data, err, fa = _post(q, {"id": disease_id, "n": top}, client, timeout)
    if err:
        return Sourced(None, prov).warn(err)
    disease = (data or {}).get("disease")
    if not disease:
        return Sourced(None, prov).warn(f"no disease for {disease_id!r}")
    rows = disease.get("associatedTargets", {}).get("rows", [])
    parsed = [
        {
            "ensembl_id": r["target"]["id"],
            "symbol": r["target"]["approvedSymbol"],
            "association_score": r["score"],
        }
        for r in rows
    ]
    return Sourced(parsed, replace(prov, query_date=fa))
