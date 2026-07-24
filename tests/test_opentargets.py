"""Hermetic tests for the Open Targets GraphQL client (mock transport, no network)."""

import json

import httpx

from pmhc_triage.opentargets import (
    associated_targets,
    resolve_disease,
    resolve_target,
    tractability,
)


def _client(payload_for):
    """Build a client whose mock transport returns payload_for(request_body)."""

    def handler(request):
        body = json.loads(request.content)
        return httpx.Response(200, json=payload_for(body))

    return httpx.Client(transport=httpx.MockTransport(handler))


def test_resolve_target_exact_match():
    c = _client(lambda b: {"data": {"search": {"hits": [
        {"id": "ENSG00000133703", "name": "KRAS", "entity": "target"},
        {"id": "ENSG00000220635", "name": "KRASP1", "entity": "target"},
    ]}}})
    s = resolve_target("KRAS", client=c)
    assert s.value == "ENSG00000133703"
    assert not s.warnings  # exact match, no warning


def test_resolve_target_no_exact_warns():
    c = _client(lambda b: {"data": {"search": {"hits": [
        {"id": "ENSG999", "name": "KRASP1", "entity": "target"},
    ]}}})
    s = resolve_target("KRAS", client=c)
    assert s.value == "ENSG999"
    assert any("no exact symbol match" in w for w in s.warnings)


def test_resolve_disease():
    c = _client(lambda b: {"data": {"search": {"hits": [
        {"id": "MONDO_0009831", "name": "malignant pancreatic neoplasm", "entity": "disease"},
    ]}}})
    s = resolve_disease("pancreatic cancer", client=c)
    assert s.value == "MONDO_0009831"


def test_tractability_returns_buckets():
    c = _client(lambda b: {"data": {"target": {
        "approvedSymbol": "KRAS",
        "tractability": [{"label": "Approved Drug", "modality": "SM", "value": True}],
    }}})
    s = tractability("ENSG00000133703", client=c)
    assert s.value[0]["modality"] == "SM"
    assert "context only" in s.provenance.method


def test_associated_targets_parsed():
    c = _client(lambda b: {"data": {"disease": {
        "name": "x",
        "associatedTargets": {"count": 2, "rows": [
            {"score": 0.81, "target": {"id": "ENSG00000133703", "approvedSymbol": "KRAS"}},
            {"score": 0.80, "target": {"id": "ENSG00000141510", "approvedSymbol": "TP53"}},
        ]},
    }}})
    s = associated_targets("MONDO_0009831", top=2, client=c)
    assert [r["symbol"] for r in s.value] == ["KRAS", "TP53"]
    assert s.value[0]["association_score"] == 0.81


def test_http_error_surfaces_missing():
    def handler(request):
        return httpx.Response(500, text="boom")

    s = resolve_target("KRAS", client=httpx.Client(transport=httpx.MockTransport(handler)))
    assert s.is_missing and any("HTTP 500" in w for w in s.warnings)


def test_graphql_errors_surface():
    c = _client(lambda b: {"errors": [{"message": "bad query"}]})
    s = tractability("ENSG_bad", client=c)
    assert s.is_missing and any("GraphQL errors" in w for w in s.warnings)
