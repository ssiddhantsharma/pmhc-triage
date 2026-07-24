"""Hermetic tests for the cBioPortal antigen client (mock transport, no network).

Mirrors the live contract verified against cBioPortal: KRAS entrez 3845, a
*_sequenced denominator, and proteinChange counting."""

import httpx
import pytest

from pmhc_triage.antigen import variant_frequency

STUDY = "paad_tcga_pan_can_atlas_2018"


def _handler(*, n_samples=179, n_g12d=49, gene_ok=True, sample_list_ok=True):
    def handler(request):
        p = request.url.path
        if "/genes/" in p:
            if not gene_ok:
                return httpx.Response(404, text="not found")
            return httpx.Response(200, json={"entrezGeneId": 3845, "hugoGeneSymbol": "KRAS"})
        if "/sample-lists/" in p:
            if not sample_list_ok:
                return httpx.Response(404, text="not found")
            return httpx.Response(200, json={"sampleIds": [f"S{i}" for i in range(n_samples)]})
        if p.endswith("/mutations/fetch"):
            recs = [{"sampleId": f"S{i}", "proteinChange": "G12D"} for i in range(n_g12d)]
            recs += [{"sampleId": "SX", "proteinChange": "G12V"}]  # noise, different change
            return httpx.Response(200, json=recs)
        return httpx.Response(404, text="unrouted")

    return httpx.Client(transport=httpx.MockTransport(handler))


def test_g12d_fraction_matches_verified_counts():
    s = variant_frequency(STUDY, "G12D", client=_handler())
    assert s.value == pytest.approx(49 / 179)
    assert "49/179" in s.provenance.method
    assert not s.is_missing
    # sample size + CI carried so the fraction isn't trusted blind
    assert s.extra["numerator"] == 49 and s.extra["denominator"] == 179
    assert s.extra["ci95_low"] < 49 / 179 < s.extra["ci95_high"]


def test_entrez_can_be_supplied_skipping_gene_lookup():
    # gene lookup returns 404, but we pass entrez explicitly -> still works
    s = variant_frequency(STUDY, "G12D", entrez=3845, client=_handler(gene_ok=False))
    assert s.value == pytest.approx(49 / 179)


def test_zero_matches_is_real_zero_not_missing():
    s = variant_frequency(STUDY, "G99Z", client=_handler(n_g12d=0))
    assert s.value == 0.0
    assert not s.is_missing
    assert any("zero samples" in w for w in s.warnings)


def test_missing_sample_list_is_missing():
    s = variant_frequency(STUDY, "G12D", client=_handler(sample_list_ok=False))
    assert s.is_missing
    assert any("HTTP 404" in w or "empty" in w for w in s.warnings)


def test_unresolvable_gene_is_missing():
    s = variant_frequency(STUDY, "G12D", gene="NOTAGENE", client=_handler(gene_ok=False))
    assert s.is_missing
    assert any("Entrez" in w for w in s.warnings)
