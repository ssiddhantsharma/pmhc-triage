"""Hermetic tests for the orchestrator + preflight.

A single mock transport routes across all three hosts (cBioPortal, Open Targets,
UniProt) by request host, so run_target/preflight can be exercised with no network."""

import httpx
import pytest

from pmhc_triage.pipeline import TargetSpec, preflight, run_target

FREQS = (
    "allele\tpopulation\tallele_frequency\n"
    "A*11:01\tEurope\t0.05\n"
    "A*03:01\tEurope\t0.13\n"
    "A*11:01\tEastAsia\t0.30\n"
    "A*03:01\tEastAsia\t0.02\n"
)
BURDEN = "disease,population,incidence\nPDAC,Europe,100000\nPDAC,EastAsia,150000\n"

# minimal real KRAS N-terminus so the UniProt WT check can pass (residue 12 = G)
KRAS_FASTA = ">sp|P01116|RASK_HUMAN\nMTEYKLVVVGAGGVGKSALTIQLIQ\n"


def _routing_client(*, study_ok=True, g12d=49, n=179, cancer_type="Pancreatic Adenocarcinoma"):
    def handler(request):
        host, path = request.url.host, request.url.path
        if host == "rest.uniprot.org":
            return httpx.Response(200, text=KRAS_FASTA)
        if host == "api.platform.opentargets.org":
            body = request.content.decode()
            if "search" in body:
                return httpx.Response(200, json={"data": {"search": {"hits": [
                    {"id": "ENSG00000133703", "name": "KRAS", "entity": "target"}]}}})
            return httpx.Response(200, json={"data": {"target": {
                "approvedSymbol": "KRAS",
                "tractability": [{"label": "Approved Drug", "modality": "SM", "value": True}]}}})
        if host == "www.cbioportal.org":
            if "/genes/" in path:
                return httpx.Response(200, json={"entrezGeneId": 3845, "hugoGeneSymbol": "KRAS"})
            if "/sample-lists/" in path:
                if not study_ok:
                    return httpx.Response(404, text="no study")
                return httpx.Response(200, json={"sampleIds": [f"S{i}" for i in range(n)]})
            if "/studies/" in path:
                return httpx.Response(200, json={"cancerType": {"name": cancer_type}})
            if path.endswith("/mutations/fetch"):
                return httpx.Response(200, json=[{"sampleId": f"S{i}", "proteinChange": "G12D"} for i in range(g12d)])
        return httpx.Response(404, text=f"unrouted {host}{path}")

    return httpx.Client(transport=httpx.MockTransport(handler))


def _spec(tmp_path, **over):
    fp = tmp_path / "freqs.tsv"; fp.write_text(FREQS)
    bp = tmp_path / "burden.csv"; bp.write_text(BURDEN)
    d = dict(gene="KRAS", variant="G12D", disease="PDAC",
             study="paad_tcga_pan_can_atlas_2018", alleles=["A*11:01", "A*03:01"],
             populations=["Europe", "EastAsia"], freqs_path=str(fp), burden_path=str(bp))
    d.update(over)
    return TargetSpec(**d)


def test_run_target_composes_end_to_end(tmp_path):
    ts = run_target(_spec(tmp_path), client=_routing_client())
    eu = ts.per_population["Europe"]
    # antigen 49/179; Europe coverage 1-(1-0.18)^2=0.3276; incidence 100000
    assert eu.effective_n == pytest.approx(100000 * (49 / 179) * 0.3276, rel=1e-4)
    assert ts.per_population["EastAsia"].computable
    # tractability context captured in the provenance log
    assert any(e["factor"] == "tractability_context" for e in ts.provenance_log)


def test_preflight_clean_passes(tmp_path):
    spec = _spec(tmp_path, uniprot="P01116")
    assert preflight(spec, client=_routing_client()) == []


def test_preflight_catches_bad_study(tmp_path):
    issues = preflight(_spec(tmp_path), client=_routing_client(study_ok=False))
    assert any("study" in i for i in issues)


def test_preflight_catches_missing_population_in_freqs(tmp_path):
    spec = _spec(tmp_path, populations=["Europe", "Africa"])  # Africa not in freqs
    issues = preflight(spec, client=_routing_client())
    assert any("Africa" in i for i in issues)


def test_preflight_catches_study_disease_mismatch(tmp_path):
    # study reports a breast cancer type, but --disease is PDAC -> must flag
    spec = _spec(tmp_path)
    issues = preflight(spec, client=_routing_client(cancer_type="Invasive Breast Carcinoma"))
    assert any("may not correspond" in i for i in issues)


def test_preflight_ok_when_study_matches_disease(tmp_path):
    spec = _spec(tmp_path, uniprot="P01116")
    issues = preflight(spec, client=_routing_client(cancer_type="Pancreatic Adenocarcinoma"))
    assert not any("correspond" in i for i in issues)


def test_preflight_catches_wt_mismatch(tmp_path):
    spec = _spec(tmp_path, variant="A12D", uniprot="P01116")  # seq has G at 12, not A
    issues = preflight(spec, client=_routing_client())
    assert any("WT mismatch" in i for i in issues)


def test_predict_alleles_mode_uses_predicted_not_manual(tmp_path, monkeypatch):
    # predict_alleles=True: pipeline should fetch seq -> peptides -> predict, and use THAT
    import pmhc_triage.pipeline as pl
    from pmhc_triage.provenance import Provenance, Sourced

    captured = {}

    def fake_predict(peptides, panel, *, threshold_percentile=2.0):
        captured["panel"] = panel
        return Sourced(["A*11:01"], Provenance(source="MHCflurry (fake)"))

    monkeypatch.setattr(pl, "predict_presenting_alleles", fake_predict)
    spec = _spec(tmp_path, alleles=[], uniprot="P01116", predict_alleles=True)
    ts = run_target(spec, client=_routing_client())
    # coverage computed from the PREDICTED allele A*11:01 (present in both pops)
    assert ts.per_population["Europe"].computable
    # the panel handed to the predictor came from the freq file
    assert "A*11:01" in captured["panel"] and "A*03:01" in captured["panel"]
    # provenance records the predicted allele source
    assert any(e["factor"] == "presenting_alleles" and "MHCflurry" in e["provenance"]["source"]
               for e in ts.provenance_log)


def test_preflight_suggests_population_alias(tmp_path):
    # freqs file uses "Europe"; request "European" -> preflight should suggest the alias
    spec = _spec(tmp_path, populations=["European"])
    issues = preflight(spec, client=_routing_client())
    assert any("did you mean 'Europe'" in i for i in issues)
