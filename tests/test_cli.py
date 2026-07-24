"""Hermetic CLI tests: monkeypatch the network-touching entry points."""

import json

from pmhc_triage import cli
from pmhc_triage.provenance import Provenance, Sourced
from pmhc_triage.score import score_target


def _fake_ts():
    def S(v):
        return Sourced(v, Provenance(source="t"))

    return score_target(
        gene="KRAS", variant="G12D", disease="PDAC",
        antigen_fraction=S(0.27),
        incidence_by_population={"Europe": S(100000)},
        coverage_by_population={"Europe": S(0.48)},
    )


_SCORE_ARGS = [
    "--gene", "KRAS", "--variant", "G12D", "--disease", "PDAC",
    "--study", "s", "--alleles", "A*11:01", "--populations", "Europe",
]


def test_score_writes_csv_and_provenance(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "run_target", lambda spec, client=None: _fake_ts())
    out = tmp_path / "r.csv"
    rc = cli.main(["score", *_SCORE_ARGS, "--out", str(out)])
    assert rc == 0
    assert out.exists() and (tmp_path / "r.provenance.json").exists()


def test_validate_exit_codes(monkeypatch):
    monkeypatch.setattr(cli, "preflight", lambda spec, client=None: [])
    assert cli.main(["validate", *_SCORE_ARGS]) == 0
    monkeypatch.setattr(cli, "preflight", lambda spec, client=None: ["study 's' not found"])
    assert cli.main(["validate", *_SCORE_ARGS]) == 1


def test_discover_writes_csv(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "resolve_disease",
                        lambda name, client=None: Sourced("MONDO_x", Provenance(source="t")))
    monkeypatch.setattr(cli, "associated_targets",
                        lambda did, top=20, client=None: Sourced(
                            [{"symbol": "KRAS", "ensembl_id": "ENSG1", "association_score": 0.8}],
                            Provenance(source="t")))
    out = tmp_path / "d.csv"
    rc = cli.main(["discover", "--disease", "x", "--out", str(out)])
    assert rc == 0
    assert "KRAS" in out.read_text()


def test_no_command_prints_help():
    assert cli.main([]) == 0
