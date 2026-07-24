import json

import pandas as pd

from pmhc_triage.provenance import Provenance, Sourced
from pmhc_triage.report import write_provenance_json, write_results_csv
from pmhc_triage.score import score_target


def _ts():
    def S(v):
        return Sourced(v, Provenance(source="t"))

    return score_target(
        gene="KRAS", variant="G12D", disease="PDAC",
        antigen_fraction=S(0.27),
        incidence_by_population={"Europe": S(100000)},
        coverage_by_population={"Europe": S(0.48)},
    )


def test_results_csv_exposes_all_factors(tmp_path):
    p = write_results_csv([_ts()], tmp_path / "r.csv")
    df = pd.read_csv(p)
    for col in ("gene", "population", "effective_n", "incidence", "antigen_fraction", "hla_coverage"):
        assert col in df.columns
    assert df.iloc[0]["gene"] == "KRAS"
    assert df.iloc[0]["effective_n"] == round(100000 * 0.27 * 0.48, 2)


def test_provenance_json_has_target_and_factors(tmp_path):
    p = write_provenance_json([_ts()], tmp_path / "p.json")
    data = json.loads(p.read_text())
    assert data[0]["target"]["gene"] == "KRAS"
    factors = {e["factor"] for e in data[0]["provenance_log"]}
    assert "antigen_fraction" in factors and "incidence" in factors and "hla_coverage" in factors
