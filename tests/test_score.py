import pytest

from pmhc_triage.provenance import Provenance, Sourced
from pmhc_triage.score import score_target


def S(value):
    return Sourced(value, Provenance(source="test"))


def missing(reason="unavailable"):
    return Sourced(None, Provenance(source="test")).warn(reason)


def test_happy_path_multiplies_per_population():
    ts = score_target(
        gene="KRAS", variant="G12D", disease="PDAC",
        antigen_fraction=S(0.27),
        incidence_by_population={"Europe": S(500_000)},
        coverage_by_population={"Europe": S(0.48)},
    )
    ps = ts.per_population["Europe"]
    assert ps.computable
    assert ps.effective_n == pytest.approx(500_000 * 0.27 * 0.48)  # 64800.0


def test_missing_antigen_makes_all_populations_none_not_zero():
    ts = score_target(
        gene="KRAS", variant="G12D", disease="PDAC",
        antigen_fraction=missing("cBioPortal 404"),
        incidence_by_population={"Europe": S(500_000)},
        coverage_by_population={"Europe": S(0.48)},
    )
    ps = ts.per_population["Europe"]
    assert ps.effective_n is None  # NOT 0
    assert any("antigen fraction missing" in r for r in ps.reasons)


def test_join_guard_refuses_cross_population():
    # incidence for India, coverage for Europe -> neither can be scored; no India*Europe product
    ts = score_target(
        gene="KRAS", variant="G12D", disease="PDAC",
        antigen_fraction=S(0.27),
        incidence_by_population={"India": S(100_000)},
        coverage_by_population={"Europe": S(0.48)},
    )
    assert ts.per_population["India"].effective_n is None
    assert ts.per_population["Europe"].effective_n is None
    assert any("no HLA coverage" in r for r in ts.per_population["India"].reasons)
    assert any("no incidence" in r for r in ts.per_population["Europe"].reasons)
    # and the label-mismatch is surfaced at the top level
    assert any("label mismatch" in w for w in ts.warnings)


def test_partial_missing_coverage_only_affects_that_population():
    ts = score_target(
        gene="KRAS", variant="G12D", disease="PDAC",
        antigen_fraction=S(0.27),
        incidence_by_population={"Europe": S(500_000), "EastAsia": S(300_000)},
        coverage_by_population={"Europe": S(0.48), "EastAsia": missing("no A*11:01 freq")},
    )
    assert ts.per_population["Europe"].computable
    assert ts.per_population["EastAsia"].effective_n is None
    assert any("coverage missing" in r for r in ts.per_population["EastAsia"].reasons)


def test_effective_n_ci_propagated_from_antigen_and_disclosed():
    af = Sourced(0.27, Provenance(source="cbio"))
    af.extra = {"ci95_low": 0.21, "ci95_high": 0.34, "denominator": 179}
    ts = score_target(
        gene="KRAS", variant="G12D", disease="PDAC",
        antigen_fraction=af,
        incidence_by_population={"Europe": S(100000)},
        coverage_by_population={"Europe": S(0.48)},
    )
    ps = ts.per_population["Europe"]
    assert ps.effective_n_low == round(100000 * 0.21 * 0.48, 2)
    assert ps.effective_n_high == round(100000 * 0.34 * 0.48, 2)
    assert ps.effective_n_low < ps.effective_n < ps.effective_n_high
    # honesty: the CI's antigen-only nature is disclosed
    assert any("antigen-fraction sampling ONLY" in w for w in ts.warnings)


def test_provenance_log_collects_every_factor():
    ts = score_target(
        gene="KRAS", variant="G12D", disease="PDAC",
        antigen_fraction=S(0.27),
        incidence_by_population={"Europe": S(500_000)},
        coverage_by_population={"Europe": S(0.48)},
    )
    factors = {entry["factor"] for entry in ts.provenance_log}
    assert {"antigen_fraction", "incidence", "hla_coverage"} <= factors


def test_rows_flat_output_exposes_factors():
    ts = score_target(
        gene="KRAS", variant="G12D", disease="PDAC",
        antigen_fraction=S(0.27),
        incidence_by_population={"Europe": S(500_000)},
        coverage_by_population={"Europe": S(0.48)},
    )
    row = ts.rows()[0]
    assert row["gene"] == "KRAS" and row["population"] == "Europe"
    # all three factors visible in the row, never collapsed into one opaque number
    assert row["incidence"] == 500_000 and row["antigen_fraction"] == 0.27 and row["hla_coverage"] == 0.48
