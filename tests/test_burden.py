import pytest

from pmhc_triage.burden import bundled_diseases, load_bundled, load_burden_table, manual_incidence


# --- manual + CSV (not bundled) --------------------------------------------

def test_manual_incidence_requires_source():
    s = manual_incidence("PDAC", "Europe", 100_000, source="my estimate")
    assert s.value == 100_000.0
    assert s.provenance.source == "my estimate"


def test_manual_rejects_negative():
    s = manual_incidence("PDAC", "Europe", -5, source="x")
    assert s.is_missing


def test_load_burden_table(tmp_path):
    p = tmp_path / "b.csv"
    p.write_text("disease,population,incidence\nPDAC,Europe,100000\nPDAC,EastAsia,150000\nCRC,Europe,50000\n")
    out = load_burden_table(p, "PDAC")
    assert set(out) == {"Europe", "EastAsia"}  # only PDAC rows
    assert out["Europe"].value == 100_000.0


# --- curated bundle (verified GLOBOCAN 2022 figures, cited) -----------------

def test_bundled_pancreatic_exact_figure():
    s = load_bundled("pancreatic cancer")
    assert s.value == 510992.0  # exact GLOBOCAN 2022 figure
    assert "GLOBOCAN 2022" in s.provenance.source


def test_bundled_case_insensitive():
    assert load_bundled("Lung Cancer").value == 2480000.0


def test_bundled_unknown_is_missing_and_lists_options():
    s = load_bundled("made up disease")
    assert s.is_missing
    assert any("available" in w for w in s.warnings)


def test_bundled_diseases_listed():
    diseases = bundled_diseases()
    assert "pancreatic cancer" in diseases
    assert len(diseases) >= 5
