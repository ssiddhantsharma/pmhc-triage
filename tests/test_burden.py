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


# --- curated bundle (exact GCO Cancer Today 2024 figures, cited) -----------

def test_bundled_pancreatic_exact_figure():
    s = load_bundled("pancreatic cancer")
    assert s.value == 531318.0  # exact GCO Cancer Today 2024 figure
    assert "2024" in s.provenance.source


def test_bundled_case_insensitive():
    assert load_bundled("Lung Cancer").value == 2637005.0


def test_bundled_surfaces_world_and_note_caveats():
    s = load_bundled("pancreatic cancer")
    assert any("World-level" in w for w in s.warnings)      # can't join per-region silently
    assert any("bundle note" in w for w in s.warnings)      # PDAC~90% caveat visible at runtime


def test_bundled_world_level_caveat_on_every_row():
    # bundle figures are now exact (2024), but all are World-level -> that join
    # limitation must still surface on any row (can't join per-region silently).
    s = load_bundled("colorectal carcinoma")
    assert any("World-level" in w for w in s.warnings)


def test_bundled_unknown_is_missing_and_lists_options():
    s = load_bundled("made up disease")
    assert s.is_missing
    assert any("available" in w for w in s.warnings)


def test_bundled_diseases_listed():
    diseases = bundled_diseases()
    assert "pancreatic cancer" in diseases
    assert len(diseases) >= 5


# --- per-region starter bundle (country-proxy figures) ---------------------

def test_bundled_per_region_exact_figures():
    # GCO 2024 country exports folded into the bundle (Germany/China/India)
    assert load_bundled("pancreatic cancer", population="Europe").value == 22941.0
    assert load_bundled("pancreatic cancer", population="EastAsia").value == 122597.0
    assert load_bundled("pancreatic cancer", population="SouthAsia").value == 20477.0


def test_bundled_per_region_surfaces_country_proxy_note():
    s = load_bundled("pancreatic cancer", population="Europe")
    # the Germany-as-Europe proxy must be surfaced at runtime, not hidden
    assert any("Germany" in w and "proxy" in w for w in s.warnings)
    # and it is NOT a World row, so no World-only caveat here
    assert not any("World-level" in w for w in s.warnings)


def test_bundled_per_region_unknown_population_missing():
    s = load_bundled("pancreatic cancer", population="Antarctica")
    assert s.is_missing
