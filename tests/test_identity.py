from pmhc_triage.identity import (
    align_populations,
    canonical_disease,
    canonical_population,
    suggest_match,
)


def test_canonical_population_synonyms():
    assert canonical_population("European") == "Europe"
    assert canonical_population("  europe ") == "Europe"
    assert canonical_population("East Asian") == "EastAsia"
    assert canonical_population("global") == "World"


def test_canonical_population_unknown_passthrough_trimmed():
    assert canonical_population("  Finland ") == "Finland"  # unknown -> trimmed, not merged


def test_conservative_no_country_to_region_merge():
    # India is a country; we must NOT silently merge it into a region
    assert canonical_population("India") == "India"
    assert canonical_population("Indian") == "Indian"


def test_canonical_disease():
    assert canonical_disease("PDAC") == "pancreatic ductal adenocarcinoma"
    assert canonical_disease("Melanoma") == "Melanoma"  # unknown passthrough


def test_suggest_match():
    assert suggest_match("Europe", ["European", "EastAsia"]) == "European"
    assert suggest_match("Europe", ["Europe"]) is None  # identical, nothing to suggest
    assert suggest_match("Mars", ["Europe"]) is None


def test_align_populations():
    out = align_populations(["Europe", "EastAsia"], ["European", "Africa"])
    assert out["matched"] == ["Europe"]
    assert "Africa" in out["only_b"]
    assert "EastAsia" in out["only_a"]
