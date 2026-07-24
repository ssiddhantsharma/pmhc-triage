import pytest

from pmhc_triage.hla import (
    coverage_by_population,
    normalize_allele,
    parse_locus,
    population_coverage,
)


# --- allele string handling -------------------------------------------------

def test_normalize_strips_hla_prefix_and_whitespace():
    assert normalize_allele(" HLA-A*02:01 ") == "A*02:01"
    assert normalize_allele("A*02:01") == "A*02:01"


def test_parse_locus():
    assert parse_locus("A*02:01") == "A"
    assert parse_locus("HLA-B*07:02") == "B"
    assert parse_locus("C*01:02") == "C"


def test_parse_locus_rejects_malformed():
    with pytest.raises(ValueError):
        parse_locus("A0201")


# --- coverage maths (analytically-derived expected values) ------------------
# Method: coverage = 1 - prod_L (1 - p_L)^2, p_L = sum of covering-allele freqs at L.

def test_single_allele_single_locus():
    # 1 - (1-0.25)^2 = 0.4375
    s = population_coverage(["A*02:01"], {"A*02:01": 0.25}, "TestPop")
    assert s.value == pytest.approx(0.4375)
    assert not s.is_missing and not s.warnings


def test_two_alleles_same_locus_sum_then_square():
    # p_A = 0.40 -> 1 - 0.6^2 = 0.64
    s = population_coverage(
        ["A*02:01", "A*01:01"], {"A*02:01": 0.25, "A*01:01": 0.15}, "TestPop"
    )
    assert s.value == pytest.approx(0.64)


def test_two_loci_combine_independently():
    # 1 - (0.75^2 * 0.90^2) = 1 - 0.455625 = 0.544375
    s = population_coverage(
        ["A*02:01", "B*07:02"], {"A*02:01": 0.25, "B*07:02": 0.10}, "TestPop"
    )
    assert s.value == pytest.approx(0.544375)


def test_hla_prefix_matches_bare_allele():
    s = population_coverage(["HLA-A*02:01"], {"A*02:01": 0.25}, "TestPop")
    assert s.value == pytest.approx(0.4375)


# --- the surfacing discipline: missing freq is NOT treated as zero ----------

def test_missing_allele_excluded_and_surfaced_not_zeroed():
    # A*99:01 has no frequency -> excluded from maths (NOT 0), coverage from A*02:01 only.
    s = population_coverage(
        ["A*02:01", "A*99:01"], {"A*02:01": 0.25}, "TestPop"
    )
    assert s.value == pytest.approx(0.4375)  # same as if A*99:01 weren't listed
    assert any("A*99:01" in w for w in s.warnings)
    assert any("NOT treated as 0" in w for w in s.warnings)


def test_all_alleles_missing_returns_missing():
    s = population_coverage(["A*99:01"], {"A*02:01": 0.25}, "TestPop")
    assert s.is_missing and s.value is None
    assert s.warnings


def test_no_covering_alleles_is_missing():
    s = population_coverage([], {"A*02:01": 0.25}, "TestPop")
    assert s.is_missing


def test_frequency_sum_over_one_is_clamped_and_warned():
    s = population_coverage(
        ["A*02:01", "A*01:01"], {"A*02:01": 0.7, "A*01:01": 0.5}, "TestPop"
    )
    assert s.value == pytest.approx(1.0)  # clamped p_A=1.0 -> full coverage
    assert any("clamped" in w for w in s.warnings)


# --- provenance travels with the number -------------------------------------

def test_provenance_recorded():
    s = population_coverage(["A*02:01"], {"A*02:01": 0.25}, "Europe", url="http://afnd")
    d = s.to_dict()
    assert "AFND" in d["provenance"]["source"]
    assert d["provenance"]["url"] == "http://afnd"
    assert "Bui" in d["provenance"]["method"]
    assert d["provenance"]["query_date"]  # stamped


# --- multi-population helper ------------------------------------------------

def test_coverage_by_population():
    out = coverage_by_population(
        ["A*02:01"],
        {"Europe": {"A*02:01": 0.28}, "EastAsia": {"A*02:01": 0.10}},
    )
    assert set(out) == {"Europe", "EastAsia"}
    assert out["Europe"].value > out["EastAsia"].value  # higher freq -> higher coverage
