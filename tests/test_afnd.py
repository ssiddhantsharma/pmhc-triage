"""Tests for the AFND parser. Uses SYNTHETIC data only -- no real AFND data is
committed to the repo (CC BY-NC; we never redistribute it)."""

import pandas as pd
import pytest

from pmhc_triage.afnd import load_afnd_frequencies
from pmhc_triage.hla import coverage_by_population


def _write(tmp_path, text, name="freqs.tsv"):
    p = tmp_path / name
    p.write_text(text)
    return p


SYNTH_TSV = (
    "allele\tpopulation\tallele_frequency\n"
    "A*02:01\tEurope\t0.28\n"
    "A*01:01\tEurope\t0.16\n"
    "A*02:01\tEastAsia\t0.10\n"
    "B*07:02\tEurope\t0.12\n"
)


def test_basic_parse(tmp_path):
    ft = load_afnd_frequencies(_write(tmp_path, SYNTH_TSV))
    assert ft.populations() == ["EastAsia", "Europe"]
    assert ft.get("Europe")["A*02:01"] == 0.28
    assert ft.get("EastAsia")["A*02:01"] == 0.10
    assert not ft.warnings


def test_hla_prefix_normalized(tmp_path):
    ft = load_afnd_frequencies(_write(tmp_path, "allele\tpopulation\tallele_frequency\nHLA-A*02:01\tX\t0.2\n"))
    assert "A*02:01" in ft.get("X")


def test_missing_column_raises(tmp_path):
    bad = "allele\tpop\tfreq\nA*02:01\tX\t0.2\n"
    with pytest.raises(ValueError, match="population"):
        load_afnd_frequencies(_write(tmp_path, bad))


def test_custom_column_names(tmp_path):
    txt = "Allele\tPopulation\tAllele Frequency\nA*02:01\tX\t0.2\n"
    ft = load_afnd_frequencies(
        _write(tmp_path, txt), freq_col="Allele Frequency"
    )
    assert ft.get("X")["A*02:01"] == 0.2


def test_percent_conversion(tmp_path):
    txt = "allele\tpopulation\tallele_frequency\nA*02:01\tX\t28\n"
    ft = load_afnd_frequencies(_write(tmp_path, txt), percent=True)
    assert ft.get("X")["A*02:01"] == pytest.approx(0.28)


def test_gt_one_without_percent_warns(tmp_path):
    txt = "allele\tpopulation\tallele_frequency\nA*02:01\tX\t28\n"
    ft = load_afnd_frequencies(_write(tmp_path, txt))
    assert any("> 1" in w for w in ft.warnings)


def test_duplicate_warns_keeps_first(tmp_path):
    txt = "allele\tpopulation\tallele_frequency\nA*02:01\tX\t0.2\nA*02:01\tX\t0.3\n"
    ft = load_afnd_frequencies(_write(tmp_path, txt))
    assert ft.get("X")["A*02:01"] == 0.2
    assert any("duplicate" in w for w in ft.warnings)


def test_integration_with_coverage(tmp_path):
    ft = load_afnd_frequencies(_write(tmp_path, SYNTH_TSV))
    cov = coverage_by_population(["A*02:01"], ft.by_population)
    # Europe freq 0.28 -> 1-(0.72)^2 = 0.4816 ; EastAsia 0.10 -> 1-(0.9)^2 = 0.19
    assert cov["Europe"].value == pytest.approx(0.4816)
    assert cov["EastAsia"].value == pytest.approx(0.19)


# --- optional sample-size column (feeds the Monte-Carlo coverage CI) --------

SYNTH_TSV_N = (
    "allele\tpopulation\tallele_frequency\tsample_size\n"
    "A*02:01\tEurope\t0.28\t1000\n"
    "A*01:01\tEurope\t0.16\t1000\n"
)


def test_sample_size_parsed_when_present(tmp_path):
    ft = load_afnd_frequencies(_write(tmp_path, SYNTH_TSV_N))
    assert ft.get_sample_sizes("Europe")["A*02:01"] == 1000
    assert not ft.warnings


def test_sample_size_absent_is_empty_not_error(tmp_path):
    # default column name asked for, but the file doesn't have it -> silently skipped
    ft = load_afnd_frequencies(_write(tmp_path, SYNTH_TSV))
    assert ft.sample_sizes == {}
    assert not any("sample" in w.lower() for w in ft.warnings)


def test_sample_size_nonpositive_warns(tmp_path):
    txt = "allele\tpopulation\tallele_frequency\tsample_size\nA*02:01\tX\t0.2\t0\n"
    ft = load_afnd_frequencies(_write(tmp_path, txt))
    assert ft.get_sample_sizes("X") == {}
    assert any("non-positive sample size" in w for w in ft.warnings)
