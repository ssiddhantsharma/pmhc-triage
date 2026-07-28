"""Tests for the threshold-sensitivity sweep. Pure -- no MHCflurry needed."""

import pandas as pd
import pytest

from pmhc_triage.presentation import sweep_presenting_alleles
from pmhc_triage.sensitivity import threshold_sensitivity


def _scored():
    # per-allele presentation percentile (lower = stronger)
    return pd.DataFrame([
        {"sample_name": "A*02:01", "presentation_percentile": 0.4},   # strong
        {"sample_name": "A*01:01", "presentation_percentile": 1.8},   # medium
        {"sample_name": "B*07:02", "presentation_percentile": 4.0},   # weak
    ])


def test_sweep_admits_more_alleles_as_threshold_loosens():
    swept = sweep_presenting_alleles(_scored(), [0.5, 2.0, 5.0])
    assert swept[0.5] == ["A*02:01"]
    assert swept[2.0] == ["A*01:01", "A*02:01"]
    assert swept[5.0] == ["A*01:01", "A*02:01", "B*07:02"]


def test_coverage_monotonic_in_threshold():
    freqs = {"Europe": {"A*02:01": 0.28, "A*01:01": 0.16, "B*07:02": 0.12}}
    rows = threshold_sensitivity(_scored(), [0.5, 2.0, 5.0], freqs)
    covs = [r.coverage_by_population["Europe"] for r in rows]
    assert covs[0] < covs[1] < covs[2]  # looser cutoff -> more alleles -> more coverage


def test_effective_n_included_when_incidence_and_antigen_given():
    freqs = {"Europe": {"A*02:01": 0.28}}
    rows = threshold_sensitivity(
        _scored(), [0.5], freqs,
        incidence_by_population={"Europe": 100000},
        antigen_fraction=0.25,
    )
    # coverage(A*02:01)=0.4816 -> 100000*0.25*0.4816 = 12040
    assert rows[0].effective_n_by_population["Europe"] == pytest.approx(12040, abs=1)


def test_missing_frequency_gives_none_coverage_not_zero():
    freqs = {"Europe": {"A*01:01": 0.16}}  # A*02:01 (the strong one) absent
    rows = threshold_sensitivity(_scored(), [0.5], freqs)
    assert rows[0].coverage_by_population["Europe"] is None  # surfaced, not 0
    assert rows[0].effective_n_by_population["Europe"] is None


def test_row_to_dict_shape():
    rows = threshold_sensitivity(_scored(), [2.0], {"Europe": {"A*02:01": 0.28}})
    d = rows[0].to_dict()
    assert d["threshold"] == 2.0 and d["n_alleles"] == 2 and "coverage_by_population" in d
