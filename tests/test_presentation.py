import importlib.util

import pandas as pd
import pytest

_HAS_MHCFLURRY = importlib.util.find_spec("mhcflurry") is not None

from pmhc_triage.presentation import (
    manual_presenting_alleles,
    predict_presenting_alleles,
    select_presenting,
)


# --- manual path (the default) ---------------------------------------------

def test_manual_normalizes_and_dedupes():
    s = manual_presenting_alleles(["HLA-A*02:01", "A*02:01", "B*07:02"])
    assert s.value == ["A*02:01", "B*07:02"]  # dedup + normalized, order preserved
    assert not s.is_missing


def test_manual_excludes_malformed_and_surfaces():
    s = manual_presenting_alleles(["A*02:01", "GARBAGE"])
    assert s.value == ["A*02:01"]
    assert any("GARBAGE" in w for w in s.warnings)


def test_manual_all_invalid_is_missing():
    s = manual_presenting_alleles(["nope", "also-nope"])
    assert s.is_missing


# --- pure selection core (tested without MHCflurry) ------------------------

def test_select_presenting_by_percentile():
    df = pd.DataFrame([
        {"allele": "A*02:01", "percentile": 0.5},   # strong -> present
        {"allele": "A*02:01", "percentile": 8.0},   # (min over peptides is 0.5)
        {"allele": "A*01:01", "percentile": 5.0},   # never below threshold
        {"allele": "B*07:02", "percentile": 1.9},   # just under threshold -> present
    ])
    assert select_presenting(df, threshold=2.0) == ["A*02:01", "B*07:02"]


def test_select_presenting_higher_is_better_mode():
    df = pd.DataFrame([
        {"allele": "A*02:01", "score": 0.9},
        {"allele": "A*01:01", "score": 0.3},
    ])
    assert select_presenting(df, score_col="score", threshold=0.7, lower_is_better=False) == ["A*02:01"]


# --- optional MHCflurry path degrades gracefully ---------------------------

def test_predict_without_mhcflurry_is_missing_not_crash():
    # No _predictor injected and (in the test env) no mhcflurry installed -> surfaced missing.
    import importlib.util

    if importlib.util.find_spec("mhcflurry") is not None:
        pytest.skip("mhcflurry installed; covered by the gated live test")
    s = predict_presenting_alleles(["VVVGADGVGK"], ["A*11:01", "A*03:01"])
    assert s.is_missing
    assert any("mhcflurry not installed" in w for w in s.warnings)


def test_predict_with_injected_predictor():
    # matches the real MHCflurry shape: per-allele 'sample_name' + 'presentation_percentile'
    class FakePredictor:
        def predict(self, peptides, alleles):
            return pd.DataFrame([
                {"sample_name": "A*11:01", "presentation_percentile": 0.3},
                {"sample_name": "A*03:01", "presentation_percentile": 9.0},
            ])

    s = predict_presenting_alleles(["VVVGADGVGK"], ["A*11:01", "A*03:01"], _predictor=FakePredictor())
    assert s.value == ["A*11:01"]  # only the strong binder


@pytest.mark.skipif(not _HAS_MHCFLURRY, reason="mhcflurry (optional extra) not installed")
def test_live_mhcflurry_kras_g12d_restriction():
    # Real end-to-end: KRAS G12D is A*11:01 / A*03:01 restricted, NOT A*02:01.
    peps = ["VVVGADGVGK", "VVGADGVGK", "GADGVGKSA"]
    s = predict_presenting_alleles(peps, ["A*11:01", "A*03:01", "A*02:01"], threshold_percentile=2.0)
    assert not s.is_missing, s.warnings
    assert "A*11:01" in s.value and "A*03:01" in s.value
    assert "A*02:01" not in s.value  # weak binder -> real biology reproduced
