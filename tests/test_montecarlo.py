"""Tests for Monte-Carlo uncertainty propagation.

These verify the *propagation logic* (draws centre on the point estimate, intervals
widen with smaller samples, missing data is surfaced not faked, runs are
reproducible) -- not that the interval matches a real addressable population.
"""

import numpy as np
import pytest

from pmhc_triage.montecarlo import (
    MCResult,
    effective_n_interval,
    sample_antigen_fraction,
    sample_coverage,
)


def _rng(seed=0):
    return np.random.default_rng(seed)


# --- antigen fraction ------------------------------------------------------

def test_antigen_draws_centre_on_point_estimate():
    draws = sample_antigen_fraction(50, 200, 50000, _rng())
    assert np.median(draws) == pytest.approx(0.25, abs=0.01)


def test_antigen_smaller_n_widens_interval():
    small = sample_antigen_fraction(5, 20, 50000, _rng())
    large = sample_antigen_fraction(500, 2000, 50000, _rng())  # same 0.25, 100x n
    width = lambda d: np.percentile(d, 97.5) - np.percentile(d, 2.5)
    assert width(small) > width(large)


def test_antigen_zero_denominator_raises():
    with pytest.raises(ValueError):
        sample_antigen_fraction(0, 0, 100, _rng())


# --- coverage --------------------------------------------------------------

def test_coverage_draws_centre_on_bui_point_estimate():
    # freq 0.28, big sample -> coverage ~ 1-(1-0.28)^2 = 0.4816
    draws, warn = sample_coverage(
        ["A*02:01"], {"A*02:01": 0.28}, {"A*02:01": 100000}, 50000, _rng()
    )
    assert np.median(draws) == pytest.approx(0.4816, abs=0.01)
    assert not warn


def test_coverage_smaller_sample_widens_interval():
    w = lambda n: (
        lambda d: float(np.percentile(d, 97.5) - np.percentile(d, 2.5))
    )(sample_coverage(["A*02:01"], {"A*02:01": 0.28}, {"A*02:01": n}, 50000, _rng())[0])
    assert w(50) > w(5000)


def test_coverage_no_sample_size_holds_fixed_and_surfaces():
    draws, warn = sample_coverage(["A*02:01"], {"A*02:01": 0.28}, {}, 1000, _rng())
    assert np.allclose(draws, 0.4816)          # no sampling error -> constant
    assert any("no AFND sample size" in w for w in warn)


def test_coverage_missing_frequency_excluded_not_zero():
    draws, warn = sample_coverage(
        ["A*02:01", "B*07:02"], {"A*02:01": 0.28}, {"A*02:01": 1000}, 1000, _rng()
    )
    assert any("B*07:02" in w and "NOT treated as 0" in w for w in warn)
    assert draws.size == 1000  # still computed from the one known allele


def test_coverage_no_known_allele_is_empty_with_warning():
    draws, warn = sample_coverage(["B*07:02"], {"A*02:01": 0.28}, {}, 1000, _rng())
    assert draws.size == 0
    assert any("NOT 0" in w for w in warn)


def test_coverage_two_loci_combine_diploid():
    # A locus p=0.28, B locus p=0.12; not-covered = (0.72^2)(0.88^2); coverage = 1 - that
    draws, _ = sample_coverage(
        ["A*02:01", "B*07:02"],
        {"A*02:01": 0.28, "B*07:02": 0.12},
        {"A*02:01": 100000, "B*07:02": 100000},
        50000,
        _rng(),
    )
    expected = 1 - (0.72**2) * (0.88**2)
    assert np.median(draws) == pytest.approx(expected, abs=0.01)


# --- effective-N interval --------------------------------------------------

def test_effective_n_median_matches_product():
    r = effective_n_interval(
        incidence=100000,
        antigen_numerator=50,
        antigen_denominator=200,   # 0.25
        covering_alleles=["A*02:01"],
        allele_freqs={"A*02:01": 0.28},
        sample_sizes={"A*02:01": 100000},
        n_draws=50000,
    )
    # 100000 * 0.25 * 0.4816 = 12040
    assert r.median == pytest.approx(12040, rel=0.03)
    assert r.ci95_low < r.median < r.ci95_high


def test_effective_n_incidence_fixed_by_default_and_labeled():
    r = effective_n_interval(
        incidence=100000, antigen_numerator=50, antigen_denominator=200,
        covering_alleles=["A*02:01"], allele_freqs={"A*02:01": 0.28},
        sample_sizes={"A*02:01": 1000}, n_draws=2000,
    )
    assert any("incidence held FIXED" in c for c in r.caveats)


def test_effective_n_user_incidence_sd_widens_and_flags():
    common = dict(
        incidence=100000, antigen_numerator=50, antigen_denominator=200,
        covering_alleles=["A*02:01"], allele_freqs={"A*02:01": 0.28},
        sample_sizes={"A*02:01": 1000}, n_draws=50000,
    )
    fixed = effective_n_interval(**common)
    noisy = effective_n_interval(**common, incidence_rel_sd=0.2)
    assert (noisy.ci95_high - noisy.ci95_low) > (fixed.ci95_high - fixed.ci95_low)
    assert any("USER-supplied" in c for c in noisy.caveats)


def test_effective_n_non_computable_when_coverage_unsampleable():
    r = effective_n_interval(
        incidence=100000, antigen_numerator=50, antigen_denominator=200,
        covering_alleles=["B*07:02"], allele_freqs={"A*02:01": 0.28},
        sample_sizes={}, n_draws=1000,
    )
    assert not r.computable
    assert r.median is None
    assert r.to_dict()["effective_n_mc_ci95"] is None


def test_reproducible_same_seed():
    kw = dict(
        incidence=100000, antigen_numerator=50, antigen_denominator=200,
        covering_alleles=["A*02:01"], allele_freqs={"A*02:01": 0.28},
        sample_sizes={"A*02:01": 1000}, n_draws=5000,
    )
    assert effective_n_interval(**kw, seed=7).median == effective_n_interval(**kw, seed=7).median


def test_different_seed_differs_slightly():
    kw = dict(
        incidence=100000, antigen_numerator=5, antigen_denominator=20,
        covering_alleles=["A*02:01"], allele_freqs={"A*02:01": 0.28},
        sample_sizes={"A*02:01": 50}, n_draws=3000,
    )
    a = effective_n_interval(**kw, seed=1).median
    b = effective_n_interval(**kw, seed=2).median
    assert a != b  # tiny cohort -> visible MC noise between seeds
