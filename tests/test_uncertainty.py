import pytest

from pmhc_triage.uncertainty import wilson_ci


def test_wilson_brackets_point_estimate():
    lo, hi = wilson_ci(49, 179)
    assert lo < 49 / 179 < hi
    assert 0.0 <= lo and hi <= 1.0


def test_small_n_wider_than_large_n_same_proportion():
    lo_s, hi_s = wilson_ci(5, 18)      # ~0.28, n=18
    lo_l, hi_l = wilson_ci(500, 1800)  # ~0.28, n=1800
    assert (hi_s - lo_s) > (hi_l - lo_l)  # small n => wider interval


def test_n_zero_is_total_ignorance():
    assert wilson_ci(0, 0) == (0.0, 1.0)


def test_known_value():
    # 49/179: Wilson 95% CI is approximately [0.212, 0.345]
    lo, hi = wilson_ci(49, 179)
    assert lo == pytest.approx(0.212, abs=0.01)
    assert hi == pytest.approx(0.345, abs=0.01)


def test_wilson_out_of_range_successes_does_not_crash():
    # defensive: successes > n must clamp, never raise a math domain error
    lo, hi = wilson_ci(200, 179)
    assert 0.0 <= lo <= 1.0 and 0.0 <= hi <= 1.0
    assert wilson_ci(-5, 100)[0] >= 0.0
