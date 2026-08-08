"""Rendering is presentation-only; these assert it runs and orders correctly.
Fully offline (constructs TargetScores directly, renders to a StringIO console)."""

import io

from rich.console import Console

from pmhc_triage.provenance import Provenance, Sourced
from pmhc_triage.render import render, render_rank
from pmhc_triage.score import score_target


def S(v):
    return Sourced(v, Provenance(source="t"))


def _score(gene, incidence, coverage=0.5, antigen=0.27):
    af = Sourced(antigen, Provenance(source="cbio"))
    af.extra = {"numerator": 49, "denominator": 179, "ci95_low": 0.21, "ci95_high": 0.34}
    cov = S(coverage) if coverage is not None else Sourced(None, Provenance(source="t")).warn("no freq")
    return score_target(
        gene=gene, variant="G12D", disease="PDAC",
        antigen_fraction=af,
        incidence_by_population={"Europe": S(incidence)},
        coverage_by_population={"Europe": cov},
    )


def _out(fn, scores):
    buf = io.StringIO()
    fn(scores, Console(file=buf, force_terminal=False, width=160))
    return buf.getvalue()


def test_render_shows_target_population_and_effective_n():
    out = _out(render, [_score("KRAS", 100000)])
    assert "KRAS" in out and "Europe" in out
    assert "13,500" in out  # 100000 * 0.27 * 0.5


def test_render_rank_orders_by_effective_n_desc():
    out = _out(render_rank, [_score("TP53", 10000), _score("KRAS", 100000)])
    assert out.index("KRAS") < out.index("TP53")  # bigger effective-N first


def test_render_rank_missing_sorts_last_not_dropped():
    good = _score("KRAS", 100000)
    bad = _score("MISSING", 100000, coverage=None)  # coverage missing -> effective_N None
    out = _out(render_rank, [bad, good])
    assert "MISSING" in out                      # kept, not dropped
    assert out.index("KRAS") < out.index("MISSING")  # computable ranks above missing


def test_render_does_not_crash_on_all_missing():
    _out(render, [_score("X", 100000, coverage=None)])  # should render without error
