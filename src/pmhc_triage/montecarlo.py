"""Monte-Carlo uncertainty propagation for effective addressable-N.

The point estimate ``effective_N = incidence x antigen_fraction x hla_coverage`` is
only as trustworthy as its softest factor, and the shipped CI historically covered
*only* the antigen-fraction sampling error. This module propagates the two factors
whose sampling uncertainty is genuinely sourceable into a real interval on
``effective_N``:

- **antigen fraction** -- a binomial proportion ``numer/denom`` from a cBioPortal
  cohort. Sampled from a Jeffreys ``Beta(numer + 1/2, denom - numer + 1/2)``.
- **HLA coverage** -- a nonlinear function of per-allele frequencies, each estimated
  from a finite AFND sample. Each covering allele's frequency is sampled from its own
  ``Beta`` posterior using the AFND sample size (2 x individuals = chromosomes typed),
  then coverage is recomputed per draw with the exact Bui-2006 diploid formula.

**Incidence is held fixed by default and said so.** GLOBOCAN publishes point
estimates without a CI, so inventing one would violate the tool's contract. The
caller may pass ``incidence_rel_sd`` to inject a *user-owned* relative uncertainty
(e.g. ``0.15`` for +/-15%); otherwise the returned interval explicitly reflects
antigen + coverage sampling only, and the ``method``/``caveats`` say so.

Reproducibility: everything is driven by a seeded ``numpy`` generator (default seed
0), so an MC interval replays byte-identically -- consistent with the rest of the
package treating numbers as sourced and reproducible, not stochastic surprises.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .hla import normalize_allele, parse_locus

# Jeffreys prior pseudo-counts -- a proper, well-behaved posterior at 0/1 and small n.
_JEFFREYS = 0.5


@dataclass
class MCResult:
    """A Monte-Carlo interval plus exactly what uncertainty went into it."""

    median: float | None
    ci95_low: float | None
    ci95_high: float | None
    n_draws: int
    method: str
    caveats: list[str] = field(default_factory=list)
    # component draws' summaries, for auditing which factor dominates the spread
    components: dict = field(default_factory=dict)

    @property
    def computable(self) -> bool:
        return self.median is not None

    def to_dict(self) -> dict:
        return {
            "effective_n_mc_median": self.median,
            "effective_n_mc_ci95": (
                [self.ci95_low, self.ci95_high] if self.median is not None else None
            ),
            "mc_n_draws": self.n_draws,
            "mc_method": self.method,
            "mc_caveats": list(self.caveats),
            "mc_components": dict(self.components),
        }


def sample_antigen_fraction(numerator: int, denominator: int, n_draws: int, rng) -> np.ndarray:
    """Draw antigen-positive fractions from the Jeffreys Beta posterior of numer/denom."""
    if denominator <= 0:
        raise ValueError("denominator must be > 0 to sample an antigen fraction")
    a = numerator + _JEFFREYS
    b = (denominator - numerator) + _JEFFREYS
    return rng.beta(a, b, size=n_draws)


def sample_coverage(
    covering_alleles,
    allele_freqs,
    sample_sizes,
    n_draws: int,
    rng,
) -> tuple[np.ndarray, list[str]]:
    """Draw HLA coverage values, propagating each allele's finite-sample frequency error.

    Parameters
    ----------
    covering_alleles
        Alleles that present the epitope (any format; normalized internally).
    allele_freqs
        ``{allele: frequency}`` for one population (per-chromosome fractions).
    sample_sizes
        ``{allele: n_individuals}`` from AFND. May be partial: an allele with no
        sample size is held at its point frequency (no sampling error) and surfaced.
    n_draws, rng
        Draw count and a seeded ``numpy`` Generator.

    Returns
    -------
    (draws, warnings)
        ``draws`` is an ``(n_draws,)`` array of coverage in ``[0, 1]``. Empty (size 0)
        with a warning if no covering allele has a known frequency -- mirroring the
        point-estimate module, a missing frequency is never treated as 0.
    """
    norm_freqs = {normalize_allele(k): float(v) for k, v in allele_freqs.items()}
    norm_n = {normalize_allele(k): v for k, v in (sample_sizes or {}).items()}
    covering = [normalize_allele(a) for a in covering_alleles]
    warnings: list[str] = []

    # Per-locus list of (n_draws,) frequency-sample arrays for covering alleles present.
    per_locus: dict[str, list[np.ndarray]] = {}
    fixed_no_n: list[str] = []
    missing: list[str] = []
    for allele in covering:
        if allele not in norm_freqs:
            missing.append(allele)
            continue
        locus = parse_locus(allele)
        p_hat = norm_freqs[allele]
        n_ind = norm_n.get(allele)
        if n_ind and n_ind > 0:
            two_n = 2 * float(n_ind)  # chromosomes
            k = p_hat * two_n
            draws = rng.beta(k + _JEFFREYS, (two_n - k) + _JEFFREYS, size=n_draws)
        else:
            # No sample size -> can't put a CI on this allele; hold it fixed and say so.
            draws = np.full(n_draws, p_hat)
            fixed_no_n.append(allele)
        per_locus.setdefault(locus, []).append(draws)

    if not per_locus:
        warnings.append("no covering allele has a known frequency; coverage not sampled (NOT 0)")
        return np.empty(0), warnings

    not_covered = np.ones(n_draws)
    clamped = False
    for locus, arrs in per_locus.items():
        p_l = np.sum(arrs, axis=0)  # summed covering-allele freq at this locus, per draw
        over = p_l > 1.0
        if over.any():
            p_l = np.minimum(p_l, 1.0)
            clamped = True
        not_covered *= (1.0 - p_l) ** 2
    coverage = 1.0 - not_covered

    if fixed_no_n:
        warnings.append(
            "no AFND sample size for "
            f"{', '.join(sorted(fixed_no_n))}; held at point frequency (no sampling error "
            "contributed -- interval is optimistically narrow for these)"
        )
    if missing:
        warnings.append(
            f"excluded (no frequency, NOT treated as 0): {', '.join(sorted(missing))}"
        )
    if clamped:
        warnings.append("covering-allele frequency sum exceeded 1.0 at a locus in some draws; clamped")
    return coverage, warnings


def _summ(a: np.ndarray) -> dict:
    return {
        "median": float(np.median(a)),
        "ci95_low": float(np.percentile(a, 2.5)),
        "ci95_high": float(np.percentile(a, 97.5)),
    }


def effective_n_interval(
    *,
    incidence: float,
    antigen_numerator: int,
    antigen_denominator: int,
    covering_alleles,
    allele_freqs,
    sample_sizes=None,
    n_draws: int = 20000,
    seed: int = 0,
    incidence_rel_sd: float | None = None,
) -> MCResult:
    """Monte-Carlo interval on ``effective_N`` for one population.

    Propagates antigen-fraction (Beta) and HLA-coverage (per-allele Beta) sampling
    error. Incidence is fixed unless ``incidence_rel_sd`` is given, in which case it
    is drawn from a truncated-at-0 normal with that relative SD -- a *user-owned*
    assumption, flagged in the caveats. Returns a surfaced non-computable ``MCResult``
    (median ``None``) if coverage cannot be sampled, never a fabricated interval.
    """
    rng = np.random.default_rng(seed)
    caveats: list[str] = []

    cov_draws, cov_warn = sample_coverage(
        covering_alleles, allele_freqs, sample_sizes, n_draws, rng
    )
    caveats.extend(cov_warn)
    method = (
        "Monte-Carlo: antigen ~ Beta(numer+0.5, denom-numer+0.5), "
        "coverage ~ per-allele Beta from AFND sample size (Bui-2006 diploid per draw)"
    )
    if cov_draws.size == 0:
        return MCResult(None, None, None, n_draws, method, caveats)

    ag_draws = sample_antigen_fraction(antigen_numerator, antigen_denominator, n_draws, rng)

    if incidence_rel_sd is not None and incidence_rel_sd > 0:
        inc_draws = rng.normal(incidence, abs(incidence) * incidence_rel_sd, size=n_draws)
        inc_draws = np.clip(inc_draws, 0.0, None)
        method += f"; incidence ~ Normal(mu, {incidence_rel_sd:g}*mu) truncated>=0 (USER-supplied)"
        caveats.append(
            f"incidence uncertainty is a USER-supplied assumption (rel. SD {incidence_rel_sd:g}), "
            "not a sourced figure"
        )
    else:
        inc_draws = float(incidence)
        caveats.append(
            "incidence held FIXED (GLOBOCAN publishes no CI); interval reflects "
            "antigen + HLA-coverage sampling only -- true interval is wider"
        )

    eff = inc_draws * ag_draws * cov_draws
    s = _summ(eff)
    return MCResult(
        median=round(s["median"], 2),
        ci95_low=round(s["ci95_low"], 2),
        ci95_high=round(s["ci95_high"], 2),
        n_draws=n_draws,
        method=method,
        caveats=caveats,
        components={
            "antigen_fraction": {k: round(v, 6) for k, v in _summ(ag_draws).items()},
            "hla_coverage": {k: round(v, 6) for k, v in _summ(cov_draws).items()},
        },
    )
