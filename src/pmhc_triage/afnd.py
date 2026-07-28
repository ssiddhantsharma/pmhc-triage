"""Load HLA allele frequencies from an AFND-format table (a parser, not a bundler).

Licensing note (load-bearing): the Allele Frequency Net Database (AFND) has been
distributed under a Creative Commons Attribution **Non-Commercial** license
(CC BY-NC) in its primary publications. This package's code is Apache-2.0, which
must not carry NC-licensed data inside it, so we deliberately **do not bundle or
redistribute AFND data**. This module instead parses a frequency table that *you*
export from AFND (or supply from another source), so the data-use terms stay
between you and AFND. Please cite AFND: http://www.allelefrequencies.net/publications.asp

The parser is format-flexible: point it at a CSV/TSV and name the allele,
population, and frequency columns. Frequencies are expected as fractions in
``[0, 1]`` (set ``percent=True`` if your column is 0-100).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .hla import normalize_allele
from .provenance import Provenance, today_iso


@dataclass
class FrequencyTable:
    """Allele frequencies grouped by population, plus provenance and any warnings.

    ``by_population`` maps population name -> {normalized allele: frequency}. Feed
    a single population's dict straight into :func:`pmhc_triage.hla.population_coverage`.

    ``sample_sizes`` maps population -> {normalized allele: n_individuals}, populated
    only when a sample-size column was supplied. It feeds the Monte-Carlo coverage CI
    (chromosomes = 2 x individuals); an empty dict means coverage uncertainty can't be
    propagated for that population (surfaced, never assumed).
    """

    by_population: dict[str, dict[str, float]]
    provenance: Provenance
    warnings: list[str] = field(default_factory=list)
    sample_sizes: dict[str, dict[str, float]] = field(default_factory=dict)

    def populations(self) -> list[str]:
        return sorted(self.by_population)

    def get(self, population: str) -> dict[str, float]:
        return self.by_population.get(population, {})

    def get_sample_sizes(self, population: str) -> dict[str, float]:
        return self.sample_sizes.get(population, {})


def load_afnd_frequencies(
    path: str | Path,
    *,
    allele_col: str = "allele",
    population_col: str = "population",
    freq_col: str = "allele_frequency",
    sample_size_col: str | None = "sample_size",
    sep: str | None = None,
    percent: bool = False,
    source: str | None = None,
) -> FrequencyTable:
    """Parse a user-supplied AFND-format frequency table.

    Column names are matched case-insensitively. ``sep=None`` sniffs the delimiter.
    Frequencies are fractions unless ``percent=True`` (then divided by 100).
    Problems are surfaced in ``FrequencyTable.warnings`` -- nothing is silently
    dropped or coerced.

    ``sample_size_col`` (number of typed individuals per allele row) is OPTIONAL: if
    that column is absent it is simply skipped (no error), and the resulting
    ``sample_sizes`` is empty so Monte-Carlo coverage CIs can't be computed -- a
    surfaced limitation, not a fabricated N.
    """
    path = Path(path)
    prov = Provenance(
        source=source or f"AFND-format export ({path.name})",
        query_date=today_iso(),
        method="user-supplied allele-frequency table (AFND CC BY-NC; not bundled)",
    )

    df = pd.read_csv(path, sep=sep, engine="python")
    lookup = {c.lower(): c for c in df.columns}

    def resolve(name: str) -> str:
        if name.lower() not in lookup:
            raise ValueError(
                f"column {name!r} not found; available columns: {list(df.columns)}"
            )
        return lookup[name.lower()]

    a_col, p_col, f_col = resolve(allele_col), resolve(population_col), resolve(freq_col)
    # sample_size is optional: resolve it only if the caller asked AND it exists.
    n_col = None
    if sample_size_col and sample_size_col.lower() in lookup:
        n_col = lookup[sample_size_col.lower()]

    warnings: list[str] = []
    by_pop: dict[str, dict[str, float]] = {}
    sizes: dict[str, dict[str, float]] = {}

    for _, row in df.iterrows():
        allele = normalize_allele(str(row[a_col]))
        pop = str(row[p_col]).strip()
        raw = str(row[f_col]).strip().rstrip("%")
        try:
            freq = float(raw)
        except ValueError:
            warnings.append(f"unparseable frequency {row[f_col]!r} for {allele}/{pop}; skipped")
            continue

        if percent:
            freq /= 100.0
        elif freq > 1.0:
            warnings.append(
                f"frequency {freq} > 1 for {allele}/{pop}; expected a fraction "
                "(pass percent=True if your column is 0-100). Kept as-is."
            )

        pop_map = by_pop.setdefault(pop, {})
        if allele in pop_map:
            warnings.append(
                f"duplicate {allele} in {pop!r}: kept first ({pop_map[allele]}), "
                f"ignored {freq} -- aggregate upstream if intended"
            )
            continue
        pop_map[allele] = freq

        if n_col is not None:
            try:
                n_ind = float(row[n_col])
                if n_ind > 0:
                    sizes.setdefault(pop, {})[allele] = n_ind
                else:
                    warnings.append(f"non-positive sample size for {allele}/{pop}; ignored")
            except (ValueError, TypeError):
                warnings.append(f"unparseable sample size {row[n_col]!r} for {allele}/{pop}; ignored")

    return FrequencyTable(by_pop, prov, warnings, sizes)
