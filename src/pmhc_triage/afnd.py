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
    """

    by_population: dict[str, dict[str, float]]
    provenance: Provenance
    warnings: list[str] = field(default_factory=list)

    def populations(self) -> list[str]:
        return sorted(self.by_population)

    def get(self, population: str) -> dict[str, float]:
        return self.by_population.get(population, {})


def load_afnd_frequencies(
    path: str | Path,
    *,
    allele_col: str = "allele",
    population_col: str = "population",
    freq_col: str = "allele_frequency",
    sep: str | None = None,
    percent: bool = False,
    source: str | None = None,
) -> FrequencyTable:
    """Parse a user-supplied AFND-format frequency table.

    Column names are matched case-insensitively. ``sep=None`` sniffs the delimiter.
    Frequencies are fractions unless ``percent=True`` (then divided by 100).
    Problems are surfaced in ``FrequencyTable.warnings`` -- nothing is silently
    dropped or coerced.
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

    warnings: list[str] = []
    by_pop: dict[str, dict[str, float]] = {}

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
        else:
            pop_map[allele] = freq

    return FrequencyTable(by_pop, prov, warnings)
