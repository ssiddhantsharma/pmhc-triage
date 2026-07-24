"""Disease burden (incidence) -- user-supplied or parsed from a table (not bundled).

There is no clean free API for cancer incidence, and GLOBOCAN/IARC and GBD carry
their own terms, so -- exactly as with AFND -- we do **not** bundle burden data.
You provide it: either a single manual value, or a CSV/TSV you export (e.g. from
GLOBOCAN) parsed here into per-population Sourced values. Cite your source.

Incidence is returned per population so it can join with HLA coverage on the same
population label in the score step. Mismatched population labels are the biggest
silent-error hazard, so keep the labels identical to those in your frequency table.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .provenance import Provenance, Sourced, today_iso


def manual_incidence(
    disease: str,
    population: str,
    cases_per_year: float,
    *,
    source: str,
    url: str | None = None,
) -> Sourced[float]:
    """Wrap a hand-entered incidence value with its (required) source."""
    prov = Provenance(
        source=source,
        url=url,
        query_date=today_iso(),
        method=f"user-supplied incidence for {disease} in {population}",
    )
    try:
        value = float(cases_per_year)
    except (TypeError, ValueError):
        return Sourced(None, prov).warn(f"non-numeric incidence {cases_per_year!r}")
    if value < 0:
        return Sourced(None, prov).warn(f"incidence must be >= 0, got {value}")
    return Sourced(value, prov)


def load_burden_table(
    path: str | Path,
    disease: str,
    *,
    disease_col: str = "disease",
    population_col: str = "population",
    incidence_col: str = "incidence",
    sep: str | None = None,
    source: str | None = None,
) -> dict[str, Sourced[float]]:
    """Parse a user-supplied incidence table, filtered to ``disease``, keyed by population.

    Column names match case-insensitively; ``sep=None`` sniffs the delimiter. The
    disease match is case-insensitive exact. Non-numeric incidences are surfaced as
    missing rather than dropped.
    """
    path = Path(path)
    df = pd.read_csv(path, sep=sep, engine="python")
    lookup = {c.lower(): c for c in df.columns}

    def resolve(name: str) -> str:
        if name.lower() not in lookup:
            raise ValueError(f"column {name!r} not found; available: {list(df.columns)}")
        return lookup[name.lower()]

    d_col, p_col, i_col = resolve(disease_col), resolve(population_col), resolve(incidence_col)
    src = source or f"burden table ({path.name})"

    out: dict[str, Sourced[float]] = {}
    for _, row in df.iterrows():
        if str(row[d_col]).strip().lower() != disease.strip().lower():
            continue
        pop = str(row[p_col]).strip()
        prov = Provenance(
            source=src,
            query_date=today_iso(),
            method=f"incidence for {disease} in {pop} (from {path.name})",
        )
        try:
            out[pop] = Sourced(float(str(row[i_col]).replace(",", "").strip()), prov)
        except ValueError:
            out[pop] = Sourced(None, prov).warn(f"non-numeric incidence {row[i_col]!r}")
    return out
