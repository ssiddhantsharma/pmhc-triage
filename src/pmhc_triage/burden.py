"""Disease burden (incidence): user-supplied value, parsed table, or the cited bundle.

Cancer incidence carries source-specific terms (GLOBOCAN/IARC, GBD), so most of it is
user-supplied; a small cited starter bundle ships for convenience. Incidence is keyed
by population so it joins HLA coverage on the same label in the score step. Keep the
labels identical to the frequency table: label mismatch is the main silent-error hazard.
"""

from __future__ import annotations

import io
import math
from importlib.resources import files
from pathlib import Path

import pandas as pd

from .provenance import Provenance, Sourced, today_iso


def _bundle_df() -> pd.DataFrame:
    text = files("pmhc_triage").joinpath("data/burden_bundle.csv").read_text()
    return pd.read_csv(io.StringIO(text))


def bundled_diseases() -> list[str]:
    """Disease labels available in the shipped GLOBOCAN-2022 starter bundle."""
    return sorted(_bundle_df()["disease"].unique().tolist())


def load_bundled(disease: str, *, population: str = "World") -> Sourced[float]:
    """Incidence from the shipped curated starter bundle (cited).

    Convenience starter values (GCO Cancer Today 2024) with the citation carried in the
    provenance. Covers ``World`` plus per-region rows (``Europe``/``EastAsia``/``SouthAsia``,
    populated from Germany/China/India country figures as regional proxies — surfaced in
    the note). Refine per-region for a real analysis. Returns a surfaced *missing* if the
    disease/population isn't in the bundle.
    """
    df = _bundle_df()
    prov = Provenance(source="pmhc-triage curated bundle", query_date=today_iso(),
                      method=f"bundled GLOBOCAN-2022 incidence for {disease!r} / {population!r}")
    hits = df[(df["disease"].str.lower() == disease.strip().lower())
              & (df["population"].str.lower() == population.strip().lower())]
    if hits.empty:
        return Sourced(None, prov).warn(
            f"{disease!r}/{population!r} not in bundle; available: {bundled_diseases()}"
        )
    row = hits.iloc[0]
    result = Sourced(
        float(row["incidence"]),
        Provenance(source=str(row["source"]), query_date=today_iso(),
                   method=f"curated GLOBOCAN-2022 bundle: {row['note']}"),
    )
    # surface the bundle's caveats at runtime, not just in the CSV
    result.warn(f"bundle note: {row['note']}")
    if str(row["population"]).strip().lower() == "world":
        result.warn(
            "World-level incidence: joins only with a 'World' HLA-coverage population. "
            "For population-specific effective-N, supply per-region incidence instead."
        )
    return result


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
    if not math.isfinite(value) or value < 0:
        return Sourced(None, prov).warn(f"incidence must be finite and >= 0, got {value}")
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
            value = float(str(row[i_col]).replace(",", "").strip())
        except ValueError:
            out[pop] = Sourced(None, prov).warn(f"non-numeric incidence {row[i_col]!r}")
            continue
        # reject non-finite/negative (a negative would become a negative effective_N)
        if not math.isfinite(value) or value < 0:
            out[pop] = Sourced(None, prov).warn(f"invalid incidence {row[i_col]!r} (must be finite and >= 0)")
            continue
        out[pop] = Sourced(value, prov)
    return out
