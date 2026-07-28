"""Fetch allele frequencies from AFND *on the user's machine* (never bundled).

The tool deliberately ships no AFND data (CC BY-NC; see :mod:`pmhc_triage.afnd`).
This module lowers the setup friction without breaking that rule: it fetches the
frequencies you need from AFND at *your* request, into a local file, so the data-use
terms stay between you and AFND -- the tool only automates the retrieval you would
otherwise do by hand. The output is exactly the ``--freqs`` format
(:func:`pmhc_triage.afnd.load_afnd_frequencies` reads it), including the
``sample_size`` column that powers the Monte-Carlo coverage CI.

Robustness notes (learned the hard way):
- AFND paginates very common alleles (e.g. ``A*02:01`` spans 500+ cohorts), so a
  bare per-allele query can silently *miss* a target cohort. We query per
  ``(allele, country)`` so the target cohort is always on the returned page.
- The results-table column order is fixed: ``[1]=allele, [3]=population,
  [5]=allele_frequency, [7]=sample_size`` (0-indexed ``<td>``). Verified against the
  live site; parsed exactly, never guessed.
- A requested ``(allele, population)`` with no AFND row is **surfaced** (listed in
  warnings), never written as a fabricated 0.
"""

from __future__ import annotations

import html as _html
import re
import urllib.parse
from dataclasses import dataclass, field

import httpx

from .hla import normalize_allele
from .provenance import today_iso

_AFND_URL = "http://www.allelefrequencies.net/hla6006a.asp"


@dataclass(frozen=True)
class ReferencePopulation:
    """A label mapped to a specific AFND cohort (chosen deliberately, not averaged).

    ``afnd_population`` must match the AFND 'Population' cell exactly; ``country`` is
    the AFND country filter used to defeat pagination on common alleles.
    """

    label: str
    afnd_population: str
    country: str


# Large, well-powered reference cohorts spanning the A*11:01 frequency gradient
# (Europe low / East Asia high / South Asia intermediate). Deliberate defaults --
# override for a different reference set. N's (2024): DKMS ~3.46M, HK BMDR ~7.6k,
# India South UCBB ~11.4k.
DEFAULT_POPULATIONS = (
    ReferencePopulation("Europe", "Germany DKMS - German donors", "Germany"),
    ReferencePopulation("EastAsia", "Hong Kong Chinese BMDR", "Hong Kong"),
    ReferencePopulation("SouthAsia", "India South UCBB", "India"),
)

# A compact common class-I panel (covers the frequent A/B/C alleles behind most
# characterized neoantigen restrictions). Not exhaustive -- pass your own list.
COMMON_CLASS_I = (
    "A*01:01", "A*02:01", "A*03:01", "A*11:01", "A*24:02", "A*26:01", "A*32:01",
    "B*07:02", "B*08:01", "B*15:01", "B*35:01", "B*40:01", "B*44:02", "B*44:03",
    "C*03:04", "C*04:01", "C*06:02", "C*07:01", "C*07:02",
)


@dataclass
class FetchResult:
    """Rows fetched plus everything that was missing or ambiguous (surfaced)."""

    rows: list[tuple]  # (allele, population_label, allele_frequency, sample_size)
    warnings: list[str] = field(default_factory=list)
    query_date: str = ""

    def write_tsv(self, path) -> int:
        """Write the ``--freqs`` TSV (allele/population/allele_frequency/sample_size)."""
        with open(path, "w") as fo:
            fo.write("allele\tpopulation\tallele_frequency\tsample_size\n")
            for allele, pop, freq, n in self.rows:
                fo.write(f"{allele}\t{pop}\t{freq}\t{n}\n")
        return len(self.rows)


def parse_afnd_table(page_html: str, allele: str) -> dict[str, tuple[float, int]]:
    """Parse one AFND results page: ``{population: (allele_frequency, sample_size)}``.

    Pure and hermetic -- unit-testable on a saved HTML snippet without the network.
    If a population appears more than once (multiple studies), the largest-N row wins.
    """
    want = normalize_allele(allele)
    best: dict[str, tuple[float, int]] = {}
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", page_html, re.S | re.I):
        cells = [
            _html.unescape(re.sub(r"<[^>]*>", "", c)).replace("\xa0", "").strip()
            for c in re.findall(r"<td[^>]*>(.*?)</td>", row, re.S | re.I)
        ]
        if len(cells) < 8:
            continue
        if normalize_allele(cells[1]) != want and not cells[1].startswith(want):
            continue
        pop = cells[3].strip()
        try:
            freq = float(cells[5])
            n = int(cells[7].replace(",", ""))
        except (ValueError, IndexError):
            continue
        if not pop or n <= 0:
            continue
        if pop not in best or n > best[pop][1]:
            best[pop] = (freq, n)
    return best


def _fetch_page(allele: str, country: str, client: httpx.Client, timeout: float) -> str | None:
    params = {"hla_locus_type": "Classical", "hla_allele1": allele, "hla_country": country}
    url = f"{_AFND_URL}?{urllib.parse.urlencode(params)}"
    try:
        r = client.get(url, headers={"User-Agent": "Mozilla/5.0 pmhc-triage"}, timeout=timeout)
    except httpx.HTTPError:
        return None
    return r.text if r.status_code == 200 else None


def fetch_frequencies(
    alleles=COMMON_CLASS_I,
    populations=DEFAULT_POPULATIONS,
    *,
    client: httpx.Client | None = None,
    timeout: float = 40.0,
) -> FetchResult:
    """Fetch ``alleles`` × ``populations`` from AFND into a ``--freqs`` table.

    Queries per ``(allele, country)`` to survive AFND pagination, keeps only the exact
    ``afnd_population`` cohort per label, and surfaces every missing ``(allele,
    population)`` pair in ``warnings`` (never writes a fabricated 0). Returns a
    :class:`FetchResult`; call ``.write_tsv(path)`` to save it.
    """
    alleles = [normalize_allele(a) for a in alleles]
    populations = list(populations)
    owns = client is None
    client = client or httpx.Client()
    rows: list[tuple] = []
    warnings: list[str] = []
    # group populations by country so each allele needs one fetch per distinct country
    by_country: dict[str, list[ReferencePopulation]] = {}
    for p in populations:
        by_country.setdefault(p.country, []).append(p)
    try:
        for allele in alleles:
            for country, pops in by_country.items():
                page = _fetch_page(allele, country, client, timeout)
                if page is None:
                    for p in pops:
                        warnings.append(f"fetch failed for {allele} in {country} -> {p.label} skipped")
                    continue
                table = parse_afnd_table(page, allele)
                for p in pops:
                    if p.afnd_population in table:
                        freq, n = table[p.afnd_population]
                        rows.append((allele, p.label, freq, n))
                    else:
                        warnings.append(
                            f"no AFND row for {allele} in {p.afnd_population!r} "
                            f"({p.label}); NOT written (surfaced, not 0)"
                        )
    finally:
        if owns:
            client.close()
    return FetchResult(rows=rows, warnings=warnings, query_date=today_iso())
