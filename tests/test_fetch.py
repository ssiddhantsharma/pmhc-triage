"""Tests for the AFND fetcher. Parsing is hermetic (synthetic HTML mirroring the real
AFND results-table column order); the network path is exercised with a fake client."""

import httpx
import pytest

from pmhc_triage.fetch import (
    DEFAULT_POPULATIONS,
    ReferencePopulation,
    fetch_frequencies,
    parse_afnd_table,
)


# One AFND-style data row: cols [0]=line [1]=allele [3]=population [5]=freq [7]=N
def _row(allele, pop, freq, n):
    return (
        f"<tr><td>1</td><td>{allele}</td><td></td><td>{pop}</td>"
        f"<td>0.0</td><td>{freq}</td><td></td><td>{n}</td>"
        f"<td>See</td><td></td><td></td><td></td></tr>"
    )


def _page(*rows):
    return "<html><table>" + "".join(rows) + "</table></html>"


def test_parse_extracts_freq_and_sample_size():
    html = _page(
        _row("A*11:01", "Germany DKMS - German donors", "0.0528", "3,456,066"),
        _row("A*11:01", "Hong Kong Chinese BMDR", "0.2968", "7595"),
    )
    t = parse_afnd_table(html, "A*11:01")
    assert t["Germany DKMS - German donors"] == (0.0528, 3456066)  # comma parsed
    assert t["Hong Kong Chinese BMDR"] == (0.2968, 7595)


def test_parse_keeps_largest_n_on_duplicate_population():
    html = _page(
        _row("A*02:01", "Germany DKMS - German donors", "0.30", "100"),
        _row("A*02:01", "Germany DKMS - German donors", "0.2839", "3456066"),
    )
    t = parse_afnd_table(html, "A*02:01")
    assert t["Germany DKMS - German donors"] == (0.2839, 3456066)  # bigger cohort wins


def test_parse_ignores_other_alleles_and_bad_rows():
    html = _page(
        _row("A*03:01", "Germany DKMS - German donors", "0.151", "3456066"),
        _row("A*11:01", "Germany DKMS - German donors", "0.0528", "3456066"),
        "<tr><td>x</td></tr>",  # malformed, too few cells
    )
    t = parse_afnd_table(html, "A*11:01")
    assert list(t) == ["Germany DKMS - German donors"]
    assert t["Germany DKMS - German donors"][0] == 0.0528


class _FakeClient:
    """Returns a page per (allele, country) from a lookup; mimics httpx.Client.get."""

    def __init__(self, pages):
        self.pages = pages

    def get(self, url, headers=None, timeout=None):
        # crude allele+country extraction from the query string
        import urllib.parse
        q = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
        allele = q["hla_allele1"][0]
        country = q["hla_country"][0]
        body = self.pages.get((allele, country), "<html></html>")
        return httpx.Response(200, text=body)

    def close(self):
        pass


def test_fetch_writes_rows_and_surfaces_missing(tmp_path):
    pops = (
        ReferencePopulation("Europe", "Germany DKMS - German donors", "Germany"),
        ReferencePopulation("EastAsia", "Hong Kong Chinese BMDR", "Hong Kong"),
    )
    pages = {
        ("A*11:01", "Germany"): _page(_row("A*11:01", "Germany DKMS - German donors", "0.0528", "3456066")),
        ("A*11:01", "Hong Kong"): _page(_row("A*11:01", "Hong Kong Chinese BMDR", "0.2968", "7595")),
        # A*02:01 present in Germany but MISSING in Hong Kong page -> must be surfaced
        ("A*02:01", "Germany"): _page(_row("A*02:01", "Germany DKMS - German donors", "0.2839", "3456066")),
        ("A*02:01", "Hong Kong"): _page(),  # empty
    }
    res = fetch_frequencies(["A*11:01", "A*02:01"], pops, client=_FakeClient(pages))
    got = {(a, p): (f, n) for a, p, f, n in res.rows}
    assert got[("A*11:01", "Europe")] == (0.0528, 3456066)
    assert got[("A*02:01", "Europe")] == (0.2839, 3456066)
    assert ("A*02:01", "EastAsia") not in got                 # missing, NOT written as 0
    assert any("A*02:01" in w and "EastAsia" in w for w in res.warnings)

    n = res.write_tsv(tmp_path / "freqs.tsv")
    assert n == 3
    # round-trips into the tool's own frequency loader
    from pmhc_triage.afnd import load_afnd_frequencies
    ft = load_afnd_frequencies(tmp_path / "freqs.tsv")
    assert ft.get("Europe")["A*11:01"] == 0.0528
    assert ft.get_sample_sizes("Europe")["A*11:01"] == 3456066


def test_default_populations_are_distinct_labels():
    labels = [p.label for p in DEFAULT_POPULATIONS]
    assert len(labels) == len(set(labels))
