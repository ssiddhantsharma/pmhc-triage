"""HLA Ligand Atlas: peptides actually eluted from benign human tissue.

The proteome self-scan asks whether a near-self mimic EXISTS in the translated proteome. Most
of the proteome is never processed and presented, so existence is weak evidence. This asks the
stronger question: was the mimic actually observed on benign tissue by mass spectrometry?

Marcu et al., "HLA Ligand Atlas: a benign reference of HLA-presented peptides to improve
T-cell-based cancer immunotherapy", J Immunother Cancer 9:e002071 (2021).
https://hla-ligand-atlas.org  --  data CC-BY 4.0, so it is fetched, not redistributed here.

WHAT THE DATA DOES AND DOES NOT SAY (from the release README, not inferred):

  - the peptide WAS measured on benign tissue. That part is evidence, not prediction.
  - `donor_alleles` is formatted `[swn]/{gene}*nn:nn`, where s/w/n is a PREDICTED strong, weak
    or non-binder call against an allele the donor carried. It is not a measured deconvolution,
    so "presented BY this allele" is not what the file states.
  - the README is explicit that tissue-to-allele pairing cannot be deduced from the aggregated
    table. A hit therefore reports tissues and allele-binding calls as two separate facts, and
    never as "presented on tissue T by allele A".
"""

from __future__ import annotations

import io
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

RELEASE = "2020.12"
ZIP_URL = "https://hla-ligand-atlas.org/rel/hla_{release}.zip"
_MEMBER = "hla_{release}/HLA_aggregated.tsv"
_HEADER = "peptide\tbinding_calls\ttissues\n"


@dataclass(frozen=True)
class BenignLigand:
    """One peptide eluted from benign tissue."""

    peptide: str
    binding: dict[str, str]   # allele -> 's' | 'w' | 'n'  (PREDICTED, see module docstring)
    tissues: tuple[str, ...]  # measured, but not resolvable to a specific allele

    def binds(self, allele: str, *, weak_counts: bool = True) -> bool:
        call = self.binding.get(allele)
        return call == "s" or (weak_counts and call == "w")


def _parse_calls(field: str) -> dict[str, str]:
    """`s/A*01:01,n/B*07:02` -> {'A*01:01': 's', 'B*07:02': 'n'}. Keeps the strongest call."""
    rank = {"s": 2, "w": 1, "n": 0}
    out: dict[str, str] = {}
    for item in field.split(","):
        item = item.strip()
        if "/" not in item:
            continue
        call, allele = item.split("/", 1)
        if call in rank and (allele not in out or rank[call] > rank[out[allele]]):
            out[allele] = call
    return out


def fetch_atlas(out_path, *, release: str = RELEASE, timeout: float = 300.0) -> int:
    """Download the Atlas and write the compact HLA-I table this package reads.

    The published artifact is a ~12 MB zip whose aggregated table is ~66 MB and mostly HLA-II;
    only the HLA-I rows are kept, which is what a class-I epitope scan can use.
    """
    import httpx

    url = ZIP_URL.format(release=release)
    with httpx.Client(follow_redirects=True, timeout=timeout) as client:
        resp = client.get(url, headers={"User-Agent": "pmhc-triage"})
        resp.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(resp.content)) as z:
        member = _MEMBER.format(release=release)
        if member not in z.namelist():
            raise RuntimeError(f"{member} not in {url}; release layout changed")
        n = 0
        with z.open(member) as fh, open(out_path, "w") as out:
            out.write(_HEADER)
            next(fh)  # header
            for raw in fh:
                parts = raw.decode("utf-8", "replace").rstrip("\n").split("\t")
                if len(parts) < 5 or parts[2] != "HLA-I":
                    continue
                out.write(f"{parts[1]}\t{parts[3]}\t{parts[4]}\n")
                n += 1
    return n


def load_atlas(path) -> Iterator[BenignLigand]:
    """Stream the compact table written by ``fetch_atlas``."""
    p = Path(path)
    with open(p) as fh:
        first = fh.readline()
        if not first.startswith("peptide\t"):
            raise ValueError(f"{p} is not a pmhc-triage ligand table (run `fetch-ligands`)")
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 3 or not parts[0]:
                continue
            tissues = tuple(t.strip() for t in parts[2].split(",") if t.strip())
            yield BenignLigand(parts[0], _parse_calls(parts[1]), tissues)
