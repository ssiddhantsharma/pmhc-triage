"""Self cross-reactivity axis: a TCR / TCR-mimic raised against the target pMHC can
recognise a near-identical SELF peptide presentable on the same HLA. A
recognition-motif definition + an opt-in proteome self-scan.

Scope, stated honestly: the scan catches CLOSE mimics -- paralogs and allelic homologs
that share the TCR-contact face (e.g. MAGE-A12 KMAELVHFL vs the MAGE-A3 target
KVAELVHFL, differing only at the A*02:01 P2 anchor). It is BLIND to distant structural
mimicry: MAGE-A3 EVDPIGHLY and cardiac titin ESDPIVAQY differ at 4 contact positions,
so no sane Hamming threshold flags the pair that killed two patients. Absence of a hit
is therefore not a safety guarantee -- this is a risk-prioritisation flag for
experimental screening, and why the curated incidents list is shipped alongside."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np

from .hla import normalize_allele, parse_locus
from .provenance import Provenance, Sourced, today_iso

_DATA = Path(__file__).parent / "data" / "agents.yaml"
_AA = "ACDEFGHIKLMNPQRSTVWY"
_CODE = np.full(256, 20, dtype=np.int8)
for _i, _c in enumerate(_AA):
    _CODE[ord(_c)] = _i

# Primary non-C-terminal anchor(s), 1-indexed (PΩ added in tcr_contact_positions).
# A*01:01 anchors at P3, not P2 -- a blanket "P2" rule mis-motifs it. Src: SYFPEITHI/IEDB.
_ANCHORS = {"A*01:01": (3,), "A*02:01": (2,), "A*03:01": (2,), "A*11:01": (2,),
            "A*24:02": (2,), "B*07:02": (2,), "B*08:01": (3, 5)}
_LOCUS_ANCHOR = {"A": (2,), "B": (2,), "C": (2,)}


def tcr_contact_positions(length: int, allele: str = "") -> list[int]:
    """0-indexed TCR-facing positions: all but the HLA anchors (allele-specific + C-terminus)."""
    a = normalize_allele(allele) if allele else ""
    anchors = set(_ANCHORS.get(a) or _LOCUS_ANCHOR.get(parse_locus(a) if a else "", (2,)))
    anchors.add(length)
    return [i for i in range(length) if (i + 1) not in anchors]


def _encode(seq: str) -> np.ndarray:
    return _CODE[np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)]


@dataclass
class SelfHit:
    accession: str
    position: int       # 1-indexed start in the protein
    peptide: str
    mismatches: int     # at the TCR-contact positions only


def read_fasta(path: str) -> Iterator[tuple[str, str]]:
    acc, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    yield acc, "".join(seq)
                h = line[1:].split()[0]
                acc = h.split("|")[1] if h[:3] in ("sp|", "tr|") else h
                seq = []
            else:
                seq.append(line.strip())
    if acc:
        yield acc, "".join(seq)


def self_scan(peptide: str, proteome: Iterable[tuple[str, str]], *, allele: str = "",
              max_mismatch: int = 3, exclude: Iterable[str] = ()) -> list[SelfHit]:
    """Self peptides matching the target at its TCR-contact positions within max_mismatch.
    Anchors are treated as wildcards (a mimic may differ there yet still present) -- e.g.
    MAGE-A12 differs from the MAGE-A3 target only at the P2 anchor, so its contact face is
    identical and it flags at mm=0. numpy Hamming per protein; the exact target is skipped."""
    L = len(peptide)
    contact = tcr_contact_positions(L, allele)
    tgt = _encode(peptide)[contact]
    excl = set(exclude)
    hits: list[SelfHit] = []
    for acc, seq in proteome:
        if acc in excl or len(seq) < L:
            continue
        w = np.lib.stride_tricks.sliding_window_view(_encode(seq), L)[:, contact]
        mm = (w != tgt).sum(axis=1)
        for idx in np.nonzero(mm <= max_mismatch)[0]:
            frag = seq[idx:idx + L]
            if frag != peptide:
                hits.append(SelfHit(acc, int(idx) + 1, frag, int(mm[idx])))
    hits.sort(key=lambda h: (h.mismatches, h.accession, h.position))
    return hits


def antigen_safety(antigen_mode: str) -> tuple[str, str]:
    if antigen_mode == "mutation":
        return ("neoantigen", "mutation-derived, absent from the normal proteome -> tumour-specific; "
                "lower central-tolerance and on-target/off-tumour risk")
    return ("shared self", "shared/overexpressed self antigen -> central tolerance may have deleted "
            "high-affinity TCRs, and normal tissue expressing it is an on-target/off-tumour risk")


@dataclass
class Agent:
    name: str
    antigen: str
    peptide: str
    allele: str
    modality: str
    status: str = ""
    note: str = ""


def load_agents(antigen: str = "") -> tuple[list[Agent], list[Agent]]:
    """(clinical agents, cross-reactivity incidents); filtered to `antigen` if given."""
    if not _DATA.exists():
        return [], []
    import yaml
    data = yaml.safe_load(_DATA.read_text()) or {}
    agents = [Agent(**a) for a in data.get("agents", [])]
    incidents = [Agent(**a) for a in data.get("incidents", [])]
    if antigen:
        key = antigen.upper()
        agents = [a for a in agents if key in a.antigen.upper()]
        incidents = [a for a in incidents if key in a.antigen.upper()]
    return agents, incidents


def scan_sourced(peptide: str, allele: str, proteome_path: str, *, max_mismatch: int = 3,
                 exclude: Iterable[str] = ()) -> Sourced[list[SelfHit]]:
    """self_scan wrapped with provenance for the report."""
    contact = tcr_contact_positions(len(peptide), allele)
    anchors = [f"P{i + 1}" for i in range(len(peptide)) if i not in contact]
    prov = Provenance(source=f"self-scan vs {Path(proteome_path).name}", query_date=today_iso(),
                      method=f"Hamming <= {max_mismatch} at TCR-contact positions "
                             f"(anchors {anchors} as wildcards)")
    try:
        hits = self_scan(peptide, read_fasta(proteome_path), allele=allele,
                         max_mismatch=max_mismatch, exclude=exclude)
    except OSError as exc:
        return Sourced(None, prov).warn(f"could not read proteome: {exc}")
    return Sourced(hits, prov)
