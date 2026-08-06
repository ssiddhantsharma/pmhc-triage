"""Neo-peptide generation from a missense variant.

Given a protein sequence and a substitution (e.g. ``G12D``), produce every k-mer of
the requested lengths that *spans the mutated residue*, computed on the **mutated**
sequence. These are the candidate peptides handed to ``presentation``.

Hard guard (the load-bearing safety check): before mutating, we assert the
sequence's residue at the variant's position equals the variant's wild-type
residue. If it doesn't -- wrong isoform, wrong accession, or an off-by-one between
1-indexed variant nomenclature and 0-indexed strings -- we **refuse** and surface
the mismatch, rather than silently mutate the wrong residue and generate confident
nonsense. Substitution nomenclature is 1-indexed (``G12D`` = position 12).
"""

from __future__ import annotations

import re
from typing import Sequence as _Seq

from .provenance import Provenance, Sourced, today_iso

_SUBSTITUTION_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

DEFAULT_LENGTHS: tuple[int, ...] = (8, 9, 10, 11)


def parse_substitution(variant: str) -> tuple[str, int, str]:
    """Parse a single-residue substitution like ``G12D`` -> ``('G', 12, 'D')``.

    Position is 1-indexed. Raises ``ValueError`` on anything that isn't a clean
    ``<WT><pos><MUT>`` (we don't guess at indels/complex variants here).
    """
    m = _SUBSTITUTION_RE.match(variant.strip())
    if not m:
        raise ValueError(
            f"unrecognized substitution {variant!r}; expected e.g. 'G12D' (1-indexed)"
        )
    wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
    if pos < 1:
        raise ValueError("position must be 1-indexed (>= 1)")
    if wt == mut:
        raise ValueError(
            f"substitution {variant!r} does not change the residue (WT == MUT); "
            "not a neoepitope-generating variant"
        )
    return wt, pos, mut


def mutant_peptides(
    sequence: str,
    variant: str,
    lengths: _Seq[int] = DEFAULT_LENGTHS,
    *,
    source: str = "derived from sequence + variant",
) -> Sourced[list[str]]:
    """All ``lengths``-mers spanning ``variant``'s position, on the mutated sequence.

    Returns a :class:`~pmhc_triage.provenance.Sourced` list. Value is ``None``
    (with a surfaced warning) if the variant is malformed, out of range, or -- the
    important case -- the sequence's residue at the position does not match the
    variant's wild-type residue.
    """
    prov = Provenance(
        source=source,
        query_date=today_iso(),
        method=f"{tuple(lengths)}-mers spanning {variant} on the mutated sequence",
    )

    try:
        wt, pos, mut = parse_substitution(variant)
    except ValueError as exc:
        return Sourced(None, prov).warn(str(exc))

    if pos > len(sequence):
        return Sourced(None, prov).warn(
            f"variant position {pos} is beyond sequence length {len(sequence)}"
        )

    observed = sequence[pos - 1]  # 1-indexed -> 0-indexed
    if observed != wt:
        return Sourced(None, prov).warn(
            f"WT residue mismatch: variant {variant} expects {wt} at position {pos}, "
            f"but the sequence has {observed} there -- wrong isoform/accession or "
            "off-by-one? Refusing to generate peptides."
        )

    mutated = sequence[: pos - 1] + mut + sequence[pos:]
    mut_idx = pos - 1  # 0-indexed mutated position

    peptides: list[str] = []
    seen: set[str] = set()
    for length in lengths:
        # windows of this length that contain mut_idx and fit fully in the sequence
        start_min = max(0, mut_idx - length + 1)
        start_max = min(len(mutated) - length, mut_idx)
        for start in range(start_min, start_max + 1):
            pep = mutated[start : start + length]
            if len(pep) == length and pep not in seen:
                seen.add(pep)
                peptides.append(pep)

    result = Sourced(peptides, prov)
    if not peptides:
        result.warn(
            f"no peptides of lengths {tuple(lengths)} fit around position {pos} "
            f"in a sequence of length {len(sequence)}"
        )
    return result
