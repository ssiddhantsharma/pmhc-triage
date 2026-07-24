"""Protein sequences from UniProt -- always the primary source, never a summary.

Text/LLM summaries of UniProt entries can misreport sequence length and invent
residue boundaries, which would silently shift every downstream peptide index. So
this module pulls the canonical sequence straight from the UniProt REST ``.fasta``
endpoint and records the exact URL + query date in the provenance. It never
guesses a sequence.
"""

from __future__ import annotations

import httpx

from .provenance import Provenance, Sourced, fetched_at_or_today, today_iso

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"


def fetch_uniprot_sequence(
    accession: str,
    *,
    client: httpx.Client | None = None,
    timeout: float = 30.0,
) -> Sourced[str]:
    """Fetch the canonical protein sequence for a UniProt ``accession`` (e.g. ``P01116``).

    Returns a :class:`~pmhc_triage.provenance.Sourced` string. On any failure
    (network error, non-200, malformed body) the value is ``None`` and the reason
    is surfaced in ``warnings`` -- never a silent empty string.

    ``client`` may be injected (e.g. an ``httpx.Client`` with a mock transport)
    for hermetic testing.
    """
    url = UNIPROT_FASTA.format(accession=accession)

    def prov_with(date: str) -> Provenance:
        return Provenance(source=f"UniProt {accession}", url=url, query_date=date,
                          method="REST .fasta canonical sequence")

    prov = prov_with(today_iso())
    owns_client = client is None
    client = client or httpx.Client(timeout=timeout, follow_redirects=True)
    try:
        resp = client.get(url)
    except httpx.HTTPError as exc:
        return Sourced(None, prov).warn(f"request to UniProt failed: {exc}")
    finally:
        if owns_client:
            client.close()

    if resp.status_code != 200:
        return Sourced(None, prov).warn(
            f"HTTP {resp.status_code} from UniProt for accession {accession!r}"
        )

    lines = resp.text.splitlines()
    if not lines or not lines[0].startswith(">"):
        return Sourced(None, prov).warn("unexpected response: no FASTA header line")

    seq = "".join(ln.strip() for ln in lines[1:] if ln and not ln.startswith(">"))
    if not seq:
        return Sourced(None, prov).warn("FASTA had a header but no sequence residues")

    return Sourced(seq, prov_with(fetched_at_or_today(resp)))
