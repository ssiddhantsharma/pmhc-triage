"""Hermetic tests for the UniProt fetcher using an httpx mock transport (no network)."""

import httpx

from pmhc_triage.sequences import fetch_uniprot_sequence

KRAS_FASTA = (
    ">sp|P01116|RASK_HUMAN GTPase KRas OS=Homo sapiens OX=9606 GN=KRAS PE=1 SV=1\n"
    "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG\n"
    "QEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDL\n"
    "PSRTVDTKQAQDLARSYGIPFIETSAKTRQRVEDAFYTLVREIRQYRLKKISKEEKTPGC\n"
    "VKIKKCIIM\n"
)


def _client(handler):
    return httpx.Client(transport=httpx.MockTransport(handler))


def test_fetch_parses_fasta():
    def handler(request):
        assert "P01116.fasta" in str(request.url)
        return httpx.Response(200, text=KRAS_FASTA)

    s = fetch_uniprot_sequence("P01116", client=_client(handler))
    assert not s.is_missing
    assert s.value.startswith("MTEYKLVVVGAGGVGKS")
    assert len(s.value) == 189
    assert s.value[11] == "G"  # residue 12, 1-indexed
    assert s.provenance.url.endswith("P01116.fasta")


def test_fetch_404_surfaces_missing():
    def handler(request):
        return httpx.Response(404, text="Not found")

    s = fetch_uniprot_sequence("NOPE", client=_client(handler))
    assert s.is_missing
    assert any("HTTP 404" in w for w in s.warnings)


def test_fetch_malformed_body_surfaces():
    def handler(request):
        return httpx.Response(200, text="garbage without a fasta header")

    s = fetch_uniprot_sequence("P01116", client=_client(handler))
    assert s.is_missing
    assert any("FASTA header" in w for w in s.warnings)


def test_network_error_surfaces_not_raises():
    def handler(request):
        raise httpx.ConnectError("boom")

    s = fetch_uniprot_sequence("P01116", client=_client(handler))
    assert s.is_missing
    assert any("failed" in w for w in s.warnings)
