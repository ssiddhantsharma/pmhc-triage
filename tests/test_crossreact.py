"""Self cross-reactivity axis: recognition motif, self-scan, curated agents.

Two primary-source-verified landmines anchor these tests, and they cut opposite ways:

- MAGE-A12 (UniProt P43365 has KMAELVHFL at 112) differs from the MAGE-A3 target
  KVAELVHFL only at P2 -- the A*02:01 anchor. Its TCR-contact face is IDENTICAL, so the
  scan flags it at zero mismatch. This is the anchor-wildcard positive control.
- Titin ESDPIVAQY differs from MAGE-A3 EVDPIGHLY at 4 TCR-contact positions on A*01:01.
  A Hamming self-scan at a sane threshold does NOT catch it -- distant structural mimicry
  is the scan's blind spot, which is exactly why curated incidents + experimental
  screening stay mandatory. Absence of a hit is not a safety guarantee.
"""

from pmhc_triage import cli
from pmhc_triage import crossreact as cr

MAGEA3_A1 = "EVDPIGHLY"   # A*01:01
TITIN = "ESDPIVAQY"       # cardiac cross-reactant, contact-distance 4
MAGEA3_A2 = "KVAELVHFL"   # A*02:01
MAGEA12 = "KMAELVHFL"     # differs only at P2 anchor


def test_a01_anchors_p3_not_p2():
    # 0-indexed; anchors excluded are P3 (idx 2) and the C-terminus (idx 8).
    contact = cr.tcr_contact_positions(9, "A*01:01")
    assert 2 not in contact and 8 not in contact
    assert 1 in contact  # P2 is a TCR-contact for A*01:01, unlike A*02:01


def test_a02_anchors_p2():
    contact = cr.tcr_contact_positions(9, "A*02:01")
    assert 1 not in contact and 8 not in contact
    assert 2 in contact  # P3 is a contact for A*02:01


def test_unknown_allele_falls_back_to_locus_default():
    assert cr.tcr_contact_positions(9, "A*99:99") == [0, 2, 3, 4, 5, 6, 7]


def test_scan_catches_anchor_only_mimic_magea12():
    # KMAELVHFL vs KVAELVHFL differ only at the A*02:01 P2 anchor -> identical contact
    # face -> flagged at zero mismatch even though the raw sequences differ.
    hits = cr.self_scan(MAGEA3_A2, [("MAGEA12", MAGEA12)], allele="A*02:01", max_mismatch=0)
    assert hits and hits[0].mismatches == 0 and hits[0].peptide == MAGEA12


def test_scan_blind_to_distant_titin_mimicry():
    # The canonical failure mode a Hamming scan cannot see: 4 contact-position differences,
    # so the default threshold does NOT flag titin. This is a documented limitation.
    assert cr.self_scan(MAGEA3_A1, [("TITIN", TITIN)], allele="A*01:01", max_mismatch=3) == []
    caught = cr.self_scan(MAGEA3_A1, [("TITIN", TITIN)], allele="A*01:01", max_mismatch=4)
    assert caught and caught[0].mismatches == 4


def test_scan_finds_embedded_window():
    prot = [("P", "MMMM" + MAGEA12 + "KKKK")]
    hits = cr.self_scan(MAGEA3_A2, prot, allele="A*02:01", max_mismatch=0)
    assert hits and hits[0].position == 5  # 1-indexed start


def test_self_scan_skips_the_exact_target():
    assert cr.self_scan(MAGEA3_A2, [("SELF", MAGEA3_A2)], allele="A*02:01", max_mismatch=3) == []


def test_hits_sorted_by_mismatch():
    prot = [("FAR", "KAAELVAFL"), ("NEAR", MAGEA12)]
    hits = cr.self_scan(MAGEA3_A2, prot, allele="A*02:01", max_mismatch=4)
    assert [h.mismatches for h in hits] == sorted(h.mismatches for h in hits)


def test_antigen_safety_modes():
    assert cr.antigen_safety("mutation")[0] == "neoantigen"
    assert cr.antigen_safety("expression")[0] == "shared self"


def test_read_fasta_parses_uniprot_headers(tmp_path):
    fa = tmp_path / "p.fasta"
    fa.write_text(">sp|P12345|NAME_HUMAN desc\nMKV\nLLL\n>bare\nAAA\n")
    assert list(cr.read_fasta(str(fa))) == [("P12345", "MKVLLL"), ("bare", "AAA")]


def test_load_agents_and_incidents():
    agents, incidents = cr.load_agents()
    assert any(a.peptide == "YLEPGPVTV" for a in agents)  # tebentafusp
    assert any("titin" in i.note.lower() for i in incidents)


def test_load_agents_filters_by_antigen():
    agents, incidents = cr.load_agents("MAGE-A3")
    assert agents == []
    assert incidents and all("MAGE-A3" in i.antigen for i in incidents)


def test_scan_sourced_carries_provenance(tmp_path):
    fa = tmp_path / "p.fasta"
    fa.write_text(f">sp|MAGEA12|T\n{MAGEA12}\n")
    s = cr.scan_sourced(MAGEA3_A2, "A*02:01", str(fa), max_mismatch=0)
    assert not s.is_missing and s.provenance.method and s.value[0].peptide == MAGEA12


def test_scan_sourced_missing_file_surfaces_not_crashes():
    s = cr.scan_sourced(MAGEA3_A2, "A*02:01", "/no/such/file.fasta")
    assert s.is_missing and s.warnings


def test_cli_crossreact_runs(capsys, tmp_path):
    fa = tmp_path / "p.fasta"
    fa.write_text(f">sp|MAGEA12|T\nMMMM{MAGEA12}KKKK\n")
    rc = cli.main(["crossreact", "--peptide", MAGEA3_A2, "--allele", "A*02:01",
                   "--proteome", str(fa), "--max-mismatch", "0"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "MAGEA12" in out and "Recognition motif" in out


def test_cli_crossreact_without_proteome(capsys):
    rc = cli.main(["crossreact", "--peptide", MAGEA3_A1, "--allele", "A*01:01"])
    assert rc == 0
    assert "--proteome" in capsys.readouterr().out  # prompts for the opt-in scan
