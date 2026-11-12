"""Benign-tissue ligand scan (HLA Ligand Atlas).

Hermetic: a handful of synthetic rows in the real column format. The 12 MB release is never
fetched here. Anchored on the same MAGE-A3 target the proteome scan uses, so the two axes are
directly comparable: KVAELVHFL on A*02:01, whose P2 is the anchor.
"""

import pytest

from pmhc_triage.crossreact import benign_sourced, scan_benign_ligands
from pmhc_triage.ligandatlas import BenignLigand, load_atlas

# real format: `[swn]/{gene}*nn:nn` comma separated, then comma-separated tissues
ROWS = """peptide\tbinding_calls\ttissues
KMAELVHFL\ts/A*02:01,n/B*07:02\tAorta, Kidney
SLAELVHAV\ts/A*02:01\tAdrenal gland, Aorta
QQQQQQQQQ\ts/A*02:01\tThymus
KVAELVHFLX\ts/A*02:01\tLung
NMAELVHFL\tn/A*02:01\tSpleen
"""


def _atlas(tmp_path):
    p = tmp_path / "benign.tsv"
    p.write_text(ROWS)
    return p


def test_parses_binding_calls_and_tissues(tmp_path):
    ligs = {lig.peptide: lig for lig in load_atlas(_atlas(tmp_path))}
    lig = ligs["KMAELVHFL"]
    assert lig.binding == {"A*02:01": "s", "B*07:02": "n"}
    assert lig.tissues == ("Aorta", "Kidney")
    assert lig.binds("A*02:01") and not lig.binds("B*07:02")
    # weak counts by default, and does not when asked not to
    weak = BenignLigand("AAA", {"A*02:01": "w"}, ())
    assert weak.binds("A*02:01") and not weak.binds("A*02:01", weak_counts=False)


def test_strongest_call_wins_when_an_allele_repeats():
    from pmhc_triage.ligandatlas import _parse_calls
    assert _parse_calls("n/A*02:01,s/A*02:01,w/A*02:01")["A*02:01"] == "s"


def test_scan_matches_contact_face_and_respects_length_and_binding(tmp_path):
    hits = scan_benign_ligands("KVAELVHFL", load_atlas(_atlas(tmp_path)),
                               allele="A*02:01", max_mismatch=1)
    found = {h.peptide for h in hits}
    # KMAELVHFL differs ONLY at P2, the A*02:01 anchor, which is wildcarded -> zero mismatches
    assert "KMAELVHFL" in found
    assert next(h for h in hits if h.peptide == "KMAELVHFL").mismatches == 0
    assert "QQQQQQQQQ" not in found, "unrelated peptide must not match"
    assert "KVAELVHFLX" not in found, "different length must not be compared"
    assert "NMAELVHFL" not in found, "predicted non-binder dropped unless --all-benign"

    everything = scan_benign_ligands("KVAELVHFL", load_atlas(_atlas(tmp_path)),
                                     allele="A*02:01", max_mismatch=1, require_binding=False)
    assert "NMAELVHFL" in {h.peptide for h in everything}


def test_hit_carries_tissue_and_the_predicted_call(tmp_path):
    hits = scan_benign_ligands("KVAELVHFL", load_atlas(_atlas(tmp_path)),
                               allele="A*02:01", max_mismatch=0)
    h = next(h for h in hits if h.peptide == "KMAELVHFL")
    assert h.binding_call == "s" and "Aorta" in h.tissues


def test_sourced_scan_states_that_tissue_cannot_be_paired_to_allele(tmp_path):
    """The release README says so explicitly; claiming 'presented on Aorta by A*02:01' would
    overstate what the aggregated table supports."""
    s = benign_sourced("KVAELVHFL", "A*02:01", str(_atlas(tmp_path)), max_mismatch=1)
    assert s.value and any("cannot be deduced" in w for w in s.warnings)
    assert "benign tissue" in s.provenance.method


def test_missing_or_wrong_file_is_surfaced_not_raised(tmp_path):
    s = benign_sourced("KVAELVHFL", "A*02:01", str(tmp_path / "nope.tsv"))
    assert s.is_missing and s.warnings

    bad = tmp_path / "bad.tsv"
    bad.write_text("something\telse\n")
    assert benign_sourced("KVAELVHFL", "A*02:01", str(bad)).is_missing

    with pytest.raises(ValueError):
        list(load_atlas(bad))
