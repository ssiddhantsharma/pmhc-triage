import pytest

from pmhc_triage.peptides import mutant_peptides, parse_substitution

# Real KRAS N-terminus (UniProt P01116, verified): residue 12 (1-indexed) is G.
KRAS_NTERM = "MTEYKLVVVGAGGVGKSALTIQLIQ"


def test_parse_substitution():
    assert parse_substitution("G12D") == ("G", 12, "D")
    assert parse_substitution(" A146T ") == ("A", 146, "T")


@pytest.mark.parametrize("bad", ["G12", "12D", "GD", "g12d", "G0D", ""])
def test_parse_substitution_rejects_malformed(bad):
    with pytest.raises(ValueError):
        parse_substitution(bad)


def test_g12d_peptides_span_mutation_and_have_right_lengths():
    s = mutant_peptides(KRAS_NTERM, "G12D")
    assert not s.is_missing and s.warnings == []
    peps = s.value
    assert peps, "expected peptides"
    # every peptide is one of the requested lengths
    assert all(len(p) in (8, 9, 10, 11) for p in peps)
    # the canonical G12D 9-mer spanning the mutation is present
    assert "VVGADGVGK" in peps
    # the 10-mer neoepitope too
    assert "VVVGADGVGK" in peps
    # no duplicates
    assert len(peps) == len(set(peps))


def test_wt_residue_mismatch_refuses():
    # sequence has G at position 12, but we claim A12D -> must refuse, not mutate blindly
    s = mutant_peptides(KRAS_NTERM, "A12D")
    assert s.is_missing
    assert any("WT residue mismatch" in w for w in s.warnings)


def test_position_out_of_range_surfaces():
    s = mutant_peptides(KRAS_NTERM, "G999D")
    assert s.is_missing
    assert any("beyond sequence length" in w for w in s.warnings)


def test_malformed_variant_surfaces_not_raises():
    s = mutant_peptides(KRAS_NTERM, "not-a-variant")
    assert s.is_missing and s.warnings


def test_provenance_records_method():
    s = mutant_peptides(KRAS_NTERM, "G12D")
    d = s.to_dict()
    assert "G12D" in d["provenance"]["method"]
    assert d["provenance"]["query_date"]


def test_synonymous_substitution_rejected():
    import pytest
    from pmhc_triage.peptides import parse_substitution, mutant_peptides
    with pytest.raises(ValueError, match="does not change"):
        parse_substitution("G12G")
    r = mutant_peptides("MTEYKLVVVGAGGVGKSALTIQ", "G12G")
    assert r.is_missing
