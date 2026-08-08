"""Validation: does the tool's output agree with facts we already know to be true?

These are LIVE checks (they hit cBioPortal / UniProt / MHCflurry), so they are
marked ``live`` and EXCLUDED from the default + CI run (which stays hermetic). Run
them explicitly:

    pytest -m live -v

Each self-skips if its service (or the optional MHCflurry extra) is unavailable, so
running them never hard-fails on a network blip -- it either confirms the fact or
skips. Ground truth is coarse (direction / range, not exact numbers) so it stays
robust as the underlying databases update.
"""

from __future__ import annotations

import importlib.util

import pytest

pytestmark = pytest.mark.live


def test_kras_g12d_antigen_fraction_matches_literature():
    # KRAS G12D is ~27% of PDAC tumors (canonical). Assert the live cBioPortal
    # count lands in a defensible range with a real denominator.
    from pmhc_triage.antigen import variant_frequency
    try:
        s = variant_frequency("paad_tcga_pan_can_atlas_2018", "G12D", gene="KRAS")
    except Exception as exc:  # network / transport
        pytest.skip(f"cBioPortal unreachable: {exc}")
    if s.is_missing:
        pytest.skip(f"cBioPortal returned no value: {s.warnings}")
    assert 0.20 <= s.value <= 0.35, f"G12D fraction {s.value} outside expected ~0.27"
    assert (s.extra or {}).get("denominator", 0) > 100


def test_paad_study_resolves_with_sequenced_samples():
    from pmhc_triage.antigen import check_study
    try:
        s = check_study("paad_tcga_pan_can_atlas_2018")
    except Exception as exc:
        pytest.skip(f"cBioPortal unreachable: {exc}")
    if s.is_missing:
        pytest.skip("cBioPortal unavailable")
    assert s.value is True


def test_uniprot_kras_sequence_has_g12():
    # The load-bearing WT guard depends on the real sequence; verify residue 12 is G.
    from pmhc_triage.sequences import fetch_uniprot_sequence
    try:
        s = fetch_uniprot_sequence("P01116")
    except Exception as exc:
        pytest.skip(f"UniProt unreachable: {exc}")
    if s.is_missing:
        pytest.skip("UniProt unavailable")
    assert s.value[11] == "G"  # 1-indexed position 12


@pytest.mark.skipif(importlib.util.find_spec("mhcflurry") is None,
                    reason="mhcflurry ([presentation] extra) not installed")
def test_kras_g12d_recovers_canonical_a3_restriction():
    # The tool should recover the textbook A3-superfamily restriction from sequence.
    from pmhc_triage.peptides import mutant_peptides
    from pmhc_triage.presentation import predict_presenting_alleles
    from pmhc_triage.sequences import fetch_uniprot_sequence
    try:
        seq = fetch_uniprot_sequence("P01116")
    except Exception as exc:
        pytest.skip(f"UniProt unreachable: {exc}")
    if seq.is_missing:
        pytest.skip("UniProt unavailable")
    peps = mutant_peptides(seq.value, "G12D").value
    s = predict_presenting_alleles(peps, ["A*02:01", "A*03:01", "A*11:01"],
                                   threshold_percentile=0.5)
    if s.is_missing:
        pytest.skip(f"MHCflurry prediction unavailable: {s.warnings}")
    assert "A*03:01" in s.value and "A*11:01" in s.value  # canonical KRAS G12D restriction
