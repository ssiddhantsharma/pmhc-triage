from pmhc_triage.provenance import Provenance, Sourced, today_iso


def test_provenance_to_dict_drops_none_fields():
    p = Provenance(source="GLOBOCAN 2022")
    assert p.to_dict() == {"source": "GLOBOCAN 2022"}
    p2 = Provenance(source="cBioPortal", url="https://cbioportal.org", method="mutation count")
    d = p2.to_dict()
    assert d["source"] == "cBioPortal" and d["url"] == "https://cbioportal.org"
    assert "query_date" not in d  # None fields omitted


def test_sourced_missing_and_warnings():
    s = Sourced(value=None, provenance=Provenance(source="x"))
    assert s.is_missing
    s.warn("could not retrieve")
    assert s.warnings == ["could not retrieve"]
    assert s.to_dict()["missing"] is True


def test_sourced_present_value():
    s = Sourced(value=0.42, provenance=Provenance(source="AFND"))
    assert not s.is_missing
    d = s.to_dict()
    assert d["value"] == 0.42 and d["missing"] is False


def test_today_iso_format():
    t = today_iso()
    assert len(t) == 10 and t[4] == "-" and t[7] == "-"
