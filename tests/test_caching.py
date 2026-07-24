import httpx

from pmhc_triage.caching import CachingTransport


def test_miss_then_hit_no_second_inner_call(tmp_path):
    calls = {"n": 0}

    def handler(request):
        calls["n"] += 1
        return httpx.Response(200, json={"hello": "world"})

    t = CachingTransport(tmp_path, inner=httpx.MockTransport(handler))
    c = httpx.Client(transport=t)

    r1 = c.get("https://example.test/api")
    assert r1.json() == {"hello": "world"}
    assert r1.headers["x-pmhc-cache"] == "MISS"

    r2 = c.get("https://example.test/api")
    assert r2.headers["x-pmhc-cache"] == "HIT"
    assert calls["n"] == 1  # inner transport not hit again


def test_offline_replay_from_disk(tmp_path):
    def ok(request):
        return httpx.Response(200, json={"v": 42})

    httpx.Client(transport=CachingTransport(tmp_path, inner=httpx.MockTransport(ok))).get(
        "https://example.test/api"
    )

    # new transport, SAME cache dir, an inner that would fail if reached
    def boom(request):
        raise httpx.ConnectError("offline")

    c2 = httpx.Client(transport=CachingTransport(tmp_path, inner=httpx.MockTransport(boom)))
    r = c2.get("https://example.test/api")
    assert r.json() == {"v": 42}  # served from disk with no network
    assert r.headers["x-pmhc-cache"] == "HIT"
    assert r.headers["x-pmhc-fetched-at"]


def test_different_body_is_a_different_key(tmp_path):
    seen = []

    def handler(request):
        seen.append(request.content)
        return httpx.Response(200, json={"n": len(seen)})

    c = httpx.Client(transport=CachingTransport(tmp_path, inner=httpx.MockTransport(handler)))
    a = c.post("https://example.test/gql", json={"query": "A"})
    b = c.post("https://example.test/gql", json={"query": "B"})
    assert a.json() != b.json()  # distinct bodies -> distinct cache entries
    assert len(seen) == 2


def test_integration_second_run_is_fully_cached(tmp_path):
    """A whole cBioPortal-style flow: second run makes zero inner calls."""
    from pmhc_triage.antigen import variant_frequency

    calls = {"n": 0}

    def handler(request):
        calls["n"] += 1
        p = request.url.path
        if "/genes/" in p:
            return httpx.Response(200, json={"entrezGeneId": 3845})
        if "/sample-lists/" in p:
            return httpx.Response(200, json={"sampleIds": [f"S{i}" for i in range(179)]})
        return httpx.Response(200, json=[{"sampleId": f"S{i}", "proteinChange": "G12D"} for i in range(49)])

    inner = httpx.MockTransport(handler)
    c = httpx.Client(transport=CachingTransport(tmp_path, inner=inner))

    s1 = variant_frequency("paad_tcga_pan_can_atlas_2018", "G12D", client=c)
    first = calls["n"]
    s2 = variant_frequency("paad_tcga_pan_can_atlas_2018", "G12D", client=c)
    assert s1.value == s2.value
    assert calls["n"] == first  # second run entirely from cache


def test_provenance_query_date_is_fetch_time_not_rerun_time(tmp_path):
    """A cached run must report the ORIGINAL fetch datetime, not today (no lying)."""
    from pmhc_triage.antigen import variant_frequency

    def handler(request):
        p = request.url.path
        if "/genes/" in p:
            return httpx.Response(200, json={"entrezGeneId": 3845})
        if "/sample-lists/" in p:
            return httpx.Response(200, json={"sampleIds": [f"S{i}" for i in range(179)]})
        return httpx.Response(200, json=[{"sampleId": f"S{i}", "proteinChange": "G12D"} for i in range(49)])

    c = httpx.Client(transport=CachingTransport(tmp_path, inner=httpx.MockTransport(handler)))
    s1 = variant_frequency("paad_tcga_pan_can_atlas_2018", "G12D", client=c)  # MISS
    s2 = variant_frequency("paad_tcga_pan_can_atlas_2018", "G12D", client=c)  # HIT
    # both stamped from the cache's fetched_at (full ISO w/ 'T'), and identical
    assert "T" in s1.provenance.query_date  # a real fetch timestamp, not a bare re-run date
    assert s1.provenance.query_date == s2.provenance.query_date
