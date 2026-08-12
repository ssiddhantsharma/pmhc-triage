"""Retry/backoff behavior. Backoff sleeps are silenced by conftest, so these are fast."""

import httpx
import pytest

from pmhc_triage.http import _retry_after, request_with_retry


def _client(script):
    """Mock client that returns/raises each item of `script` on successive calls."""
    calls = {"n": 0}

    def handler(request):
        i = min(calls["n"], len(script) - 1)
        calls["n"] += 1
        item = script[i]
        if isinstance(item, Exception):
            raise item
        return httpx.Response(item, text="ok" if item == 200 else "err")

    return httpx.Client(transport=httpx.MockTransport(handler)), calls


def test_retries_transient_5xx_then_succeeds():
    client, calls = _client([503, 503, 200])
    resp = request_with_retry(client, "GET", "http://x/", retries=3)
    assert resp.status_code == 200
    assert calls["n"] == 3


def test_exhausts_retries_returns_final_error():
    client, calls = _client([500])  # always 500
    resp = request_with_retry(client, "GET", "http://x/", retries=2)
    assert resp.status_code == 500  # surfaced to caller, not swallowed
    assert calls["n"] == 3          # 1 attempt + 2 retries


def test_transport_error_retried_then_reraised():
    client, calls = _client([httpx.ConnectError("boom")])
    with pytest.raises(httpx.TransportError):
        request_with_retry(client, "GET", "http://x/", retries=2)
    assert calls["n"] == 3


def test_transport_error_then_success():
    client, calls = _client([httpx.ConnectError("boom"), 200])
    resp = request_with_retry(client, "GET", "http://x/", retries=3)
    assert resp.status_code == 200
    assert calls["n"] == 2


def test_non_retryable_status_returns_immediately():
    client, calls = _client([404])
    resp = request_with_retry(client, "GET", "http://x/", retries=3)
    assert resp.status_code == 404
    assert calls["n"] == 1  # 404 is a real answer, not retried


def test_retry_after_header_honored(monkeypatch):
    import pmhc_triage.http as h
    slept = []
    monkeypatch.setattr(h, "_sleep", lambda s: slept.append(s))
    state = {"first": True}

    def handler(request):
        if state["first"]:
            state["first"] = False
            return httpx.Response(429, headers={"Retry-After": "2"})
        return httpx.Response(200)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    resp = request_with_retry(client, "GET", "http://x/", retries=2)
    assert resp.status_code == 200
    assert slept == [2.0]  # used the header, not exponential backoff


def test_retry_after_parsing():
    def r(headers):
        return httpx.Response(429, headers=headers)
    assert _retry_after(r({"retry-after": "5"})) == 5.0
    assert _retry_after(r({"retry-after": "9999"})) == 30.0            # capped
    assert _retry_after(r({})) is None                                # absent
    assert _retry_after(r({"retry-after": "Wed, 21 Oct 2099 07:28:00 GMT"})) is None  # date form
