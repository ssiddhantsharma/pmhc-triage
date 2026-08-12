"""Send an HTTP request, retrying transient failures with exponential backoff.

External services (cBioPortal, AFND, UniProt, Open Targets) blip. A single-attempt call
loses a factor to a momentary 429/5xx or timeout and forces a re-run/re-fetch, so this
retries transient failures (honoring Retry-After) then gives up and lets the caller
surface it as missing. Non-transient responses (2xx, 4xx) return immediately; a
non-transient exception propagates unchanged.
"""

from __future__ import annotations

import time

import httpx

_RETRYABLE_STATUS = frozenset({429, 500, 502, 503, 504})
_MAX_RETRY_AFTER = 30.0  # cap a Retry-After header so a hostile value can't hang a run

# Module-level indirection so tests can silence backoff without real sleeps.
_sleep = time.sleep


def request_with_retry(
    client: httpx.Client,
    method: str,
    url: str,
    *,
    retries: int = 3,
    base_delay: float = 0.5,
    max_delay: float = 8.0,
    timeout: float | None = None,
    **kwargs,
) -> httpx.Response:
    """client.request(method, url, **kwargs) with retry on 429/5xx and transport errors.

    Retries up to ``retries`` extra times with exponential backoff (capped at
    ``max_delay``, honoring a numeric Retry-After). Returns the final response (the
    caller checks status as usual); re-raises the last transport error if every attempt
    fails on a connection/timeout.
    """
    if timeout is not None:
        kwargs.setdefault("timeout", timeout)
    for attempt in range(retries + 1):
        try:
            resp = client.request(method, url, **kwargs)
        except httpx.TransportError:
            if attempt >= retries:
                raise
            _sleep(_backoff(attempt, base_delay, max_delay))
            continue
        if resp.status_code in _RETRYABLE_STATUS and attempt < retries:
            _sleep(_retry_after(resp) or _backoff(attempt, base_delay, max_delay))
            continue
        return resp


def _backoff(attempt: int, base: float, cap: float) -> float:
    return min(base * (2 ** attempt), cap)


def _retry_after(resp: httpx.Response) -> float | None:
    """Seconds from a numeric Retry-After header, capped; None if absent or a date form."""
    v = resp.headers.get("retry-after")
    if not v:
        return None
    try:
        return min(max(float(v), 0.0), _MAX_RETRY_AFTER)
    except ValueError:
        return None
