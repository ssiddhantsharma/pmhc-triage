"""On-disk HTTP cache for reproducible, offline re-runs.

Every remote factor is fetched through httpx clients that accept an injected
transport. ``cached_client(dir)`` returns a client whose responses are cached to
disk keyed by (method, URL, request body), each stamped with the UTC datetime it
was fetched. A re-run served from cache returns byte-identical responses -- so a
provenance run reproduces exactly, and works with the network unplugged.

Cache entries carry ``x-pmhc-cache: HIT|MISS`` and ``x-pmhc-fetched-at`` response
headers so a caller can tell fresh from replayed data.
"""

from __future__ import annotations

import base64
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import httpx


def _fetched_at() -> str:
    return datetime.now(timezone.utc).isoformat()


class CachingTransport(httpx.BaseTransport):
    """httpx transport that caches responses to ``cache_dir`` and replays them."""

    def __init__(self, cache_dir: str | Path, inner: httpx.BaseTransport | None = None):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.inner = inner or httpx.HTTPTransport()

    def _key(self, request: httpx.Request) -> str:
        h = hashlib.sha256()
        h.update(request.method.encode())
        h.update(b"\0")
        h.update(str(request.url).encode())
        h.update(b"\0")
        h.update(request.content or b"")
        return h.hexdigest()

    # The inner transport's .read() returns ALREADY-DECODED bytes, but leaves the
    # content-encoding header in place. So we must drop content-encoding (else the
    # client tries to gunzip plain data -> "incorrect header check") and
    # content-length (httpx recomputes it from the returned content).
    _DROP_HEADERS = {"content-encoding", "content-length", "transfer-encoding"}

    @classmethod
    def _replay_headers(cls, stored: dict) -> dict:
        return {k: v for k, v in stored.items() if k.lower() not in cls._DROP_HEADERS}

    def handle_request(self, request: httpx.Request) -> httpx.Response:
        path = self.cache_dir / f"{self._key(request)}.json"
        if path.exists():
            rec = json.loads(path.read_text())
            content = base64.b64decode(rec["content_b64"])
            headers = self._replay_headers(rec.get("headers", {}))
            headers.update({"x-pmhc-cache": "HIT", "x-pmhc-fetched-at": rec["fetched_at"]})
            return httpx.Response(rec["status"], content=content, headers=headers, request=request)

        resp = self.inner.handle_request(request)
        body = resp.read()
        fetched_at = _fetched_at()
        stored_headers = self._replay_headers(dict(resp.headers))
        rec = {
            "method": request.method,
            "url": str(request.url),
            "status": resp.status_code,
            "content_b64": base64.b64encode(body).decode(),
            "headers": stored_headers,
            "fetched_at": fetched_at,
        }
        path.write_text(json.dumps(rec))
        headers = {**stored_headers, "x-pmhc-cache": "MISS", "x-pmhc-fetched-at": fetched_at}
        return httpx.Response(resp.status_code, content=body, headers=headers, request=request)


def cached_client(cache_dir: str | Path, *, inner: httpx.BaseTransport | None = None, **kwargs) -> httpx.Client:
    """An ``httpx.Client`` backed by an on-disk cache at ``cache_dir``."""
    return httpx.Client(transport=CachingTransport(cache_dir, inner=inner), **kwargs)
