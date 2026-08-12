import pytest

import pmhc_triage.http as _http


@pytest.fixture(autouse=True)
def _no_backoff_sleep(monkeypatch):
    """Silence retry backoff in tests: the retry logic still runs, it just doesn't sleep."""
    monkeypatch.setattr(_http, "_sleep", lambda *a, **k: None)
