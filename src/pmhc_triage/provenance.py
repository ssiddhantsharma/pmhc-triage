"""Provenance primitives: every emitted value carries where it came from.

Provenance: source, optional url, query_date, derivation method.
Sourced[T]: a value bound to its Provenance plus warnings; value None means missing
(a factor that could not be computed), never a laundered zero.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any, Generic, TypeVar

T = TypeVar("T")


def today_iso() -> str:
    """UTC date as YYYY-MM-DD."""
    return datetime.now(timezone.utc).date().isoformat()


# Header the caching transport stamps with the original fetch datetime.
CACHE_FETCHED_AT_HEADER = "x-pmhc-fetched-at"


def fetched_at_or_today(response) -> str:
    """Original fetch datetime if the response came from cache, else today's date."""
    try:
        value = response.headers.get(CACHE_FETCHED_AT_HEADER)
    except Exception:
        value = None
    return value or today_iso()


@dataclass(frozen=True)
class Provenance:
    """Where one datum came from. ``source`` is required; the rest are optional."""

    source: str
    url: str | None = None
    query_date: str | None = None
    method: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {k: v for k, v in asdict(self).items() if v is not None}


@dataclass
class Sourced(Generic[T]):
    """A value bound to its provenance; value None means missing.

    extra carries optional structured metadata (e.g. sample size and CI behind a
    fraction) so a number from n=23 is not trusted like one from n=1500.
    """

    value: T | None
    provenance: Provenance
    warnings: list[str] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def is_missing(self) -> bool:
        return self.value is None

    def warn(self, message: str) -> "Sourced[T]":
        """Attach a warning and return ``self`` (chainable). Warnings are never fatal."""
        self.warnings.append(message)
        return self

    def to_dict(self) -> dict[str, Any]:
        d = {
            "value": self.value,
            "missing": self.is_missing,
            "warnings": list(self.warnings),
            "provenance": self.provenance.to_dict(),
        }
        if self.extra:
            d["extra"] = dict(self.extra)
        return d
