"""Provenance primitive: every value the tool emits carries where it came from.

The whole point of ``pmhc-triage`` is *sourced, reproducible* numbers. A language
model will happily invent an incidence figure or an allele frequency; this module
makes that impossible here by forcing every number to travel with its source, URL,
query date, and derivation method -- and by *surfacing* (never silently dropping)
missing inputs.

Two types:

``Provenance``
    Where a single datum came from: a human-readable ``source``, an optional
    ``url``, the ISO ``query_date`` it was retrieved/computed, and the ``method``
    used to derive it.

``Sourced[T]``
    A value bound to its ``Provenance``, plus a list of ``warnings``. A value of
    ``None`` means *missing* -- the caller must be able to see that a factor could
    not be computed, rather than have a zero laundered into the final number.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import date, timezone, datetime
from typing import Any, Generic, TypeVar

T = TypeVar("T")


def today_iso() -> str:
    """UTC date as ``YYYY-MM-DD`` -- used to stamp ``query_date`` on retrievals."""
    return datetime.now(timezone.utc).date().isoformat()


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
    """A value bound to its provenance. ``value is None`` means the datum is missing.

    ``extra`` carries optional structured metadata about the value -- e.g. the
    sample size ``n`` and confidence interval behind a fraction -- so a number from
    n=23 is never silently trusted like one from n=1500.
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
