"""Writers: results CSV (flat, factors visible) + provenance JSON (every number sourced)."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .score import TargetScore


def write_results_csv(scores: list[TargetScore], path: str | Path) -> Path:
    """One row per (target, population); all factors present as columns, never hidden."""
    rows = [row for ts in scores for row in ts.rows()]
    df = pd.DataFrame(rows)
    path = Path(path)
    df.to_csv(path, index=False)
    return path


def write_provenance_json(scores: list[TargetScore], path: str | Path) -> Path:
    """The audit trail: every factor's value + source + url + query_date + method."""
    payload = [ts.to_dict() for ts in scores]
    path = Path(path)
    path.write_text(json.dumps(payload, indent=2))
    return path
