"""Wilson score interval for a binomial proportion.

49/179 and 4900/17900 are both ~0.27 but not the same evidence; the Wilson interval
(better than the normal approximation at small n and near 0/1) shows how soft it is.
"""

from __future__ import annotations

import math


def wilson_ci(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score confidence interval for a binomial proportion (default 95%).

    Returns ``(low, high)`` clamped to ``[0, 1]``. For ``n == 0`` returns
    ``(0.0, 1.0)`` (total ignorance) rather than raising.
    """
    if n <= 0:
        return (0.0, 1.0)
    # Defensive: clamp successes into [0, n] so a caller passing an out-of-range count
    # can't drive p out of [0,1] and trigger a sqrt-of-negative (math domain error).
    successes = min(max(successes, 0), n)
    p = successes / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (max(0.0, center - half), min(1.0, center + half))
