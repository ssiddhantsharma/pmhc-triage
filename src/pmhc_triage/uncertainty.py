"""Small-sample uncertainty helpers.

A proportion estimated from 49/179 samples is not the same evidence as 4900/17900,
even though both are ~0.27. We attach a Wilson score interval (better than the
normal approximation at small n and near 0/1) so downstream consumers can see how
soft a fraction is.
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
    p = successes / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (max(0.0, center - half), min(1.0, center + half))
