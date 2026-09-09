"""Numerically stable primitives shared by the transfer equations.

Every burn magnitude in this project has the algebraic shape ``sqrt(1 + x) - 1``
for some ``x`` that tends to zero at a physically interesting boundary
(``R -> 1+`` or ``B -> R+``). Evaluating that difference directly loses roughly
half the available significant digits once ``|x|`` falls below about ``1e-8``.
Rewriting it in the equivalent, cancellation-free form removes the problem
exactly -- this is an algebraic identity, not an approximation.
"""

from __future__ import annotations

import math


def sqrt1pm1(x: float) -> float:
    """Return ``sqrt(1 + x) - 1`` without cancellation for small ``|x|``.

    Uses the identity ``sqrt(1+x) - 1 = x / (sqrt(1+x) + 1)``, whose denominator
    is near 2 rather than near 0, so no significant digits are lost. Exact at
    ``x = 0``, which is what makes the ``B = R`` bi-elliptic/Hohmann identity
    hold to the last bit rather than merely to a tolerance.

    Requires ``x >= -1``.
    """
    if x < -1.0:
        raise ValueError(f"sqrt1pm1 requires x >= -1, got {x!r}")
    return x / (math.sqrt(1.0 + x) + 1.0)
