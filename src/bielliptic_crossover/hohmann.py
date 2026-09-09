"""Two-impulse Hohmann transfer between coplanar circular orbits.

Normalization (see DESIGN.md Sections 1-2)::

    R      = r2 / r1            final-to-initial RADIUS ratio, R >= 1
    v1     = sqrt(mu / r1)      initial circular speed
    dv_bar = dv / v1            normalized delta-v

Both impulses are prograde. The closed forms are

    dv_bar_H1 = sqrt(2R/(1+R)) - 1
    dv_bar_H2 = 1/sqrt(R) - sqrt(2/(R*(1+R)))

evaluated here in the algebraically equivalent, cancellation-free ``sqrt1pm1``
form so that the ``R -> 1+`` limit stays accurate to full precision.
"""

from __future__ import annotations

import math

from ._stable import sqrt1pm1

__all__ = [
    "hohmann_burns_normalized",
    "hohmann_total_normalized",
    "hohmann_total_dv",
]


def _validate_R(R: float) -> float:
    R = float(R)
    if math.isnan(R):
        raise ValueError("R must not be NaN")
    if R < 1.0:
        raise ValueError(f"R = r2/r1 must satisfy R >= 1, got {R!r}")
    if math.isinf(R):
        raise ValueError("R must be finite")
    return R


def hohmann_burns_normalized(R: float) -> tuple[float, float]:
    """Return the two normalized Hohmann impulse magnitudes ``(dv1_bar, dv2_bar)``.

    Both are positive for ``R > 1`` and both burns are **prograde**. At ``R = 1``
    both are exactly ``0.0``.

    Raises ``ValueError`` for ``R < 1``, NaN or infinite ``R``.
    """
    R = _validate_R(R)

    # dv1: sqrt(2R/(1+R)) - 1, with 2R/(1+R) = 1 + (R-1)/(1+R).
    dv1 = sqrt1pm1((R - 1.0) / (1.0 + R))

    # dv2: (1/sqrt(R)) * (1 - sqrt(2/(1+R))), with 2/(1+R) = 1 + (1-R)/(1+R).
    dv2 = -sqrt1pm1((1.0 - R) / (1.0 + R)) / math.sqrt(R)

    return dv1, dv2


def hohmann_total_normalized(R: float) -> float:
    """Total normalized Hohmann delta-v ``dv/v1``. Exactly ``0.0`` at ``R = 1``."""
    dv1, dv2 = hohmann_burns_normalized(R)
    return dv1 + dv2


def hohmann_total_dv(mu: float, r1: float, r2: float) -> float:
    """Total dimensional Hohmann delta-v [same length/time units as inputs].

    Computed by rescaling the normalized result. For an implementation that
    instead evaluates every apsis speed directly from vis-viva, see
    :func:`bielliptic_crossover.dimensional.hohmann_total_dv_direct`; the two
    agree to floating-point tolerance and are cross-checked in the test suite.
    """
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu!r}")
    if r1 <= 0.0:
        raise ValueError(f"r1 must be positive, got {r1!r}")
    if r2 < r1:
        raise ValueError(f"r2 must satisfy r2 >= r1, got r2={r2!r}, r1={r1!r}")
    return hohmann_total_normalized(r2 / r1) * math.sqrt(mu / r1)
