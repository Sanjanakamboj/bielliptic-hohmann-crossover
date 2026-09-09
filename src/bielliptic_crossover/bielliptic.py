"""Three-impulse bi-elliptic transfer between coplanar circular orbits.

Normalization (see DESIGN.md Sections 1 and 3)::

    R = r2 / r1     final-to-initial RADIUS ratio,        R >= 1
    B = rb / r1     intermediate apoapsis RADIUS ratio,   B >= R

The manoeuvre is:

1. at ``r1``: circular -> ellipse 1 (periapsis ``r1``, apoapsis ``rb``) -- PROGRADE
2. at ``rb``: ellipse 1 -> ellipse 2 (periapsis ``r2``, apoapsis ``rb``) -- PROGRADE
3. at ``r2``: ellipse 2 -> circular                                     -- RETROGRADE

Burn 3 is retrograde because arrival at ``r2`` (the periapsis of ellipse 2) is
*faster* than circular whenever ``B > R``. This is the qualitative difference
from the Hohmann transfer, whose second burn is prograde. All three functions
below return positive **magnitudes**; direction is documented, never encoded as
a sign. Use :func:`burn_directions` to query direction explicitly.

At ``B = R`` the transfer degenerates exactly to the Hohmann transfer: burn 3
vanishes and burns 1 and 2 reduce to the Hohmann impulses. With the
cancellation-free formulation used here that identity holds to the last bit.
"""

from __future__ import annotations

import math

from ._stable import sqrt1pm1

__all__ = [
    "bielliptic_burns_normalized",
    "bielliptic_total_normalized",
    "bielliptic_infinite_limit_normalized",
    "bielliptic_total_derivative_wrt_B",
    "burn_directions",
]

#: Human-readable direction of each impulse, for ``B > R > 1``.
_BURN_DIRECTIONS = ("prograde", "prograde", "retrograde")


def burn_directions() -> tuple[str, str, str]:
    """Physical direction of each bi-elliptic impulse for ``B > R > 1``."""
    return _BURN_DIRECTIONS


def _validate_R_B(R: float, B: float) -> tuple[float, float]:
    R = float(R)
    B = float(B)
    if math.isnan(R) or math.isnan(B):
        raise ValueError("R and B must not be NaN")
    if R < 1.0:
        raise ValueError(f"R = r2/r1 must satisfy R >= 1, got {R!r}")
    if math.isinf(R):
        raise ValueError("R must be finite")
    if math.isinf(B):
        raise ValueError(
            "B must be finite; use bielliptic_infinite_limit_normalized(R) "
            "for the analytical B -> infinity limit"
        )
    if B < R:
        raise ValueError(f"B = rb/r1 must satisfy B >= R, got B={B!r}, R={R!r}")
    return R, B


def bielliptic_burns_normalized(R: float, B: float) -> tuple[float, float, float]:
    """Return the three normalized impulse MAGNITUDES ``(dv1_bar, dv2_bar, dv3_bar)``.

    All three are non-negative on the admissible domain ``B >= R >= 1``.
    ``dv3_bar`` is exactly ``0.0`` at ``B = R``.

    Directions (see :func:`burn_directions`): prograde, prograde, **retrograde**.
    """
    R, B = _validate_R_B(R, B)

    # dv1 = sqrt(2B/(1+B)) - 1,  with 2B/(1+B) = 1 + (B-1)/(1+B)
    dv1 = sqrt1pm1((B - 1.0) / (1.0 + B))

    # dv2 = sqrt(2R/(B(R+B))) - sqrt(2/(B(1+B)))
    #     = sqrt(2/(B(1+B))) * (sqrt(R(1+B)/(R+B)) - 1),
    #   with R(1+B)/(R+B) = 1 + B(R-1)/(R+B)
    dv2 = math.sqrt(2.0 / (B * (1.0 + B))) * sqrt1pm1(B * (R - 1.0) / (R + B))

    # dv3 = sqrt(2B/(R(R+B))) - 1/sqrt(R)
    #     = (1/sqrt(R)) * (sqrt(2B/(R+B)) - 1),  with 2B/(R+B) = 1 + (B-R)/(R+B)
    dv3 = sqrt1pm1((B - R) / (R + B)) / math.sqrt(R)

    return dv1, dv2, dv3


def bielliptic_total_normalized(R: float, B: float) -> float:
    """Total normalized bi-elliptic delta-v ``dv/v1`` (sum of the three magnitudes)."""
    dv1, dv2, dv3 = bielliptic_burns_normalized(R, B)
    return dv1 + dv2 + dv3


def bielliptic_infinite_limit_normalized(R: float) -> float:
    """Exact analytical ``B -> infinity`` limit ``(sqrt(2) - 1) * (1 + 1/sqrt(R))``.

    Derived term by term in DESIGN.md Section 4.2 -- **not** obtained by
    evaluating the finite-``B`` expression at some large ``B``.

    Physically: burn 1 becomes a parabolic escape increment from ``r1``, burn 2
    becomes free, and burn 3 becomes the mirror-image parabolic capture
    increment at ``r2``.

    This limit costs unbounded transfer time and is a mathematical asymptote,
    never a flight recommendation.
    """
    R = float(R)
    if math.isnan(R):
        raise ValueError("R must not be NaN")
    if R < 1.0:
        raise ValueError(f"R = r2/r1 must satisfy R >= 1, got {R!r}")
    if math.isinf(R):
        raise ValueError("R must be finite")
    return (math.sqrt(2.0) - 1.0) * (1.0 + 1.0 / math.sqrt(R))


def bielliptic_total_derivative_wrt_B(R: float, B: float) -> float:
    """Analytical ``d/dB`` of :func:`bielliptic_total_normalized` at fixed ``R``.

    Closed-form term-by-term derivative of the four ``B``-dependent terms
    (DESIGN.md Section 5.2). Valid on ``B >= R >= 1``; at ``B = R`` it gives the
    one-sided slope that defines threshold ``R2*``.

    Cross-checked in the test suite against complex-step differentiation, which
    shares none of this algebra.
    """
    R, B = _validate_R_B(R, B)

    sqrt2 = math.sqrt(2.0)

    # d/dB [ sqrt(2B/(1+B)) ]
    d1 = (1.0 / sqrt2) * B**-0.5 * (1.0 + B) ** -1.5

    # d/dB [ sqrt(2R/(B(R+B))) ]
    d2 = -(math.sqrt(2.0 * R) / 2.0) * B**-1.5 * (R + B) ** -1.5 * (R + 2.0 * B)

    # d/dB [ -sqrt(2/(B(1+B))) ]
    d3 = (sqrt2 / 2.0) * B**-1.5 * (1.0 + B) ** -1.5 * (1.0 + 2.0 * B)

    # d/dB [ sqrt(2B/(R(R+B))) ]
    d4 = 0.5 * math.sqrt(2.0 / R) * R * B**-0.5 * (R + B) ** -1.5

    return d1 + d2 + d3 + d4
