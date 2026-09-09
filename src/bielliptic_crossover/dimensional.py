"""Independent dimensional transfer computations, direct from vis-viva.

This module is deliberately written **without reference to the normalized
closed forms** in :mod:`bielliptic_crossover.hohmann` and
:mod:`bielliptic_crossover.bielliptic`. It evaluates the actual speed at every
apsis from

    v^2 = mu * (2/r - 1/a)

and differences those speeds. It therefore shares no algebra with the
normalized path beyond the vis-viva equation itself, which makes agreement
between the two a genuine check against an algebra transcription bug rather
than a restatement of the same expression.

Speeds and signed velocity changes are exposed as well as magnitudes, so that
burn *direction* can be verified physically rather than assumed.
"""

from __future__ import annotations

import math
from typing import NamedTuple

__all__ = [
    "vis_viva_speed",
    "circular_speed",
    "HohmannBurnsDimensional",
    "BiellipticBurnsDimensional",
    "hohmann_burns_direct",
    "hohmann_total_dv_direct",
    "bielliptic_burns_direct",
    "bielliptic_total_dv",
]


def vis_viva_speed(mu: float, r: float, a: float) -> float:
    """Speed on a Keplerian orbit of semi-major axis ``a`` at radius ``r``."""
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu!r}")
    if r <= 0.0:
        raise ValueError(f"r must be positive, got {r!r}")
    if a <= 0.0:
        raise ValueError(f"a must be positive (bound orbits only), got {a!r}")
    v_squared = mu * (2.0 / r - 1.0 / a)
    if v_squared < 0.0:
        raise ValueError(
            f"vis-viva gave v^2 = {v_squared!r} < 0: radius r={r!r} is outside "
            f"the orbit with a={a!r}"
        )
    return math.sqrt(v_squared)


def circular_speed(mu: float, r: float) -> float:
    """Circular orbit speed ``sqrt(mu/r)``."""
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu!r}")
    if r <= 0.0:
        raise ValueError(f"r must be positive, got {r!r}")
    return math.sqrt(mu / r)


class HohmannBurnsDimensional(NamedTuple):
    """Dimensional Hohmann burn breakdown. ``dv*_signed`` is ``v_after - v_before``."""

    dv1: float
    dv2: float
    total: float
    dv1_signed: float
    dv2_signed: float
    v_circ_1: float
    v_peri_transfer: float
    v_apo_transfer: float
    v_circ_2: float


class BiellipticBurnsDimensional(NamedTuple):
    """Dimensional bi-elliptic burn breakdown. ``dv*_signed`` is ``v_after - v_before``."""

    dv1: float
    dv2: float
    dv3: float
    total: float
    dv1_signed: float
    dv2_signed: float
    dv3_signed: float
    v_circ_1: float
    v_peri_ellipse1: float
    v_apo_ellipse1: float
    v_apo_ellipse2: float
    v_peri_ellipse2: float
    v_circ_2: float


def _validate_radii(mu: float, r1: float, r2: float) -> None:
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu!r}")
    if r1 <= 0.0:
        raise ValueError(f"r1 must be positive, got {r1!r}")
    if r2 < r1:
        raise ValueError(f"r2 must satisfy r2 >= r1, got r2={r2!r}, r1={r1!r}")


def hohmann_burns_direct(mu: float, r1: float, r2: float) -> HohmannBurnsDimensional:
    """Hohmann burns computed directly from vis-viva speeds at each apsis."""
    _validate_radii(mu, r1, r2)
    a_h = (r1 + r2) / 2.0

    v_circ_1 = circular_speed(mu, r1)
    v_peri = vis_viva_speed(mu, r1, a_h)
    v_apo = vis_viva_speed(mu, r2, a_h)
    v_circ_2 = circular_speed(mu, r2)

    dv1_signed = v_peri - v_circ_1  # prograde: speed up onto the transfer ellipse
    dv2_signed = v_circ_2 - v_apo  # prograde: speed up to circularize

    return HohmannBurnsDimensional(
        dv1=abs(dv1_signed),
        dv2=abs(dv2_signed),
        total=abs(dv1_signed) + abs(dv2_signed),
        dv1_signed=dv1_signed,
        dv2_signed=dv2_signed,
        v_circ_1=v_circ_1,
        v_peri_transfer=v_peri,
        v_apo_transfer=v_apo,
        v_circ_2=v_circ_2,
    )


def hohmann_total_dv_direct(mu: float, r1: float, r2: float) -> float:
    """Total dimensional Hohmann delta-v, direct vis-viva path."""
    return hohmann_burns_direct(mu, r1, r2).total


def bielliptic_burns_direct(
    mu: float, r1: float, r2: float, rb: float
) -> BiellipticBurnsDimensional:
    """Bi-elliptic burns computed directly from vis-viva speeds at each apsis.

    Note the sign structure, which is checked in the test suite rather than
    assumed: ``dv1_signed > 0`` and ``dv2_signed > 0`` (prograde), while
    ``dv3_signed < 0`` for ``rb > r2`` -- the spacecraft arrives at ``r2``
    faster than circular and must brake.
    """
    _validate_radii(mu, r1, r2)
    if rb < r2:
        raise ValueError(f"rb must satisfy rb >= r2, got rb={rb!r}, r2={r2!r}")

    a_1 = (r1 + rb) / 2.0
    a_2 = (r2 + rb) / 2.0

    v_circ_1 = circular_speed(mu, r1)
    v_peri_1 = vis_viva_speed(mu, r1, a_1)
    v_apo_1 = vis_viva_speed(mu, rb, a_1)
    v_apo_2 = vis_viva_speed(mu, rb, a_2)
    v_peri_2 = vis_viva_speed(mu, r2, a_2)
    v_circ_2 = circular_speed(mu, r2)

    dv1_signed = v_peri_1 - v_circ_1  # prograde
    dv2_signed = v_apo_2 - v_apo_1  # prograde
    dv3_signed = v_circ_2 - v_peri_2  # RETROGRADE (negative for rb > r2)

    total = abs(dv1_signed) + abs(dv2_signed) + abs(dv3_signed)

    return BiellipticBurnsDimensional(
        dv1=abs(dv1_signed),
        dv2=abs(dv2_signed),
        dv3=abs(dv3_signed),
        total=total,
        dv1_signed=dv1_signed,
        dv2_signed=dv2_signed,
        dv3_signed=dv3_signed,
        v_circ_1=v_circ_1,
        v_peri_ellipse1=v_peri_1,
        v_apo_ellipse1=v_apo_1,
        v_apo_ellipse2=v_apo_2,
        v_peri_ellipse2=v_peri_2,
        v_circ_2=v_circ_2,
    )


def bielliptic_total_dv(mu: float, r1: float, r2: float, rb: float) -> float:
    """Total dimensional bi-elliptic delta-v, direct vis-viva path."""
    return bielliptic_burns_direct(mu, r1, r2, rb).total
