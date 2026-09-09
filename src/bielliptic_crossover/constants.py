"""Physical constants and the Earth reference case.

Deliberately contains **no crossover constants**. The radius-ratio thresholds
``R1*`` and ``R2*`` are numerically recomputed by :mod:`bielliptic_crossover.crossover`
from the transfer equations; hardcoding them here would defeat the verification.
"""

from __future__ import annotations

import math

#: Earth gravitational parameter [km^3/s^2].
MU_EARTH = 398600.4418

#: Earth mean equatorial radius [km].
R_EARTH = 6378.1363

#: Reference initial orbit ALTITUDE [km] (measured from the surface).
H1_REFERENCE = 300.0

#: Reference initial orbit RADIUS [km] (measured from the centre of the Earth).
R1_REFERENCE = R_EARTH + H1_REFERENCE

#: Reference initial circular speed [km/s].
V1_REFERENCE = math.sqrt(MU_EARTH / R1_REFERENCE)

#: Mars gravitational parameter [km^3/s^2], used only for scale-invariance tests.
MU_MARS = 42828.375214

#: Mars mean equatorial radius [km], used only for scale-invariance tests.
R_MARS = 3396.19


def altitude_from_radius(r: float, body_radius: float = R_EARTH) -> float:
    """Convert an orbital RADIUS to an ALTITUDE above the body's surface."""
    return r - body_radius


def radius_from_altitude(h: float, body_radius: float = R_EARTH) -> float:
    """Convert an ALTITUDE above the body's surface to an orbital RADIUS."""
    return h + body_radius
