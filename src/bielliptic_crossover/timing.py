"""Transfer times for the Hohmann and bi-elliptic manoeuvres.

Each transfer arc is half of an ellipse, so its duration is half the orbital
period of that ellipse::

    t_H = pi * sqrt(a_H^3 / mu),                a_H = (r1 + r2)/2

    t_B = pi * sqrt(a_1^3 / mu)
        + pi * sqrt(a_2^3 / mu),                a_1 = (r1 + rb)/2
                                                a_2 = (r2 + rb)/2

Nondimensional time scale::

    t_star = sqrt(r1^3 / mu)
    t_bar  = t / t_star

so that ``t_bar_H = pi * ((1+R)/2)^(3/2)`` and
``t_bar_B = pi * [((1+B)/2)^(3/2) + ((R+B)/2)^(3/2)]`` depend only on ``R``
and ``B``.

.. warning::
   **The ``B = R`` limit does NOT reduce to the Hohmann transfer time.**

   Delta-v and time degenerate differently at ``B = R``. The delta-v does
   reduce exactly (burn 3 vanishes), but ellipse 2 degenerates into the
   *circular target orbit*, and the bi-elliptic path still coasts a half
   revolution of it before reaching the nominal periapsis point. Hence

       t_bar_B(R, R) = t_bar_H(R) + pi * R^(3/2)

   i.e. the Hohmann time plus half the period of the final circular orbit.
   The extra term is the half-lap of loitering at ``r2`` that the Hohmann
   transfer never flies. See :func:`bielliptic_time_excess_at_B_equals_R`.

   This is a genuine property of the manoeuvre geometry, not a bug and not a
   discontinuity: ``t_bar_B`` is continuous in ``B`` down to ``B = R``; it
   simply does not converge to ``t_bar_H`` there.
"""

from __future__ import annotations

import math

__all__ = [
    "time_scale",
    "hohmann_transfer_time_normalized",
    "bielliptic_transfer_time_normalized",
    "hohmann_transfer_time",
    "bielliptic_transfer_time",
    "bielliptic_time_excess_at_B_equals_R",
]


def time_scale(mu: float, r1: float) -> float:
    """Nondimensional time scale ``t_star = sqrt(r1^3 / mu)``."""
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu!r}")
    if r1 <= 0.0:
        raise ValueError(f"r1 must be positive, got {r1!r}")
    return math.sqrt(r1**3 / mu)


def hohmann_transfer_time_normalized(R: float) -> float:
    """Hohmann transfer time in units of ``t_star``: ``pi * ((1+R)/2)^(3/2)``."""
    R = float(R)
    if math.isnan(R) or R < 1.0 or math.isinf(R):
        raise ValueError(f"R must be finite and satisfy R >= 1, got {R!r}")
    return math.pi * ((1.0 + R) / 2.0) ** 1.5


def bielliptic_transfer_time_normalized(R: float, B: float) -> float:
    """Bi-elliptic transfer time in units of ``t_star``.

    ``pi * [((1+B)/2)^(3/2) + ((R+B)/2)^(3/2)]``.

    Grows as ``B^(3/2)``, so ``B -> infinity`` implies unbounded duration.
    See the module docstring for why ``B = R`` does **not** give the Hohmann time.
    """
    R = float(R)
    B = float(B)
    if math.isnan(R) or math.isnan(B):
        raise ValueError("R and B must not be NaN")
    if R < 1.0 or math.isinf(R):
        raise ValueError(f"R must be finite and satisfy R >= 1, got {R!r}")
    if math.isinf(B):
        raise ValueError("B must be finite; B -> infinity implies unbounded transfer time")
    if B < R:
        raise ValueError(f"B must satisfy B >= R, got B={B!r}, R={R!r}")
    return math.pi * (((1.0 + B) / 2.0) ** 1.5 + ((R + B) / 2.0) ** 1.5)


def bielliptic_time_excess_at_B_equals_R(R: float) -> float:
    """Normalized time by which ``t_bar_B(R, R)`` exceeds ``t_bar_H(R)``.

    Equals ``pi * R^(3/2)`` -- half the period of the final circular orbit.
    Provided so the degeneracy documented in the module docstring is an
    explicit, testable quantity rather than a footnote.
    """
    R = float(R)
    if math.isnan(R) or R < 1.0 or math.isinf(R):
        raise ValueError(f"R must be finite and satisfy R >= 1, got {R!r}")
    return math.pi * R**1.5


def hohmann_transfer_time(mu: float, r1: float, r2: float) -> float:
    """Dimensional Hohmann transfer time, computed directly from ``a_H``."""
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu!r}")
    if r1 <= 0.0:
        raise ValueError(f"r1 must be positive, got {r1!r}")
    if r2 < r1:
        raise ValueError(f"r2 must satisfy r2 >= r1, got r2={r2!r}, r1={r1!r}")
    a_h = (r1 + r2) / 2.0
    return math.pi * math.sqrt(a_h**3 / mu)


def bielliptic_transfer_time(mu: float, r1: float, r2: float, rb: float) -> float:
    """Dimensional bi-elliptic transfer time, computed directly from ``a_1`` and ``a_2``."""
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu!r}")
    if r1 <= 0.0:
        raise ValueError(f"r1 must be positive, got {r1!r}")
    if r2 < r1:
        raise ValueError(f"r2 must satisfy r2 >= r1, got r2={r2!r}, r1={r1!r}")
    if rb < r2:
        raise ValueError(f"rb must satisfy rb >= r2, got rb={rb!r}, r2={r2!r}")
    a_1 = (r1 + rb) / 2.0
    a_2 = (r2 + rb) / 2.0
    return math.pi * (math.sqrt(a_1**3 / mu) + math.sqrt(a_2**3 / mu))
