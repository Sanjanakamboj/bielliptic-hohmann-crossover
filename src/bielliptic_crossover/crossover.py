"""Crossover thresholds, break-even apoapsis ratio, and structural verification.

Two distinct radius-ratio thresholds govern this problem. They answer different
questions and are computed from **opposite ends** of the ``B`` domain
(DESIGN.md Section 5):

``R1*`` -- condition ``dv_bar_H(R) = dv_bar_B_inf(R)``, examined at ``B -> infinity``.
          Asks: can **any** bi-elliptic transfer beat Hohmann?  (exists B)

``R2*`` -- condition ``d/dB[dv_bar_B] = 0`` at ``B = R+``, examined at ``B -> R+``.
          Asks: do **all** bi-elliptic transfers beat Hohmann?  (for all B)

Nothing in this module hardcodes either value: both are recovered by root-finding
on the transfer equations, and each is independently cross-checked against an
exact polynomial characterization derived in M1.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import NamedTuple

import numpy as np
from scipy.optimize import brentq

from .bielliptic import (
    bielliptic_infinite_limit_normalized,
    bielliptic_total_derivative_wrt_B,
    bielliptic_total_normalized,
)
from .hohmann import hohmann_total_normalized

__all__ = [
    "threshold_R1",
    "threshold_R1_from_polynomial",
    "threshold_R2",
    "threshold_R2_from_cubic",
    "BreakEvenResult",
    "break_even_B",
    "classify_radius_ratio",
    "BStructure",
    "scan_B_structure",
    "B_RATIO_RESOLVABLE_MAX",
    "bielliptic_infimum_normalized",
    "REGIME_HOHMANN_ONLY",
    "REGIME_LARGE_B_BIELLIPTIC",
    "REGIME_ALL_BIELLIPTIC",
]

#: Region A -- Hohmann beats every admissible bi-elliptic transfer.
REGIME_HOHMANN_ONLY = "hohmann_only"
#: Region B -- only ``B > B_crit(R)`` beats Hohmann.
REGIME_LARGE_B_BIELLIPTIC = "large_B_bielliptic"
#: Region C -- every ``B > R`` beats Hohmann.
REGIME_ALL_BIELLIPTIC = "all_bielliptic"

# Largest B the break-even search will consider. Far beyond any physically
# meaningful apoapsis; bounded so that failure is reported rather than hidden.
_B_SEARCH_MAX = 1.0e30

# Relative offset used to step off the trivial B = R root without landing in
# its floating-point noise floor.
_B_PROBE_OFFSET = 1.0e-6


def _f1(R: float) -> float:
    """``dv_bar_H(R) - dv_bar_B_inf(R)``; its positive root is ``R1*``."""
    return hohmann_total_normalized(R) - bielliptic_infinite_limit_normalized(R)


@lru_cache(maxsize=1)
def threshold_R1() -> float:
    """Threshold ``R1*``, from the transfer equations.

    Root of ``dv_bar_H(R) = dv_bar_B_inf(R)``. Below it no bi-elliptic transfer
    of any ``B`` can beat Hohmann; above it, sufficiently distant ones can.

    Bracket: ``f1(1) = -2(sqrt(2)-1) < 0`` and ``f1`` is positive for large ``R``
    (it decays as ``(2-sqrt(2))/sqrt(R)``), so ``[1, 1e4]`` brackets the root.
    """
    return float(brentq(_f1, 1.0, 1.0e4, xtol=1e-15, rtol=8.9e-16, maxiter=200))


@lru_cache(maxsize=1)
def threshold_R1_from_polynomial() -> float:
    """Threshold ``R1*`` from the exact polynomial ``u^3 - (1+2*sqrt(2))*u^2 + u + 1 = 0``.

    With ``u = sqrt(R)``, squaring the reduced condition
    ``(R-1)/sqrt(R+1) = sqrt(R) + 1 - sqrt(2)`` yields this cubic (DESIGN.md
    Section 5.1). The physically relevant root is the unique real root with
    ``u > 1``. Shares no code path with :func:`threshold_R1`.
    """
    coefficients = [1.0, -(1.0 + 2.0 * math.sqrt(2.0)), 1.0, 1.0]
    roots = np.roots(coefficients)
    real_roots = [r.real for r in roots if abs(r.imag) < 1e-12 and r.real > 1.0]
    if len(real_roots) != 1:
        raise RuntimeError(f"expected exactly one root u > 1, got {real_roots!r}")
    return float(real_roots[0] ** 2)


@lru_cache(maxsize=1)
def threshold_R2() -> float:
    """Threshold ``R2*``, from the ``B -> R+`` slope condition.

    Root of ``d/dB[dv_bar_B(R,B)]|_(B=R) = 0``, using the analytical derivative
    in :func:`~bielliptic_crossover.bielliptic.bielliptic_total_derivative_wrt_B`.
    The slope is positive below ``R2*`` (a marginally bi-elliptic transfer is
    worse than Hohmann) and negative above it (it is immediately better).
    """

    def slope_at_B_equals_R(R: float) -> float:
        return bielliptic_total_derivative_wrt_B(R, R)

    return float(
        brentq(slope_at_B_equals_R, 1.0, 1.0e3, xtol=1e-15, rtol=8.9e-16, maxiter=200)
    )


@lru_cache(maxsize=1)
def threshold_R2_from_cubic() -> float:
    """Threshold ``R2*`` from the exact cubic ``R^3 - 15*R^2 - 9*R - 1 = 0``.

    Obtained by clearing radicals in the slope condition (DESIGN.md Section 5.2).
    The cubic has exactly one positive real root. Shares no code path with
    :func:`threshold_R2`.
    """
    roots = np.roots([1.0, -15.0, -9.0, -1.0])
    real_positive = [r.real for r in roots if abs(r.imag) < 1e-12 and r.real > 0.0]
    if len(real_positive) != 1:
        raise RuntimeError(f"expected exactly one positive root, got {real_positive!r}")
    return float(real_positive[0])


class BreakEvenResult(NamedTuple):
    """Outcome of a break-even apoapsis-ratio search at fixed ``R``.

    ``regime``
        One of :data:`REGIME_HOHMANN_ONLY`, :data:`REGIME_LARGE_B_BIELLIPTIC`,
        :data:`REGIME_ALL_BIELLIPTIC`.
    ``B_crit``
        The **nontrivial** break-even ratio, i.e. the unique ``B > R`` at which
        ``dv_bar_B(R,B) = dv_bar_H(R)``. Only meaningful in Region B; ``None``
        otherwise. The trivial root at ``B = R`` (where equality holds by
        construction) is never returned.
    ``winning_B_infimum``
        Infimum of the set of ``B`` that strictly beat Hohmann: ``B_crit`` in
        Region B, ``R`` in Region C (an open boundary, not attained), ``None``
        in Region A.
    """

    R: float
    regime: str
    B_crit: float | None
    winning_B_infimum: float | None
    description: str


def break_even_B(R: float) -> BreakEvenResult:
    """Find the nontrivial break-even ``B`` where bi-elliptic ties Hohmann.

    Solves ``dv_bar_B(R,B) - dv_bar_H(R) = 0`` for ``B > R``, explicitly
    stepping off the trivial root at ``B = R``.

    Region A returns ``regime = "hohmann_only"`` with ``B_crit = None``.
    Region C returns ``regime = "all_bielliptic"`` with ``B_crit = None`` and
    ``winning_B_infimum = R``: every ``B > R`` already wins, so the meaningful
    threshold collapses to the open ``B -> R+`` boundary and there is no finite
    break-even point to report. Manufacturing one there would be fiction.

    Raises ``RuntimeError`` if a root must exist by the infinite-``B`` limit but
    lies beyond the search ceiling (``R`` extremely close to ``R1*``).
    """
    R = float(R)
    if math.isnan(R) or math.isinf(R) or R <= 1.0:
        raise ValueError(f"break_even_B requires a finite R > 1, got {R!r}")

    dv_hohmann = hohmann_total_normalized(R)

    def excess(B: float) -> float:
        """Bi-elliptic minus Hohmann. Zero at B = R by construction."""
        return bielliptic_total_normalized(R, B) - dv_hohmann

    b_probe = R * (1.0 + _B_PROBE_OFFSET)
    if excess(b_probe) < 0.0:
        return BreakEvenResult(
            R=R,
            regime=REGIME_ALL_BIELLIPTIC,
            B_crit=None,
            winning_B_infimum=R,
            description=(
                "Region C: every admissible B > R beats Hohmann, so the winning "
                "set is the open interval (R, inf) and its infimum B = R is not "
                "attained. No finite nontrivial break-even B exists."
            ),
        )

    # Expand geometrically until the excess changes sign.
    b_lo = b_probe
    b_hi = b_probe
    while b_hi < _B_SEARCH_MAX:
        b_next = min(b_hi * 4.0, _B_SEARCH_MAX)
        if excess(b_next) < 0.0:
            b_lo, b_hi = b_hi, b_next
            break
        b_lo, b_hi = b_next, b_next
    else:  # pragma: no cover - loop always terminates via break or exhaustion
        pass

    if excess(b_hi) >= 0.0:
        # No sign change found. Either genuinely Region A, or R sits so close to
        # R1* that B_crit exceeds the search ceiling. Distinguish honestly.
        if bielliptic_infinite_limit_normalized(R) < dv_hohmann:
            raise RuntimeError(
                f"R = {R!r} lies in Region B (the infinite-B limit does beat "
                f"Hohmann) but the break-even B exceeds the search ceiling "
                f"{_B_SEARCH_MAX:g}. R is extremely close to R1*."
            )
        return BreakEvenResult(
            R=R,
            regime=REGIME_HOHMANN_ONLY,
            B_crit=None,
            winning_B_infimum=None,
            description=(
                "Region A: Hohmann beats every admissible bi-elliptic transfer. "
                "No break-even B exists."
            ),
        )

    # Solve in log-B: the root can span many orders of magnitude.
    log_root = brentq(
        lambda log_b: excess(math.exp(log_b)),
        math.log(b_lo),
        math.log(b_hi),
        xtol=1e-15,
        rtol=8.9e-16,
        maxiter=300,
    )
    b_crit = float(math.exp(log_root))

    if not b_crit > R:
        raise RuntimeError(
            f"break-even solver returned B = {b_crit!r}, which is not > R = {R!r}; "
            "this would be the trivial root"
        )

    return BreakEvenResult(
        R=R,
        regime=REGIME_LARGE_B_BIELLIPTIC,
        B_crit=b_crit,
        winning_B_infimum=b_crit,
        description=(
            "Region B: bi-elliptic beats Hohmann only for B > B_crit. Smaller "
            "apoapsis ratios near B = R are worse than Hohmann."
        ),
    )


def classify_radius_ratio(R: float, tol: float = 1e-12) -> str:
    """Classify ``R`` into the Region A / B / C structure.

    Returns :data:`REGIME_HOHMANN_ONLY`, :data:`REGIME_LARGE_B_BIELLIPTIC` or
    :data:`REGIME_ALL_BIELLIPTIC`, using the numerically recomputed thresholds
    rather than duplicated constants.

    ``tol`` is a relative tolerance applied at each threshold; a value within
    ``tol`` of a threshold is assigned to the lower region, so classification is
    well defined exactly at the boundary.
    """
    R = float(R)
    if math.isnan(R) or math.isinf(R) or R <= 1.0:
        raise ValueError(f"classify_radius_ratio requires a finite R > 1, got {R!r}")

    r1_star = threshold_R1()
    r2_star = threshold_R2()

    if R <= r1_star * (1.0 + tol):
        return REGIME_HOHMANN_ONLY
    if R <= r2_star * (1.0 + tol):
        return REGIME_LARGE_B_BIELLIPTIC
    return REGIME_ALL_BIELLIPTIC


def bielliptic_infimum_normalized(R: float) -> tuple[float, str]:
    """Infimum of ``dv_bar_B(R, .)`` over ``B > R``, and which endpoint attains it.

    M1 established -- and :func:`scan_B_structure` re-verifies numerically --
    that ``dv_bar_B(R, .)`` has no interior local minimum, so the infimum is
    always an endpoint value::

        min( dv_bar_H(R), dv_bar_B_inf(R) )

    Returns ``(value, endpoint)`` where ``endpoint`` is ``"B->R+"`` or
    ``"B->inf"``. Neither endpoint is attained on the open domain ``B > R``.
    """
    dv_hohmann = hohmann_total_normalized(R)
    dv_infinite = bielliptic_infinite_limit_normalized(R)
    if dv_hohmann <= dv_infinite:
        return dv_hohmann, "B->R+"
    return dv_infinite, "B->inf"


class BStructure(NamedTuple):
    """Result of a dense structural sweep of ``dv_bar_B(R, .)`` in ``B``."""

    R: float
    n_samples: int
    b_ratio_max: float
    n_interior_minima: int
    n_interior_maxima: int
    b_at_interior_maximum: float | None
    dv_at_B_near_R: float
    dv_at_B_max: float
    dv_hohmann: float
    dv_infinite_limit: float
    infimum_endpoint: str


#: Largest ``B/R`` at which the monotone large-``B`` tail is still resolvable in
#: double precision. Beyond roughly this ratio the residual
#: ``dv_bar_B - dv_bar_B_inf ~ (sqrt(R)-3)/(sqrt(2)B)`` falls below a few ulp of
#: ``dv_bar_B`` itself and the curve is flat to machine precision. Sweeping past
#: it does not reveal new structure, only rounding noise.
B_RATIO_RESOLVABLE_MAX = 1.0e10


def scan_B_structure(
    R: float,
    n_samples: int = 4001,
    b_ratio_max: float = B_RATIO_RESOLVABLE_MAX,
    noise_ulps: float = 8.0,
) -> BStructure:
    """Densely sweep ``B`` and report the interior stationary structure.

    Samples ``B/R`` logarithmically from just above 1 to ``b_ratio_max`` and
    counts interior local minima and maxima by three-point comparison.

    The expected -- and separately proved -- result is
    ``n_interior_minima == 0`` for every ``R``: the curve is monotone increasing,
    monotone decreasing, or up-then-down. A nonzero count means either a genuine
    discovery or, far more likely, a bug; it must be investigated, not accepted.

    ``noise_ulps`` sets how far a sample must sit below (or above) both
    neighbours to count as an extremum, measured in units in the last place. A
    turning point smaller than the rounding noise is not detectable, and
    counting one would report float noise as structure. Reducing it to ``0``
    recovers the raw comparison, which is dominated by noise once
    ``b_ratio_max`` exceeds :data:`B_RATIO_RESOLVABLE_MAX`.
    """
    R = float(R)
    if math.isnan(R) or math.isinf(R) or R <= 1.0:
        raise ValueError(f"scan_B_structure requires a finite R > 1, got {R!r}")
    if n_samples < 3:
        raise ValueError(f"n_samples must be at least 3, got {n_samples!r}")
    if b_ratio_max <= 1.0:
        raise ValueError(f"b_ratio_max must exceed 1, got {b_ratio_max!r}")

    ratios = np.logspace(math.log10(1.0 + 1e-10), math.log10(b_ratio_max), n_samples)
    b_values = R * ratios
    dv = np.array([bielliptic_total_normalized(R, float(b)) for b in b_values])

    left = dv[:-2]
    middle = dv[1:-1]
    right = dv[2:]
    # An extremum must clear the rounding noise floor to be real.
    tol = noise_ulps * np.spacing(np.abs(middle))
    minima_mask = (middle < left - tol) & (middle < right - tol)
    maxima_mask = (middle > left + tol) & (middle > right + tol)

    n_min = int(np.count_nonzero(minima_mask))
    n_max = int(np.count_nonzero(maxima_mask))

    b_at_max: float | None = None
    if n_max > 0:
        b_at_max = float(b_values[1:-1][maxima_mask][int(np.argmax(middle[maxima_mask]))])

    _, endpoint = bielliptic_infimum_normalized(R)

    return BStructure(
        R=R,
        n_samples=n_samples,
        b_ratio_max=b_ratio_max,
        n_interior_minima=n_min,
        n_interior_maxima=n_max,
        b_at_interior_maximum=b_at_max,
        dv_at_B_near_R=float(dv[0]),
        dv_at_B_max=float(dv[-1]),
        dv_hohmann=hohmann_total_normalized(R),
        dv_infinite_limit=bielliptic_infinite_limit_normalized(R),
        infimum_endpoint=endpoint,
    )
