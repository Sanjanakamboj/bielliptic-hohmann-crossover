"""Milestone 3 trade analysis: crossover map, break-even locus, delta-v/time trade.

This module composes the verified M2 primitives into the engineering trade study.
It contains **no transfer equations of its own** -- every delta-v and every
transfer time is obtained by calling the M2 production functions, so there is a
single source of truth for the physics.

Three ideas are kept rigorously separate throughout, because conflating them is
the standard error in this problem:

1. **Mathematical delta-v dominance** -- which family contains a cheaper transfer
   at all. Decided by the two thresholds ``R1*`` and ``R2*``.
2. **Finite-``B`` practical saving** -- what a real, bounded apoapsis actually
   buys. Always smaller than the mathematical bound, sometimes negative.
3. **Transfer-time and model-validity penalty** -- what that saving costs in
   duration, and whether the two-body model still applies at the required
   apoapsis.

The ``B -> infinity`` branch is referred to throughout as the *mathematical
infimum*, never as an optimum or a recommendation: it costs unbounded transfer
time and is not a realizable transfer.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import NamedTuple

import numpy as np

from .bielliptic import (
    bielliptic_infinite_limit_normalized,
    bielliptic_total_normalized,
)
from .constants import MU_EARTH, R1_REFERENCE, R_EARTH
from .crossover import (
    REGIME_ALL_BIELLIPTIC,
    REGIME_HOHMANN_ONLY,
    REGIME_LARGE_B_BIELLIPTIC,
    bielliptic_infimum_normalized,
    break_even_B,
    classify_radius_ratio,
    threshold_R1,
    threshold_R2,
)
from .hohmann import hohmann_total_normalized
from .timing import bielliptic_transfer_time, hohmann_transfer_time

__all__ = [
    "LUNAR_DISTANCE_KM",
    "EARTH_HILL_RADIUS_KM",
    "MODEL_OK",
    "MODEL_CAUTION",
    "MODEL_INVALID",
    "GEO_RADIUS_KM",
    "FINITE_B_OVER_R_FAMILY",
    "TradePoint",
    "RadiusRatioTrade",
    "EndpointInfimum",
    "model_validity",
    "minimum_endpoint_transfer",
    "practical_trade_metrics",
    "break_even_trade",
    "evaluate_trade",
    "authoritative_R_grid",
    "break_even_locus_grid",
    "break_even_locus",
    "sweep_radius_ratios",
    "recommendation_for",
    "recommendation_summary",
]

#: Mean Earth-Moon distance [km]. Used only to flag where an Earth-centred
#: two-body model stops being physically appropriate.
LUNAR_DISTANCE_KM = 384400.0

#: Approximate Earth Hill-sphere radius [km]. Beyond this the Earth-centred
#: two-body model is not merely inaccurate, it is the wrong model.
EARTH_HILL_RADIUS_KM = 1.5e6

#: Geostationary orbit RADIUS [km].
GEO_RADIUS_KM = 42164.0

MODEL_OK = "two_body_reasonable"
MODEL_CAUTION = "lunar_third_body_caution"
MODEL_INVALID = "two_body_invalid"

#: Representative finite apoapsis strategies. These are illustrative bounded
#: designs, NOT optima -- M1/M2 established that no finite interior optimum
#: exists, so no value here can be called optimal.
FINITE_B_OVER_R_FAMILY = (1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0)

_SECONDS_PER_DAY = 86400.0


def model_validity(rb_km: float) -> tuple[str, str]:
    """Classify whether an Earth-centred two-body model still applies at ``rb_km``.

    Returns ``(flag, note)``. The thresholds are the lunar distance (beyond which
    lunar third-body acceleration is not a perturbation but a dominant term) and
    the Earth Hill radius (beyond which the spacecraft is not bound to Earth in
    any useful sense).
    """
    if rb_km >= EARTH_HILL_RADIUS_KM:
        return (
            MODEL_INVALID,
            "apoapsis exceeds the Earth Hill radius; the spacecraft is not "
            "meaningfully Earth-bound and this two-body result is a mathematical "
            "extrapolation only",
        )
    if rb_km >= LUNAR_DISTANCE_KM:
        return (
            MODEL_INVALID,
            "apoapsis exceeds the lunar distance; an Earth-only two-body model "
            "is no longer physically appropriate here",
        )
    if rb_km >= 0.1 * LUNAR_DISTANCE_KM:
        return (
            MODEL_CAUTION,
            "apoapsis exceeds 10 percent of the lunar distance; lunar and solar "
            "third-body effects would be significant in reality",
        )
    return (MODEL_OK, "")


class EndpointInfimum(NamedTuple):
    """The best delta-v achievable over ``B``, which is always an endpoint value.

    M1/M2 established that ``dv_bar_B(R, .)`` has no finite interior minimum, so

        inf over B of dv_bar_B(R,B) = min( dv_bar_H(R), dv_bar_B_inf(R) )

    ``endpoint`` is ``"B->R+"`` (i.e. the Hohmann transfer) or ``"B->inf"``.
    Neither endpoint is attained on the open domain ``B > R``: the ``B->inf``
    branch in particular is a **mathematical infimum with unbounded transfer
    time**, never a realizable optimum.
    """

    R: float
    dv_hohmann_bar: float
    dv_infinite_bar: float
    dv_best_bar: float
    best_endpoint: str
    saving_bar: float
    saving_m_s: float
    attainable: bool


class TradePoint(NamedTuple):
    """Full delta-v / time / geometry metrics for one ``(R, B)`` pair."""

    R: float
    B: float
    B_over_R: float
    label: str
    dv_hohmann_bar: float
    dv_bielliptic_bar: float
    saving_bar: float
    dv_hohmann_m_s: float
    dv_bielliptic_m_s: float
    dv_saving_m_s: float
    t_hohmann_days: float
    t_bielliptic_days: float
    time_ratio: float
    extra_time_days: float
    rb_km: float
    rb_altitude_km: float
    rb_over_lunar_distance: float
    model_validity_flag: str
    model_validity_note: str


class RadiusRatioTrade(NamedTuple):
    """Everything M3 computes at a single radius ratio ``R``."""

    R: float
    regime: str
    r2_km: float
    r2_altitude_km: float
    infimum: EndpointInfimum
    B_crit: float | None
    points: tuple[TradePoint, ...]


def minimum_endpoint_transfer(
    R: float, mu: float = MU_EARTH, r1: float = R1_REFERENCE
) -> EndpointInfimum:
    """Mathematical best-possible normalized delta-v at ``R``, and where it sits.

    ``saving_bar = dv_bar_H(R) - dv_best_bar(R)`` is zero throughout Region A and
    positive in Regions B and C. In Regions B and C it is a **lower bound
    requiring B -> infinity**, hence ``attainable = False`` there.
    """
    dv_hohmann = hohmann_total_normalized(R)
    dv_infinite = bielliptic_infinite_limit_normalized(R)
    dv_best, endpoint = bielliptic_infimum_normalized(R)
    saving_bar = dv_hohmann - dv_best
    v1 = math.sqrt(mu / r1)
    return EndpointInfimum(
        R=R,
        dv_hohmann_bar=dv_hohmann,
        dv_infinite_bar=dv_infinite,
        dv_best_bar=dv_best,
        best_endpoint=endpoint,
        saving_bar=saving_bar,
        saving_m_s=saving_bar * v1 * 1000.0,
        # The Hohmann endpoint IS a real transfer; the B->inf endpoint is not.
        attainable=(endpoint == "B->R+"),
    )


def practical_trade_metrics(
    R: float,
    B: float,
    mu: float = MU_EARTH,
    r1: float = R1_REFERENCE,
    body_radius: float = R_EARTH,
    label: str = "",
) -> TradePoint:
    """Delta-v, time and geometry metrics for one finite bi-elliptic design.

    Every physical quantity is delegated to the M2 production functions.
    ``dv_saving_m_s`` is positive when the bi-elliptic transfer is cheaper.
    """
    r2 = R * r1
    rb = B * r1
    v1 = math.sqrt(mu / r1)

    dv_hohmann_bar = hohmann_total_normalized(R)
    dv_bielliptic_bar = bielliptic_total_normalized(R, B)
    saving_bar = dv_hohmann_bar - dv_bielliptic_bar

    t_hohmann = hohmann_transfer_time(mu, r1, r2)
    t_bielliptic = bielliptic_transfer_time(mu, r1, r2, rb)

    flag, note = model_validity(rb)

    return TradePoint(
        R=R,
        B=B,
        B_over_R=B / R,
        label=label,
        dv_hohmann_bar=dv_hohmann_bar,
        dv_bielliptic_bar=dv_bielliptic_bar,
        saving_bar=saving_bar,
        dv_hohmann_m_s=dv_hohmann_bar * v1 * 1000.0,
        dv_bielliptic_m_s=dv_bielliptic_bar * v1 * 1000.0,
        dv_saving_m_s=saving_bar * v1 * 1000.0,
        t_hohmann_days=t_hohmann / _SECONDS_PER_DAY,
        t_bielliptic_days=t_bielliptic / _SECONDS_PER_DAY,
        time_ratio=t_bielliptic / t_hohmann,
        extra_time_days=(t_bielliptic - t_hohmann) / _SECONDS_PER_DAY,
        rb_km=rb,
        rb_altitude_km=rb - body_radius,
        rb_over_lunar_distance=rb / LUNAR_DISTANCE_KM,
        model_validity_flag=flag,
        model_validity_note=note,
    )


def break_even_trade(
    R: float, mu: float = MU_EARTH, r1: float = R1_REFERENCE
) -> TradePoint | None:
    """Metrics exactly at the Region-B break-even apoapsis, or ``None``.

    Returns ``None`` outside Region B: in Region A no break-even exists, and in
    Region C every ``B > R`` already wins so the boundary is open and there is no
    finite break-even point.

    At the returned point the delta-v saving is **zero by definition** -- the
    value of this point is the transfer time and apoapsis it already demands.
    """
    result = break_even_B(R)
    if result.B_crit is None:
        return None
    return practical_trade_metrics(R, result.B_crit, mu=mu, r1=r1, label="B_crit")


def evaluate_trade(
    R: float,
    B_values: Sequence[float] | None = None,
    mu: float = MU_EARTH,
    r1: float = R1_REFERENCE,
    body_radius: float = R_EARTH,
    include_break_even: bool = True,
) -> RadiusRatioTrade:
    """Evaluate the full finite-``B`` trade family at one radius ratio.

    ``B_values`` defaults to :data:`FINITE_B_OVER_R_FAMILY` scaled by ``R``.
    When ``include_break_even`` is set and ``R`` lies in Region B, the break-even
    apoapsis and two points beyond it (``1.1 B_crit`` and ``2 B_crit``) are
    appended -- these quantify how *slowly* savings actually appear past
    break-even, which matters most near ``R1*``.
    """
    R = float(R)
    if math.isnan(R) or math.isinf(R) or R <= 1.0:
        raise ValueError(f"evaluate_trade requires a finite R > 1, got {R!r}")

    result = break_even_B(R)
    points: list[TradePoint] = []

    if B_values is None:
        for ratio in FINITE_B_OVER_R_FAMILY:
            points.append(
                practical_trade_metrics(
                    R, R * ratio, mu=mu, r1=r1, body_radius=body_radius,
                    label=f"B/R={ratio:g}",
                )
            )
    else:
        for B in B_values:
            points.append(
                practical_trade_metrics(
                    R, float(B), mu=mu, r1=r1, body_radius=body_radius,
                    label=f"B={float(B):g}",
                )
            )

    if include_break_even and result.B_crit is not None:
        for multiple, name in ((1.0, "B_crit"), (1.1, "1.1*B_crit"), (2.0, "2*B_crit")):
            points.append(
                practical_trade_metrics(
                    R, multiple * result.B_crit, mu=mu, r1=r1,
                    body_radius=body_radius, label=name,
                )
            )

    points.sort(key=lambda point: point.B)

    return RadiusRatioTrade(
        R=R,
        regime=classify_radius_ratio(R),
        r2_km=R * r1,
        r2_altitude_km=R * r1 - body_radius,
        infimum=minimum_endpoint_transfer(R, mu=mu, r1=r1),
        B_crit=result.B_crit,
        points=tuple(points),
    )


def authoritative_R_grid(
    r_min: float = 1.01,
    r_max: float = 100.0,
    n_global: int = 600,
    n_window: int = 241,
    window_half_width: float = 0.5,
) -> np.ndarray:
    """The authoritative radius-ratio grid used by every M3 sweep and figure.

    Construction (documented in DESIGN.md M3.2):

    - ``n_global`` logarithmically spaced points over ``[r_min, r_max]``;
    - ``n_window`` linearly spaced points in a window of half-width
      ``window_half_width`` about each threshold, to resolve them locally;
    - the two thresholds themselves, and a few points immediately either side;
    - deduplicated and sorted.

    The grid is for *plotting and tabulation only*. The thresholds are always
    obtained from the M2 root solvers, never read off this grid.
    """
    if not 1.0 < r_min < r_max:
        raise ValueError(f"require 1 < r_min < r_max, got {r_min!r}, {r_max!r}")

    r1_star = threshold_R1()
    r2_star = threshold_R2()

    pieces = [np.logspace(math.log10(r_min), math.log10(r_max), n_global)]
    for centre in (r1_star, r2_star):
        pieces.append(
            np.linspace(
                max(r_min, centre - window_half_width),
                min(r_max, centre + window_half_width),
                n_window,
            )
        )
        pieces.append(
            np.array(
                [
                    centre * (1.0 - 1e-9),
                    centre,
                    centre * (1.0 + 1e-9),
                    centre * (1.0 + 1e-6),
                    centre * (1.0 + 1e-3),
                ]
            )
        )
    grid = np.unique(np.concatenate(pieces))
    return grid[(grid >= r_min) & (grid <= r_max)]


def break_even_locus_grid(n_points: int = 220) -> np.ndarray:
    """Radius-ratio grid spanning Region B, dense at both ends.

    Parameterised as ``R = R1* + (R2* - R1*) * s`` with ``s`` logarithmically
    spaced. That resolves the divergence of ``B_crit`` as ``R -> R1*+`` while
    still reaching the far end where ``B_crit -> R``.
    """
    r1_star = threshold_R1()
    r2_star = threshold_R2()
    span = r2_star - r1_star
    s_low = np.logspace(-5.0, math.log10(0.5), n_points // 2)
    s_high = 1.0 - np.logspace(-6.0, math.log10(0.5), n_points // 2)
    s = np.unique(np.concatenate([s_low, s_high]))
    return r1_star + span * s


class BreakEvenLocusPoint(NamedTuple):
    """One point on the Region-B break-even locus."""

    R: float
    B_crit: float
    B_crit_over_R: float
    rb_km: float
    rb_altitude_km: float
    rb_over_lunar_distance: float
    time_ratio_at_break_even: float
    model_validity_flag: str


def break_even_locus(
    R_values: Sequence[float] | np.ndarray | None = None,
    mu: float = MU_EARTH,
    r1: float = R1_REFERENCE,
    body_radius: float = R_EARTH,
) -> tuple[BreakEvenLocusPoint, ...]:
    """Compute ``B_crit(R)`` across Region B by root solving at every ``R``.

    Each point calls the M2 break-even solver directly. Nothing is interpolated:
    the locus is a set of solved roots, not a fitted curve.
    """
    if R_values is None:
        R_values = break_even_locus_grid()

    locus: list[BreakEvenLocusPoint] = []
    for R in R_values:
        R = float(R)
        result = break_even_B(R)
        if result.B_crit is None:
            continue
        point = practical_trade_metrics(
            R, result.B_crit, mu=mu, r1=r1, body_radius=body_radius, label="B_crit"
        )
        locus.append(
            BreakEvenLocusPoint(
                R=R,
                B_crit=result.B_crit,
                B_crit_over_R=result.B_crit / R,
                rb_km=point.rb_km,
                rb_altitude_km=point.rb_altitude_km,
                rb_over_lunar_distance=point.rb_over_lunar_distance,
                time_ratio_at_break_even=point.time_ratio,
                model_validity_flag=point.model_validity_flag,
            )
        )
    return tuple(locus)


def sweep_radius_ratios(
    R_values: Sequence[float] | np.ndarray | None = None,
    mu: float = MU_EARTH,
    r1: float = R1_REFERENCE,
) -> tuple[EndpointInfimum, ...]:
    """Endpoint-infimum envelope over a radius-ratio grid.

    Returns one :class:`EndpointInfimum` per ``R``. This is the data behind the
    top panel of the main crossover figure.
    """
    if R_values is None:
        R_values = authoritative_R_grid()
    return tuple(minimum_endpoint_transfer(float(R), mu=mu, r1=r1) for R in R_values)


class Recommendation(NamedTuple):
    """Scoped engineering guidance at one radius ratio."""

    R: float
    regime: str
    headline: str
    delta_v_verdict: str
    practical_caveat: str
    model_validity_flag: str


def recommendation_for(
    R: float, mu: float = MU_EARTH, r1: float = R1_REFERENCE
) -> Recommendation:
    """Conditional recommendation at ``R``, derived entirely from computed values.

    Deliberately separates *mathematical* delta-v dominance from *practical*
    benefit. It never asserts that a bi-elliptic transfer is simply "better"
    above a threshold, and it never presents ``B -> infinity`` as a design.
    """
    regime = classify_radius_ratio(R)
    infimum = minimum_endpoint_transfer(R, mu=mu, r1=r1)
    trade = evaluate_trade(R, mu=mu, r1=r1)

    # Best finite strategy actually evaluated, and whether it saves anything.
    best_point = max(trade.points, key=lambda point: point.dv_saving_m_s)

    if regime == REGIME_HOHMANN_ONLY:
        headline = "Use the Hohmann transfer."
        verdict = (
            f"Hohmann is delta-v superior to EVERY admissible bi-elliptic transfer "
            f"at R = {R:.6g}. No choice of intermediate apoapsis beats it."
        )
        caveat = (
            "The bi-elliptic family offers no delta-v benefit here at any B, so "
            "its longer transfer time buys nothing."
        )
    elif regime == REGIME_LARGE_B_BIELLIPTIC:
        assert trade.B_crit is not None
        break_even = break_even_trade(R, mu=mu, r1=r1)
        assert break_even is not None
        headline = "Hohmann unless a very distant apoapsis is acceptable."
        verdict = (
            f"A bi-elliptic transfer beats Hohmann only if B exceeds "
            f"B_crit = {trade.B_crit:.6g} (r_b = {break_even.rb_km:.4g} km, "
            f"{break_even.rb_over_lunar_distance:.3g} lunar distances). The entire "
            f"mathematical prize, requiring B -> infinity, is "
            f"{infimum.saving_m_s:.3g} m/s."
        )
        caveat = (
            f"Merely reaching break-even already costs "
            f"{break_even.time_ratio:.4g}x the Hohmann transfer time "
            f"({break_even.t_bielliptic_days:.4g} d vs "
            f"{break_even.t_hohmann_days:.4g} d) for zero saving. The best finite "
            f"strategy evaluated saves {best_point.dv_saving_m_s:.3g} m/s at "
            f"{best_point.time_ratio:.4g}x the Hohmann time."
        )
    else:
        headline = "Bi-elliptic is delta-v favourable; weigh the saving against time."
        verdict = (
            f"Every admissible B > R beats Hohmann at R = {R:.6g}. The mathematical "
            f"infimum, requiring B -> infinity, saves {infimum.saving_m_s:.3g} m/s."
        )
        caveat = (
            f"Delta-v dominance does not imply a large practical benefit: the best "
            f"finite strategy evaluated saves {best_point.dv_saving_m_s:.3g} m/s at "
            f"{best_point.time_ratio:.4g}x the Hohmann transfer time "
            f"({best_point.extra_time_days:.4g} extra days)."
        )

    flag, _ = model_validity(best_point.rb_km)
    return Recommendation(
        R=R,
        regime=regime,
        headline=headline,
        delta_v_verdict=verdict,
        practical_caveat=caveat,
        model_validity_flag=flag,
    )


def recommendation_summary(
    mu: float = MU_EARTH, r1: float = R1_REFERENCE
) -> dict[str, object]:
    """Banded recommendation for the whole radius-ratio domain.

    Built from the recomputed thresholds, not from literals.
    """
    r1_star = threshold_R1()
    r2_star = threshold_R2()
    geo_R = GEO_RADIUS_KM / r1

    return {
        "band_A": {
            "condition": f"1 < R < {r1_star:.14f}",
            "regime": REGIME_HOHMANN_ONLY,
            "guidance": (
                "Hohmann is delta-v superior to every admissible bi-elliptic "
                "transfer. No intermediate apoapsis helps."
            ),
        },
        "band_B": {
            "condition": f"{r1_star:.14f} < R < {r2_star:.14f}",
            "regime": REGIME_LARGE_B_BIELLIPTIC,
            "guidance": (
                "A bi-elliptic transfer can beat Hohmann ONLY if the intermediate "
                "apoapsis exceeds B_crit(R). Near the lower threshold the required "
                "apoapsis and transfer time are absurdly large for a vanishing "
                "delta-v saving."
            ),
        },
        "band_C": {
            "condition": f"R > {r2_star:.14f}",
            "regime": REGIME_ALL_BIELLIPTIC,
            "guidance": (
                "Every admissible B > R is delta-v better than Hohmann, but the "
                "practical size of the saving must still be weighed against the "
                "transfer-time penalty. Delta-v dominance is not the same as a "
                "material benefit."
            ),
        },
        "infinite_B_status": (
            "B -> infinity is a MATHEMATICAL INFIMUM with unbounded transfer time. "
            "It is a lower bound on achievable delta-v, never an engineering "
            "recommendation and never a realizable transfer."
        ),
        "earth_leo_to_geo": {
            "R": geo_R,
            "regime": classify_radius_ratio(geo_R),
            "guidance": (
                "A 300 km LEO to GEO transfer has R = "
                f"{geo_R:.9f}, firmly inside Region A. Hohmann remains preferred "
                "within this idealized coplanar impulsive two-body model."
            ),
        },
        "model_validity_warning": (
            f"Apoapsis radii approaching or exceeding the lunar distance "
            f"({LUNAR_DISTANCE_KM:.0f} km) invalidate the Earth-only two-body "
            "model used throughout. Results quoted beyond that radius are "
            "mathematical extrapolations of the idealized model, not trajectory "
            "designs."
        ),
    }
