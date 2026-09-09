"""Milestone 3 verification: trade sweep, break-even locus, envelope, recommendation.

Every expectation is derived from production code rather than transcribed, so
these tests fail if the M3 layer ever stops agreeing with the verified M2
primitives it is supposed to be composing.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from bielliptic_crossover import (
    MU_EARTH,
    R1_REFERENCE,
    R_EARTH,
    V1_REFERENCE,
    bielliptic_infinite_limit_normalized,
    bielliptic_total_normalized,
    bielliptic_transfer_time,
    break_even_B,
    classify_radius_ratio,
    hohmann_total_normalized,
    hohmann_transfer_time,
    threshold_R1,
    threshold_R2,
)
from bielliptic_crossover.crossover import (
    REGIME_ALL_BIELLIPTIC,
    REGIME_HOHMANN_ONLY,
    REGIME_LARGE_B_BIELLIPTIC,
    bielliptic_infimum_normalized,
)
from bielliptic_crossover.trade import (
    EARTH_HILL_RADIUS_KM,
    FINITE_B_OVER_R_FAMILY,
    GEO_RADIUS_KM,
    LUNAR_DISTANCE_KM,
    MODEL_CAUTION,
    MODEL_INVALID,
    MODEL_OK,
    authoritative_R_grid,
    break_even_locus,
    break_even_locus_grid,
    break_even_trade,
    evaluate_trade,
    minimum_endpoint_transfer,
    model_validity,
    practical_trade_metrics,
    recommendation_for,
    recommendation_summary,
    sweep_radius_ratios,
)

DAY = 86400.0


# ---------------------------------------------------------------- A. resolution


def test_thresholds_are_independent_of_sweep_resolution() -> None:
    """The grid is for plotting; thresholds come from the root solvers."""
    reference_R1, reference_R2 = threshold_R1(), threshold_R2()
    for n_global, n_window in ((80, 21), (600, 241), (1500, 601)):
        grid = authoritative_R_grid(n_global=n_global, n_window=n_window)
        assert grid.size > 0
        # The solvers are untouched by the grid used around them.
        assert threshold_R1() == reference_R1
        assert threshold_R2() == reference_R2


def test_envelope_is_independent_of_grid_density() -> None:
    """Sampling more finely changes only where we sample, never the values."""
    coarse = authoritative_R_grid(n_global=120, n_window=31)
    for R in coarse[::7]:
        item = minimum_endpoint_transfer(float(R))
        assert item.dv_best_bar == pytest.approx(
            min(hohmann_total_normalized(float(R)),
                bielliptic_infinite_limit_normalized(float(R))),
            rel=1e-15,
        )


def test_authoritative_grid_is_sorted_unique_and_bounded() -> None:
    grid = authoritative_R_grid()
    assert np.all(np.diff(grid) > 0)
    assert grid.min() >= 1.01
    assert grid.max() <= 100.0
    # Both thresholds are resolved by a dense local window.
    for threshold in (threshold_R1(), threshold_R2()):
        nearby = grid[np.abs(grid - threshold) < 0.5]
        assert nearby.size >= 200


def test_grid_rejects_bad_bounds() -> None:
    with pytest.raises(ValueError, match="r_min"):
        authoritative_R_grid(r_min=0.5)


# ------------------------------------------------- B. locus is solved, not fitted


def test_break_even_locus_points_are_solved_roots_not_interpolants() -> None:
    """Each locus point must satisfy the defining equality to machine precision.

    An interpolated or fitted curve would not.
    """
    for point in break_even_locus():
        residual = bielliptic_total_normalized(point.R, point.B_crit) - (
            hohmann_total_normalized(point.R)
        )
        assert abs(residual) < 1e-15
        assert point.B_crit == break_even_B(point.R).B_crit


def test_locus_matches_independent_bisection() -> None:
    """Cross-check a few locus points against a bisection written here."""
    for point in break_even_locus()[::40]:
        R = point.R
        dv_hohmann = hohmann_total_normalized(R)
        low, high = R * (1.0 + 1e-9), 1.0e20
        for _ in range(300):
            middle = math.sqrt(low * high)
            if bielliptic_total_normalized(R, middle) - dv_hohmann <= 0.0:
                high = middle
            else:
                low = middle
        assert math.sqrt(low * high) == pytest.approx(point.B_crit, rel=1e-6)


# ---------------------------------------------------------------- C, D, E


def test_break_even_B_decreases_strictly_across_region_B() -> None:
    locus = break_even_locus()
    ratios = np.array([point.B_crit_over_R for point in locus])
    absolute = np.array([point.B_crit for point in locus])
    assert np.all(np.diff(ratios) < 0.0)
    assert np.all(np.diff(absolute) < 0.0)


def test_break_even_B_diverges_toward_R1_star() -> None:
    r1_star = threshold_R1()
    previous = None
    for offset in (1e-2, 1e-3, 1e-4, 1e-5):
        result = break_even_B(r1_star + offset)
        assert result.B_crit is not None
        if previous is not None:
            # Each decade closer to R1* multiplies B_crit by roughly ten.
            assert result.B_crit > 5.0 * previous
        previous = result.B_crit
    assert previous > 1e5


def test_break_even_B_diverges_as_inverse_first_power() -> None:
    """Measured asymptotic law: B_crit/R ~ (R - R1*)^-1."""
    r1_star = threshold_R1()
    offsets = np.array([1e-5, 1e-4, 1e-3, 1e-2])
    ratios = np.array(
        [break_even_B(r1_star + float(d)).B_crit / (r1_star + float(d)) for d in offsets]
    )
    slope, _ = np.polyfit(np.log10(offsets), np.log10(ratios), 1)
    assert slope == pytest.approx(-1.0, abs=5e-3)


def test_break_even_ratio_tends_to_one_toward_R2_star() -> None:
    """B_crit/R -> 1 as R -> R2*-, measured within the resolvable band."""
    r2_star = threshold_R2()
    previous = None
    for offset in (1e-2, 1e-3, 1e-4, 1e-5):
        result = break_even_B(r2_star - offset)
        assert result.B_crit is not None
        ratio = result.B_crit / (r2_star - offset)
        assert ratio > 1.0
        if previous is not None:
            assert ratio < previous
        previous = ratio
    assert previous == pytest.approx(1.0, abs=1e-5)


def test_break_even_excess_vanishes_linearly_at_R2_star() -> None:
    """B_crit/R - 1 ~ 0.306 * (R2* - R): linear, so it passes below any fixed
    probe offset at some finite distance from the threshold."""
    r2_star = threshold_R2()
    offsets = np.array([1e-5, 1e-4, 1e-3, 1e-2])
    excess = np.array(
        [break_even_B(r2_star - float(d)).B_crit / (r2_star - float(d)) - 1.0
         for d in offsets]
    )
    slope, intercept = np.polyfit(np.log10(offsets), np.log10(excess), 1)
    assert slope == pytest.approx(1.0, abs=5e-3)
    assert 10.0**intercept == pytest.approx(0.306, rel=0.05)


def test_threshold_neighbourhoods_are_resolution_limited_not_wrong() -> None:
    """Documented precision boundaries at BOTH thresholds.

    These are resolution limits of double-precision arithmetic, not errors, and
    they are recorded so nobody later "fixes" them into a false certainty.

    Near R2*: ``break_even_B`` steps off the trivial B = R root by a fixed
    relative probe of 1e-6. Because ``B_crit/R - 1`` shrinks linearly to zero at
    R2*, it drops below that probe once ``R2* - R`` is under about 3.3e-6, and
    the solver then reports Region C while the threshold-based classifier still
    reports Region B. Inside that band the break-even apoapsis is less than one
    part in a million above the target radius, so no physically meaningful
    distinction is being lost.

    Near R1*: the double-rounded ``threshold_R1()`` sits about 2.2e-14 above the
    true root, so at exactly that R a finite break-even does exist mathematically
    (about 4e15) even though the classifier assigns Region A. Both readings are
    defensible; the value is a mathematical extrapolation with no physical
    meaning.
    """
    r1_star, r2_star = threshold_R1(), threshold_R2()

    # Outside the band the two agree.
    for offset in (1e-5, 1e-4, 1e-3):
        R = r2_star - offset
        assert break_even_B(R).regime == classify_radius_ratio(R)

    # Inside it they diverge, in the documented direction.
    for offset in (1e-6, 1e-7):
        R = r2_star - offset
        assert classify_radius_ratio(R) == REGIME_LARGE_B_BIELLIPTIC
        assert break_even_B(R).regime == REGIME_ALL_BIELLIPTIC
        assert break_even_B(R).B_crit is None

    # At R1* exactly: classifier says Region A, solver finds a huge finite root.
    assert classify_radius_ratio(r1_star) == REGIME_HOHMANN_ONLY
    at_threshold = break_even_B(r1_star)
    assert at_threshold.B_crit is not None
    assert at_threshold.B_crit > 1e12
    assert break_even_trade(r1_star).model_validity_flag == MODEL_INVALID


def test_break_even_time_penalty_diverges_as_three_halves_power() -> None:
    """t_B/t_H at break-even ~ (R - R1*)^(-3/2), since t ~ B^(3/2)."""
    r1_star = threshold_R1()
    offsets = np.array([1e-5, 1e-4, 1e-3, 1e-2])
    ratios = np.array(
        [break_even_trade(r1_star + float(d)).time_ratio for d in offsets]
    )
    slope, _ = np.polyfit(np.log10(offsets), np.log10(ratios), 1)
    assert slope == pytest.approx(-1.5, abs=1e-2)


def test_locus_grid_stays_inside_region_B() -> None:
    r1_star, r2_star = threshold_R1(), threshold_R2()
    grid = break_even_locus_grid()
    assert np.all(grid > r1_star)
    assert np.all(grid < r2_star)


# ---------------------------------------------------------------- F, G


@pytest.mark.parametrize("R", [1.5, 5.0, 10.0, 11.9, 12.0, 14.0, 16.0, 20.0, 50.0, 99.0])
def test_envelope_equals_endpoint_minimum(R: float) -> None:
    item = minimum_endpoint_transfer(R)
    value, endpoint = bielliptic_infimum_normalized(R)
    assert item.dv_best_bar == pytest.approx(value, rel=1e-15)
    assert item.best_endpoint == endpoint
    assert item.dv_best_bar == pytest.approx(
        min(item.dv_hohmann_bar, item.dv_infinite_bar), rel=1e-15
    )


@pytest.mark.parametrize("R", [1.5, 5.0, 11.0])
def test_saving_is_exactly_zero_in_region_A(R: float) -> None:
    item = minimum_endpoint_transfer(R)
    assert item.saving_bar == 0.0
    assert item.saving_m_s == 0.0
    assert item.best_endpoint == "B->R+"
    assert item.attainable is True


@pytest.mark.parametrize("R", [12.0, 14.0, 16.0, 30.0, 99.0])
def test_saving_is_positive_and_unattainable_in_regions_B_and_C(R: float) -> None:
    item = minimum_endpoint_transfer(R)
    assert item.saving_bar > 0.0
    assert item.best_endpoint == "B->inf"
    # The B->infinity endpoint is NOT a realizable transfer.
    assert item.attainable is False


@pytest.mark.parametrize("R", [2.0, 10.0, 12.0, 15.0, 16.0, 20.0, 50.0])
def test_M3_introduces_no_finite_interior_optimum(R: float) -> None:
    """No finite B evaluated anywhere in the M3 layer may beat the endpoint infimum."""
    infimum, _ = bielliptic_infimum_normalized(R)
    trade = evaluate_trade(R)
    for point in trade.points:
        assert point.dv_bielliptic_bar >= infimum - 1e-15
    dense = [
        bielliptic_total_normalized(R, R * (1.0 + 1e-9) * 10.0 ** (10.0 * i / 500))
        for i in range(501)
    ]
    assert min(dense) >= infimum - 1e-14


# ---------------------------------------------------------------- H. delegation


@pytest.mark.parametrize("R", [2.0, 12.0, 16.0, 50.0])
@pytest.mark.parametrize("ratio", FINITE_B_OVER_R_FAMILY)
def test_finite_strategies_reproduce_direct_M2_calls(R: float, ratio: float) -> None:
    """The trade layer must delegate, never re-derive."""
    B = R * ratio
    point = practical_trade_metrics(R, B)
    assert point.dv_bielliptic_bar == bielliptic_total_normalized(R, B)
    assert point.dv_hohmann_bar == hohmann_total_normalized(R)
    assert point.saving_bar == pytest.approx(
        hohmann_total_normalized(R) - bielliptic_total_normalized(R, B), abs=1e-18
    )
    assert point.t_bielliptic_days == pytest.approx(
        bielliptic_transfer_time(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, B * R1_REFERENCE)
        / DAY,
        rel=1e-15,
    )
    assert point.t_hohmann_days == pytest.approx(
        hohmann_transfer_time(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE) / DAY, rel=1e-15
    )
    assert point.rb_km == pytest.approx(B * R1_REFERENCE, rel=1e-15)
    assert point.rb_altitude_km == pytest.approx(B * R1_REFERENCE - R_EARTH, rel=1e-12)


def test_evaluate_trade_includes_break_even_points_only_in_region_B() -> None:
    labels_B = {point.label for point in evaluate_trade(13.0).points}
    assert {"B_crit", "1.1*B_crit", "2*B_crit"} <= labels_B
    for R in (5.0, 20.0):
        labels = {point.label for point in evaluate_trade(R).points}
        assert "B_crit" not in labels


def test_evaluate_trade_points_are_sorted_by_B() -> None:
    for R in (12.0, 16.0):
        values = [point.B for point in evaluate_trade(R).points]
        assert values == sorted(values)


def test_evaluate_trade_rejects_R_at_or_below_one() -> None:
    with pytest.raises(ValueError, match="R > 1"):
        evaluate_trade(1.0)


# ---------------------------------------------------------------- I. scaling


@pytest.mark.parametrize("mu, r1", [(MU_EARTH, R1_REFERENCE), (42828.375214, 3796.19),
                                    (1.32712440018e11, 1.496e8), (1.0e3, 1.0e2)])
def test_dimensional_scaling_is_exact(mu: float, r1: float) -> None:
    """Normalized quantities are invariant; dimensional ones scale as sqrt(mu/r1)."""
    v1 = math.sqrt(mu / r1)
    for R in (2.0, 12.0, 50.0):
        for ratio in (1.5, 5.0):
            point = practical_trade_metrics(R, R * ratio, mu=mu, r1=r1)
            assert point.dv_bielliptic_bar == pytest.approx(
                bielliptic_total_normalized(R, R * ratio), rel=1e-15
            )
            assert point.dv_bielliptic_m_s == pytest.approx(
                point.dv_bielliptic_bar * v1 * 1000.0, rel=1e-15
            )
            # Time ratio is dimensionless and therefore body-independent.
            reference = practical_trade_metrics(R, R * ratio)
            assert point.time_ratio == pytest.approx(reference.time_ratio, rel=1e-13)


def test_saving_scales_linearly_with_v1() -> None:
    base = minimum_endpoint_transfer(20.0)
    scaled = minimum_endpoint_transfer(20.0, mu=4.0 * MU_EARTH, r1=R1_REFERENCE)
    assert scaled.saving_m_s == pytest.approx(2.0 * base.saving_m_s, rel=1e-13)
    assert scaled.saving_bar == pytest.approx(base.saving_bar, rel=1e-15)


# ---------------------------------------------------------------- J. GEO


def test_geo_is_firmly_region_A() -> None:
    geo_R = GEO_RADIUS_KM / R1_REFERENCE
    assert geo_R == pytest.approx(6.313737562, abs=1e-8)
    assert classify_radius_ratio(geo_R) == REGIME_HOHMANN_ONLY
    assert geo_R < threshold_R1()
    item = minimum_endpoint_transfer(geo_R)
    assert item.saving_bar == 0.0
    assert break_even_B(geo_R).B_crit is None
    # Every finite strategy is strictly worse than Hohmann.
    for point in evaluate_trade(geo_R).points:
        assert point.dv_saving_m_s < 0.0


def test_geo_recommendation_says_use_hohmann() -> None:
    recommendation = recommendation_for(GEO_RADIUS_KM / R1_REFERENCE)
    assert recommendation.regime == REGIME_HOHMANN_ONLY
    assert "Hohmann" in recommendation.headline


# ---------------------------------------------------------------- K, L. time


@pytest.mark.parametrize("R", [2.0, 12.0, 16.0, 50.0])
@pytest.mark.parametrize("ratio", [1.0, 1.25, 2.0, 10.0])
def test_time_ratio_always_exceeds_one(R: float, ratio: float) -> None:
    """Even at B = R the bi-elliptic path is slower (the M2 half-lap degeneracy)."""
    point = practical_trade_metrics(R, R * ratio)
    assert point.time_ratio > 1.0
    assert point.extra_time_days > 0.0


def test_time_ratio_at_B_equals_R_reproduces_the_M2_degeneracy() -> None:
    """t_B(R,R) - t_H(R) = pi * R^(3/2) * t_star, not zero."""
    t_star = math.sqrt(R1_REFERENCE**3 / MU_EARTH)
    for R in (2.0, 12.0, 16.0):
        point = practical_trade_metrics(R, R)
        excess_days = point.extra_time_days
        assert excess_days == pytest.approx(math.pi * R**1.5 * t_star / DAY, rel=1e-13)


@pytest.mark.parametrize("R", [12.0, 16.0])
def test_large_B_time_growth_is_three_halves_power(R: float) -> None:
    big = 1e6 * R
    first = practical_trade_metrics(R, big).t_bielliptic_days
    second = practical_trade_metrics(R, 2.0 * big).t_bielliptic_days
    assert second / first == pytest.approx(2.0**1.5, rel=1e-6)


# ---------------------------------------------------------------- M. R=12


def test_R12_practical_penalty_reproduces_M2() -> None:
    """The M2 headline trade must survive unchanged through the M3 layer."""
    trade = evaluate_trade(12.0)
    assert trade.regime == REGIME_LARGE_B_BIELLIPTIC
    assert trade.B_crit == 815.8202504753092
    assert trade.infimum.saving_m_s == pytest.approx(3.0374130154, abs=1e-9)

    break_even = break_even_trade(12.0)
    assert break_even is not None
    assert break_even.dv_saving_m_s == pytest.approx(0.0, abs=1e-12)
    assert break_even.time_ratio == pytest.approx(1006.1987, rel=1e-6)
    assert break_even.t_hohmann_days == pytest.approx(0.520859, rel=1e-5)
    assert break_even.t_bielliptic_days == pytest.approx(524.0877, rel=1e-6)
    assert break_even.rb_over_lunar_distance == pytest.approx(14.17, rel=1e-3)
    assert break_even.model_validity_flag == MODEL_INVALID


def test_R12_moderate_apoapsis_is_worse_than_hohmann() -> None:
    """In Region B a modest B is strictly WORSE, which the trade must expose."""
    for ratio in (1.25, 1.5, 2.0, 5.0, 10.0):
        point = practical_trade_metrics(12.0, 12.0 * ratio)
        assert point.dv_saving_m_s < 0.0


def test_R16_trade_shows_immediate_but_small_benefit() -> None:
    trade = evaluate_trade(16.0)
    assert trade.regime == REGIME_ALL_BIELLIPTIC
    assert trade.B_crit is None
    # Every B > R already wins, but the smallest strategy wins by very little.
    smallest = min(trade.points, key=lambda point: point.B)
    assert 0.0 < smallest.dv_saving_m_s < 20.0
    assert smallest.time_ratio > 1.0
    assert trade.infimum.saving_m_s == pytest.approx(142.7136, rel=1e-5)


def test_just_above_R2_star_delta_v_dominance_is_nearly_worthless() -> None:
    """'Every B wins' does not mean the win is worth anything."""
    R = threshold_R2() * (1.0 + 1e-3)
    point = practical_trade_metrics(R, R * 1.001)
    assert point.dv_saving_m_s > 0.0
    assert point.dv_saving_m_s < 0.01  # sub-centimetre-per-second
    assert point.time_ratio > 3.5


# ---------------------------------------------------------------- model validity


def test_model_validity_thresholds() -> None:
    assert model_validity(0.01 * LUNAR_DISTANCE_KM)[0] == MODEL_OK
    assert model_validity(0.5 * LUNAR_DISTANCE_KM)[0] == MODEL_CAUTION
    assert model_validity(2.0 * LUNAR_DISTANCE_KM)[0] == MODEL_INVALID
    assert model_validity(2.0 * EARTH_HILL_RADIUS_KM)[0] == MODEL_INVALID
    assert "lunar" in model_validity(2.0 * LUNAR_DISTANCE_KM)[1]
    assert "Hill" in model_validity(2.0 * EARTH_HILL_RADIUS_KM)[1]
    assert model_validity(0.01 * LUNAR_DISTANCE_KM)[1] == ""


def test_region_B_break_even_is_beyond_the_moon_at_low_R() -> None:
    """Across the lower part of Region B, break-even apoapsis exceeds lunar distance."""
    for R in (12.0, 12.5, 12.8):
        break_even = break_even_trade(R)
        assert break_even is not None
        assert break_even.rb_over_lunar_distance > 1.0
        assert break_even.model_validity_flag == MODEL_INVALID


# ---------------------------------------------------------------- N. recommendation


def test_recommendation_summary_uses_recomputed_thresholds() -> None:
    summary = recommendation_summary()
    assert f"{threshold_R1():.14f}" in summary["band_A"]["condition"]
    assert f"{threshold_R2():.14f}" in summary["band_C"]["condition"]
    assert summary["band_A"]["regime"] == REGIME_HOHMANN_ONLY
    assert summary["band_B"]["regime"] == REGIME_LARGE_B_BIELLIPTIC
    assert summary["band_C"]["regime"] == REGIME_ALL_BIELLIPTIC


def test_recommendation_never_calls_infinite_B_optimal() -> None:
    summary = recommendation_summary()
    text = " ".join(str(value) for value in summary.values()).lower()
    assert "infimum" in text
    assert "never an engineering recommendation" in text
    for R in (5.0, 12.0, 16.0, 50.0):
        recommendation = recommendation_for(R)
        combined = (
            recommendation.headline
            + recommendation.delta_v_verdict
            + recommendation.practical_caveat
        ).lower()
        assert "b -> infinity is optimal" not in combined
        assert "optimal transfer" not in combined


def test_recommendation_regimes_match_the_classifier() -> None:
    for R in (2.0, 11.0, 12.0, 15.0, 16.0, 40.0):
        assert recommendation_for(R).regime == classify_radius_ratio(R)


def test_region_B_recommendation_quantifies_the_break_even_cost() -> None:
    recommendation = recommendation_for(12.0)
    assert "B_crit" in recommendation.delta_v_verdict
    assert "815" in recommendation.delta_v_verdict
    assert "Hohmann transfer time" in recommendation.practical_caveat


def test_earth_leo_to_geo_band_is_region_A() -> None:
    summary = recommendation_summary()
    assert summary["earth_leo_to_geo"]["regime"] == REGIME_HOHMANN_ONLY


# ---------------------------------------------------------------- sweep helpers


def test_sweep_radius_ratios_matches_pointwise_evaluation() -> None:
    grid = authoritative_R_grid(n_global=60, n_window=11)
    swept = sweep_radius_ratios(grid)
    assert len(swept) == grid.size
    for item, R in zip(swept, grid):
        assert item.R == float(R)
        assert item.dv_best_bar == minimum_endpoint_transfer(float(R)).dv_best_bar


def test_evaluate_trade_reports_target_geometry() -> None:
    trade = evaluate_trade(20.0)
    assert trade.r2_km == pytest.approx(20.0 * R1_REFERENCE, rel=1e-15)
    assert trade.r2_altitude_km == pytest.approx(20.0 * R1_REFERENCE - R_EARTH, rel=1e-12)


def test_v1_reference_is_unchanged_from_M2() -> None:
    assert V1_REFERENCE == pytest.approx(7.72576063698292, abs=1e-11)
