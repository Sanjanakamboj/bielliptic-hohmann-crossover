"""Verification D, E and F: thresholds, region classifier, break-even B."""

from __future__ import annotations

import itertools
import math

import pytest

from bielliptic_crossover.bielliptic import (
    bielliptic_infinite_limit_normalized,
    bielliptic_total_normalized,
)
from bielliptic_crossover.crossover import (
    REGIME_ALL_BIELLIPTIC,
    REGIME_HOHMANN_ONLY,
    REGIME_LARGE_B_BIELLIPTIC,
    bielliptic_infimum_normalized,
    break_even_B,
    classify_radius_ratio,
    threshold_R1,
    threshold_R1_from_polynomial,
    threshold_R2,
    threshold_R2_from_cubic,
)
from bielliptic_crossover.hohmann import hohmann_total_normalized

# M1 accepted values, used as regression targets only -- never as inputs.
M1_R1 = 11.938765472645870716
M1_R2 = 15.581718738763179213


# --------------------------------------------------------------- threshold R1


def test_R1_from_transfer_equations() -> None:
    assert threshold_R1() == pytest.approx(M1_R1, rel=1e-12)


def test_R1_from_polynomial() -> None:
    assert threshold_R1_from_polynomial() == pytest.approx(M1_R1, rel=1e-12)


def test_R1_two_independent_paths_agree() -> None:
    assert threshold_R1() == pytest.approx(threshold_R1_from_polynomial(), rel=1e-12)


def test_R1_actually_satisfies_its_defining_condition() -> None:
    r1_star = threshold_R1()
    assert hohmann_total_normalized(r1_star) == pytest.approx(
        bielliptic_infinite_limit_normalized(r1_star), abs=1e-15
    )


def test_R1_satisfies_the_reduced_algebraic_form() -> None:
    """(R-1)/sqrt(R+1) = sqrt(R) + 1 - sqrt(2)."""
    R = threshold_R1()
    assert (R - 1.0) / math.sqrt(R + 1.0) == pytest.approx(
        math.sqrt(R) + 1.0 - math.sqrt(2.0), abs=1e-14
    )


def test_R1_polynomial_root_satisfies_its_cubic() -> None:
    u = math.sqrt(threshold_R1_from_polynomial())
    assert u**3 - (1.0 + 2.0 * math.sqrt(2.0)) * u**2 + u + 1.0 == pytest.approx(
        0.0, abs=1e-12
    )


def test_hohmann_beats_infinite_limit_below_R1_and_loses_above() -> None:
    r1_star = threshold_R1()
    below, above = r1_star * 0.99, r1_star * 1.01
    assert hohmann_total_normalized(below) < bielliptic_infinite_limit_normalized(below)
    assert hohmann_total_normalized(above) > bielliptic_infinite_limit_normalized(above)


# --------------------------------------------------------------- threshold R2


def test_R2_from_slope_condition() -> None:
    assert threshold_R2() == pytest.approx(M1_R2, rel=1e-12)


def test_R2_from_cubic() -> None:
    assert threshold_R2_from_cubic() == pytest.approx(M1_R2, rel=1e-12)


def test_R2_two_independent_paths_agree() -> None:
    assert threshold_R2() == pytest.approx(threshold_R2_from_cubic(), rel=1e-12)


def test_R2_satisfies_its_cubic() -> None:
    R = threshold_R2()
    assert R**3 - 15.0 * R**2 - 9.0 * R - 1.0 == pytest.approx(0.0, abs=1e-11)


def test_R2_satisfies_the_pre_squaring_form() -> None:
    """(1+R)^(3/2) = sqrt(2) * (1+3R), the form before radicals were cleared."""
    R = threshold_R2()
    assert (1.0 + R) ** 1.5 == pytest.approx(math.sqrt(2.0) * (1.0 + 3.0 * R), rel=1e-13)


# --------------------------------------------------- the two are distinct


def test_the_two_thresholds_are_distinct_and_ordered() -> None:
    assert threshold_R1() < threshold_R2()
    assert threshold_R2() - threshold_R1() == pytest.approx(3.6429532661, abs=1e-8)


def test_between_thresholds_infinite_limit_wins_but_small_B_does_not() -> None:
    """The defining difference between the two thresholds, stated as a test."""
    R = 0.5 * (threshold_R1() + threshold_R2())  # strictly inside Region B
    assert bielliptic_infinite_limit_normalized(R) < hohmann_total_normalized(R)
    assert bielliptic_total_normalized(R, R * 1.001) > hohmann_total_normalized(R)


# --------------------------------------------------------------- classifier


@pytest.mark.parametrize(
    "R, expected",
    [
        (1.5, REGIME_HOHMANN_ONLY),
        (5.0, REGIME_HOHMANN_ONLY),
        (11.0, REGIME_HOHMANN_ONLY),
        (11.93, REGIME_HOHMANN_ONLY),
        (11.95, REGIME_LARGE_B_BIELLIPTIC),
        (13.0, REGIME_LARGE_B_BIELLIPTIC),
        (15.5, REGIME_LARGE_B_BIELLIPTIC),
        (15.6, REGIME_ALL_BIELLIPTIC),
        (20.0, REGIME_ALL_BIELLIPTIC),
        (1000.0, REGIME_ALL_BIELLIPTIC),
    ],
)
def test_classifier(R: float, expected: str) -> None:
    assert classify_radius_ratio(R) == expected


def test_classifier_on_both_sides_of_each_threshold() -> None:
    r1_star, r2_star = threshold_R1(), threshold_R2()
    assert classify_radius_ratio(r1_star * (1.0 - 1e-6)) == REGIME_HOHMANN_ONLY
    assert classify_radius_ratio(r1_star * (1.0 + 1e-6)) == REGIME_LARGE_B_BIELLIPTIC
    assert classify_radius_ratio(r2_star * (1.0 - 1e-6)) == REGIME_LARGE_B_BIELLIPTIC
    assert classify_radius_ratio(r2_star * (1.0 + 1e-6)) == REGIME_ALL_BIELLIPTIC


def test_classifier_exactly_at_thresholds_is_well_defined() -> None:
    assert classify_radius_ratio(threshold_R1()) == REGIME_HOHMANN_ONLY
    assert classify_radius_ratio(threshold_R2()) == REGIME_LARGE_B_BIELLIPTIC


def test_classifier_rejects_R_at_or_below_one() -> None:
    for bad in (1.0, 0.5, -3.0):
        with pytest.raises(ValueError, match="R > 1"):
            classify_radius_ratio(bad)


# --------------------------------------------------------------- break-even B

# Production values, recomputed by the solver itself (M1 quoted Bcrit(12) ~ 815.8).
EXPECTED_BCRIT = {
    12.0: 815.8202504753092,
    13.0: 48.90484332838889,
    14.0: 26.10461128235042,
    15.0: 18.19028151222189,
}


@pytest.mark.parametrize("R, expected", sorted(EXPECTED_BCRIT.items()))
def test_break_even_B_regression(R: float, expected: float) -> None:
    result = break_even_B(R)
    assert result.regime == REGIME_LARGE_B_BIELLIPTIC
    assert result.B_crit == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize("R", sorted(EXPECTED_BCRIT))
def test_equality_holds_at_break_even_B(R: float) -> None:
    b_crit = break_even_B(R).B_crit
    assert b_crit is not None
    assert bielliptic_total_normalized(R, b_crit) == pytest.approx(
        hohmann_total_normalized(R), abs=1e-15
    )


@pytest.mark.parametrize("R", sorted(EXPECTED_BCRIT))
def test_sign_behaviour_around_break_even_B(R: float) -> None:
    """Below B_crit Hohmann wins; above it the bi-elliptic wins."""
    b_crit = break_even_B(R).B_crit
    assert b_crit is not None
    dv_h = hohmann_total_normalized(R)
    assert bielliptic_total_normalized(R, b_crit * 0.99) > dv_h
    assert bielliptic_total_normalized(R, b_crit * 1.01) < dv_h


@pytest.mark.parametrize("R", sorted(EXPECTED_BCRIT))
def test_break_even_B_is_not_the_trivial_root(R: float) -> None:
    """B = R always satisfies equality by construction and must never be returned."""
    result = break_even_B(R)
    assert result.B_crit is not None
    assert result.B_crit > R * 1.01
    assert bielliptic_total_normalized(R, R) == pytest.approx(
        hohmann_total_normalized(R), abs=1e-15
    )


@pytest.mark.parametrize("R", [1.5, 2.0, 5.0, 10.0, 11.5, 11.9])
def test_region_A_reports_no_bielliptic_win(R: float) -> None:
    result = break_even_B(R)
    assert result.regime == REGIME_HOHMANN_ONLY
    assert result.B_crit is None
    assert result.winning_B_infimum is None


@pytest.mark.parametrize("R", [16.0, 20.0, 50.0, 100.0, 1000.0])
def test_region_C_collapses_to_the_B_to_R_boundary(R: float) -> None:
    """No finite break-even root is manufactured in Region C."""
    result = break_even_B(R)
    assert result.regime == REGIME_ALL_BIELLIPTIC
    assert result.B_crit is None
    assert result.winning_B_infimum == R
    # every B > R really does win
    dv_h = hohmann_total_normalized(R)
    for factor in (1.000001, 1.01, 2.0, 100.0):
        assert bielliptic_total_normalized(R, R * factor) < dv_h


def test_break_even_regime_matches_classifier() -> None:
    for R in (2.0, 8.0, 11.9, 12.0, 13.0, 15.0, 16.0, 30.0, 200.0):
        assert break_even_B(R).regime == classify_radius_ratio(R)


def test_break_even_B_rejects_R_at_or_below_one() -> None:
    with pytest.raises(ValueError, match="R > 1"):
        break_even_B(1.0)


def test_break_even_B_diverges_approaching_R1_from_above() -> None:
    """B_crit -> infinity as R -> R1*+, and -> R as R -> R2*-."""
    r1_star, r2_star = threshold_R1(), threshold_R2()
    near_r1 = break_even_B(r1_star * (1.0 + 1e-4)).B_crit
    mid = break_even_B(0.5 * (r1_star + r2_star)).B_crit
    near_r2 = break_even_B(r2_star * (1.0 - 1e-6)).B_crit
    assert near_r1 is not None and mid is not None and near_r2 is not None
    assert near_r1 > mid > near_r2
    assert near_r2 == pytest.approx(r2_star, rel=1e-3)


# --------------------------------------------------------------- infimum


@pytest.mark.parametrize("R", [2.0, 5.0, 11.0, 12.0, 15.0, 16.0, 50.0])
def test_infimum_is_an_endpoint_value(R: float) -> None:
    value, endpoint = bielliptic_infimum_normalized(R)
    dv_h = hohmann_total_normalized(R)
    dv_inf = bielliptic_infinite_limit_normalized(R)
    assert value == pytest.approx(min(dv_h, dv_inf), rel=1e-15)
    assert endpoint == ("B->R+" if dv_h <= dv_inf else "B->inf")


@pytest.mark.parametrize("R", [2.0, 5.0, 11.0, 12.0, 15.0, 16.0, 50.0])
def test_no_finite_B_beats_the_infimum(R: float) -> None:
    infimum, _ = bielliptic_infimum_normalized(R)
    for exponent in range(0, 13):
        B = R * (1.0 + 10.0**-6) * 10.0**exponent
        assert bielliptic_total_normalized(R, B) >= infimum - 1e-15


# --------------------------------------------------------------- continuity


@pytest.mark.parametrize("R", [2.0, 12.0, 16.0, 50.0])
def test_curve_is_finite_over_twelve_decades_of_B(R: float) -> None:
    """No NaNs, infinities or sign pathologies across a huge dynamic range."""
    for exponent in range(0, 121):
        value = bielliptic_total_normalized(R, R * 10.0 ** (exponent / 10.0))
        assert math.isfinite(value)
        assert 0.0 < value < 1.0


@pytest.mark.parametrize("R", [2.0, 12.0, 16.0, 50.0])
def test_curve_is_continuous_in_B(R: float) -> None:
    """Continuity via refinement: halving the grid step must roughly halve the
    largest jump between neighbouring samples. A genuine discontinuity would
    leave the largest jump stuck at a finite size instead."""

    def largest_jump(n_samples: int) -> float:
        values = [
            bielliptic_total_normalized(R, R * (1.0 + 1e-9) * 10.0 ** (12.0 * i / n_samples))
            for i in range(n_samples + 1)
        ]
        return max(abs(b - a) for a, b in itertools.pairwise(values))

    coarse = largest_jump(500)
    fine = largest_jump(1000)
    finer = largest_jump(2000)
    assert fine < 0.6 * coarse
    assert finer < 0.6 * fine
    assert finer < 1e-2


def test_hohmann_curve_is_continuous_in_R() -> None:
    """Same refinement argument for the Hohmann curve in R."""

    def largest_jump(n_samples: int) -> float:
        values = [
            hohmann_total_normalized(1.0 + 10.0 ** (6.0 * i / n_samples - 3.0))
            for i in range(n_samples + 1)
        ]
        return max(abs(b - a) for a, b in itertools.pairwise(values))

    coarse = largest_jump(500)
    fine = largest_jump(1000)
    assert fine < 0.6 * coarse
    assert fine < 1e-2
    assert all(
        math.isfinite(hohmann_total_normalized(1.0 + 10.0 ** (e / 10.0 - 4.0)))
        for e in range(0, 81)
    )
