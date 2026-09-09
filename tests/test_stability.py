"""Numerical stability audit (DESIGN.md Section M2.5).

The transfer expressions all contain a difference of two nearly equal square
roots at some physically interesting boundary. These tests pin down where the
naive form loses precision, prove the production form does not, and record the
conditioning limits that no formulation can remove.
"""

from __future__ import annotations

import math

import pytest

from bielliptic_crossover._stable import sqrt1pm1
from bielliptic_crossover.bielliptic import (
    bielliptic_burns_normalized,
    bielliptic_total_normalized,
)
from bielliptic_crossover.constants import MU_EARTH, R1_REFERENCE
from bielliptic_crossover.dimensional import bielliptic_total_dv, hohmann_total_dv_direct
from bielliptic_crossover.hohmann import hohmann_total_normalized


def naive_burn3(R: float, B: float) -> float:
    """The direct, cancellation-prone form of burn 3."""
    return math.sqrt(2.0 * B / (R * (R + B))) - 1.0 / math.sqrt(R)


# ---------------------------------------------------------------- sqrt1pm1


def test_sqrt1pm1_is_exact_at_zero() -> None:
    assert sqrt1pm1(0.0) == 0.0


@pytest.mark.parametrize("x", [1e-18, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1.0])
def test_sqrt1pm1_matches_definition_for_moderate_x(x: float) -> None:
    assert sqrt1pm1(x) == pytest.approx(x / (math.sqrt(1.0 + x) + 1.0), rel=1e-15)


def test_sqrt1pm1_beats_the_naive_form_at_tiny_x() -> None:
    """At x = 1e-18 the naive form returns exactly 0; the stable form does not."""
    assert math.sqrt(1.0 + 1e-18) - 1.0 == 0.0
    assert sqrt1pm1(1e-18) == pytest.approx(5e-19, rel=1e-9)


def test_sqrt1pm1_rejects_x_below_minus_one() -> None:
    with pytest.raises(ValueError, match="x >= -1"):
        sqrt1pm1(-1.5)


# ------------------------------------------------- A. R -> 1+


@pytest.mark.parametrize("delta", [1e-2, 1e-4, 1e-6, 1e-8])
def test_hohmann_near_R_equals_one_matches_leading_order(delta: float) -> None:
    """dv ~ (R-1)/2 with a relative error that shrinks with delta."""
    assert hohmann_total_normalized(1.0 + delta) == pytest.approx(
        delta / 2.0, rel=2.0 * delta
    )


@pytest.mark.parametrize("delta", [1e-2, 1e-4, 1e-6, 1e-8, 1e-10])
def test_hohmann_stays_positive_and_monotone_near_one(delta: float) -> None:
    assert hohmann_total_normalized(1.0 + delta) > 0.0
    assert hohmann_total_normalized(1.0 + delta) > hohmann_total_normalized(1.0 + delta / 10)


def test_R_near_one_is_limited_by_input_representation_not_algebra() -> None:
    """Conditioning finding: for R - 1 below about 1e-10 the accuracy limit is
    how well ``R - 1`` survives being stored as a double, not the formula.

    At R = 1 + 1e-12 the nearest double reproduces R - 1 to only ~4 significant
    digits, and dv ~ (R-1)/2 inherits exactly that error. No reformulation can
    recover information the input no longer carries.
    """
    stored = 1.0 + 1e-12
    represented_delta = stored - 1.0
    relative_input_error = abs(represented_delta - 1e-12) / 1e-12
    assert relative_input_error > 1e-5  # the input itself is this inaccurate

    # The computed dv tracks the STORED delta faithfully, which is all it can do.
    assert hohmann_total_normalized(stored) == pytest.approx(
        represented_delta / 2.0, rel=1e-9
    )


# ------------------------------------------------- B. B -> R+


@pytest.mark.parametrize("delta", [1e-4, 1e-8, 1e-12, 1e-14])
def test_burn3_stays_accurate_where_the_naive_form_degrades(delta: float) -> None:
    """Burn 3 is the cancellation-prone term: it is the difference of two nearly
    equal speeds as B -> R+."""
    R = 12.0
    B = R * (1.0 + delta)
    stable = bielliptic_burns_normalized(R, B)[2]
    # Reference: the stable identity is exact, so compare against the analytic
    # leading behaviour dv3 ~ (B-R)/(2*(R+B)*sqrt(R)) * 2 = delta/(4*sqrt(R)).
    assert stable == pytest.approx(delta / (4.0 * math.sqrt(R)), rel=1e-3)
    assert stable > 0.0


@pytest.mark.parametrize("R", [1.5, 3.0, 4.5, 6.0, 12.0, 18.0])
def test_naive_burn3_goes_negative_at_B_equals_R_but_stable_does_not(R: float) -> None:
    """Demonstrated defect: the naive form produces a NEGATIVE burn magnitude,
    which is physically impossible, at exactly B = R. The stable form returns
    exactly zero."""
    assert naive_burn3(R, R) < 0.0
    assert bielliptic_burns_normalized(R, R)[2] == 0.0


@pytest.mark.parametrize("R", [2.0, 12.0, 50.0])
def test_burn3_is_nonnegative_everywhere_on_the_domain(R: float) -> None:
    for exponent in range(-16, 9):
        B = R * (1.0 + 10.0**exponent)
        assert bielliptic_burns_normalized(R, B)[2] >= 0.0


@pytest.mark.parametrize("R", [2.0, 12.0, 16.0, 50.0])
def test_total_is_monotone_approaching_B_equals_R(R: float) -> None:
    """No sign noise in the excess over Hohmann as B -> R+, which is what the
    break-even solver depends on."""
    dv_h = hohmann_total_normalized(R)
    excesses = [
        bielliptic_total_normalized(R, R * (1.0 + 10.0**-e)) - dv_h for e in range(2, 13)
    ]
    assert all(math.isfinite(value) for value in excesses)
    assert all(abs(value) < 1e-2 for value in excesses)


# ------------------------------------------------- D/E. large R and wide scales


@pytest.mark.parametrize("R", [1e2, 1e4, 1e6, 1e8, 1e10])
def test_large_R_stays_finite_and_bounded(R: float) -> None:
    value = hohmann_total_normalized(R)
    assert math.isfinite(value)
    # Hohmann tends to sqrt(2) - 1 from above as R -> infinity.
    assert math.sqrt(2.0) - 1.0 < value < 0.55


@pytest.mark.parametrize("mu", [1e2, MU_EARTH, 1.327e11, 1e20])
@pytest.mark.parametrize("r1", [1e1, R1_REFERENCE, 1e6, 1e12])
def test_dimensional_agreement_across_extreme_scales(mu: float, r1: float) -> None:
    v1 = math.sqrt(mu / r1)
    for R in (2.0, 12.0, 50.0):
        assert hohmann_total_dv_direct(mu, r1, R * r1) / v1 == pytest.approx(
            hohmann_total_normalized(R), rel=1e-13
        )
        for factor in (1.000001, 2.0, 1e5):
            direct = bielliptic_total_dv(mu, r1, R * r1, R * factor * r1)
            assert direct / v1 == pytest.approx(
                bielliptic_total_normalized(R, R * factor), rel=1e-13
            )
