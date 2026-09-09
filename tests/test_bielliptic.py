"""Verification B and C: bi-elliptic transfer and its infinite-B limit."""

from __future__ import annotations

import cmath
import math

import pytest

from bielliptic_crossover.bielliptic import (
    bielliptic_burns_normalized,
    bielliptic_infinite_limit_normalized,
    bielliptic_total_derivative_wrt_B,
    bielliptic_total_normalized,
    burn_directions,
)
from bielliptic_crossover.constants import MU_EARTH, R1_REFERENCE
from bielliptic_crossover.dimensional import bielliptic_burns_direct, bielliptic_total_dv
from bielliptic_crossover.hohmann import hohmann_burns_normalized, hohmann_total_normalized

R_GRID = [1.5, 2.0, 5.0, 10.0, 12.0, 15.0, 16.0, 20.0, 50.0, 200.0]
B_FACTORS = [1.0 + 1e-9, 1.000001, 1.5, 2.0, 5.0, 50.0, 1000.0, 1e6]


# ---------------------------------------------------------------- B = R identity


@pytest.mark.parametrize("R", R_GRID)
def test_burn3_is_exactly_zero_at_B_equals_R(R: float) -> None:
    assert bielliptic_burns_normalized(R, R)[2] == 0.0


@pytest.mark.parametrize("R", R_GRID)
def test_burn1_is_bitwise_identical_to_hohmann_at_B_equals_R(R: float) -> None:
    assert bielliptic_burns_normalized(R, R)[0] == hohmann_burns_normalized(R)[0]


@pytest.mark.parametrize("R", R_GRID)
def test_total_reduces_to_hohmann_at_B_equals_R(R: float) -> None:
    """Exact mathematically; in float, burn 2 is reached by two different but
    algebraically equal expressions, so agreement is to within ~1 ulp."""
    total = bielliptic_total_normalized(R, R)
    hohmann = hohmann_total_normalized(R)
    assert abs(total - hohmann) <= 2.0 * math.ulp(hohmann)


# ---------------------------------------------------------------- burn magnitudes


@pytest.mark.parametrize("R", R_GRID)
@pytest.mark.parametrize("factor", B_FACTORS)
def test_all_burn_magnitudes_positive(R: float, factor: float) -> None:
    dv1, dv2, dv3 = bielliptic_burns_normalized(R, R * factor)
    assert dv1 > 0.0
    assert dv2 > 0.0
    assert dv3 > 0.0


@pytest.mark.parametrize("R", R_GRID)
@pytest.mark.parametrize("factor", B_FACTORS)
def test_total_equals_sum_of_magnitudes(R: float, factor: float) -> None:
    """The total is the sum of the three magnitudes to within rounding.

    Compared at a few ulp rather than bitwise: CPython's ``sum`` applies
    compensated (Neumaier) summation to floats, so it can differ from the
    plain left-to-right addition used in the production function by 1 ulp.
    """
    burns = bielliptic_burns_normalized(R, R * factor)
    total = bielliptic_total_normalized(R, R * factor)
    assert abs(total - sum(burns)) <= 4.0 * math.ulp(total)


def test_burn_directions_documented() -> None:
    assert burn_directions() == ("prograde", "prograde", "retrograde")


@pytest.mark.parametrize("R", R_GRID)
@pytest.mark.parametrize("factor", [1.000001, 2.0, 50.0])
def test_burn3_is_physically_retrograde(R: float, factor: float) -> None:
    """Arrival at r2 (periapsis of ellipse 2) is FASTER than circular, so the
    third burn must decelerate. Checked from actual speeds, not assumed."""
    burns = bielliptic_burns_direct(
        MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, R * factor * R1_REFERENCE
    )
    assert burns.v_peri_ellipse2 > burns.v_circ_2
    assert burns.dv3_signed < 0.0
    # ... while the first two are prograde
    assert burns.dv1_signed > 0.0
    assert burns.dv2_signed > 0.0


# ---------------------------------------------------------------- dimensional path


@pytest.mark.parametrize("R", R_GRID)
@pytest.mark.parametrize("factor", B_FACTORS)
def test_normalized_matches_independent_direct_vis_viva(R: float, factor: float) -> None:
    v1 = math.sqrt(MU_EARTH / R1_REFERENCE)
    direct = bielliptic_total_dv(
        MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, R * factor * R1_REFERENCE
    )
    assert direct / v1 == pytest.approx(
        bielliptic_total_normalized(R, R * factor), rel=1e-13
    )


# ---------------------------------------------------------------- infinite limit


@pytest.mark.parametrize("R", R_GRID)
def test_infinite_limit_closed_form(R: float) -> None:
    assert bielliptic_infinite_limit_normalized(R) == pytest.approx(
        (math.sqrt(2.0) - 1.0) * (1.0 + 1.0 / math.sqrt(R)), rel=1e-15
    )


@pytest.mark.parametrize("R", [2.0, 12.0, 20.0, 50.0])
def test_finite_B_converges_to_infinite_limit(R: float) -> None:
    limit = bielliptic_infinite_limit_normalized(R)
    previous = abs(bielliptic_total_normalized(R, R * 1e3) - limit)
    for exponent in (5, 7, 9, 11):
        residual = abs(bielliptic_total_normalized(R, R * 10.0**exponent) - limit)
        assert residual < previous
        previous = residual
    assert previous < 1e-11


@pytest.mark.parametrize("R", [2.0, 4.0, 12.0, 16.0, 25.0])
def test_leading_order_asymptotic_correction(R: float) -> None:
    """M1 Section 4.3: dv_B - dv_B_inf ~ (sqrt(R) - 3) / (sqrt(2) * B)."""
    limit = bielliptic_infinite_limit_normalized(R)
    for B in (1e6, 1e8):
        actual = bielliptic_total_normalized(R, B) - limit
        predicted = (math.sqrt(R) - 3.0) / (math.sqrt(2.0) * B)
        assert actual == pytest.approx(predicted, rel=1e-4)


def test_asymptotic_sign_flips_at_R_equals_9() -> None:
    """Approach is from below for R < 9 and from above for R > 9."""
    limit_below = bielliptic_infinite_limit_normalized(4.0)
    limit_above = bielliptic_infinite_limit_normalized(16.0)
    assert bielliptic_total_normalized(4.0, 1e8) - limit_below < 0.0
    assert bielliptic_total_normalized(16.0, 1e8) - limit_above > 0.0


def test_infinite_limit_rejects_infinite_B_on_finite_function() -> None:
    with pytest.raises(ValueError, match="B must be finite"):
        bielliptic_total_normalized(12.0, float("inf"))


# ---------------------------------------------------------------- derivative


@pytest.mark.parametrize("R", [2.0, 10.0, 12.0, 15.0, 16.0, 20.0, 50.0])
@pytest.mark.parametrize("factor", [1.0, 1.5, 3.0, 20.0])
def test_analytic_derivative_matches_complex_step(R: float, factor: float) -> None:
    """Complex-step differentiation shares none of the analytic algebra and is
    free of subtractive cancellation."""

    def total_complex(B: complex) -> complex:
        return (
            cmath.sqrt(2.0 * B / (1.0 + B))
            - 1.0
            + cmath.sqrt(2.0 * R / (B * (R + B)))
            - cmath.sqrt(2.0 / (B * (1.0 + B)))
            + cmath.sqrt(2.0 * B / (R * (R + B)))
            - 1.0 / cmath.sqrt(R)
        )

    B = R * factor
    step = 1e-30
    complex_step = total_complex(complex(B, step)).imag / step
    analytic = bielliptic_total_derivative_wrt_B(R, B)
    assert analytic == pytest.approx(complex_step, rel=1e-12)


@pytest.mark.parametrize(
    "R, expected_sign",
    [(2.0, 1), (10.0, 1), (15.0, 1), (16.0, -1), (20.0, -1), (50.0, -1)],
)
def test_slope_sign_at_B_equals_R_flips_across_R2(R: float, expected_sign: int) -> None:
    slope = bielliptic_total_derivative_wrt_B(R, R)
    assert math.copysign(1.0, slope) == expected_sign


# ---------------------------------------------------------------- validation


def test_B_below_R_raises() -> None:
    with pytest.raises(ValueError, match="B >= R"):
        bielliptic_total_normalized(12.0, 11.9)


def test_R_below_one_raises() -> None:
    with pytest.raises(ValueError, match="R >= 1"):
        bielliptic_total_normalized(0.5, 10.0)


def test_M1_table_values_reproduced() -> None:
    """Verification J: M1 DESIGN.md Table 1 bi-elliptic columns."""
    expected = {
        2.0: (0.46632139, 0.60189809, 0.70710566, 0.70710678),
        5.0: (0.54094283, 0.57688720, 0.59945496, 0.59945550),
        10.0: (0.54261935, 0.54594466, 0.54519951, 0.54519939),
        12.0: (0.53923048, 0.53773598, 0.53378705, 0.53378672),
        15.0: (0.53385750, 0.52794811, 0.52116366, 0.52116304),
        16.0: (0.53211454, 0.52518686, 0.51776766, 0.51776695),
        20.0: (0.52563061, 0.51592650, 0.50683557, 0.50683453),
        50.0: (0.49665076, 0.48341568, 0.47279508, 0.47279221),
    }
    for R, (b2, b5, b_large, b_inf) in expected.items():
        assert bielliptic_total_normalized(R, 2 * R) == pytest.approx(b2, abs=6e-9)
        assert bielliptic_total_normalized(R, 5 * R) == pytest.approx(b5, abs=6e-9)
        assert bielliptic_total_normalized(R, 1e6) == pytest.approx(b_large, abs=6e-9)
        assert bielliptic_infinite_limit_normalized(R) == pytest.approx(b_inf, abs=6e-9)
