"""Verification I: the no-interior-minimum structure.

M1 established that ``dv_bar_B(R, .)`` has no interior local minimum in ``B``,
so the infimum over admissible ``B`` is always an endpoint. These tests
re-derive that numerically rather than assuming it, and are written so that they
would FAIL if a future implementation ever introduced a spurious finite interior
optimum.
"""

from __future__ import annotations

import math

import pytest

from bielliptic_crossover.bielliptic import (
    bielliptic_infinite_limit_normalized,
    bielliptic_total_derivative_wrt_B,
    bielliptic_total_normalized,
)
from bielliptic_crossover.crossover import (
    B_RATIO_RESOLVABLE_MAX,
    bielliptic_infimum_normalized,
    scan_B_structure,
    threshold_R2,
)
from bielliptic_crossover.hohmann import hohmann_total_normalized

# Spanning Region A, B and C, plus the R = 9 shape change.
R_SPAN = [2.0, 5.0, 8.0, 9.0, 10.0, 12.0, 14.0, 15.0, 16.0, 20.0, 50.0, 100.0]


@pytest.mark.parametrize("R", R_SPAN)
def test_no_interior_local_minimum(R: float) -> None:
    structure = scan_B_structure(R, n_samples=4001)
    assert structure.n_interior_minima == 0


@pytest.mark.parametrize("R", R_SPAN)
def test_at_most_one_interior_stationary_point_and_it_is_a_maximum(R: float) -> None:
    structure = scan_B_structure(R, n_samples=4001)
    assert structure.n_interior_maxima <= 1
    if structure.n_interior_maxima == 1:
        b_max = structure.b_at_interior_maximum
        assert b_max is not None
        peak = bielliptic_total_normalized(R, b_max)
        # A genuine maximum: strictly above neighbours on both sides. The lower
        # probe is placed midway back to B = R so that it stays in the domain
        # even when the peak sits close to R (as it does near R2*).
        lower = R + 0.5 * (b_max - R)
        upper = b_max * 1.1
        assert lower > R
        assert peak > bielliptic_total_normalized(R, lower)
        assert peak > bielliptic_total_normalized(R, upper)


@pytest.mark.parametrize("R", R_SPAN)
def test_infimum_is_an_endpoint_not_an_interior_point(R: float) -> None:
    infimum, _ = bielliptic_infimum_normalized(R)
    structure = scan_B_structure(R, n_samples=4001)
    assert infimum == pytest.approx(
        min(structure.dv_hohmann, structure.dv_infinite_limit), rel=1e-15
    )
    # No sampled interior point beats the endpoint infimum.
    assert structure.dv_at_B_near_R >= infimum - 1e-14
    assert structure.dv_at_B_max >= infimum - 1e-14


@pytest.mark.parametrize("R", R_SPAN)
def test_brute_force_minimum_sits_at_a_domain_endpoint(R: float) -> None:
    """Independent dense sweep: the smallest sampled value must be at one end."""
    n = 3000
    ratios = [
        10.0 ** (math.log10(B_RATIO_RESOLVABLE_MAX) * i / n) for i in range(n + 1)
    ]
    values = [bielliptic_total_normalized(R, R * (1.0 + 1e-10) * q) for q in ratios]
    argmin = min(range(len(values)), key=values.__getitem__)
    assert argmin in (0, len(values) - 1)


@pytest.mark.parametrize(
    "R, expected_shape",
    [
        (2.0, "increasing"),
        (5.0, "increasing"),
        (8.0, "increasing"),
        (10.0, "up_then_down"),
        (12.0, "up_then_down"),
        (15.0, "up_then_down"),
        (16.0, "decreasing"),
        (20.0, "decreasing"),
        (100.0, "decreasing"),
    ],
)
def test_shape_matches_the_M1_endpoint_slope_table(R: float, expected_shape: str) -> None:
    """M1 Section 6.1: the shape is set by the slope at B=R+ (flips at R2*) and
    the slope at large B (flips at R=9)."""
    slope_near = bielliptic_total_derivative_wrt_B(R, R)
    slope_far = bielliptic_total_derivative_wrt_B(R, R * 1e10)
    shape = {
        (True, True): "increasing",
        (True, False): "up_then_down",
        (False, False): "decreasing",
    }[(slope_near > 0.0, slope_far > 0.0)]
    assert shape == expected_shape


def test_large_B_slope_sign_flips_exactly_at_R_equals_9() -> None:
    assert bielliptic_total_derivative_wrt_B(8.9, 8.9 * 1e10) > 0.0
    assert bielliptic_total_derivative_wrt_B(9.1, 9.1 * 1e10) < 0.0
    assert abs(bielliptic_total_derivative_wrt_B(9.0, 9.0 * 1e12)) < 1e-24


def test_no_interior_minimum_over_a_wide_R_sweep() -> None:
    """Coarser per-R sampling, far wider R coverage."""
    R = 1.5
    while R <= 200.0:
        structure = scan_B_structure(R, n_samples=801)
        assert structure.n_interior_minima == 0, f"interior minimum found at R={R}"
        R = round(R + 2.5, 10)


def test_a_finite_B_win_forces_an_infinite_B_win() -> None:
    """M1 Section 6.3: the 'finite wins while infinite loses' region is empty.

    Wherever the infinite-B limit does NOT beat Hohmann, no finite B may either.
    """
    R = 1.5
    while R <= 60.0:
        dv_h = hohmann_total_normalized(R)
        if bielliptic_infinite_limit_normalized(R) >= dv_h:
            best = min(
                bielliptic_total_normalized(R, R * (1.0 + 1e-9) * 10.0 ** (10.0 * i / 400))
                for i in range(401)
            )
            assert best >= dv_h - 1e-15, f"finite-B win without infinite-B win at R={R}"
        R = round(R + 0.5, 10)


def test_large_B_tail_is_unresolvable_beyond_the_documented_limit() -> None:
    """Numerical-conditioning finding, recorded so it cannot be rediscovered as
    a 'structural' result.

    The large-B residual decays as ``(sqrt(R)-3)/(sqrt(2)B)``, so it eventually
    drops below a few ulp of the value itself and the curve goes flat to machine
    precision. With the noise floor switched off, a sweep past that point
    produces many spurious 'extrema' that are pure rounding noise.

    R = 9 is the worst case: the leading coefficient ``sqrt(R) - 3`` vanishes
    identically there, so the residual decays as ``1/B^2`` instead of ``1/B``
    and reaches the noise floor around ``B/R ~ 1e8`` -- well inside the range
    that is comfortably resolvable for a generic R. This is exactly why the
    detector needs a noise floor rather than a raw comparison.
    """
    noisy = scan_B_structure(9.0, n_samples=2001, b_ratio_max=1e14, noise_ulps=0.0)
    assert noisy.n_interior_minima > 50  # noise, not structure

    clean = scan_B_structure(9.0, n_samples=2001, b_ratio_max=1e14, noise_ulps=8.0)
    assert clean.n_interior_minima == 0

    # For a generic R the raw comparison is already clean within the documented
    # resolvable range, confirming the limit is about resolution and not about
    # the shape of the curve.
    for R in (12.0, 20.0):
        raw = scan_B_structure(
            R, n_samples=2001, b_ratio_max=B_RATIO_RESOLVABLE_MAX, noise_ulps=0.0
        )
        assert raw.n_interior_minima == 0


def test_R_equals_9_tail_decays_quadratically() -> None:
    """The 1/B term vanishes at R = 9, leaving a 1/B^2 tail."""
    limit = bielliptic_infinite_limit_normalized(9.0)
    previous = abs(bielliptic_total_normalized(9.0, 9.0 * 1e2) - limit)
    for exponent in (3, 4, 5, 6):
        residual = abs(bielliptic_total_normalized(9.0, 9.0 * 10.0**exponent) - limit)
        assert math.log10(previous / residual) == pytest.approx(2.0, abs=0.01)
        previous = residual


def test_residual_is_above_the_noise_floor_within_the_resolvable_range() -> None:
    for R in (12.0, 20.0, 100.0):
        limit = bielliptic_infinite_limit_normalized(R)
        residual = abs(
            bielliptic_total_normalized(R, R * B_RATIO_RESOLVABLE_MAX) - limit
        )
        assert residual > 100.0 * math.ulp(limit)


def test_scan_validates_its_arguments() -> None:
    with pytest.raises(ValueError, match="R > 1"):
        scan_B_structure(1.0)
    with pytest.raises(ValueError, match="n_samples"):
        scan_B_structure(12.0, n_samples=2)
    with pytest.raises(ValueError, match="b_ratio_max"):
        scan_B_structure(12.0, b_ratio_max=1.0)


def test_interior_maximum_location_regression() -> None:
    """M1 quoted B_max ~ 66.47 at R=10 and ~26.47 at R=12."""
    for R, expected in ((10.0, 66.4664), (12.0, 26.4660)):
        structure = scan_B_structure(R, n_samples=40001, b_ratio_max=1e6)
        assert structure.b_at_interior_maximum is not None
        assert structure.b_at_interior_maximum == pytest.approx(expected, rel=2e-3)


def test_shape_classification_is_exhaustive() -> None:
    """The (near-slope, far-slope) = (negative, positive) combination cannot occur."""
    R = 1.5
    while R <= 100.0:
        slope_near = bielliptic_total_derivative_wrt_B(R, R)
        slope_far = bielliptic_total_derivative_wrt_B(R, R * 1e10)
        assert not (slope_near < 0.0 and slope_far > 0.0), R
        R = round(R + 0.5, 10)


def test_R2_is_exactly_where_the_shape_changes() -> None:
    r2_star = threshold_R2()
    assert bielliptic_total_derivative_wrt_B(r2_star * (1 - 1e-9), r2_star * (1 - 1e-9)) > 0
    assert bielliptic_total_derivative_wrt_B(r2_star * (1 + 1e-9), r2_star * (1 + 1e-9)) < 0
    assert math.isclose(
        bielliptic_total_derivative_wrt_B(r2_star, r2_star), 0.0, abs_tol=1e-15
    )
