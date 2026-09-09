"""Verification G and H: scaling, invariance, and dimensional round trips."""

from __future__ import annotations

import math

import pytest

from bielliptic_crossover.bielliptic import bielliptic_total_normalized
from bielliptic_crossover.constants import (
    MU_EARTH,
    MU_MARS,
    R1_REFERENCE,
    R_EARTH,
    R_MARS,
    V1_REFERENCE,
    altitude_from_radius,
    radius_from_altitude,
)
from bielliptic_crossover.crossover import threshold_R1, threshold_R2
from bielliptic_crossover.dimensional import (
    bielliptic_total_dv,
    circular_speed,
    hohmann_total_dv_direct,
    vis_viva_speed,
)
from bielliptic_crossover.hohmann import hohmann_total_normalized

# mu and r1 spanning many orders of magnitude.
BODIES = [
    (MU_EARTH, R1_REFERENCE),
    (MU_MARS, R_MARS + 400.0),
    (1.32712440018e11, 1.496e8),  # heliocentric
    (1.0e2, 1.0e1),
    (1.0e20, 1.0e12),
]
R_GRID = [2.0, 12.0, 15.0, 20.0, 50.0]


def test_earth_reference_constants() -> None:
    assert R1_REFERENCE == pytest.approx(6678.1363, abs=1e-9)
    assert V1_REFERENCE == pytest.approx(7.72576063698292, abs=1e-11)


def test_radius_altitude_round_trip() -> None:
    for h in (300.0, 35786.0, 100000.0):
        assert altitude_from_radius(radius_from_altitude(h)) == pytest.approx(h, abs=1e-9)
    assert radius_from_altitude(300.0) == pytest.approx(R_EARTH + 300.0, abs=1e-12)


def test_vis_viva_reduces_to_circular_speed() -> None:
    for r in (1.0e3, 6678.1363, 1.0e6):
        assert vis_viva_speed(MU_EARTH, r, r) == pytest.approx(
            circular_speed(MU_EARTH, r), rel=1e-15
        )


@pytest.mark.parametrize("mu, r1", BODIES)
@pytest.mark.parametrize("R", R_GRID)
def test_dimensional_dv_scales_as_sqrt_mu_over_r1(mu: float, r1: float, R: float) -> None:
    v1 = math.sqrt(mu / r1)
    assert hohmann_total_dv_direct(mu, r1, R * r1) / v1 == pytest.approx(
        hohmann_total_normalized(R), rel=1e-13
    )


@pytest.mark.parametrize("mu, r1", BODIES)
@pytest.mark.parametrize("R", R_GRID)
@pytest.mark.parametrize("factor", [1.000001, 2.0, 1e5])
def test_bielliptic_dimensional_round_trip(
    mu: float, r1: float, R: float, factor: float
) -> None:
    v1 = math.sqrt(mu / r1)
    direct = bielliptic_total_dv(mu, r1, R * r1, R * factor * r1)
    assert direct / v1 == pytest.approx(bielliptic_total_normalized(R, R * factor), rel=1e-13)


def test_explicit_sqrt_mu_over_r1_scaling_law() -> None:
    """Doubling mu multiplies dv by sqrt(2); quadrupling r1 halves it."""
    base = hohmann_total_dv_direct(MU_EARTH, R1_REFERENCE, 12.0 * R1_REFERENCE)
    scaled_mu = hohmann_total_dv_direct(2.0 * MU_EARTH, R1_REFERENCE, 12.0 * R1_REFERENCE)
    scaled_r1 = hohmann_total_dv_direct(
        MU_EARTH, 4.0 * R1_REFERENCE, 12.0 * 4.0 * R1_REFERENCE
    )
    assert scaled_mu / base == pytest.approx(math.sqrt(2.0), rel=1e-13)
    assert scaled_r1 / base == pytest.approx(0.5, rel=1e-13)


def test_crossover_thresholds_are_invariant_to_mu_and_r1() -> None:
    """The thresholds are pure radius ratios: they contain no mu or r1 at all.

    Recomputing them from dimensional transfers for wildly different bodies must
    give the same numbers.
    """
    from scipy.optimize import brentq

    from bielliptic_crossover.bielliptic import bielliptic_infinite_limit_normalized

    reference_R1 = threshold_R1()
    for mu, r1 in BODIES:
        v1 = math.sqrt(mu / r1)

        def excess(R: float, mu: float = mu, r1: float = r1, v1: float = v1) -> float:
            return (
                hohmann_total_dv_direct(mu, r1, R * r1) / v1
                - bielliptic_infinite_limit_normalized(R)
            )

        root = brentq(excess, 2.0, 100.0, xtol=1e-13, rtol=8.9e-16)
        assert root == pytest.approx(reference_R1, rel=1e-11)


def test_threshold_R2_is_dimensionless_by_construction() -> None:
    """R2* comes from a cubic with integer coefficients: no physical scale enters."""
    r2_star = threshold_R2()
    assert r2_star**3 - 15.0 * r2_star**2 - 9.0 * r2_star - 1.0 == pytest.approx(
        0.0, abs=1e-11
    )


def test_M1_earth_reference_threshold_radii() -> None:
    """M1 DESIGN.md Section 10 Earth reference conversions."""
    r2_at_R1 = threshold_R1() * R1_REFERENCE
    r2_at_R2 = threshold_R2() * R1_REFERENCE
    assert r2_at_R1 == pytest.approx(79728.7030801, abs=1e-4)
    assert altitude_from_radius(r2_at_R1) == pytest.approx(73350.5667801, abs=1e-4)
    assert r2_at_R2 == pytest.approx(104056.841526, abs=1e-3)
    assert altitude_from_radius(r2_at_R2) == pytest.approx(97678.7052257, abs=1e-3)


def test_geo_from_reference_parking_orbit_is_region_A() -> None:
    """GEO from 300 km is R ~ 6.31, far below R1*."""
    r_geo = 42164.0
    assert r_geo / R1_REFERENCE == pytest.approx(6.313737562, abs=1e-8)
    assert r_geo / R1_REFERENCE < threshold_R1()


def test_vis_viva_rejects_radius_outside_orbit() -> None:
    with pytest.raises(ValueError, match="outside"):
        vis_viva_speed(MU_EARTH, 1.0e6, 1.0e3)


def test_bielliptic_direct_rejects_rb_below_r2() -> None:
    with pytest.raises(ValueError, match="rb >= r2"):
        bielliptic_total_dv(MU_EARTH, R1_REFERENCE, 10.0 * R1_REFERENCE, 5.0 * R1_REFERENCE)
