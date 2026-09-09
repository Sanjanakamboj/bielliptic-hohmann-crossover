"""Verification H: transfer times, including the B = R degeneracy."""

from __future__ import annotations

import math

import pytest

from bielliptic_crossover.constants import MU_EARTH, R1_REFERENCE
from bielliptic_crossover.timing import (
    bielliptic_time_excess_at_B_equals_R,
    bielliptic_transfer_time,
    bielliptic_transfer_time_normalized,
    hohmann_transfer_time,
    hohmann_transfer_time_normalized,
    time_scale,
)

R_GRID = [2.0, 5.0, 10.0, 12.0, 16.0, 20.0, 50.0]
DAY = 86400.0


def test_time_scale_definition() -> None:
    assert time_scale(MU_EARTH, R1_REFERENCE) == pytest.approx(
        math.sqrt(R1_REFERENCE**3 / MU_EARTH), rel=1e-15
    )


@pytest.mark.parametrize("R", R_GRID)
def test_hohmann_time_normalized_matches_dimensional(R: float) -> None:
    t_star = time_scale(MU_EARTH, R1_REFERENCE)
    dimensional = hohmann_transfer_time(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE)
    assert dimensional / t_star == pytest.approx(
        hohmann_transfer_time_normalized(R), rel=1e-14
    )


@pytest.mark.parametrize("R", R_GRID)
@pytest.mark.parametrize("factor", [1.0, 1.5, 3.0, 100.0])
def test_bielliptic_time_normalized_matches_dimensional(R: float, factor: float) -> None:
    t_star = time_scale(MU_EARTH, R1_REFERENCE)
    dimensional = bielliptic_transfer_time(
        MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, R * factor * R1_REFERENCE
    )
    assert dimensional / t_star == pytest.approx(
        bielliptic_transfer_time_normalized(R, R * factor), rel=1e-14
    )


@pytest.mark.parametrize("R", R_GRID)
def test_hohmann_time_is_half_the_transfer_ellipse_period(R: float) -> None:
    """Independent formulation: half of 2*pi*sqrt(a^3/mu)."""
    a_h = (R1_REFERENCE + R * R1_REFERENCE) / 2.0
    half_period = 0.5 * 2.0 * math.pi * math.sqrt(a_h**3 / MU_EARTH)
    assert hohmann_transfer_time(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE) == pytest.approx(
        half_period, rel=1e-15
    )


@pytest.mark.parametrize("R", R_GRID)
def test_bielliptic_time_is_sum_of_two_half_periods(R: float) -> None:
    rb = 4.0 * R * R1_REFERENCE
    a_1 = (R1_REFERENCE + rb) / 2.0
    a_2 = (R * R1_REFERENCE + rb) / 2.0
    expected = 0.5 * 2.0 * math.pi * (
        math.sqrt(a_1**3 / MU_EARTH) + math.sqrt(a_2**3 / MU_EARTH)
    )
    assert bielliptic_transfer_time(
        MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, rb
    ) == pytest.approx(expected, rel=1e-15)


@pytest.mark.parametrize("R", R_GRID)
def test_time_increases_monotonically_with_B(R: float) -> None:
    previous = bielliptic_transfer_time_normalized(R, R)
    for factor in (1.0001, 1.1, 2.0, 10.0, 1e3, 1e6):
        current = bielliptic_transfer_time_normalized(R, R * factor)
        assert current > previous
        previous = current


@pytest.mark.parametrize("R", [2.0, 12.0, 50.0])
def test_time_grows_as_B_to_the_three_halves(R: float) -> None:
    """Doubling B multiplies the time by 2^(3/2) in the large-B limit."""
    big = 1e8 * R
    ratio = bielliptic_transfer_time_normalized(
        R, 2.0 * big
    ) / bielliptic_transfer_time_normalized(R, big)
    assert ratio == pytest.approx(2.0**1.5, rel=1e-6)


@pytest.mark.parametrize("R", R_GRID)
def test_time_is_unbounded_as_B_grows(R: float) -> None:
    assert bielliptic_transfer_time_normalized(R, R * 1e12) > 1e15


def test_infinite_B_is_rejected_rather_than_faked() -> None:
    with pytest.raises(ValueError, match="unbounded transfer time"):
        bielliptic_transfer_time_normalized(12.0, float("inf"))


# ------------------------------------------------- the B = R degeneracy


@pytest.mark.parametrize("R", R_GRID)
def test_B_equals_R_does_NOT_reduce_to_hohmann_time(R: float) -> None:
    """Delta-v degenerates to Hohmann at B = R, but time does NOT.

    Ellipse 2 becomes the circular target orbit and the path still coasts a
    half revolution of it, so the bi-elliptic time exceeds the Hohmann time by
    half the final circular period. This test exists to stop that subtlety from
    being silently 'fixed' into a false equality.
    """
    t_b = bielliptic_transfer_time_normalized(R, R)
    t_h = hohmann_transfer_time_normalized(R)
    assert t_b > t_h
    assert t_b - t_h == pytest.approx(math.pi * R**1.5, rel=1e-14)
    assert t_b - t_h == pytest.approx(bielliptic_time_excess_at_B_equals_R(R), rel=1e-14)


@pytest.mark.parametrize("R", R_GRID)
def test_time_excess_equals_half_final_circular_period(R: float) -> None:
    r2 = R * R1_REFERENCE
    half_period_at_r2 = math.pi * math.sqrt(r2**3 / MU_EARTH)
    t_star = time_scale(MU_EARTH, R1_REFERENCE)
    assert bielliptic_time_excess_at_B_equals_R(R) * t_star == pytest.approx(
        half_period_at_r2, rel=1e-14
    )


def test_M1_transfer_time_table_reproduced() -> None:
    """M1 DESIGN.md Table 3, days, to its printed precision."""
    expected_hohmann = {
        2.0: 0.0577414,
        5.0: 0.163317,
        10.0: 0.40541,
        12.0: 0.520859,
        16.0: 0.778894,
        20.0: 1.06939,
        50.0: 4.04725,
    }
    for R, days in expected_hohmann.items():
        got = hohmann_transfer_time(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE) / DAY
        assert got == pytest.approx(days, rel=1e-5)

    expected_2R = {2.0: 0.287557, 12.0: 3.7893, 16.0: 5.80202, 50.0: 31.6941}
    for R, days in expected_2R.items():
        got = (
            bielliptic_transfer_time(
                MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, 2.0 * R * R1_REFERENCE
            )
            / DAY
        )
        assert got == pytest.approx(days, rel=1e-5)
