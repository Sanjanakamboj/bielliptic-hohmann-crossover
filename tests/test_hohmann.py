"""Verification A: Hohmann transfer."""

from __future__ import annotations

import math

import pytest

from bielliptic_crossover.constants import MU_EARTH, R1_REFERENCE
from bielliptic_crossover.dimensional import hohmann_burns_direct, hohmann_total_dv_direct
from bielliptic_crossover.hohmann import (
    hohmann_burns_normalized,
    hohmann_total_normalized,
    hohmann_total_dv,
)

R_GRID = [1.0001, 1.5, 2.0, 5.0, 10.0, 11.9, 12.0, 15.0, 15.6, 16.0, 20.0, 50.0, 200.0, 1e4]


def test_R_equals_one_gives_exactly_zero() -> None:
    assert hohmann_total_normalized(1.0) == 0.0
    assert hohmann_burns_normalized(1.0) == (0.0, 0.0)


def test_R_approaching_one_tends_to_zero_like_half_delta() -> None:
    """Leading behaviour is dv ~ (R-1)/2."""
    for d in (1e-3, 1e-5, 1e-7):
        dv = hohmann_total_normalized(1.0 + d)
        assert dv == pytest.approx(d / 2.0, rel=1e-3)


@pytest.mark.parametrize("R", R_GRID)
def test_both_burns_positive_for_R_above_one(R: float) -> None:
    dv1, dv2 = hohmann_burns_normalized(R)
    assert dv1 > 0.0
    assert dv2 > 0.0
    assert hohmann_total_normalized(R) == pytest.approx(dv1 + dv2, rel=0, abs=1e-18)


@pytest.mark.parametrize("R", R_GRID)
def test_both_burns_are_prograde(R: float) -> None:
    """Physically: transfer periapsis is faster than circular at r1, and
    transfer apoapsis is slower than circular at r2."""
    burns = hohmann_burns_direct(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE)
    assert burns.v_peri_transfer > burns.v_circ_1
    assert burns.v_apo_transfer < burns.v_circ_2
    assert burns.dv1_signed > 0.0
    assert burns.dv2_signed > 0.0


@pytest.mark.parametrize("R", R_GRID)
def test_normalized_matches_direct_vis_viva(R: float) -> None:
    v1 = math.sqrt(MU_EARTH / R1_REFERENCE)
    direct = hohmann_total_dv_direct(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE)
    assert direct / v1 == pytest.approx(hohmann_total_normalized(R), rel=1e-14)


@pytest.mark.parametrize("R", R_GRID)
def test_dimensional_wrapper_matches_normalized_times_v1(R: float) -> None:
    v1 = math.sqrt(MU_EARTH / R1_REFERENCE)
    wrapper = hohmann_total_dv(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE)
    assert wrapper == pytest.approx(hohmann_total_normalized(R) * v1, rel=1e-15)


@pytest.mark.parametrize("bad_R", [0.999, 0.0, -1.0, -0.5])
def test_invalid_R_below_one_raises(bad_R: float) -> None:
    with pytest.raises(ValueError, match="R >= 1"):
        hohmann_total_normalized(bad_R)


@pytest.mark.parametrize("bad_R", [float("nan"), float("inf")])
def test_nan_and_inf_R_raise(bad_R: float) -> None:
    with pytest.raises(ValueError):
        hohmann_total_normalized(bad_R)


def test_dimensional_wrapper_validates_inputs() -> None:
    with pytest.raises(ValueError, match="mu must be positive"):
        hohmann_total_dv(-1.0, 1.0, 2.0)
    with pytest.raises(ValueError, match="r1 must be positive"):
        hohmann_total_dv(1.0, 0.0, 2.0)
    with pytest.raises(ValueError, match="r2 >= r1"):
        hohmann_total_dv(1.0, 2.0, 1.0)


def test_M1_table_values_reproduced() -> None:
    """Verification J: M1 DESIGN.md Table 1, to its printed precision."""
    expected = {
        2.0: 0.28445705,
        5.0: 0.48000915,
        10.0: 0.52978752,
        12.0: 0.53417987,
        15.0: 0.53621819,
        16.0: 0.53623939,
        20.0: 0.53473136,
        50.0: 0.51369584,
    }
    for R, value in expected.items():
        assert hohmann_total_normalized(R) == pytest.approx(value, abs=6e-9)


def test_M1_earth_dimensional_values_reproduced() -> None:
    """M1 DESIGN.md Table 2, dv_H in km/s to printed precision."""
    expected = {2.0: 2.197647, 12.0: 4.126946, 16.0: 4.142857, 50.0: 3.968691}
    for R, value in expected.items():
        got = hohmann_total_dv(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE)
        assert got == pytest.approx(value, abs=5e-7)
