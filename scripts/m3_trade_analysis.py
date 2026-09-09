#!/usr/bin/env python3
"""Regenerate the Milestone 3 trade-study result artifacts.

Writes, all from production code in ``bielliptic_crossover.trade``:

    results/m3_radius_trade.csv         endpoint-infimum envelope over the grid
    results/m3_finite_b_strategies.csv  finite-B strategy metrics per R
    results/m3_break_even_curve.csv     B_crit(R) locus across Region B
    results/m3_earth_examples.csv       Earth dimensional reference cases
    results/m3_summary.json             headline numbers and recommendation

All outputs are deterministic: no timestamps, no wall-clock values, no
machine-specific paths, sorted JSON keys, and full-precision ``repr`` floats in
the CSVs so figures can be reproduced exactly.

Usage (from the repository root)::

    python scripts/m3_trade_analysis.py
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from bielliptic_crossover import (  # noqa: E402
    MU_EARTH,
    R1_REFERENCE,
    R_EARTH,
    V1_REFERENCE,
    bielliptic_infinite_limit_normalized,
    break_even_B,
    classify_radius_ratio,
    hohmann_total_normalized,
    hohmann_transfer_time,
    threshold_R1,
    threshold_R2,
)
from bielliptic_crossover.trade import (  # noqa: E402
    GEO_RADIUS_KM,
    LUNAR_DISTANCE_KM,
    authoritative_R_grid,
    break_even_locus,
    break_even_trade,
    evaluate_trade,
    minimum_endpoint_transfer,
    model_validity,
    practical_trade_metrics,
    recommendation_for,
    recommendation_summary,
    sweep_radius_ratios,
)

RESULTS = pathlib.Path("results")
CSV_RADIUS_TRADE = RESULTS / "m3_radius_trade.csv"
CSV_FINITE_B = RESULTS / "m3_finite_b_strategies.csv"
CSV_BREAK_EVEN = RESULTS / "m3_break_even_curve.csv"
CSV_EARTH = RESULTS / "m3_earth_examples.csv"
JSON_SUMMARY = RESULTS / "m3_summary.json"

DAY = 86400.0

#: Radius ratios used for the finite-B strategy table and the Earth examples.
STRATEGY_R_VALUES = (2.0, 5.0, 10.0, 11.5, 11.9, 12.0, 12.5, 13.0, 14.0,
                     15.0, 15.5, 16.0, 17.0, 20.0, 50.0, 100.0)

#: Earth dimensional reference cases (GEO is inserted separately).
EARTH_R_VALUES = (10.0, 12.0, 13.0, 15.0, 16.0, 20.0, 50.0)


def _write_csv(path: pathlib.Path, header: list[str], rows: list[list[object]]) -> None:
    target = REPO_ROOT / path
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _f(value: float) -> str:
    """Full-precision, round-trippable float formatting."""
    return repr(float(value))


def write_radius_trade() -> int:
    """Endpoint-infimum envelope over the authoritative grid."""
    grid = authoritative_R_grid()
    envelope = sweep_radius_ratios(grid)
    rows = [
        [
            _f(item.R),
            classify_radius_ratio(item.R),
            _f(item.dv_hohmann_bar),
            _f(item.dv_infinite_bar),
            _f(item.dv_best_bar),
            item.best_endpoint,
            _f(item.saving_bar),
            _f(item.saving_m_s),
            "yes" if item.attainable else "no",
        ]
        for item in envelope
    ]
    _write_csv(
        CSV_RADIUS_TRADE,
        [
            "R",
            "regime",
            "dv_hohmann_bar",
            "dv_bielliptic_infinite_bar",
            "dv_best_bar",
            "best_endpoint",
            "saving_bar",
            "saving_m_s",
            "best_endpoint_attainable",
        ],
        rows,
    )
    return len(rows)


def write_finite_b_strategies() -> int:
    rows: list[list[object]] = []
    for R in STRATEGY_R_VALUES:
        trade = evaluate_trade(R)
        for point in trade.points:
            rows.append(
                [
                    _f(point.R),
                    trade.regime,
                    point.label,
                    _f(point.B),
                    _f(point.B_over_R),
                    _f(point.dv_hohmann_bar),
                    _f(point.dv_bielliptic_bar),
                    _f(point.saving_bar),
                    _f(point.dv_saving_m_s),
                    _f(point.t_hohmann_days),
                    _f(point.t_bielliptic_days),
                    _f(point.time_ratio),
                    _f(point.extra_time_days),
                    _f(point.rb_km),
                    _f(point.rb_altitude_km),
                    _f(point.rb_over_lunar_distance),
                    point.model_validity_flag,
                ]
            )
    _write_csv(
        CSV_FINITE_B,
        [
            "R",
            "regime",
            "strategy",
            "B",
            "B_over_R",
            "dv_hohmann_bar",
            "dv_bielliptic_bar",
            "saving_bar",
            "dv_saving_m_s",
            "t_hohmann_days",
            "t_bielliptic_days",
            "time_ratio",
            "extra_time_days",
            "rb_km",
            "rb_altitude_km",
            "rb_over_lunar_distance",
            "model_validity_flag",
        ],
        rows,
    )
    return len(rows)


def write_break_even_curve() -> int:
    locus = break_even_locus()
    rows = [
        [
            _f(point.R),
            _f(point.B_crit),
            _f(point.B_crit_over_R),
            _f(point.rb_km),
            _f(point.rb_altitude_km),
            _f(point.rb_over_lunar_distance),
            _f(point.time_ratio_at_break_even),
            point.model_validity_flag,
        ]
        for point in locus
    ]
    _write_csv(
        CSV_BREAK_EVEN,
        [
            "R",
            "B_crit",
            "B_crit_over_R",
            "rb_km",
            "rb_altitude_km",
            "rb_over_lunar_distance",
            "time_ratio_at_break_even",
            "model_validity_flag",
        ],
        rows,
    )
    return len(rows)


def _earth_row(R: float, name: str) -> list[object]:
    trade = evaluate_trade(R)
    infimum = trade.infimum
    # Representative bounded design: B = 5R, a deliberately arbitrary finite
    # strategy, NOT an optimum (no finite optimum exists).
    representative = practical_trade_metrics(R, 5.0 * R, label="B/R=5")
    break_even = break_even_trade(R)
    flag, _ = model_validity(representative.rb_km)
    return [
        name,
        _f(R),
        trade.regime,
        _f(trade.r2_km),
        _f(trade.r2_altitude_km),
        _f(infimum.dv_hohmann_bar * V1_REFERENCE),
        _f(infimum.dv_infinite_bar * V1_REFERENCE),
        _f(infimum.saving_m_s),
        _f(representative.dv_bielliptic_m_s / 1000.0),
        _f(representative.dv_saving_m_s),
        _f(representative.rb_km),
        _f(representative.rb_over_lunar_distance),
        flag,
        _f(hohmann_transfer_time(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE) / DAY),
        _f(representative.t_bielliptic_days),
        _f(representative.time_ratio),
        "" if trade.B_crit is None else _f(trade.B_crit),
        "" if break_even is None else _f(break_even.t_bielliptic_days),
        "" if break_even is None else _f(break_even.time_ratio),
    ]


def write_earth_examples() -> int:
    geo_R = GEO_RADIUS_KM / R1_REFERENCE
    rows = [_earth_row(geo_R, "GEO")]
    rows.extend(_earth_row(R, f"R={R:g}") for R in EARTH_R_VALUES)
    _write_csv(
        CSV_EARTH,
        [
            "case",
            "R",
            "regime",
            "r2_km",
            "r2_altitude_km",
            "dv_hohmann_km_s",
            "dv_infinite_infimum_km_s",
            "infimum_saving_m_s",
            "dv_representative_B5R_km_s",
            "representative_B5R_saving_m_s",
            "representative_rb_km",
            "representative_rb_over_lunar_distance",
            "representative_model_validity_flag",
            "t_hohmann_days",
            "t_representative_B5R_days",
            "time_ratio_B5R",
            "B_crit",
            "t_at_break_even_days",
            "time_ratio_at_break_even",
        ],
        rows,
    )
    return len(rows)


def _trade_block(R: float) -> dict[str, object]:
    trade = evaluate_trade(R)
    break_even = break_even_trade(R)
    recommendation = recommendation_for(R)
    block: dict[str, object] = {
        "R": R,
        "regime": trade.regime,
        "r2_km": trade.r2_km,
        "r2_altitude_km": trade.r2_altitude_km,
        "dv_hohmann_bar": trade.infimum.dv_hohmann_bar,
        "dv_hohmann_km_s": trade.infimum.dv_hohmann_bar * V1_REFERENCE,
        "dv_infinite_infimum_bar": trade.infimum.dv_infinite_bar,
        "infimum_saving_m_s": trade.infimum.saving_m_s,
        "infimum_endpoint": trade.infimum.best_endpoint,
        "infimum_attainable": trade.infimum.attainable,
        "t_hohmann_days": hohmann_transfer_time(
            MU_EARTH, R1_REFERENCE, R * R1_REFERENCE
        )
        / DAY,
        "B_crit": trade.B_crit,
        "finite_strategies": [
            {
                "strategy": point.label,
                "B": point.B,
                "B_over_R": point.B_over_R,
                "dv_saving_m_s": point.dv_saving_m_s,
                "time_ratio": point.time_ratio,
                "extra_time_days": point.extra_time_days,
                "rb_km": point.rb_km,
                "rb_over_lunar_distance": point.rb_over_lunar_distance,
                "model_validity_flag": point.model_validity_flag,
            }
            for point in trade.points
        ],
        "recommendation": {
            "headline": recommendation.headline,
            "delta_v_verdict": recommendation.delta_v_verdict,
            "practical_caveat": recommendation.practical_caveat,
        },
    }
    if break_even is not None:
        block["break_even"] = {
            "B_crit": break_even.B,
            "B_crit_over_R": break_even.B_over_R,
            "rb_km": break_even.rb_km,
            "rb_over_lunar_distance": break_even.rb_over_lunar_distance,
            "t_bielliptic_days": break_even.t_bielliptic_days,
            "time_ratio": break_even.time_ratio,
            "dv_saving_m_s": break_even.dv_saving_m_s,
            "model_validity_flag": break_even.model_validity_flag,
            "note": "saving is zero at break-even by definition",
        }
    return block


def build_summary() -> dict[str, object]:
    r1_star = threshold_R1()
    r2_star = threshold_R2()
    geo_R = GEO_RADIUS_KM / R1_REFERENCE

    return {
        "R1_star": r1_star,
        "R1_star_definition": "dv_bar_H(R) = dv_bar_B_inf(R); first R at which ANY bi-elliptic can beat Hohmann",
        "R2_star": r2_star,
        "R2_star_definition": "d(dv_bar_B)/dB = 0 at B = R+; R above which EVERY admissible B > R beats Hohmann",
        "Earth_reference": {
            "mu_km3_s2": MU_EARTH,
            "R_earth_km": R_EARTH,
            "h1_km": R1_REFERENCE - R_EARTH,
            "r1_km": R1_REFERENCE,
            "v1_km_s": V1_REFERENCE,
            "lunar_distance_km": LUNAR_DISTANCE_KM,
            "r2_at_R1_star_km": r1_star * R1_REFERENCE,
            "r2_at_R2_star_km": r2_star * R1_REFERENCE,
        },
        "Bcrit_examples": {
            f"{R:g}": break_even_B(R).B_crit
            for R in (12.0, 12.5, 13.0, 14.0, 15.0, 15.5)
        },
        "Bcrit_asymptotics": {
            "near_R1_star": "B_crit -> infinity as R -> R1*+",
            "near_R2_star": "B_crit -> R as R -> R2*-",
            "monotonicity": "B_crit/R decreases strictly across Region B",
        },
        "GEO_classification": {
            "r_geo_km": GEO_RADIUS_KM,
            "R": geo_R,
            "regime": classify_radius_ratio(geo_R),
            "dv_hohmann_km_s": hohmann_total_normalized(geo_R) * V1_REFERENCE,
            "dv_infinite_infimum_km_s": bielliptic_infinite_limit_normalized(geo_R)
            * V1_REFERENCE,
            "infimum_saving_m_s": minimum_endpoint_transfer(geo_R).saving_m_s,
            "verdict": "Region A: Hohmann beats every admissible bi-elliptic transfer",
        },
        "R12_trade": _trade_block(12.0),
        "R16_trade": _trade_block(16.0),
        "recommendation_summary": recommendation_summary(),
        "scope": {
            "model": "idealized coplanar two-body impulsive transfer",
            "excluded": [
                "J2 and higher gravity harmonics",
                "third-body (lunar/solar) perturbations",
                "atmospheric drag",
                "finite-burn losses",
                "plane changes",
                "low thrust",
                "radiation and environment models",
                "launch and operational constraints",
                "mission-specific optimization",
            ],
        },
    }


def main() -> int:
    n_envelope = write_radius_trade()
    n_strategies = write_finite_b_strategies()
    n_locus = write_break_even_curve()
    n_earth = write_earth_examples()

    summary_path = REPO_ROOT / JSON_SUMMARY
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(build_summary(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    # Relative paths only: never print a machine-specific absolute path.
    print(f"wrote {CSV_RADIUS_TRADE.as_posix()} ({n_envelope} rows)")
    print(f"wrote {CSV_FINITE_B.as_posix()} ({n_strategies} rows)")
    print(f"wrote {CSV_BREAK_EVEN.as_posix()} ({n_locus} rows)")
    print(f"wrote {CSV_EARTH.as_posix()} ({n_earth} rows)")
    print(f"wrote {JSON_SUMMARY.as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
