#!/usr/bin/env python3
"""Generate the final numerical sanity table from production code.

Writes:

    results/m4_final_summary.md    markdown table, embedded verbatim in the docs
    results/m4_final_summary.csv   the same rows, machine readable

The markdown table is the single source of truth for the headline numbers quoted
in ``README.md`` and ``DESIGN.md``. It is generated, never hand-maintained, and
``tests/test_m4_docs.py`` asserts the documented values still match production --
this exists because an M4 audit found hand-transcribed prose tables had drifted
from the production output (see DESIGN.md, M4 audit section).

Deterministic: no timestamps, no wall-clock values, no machine-specific paths.

Usage (from the repository root)::

    python scripts/m4_final_summary.py
"""

from __future__ import annotations

import csv
import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from bielliptic_crossover import (
    R1_REFERENCE,
    R_EARTH,
    V1_REFERENCE,
    classify_radius_ratio,
)
from bielliptic_crossover.trade import (
    GEO_RADIUS_KM,
    break_even_trade,
    minimum_endpoint_transfer,
    practical_trade_metrics,
)

OUT_MD = pathlib.Path("results/m4_final_summary.md")
OUT_CSV = pathlib.Path("results/m4_final_summary.csv")

#: Representative finite design used in the headline table. Deliberately an
#: arbitrary bounded choice, NOT an optimum -- no finite optimum in B exists.
REPRESENTATIVE_B_OVER_R = 2.0

CASES: tuple[tuple[str, float | None], ...] = (
    ("GEO", None),  # R derived from the GEO radius
    ("R=10", 10.0),
    ("R=12", 12.0),
    ("R=13", 13.0),
    ("R=15", 15.0),
    ("R=16", 16.0),
    ("R=20", 20.0),
    ("R=50", 50.0),
)

REGION_LABEL = {
    "hohmann_only": "A",
    "large_B_bielliptic": "B",
    "all_bielliptic": "C",
}

FLAG_LABEL = {
    "two_body_reasonable": "ok",
    "lunar_third_body_caution": "caution",
    "two_body_invalid": "INVALID",
}

HEADER = [
    "case",
    "R",
    "region",
    "dv_hohmann_km_s",
    "dv_infimum_km_s",
    "math_saving_m_s",
    "representative_B_over_R",
    "finite_B_saving_m_s",
    "time_ratio",
    "model_validity",
]


def rows() -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for name, value in CASES:
        R = GEO_RADIUS_KM / R1_REFERENCE if value is None else value
        infimum = minimum_endpoint_transfer(R)
        point = practical_trade_metrics(R, REPRESENTATIVE_B_OVER_R * R)
        out.append(
            {
                "case": name,
                "R": R,
                "region": REGION_LABEL[classify_radius_ratio(R)],
                "dv_hohmann_km_s": infimum.dv_hohmann_bar * V1_REFERENCE,
                "dv_infimum_km_s": infimum.dv_infinite_bar * V1_REFERENCE,
                "math_saving_m_s": infimum.saving_m_s,
                "representative_B_over_R": REPRESENTATIVE_B_OVER_R,
                "finite_B_saving_m_s": point.dv_saving_m_s,
                "time_ratio": point.time_ratio,
                "model_validity": FLAG_LABEL[point.model_validity_flag],
            }
        )
    return out


def markdown_table() -> str:
    data = rows()
    lines = [
        "| case | `R` | region | Hohmann Δv [km/s] | `B→∞` infimum [km/s] | math. saving | `B/R` | finite-`B` saving | `t_B/t_H` | model |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in data:
        lines.append(
            f"| {row['case']} | {row['R']:.4f} | {row['region']} | "
            f"{row['dv_hohmann_km_s']:.4f} | {row['dv_infimum_km_s']:.4f} | "
            f"{row['math_saving_m_s']:+.2f} m/s | {row['representative_B_over_R']:g} | "
            f"{row['finite_B_saving_m_s']:+.2f} m/s | {row['time_ratio']:.2f} | "
            f"{row['model_validity']} |"
        )
    return "\n".join(lines)


def detail_tables() -> str:
    """Per-strategy tables for the two headline radius ratios."""
    blocks: list[str] = []
    for R in (12.0, 16.0):
        region = REGION_LABEL[classify_radius_ratio(R)]
        infimum = minimum_endpoint_transfer(R)
        blocks.append(
            f"\n### `R = {R:g}` (Region {region}) — "
            f"mathematical infimum saving {infimum.saving_m_s:+.2f} m/s "
            f"(requires `B → ∞`)\n"
        )
        blocks.append(
            "| strategy | `B` | saving [m/s] | `t_B/t_H` | `r_b` / lunar dist. | model |"
        )
        blocks.append("|---|---|---|---|---|---|")
        strategies: list[tuple[str, float]] = [
            (f"B/R={ratio:g}", R * ratio) for ratio in (1.25, 2.0, 5.0, 10.0, 100.0)
        ]
        break_even = break_even_trade(R)
        if break_even is not None:
            strategies.extend(
                [
                    ("B_crit", break_even.B),
                    ("2·B_crit", 2.0 * break_even.B),
                ]
            )
        for label, B in sorted(strategies, key=lambda item: item[1]):
            point = practical_trade_metrics(R, B)
            blocks.append(
                f"| `{label}` | {point.B:.2f} | {point.dv_saving_m_s:+.2f} | "
                f"{point.time_ratio:.2f} | {point.rb_over_lunar_distance:.2f} | "
                f"{FLAG_LABEL[point.model_validity_flag]} |"
            )
    return "\n".join(blocks)


def build_document() -> str:
    geo_R = GEO_RADIUS_KM / R1_REFERENCE
    parts = [
        "# Final numerical sanity table",
        "",
        "Generated by `scripts/m4_final_summary.py` from production code. Do not edit by hand.",
        "",
        f"Earth reference: `r1 = R_Earth + 300 km = {R1_REFERENCE:.4f} km`, "
        f"`v1 = {V1_REFERENCE:.14f} km/s`, `R_Earth = {R_EARTH:.4f} km`.",
        "",
        "`B/R = 2` is an arbitrary bounded design used for a like-for-like comparison,",
        "**not** an optimum: no finite optimum in `B` exists. The `B→∞` column is a",
        "**mathematical infimum** requiring unbounded transfer time, never a realizable",
        "transfer.",
        "",
        markdown_table(),
        "",
        f"GEO corresponds to `R = {geo_R:.9f}` and sits in Region A.",
        "",
        "## Per-strategy detail at the two headline radius ratios",
        detail_tables(),
        "",
    ]
    return "\n".join(parts)


def main() -> int:
    document = build_document()
    (REPO_ROOT / OUT_MD).parent.mkdir(parents=True, exist_ok=True)
    (REPO_ROOT / OUT_MD).write_text(document, encoding="utf-8")

    with (REPO_ROOT / OUT_CSV).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER, lineterminator="\n")
        writer.writeheader()
        for row in rows():
            writer.writerow(
                {
                    key: (repr(value) if isinstance(value, float) else value)
                    for key, value in row.items()
                }
            )

    print(f"wrote {OUT_MD.as_posix()}")
    print(f"wrote {OUT_CSV.as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
