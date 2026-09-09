#!/usr/bin/env python3
"""Generate the two Milestone 2 diagnostic figures.

These are VERIFICATION figures, not the final portfolio crossover plot (which
is Milestone 4). They exist to make the M2 structural results visually checkable.

Figure 1  figures/m2_fig1_excess_vs_B.png
    dv_B(R,B) - dv_H(R) against B/R on a logarithmic axis, for R spanning
    Regions A, B and C. Makes the region structure and the absence of any
    interior minimum directly visible.

Figure 2  figures/m2_fig2_time_trade_R12.png
    Delta-v and transfer time together for the Region-B case R = 12, showing
    how a few m/s of saving demands an enormous transfer duration.

Usage (from the repository root)::

    python scripts/m2_diagnostic_figures.py
"""

from __future__ import annotations

import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from bielliptic_crossover import (  # noqa: E402
    MU_EARTH,
    R1_REFERENCE,
    bielliptic_infinite_limit_normalized,
    bielliptic_total_normalized,
    bielliptic_transfer_time,
    break_even_B,
    hohmann_total_normalized,
    hohmann_transfer_time,
    threshold_R1,
    threshold_R2,
)

FIGURE_DIR = pathlib.Path("figures")
FIG1 = FIGURE_DIR / "m2_fig1_excess_vs_B.png"
FIG2 = FIGURE_DIR / "m2_fig2_time_trade_R12.png"
DAY = 86400.0

# B/R range: the upper end is kept inside the range where the tail is resolvable
# in double precision, so the plot shows structure rather than rounding noise.
B_RATIO_MIN = 1.0 + 1e-6
B_RATIO_MAX = 1.0e8


def figure_1(output_path: pathlib.Path) -> None:
    r1_star, r2_star = threshold_R1(), threshold_R2()
    cases = [
        (5.0, "A"),
        (10.0, "A"),
        (12.0, "B"),
        (13.0, "B"),
        (15.0, "B"),
        (16.0, "C"),
        (20.0, "C"),
        (50.0, "C"),
    ]
    colours = plt.get_cmap("viridis")(np.linspace(0.05, 0.92, len(cases)))
    ratio_max = 1.0e6
    ratios = np.logspace(np.log10(B_RATIO_MIN), np.log10(ratio_max), 1600)

    fig, ax = plt.subplots(figsize=(10.6, 7.4))
    extremes = []
    for (R, region), colour in zip(cases, colours):
        dv_hohmann = hohmann_total_normalized(R)
        excess = np.array(
            [bielliptic_total_normalized(R, R * q) - dv_hohmann for q in ratios]
        )
        extremes.append((float(excess.min()), float(excess.max())))
        ax.plot(ratios, excess, color=colour, lw=1.9, label=f"R = {R:g}  (Region {region})")

        result = break_even_B(R)
        if result.B_crit is not None:
            ax.plot(
                result.B_crit / R, 0.0, marker="o", ms=8, mfc="white",
                mec=colour, mew=2.0, zorder=6,
            )

    ax.axhline(0.0, color="0.25", lw=1.2, ls="--", zorder=1)

    # Symmetric-log y axis: no curve is clipped, and the near-zero structure
    # (where the crossings live) stays legible. Linear within +/- 1e-3.
    linear_threshold = 1.0e-3
    ax.set_yscale("symlog", linthresh=linear_threshold, linscale=1.1)
    lowest = min(low for low, _ in extremes)
    highest = max(high for _, high in extremes)
    ax.set_ylim(1.35 * lowest, 1.35 * highest)
    ax.set_xscale("log")
    ax.set_xlim(B_RATIO_MIN, ratio_max)

    ax.axhspan(-linear_threshold, linear_threshold, color="0.85", alpha=0.45, zorder=0)

    ax.set_xlabel(
        r"apoapsis ratio  $B/R = r_b / r_2$   (logarithmic; $B/R = 1$ means $r_b = r_2$)"
    )
    ax.set_ylabel(
        r"$\overline{\Delta v}_B(R,B) - \overline{\Delta v}_H(R)$   [units of $v_1$]"
        "\nsymlog scale, linear inside the shaded band"
    )
    ax.set_title(
        "M2 diagnostic 1: bi-elliptic excess over Hohmann vs. intermediate apoapsis\n"
        "above zero = Hohmann better;  below zero = bi-elliptic better    "
        f"($R_1^*$ = {r1_star:.4f},  $R_2^*$ = {r2_star:.4f})",
        fontsize=11,
    )
    ax.text(
        0.985,
        0.965,
        "Region A: never crosses zero\n"
        "Region B: crosses zero at a finite $B_{crit}$ (open circles)\n"
        "Region C: below zero immediately above $B = R$\n"
        "every interior turning point is a MAXIMUM, never a minimum",
        transform=ax.transAxes,
        fontsize=9,
        va="top",
        ha="right",
        bbox={"boxstyle": "round,pad=0.45", "fc": "white", "ec": "0.6", "alpha": 0.94},
    )
    ax.grid(True, which="major", alpha=0.30)
    ax.grid(True, which="minor", alpha=0.12)
    # Legend placed below the axes so it never occludes the flat Region-C tails.
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        fontsize=9,
        framealpha=0.94,
        ncol=4,
    )

    fig.tight_layout()
    fig.savefig(REPO_ROOT / output_path, dpi=170)
    plt.close(fig)


def figure_2(output_path: pathlib.Path) -> None:
    R = 12.0
    dv_hohmann = hohmann_total_normalized(R)
    result = break_even_B(R)
    b_crit = result.B_crit
    assert b_crit is not None
    limit = bielliptic_infinite_limit_normalized(R)

    ratios = np.logspace(np.log10(B_RATIO_MIN), np.log10(1.0e6), 1200)
    b_values = R * ratios

    v1_km_s = (MU_EARTH / R1_REFERENCE) ** 0.5
    saving_ms = np.array(
        [(dv_hohmann - bielliptic_total_normalized(R, b)) * v1_km_s * 1000.0 for b in b_values]
    )
    times_days = np.array(
        [
            bielliptic_transfer_time(
                MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, b * R1_REFERENCE
            )
            / DAY
            for b in b_values
        ]
    )
    t_hohmann = hohmann_transfer_time(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE) / DAY
    max_saving_ms = (dv_hohmann - limit) * v1_km_s * 1000.0

    fig, (ax_dv, ax_t) = plt.subplots(
        2, 1, figsize=(10.0, 7.8), sharex=True, layout="constrained"
    )

    ax_dv.plot(ratios, saving_ms, color="#1f5fa9", lw=2.0, label=r"$\Delta v$ saving vs. Hohmann")
    ax_dv.axhline(0.0, color="0.25", lw=1.2, ls="--")
    ax_dv.axhline(
        max_saving_ms,
        color="#b03030",
        lw=1.4,
        ls=":",
        label=f"$B\\to\\infty$ asymptote = {max_saving_ms:.2f} m/s (infinite time)",
    )
    ax_dv.axvline(b_crit / R, color="0.45", lw=1.3, ls="-.")
    ax_dv.set_ylabel("$\\Delta v$ saving  [m/s]\n(positive = bi-elliptic better)")
    # Full data range with padding: the curve dips well below zero before it ever
    # turns positive, and clipping that away would hide the real cost of a
    # moderate apoapsis in Region B.
    span = float(saving_ms.max() - saving_ms.min())
    ax_dv.set_ylim(float(saving_ms.min()) - 0.08 * span, float(saving_ms.max()) + 0.22 * span)
    ax_dv.annotate(
        f"the ENTIRE prize at this R is {max_saving_ms:.2f} m/s",
        xy=(2.0e4, max_saving_ms),
        xytext=(30.0, max_saving_ms - 0.30 * span),
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "#b03030", "lw": 1.2},
        bbox={"boxstyle": "round,pad=0.4", "fc": "white", "ec": "#b03030", "alpha": 0.94},
    )
    ax_dv.annotate(
        f"worst case {saving_ms.min():.1f} m/s\nat B/R = {ratios[int(saving_ms.argmin())]:.1f}",
        xy=(ratios[int(saving_ms.argmin())], float(saving_ms.min())),
        xytext=(1.6, float(saving_ms.min()) + 0.18 * span),
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "0.35", "lw": 1.2},
        bbox={"boxstyle": "round,pad=0.4", "fc": "white", "ec": "0.6", "alpha": 0.94},
    )
    ax_dv.set_title(
        f"M2 diagnostic 2: delta-v vs. time trade for R = {R:g} (Region B), "
        f"Earth reference $r_1$ = {R1_REFERENCE:.4f} km",
        fontsize=11,
    )
    ax_dv.grid(True, which="major", alpha=0.30)
    ax_dv.grid(True, which="minor", alpha=0.12)
    ax_dv.legend(loc="center right", fontsize=9, framealpha=0.92)

    ax_t.plot(ratios, times_days, color="#1f5fa9", lw=2.0, label="bi-elliptic transfer time")
    ax_t.axhline(
        t_hohmann, color="#2e7d32", lw=1.4, ls=":", label=f"Hohmann = {t_hohmann:.3f} d"
    )
    ax_t.axvline(b_crit / R, color="0.45", lw=1.3, ls="-.")
    ax_t.set_yscale("log")
    ax_t.set_xscale("log")
    ax_t.set_xlim(B_RATIO_MIN, 1.0e6)
    ax_t.set_xlabel(r"apoapsis ratio  $B/R = r_b/r_2$   (logarithmic)")
    ax_t.set_ylabel("transfer time  [days]\n(logarithmic)")
    ax_t.grid(True, which="major", alpha=0.30)
    ax_t.grid(True, which="minor", alpha=0.12)
    ax_t.legend(loc="upper left", fontsize=9, framealpha=0.92)

    t_at_crit = (
        bielliptic_transfer_time(
            MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, b_crit * R1_REFERENCE
        )
        / DAY
    )
    ax_t.annotate(
        f"break-even $B_{{crit}}/R$ = {b_crit / R:.1f}\n"
        f"zero saving, {t_at_crit:.0f} d\n"
        f"= {t_at_crit / t_hohmann:.0f}x the Hohmann time",
        xy=(b_crit / R, t_at_crit),
        xytext=(3.0, 4.0e3),
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "0.35", "lw": 1.2},
        bbox={"boxstyle": "round,pad=0.45", "fc": "white", "ec": "0.6", "alpha": 0.92},
    )

    fig.savefig(REPO_ROOT / output_path, dpi=170)
    plt.close(fig)


def main() -> int:
    (REPO_ROOT / FIGURE_DIR).mkdir(parents=True, exist_ok=True)
    figure_1(FIG1)
    figure_2(FIG2)
    # Relative paths only: never print a machine-specific absolute path.
    print(f"wrote {FIG1.as_posix()}")
    print(f"wrote {FIG2.as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
