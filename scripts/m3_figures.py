#!/usr/bin/env python3
"""Generate the Milestone 3 portfolio figures.

    figures/m3_crossover_map.png            PRIMARY project figure
    figures/m3_dv_time_trade.png            practical delta-v / time story
    figures/m3_break_even_time_penalty.png  Region-B divergence detail

All plotted quantities come from the production trade module, which in turn
calls only the verified M2 transfer functions. Thresholds are taken from the M2
root solvers, never read off a plotting grid.

``B -> infinity`` is drawn only as a labelled asymptote or an arrow direction.
It is never plotted as a finite point, because it is a mathematical infimum with
unbounded transfer time, not a realizable transfer.

Usage (from the repository root)::

    python scripts/m3_figures.py
"""

from __future__ import annotations

import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402

from bielliptic_crossover import (  # noqa: E402
    R1_REFERENCE,
    V1_REFERENCE,
    bielliptic_infinite_limit_normalized,
    bielliptic_total_normalized,
    break_even_B,
    hohmann_total_normalized,
    threshold_R1,
    threshold_R2,
)
from bielliptic_crossover.trade import (  # noqa: E402
    FINITE_B_OVER_R_FAMILY,
    LUNAR_DISTANCE_KM,
    authoritative_R_grid,
    break_even_locus,
    break_even_trade,
    minimum_endpoint_transfer,
    practical_trade_metrics,
)

FIGURES = pathlib.Path("figures")
FIG_MAP = FIGURES / "m3_crossover_map.png"
FIG_TRADE = FIGURES / "m3_dv_time_trade.png"
FIG_PENALTY = FIGURES / "m3_break_even_time_penalty.png"

COLOUR_HOHMANN = "#1f5fa9"
COLOUR_INFIMUM = "#b03030"
COLOUR_SAVING = "#2e7d32"
COLOUR_BCRIT = "#7b3fa0"

SHADE_A = "#dfe9f5"
SHADE_B = "#fdf0d5"
SHADE_C = "#e3f2e3"


def _shade_regions(ax, r1_star: float, r2_star: float, r_min: float, r_max: float) -> None:
    ax.axvspan(r_min, r1_star, color=SHADE_A, zorder=0)
    ax.axvspan(r1_star, r2_star, color=SHADE_B, zorder=0)
    ax.axvspan(r2_star, r_max, color=SHADE_C, zorder=0)
    ax.axvline(r1_star, color="0.35", lw=1.3, ls="--", zorder=2)
    ax.axvline(r2_star, color="0.35", lw=1.3, ls="-.", zorder=2)


def figure_crossover_map(output_path: pathlib.Path) -> None:
    r1_star, r2_star = threshold_R1(), threshold_R2()
    grid = authoritative_R_grid()
    r_min, r_max = float(grid.min()), float(grid.max())

    dv_hohmann = np.array([hohmann_total_normalized(float(R)) for R in grid])
    dv_infimum = np.array([bielliptic_infinite_limit_normalized(float(R)) for R in grid])
    saving_ms = np.array([minimum_endpoint_transfer(float(R)).saving_m_s for R in grid])

    locus = break_even_locus()
    locus_R = np.array([point.R for point in locus])
    locus_ratio = np.array([point.B_crit_over_R for point in locus])

    fig, (ax_dv, ax_save, ax_bcrit) = plt.subplots(
        3, 1, figsize=(11.0, 12.2), sharex=True, layout="constrained",
        gridspec_kw={"height_ratios": [1.15, 1.0, 1.0]},
    )

    # ---------------- panel 1: normalized delta-v -----------------------------
    _shade_regions(ax_dv, r1_star, r2_star, r_min, r_max)
    ax_dv.plot(grid, dv_hohmann, color=COLOUR_HOHMANN, lw=2.2,
               label=r"Hohmann  $\overline{\Delta v}_H(R)$  (a real transfer)")
    ax_dv.plot(grid, dv_infimum, color=COLOUR_INFIMUM, lw=2.2, ls=(0, (5, 2)),
               label=r"bi-elliptic $B\to\infty$ MATHEMATICAL INFIMUM  "
                     r"$(\sqrt{2}-1)(1+1/\sqrt{R})$")
    ax_dv.set_xscale("log")
    ax_dv.set_ylabel(r"normalized $\Delta v$   [units of $v_1=\sqrt{\mu/r_1}$]")
    # Full data range: the infimum curve reaches ~0.826 as R -> 1, and
    # clipping it would misrepresent the Region-A gap.
    ax_dv.set_ylim(0.0, 0.98)
    ax_dv.grid(True, which="major", alpha=0.28)
    ax_dv.grid(True, which="minor", alpha=0.10)
    curve_legend = ax_dv.legend(loc="upper right", fontsize=9.5, framealpha=0.96)
    ax_dv.add_artist(curve_legend)
    # Region key as a separate legend so it can never be occluded by the curves
    # or by the curve legend.
    ax_dv.legend(
        handles=[
            Patch(facecolor=SHADE_A, edgecolor="0.5",
                  label="REGION A  ($R<R_1^*$): Hohmann beats EVERY $B$"),
            Patch(facecolor=SHADE_B, edgecolor="0.5",
                  label="REGION B  ($R_1^*<R<R_2^*$): only $B>B_{crit}(R)$ beats Hohmann"),
            Patch(facecolor=SHADE_C, edgecolor="0.5",
                  label="REGION C  ($R>R_2^*$): every $B>R$ beats Hohmann"),
        ],
        loc="lower right", fontsize=9.5, framealpha=0.96,
    )
    ax_dv.set_title(
        "Bi-elliptic vs. Hohmann crossover map\n"
        "coplanar two-body impulsive transfer between circular orbits; "
        "result depends only on the RADIUS ratio $R=r_2/r_1$",
        fontsize=12.5,
    )

    # Inset: the curves cross in a narrow band, invisible at full scale.
    inset = ax_dv.inset_axes([0.09, 0.50, 0.40, 0.44])
    zoom = (grid >= 9.0) & (grid <= 22.0)
    _shade_regions(inset, r1_star, r2_star, 9.0, 22.0)
    inset.plot(grid[zoom], dv_hohmann[zoom], color=COLOUR_HOHMANN, lw=1.9)
    inset.plot(grid[zoom], dv_infimum[zoom], color=COLOUR_INFIMUM, lw=1.9, ls=(0, (5, 2)))
    inset.set_xlim(9.0, 22.0)
    inset.set_ylim(0.505, 0.556)
    inset.tick_params(labelsize=8)
    inset.set_title("zoom: where the curves cross", fontsize=9)
    inset.grid(True, alpha=0.25)
    ax_dv.indicate_inset_zoom(inset, edgecolor="0.45")

    # ---------------- panel 2: best achievable saving -------------------------
    _shade_regions(ax_save, r1_star, r2_star, r_min, r_max)
    ax_save.plot(grid, saving_ms, color=COLOUR_SAVING, lw=2.2,
                 label=r"$\overline{\Delta v}_H-\min(\overline{\Delta v}_H,"
                       r"\overline{\Delta v}_{B\infty})$, Earth ref.")
    ax_save.axhline(0.0, color="0.35", lw=1.1, ls=":")
    ax_save.set_ylabel("best possible $\\Delta v$ saving  [m/s]\n"
                       "(requires $B\\to\\infty$: infinite transfer time)")
    ax_save.set_ylim(-15.0, 360.0)
    ax_save.grid(True, which="major", alpha=0.28)
    ax_save.grid(True, which="minor", alpha=0.10)
    ax_save.legend(loc="upper left", fontsize=9.5, framealpha=0.95)
    ax_save.annotate(
        "exactly zero throughout Region A:\nno bi-elliptic transfer helps at all",
        xy=(4.0, 0.0), xytext=(1.35, 105.0), fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "0.35", "lw": 1.1},
        bbox={"boxstyle": "round,pad=0.4", "fc": "white", "ec": "0.6", "alpha": 0.95},
    )
    ax_save.annotate(
        f"only {minimum_endpoint_transfer(12.0).saving_m_s:.2f} m/s at $R=12$,\n"
        "and only as $B\\to\\infty$",
        xy=(12.0, minimum_endpoint_transfer(12.0).saving_m_s),
        xytext=(17.0, 40.0), fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "0.35", "lw": 1.1},
        bbox={"boxstyle": "round,pad=0.4", "fc": "white", "ec": "0.6", "alpha": 0.95},
    )

    # ---------------- panel 3: break-even locus -------------------------------
    _shade_regions(ax_bcrit, r1_star, r2_star, r_min, r_max)
    ax_bcrit.plot(locus_R, locus_ratio, color=COLOUR_BCRIT, lw=2.4,
                  label=r"$B_{crit}(R)/R = r_{b,crit}/r_2$  (solved root at each $R$)")
    ax_bcrit.axhline(1.0, color="0.35", lw=1.1, ls=":")
    ax_bcrit.set_yscale("log")
    ax_bcrit.set_xscale("log")
    ax_bcrit.set_xlim(r_min, r_max)
    ax_bcrit.set_ylim(0.6, 5e5)
    ax_bcrit.set_xlabel(r"final-to-initial RADIUS ratio  $R=r_2/r_1$   (logarithmic)")
    ax_bcrit.set_ylabel("break-even apoapsis ratio\n$B_{crit}/R$   (logarithmic)")
    ax_bcrit.grid(True, which="major", alpha=0.28)
    ax_bcrit.grid(True, which="minor", alpha=0.10)
    ax_bcrit.legend(loc="upper right", fontsize=9.5, framealpha=0.95)
    ax_bcrit.annotate(
        f"$B_{{crit}}\\to\\infty$ as $R\\to R_1^{{*+}}$\n"
        f"($R_1^*$ = {r1_star:.6f}: first $R$ where\nANY bi-elliptic can win)\n"
        "curve runs off the top; it is not clipped to a value",
        xy=(locus_R[0], locus_ratio[0]), xytext=(1.30, 1.1e4), fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "0.35", "lw": 1.1},
        bbox={"boxstyle": "round,pad=0.4", "fc": "white", "ec": "0.6", "alpha": 0.95},
    )
    ax_bcrit.annotate(
        f"$B_{{crit}}\\to R$ as $R\\to R_2^{{*-}}$\n"
        f"($R_2^*$ = {r2_star:.6f})\n"
        "beyond it no $B_{crit}$ is defined:\nevery $B>R$ already wins",
        xy=(locus_R[-1], locus_ratio[-1]), xytext=(20.0, 30.0), fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "0.35", "lw": 1.1},
        bbox={"boxstyle": "round,pad=0.4", "fc": "white", "ec": "0.6", "alpha": 0.95},
    )

    # ---------------- region labels and threshold callouts --------------------
    for ax in (ax_dv, ax_save, ax_bcrit):
        ax.set_xlim(r_min, r_max)
    fig.savefig(REPO_ROOT / output_path, dpi=165)
    plt.close(fig)


def _trade_curve(R: float, n: int = 900) -> tuple[np.ndarray, np.ndarray]:
    """(transfer time [days], delta-v saving [m/s]) along a dense B sweep."""
    ratios = np.logspace(np.log10(1.0 + 1e-6), np.log10(1.0e5), n)
    times = np.empty(n)
    savings = np.empty(n)
    for index, ratio in enumerate(ratios):
        point = practical_trade_metrics(R, R * ratio)
        times[index] = point.t_bielliptic_days
        savings[index] = point.dv_saving_m_s
    return times, savings


def figure_dv_time_trade(output_path: pathlib.Path) -> None:
    cases = (12.0, 16.0)
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 6.6), layout="constrained")
    y_limits = {12.0: (-60.0, 12.0), 16.0: (-30.0, 165.0)}

    for ax, R in zip(axes, cases):
        ax.set_ylim(*y_limits[R])
        times, savings = _trade_curve(R)
        infimum = minimum_endpoint_transfer(R)
        hohmann_days = practical_trade_metrics(R, R * 2.0).t_hohmann_days

        ax.plot(times, savings, color=COLOUR_HOHMANN, lw=2.0, zorder=3,
                label="bi-elliptic family (dense $B$ sweep)")
        ax.axhline(0.0, color="0.3", lw=1.2, ls="--", zorder=2)
        ax.axhline(infimum.saving_m_s, color=COLOUR_INFIMUM, lw=1.5, ls=":", zorder=2,
                   label=f"$B\\to\\infty$ infimum = {infimum.saving_m_s:.2f} m/s")
        ax.plot([hohmann_days], [0.0], marker="*", ms=17, color=COLOUR_SAVING,
                mec="black", mew=0.7, zorder=6,
                label=f"Hohmann ({hohmann_days:.2f} d, 0 m/s)")

        for ratio in FINITE_B_OVER_R_FAMILY:
            point = practical_trade_metrics(R, R * ratio)
            ax.plot(point.t_bielliptic_days, point.dv_saving_m_s, marker="o", ms=6.5,
                    mfc="white", mec=COLOUR_HOHMANN, mew=1.6, zorder=5)
            if ratio in (1.25, 2.0, 10.0, 100.0):
                ax.annotate(f"$B/R$={ratio:g}",
                            xy=(point.t_bielliptic_days, point.dv_saving_m_s),
                            xytext=(6, -13), textcoords="offset points", fontsize=8.5)

        break_even = break_even_trade(R)
        if break_even is not None:
            ax.plot(break_even.t_bielliptic_days, 0.0, marker="D", ms=9,
                    color=COLOUR_BCRIT, mec="black", mew=0.7, zorder=6,
                    label=f"$B_{{crit}}$ = {break_even.B:.4g} ({break_even.t_bielliptic_days:.0f} d, 0 m/s)")

        # Asymptotic DIRECTION only: never a finite plotted point for B -> infinity.
        ax.annotate(
            "", xy=(times[-1] * 0.85, infimum.saving_m_s),
            xytext=(times[-1] * 0.06, infimum.saving_m_s),
            arrowprops={"arrowstyle": "->", "color": COLOUR_INFIMUM, "lw": 1.5,
                        "ls": ":"}, zorder=4,
        )
        ax.text(times[-1] * 0.1, infimum.saving_m_s * 1.02 + 1.0,
                "$B\\to\\infty$: approached only as\ntransfer time $\\to\\infty$",
                fontsize=8.5, color=COLOUR_INFIMUM, va="bottom")

        ax.set_xscale("log")
        ax.set_xlabel("bi-elliptic transfer time  [days]   (logarithmic)")
        ax.grid(True, which="major", alpha=0.28)
        ax.grid(True, which="minor", alpha=0.10)
        ax.legend(loc="lower right", fontsize=8.5, framealpha=0.95)

        # Model-validity marker: B at which the apoapsis RADIUS equals the lunar
        # distance. Beyond it an Earth-only two-body model is not appropriate.
        b_lunar = LUNAR_DISTANCE_KM / R1_REFERENCE
        lunar_point = practical_trade_metrics(R, b_lunar)
        ax.axvline(lunar_point.t_bielliptic_days, color="0.4", lw=1.4, ls="-.", zorder=2)
        ax.axvspan(lunar_point.t_bielliptic_days, times[-1] * 3.0,
                   color="#f3d9d9", alpha=0.35, zorder=0)
        ax.text(
            lunar_point.t_bielliptic_days * 1.35, ax.get_ylim()[0],
            "$r_b$ = lunar distance\n→ Earth-only two-body\nmodel INVALID to the right",
            fontsize=8.5, color="0.2", va="bottom",
            bbox={"boxstyle": "round,pad=0.35", "fc": "white", "ec": "0.55", "alpha": 0.95},
        )

        regime = "Region B" if break_even is not None else "Region C"
        ax.set_title(f"$R$ = {R:g}  ({regime})", fontsize=12)

    axes[0].set_ylabel("$\\Delta v$ saving vs. Hohmann  [m/s]\n(positive = bi-elliptic cheaper)")

    fig.suptitle(
        "Delta-v vs. transfer-time trade   (Earth reference: $r_1$ = "
        f"{R1_REFERENCE:.4f} km, $v_1$ = {V1_REFERENCE:.4f} km/s)\n"
        "circles are representative FINITE apoapsis designs, not optima — "
        "no finite optimum in $B$ exists",
        fontsize=12,
    )
    fig.savefig(REPO_ROOT / output_path, dpi=165)
    plt.close(fig)


def figure_break_even_penalty(output_path: pathlib.Path) -> None:
    r1_star, r2_star = threshold_R1(), threshold_R2()
    locus = break_even_locus()
    R_values = np.array([point.R for point in locus])
    ratio = np.array([point.B_crit_over_R for point in locus])
    time_ratio = np.array([point.time_ratio_at_break_even for point in locus])
    lunar = np.array([point.rb_over_lunar_distance for point in locus])

    fig, ax = plt.subplots(figsize=(11.0, 6.4), layout="constrained")
    ax.plot(R_values - r1_star, ratio, color=COLOUR_BCRIT, lw=2.3,
            label=r"$B_{crit}/R$  (required apoapsis, relative to the target)")
    ax.plot(R_values - r1_star, time_ratio, color=COLOUR_HOHMANN, lw=2.3, ls=(0, (5, 2)),
            label=r"$t_B/t_H$ at break-even  (time cost for ZERO saving)")
    ax.plot(R_values - r1_star, lunar, color="#c07000", lw=1.9, ls=(0, (1, 1.4)),
            label=r"$r_{b,crit}$ / lunar distance  (Earth reference)")
    ax.axhline(1.0, color="0.35", lw=1.1, ls=":")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"distance above the lower threshold,  $R - R_1^*$   (logarithmic)")
    ax.set_ylabel("ratio   (logarithmic)")
    ax.grid(True, which="major", alpha=0.28)
    ax.grid(True, which="minor", alpha=0.10)
    ax.legend(loc="upper right", fontsize=9.5, framealpha=0.95)
    ax.set_title(
        "Region B: the cost of merely reaching break-even\n"
        f"$R_1^*$ = {r1_star:.9f},  $R_2^*$ = {r2_star:.9f};  "
        "all three ratios diverge as $R\\to R_1^{*+}$",
        fontsize=12,
    )
    # Where the break-even apoapsis first falls inside the lunar distance.
    from scipy.optimize import brentq

    r_lunar = brentq(
        lambda R: break_even_trade(R).rb_over_lunar_distance - 1.0,
        r1_star * 1.0001, r2_star * 0.9999, xtol=1e-12,
    )
    ax.axvline(r_lunar - r1_star, color="#c07000", lw=1.4, ls="-")
    ax.annotate(
        f"$r_{{b,crit}}$ reaches the lunar distance\nonly at $R$ = {r_lunar:.4f}.\n"
        "BELOW that, break-even apoapsis lies\nbeyond the Moon and the Earth-only\ntwo-body model does not apply.",
        xy=(r_lunar - r1_star, 1.0), xytext=(3.0e-3, 2.2e2), fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "#c07000", "lw": 1.3},
        bbox={"boxstyle": "round,pad=0.4", "fc": "white", "ec": "#c07000", "alpha": 0.96},
    )
    ax.text(
        0.015, 0.055,
        "Measured power laws on the asymptotic branch:\n"
        r"    $B_{crit}/R \propto (R-R_1^*)^{-1}$        (fitted slope $-1.00005$)"
        "\n"
        r"    $t_B/t_H \propto (R-R_1^*)^{-3/2}$      (fitted slope $-1.49983$)"
        "\n"
        "At break-even the $\\Delta v$ saving is exactly ZERO by definition; these are\n"
        "the apoapsis and duration required before any saving begins at all.",
        transform=ax.transAxes, ha="left", fontsize=9,
        bbox={"boxstyle": "round,pad=0.45", "fc": "white", "ec": "0.6", "alpha": 0.96},
    )
    fig.savefig(REPO_ROOT / output_path, dpi=165)
    plt.close(fig)


def main() -> int:
    (REPO_ROOT / FIGURES).mkdir(parents=True, exist_ok=True)
    figure_crossover_map(FIG_MAP)
    figure_dv_time_trade(FIG_TRADE)
    figure_break_even_penalty(FIG_PENALTY)
    print(f"wrote {FIG_MAP.as_posix()}")
    print(f"wrote {FIG_TRADE.as_posix()}")
    print(f"wrote {FIG_PENALTY.as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
