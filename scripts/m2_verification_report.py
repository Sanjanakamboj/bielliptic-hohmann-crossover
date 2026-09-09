#!/usr/bin/env python3
"""Regenerate the Milestone 2 verification report.

Recomputes every headline number from the production package and writes
``results/m2_verification_report.txt``.

The report is deterministic: it contains no timestamps, no wall-clock timings
and no machine-specific paths, so regenerating it twice must produce byte
identical output. That is checked by ``tests/test_report.py``.

Usage (from the repository root)::

    python scripts/m2_verification_report.py
"""

from __future__ import annotations

import math
import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from bielliptic_crossover import (  # noqa: E402
    MU_EARTH,
    R1_REFERENCE,
    R_EARTH,
    V1_REFERENCE,
    altitude_from_radius,
    bielliptic_burns_normalized,
    bielliptic_infinite_limit_normalized,
    bielliptic_total_derivative_wrt_B,
    bielliptic_total_dv,
    bielliptic_total_normalized,
    bielliptic_transfer_time,
    bielliptic_transfer_time_normalized,
    break_even_B,
    classify_radius_ratio,
    hohmann_burns_normalized,
    hohmann_total_dv_direct,
    hohmann_total_normalized,
    hohmann_transfer_time,
    hohmann_transfer_time_normalized,
    scan_B_structure,
    threshold_R1,
    threshold_R1_from_polynomial,
    threshold_R2,
    threshold_R2_from_cubic,
)

OUTPUT_PATH = pathlib.Path("results/m2_verification_report.txt")
DAY = 86400.0
R_TABLE = [2.0, 5.0, 10.0, 12.0, 15.0, 16.0, 20.0, 50.0]
R_BREAK_EVEN = [12.0, 13.0, 14.0, 15.0]
R_STRUCTURE = [2.0, 5.0, 10.0, 12.0, 14.0, 15.0, 16.0, 20.0, 50.0, 100.0]


def rule(title: str) -> str:
    return f"\n{'=' * 78}\n{title}\n{'=' * 78}"


def build_report() -> str:
    out: list[str] = []
    add = out.append

    add("MILESTONE 2 VERIFICATION REPORT")
    add("Bi-elliptic / Hohmann crossover")
    add("")
    add("Every value below is recomputed from the production package at run time.")
    add("No crossover constant is hardcoded. Regenerating this report must produce")
    add("byte-identical output.")

    # ---------------------------------------------------------------- thresholds
    add(rule("1. CROSSOVER THRESHOLDS (two independent formulations each)"))
    r1_transfer = threshold_R1()
    r1_poly = threshold_R1_from_polynomial()
    r2_slope = threshold_R2()
    r2_cubic = threshold_R2_from_cubic()
    add("")
    add("Threshold R1*  --  condition dv_H(R) = dv_B_inf(R), examined at B -> infinity")
    add("  asks: can ANY bi-elliptic transfer beat Hohmann?   (exists B)")
    add(f"  from transfer equations (brentq) : {r1_transfer:.16f}")
    add(f"  from cubic u^3-(1+2sqrt2)u^2+u+1=0  : {r1_poly:.16f}   (u = sqrt(R))")
    add(f"  absolute difference              : {abs(r1_transfer - r1_poly):.3e}")
    add("")
    add("Threshold R2*  --  condition d(dv_B)/dB = 0 at B = R+, examined at B -> R+")
    add("  asks: do ALL bi-elliptic transfers beat Hohmann?   (for all B)")
    add(f"  from slope condition (brentq)    : {r2_slope:.16f}")
    add(f"  from cubic R^3-15R^2-9R-1=0      : {r2_cubic:.16f}")
    add(f"  absolute difference              : {abs(r2_slope - r2_cubic):.3e}")
    add("")
    add(f"  separation R2* - R1*             : {r2_slope - r1_transfer:.16f}")
    add("")
    add("Residual checks at the recomputed roots:")
    add(
        f"  dv_H(R1*) - dv_B_inf(R1*)        : "
        f"{hohmann_total_normalized(r1_transfer) - bielliptic_infinite_limit_normalized(r1_transfer):+.3e}"
    )
    add(
        f"  d(dv_B)/dB at (R2*, R2*)         : "
        f"{bielliptic_total_derivative_wrt_B(r2_slope, r2_slope):+.3e}"
    )
    add(
        f"  R2*^3 - 15 R2*^2 - 9 R2* - 1     : "
        f"{r2_slope**3 - 15 * r2_slope**2 - 9 * r2_slope - 1:+.3e}"
    )

    # ---------------------------------------------------------------- earth ref
    add(rule("2. EARTH REFERENCE CASE"))
    add("")
    add(f"  mu_Earth                         : {MU_EARTH:.7f} km^3/s^2")
    add(f"  R_Earth (mean equatorial RADIUS) : {R_EARTH:.7f} km")
    add(f"  initial ALTITUDE h1              : {R1_REFERENCE - R_EARTH:.7f} km")
    add(f"  initial RADIUS r1                : {R1_REFERENCE:.7f} km")
    add(f"  initial circular speed v1        : {V1_REFERENCE:.14f} km/s")
    add("")
    add("  Threshold radii (RADIUS ratio converted to this reference case):")
    for label, ratio in (("R1*", r1_transfer), ("R2*", r2_slope)):
        r2_km = ratio * R1_REFERENCE
        add(
            f"    {label}: R = {ratio:.12f}  ->  r2 = {r2_km:.6f} km"
            f"  (altitude {altitude_from_radius(r2_km):.6f} km)"
        )
    add("")
    add(f"  GEO (r2 = 42164 km) corresponds to R = {42164.0 / R1_REFERENCE:.9f}")
    add(f"    classification: {classify_radius_ratio(42164.0 / R1_REFERENCE)}")

    # ---------------------------------------------------------------- normalized
    add(rule("3. NORMALIZED DELTA-V COMPARISON (units of v1)"))
    add("")
    add(
        f"{'R':>6} {'dv_H':>13} {'dv_B(B=2R)':>13} {'dv_B(B=5R)':>13} "
        f"{'dv_B(B=1e6)':>13} {'dv_B_inf':>13} {'better':>9} {'region':>20}"
    )
    add("-" * 110)
    for R in R_TABLE:
        dv_h = hohmann_total_normalized(R)
        dv_inf = bielliptic_infinite_limit_normalized(R)
        better = "Hohmann" if dv_h <= dv_inf else "bi-ell"
        add(
            f"{R:>6.1f} {dv_h:>13.8f} {bielliptic_total_normalized(R, 2 * R):>13.8f} "
            f"{bielliptic_total_normalized(R, 5 * R):>13.8f} "
            f"{bielliptic_total_normalized(R, 1e6):>13.8f} {dv_inf:>13.8f} "
            f"{better:>9} {classify_radius_ratio(R):>20}"
        )
    add("")
    add("Note: B = 2R and B = 5R are arbitrary samples, NOT optima. There is no")
    add("finite interior optimum in B (see section 6).")

    # ---------------------------------------------------------------- burns
    add(rule("4. BURN BREAKDOWN AND DIRECTIONS"))
    add("")
    add("Hohmann (both impulses prograde):")
    add(f"{'R':>6} {'dv1_bar':>13} {'dv2_bar':>13} {'total':>13}")
    add("-" * 48)
    for R in R_TABLE:
        dv1, dv2 = hohmann_burns_normalized(R)
        add(f"{R:>6.1f} {dv1:>13.8f} {dv2:>13.8f} {dv1 + dv2:>13.8f}")
    add("")
    add("Bi-elliptic at B = 5R (burns 1 and 2 prograde, burn 3 RETROGRADE;")
    add("magnitudes are positive by convention, direction is stated not signed):")
    add(f"{'R':>6} {'dv1_bar':>13} {'dv2_bar':>13} {'dv3_bar':>13} {'total':>13}")
    add("-" * 62)
    for R in R_TABLE:
        dv1, dv2, dv3 = bielliptic_burns_normalized(R, 5 * R)
        add(f"{R:>6.1f} {dv1:>13.8f} {dv2:>13.8f} {dv3:>13.8f} {dv1 + dv2 + dv3:>13.8f}")

    # ---------------------------------------------------------------- break-even
    add(rule("5. BREAK-EVEN APOAPSIS RATIO B_crit(R)"))
    add("")
    add("The trivial root at B = R (where equality holds by construction) is")
    add("explicitly excluded. Region A has no break-even B; Region C has none")
    add("either, because every B > R already wins and the boundary B = R is open.")
    add("")
    add(
        f"{'R':>6} {'regime':>20} {'B_crit':>20} {'rb_crit/r2':>12} "
        f"{'rb_crit alt [km]':>18} {'residual':>11}"
    )
    add("-" * 92)
    for R in [2.0, 5.0, 10.0, 11.9] + R_BREAK_EVEN + [15.5, 16.0, 20.0, 50.0]:
        result = break_even_B(R)
        if result.B_crit is None:
            add(f"{R:>6.1f} {result.regime:>20} {'-':>20} {'-':>12} {'-':>18} {'-':>11}")
        else:
            b_crit = result.B_crit
            residual = bielliptic_total_normalized(R, b_crit) - hohmann_total_normalized(R)
            altitude = altitude_from_radius(b_crit * R1_REFERENCE)
            add(
                f"{R:>6.1f} {result.regime:>20} {b_crit:>20.10f} {b_crit / R:>12.6f} "
                f"{altitude:>18.4e} {residual:>+11.2e}"
            )
    add("")
    add("Required regression targets (R = 12, 13, 14, 15), full production precision:")
    for R in R_BREAK_EVEN:
        add(f"  B_crit({R:.0f}) = {break_even_B(R).B_crit!r}")

    # ---------------------------------------------------------------- structure
    add(rule("6. STRUCTURAL VERIFICATION: NO FINITE INTERIOR MINIMUM"))
    add("")
    add("Dense logarithmic sweep in B for each R. Any interior stationary point")
    add("must be a MAXIMUM. The infimum is always an endpoint:")
    add("    inf over B of dv_B(R,B) = min( dv_H(R), dv_B_inf(R) )")
    add("")
    add(
        f"{'R':>6} {'samples':>9} {'B/R max':>10} {'interior min':>13} "
        f"{'interior max':>13} {'B at max':>13} {'infimum at':>12}"
    )
    add("-" * 82)
    for R in R_STRUCTURE:
        structure = scan_B_structure(R, n_samples=4001)
        b_at_max = (
            f"{structure.b_at_interior_maximum:.6f}"
            if structure.b_at_interior_maximum is not None
            else "-"
        )
        add(
            f"{R:>6.1f} {structure.n_samples:>9} {structure.b_ratio_max:>10.1e} "
            f"{structure.n_interior_minima:>13} {structure.n_interior_maxima:>13} "
            f"{b_at_max:>13} {structure.infimum_endpoint:>12}"
        )
    add("")
    add("Total interior minima found across all R above: ")
    total_minima = sum(
        scan_B_structure(R, n_samples=4001).n_interior_minima for R in R_STRUCTURE
    )
    add(f"    {total_minima}   (expected 0)")

    # ---------------------------------------------------------------- residuals
    add(rule("7. NORMALIZED vs INDEPENDENT DIMENSIONAL VIS-VIVA RESIDUALS"))
    add("")
    add("The dimensional path evaluates every apsis speed directly from vis-viva")
    add("and shares no algebra with the normalized closed forms.")
    add("")
    add(f"{'R':>6} {'B/R':>10} {'|direct/v1 - normalized|':>26} {'relative':>12}")
    add("-" * 58)
    worst_absolute = 0.0
    worst_relative = 0.0
    for R in R_TABLE:
        for factor in (1.000001, 2.0, 5.0, 1000.0):
            B = R * factor
            direct = bielliptic_total_dv(
                MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, B * R1_REFERENCE
            )
            normalized = bielliptic_total_normalized(R, B)
            absolute = abs(direct / V1_REFERENCE - normalized)
            worst_absolute = max(worst_absolute, absolute)
            worst_relative = max(worst_relative, absolute / normalized)
    for R in (2.0, 12.0, 50.0):
        for factor in (1.000001, 2.0, 1000.0):
            B = R * factor
            direct = bielliptic_total_dv(
                MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, B * R1_REFERENCE
            )
            normalized = bielliptic_total_normalized(R, B)
            absolute = abs(direct / V1_REFERENCE - normalized)
            add(f"{R:>6.1f} {factor:>10.6f} {absolute:>26.3e} {absolute / normalized:>12.3e}")
    add("")
    add(f"  worst absolute residual over the full grid : {worst_absolute:.3e}")
    add(f"  worst relative residual over the full grid : {worst_relative:.3e}")
    add("")
    worst_hohmann = max(
        abs(hohmann_total_dv_direct(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE) / V1_REFERENCE
            - hohmann_total_normalized(R))
        for R in R_TABLE
    )
    add(f"  worst Hohmann absolute residual            : {worst_hohmann:.3e}")

    # ---------------------------------------------------------------- dimensional
    add(rule("8. DIMENSIONAL EARTH EXAMPLES"))
    add("")
    add("Final ALTITUDE is a conversion from the universal RADIUS RATIO result:")
    add("    h2 = R * r1 - R_Earth")
    add("")
    add(
        f"{'R':>6} {'h2 [km]':>14} {'dv_H [km/s]':>13} {'dv_B_inf [km/s]':>16} "
        f"{'saving [m/s]':>13} {'saving [%]':>11}"
    )
    add("-" * 78)
    for R in R_TABLE:
        dv_h = hohmann_total_normalized(R) * V1_REFERENCE
        dv_inf = bielliptic_infinite_limit_normalized(R) * V1_REFERENCE
        add(
            f"{R:>6.1f} {R * R1_REFERENCE - R_EARTH:>14.4f} {dv_h:>13.6f} "
            f"{dv_inf:>16.6f} {(dv_h - dv_inf) * 1000.0:>13.4f} "
            f"{(dv_h - dv_inf) / dv_h * 100.0:>11.4f}"
        )
    add("")
    add("The saving column is the MAXIMUM theoretically available prize, reached")
    add("only as B -> infinity, i.e. only in infinite transfer time.")

    # ---------------------------------------------------------------- times
    add(rule("9. TRANSFER TIMES"))
    add("")
    add(
        f"{'R':>6} {'t_H [d]':>11} {'t_B(2R) [d]':>13} {'t_B(5R) [d]':>13} "
        f"{'B_crit':>14} {'t_B(B_crit) [d]':>17} {'t_B/t_H':>11}"
    )
    add("-" * 90)
    for R in R_TABLE:
        t_h = hohmann_transfer_time(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE) / DAY
        t_2r = (
            bielliptic_transfer_time(
                MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, 2 * R * R1_REFERENCE
            )
            / DAY
        )
        t_5r = (
            bielliptic_transfer_time(
                MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, 5 * R * R1_REFERENCE
            )
            / DAY
        )
        result = break_even_B(R)
        b_ref = result.B_crit if result.B_crit is not None else result.winning_B_infimum
        if b_ref is None:
            add(
                f"{R:>6.1f} {t_h:>11.6f} {t_2r:>13.6f} {t_5r:>13.6f} "
                f"{'none':>14} {'-':>17} {'-':>11}"
            )
        else:
            t_crit = (
                bielliptic_transfer_time(
                    MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, b_ref * R1_REFERENCE
                )
                / DAY
            )
            add(
                f"{R:>6.1f} {t_h:>11.6f} {t_2r:>13.6f} {t_5r:>13.6f} "
                f"{b_ref:>14.6f} {t_crit:>17.6f} {t_crit / t_h:>11.4f}"
            )
    add("")
    add("The B = R degeneracy (delta-v and time degenerate DIFFERENTLY):")
    add("  dv_B(R,R) = dv_H(R) exactly, but")
    add("  t_B(R,R) = t_H(R) + pi*R^(3/2) in units of t_star = sqrt(r1^3/mu),")
    add("  i.e. the Hohmann time plus HALF THE PERIOD of the final circular orbit,")
    add("  because ellipse 2 degenerates to that circular orbit and the path still")
    add("  coasts a half revolution of it.")
    add("")
    add(f"{'R':>6} {'t_bar_H':>14} {'t_bar_B(B=R)':>15} {'excess':>14} {'pi*R^1.5':>14}")
    add("-" * 66)
    for R in (2.0, 12.0, 16.0, 50.0):
        t_h_bar = hohmann_transfer_time_normalized(R)
        t_b_bar = bielliptic_transfer_time_normalized(R, R)
        add(
            f"{R:>6.1f} {t_h_bar:>14.8f} {t_b_bar:>15.8f} "
            f"{t_b_bar - t_h_bar:>14.8f} {math.pi * R**1.5:>14.8f}"
        )

    # ---------------------------------------------------------------- trade
    add(rule("10. THE REGION-B TRADE, QUANTIFIED (R = 12)"))
    add("")
    R = 12.0
    dv_h = hohmann_total_normalized(R) * V1_REFERENCE
    dv_inf = bielliptic_infinite_limit_normalized(R) * V1_REFERENCE
    b_crit = break_even_B(R).B_crit
    assert b_crit is not None
    t_h = hohmann_transfer_time(MU_EARTH, R1_REFERENCE, R * R1_REFERENCE) / DAY
    t_crit = (
        bielliptic_transfer_time(
            MU_EARTH, R1_REFERENCE, R * R1_REFERENCE, b_crit * R1_REFERENCE
        )
        / DAY
    )
    add(f"  region                                : {classify_radius_ratio(R)}")
    add(f"  entire theoretical prize (B -> inf)   : {(dv_h - dv_inf) * 1000.0:.4f} m/s")
    add(f"  break-even B_crit                     : {b_crit:.6f}")
    add(f"  break-even apoapsis RADIUS            : {b_crit * R1_REFERENCE:.4e} km")
    add(
        f"  break-even apoapsis ALTITUDE          : "
        f"{altitude_from_radius(b_crit * R1_REFERENCE):.4e} km"
    )
    add(f"  ... as a multiple of lunar distance   : {b_crit * R1_REFERENCE / 384400.0:.2f}")
    add(f"  Hohmann transfer time                 : {t_h:.6f} d")
    add(f"  transfer time merely to BREAK EVEN    : {t_crit:.4f} d")
    add(f"  time penalty factor                   : {t_crit / t_h:.2f}x")
    add("")
    add("  At break-even the delta-v saving is exactly zero by definition; the")
    add("  full prize above requires B -> infinity and therefore infinite time.")
    add("  This is why the deliverable must be an engineering trade, not a")
    add("  minimization. The recommendation itself is Milestone 6, not M2.")

    add("")
    add("=" * 78)
    add("END OF REPORT")
    add("=" * 78)
    return "\n".join(out) + "\n"


def main() -> int:
    report = build_report()
    output_path = REPO_ROOT / OUTPUT_PATH
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(report, encoding="utf-8")
    # Relative path only: never print a machine-specific absolute path.
    print(f"wrote {OUTPUT_PATH.as_posix()} ({len(report)} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
