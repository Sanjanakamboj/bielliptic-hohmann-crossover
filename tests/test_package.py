"""Package-level scaffold and API-surface tests.

At M1 this file asserted that no solver API existed yet. M2 implements that API,
so those assertions have been replaced by their M2 counterparts: the public
surface is now checked to exist and to be complete, and the *M3* surface is
checked to be absent. The M1 version remains in git history unchanged.
"""

from __future__ import annotations

import pathlib

import pytest

import bielliptic_crossover

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]


def test_package_imports() -> None:
    assert bielliptic_crossover is not None


def test_version_string() -> None:
    assert bielliptic_crossover.__version__ == "1.0.0"


@pytest.mark.parametrize(
    "relative_path",
    [
        "README.md",
        "DESIGN.md",
        ".gitignore",
        "pyproject.toml",
        "src/bielliptic_crossover/__init__.py",
        "src/bielliptic_crossover/constants.py",
        "src/bielliptic_crossover/_stable.py",
        "src/bielliptic_crossover/hohmann.py",
        "src/bielliptic_crossover/bielliptic.py",
        "src/bielliptic_crossover/timing.py",
        "src/bielliptic_crossover/dimensional.py",
        "src/bielliptic_crossover/crossover.py",
        "scripts/m2_verification_report.py",
        "scripts/m2_diagnostic_figures.py",
        "scripts/m3_trade_analysis.py",
        "scripts/m3_figures.py",
        "src/bielliptic_crossover/trade.py",
        "LICENSE",
        ".github/workflows/tests.yml",
        "scripts/m4_final_summary.py",
        "scripts/m4_verify_artifacts.py",
        "results/m2_verification_report.txt",
        "results/m4_final_summary.md",
        "results/m3_summary.json",
        "figures/m3_crossover_map.png",
        "figures/m3_dv_time_trade.png",
        "figures/m2_fig1_excess_vs_B.png",
        "figures/m2_fig2_time_trade_R12.png",
    ],
)
def test_scaffold_file_exists(relative_path: str) -> None:
    assert (REPO_ROOT / relative_path).is_file()


def test_src_layout_is_used() -> None:
    assert (REPO_ROOT / "src" / "bielliptic_crossover").is_dir()
    assert not (REPO_ROOT / "bielliptic_crossover").exists()


@pytest.mark.parametrize(
    "attribute_name",
    [
        "hohmann_total_normalized",
        "bielliptic_total_normalized",
        "bielliptic_infinite_limit_normalized",
        "threshold_R1",
        "threshold_R2",
        "break_even_B",
        "classify_radius_ratio",
        "scan_B_structure",
        "hohmann_transfer_time",
        "bielliptic_transfer_time",
        "bielliptic_total_dv",
    ],
)
def test_m2_public_api_exists(attribute_name: str) -> None:
    assert hasattr(bielliptic_crossover, attribute_name)
    assert attribute_name in bielliptic_crossover.__all__


@pytest.mark.parametrize(
    "attribute_name",
    [
        # Physics deliberately excluded from this idealized model, in EVERY
        # milestone. M3 delivered the sweeps, figures and recommendation (they
        # live in bielliptic_crossover.trade), but none of the following may
        # appear without an explicit scope change.
        "j2_correction",
        "j2_perturbation",
        "third_body_acceleration",
        "lunar_perturbation",
        "finite_burn_loss",
        "gravity_loss",
        "plane_change_dv",
        "inclination_change",
        "low_thrust_transfer",
        "lambert_solve",
        "radiation_dose",
        "mission_optimizer",
        "atmospheric_drag",
    ],
)
def test_excluded_physics_is_absent(attribute_name: str) -> None:
    """The model stays coplanar, two-body and impulsive."""
    import bielliptic_crossover.trade as trade_module

    assert not hasattr(bielliptic_crossover, attribute_name)
    assert not hasattr(trade_module, attribute_name)


def test_trade_layer_is_not_exported_from_the_M2_surface() -> None:
    """M3 is additive: the verified M2 public API in __init__ is unchanged."""
    assert "trade" not in bielliptic_crossover.__all__
    import bielliptic_crossover.trade as trade_module

    assert hasattr(trade_module, "evaluate_trade")
    assert hasattr(trade_module, "break_even_locus")
    assert hasattr(trade_module, "recommendation_summary")


def test_all_exports_resolve() -> None:
    for name in bielliptic_crossover.__all__:
        assert hasattr(bielliptic_crossover, name), name


def test_constants_module_holds_no_crossover_constants() -> None:
    """Thresholds must be recomputed, never stored as magic numbers."""
    source = (REPO_ROOT / "src" / "bielliptic_crossover" / "constants.py").read_text()
    for forbidden in ("11.93", "15.58", "11.9387", "15.5817"):
        assert forbidden not in source
