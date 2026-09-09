"""M4 audit tests: documentation must not drift from production output.

The M4 audit found that hand-transcribed prose tables in DESIGN.md had drifted
from the production numbers (production code, CSVs and figures were correct; five
table cells were not). These tests exist so that can never recur silently.
"""

from __future__ import annotations

import pathlib
import re

import pytest

import bielliptic_crossover
from bielliptic_crossover import (
    R1_REFERENCE,
    V1_REFERENCE,
    break_even_B,
    classify_radius_ratio,
    threshold_R1,
    threshold_R2,
)
from bielliptic_crossover.trade import (
    GEO_RADIUS_KM,
    break_even_trade,
    minimum_endpoint_transfer,
    practical_trade_metrics,
    recommendation_for,
)

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
README = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
DESIGN = (REPO_ROOT / "DESIGN.md").read_text(encoding="utf-8")
SUMMARY_MD = (REPO_ROOT / "results" / "m4_final_summary.md").read_text(encoding="utf-8")


# ------------------------------------------------------- version metadata


PYPROJECT = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")


def _pyproject_field(name: str) -> str:
    """Read a top-level [project] string field.

    Deliberately regex-based rather than using ``tomllib``: that module is
    standard library only from Python 3.11, and this package supports 3.10.
    """
    match = re.search(rf'^{name} = "([^"]+)"$', PYPROJECT, re.MULTILINE)
    assert match is not None, f"{name} not found in pyproject.toml"
    return match.group(1)


def test_version_metadata_is_consistent() -> None:
    assert _pyproject_field("version") == bielliptic_crossover.__version__
    assert bielliptic_crossover.__version__ == "1.0.0"
    assert f"Version `{bielliptic_crossover.__version__}`" in README


def test_license_is_present_and_declared() -> None:
    licence = (REPO_ROOT / "LICENSE").read_text(encoding="utf-8")
    assert "MIT License" in licence
    assert _pyproject_field("license") == "MIT"


def test_declared_python_floor_is_actually_supported() -> None:
    """Nothing in the package or its tests may need a newer stdlib than declared."""
    assert _pyproject_field("requires-python") == ">=3.10"
    for path in sorted((REPO_ROOT / "src").rglob("*.py")) + sorted(
        (REPO_ROOT / "tests").glob("*.py")
    ):
        text = path.read_text(encoding="utf-8")
        # Match a real import statement, not a mention of one in a string.
        assert not re.search(r"^\s*import tomllib\b", text, re.MULTILINE), (
            f"{path.name} imports tomllib, which needs Python 3.11+"
        )


def test_ci_workflow_runs_the_suite_on_both_python_versions() -> None:
    workflow = (REPO_ROOT / ".github" / "workflows" / "tests.yml").read_text(encoding="utf-8")
    assert "pytest -W error" in workflow
    assert '"3.11"' in workflow
    assert '"3.12"' in workflow


# ------------------------------------------------- README matches production


def test_readme_embeds_the_generated_summary_table() -> None:
    """The README's headline table is the generated one, verbatim."""
    table = SUMMARY_MD[SUMMARY_MD.index("| case |") : SUMMARY_MD.index("\n\nGEO corresponds")]
    assert table in README


def test_readme_thresholds_match_production() -> None:
    assert repr(threshold_R1()) in README
    assert repr(threshold_R2()) in README


def test_readme_earth_reference_matches_production() -> None:
    assert f"{R1_REFERENCE:.4f} km" in README
    assert f"{V1_REFERENCE:.14f} km/s" in README


def test_readme_embedded_figures_exist() -> None:
    for match in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", README):
        if match.startswith("http"):
            continue
        assert (REPO_ROOT / match).is_file(), match


def test_readme_relative_links_resolve() -> None:
    for target in re.findall(r"\]\((?!http)([^)#]+)\)", README):
        assert (REPO_ROOT / target).exists(), target


# ------------------------------------------- DESIGN prose tables match production

# (R, B, saving [m/s], time ratio) as printed in the DESIGN.md M3.8 tables.
DOCUMENTED_ROWS = [
    (12.0, 15.0, -22.30, 4.36),
    (12.0, 24.0, -39.02, 7.28),
    (12.0, 120.0, -14.99, 60.75),
    (16.0, 20.0, 5.97, 4.45),
    (16.0, 32.0, 31.87, 7.45),
    (16.0, 80.0, 85.39, 23.82),
    (16.0, 1600.0, 139.33, 1840.74),
]


@pytest.mark.parametrize("R, B, saving, time_ratio", DOCUMENTED_ROWS)
def test_documented_strategy_rows_match_production(
    R: float, B: float, saving: float, time_ratio: float
) -> None:
    point = practical_trade_metrics(R, B)
    assert point.dv_saving_m_s == pytest.approx(saving, abs=5e-3)
    assert point.time_ratio == pytest.approx(time_ratio, abs=5e-3)
    # ... and the printed value really is in DESIGN.md.
    assert f"{saving:+.2f}".lstrip("+") in DESIGN or f"{abs(saving):.2f}" in DESIGN


def test_documented_break_even_rows_match_production() -> None:
    break_even = break_even_trade(12.0)
    assert break_even is not None
    assert break_even.time_ratio == pytest.approx(1006.20, abs=5e-3)
    assert break_even.rb_over_lunar_distance == pytest.approx(14.17, abs=5e-3)
    doubled = practical_trade_metrics(12.0, 2.0 * break_even.B)
    assert doubled.dv_saving_m_s == pytest.approx(1.50, abs=5e-3)
    assert doubled.time_ratio == pytest.approx(2829.07, abs=5e-3)


def test_design_thresholds_and_headline_values_present() -> None:
    assert repr(threshold_R1()) in DESIGN
    assert repr(threshold_R2()) in DESIGN
    assert "12.834879147" in DESIGN  # lunar-distance crossover
    for value in ("815.8202504753092", "48.90484332838889",
                  "26.10461128235042", "18.19028151222189"):
        assert value in DESIGN


def test_design_has_a_status_banner_naming_M3_as_authoritative() -> None:
    head = DESIGN[:2000]
    assert "DESIGN — status" in head
    assert "authoritative" in head
    assert "M4" in head


# ------------------------------------------------------- claim / language audit

FORBIDDEN = [
    "bi-elliptic is always better above",
    "is where bi-elliptic first wins",
    "B→∞ is the practical optimum",
    "B -> infinity is optimal",
    "mission optimal",
    "flight optimal",
    "operationally optimal",
]


@pytest.mark.parametrize("phrase", FORBIDDEN)
@pytest.mark.parametrize("document", ["README.md", "DESIGN.md"])
def test_no_misleading_claims(document: str, phrase: str) -> None:
    text = (REPO_ROOT / document).read_text(encoding="utf-8").lower()
    assert phrase.lower() not in text


@pytest.mark.parametrize("document", ["README.md", "DESIGN.md"])
def test_infinite_B_is_always_qualified(document: str) -> None:
    text = (REPO_ROOT / document).read_text(encoding="utf-8").lower()
    assert "mathematical infimum" in text
    assert "never a realizable transfer" in text or "not a realizable transfer" in text


@pytest.mark.parametrize("document", ["README.md", "DESIGN.md"])
def test_both_thresholds_are_distinguished(document: str) -> None:
    text = (REPO_ROOT / document).read_text(encoding="utf-8")
    assert "11.9388" in text and "15.5817" in text
    assert "any" in text.lower() and "all" in text.lower()


# ------------------------------------------- Region C has no finite B_crit

M2_REPORT = (REPO_ROOT / "results" / "m2_verification_report.txt").read_text(
    encoding="utf-8"
)


def _transfer_time_rows() -> list[list[str]]:
    """Data rows of the generated report's TRANSFER TIMES table."""
    section = M2_REPORT[M2_REPORT.index("9. TRANSFER TIMES") :]
    section = section[: section.index("basis legend:")]
    rows = []
    for line in section.splitlines():
        fields = line.split()
        if len(fields) >= 7 and re.fullmatch(r"\d+\.\d+", fields[0]):
            rows.append(fields)
    return rows


def test_report_transfer_time_table_was_parsed() -> None:
    rows = _transfer_time_rows()
    assert len(rows) == 8
    assert [row[0] for row in rows] == [
        "2.0", "5.0", "10.0", "12.0", "15.0", "16.0", "20.0", "50.0"
    ]


def test_region_C_rows_are_never_labelled_as_break_even_roots() -> None:
    """Regression guard: in Region C no finite B_crit exists.

    The value tabulated there is the OPEN ``B -> R+`` boundary, evaluated only for
    the ``B = R`` timing degeneracy. Presenting it under a ``B_crit`` heading
    reads as 'the break-even apoapsis equals the target radius', which is false.
    """
    for row in _transfer_time_rows():
        R = float(row[0])
        region, basis = row[1], row[6]
        assert region == {"hohmann_only": "A", "large_B_bielliptic": "B",
                          "all_bielliptic": "C"}[classify_radius_ratio(R)]
        if region == "C":
            assert basis == "B=R+", f"R={R} is Region C but basis is {basis!r}"
            assert break_even_B(R).B_crit is None
            # The tabulated value is the open boundary, i.e. R itself.
            assert float(row[5]) == pytest.approx(R, rel=1e-12)
        elif region == "B":
            assert basis == "B_crit"
            b_crit = break_even_B(R).B_crit
            assert b_crit is not None
            assert float(row[5]) == pytest.approx(b_crit, rel=1e-6)
        else:
            assert basis == "-"
            assert row[5] == "none"
            assert break_even_B(R).B_crit is None


def test_report_states_region_C_has_no_finite_break_even() -> None:
    legend = M2_REPORT[M2_REPORT.index("basis legend:") :][:900]
    assert "NO finite B_crit exists" in legend
    assert "NOT a break-even root" in legend
    assert "open interval (R, inf)" in legend


def test_design_states_region_C_has_no_finite_break_even() -> None:
    assert "**none exists**" in DESIGN
    assert "no finite `B_crit` exists" in DESIGN
    # ... and never claims B_crit equals R beyond R2*.
    assert "| ≥ `R2*` | `= R` |" not in DESIGN


@pytest.mark.parametrize("R", [16.0, 20.0, 50.0, 100.0])
def test_production_reports_no_break_even_root_in_region_C(R: float) -> None:
    result = break_even_B(R)
    assert result.regime == "all_bielliptic"
    assert result.B_crit is None
    assert result.winning_B_infimum == R


# ------------------------------------------------------- artifacts


def test_no_absolute_paths_in_any_result_artifact() -> None:
    for path in sorted((REPO_ROOT / "results").glob("*")):
        text = path.read_text(encoding="utf-8")
        for fragment in ("/Users/", "/home/", "/private/", "C:\\"):
            assert fragment not in text, f"{path.name} contains {fragment}"


def test_final_summary_regenerates_deterministically() -> None:
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "m4_final_summary", REPO_ROOT / "scripts" / "m4_final_summary.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.build_document() == module.build_document()
    assert module.build_document() == SUMMARY_MD


# ------------------------------------------------------- recommendations


@pytest.mark.parametrize(
    "R, expected_regime",
    [
        (GEO_RADIUS_KM / R1_REFERENCE, "hohmann_only"),
        (5.0, "hohmann_only"),
        (12.0, "large_B_bielliptic"),
        (15.0, "large_B_bielliptic"),
        (16.0, "all_bielliptic"),
        (50.0, "all_bielliptic"),
    ],
)
def test_recommendation_regimes(R: float, expected_regime: str) -> None:
    recommendation = recommendation_for(R)
    assert recommendation.regime == expected_regime
    assert recommendation.regime == classify_radius_ratio(R)


def test_geo_recommendation_is_hohmann() -> None:
    geo_R = GEO_RADIUS_KM / R1_REFERENCE
    assert "Hohmann" in recommendation_for(geo_R).headline
    assert minimum_endpoint_transfer(geo_R).saving_bar == 0.0
