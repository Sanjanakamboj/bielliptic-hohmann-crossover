"""M3 result artifacts must be reproducible and free of drift or absolute paths."""

from __future__ import annotations

import csv
import importlib.util
import itertools
import json
import pathlib

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "m3_trade_analysis.py"
RESULTS = REPO_ROOT / "results"

CSV_FILES = [
    "m3_radius_trade.csv",
    "m3_finite_b_strategies.csv",
    "m3_break_even_curve.csv",
    "m3_earth_examples.csv",
]


def _load_module():
    spec = importlib.util.spec_from_file_location("m3_trade_analysis", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_summary_build_is_deterministic() -> None:
    module = _load_module()
    first = json.dumps(module.build_summary(), indent=2, sort_keys=True)
    second = json.dumps(module.build_summary(), indent=2, sort_keys=True)
    assert first == second


def test_committed_summary_matches_regeneration() -> None:
    module = _load_module()
    expected = json.dumps(module.build_summary(), indent=2, sort_keys=True) + "\n"
    assert (RESULTS / "m3_summary.json").read_text(encoding="utf-8") == expected


@pytest.mark.parametrize("name", CSV_FILES + ["m3_summary.json"])
def test_artifact_exists_and_has_no_absolute_paths(name: str) -> None:
    text = (RESULTS / name).read_text(encoding="utf-8")
    assert text
    for fragment in ("/Users/", "/home/", "/private/", "C:\\"):
        assert fragment not in text


@pytest.mark.parametrize("name", CSV_FILES)
def test_csv_rows_are_round_trippable_floats(name: str) -> None:
    with (RESULTS / name).open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows
    for row in rows[:40]:
        for key, value in row.items():
            if key in {"regime", "strategy", "case", "model_validity_flag",
                       "best_endpoint", "best_endpoint_attainable",
                       "representative_model_validity_flag"}:
                continue
            if value == "":
                continue
            assert repr(float(value)) == value


def test_summary_contains_required_keys() -> None:
    summary = json.loads((RESULTS / "m3_summary.json").read_text(encoding="utf-8"))
    for key in (
        "R1_star",
        "R2_star",
        "Earth_reference",
        "Bcrit_examples",
        "GEO_classification",
        "R12_trade",
        "R16_trade",
        "recommendation_summary",
    ):
        assert key in summary


def test_summary_headline_values() -> None:
    summary = json.loads((RESULTS / "m3_summary.json").read_text(encoding="utf-8"))
    assert summary["R1_star"] == pytest.approx(11.93876547264589, rel=1e-13)
    assert summary["R2_star"] == pytest.approx(15.58171873876318, rel=1e-13)
    assert summary["Earth_reference"]["v1_km_s"] == pytest.approx(
        7.72576063698292, abs=1e-11
    )
    assert summary["Bcrit_examples"]["12"] == 815.8202504753092
    assert summary["Bcrit_examples"]["13"] == 48.90484332838889
    assert summary["Bcrit_examples"]["14"] == 26.10461128235042
    assert summary["Bcrit_examples"]["15"] == 18.19028151222189
    assert summary["GEO_classification"]["regime"] == "hohmann_only"
    assert summary["R12_trade"]["B_crit"] == 815.8202504753092
    assert summary["R16_trade"]["B_crit"] is None


def test_break_even_curve_is_monotone_and_inside_region_B() -> None:
    with (RESULTS / "m3_break_even_curve.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    ratios = [float(row["B_crit_over_R"]) for row in rows]
    assert all(later < earlier for earlier, later in itertools.pairwise(ratios))
    assert ratios[0] > 1e4
    assert ratios[-1] == pytest.approx(1.0, abs=1e-5)


def test_radius_trade_saving_is_zero_exactly_in_region_A() -> None:
    with (RESULTS / "m3_radius_trade.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    region_a = [row for row in rows if row["regime"] == "hohmann_only"]
    assert region_a
    assert all(float(row["saving_m_s"]) == 0.0 for row in region_a)
    assert all(row["best_endpoint"] == "B->R+" for row in region_a)
    others = [row for row in rows if row["regime"] != "hohmann_only"]
    assert all(float(row["saving_m_s"]) > 0.0 for row in others)
    assert all(row["best_endpoint_attainable"] == "no" for row in others)


def test_earth_examples_flag_geo_as_region_A() -> None:
    with (RESULTS / "m3_earth_examples.csv").open(encoding="utf-8") as handle:
        rows = {row["case"]: row for row in csv.DictReader(handle)}
    assert rows["GEO"]["regime"] == "hohmann_only"
    assert float(rows["GEO"]["infimum_saving_m_s"]) == 0.0
    assert float(rows["GEO"]["representative_B5R_saving_m_s"]) < 0.0
