"""The M2 verification report must be reproducible and free of drift."""

from __future__ import annotations

import importlib.util
import pathlib

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "m2_verification_report.py"
REPORT = REPO_ROOT / "results" / "m2_verification_report.txt"


def _load_report_module():
    spec = importlib.util.spec_from_file_location("m2_verification_report", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_report_generation_is_deterministic() -> None:
    module = _load_report_module()
    assert module.build_report() == module.build_report()


def test_committed_report_matches_regeneration() -> None:
    """Guards against a stale committed artifact."""
    module = _load_report_module()
    assert REPORT.read_text(encoding="utf-8") == module.build_report()


def test_report_contains_no_absolute_paths() -> None:
    text = REPORT.read_text(encoding="utf-8")
    for fragment in ("/Users/", "/home/", "/private/", "C:\\"):
        assert fragment not in text


@pytest.mark.parametrize(
    "expected",
    [
        "11.9387654726",
        "15.5817187387",
        "7.72576063698292",
        "815.8202504753092",
        "48.90484332838889",
        "26.10461128235042",
        "18.19028151222189",
        "NO FINITE INTERIOR MINIMUM",
        "hohmann_only",
        "large_B_bielliptic",
        "all_bielliptic",
    ],
)
def test_report_contains_headline_numbers(expected: str) -> None:
    assert expected in REPORT.read_text(encoding="utf-8")


def test_report_records_zero_interior_minima() -> None:
    text = REPORT.read_text(encoding="utf-8")
    assert "0   (expected 0)" in text
