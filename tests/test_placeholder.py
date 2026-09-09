"""Milestone 1 scaffold tests.

These tests deliberately assert only what legitimately exists at M1: that the
package imports, that it reports the expected version, that the repository
scaffold is in place, and that no Milestone 2 production API has been
introduced ahead of schedule.

There are no numerical physics tests here, because no numerical physics code
has been implemented yet. Adding assertions about transfer delta-v now would
be testing constants transcribed from a document rather than testing code.
"""

from __future__ import annotations

import pathlib

import pytest

import bielliptic_crossover

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]


def test_package_imports() -> None:
    assert bielliptic_crossover is not None


def test_version_string() -> None:
    assert bielliptic_crossover.__version__ == "0.1.0"


@pytest.mark.parametrize(
    "relative_path",
    [
        "README.md",
        "DESIGN.md",
        ".gitignore",
        "pyproject.toml",
        "src/bielliptic_crossover/__init__.py",
        "tests/test_placeholder.py",
        "scripts/.gitkeep",
        "figures/.gitkeep",
        "results/.gitkeep",
    ],
)
def test_scaffold_file_exists(relative_path: str) -> None:
    assert (REPO_ROOT / relative_path).is_file()


def test_src_layout_is_used() -> None:
    """The package must be importable from src/, not from the repo root."""
    assert (REPO_ROOT / "src" / "bielliptic_crossover").is_dir()
    assert not (REPO_ROOT / "bielliptic_crossover").exists()


@pytest.mark.parametrize(
    "attribute_name",
    [
        # Milestone 2+ transfer solver / sweep / optimizer / plotting surface.
        "hohmann_delta_v",
        "bielliptic_delta_v",
        "delta_v_hohmann",
        "delta_v_bielliptic",
        "crossover_ratio",
        "find_crossover",
        "optimize_b",
        "optimal_b",
        "sweep",
        "sweep_ratio",
        "transfer_time",
        "plot_crossover",
    ],
)
def test_no_milestone_2_api_yet(attribute_name: str) -> None:
    """M1 is derivation only; no production solver API should exist yet."""
    assert not hasattr(bielliptic_crossover, attribute_name)


def test_public_api_is_version_only() -> None:
    assert bielliptic_crossover.__all__ == ["__version__"]


def test_no_production_modules_yet() -> None:
    package_dir = REPO_ROOT / "src" / "bielliptic_crossover"
    modules = sorted(p.name for p in package_dir.glob("*.py"))
    assert modules == ["__init__.py"]
