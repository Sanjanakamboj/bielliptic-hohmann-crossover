#!/usr/bin/env python3
"""Verify the committed result artifacts are current and reproducible.

Two distinct properties are checked, because they are not the same thing:

1. **Determinism within an environment.** Every generator is run twice into
   separate directories and the outputs must be **byte-identical**. This is the
   property that actually guarantees the artifacts are a pure function of the
   code.

2. **The committed artifacts are current.** Regenerated output is compared with
   the committed files numerically, cell by cell, to a tight relative tolerance.

Byte equality against the committed files is reported but **not required**, and
deliberately so. The generators call ``math.exp``/``math.log`` and
``numpy.logspace``; those are libm/BLAS-backed and are not bit-identical across
C libraries, so full-precision ``repr`` output can differ in the last unit in the
last place between macOS and Linux. Demanding byte equality across platforms
would be an overclaim; demanding it *within* a platform, plus numerical agreement
everywhere, is the honest and still-strict guarantee.

Usage (from the repository root)::

    python scripts/m4_verify_artifacts.py
"""

from __future__ import annotations

import csv
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]

GENERATORS = (
    "m2_verification_report.py",
    "m3_trade_analysis.py",
    "m4_final_summary.py",
)

#: Relative tolerance for the committed-vs-regenerated numeric comparison. Wide
#: enough to absorb a few ulp of libm difference, far tighter than any physical
#: significance.
TOLERANCE = 1e-12


def _run_generators(target: pathlib.Path) -> None:
    """Run every generator with ``results/`` redirected to ``target``."""
    results = REPO_ROOT / "results"
    backup = target / "_backup"
    shutil.copytree(results, backup)
    try:
        for script in GENERATORS:
            subprocess.run(
                [sys.executable, str(REPO_ROOT / "scripts" / script)],
                cwd=REPO_ROOT, check=True, capture_output=True,
            )
        shutil.copytree(results, target / "generated")
    finally:
        shutil.rmtree(results)
        shutil.copytree(backup, results)
        shutil.rmtree(backup)


def _numbers(path: pathlib.Path) -> list[float]:
    """Every float in a result file, in document order."""
    values: list[float] = []
    if path.suffix == ".json":

        def walk(node: object) -> None:
            if isinstance(node, bool):
                return
            if isinstance(node, (int, float)):
                values.append(float(node))
            elif isinstance(node, dict):
                for key in sorted(node):
                    walk(node[key])
            elif isinstance(node, list):
                for item in node:
                    walk(item)

        walk(json.loads(path.read_text(encoding="utf-8")))
    elif path.suffix == ".csv":
        with path.open(encoding="utf-8") as handle:
            for row in csv.reader(handle):
                for cell in row:
                    try:
                        values.append(float(cell))
                    except ValueError:
                        continue
    else:
        for token in path.read_text(encoding="utf-8").replace(",", " ").split():
            try:
                values.append(float(token.strip("|()")))
            except ValueError:
                continue
    return values


def main() -> int:
    failures: list[str] = []
    with tempfile.TemporaryDirectory() as first_dir, tempfile.TemporaryDirectory() as second_dir:
        first, second = pathlib.Path(first_dir), pathlib.Path(second_dir)
        _run_generators(first)
        _run_generators(second)

        print("1. Determinism within this environment (must be byte-identical)")
        for path in sorted((first / "generated").iterdir()):
            if path.name.startswith("."):
                continue
            twin = second / "generated" / path.name
            identical = path.read_bytes() == twin.read_bytes()
            print(f"   [{'OK ' if identical else 'FAIL'}] {path.name}")
            if not identical:
                failures.append(f"{path.name} is not deterministic within one environment")

        print()
        print(f"2. Committed artifacts are current (numeric, rtol={TOLERANCE:g})")
        for path in sorted((first / "generated").iterdir()):
            if path.name.startswith("."):
                continue
            committed = REPO_ROOT / "results" / path.name
            if not committed.is_file():
                failures.append(f"{path.name} is not committed")
                print(f"   [FAIL] {path.name}: missing from results/")
                continue
            fresh_values, committed_values = _numbers(path), _numbers(committed)
            byte_equal = path.read_bytes() == committed.read_bytes()
            if len(fresh_values) != len(committed_values):
                failures.append(f"{path.name} has a different number of values")
                print(f"   [FAIL] {path.name}: {len(committed_values)} vs {len(fresh_values)} values")
                continue
            worst = 0.0
            for a, b in zip(committed_values, fresh_values, strict=True):
                scale = max(abs(a), abs(b), 1.0)
                worst = max(worst, abs(a - b) / scale)
            ok = worst <= TOLERANCE
            if not ok:
                failures.append(f"{path.name} differs by {worst:.3e} (stale?)")
            note = "byte-identical" if byte_equal else f"max rel diff {worst:.3e}"
            print(f"   [{'OK ' if ok else 'FAIL'}] {path.name}: {note}")

    print()
    if failures:
        print("ARTIFACT VERIFICATION FAILED:")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print("ARTIFACT VERIFICATION PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
