#!/usr/bin/env python3
"""
scripts/tests/run_audit.py
==========================
Centralized test runner and QA audit script for Packera dubia morphometrics pipeline.
Executes all active unit and integration tests via unittest (or fast Python-only
checks when --fast / --python-only is supplied), logs execution results, tears down
any warm Rscript worker sessions cleanly, and asserts a 100% pass rate.
"""

from __future__ import annotations

import argparse
import os
import sys
import unittest
from pathlib import Path

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests._r_session import close_all_warm_r_sessions


def run_audit(
    verbose: bool = True,
    log_file: str = "/tmp/lm3_audit_results.txt",
    fast: bool = False,
    pattern: str = "test_*.py",
) -> unittest.TestResult:
    """Discovers and executes active unit tests in scripts/tests, asserting 100% pass rate."""
    prev_fast = os.environ.get("PACKERA_FAST_AUDIT")
    if fast:
        os.environ["PACKERA_FAST_AUDIT"] = "1"

    try:
        test_dir = PROJECT_ROOT / "scripts" / "tests"
        loader = unittest.TestLoader()
        suite = loader.discover(str(test_dir), pattern=pattern)
        mode_tag = " [FAST PYTHON-ONLY MODE]" if fast else ""

        header = (
            f"==================================================================\n"
            f"PACKERA DUBIA MORPHOMETRICS - CENTRAL TEST SUITE & QA AUDIT RUNNER{mode_tag}\n"
            f"Discovered {suite.countTestCases()} unit and integration tests\n"
            f"=================================================================="
        )
        print(header, flush=True)

        runner = unittest.TextTestRunner(verbosity=2 if verbose else 1)
        result = runner.run(suite)
    finally:
        close_all_warm_r_sessions()
        if fast:
            if prev_fast is None:
                os.environ.pop("PACKERA_FAST_AUDIT", None)
            else:
                os.environ["PACKERA_FAST_AUDIT"] = prev_fast

    total_run = result.testsRun
    num_failures = len(result.failures)
    num_errors = len(result.errors)
    num_skipped = len(result.skipped)
    passed = total_run - num_failures - num_errors
    pass_rate = (passed / total_run * 100.0) if total_run > 0 else 0.0

    summary = (
        f"\n==================================================================\n"
        f"AUDIT SUMMARY:\n"
        f"  Total tests executed: {total_run}\n"
        f"  Passed:               {passed}\n"
        f"  Failures:             {num_failures}\n"
        f"  Errors:               {num_errors}\n"
        f"  Skipped:              {num_skipped}\n"
        f"  Pass Rate:            {pass_rate:.1f}%\n"
        f"=================================================================="
    )
    print(summary, flush=True)

    try:
        with open(log_file, "w", encoding="utf-8") as f:
            f.write(header + "\n\n")
            f.write(summary + "\n")
            if result.failures:
                f.write("\nFAILURES:\n")
                for test, err in result.failures:
                    f.write(f"- {test}: {err}\n")
            if result.errors:
                f.write("\nERRORS:\n")
                for test, err in result.errors:
                    f.write(f"- {test}: {err}\n")
    except Exception as e:
        print(f"Warning: Could not write to log file {log_file}: {e}", file=sys.stderr)

    assert result.wasSuccessful(), (
        f"Central QA Audit Failed: {num_failures} failures and {num_errors} errors detected (Pass Rate: {pass_rate:.1f}%)."
    )

    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Centralized QA test runner for the Packera morphometrics pipeline."
    )
    parser.add_argument(
        "--fast",
        "--python-only",
        dest="fast",
        action="store_true",
        default=False,
        help="Skip Rscript statistical/GIS integration suites for fast (~20s) Python-only verification.",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        default=False,
        help="Use compact single-dot unittest output instead of verbose per-test lines.",
    )
    parser.add_argument(
        "--pattern",
        default="test_*.py",
        help="Test filename discovery glob pattern (default: test_*.py).",
    )
    parser.add_argument(
        "--log-file",
        default="/tmp/lm3_audit_results.txt",
        help="Path to write QA audit summary log.",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    try:
        res = run_audit(
            verbose=not args.quiet,
            log_file=args.log_file,
            fast=args.fast,
            pattern=args.pattern,
        )
        sys.exit(0 if res.wasSuccessful() else 1)
    except AssertionError as err:
        print(f"\n[FATAL AUDIT FAILURE] {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
