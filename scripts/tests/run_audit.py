#!/usr/bin/env python3
"""
scripts/tests/run_audit.py
==========================
Centralized test runner and QA audit script for Packera dubia morphometrics pipeline.
Executes all active unit and integration tests via unittest, logs execution results,
and asserts a 100% pass rate.
"""

import os
import sys
import unittest
from pathlib import Path

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def run_audit(verbose: bool = True, log_file: str = "/tmp/lm2_audit_results.txt") -> unittest.TestResult:
    """Discovers and executes all active unit tests in scripts/tests, asserting 100% pass rate."""
    test_dir = PROJECT_ROOT / "scripts" / "tests"
    loader = unittest.TestLoader()
    suite = loader.discover(str(test_dir), pattern="test_*.py")

    header = (
        f"==================================================================\n"
        f"PACKERA DUBIA MORPHOMETRICS - CENTRAL TEST SUITE & QA AUDIT RUNNER\n"
        f"Discovered {suite.countTestCases()} unit and integration tests\n"
        f"=================================================================="
    )
    print(header, flush=True)

    runner = unittest.TextTestRunner(verbosity=2 if verbose else 1)
    result = runner.run(suite)

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

    # Write results to audit log
    try:
        with open(log_file, "w") as f:
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


def main():
    try:
        res = run_audit(verbose=True)
        sys.exit(0 if res.wasSuccessful() else 1)
    except AssertionError as err:
        print(f"\n[FATAL AUDIT FAILURE] {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
