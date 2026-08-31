import os
import subprocess
import sys

sys.path.insert(0, os.path.abspath("."))
from scripts.tests.test_audit_all_modules import CORE_MODULES

def main():
    print(f"Testing {len(CORE_MODULES)} modules in isolated subprocesses...\n")
    passed = []
    failed = []

    for mod in CORE_MODULES:
        cmd = [
            sys.executable, "-c",
            f'import sys, os; sys.path.insert(0, "."); sys.path.insert(0, "LeafMachine2"); sys.path.insert(0, "LeafMachine2/leafmachine2"); import {mod}'
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0:
            passed.append(mod)
            print(f"  [PASS] {mod}")
        else:
            err = res.stderr.strip().splitlines()[-1] if res.stderr.strip() else "Non-zero exit"
            failed.append((mod, err))
            print(f"  [FAIL] {mod} -> {err}")

    print(f"\n==========================================")
    print(f"SUMMARY: {len(passed)}/{len(CORE_MODULES)} PASSED, {len(failed)} FAILED")
    print(f"==========================================")
    if failed:
        sys.exit(1)

if __name__ == "__main__":
    main()
