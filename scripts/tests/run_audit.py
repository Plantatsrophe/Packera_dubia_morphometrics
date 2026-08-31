import os
import subprocess
import sys

sys.path.insert(0, os.path.abspath("."))
from scripts.tests.test_audit_all_modules import CORE_MODULES

def main():
    log_file = "/tmp/lm2_audit_results.txt"
    with open(log_file, "w") as out:
        def log(msg):
            print(msg, flush=True)
            out.write(msg + "\n")
            out.flush()

        log(f"Testing {len(CORE_MODULES)} modules in isolated subprocesses...\n")
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
                log(f"  [PASS] {mod}")
            else:
                err = res.stderr.strip().splitlines()[-1] if res.stderr.strip() else "Non-zero exit"
                failed.append((mod, err))
                log(f"  [FAIL] {mod} -> {err}")

        log(f"\n==========================================")
        log(f"SUMMARY: {len(passed)}/{len(CORE_MODULES)} PASSED, {len(failed)} FAILED")
        log(f"==========================================")
        if failed:
            sys.exit(1)

if __name__ == "__main__":
    main()
