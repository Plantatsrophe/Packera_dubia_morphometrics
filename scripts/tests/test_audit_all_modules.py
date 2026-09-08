#!/usr/bin/env python3
"""
scripts/tests/test_audit_all_modules.py
=======================================
Audit test suite verifying importability, structural integrity, and contract
preservation across all active botanical morphometrics pipeline modules.
"""

import importlib
import sys
import unittest
from pathlib import Path

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

ACTIVE_MODULES = [
    "config",
    "main",
    "scripts.core.config",
    "scripts.core.harvester",
    "scripts.data_prep.01_voucher_harvester",
    "scripts.data_prep.audit_and_prune_vouchers",
    "scripts.annotation_and_training.annotate_with_sam2",
    "scripts.annotation_and_training.sam2_annotator_utils",
    "scripts.vision.lm2_geometry_utils",
    "scripts.vision.configure_leafmachine2",
    "scripts.pipeline.02_segment_and_extract",
    "scripts.analysis.05_cleanlab_vision_xai",
]

# Backwards compatibility alias
CORE_MODULES = ACTIVE_MODULES


class TestAuditAllModules(unittest.TestCase):
    """Audit test suite verifying all active modules import cleanly and adhere to contract."""

    def test_active_modules_importable(self):
        """Verify that every active module in the pipeline can be imported without error."""
        for mod_name in ACTIVE_MODULES:
            with self.subTest(module=mod_name):
                mod = importlib.import_module(mod_name)
                self.assertIsNotNone(mod, f"Module '{mod_name}' should import cleanly.")

    def test_pipeline_config_contract_across_modules(self):
        """Verify that entry points and modules expose the standard PipelineConfig."""
        import config
        import main
        import scripts.core.config as core_cfg
        seg_mod = importlib.import_module("scripts.pipeline.02_segment_and_extract")
        harv_mod = importlib.import_module("scripts.data_prep.01_voucher_harvester")

        self.assertTrue(hasattr(config, "load_config"))
        self.assertTrue(hasattr(core_cfg, "PipelineConfig"))
        self.assertTrue(hasattr(main, "PipelineConfig"))
        self.assertTrue(hasattr(seg_mod, "PipelineConfig"))
        self.assertTrue(hasattr(harv_mod, "PipelineConfig"))


def run_all_tests() -> unittest.TestResult:
    """Discovers and runs all active unit tests in scripts/tests, asserting 100% pass rate."""
    test_dir = Path(__file__).resolve().parent
    loader = unittest.TestLoader()
    suite = loader.discover(str(test_dir), pattern="test_*.py")
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    assert result.wasSuccessful(), f"Audit failed with {len(result.failures)} failures and {len(result.errors)} errors."
    return result


if __name__ == "__main__":
    result = run_all_tests()
    sys.exit(0 if result.wasSuccessful() else 1)
