#!/usr/bin/env python3
"""
Unit tests for centralized configuration and root CLI runner (main.py).
Tests argument parsing, preflight checks, dispatch, and dry-run execution.
"""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import main
from config import load_config, DEFAULT_CONFIG_PATH


class TestMainPipelineCLI(unittest.TestCase):
    """Test suite for main.py entrypoint and CLI argument parsing."""

    def test_load_config_structure(self):
        """Verify config.yaml loads with all required pipeline sections."""
        cfg = load_config(DEFAULT_CONFIG_PATH)
        self.assertIn("paths", cfg)
        self.assertIn("taxonomy", cfg)
        self.assertIn("harvesting", cfg)
        self.assertIn("segmentation", cfg)
        self.assertIn("morphometrics", cfg)

        # Verify resolved paths
        self.assertIn("resolved_paths", cfg)
        self.assertTrue(cfg["resolved_paths"]["raw_vouchers_dir"].is_absolute())
        self.assertTrue(cfg["resolved_paths"]["curated_vouchers_csv"].is_absolute())

        # Verify key hyperparameters
        self.assertEqual(cfg["harvesting"]["max_records_per_taxon"], 5000)
        self.assertEqual(cfg["segmentation"]["score_thresh"], 0.40)
        self.assertEqual(cfg["morphometrics"]["harmonics"], 12)

    def test_build_parser_subcommands(self):
        """Verify CLI parser parses all subcommands correctly."""
        parser = main.build_parser()

        # harvest
        args_harvest = parser.parse_args(["harvest", "--max-records", "100", "--download-images"])
        self.assertEqual(args_harvest.subcommand, "harvest")
        self.assertEqual(args_harvest.max_records, 100)
        self.assertTrue(args_harvest.download_images)

        # segment
        args_seg = parser.parse_args(["segment", "--device", "cpu", "--limit", "10"])
        self.assertEqual(args_seg.subcommand, "segment")
        self.assertEqual(args_seg.device, "cpu")
        self.assertEqual(args_seg.limit, 10)

        # morphometrics
        args_morph = parser.parse_args(["morphometrics", "--harmonics", "15", "--num-pcs", "6"])
        self.assertEqual(args_morph.subcommand, "morphometrics")
        self.assertEqual(args_morph.harmonics, 15)
        self.assertEqual(args_morph.num_pcs, 6)

        # synthesis
        args_synth = parser.parse_args([
            "synthesis",
            "--cleanlab-threshold", "0.90",
            "--no-export-figures",
            "--permutations", "50",
        ])
        self.assertEqual(args_synth.subcommand, "synthesis")
        self.assertEqual(args_synth.cleanlab_threshold, 0.90)
        self.assertFalse(args_synth.export_figures)
        self.assertEqual(args_synth.permutations, 50)

        # run-all
        args_all = parser.parse_args(["run-all", "--download-images"])
        self.assertEqual(args_all.subcommand, "run-all")
        self.assertTrue(args_all.download_images)

    def test_gpu_availability_check(self):
        """Verify GPU availability check runs and returns boolean without unhandled exception."""
        result = main.check_gpu_availability(warn_only=True)
        self.assertIsInstance(result, bool)

    def test_verify_file_exists_fails_cleanly(self):
        """Verify preflight check raises SystemExit on missing upstream file."""
        with self.assertRaises(SystemExit):
            main.verify_file_exists(
                Path("/non/existent/path/voucher_table.csv"),
                description="test table",
                hint_cmd="python main.py harvest",
            )

    def test_verify_dir_has_files_fails_cleanly(self):
        """Verify preflight check raises SystemExit on empty directory."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            empty_dir = Path(tmp_dir) / "empty_dir"
            empty_dir.mkdir()
            with self.assertRaises(SystemExit):
                main.verify_dir_has_files(
                    empty_dir,
                    pattern="*.csv",
                    description="test empty dir",
                    hint_cmd="python main.py segment",
                )

    def test_main_dispatch_harvest(self):
        """Verify main() entrypoint dispatches harvest subcommand to handler."""
        with patch.object(sys, "argv", ["main.py", "harvest", "--max-records", "10"]), \
             patch("main.run_harvest") as mock_harvest:
            main.main()
            mock_harvest.assert_called_once()
            args, cfg = mock_harvest.call_args[0]
            self.assertEqual(args.subcommand, "harvest")
            self.assertEqual(args.max_records, 10)

    def test_main_dispatch_segment(self):
        """Verify main() entrypoint dispatches segment subcommand to handler."""
        with patch.object(sys, "argv", ["main.py", "segment", "--device", "cpu"]), \
             patch("main.run_segment") as mock_segment:
            main.main()
            mock_segment.assert_called_once()
            args, cfg = mock_segment.call_args[0]
            self.assertEqual(args.subcommand, "segment")
            self.assertEqual(args.device, "cpu")

    def test_main_dispatch_run_all(self):
        """Verify main() entrypoint dispatches run-all subcommand to handler."""
        with patch.object(sys, "argv", ["main.py", "run-all"]), \
             patch("main.run_all") as mock_run_all:
            main.main()
            mock_run_all.assert_called_once()


if __name__ == "__main__":
    unittest.main()
