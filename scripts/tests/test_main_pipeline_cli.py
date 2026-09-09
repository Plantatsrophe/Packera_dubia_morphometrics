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
from config import DEFAULT_CONFIG_PATH, load_config
from scripts.core.artifact_manager import (
    compute_sha256,
    download_file_with_progress,
    ensure_model_weights,
)


class TestMainPipelineCLI(unittest.TestCase):
    """Test suite for main.py entrypoint and CLI argument parsing."""

    def test_load_config_structure(self):
        """Verify config.yaml loads with all required pipeline sections."""
        cfg = load_config(DEFAULT_CONFIG_PATH)
        self.assertIn("paths", cfg)
        self.assertIn("models", cfg)
        self.assertIn("taxonomy", cfg)
        self.assertIn("harvesting", cfg)
        self.assertIn("segmentation", cfg)
        self.assertIn("morphometrics", cfg)
        self.assertIn("pcd_weights_path", cfg["models"])
        self.assertIn("pcd_weights_url", cfg["models"])

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

        # check-env
        args_check = parser.parse_args(["check-env", "--strict"])
        self.assertEqual(args_check.subcommand, "check-env")
        self.assertTrue(args_check.strict)

        # download-weights
        args_dl = parser.parse_args(["download-weights", "--force"])
        self.assertEqual(args_dl.subcommand, "download-weights")
        self.assertTrue(args_dl.force)

        # calibrate-geometry
        args_calib = parser.parse_args([
            "calibrate-geometry",
            "--annotations", "data/annotations/packera_train_coco.json",
            "--update-config",
            "--min-area", "200.0"
        ])
        self.assertEqual(args_calib.subcommand, "calibrate-geometry")
        self.assertEqual(str(args_calib.annotations), "data/annotations/packera_train_coco.json")
        self.assertTrue(args_calib.update_config)
        self.assertEqual(args_calib.min_area, 200.0)

        # force / overwrite flags across subparsers
        self.assertTrue(parser.parse_args(["harvest", "--force"]).force)
        self.assertTrue(parser.parse_args(["harvest", "--overwrite"]).force)
        self.assertTrue(parser.parse_args(["segment", "--force"]).force)
        self.assertTrue(parser.parse_args(["segment", "--overwrite"]).force)
        self.assertTrue(parser.parse_args(["download-weights", "--force"]).force)
        self.assertTrue(parser.parse_args(["download-weights", "--overwrite"]).force)
        self.assertTrue(parser.parse_args(["run-all", "--force"]).force)
        self.assertTrue(parser.parse_args(["run-all", "--overwrite"]).force)

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

    def test_run_harvest_forwards_force_flag(self):
        """Verify run_harvest appends --force flag to child subprocess command."""
        from unittest.mock import MagicMock
        parser = main.build_parser()
        args = parser.parse_args(["harvest", "--force"])
        cfg = {
            "paths": {
                "raw_vouchers_dir": "/tmp/dummy_raw",
                "curated_vouchers_csv": "/tmp/dummy_curated.csv",
            },
            "harvesting": {
                "max_records_per_taxon": 100,
                "max_uncertainty_meters": 5000,
                "min_megapixels": 8.0,
                "min_file_size_kb": 500,
                "download_concurrency": 10,
            },
        }
        with patch("subprocess.run") as mock_subp, \
             patch("main.verify_file_exists"):
            main.run_harvest(args, cfg)
            mock_subp.assert_called_once()
            cmd_args = mock_subp.call_args[0][0]
            self.assertIn("--force", cmd_args)

    def test_run_segment_forwards_force_flag(self):
        """Verify run_segment appends --force flag to child subprocess command."""
        parser = main.build_parser()
        args = parser.parse_args(["segment", "--force"])
        cfg = {
            "paths": {
                "curated_vouchers_csv": "/tmp/dummy_curated.csv",
                "model_weights": "models/dummy.pth",
                "workspace_root": "/tmp/dummy_root",
                "contours_dir": "/tmp/dummy_contours",
            },
            "segmentation": {
                "device": "cpu",
                "min_solidity": 0.72,
                "min_ucs": 0.85,
                "score_thresh": 0.40,
            },
        }
        with patch("subprocess.run") as mock_subp, \
             patch("main.verify_file_exists"), \
             patch("main.verify_dir_has_files"):
            main.run_segment(args, cfg)
            mock_subp.assert_called_once()
            cmd_args = mock_subp.call_args[0][0]
            self.assertIn("--force", cmd_args)

    def test_main_dispatch_run_all(self):
        """Verify main() entrypoint dispatches run-all subcommand to handler."""
        with patch.object(sys, "argv", ["main.py", "run-all"]), \
             patch("main.run_all") as mock_run_all:
            main.main()
            mock_run_all.assert_called_once()

    def test_check_environment_clean_pass(self):
        """Verify check_environment() returns True and passes cleanly when dependencies are satisfied."""
        dummy_completed_probe = MagicMock(returncode=0, stdout="ALL_INSTALLED\n", stderr="")
        dummy_ver_probe = MagicMock(returncode=0, stdout="4.3.2\n", stderr="")

        def mock_subp_run(cmd, *args, **kwargs):
            if isinstance(cmd, list) and len(cmd) > 2 and "getRversion" in str(cmd[2]):
                return dummy_ver_probe
            return dummy_completed_probe

        with patch("shutil.which", return_value="/usr/bin/Rscript"), \
             patch("subprocess.run", side_effect=mock_subp_run), \
             patch("importlib.import_module", return_value=MagicMock(__version__="1.0.0")):
            result = main.check_environment()
            self.assertTrue(result)

    def test_check_env_cli_exits_cleanly(self):
        """Verify 'main.py check-env' CLI subcommand exits with status code 0 when diagnostics pass."""
        with patch.object(sys, "argv", ["main.py", "check-env"]), \
             patch("main.check_environment", return_value=True):
            with self.assertRaises(SystemExit) as cm:
                main.main()
            self.assertEqual(cm.exception.code, 0)

    def test_check_env_cli_exits_code_1_on_failure(self):
        """Verify 'main.py check-env' CLI subcommand exits with status code 1 when diagnostics fail."""
        with patch.object(sys, "argv", ["main.py", "check-env"]), \
             patch("main.check_environment", return_value=False):
            with self.assertRaises(SystemExit) as cm:
                main.main()
            self.assertEqual(cm.exception.code, 1)

    def test_check_environment_missing_rscript(self):
        """Verify check_environment() reports failure when Rscript is missing from PATH."""
        with patch("shutil.which", return_value=None):
            result = main.check_environment()
            self.assertFalse(result)

    def test_check_environment_missing_r_packages_remediation(self):
        """Verify check_environment() identifies missing R packages and prints remediation command."""
        missing_probe = MagicMock(
            returncode=1,
            stdout="MISSING: Momocs,MorphoTools2\n",
            stderr="",
        )
        with patch("shutil.which", return_value="/usr/bin/Rscript"), \
             patch("subprocess.run", return_value=missing_probe):
            result = main.check_environment()
            self.assertFalse(result)


    def test_get_lm2_python_executable_unix(self):
        """Verify get_lm2_python_executable resolves .venv_LM2/bin/python on Unix systems."""
        with patch.object(Path, "exists", side_effect=[True]):
            resolved = main.get_lm2_python_executable()
            self.assertEqual(resolved, main.PROJECT_ROOT / ".venv_LM2" / "bin" / "python")

    def test_get_lm2_python_executable_windows(self):
        """Verify get_lm2_python_executable resolves .venv_LM2/Scripts/python.exe on Windows."""
        with patch.object(Path, "exists", side_effect=[False, True]):
            resolved = main.get_lm2_python_executable()
            self.assertEqual(resolved, main.PROJECT_ROOT / ".venv_LM2" / "Scripts" / "python.exe")

    def test_get_lm2_python_executable_fallback_warning(self):
        """Verify get_lm2_python_executable falls back to sys.executable and warns if .venv_LM2 is absent."""
        with patch.object(Path, "exists", return_value=False), \
             patch.object(main.logger, "warning") as mock_warn:
            resolved = main.get_lm2_python_executable()
            self.assertEqual(resolved, Path(sys.executable))
            mock_warn.assert_called_once()
            self.assertIn("LeafMachine2 dedicated virtual environment", mock_warn.call_args[0][0])

    def test_run_segment_uses_resolved_lm2_interpreter(self):
        """Verify run_segment invokes 02_segment_and_extract.py with resolved LM2 interpreter and logs diagnostic."""
        parser = main.build_parser()
        args = parser.parse_args(["segment"])
        cfg = {
            "paths": {
                "curated_vouchers_csv": "/tmp/dummy_curated.csv",
                "model_weights": "models/dummy.pth",
                "contours_dir": "/tmp/dummy_contours",
            },
            "segmentation": {
                "device": "cpu",
                "min_solidity": 0.72,
                "min_ucs": 0.85,
                "score_thresh": 0.40,
            },
        }
        dummy_lm2_python = main.PROJECT_ROOT / ".venv_LM2" / "bin" / "python"
        with patch("main.get_lm2_python_executable", return_value=dummy_lm2_python), \
             patch("subprocess.run") as mock_subp, \
             patch("main.verify_file_exists"), \
             patch("main.verify_dir_has_files"), \
             patch.object(main.logger, "info") as mock_info:
            main.run_segment(args, cfg)
            mock_subp.assert_called_once()
            cmd_called = mock_subp.call_args[0][0]
            self.assertEqual(cmd_called[0], str(dummy_lm2_python))
            # Verify diagnostic log entry
            log_messages = [call[0][0] for call in mock_info.call_args_list if call[0]]
            self.assertTrue(
                any("Executing LeafMachine2 segmentation via: .venv_LM2/bin/python" in msg for msg in log_messages),
                f"Expected diagnostic log missing from: {log_messages}",
            )

    def test_run_segment_propagates_non_zero_exit_code(self):
        """Verify run_segment catches CalledProcessError and immediately exits with non-zero status."""
        parser = main.build_parser()
        args = parser.parse_args(["segment"])
        cfg = {
            "paths": {
                "curated_vouchers_csv": "/tmp/dummy_curated.csv",
                "model_weights": "models/dummy.pth",
                "contours_dir": "/tmp/dummy_contours",
            },
            "segmentation": {
                "device": "cpu",
                "min_solidity": 0.72,
                "min_ucs": 0.85,
                "score_thresh": 0.40,
            },
        }
        import subprocess as sp
        error = sp.CalledProcessError(returncode=2, cmd=["dummy_cmd"])
        with patch("subprocess.run", side_effect=error), \
             patch("main.verify_file_exists"):
            with self.assertRaises(SystemExit) as cm:
                main.run_segment(args, cfg)
            self.assertEqual(cm.exception.code, 2)

    def test_ensure_model_weights_preexisting(self):
        """Verify ensure_model_weights returns existing non-empty weights without downloading."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            dummy_pth = Path(tmp_dir) / "test_weights.pth"
            dummy_pth.write_bytes(b"dummy model checkpoint data")

            cfg = {
                "models": {
                    "pcd_weights_path": str(dummy_pth),
                    "pcd_weights_url": "https://example.com/weights.pth",
                    "pcd_weights_sha256": "",
                }
            }
            with patch("scripts.core.artifact_manager.download_file_with_progress") as mock_dl:
                result = ensure_model_weights(cfg, force=False)
                self.assertEqual(result, dummy_pth)
                mock_dl.assert_not_called()

    def test_ensure_model_weights_missing_url_graceful_error(self):
        """Verify ensure_model_weights raises ValueError when weights are missing and URL is empty."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            non_existent_pth = Path(tmp_dir) / "missing_weights.pth"
            cfg = {
                "models": {
                    "pcd_weights_path": str(non_existent_pth),
                    "pcd_weights_url": "",
                    "pcd_weights_sha256": "",
                }
            }
            with self.assertRaises(ValueError) as ctx:
                ensure_model_weights(cfg, force=False)
            self.assertIn("no remote download URL", str(ctx.exception))

    def test_ensure_model_weights_download_success(self):
        """Verify ensure_model_weights triggers download and atomic staging when file is missing."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            target_pth = Path(tmp_dir) / "downloaded_weights.pth"
            cfg = {
                "models": {
                    "pcd_weights_path": str(target_pth),
                    "pcd_weights_url": "https://github.com/Plantatsrophe/Packera_dubia_morphometrics/releases/download/v1.0-weights/lm2_packera_pcd_finetuned.pth",
                    "pcd_weights_sha256": "",
                }
            }

            def fake_download(url, dest, expected_sha256=None, chunk_size=None):
                dest = Path(dest)
                dest.write_bytes(b"downloaded simulated checkpoint")
                return dest

            with patch("scripts.core.artifact_manager.download_file_with_progress", side_effect=fake_download) as mock_dl:
                res = ensure_model_weights(cfg, force=False)
                self.assertEqual(res, target_pth)
                mock_dl.assert_called_once()
                self.assertTrue(target_pth.exists())

    def test_ensure_model_weights_checksum_verification(self):
        """Verify ensure_model_weights and download_file_with_progress enforce SHA256 integrity."""
        import hashlib
        data = b"synthetic torch model payload"
        correct_hash = hashlib.sha256(data).hexdigest()
        bad_hash = "0123456789abcdef" * 4

        with tempfile.TemporaryDirectory() as tmp_dir:
            target_pth = Path(tmp_dir) / "payload.pth"

            class DummyResponse:
                def __init__(self):
                    self.headers = {"Content-Length": str(len(data))}
                def read(self, chunk_size):
                    if not hasattr(self, "_read_done"):
                        self._read_done = True
                        return data
                    return b""
                def __enter__(self):
                    return self
                def __exit__(self, *args):
                    pass

            # Test successful checksum
            with patch("urllib.request.urlopen", return_value=DummyResponse()):
                download_file_with_progress("https://example.com/model.pth", target_pth, expected_sha256=correct_hash)
                self.assertTrue(target_pth.exists())
                self.assertEqual(compute_sha256(target_pth), correct_hash)

            # Test checksum mismatch raises ValueError and removes partial file
            with patch("urllib.request.urlopen", return_value=DummyResponse()):
                with self.assertRaises(ValueError) as ctx:
                    download_file_with_progress("https://example.com/model.pth", target_pth, expected_sha256=bad_hash)
                self.assertIn("Checksum verification failed", str(ctx.exception))

    def test_main_dispatch_download_weights(self):
        """Verify main() entrypoint dispatches download-weights subcommand."""
        with patch.object(sys, "argv", ["main.py", "download-weights", "--force"]), \
             patch("main.run_download_weights") as mock_dl:
            main.main()
            mock_dl.assert_called_once()
            args, cfg = mock_dl.call_args[0]
            self.assertEqual(args.subcommand, "download-weights")
            self.assertTrue(args.force)


if __name__ == "__main__":
    unittest.main()

