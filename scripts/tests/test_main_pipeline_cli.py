#!/usr/bin/env python3
"""
Unit tests for centralized configuration and root CLI runner (main.py).
Tests argument parsing, preflight checks, dispatch, and dry-run execution.
"""

import argparse
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
        self.assertEqual(cfg["segmentation"]["score_thresh"], 0.35)
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
            "--no-export-figures",
            "--permutations", "50",
        ])
        self.assertEqual(args_synth.subcommand, "synthesis")
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

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            dummy_weights = tmp_path / "lm3_leaf_segmenter__yolo26x_seg_1024.onnx"
            dummy_weights.write_bytes(b"MOCK_WEIGHTS_BYTES" * 64)
            dummy_cfg = tmp_path / "config.yaml"
            dummy_cfg.write_text(
                f"paths:\n  model_weights: {dummy_weights}\n",
                encoding="utf-8",
            )
            args = argparse.Namespace(strict=False, config=str(dummy_cfg))

            with patch("shutil.which", return_value="/usr/bin/Rscript"), \
                 patch("subprocess.run", side_effect=mock_subp_run), \
                 patch("importlib.import_module", return_value=MagicMock(__version__="1.0.0")), \
                 patch("main.ensure_model_weights", return_value=dummy_weights):
                result = main.check_environment(args)
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
        with patch("shutil.which", return_value=None), \
             patch("main.ensure_model_weights", return_value=Path("/tmp/mock.pth")):
            result = main.check_environment()
            self.assertFalse(result)

    def test_check_environment_missing_r_packages_remediation(self):
        """Verify check_environment() identifies missing R packages and prints remediation command."""
        try:
            import pygbif  # noqa: F401
        except Exception:
            pass
        missing_probe = MagicMock(
            returncode=1,
            stdout="MISSING: Momocs,MorphoTools2\n",
            stderr="",
        )
        with patch("shutil.which", return_value="/usr/bin/Rscript"), \
             patch("subprocess.run", return_value=missing_probe), \
             patch("main.ensure_model_weights", return_value=Path("/tmp/mock.pth")):
            result = main.check_environment()
            self.assertFalse(result)
        main_src = (PROJECT_ROOT / "main.py").read_text(encoding="utf-8")
        self.assertIn("remotes::install_github", main_src)
        self.assertIn("MomX/Momocs", main_src)


    def test_get_lm3_python_executable_unix(self):
        """Verify get_lm3_python_executable resolves LeafMachine3/.venv/bin/python and .venv_LM3/bin/python on Unix systems."""
        with patch.object(Path, "exists", side_effect=[True]):
            resolved = main.get_lm3_python_executable()
            self.assertEqual(resolved, main.PROJECT_ROOT / "LeafMachine3" / ".venv" / "bin" / "python")
        with patch.object(Path, "exists", side_effect=[False, True]):
            resolved_alt = main.get_lm3_python_executable()
            self.assertEqual(resolved_alt, main.PROJECT_ROOT / ".venv_LM3" / "bin" / "python")

    def test_get_lm3_python_executable_windows(self):
        """Verify get_lm3_python_executable resolves LeafMachine3/.venv/Scripts/python.exe and .venv_LM3/Scripts/python.exe on Windows."""
        with patch.object(Path, "exists", side_effect=[False, False, True]):
            resolved = main.get_lm3_python_executable()
            self.assertEqual(resolved, main.PROJECT_ROOT / "LeafMachine3" / ".venv" / "Scripts" / "python.exe")
        with patch.object(Path, "exists", side_effect=[False, False, False, True]):
            resolved_alt = main.get_lm3_python_executable()
            self.assertEqual(resolved_alt, main.PROJECT_ROOT / ".venv_LM3" / "Scripts" / "python.exe")

    def test_get_lm3_python_executable_fallback_warning(self):
        """Verify get_lm3_python_executable falls back to sys.executable and warns if LeafMachine3/.venv and .venv_LM3 are absent."""
        with patch.object(Path, "exists", return_value=False), \
             patch.object(main.logger, "warning") as mock_warn:
            resolved = main.get_lm3_python_executable()
            self.assertEqual(resolved, Path(sys.executable))
            mock_warn.assert_called_once()
            self.assertIn("LeafMachine3 dedicated virtual environment", mock_warn.call_args[0][0])

    def test_run_segment_uses_resolved_lm3_interpreter(self):
        """Verify run_segment invokes 02_segment_and_extract.py with resolved LM3 interpreter and logs diagnostic."""
        parser = main.build_parser()
        args = parser.parse_args(["segment"])
        cfg = {
            "paths": {
                "curated_vouchers_csv": "/tmp/dummy_curated.csv",
                "model_weights": "models/dummy.onnx",
                "contours_dir": "/tmp/dummy_contours",
            },
            "segmentation": {
                "device": "cpu",
                "min_solidity": 0.72,
                "min_ucs": 0.85,
                "score_thresh": 0.40,
            },
        }
        dummy_lm3_python = main.PROJECT_ROOT / "LeafMachine3" / ".venv" / "bin" / "python"
        with patch("main.get_lm3_python_executable", return_value=dummy_lm3_python), \
             patch("subprocess.run") as mock_subp, \
             patch("main.verify_file_exists"), \
             patch("main.verify_dir_has_files"), \
             patch.object(main.logger, "info") as mock_info:
            main.run_segment(args, cfg)
            mock_subp.assert_called_once()
            cmd_called = mock_subp.call_args[0][0]
            self.assertEqual(cmd_called[0], str(dummy_lm3_python))
            # Verify diagnostic log entry
            log_messages = [call[0][0] for call in mock_info.call_args_list if call[0]]
            self.assertTrue(
                any("Executing LeafMachine3 segmentation via: LeafMachine3/.venv/bin/python" in msg for msg in log_messages),
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
                    "pcd_weights_url": "https://github.com/Plantatsrophe/Packera_dubia_morphometrics/releases/download/v2.0-weights/lm3_leaf_segmenter__yolo26x_seg_1024.onnx",
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

    def test_run_pipeline_sh_gold_only_and_determiner_tier_filtering(self):
        """Verify run_pipeline.sh --gold-only and --determiner-tier filter curated_vouchers.csv, Step 02 manifests, and cropped_patches strictly to Tier_1_Gold before Steps 02-07."""
        import subprocess
        import pandas as pd

        sh_path = PROJECT_ROOT / "run_pipeline.sh"
        sh_src = sh_path.read_text(encoding="utf-8")
        self.assertIn("--gold-only)", sh_src)
        self.assertIn("--determiner-tier)", sh_src)
        self.assertIn("--determiner-tier=*)", sh_src)
        # Verify preflight filter runs before Step 02 segmentation
        preflight_pos = sh_src.index("\nrun_preflight_subset_filter\n")
        step02_pos = sh_src.index("02_segment_and_extract.py")
        self.assertLess(preflight_pos, step02_pos)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            root_data = tmp_path / "root_data"
            run_dir = tmp_path / "runs" / "test_gold_only"
            tables_dir = run_dir / "data" / "tables"
            tables_dir.mkdir(parents=True, exist_ok=True)

            for sub in ("capitula",):
                pdir = root_data / "cropped_patches" / sub
                pdir.mkdir(parents=True, exist_ok=True)

            vouchers_df = pd.DataFrame({
                "catalogNumber": ["NCU001", "UNCC_49308", "US002", "NY003"],
                "scientificName": ["Packera dubia", "Packera anonyma", "Packera paupercula", "Packera plattensis"],
                "determiner_tier": ["Tier_1_Gold", "Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze"],
                "authority_weight": [3.0, 3.0, 1.5, 0.5],
            })
            vouchers_df.to_csv(tables_dir / "curated_vouchers.csv", index=False)

            for tbl_name in (
                "extracted_leaf_manifest.csv",
                "extracted_leaves_manifest.csv",
                "failed_qc_vouchers.csv",
                "voucher_reproductive_metrics.csv",
                "voucher_phenological_states.csv",
            ):
                pd.DataFrame({
                    "catalogNumber": ["NCU001", "UNCC_49308", "US002", "NY003"],
                    "plant_individual_id": [0, 0, 0, 0],
                }).to_csv(tables_dir / tbl_name, index=False)

            # Create cropped capitula patches for Gold and Silver vouchers
            (root_data / "cropped_patches" / "capitula" / "NCU001_cap1.jpg").write_bytes(b"jpg")
            (root_data / "cropped_patches" / "capitula" / "US002_cap1.jpg").write_bytes(b"jpg")

            # Extract the run_preflight_subset_filter function from run_pipeline.sh and run it
            fn_start = sh_src.index("run_preflight_subset_filter() {")
            fn_end = sh_src.index("\nrun_preflight_subset_filter\n", fn_start)
            fn_body = sh_src[fn_start:fn_end]

            test_bash = f"""
set -euo pipefail
PROJECT_ROOT="{PROJECT_ROOT}"
RUN_DIR="{run_dir}"
PYTHON_MAIN="{sys.executable}"
FOCAL_TAXON=""
EXCLUDE_TAXON_CSV=""
INTRASPECIFIC_MODE="false"
GOLD_ANCHORS=30
GOLD_ONLY=true
DETERMINER_TIER="Tier_1_Gold"
REUSE_CANONICAL_DATA="true"
{fn_body}
"""
            test_bash = test_bash.replace('"${PROJECT_ROOT}/data"', f'"{root_data}"')
            test_bash += "\nrun_preflight_subset_filter\n"
            res = subprocess.run(["bash", "-c", test_bash], cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=True)
            self.assertIn("Tier_1_Gold", res.stdout)

            out_vouchers = pd.read_csv(tables_dir / "curated_vouchers.csv")
            self.assertEqual(set(out_vouchers["determiner_tier"].unique()), {"Tier_1_Gold"})
            self.assertEqual(set(out_vouchers["catalogNumber"]), {"NCU001", "UNCC_49308"})

            for tbl_name in (
                "extracted_leaf_manifest.csv",
                "extracted_leaves_manifest.csv",
                "failed_qc_vouchers.csv",
                "voucher_reproductive_metrics.csv",
                "voucher_phenological_states.csv",
            ):
                tdf = pd.read_csv(tables_dir / tbl_name)
                self.assertEqual(set(tdf["catalogNumber"]), {"NCU001", "UNCC_49308"}, f"Failed filter on {tbl_name}")

            for sub in ("capitula",):
                link_p = run_dir / "data" / "cropped_patches" / sub
                self.assertTrue(link_p.is_symlink(), f"Expected directory-level symlink for {sub}")
                self.assertEqual(link_p.resolve(), (root_data / "cropped_patches" / sub).resolve())

            import sqlite3
            db_file = root_data / "packera_pipeline.sqlite"
            self.assertTrue(db_file.is_file())
            with sqlite3.connect(db_file) as conn:
                cohort_rows = conn.execute(
                    "SELECT catalogNumber, determiner_tier FROM run_active_cohort WHERE run_id = ?;",
                    (run_dir.name,),
                ).fetchall()
            self.assertEqual({r[0] for r in cohort_rows}, {"NCU001", "UNCC_49308"})
            self.assertEqual({r[1] for r in cohort_rows}, {"Tier_1_Gold"})

    def test_reconcile_curated_vouchers_csv_with_db_restores_requested_taxa(self):
        """Verify SQLite vouchers missing from curated_vouchers.csv are appended (taxon-restricted), preserving existing rows byte-for-byte."""
        from scripts.core.db import get_db_connection, reconcile_curated_vouchers_csv_with_db, upsert_vouchers_records

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            csv_path = tmp_path / "curated_vouchers.csv"
            original = (
                "catalogNumber,scientificName,determiner_tier,is_primary_duplicate,duplicate_count\n"
                "A1,Packera anonyma,Tier_1_Gold,True,1\n"
            )
            csv_path.write_text(original, encoding="utf-8")

            conn = get_db_connection(db_path=tmp_path / "test.sqlite", init_schema=True)
            try:
                upsert_vouchers_records(conn, [
                    {"catalogNumber": "A1", "scientificName": "Packera anonyma", "determiner_tier": "Tier_1_Gold"},
                    {"catalogNumber": "O1", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold",
                     "species_standardized": "Packera obovata", "is_primary_duplicate": 1, "duplicate_count": 2},
                    {"catalogNumber": "P1", "scientificName": "Packera paupercula", "determiner_tier": "Tier_3_Bronze",
                     "species_standardized": "Packera paupercula"},
                ])
                n_added = reconcile_curated_vouchers_csv_with_db(conn, csv_path, taxa={"Packera obovata"})
                self.assertEqual(n_added, 1)
                # Idempotent: a second pass appends nothing
                self.assertEqual(reconcile_curated_vouchers_csv_with_db(conn, csv_path, taxa={"Packera obovata"}), 0)
                # Empty taxa set restores nothing
                self.assertEqual(reconcile_curated_vouchers_csv_with_db(conn, csv_path, taxa=set()), 0)
            finally:
                conn.close()

            text = csv_path.read_text(encoding="utf-8")
            self.assertTrue(text.startswith(original))
            self.assertEqual(text[len(original):], "O1,Packera obovata,Tier_1_Gold,True,2\n")

    def test_run_pipeline_sh_preflight_restores_sqlite_taxa_and_fails_on_missing_taxa(self):
        """Verify --taxon restores requested-taxon vouchers present in SQLite but missing from curated_vouchers.csv, and fails fast on taxa with zero vouchers."""
        import subprocess
        import pandas as pd
        from scripts.core.db import get_db_connection, upsert_vouchers_records

        sh_src = (PROJECT_ROOT / "run_pipeline.sh").read_text(encoding="utf-8")
        fn_start = sh_src.index("run_preflight_subset_filter() {")
        fn_end = sh_src.index("\nrun_preflight_subset_filter\n", fn_start)
        fn_body = sh_src[fn_start:fn_end]

        def _run_preflight(run_dir, root_data, custom_taxa_csv):
            test_bash = f"""
set -euo pipefail
PROJECT_ROOT="{PROJECT_ROOT}"
RUN_DIR="{run_dir}"
PYTHON_MAIN="{sys.executable}"
FOCAL_TAXON=""
EXCLUDE_TAXON_CSV=""
INTRASPECIFIC_MODE="false"
GOLD_ANCHORS=30
DETERMINER_TIER=""
REUSE_CANONICAL_DATA="true"
CUSTOM_TAXA_CSV="{custom_taxa_csv}"
SEED_DATA_DIR="{root_data}"
START_STEP=3
{fn_body}
"""
            test_bash = test_bash.replace('"${PROJECT_ROOT}/data"', f'"{root_data}"')
            test_bash += "\nrun_preflight_subset_filter\n"
            env = {k: v for k, v in __import__("os").environ.items() if k not in ("PACKERA_DB_PATH", "PACKERA_DEMO_MODE")}
            return subprocess.run(["bash", "-c", test_bash], cwd=str(PROJECT_ROOT), capture_output=True, text=True, env=env)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            root_data = tmp_path / "root_data"
            run_dir = tmp_path / "runs" / "test_sqlite_taxa_sync"
            tables_dir = run_dir / "data" / "tables"
            tables_dir.mkdir(parents=True, exist_ok=True)
            for sub in ("capitula",):
                (root_data / "cropped_patches" / sub).mkdir(parents=True, exist_ok=True)

            # Stale CSV snapshot lacking Packera obovata
            pd.DataFrame({
                "catalogNumber": ["A1", "D1"],
                "scientificName": ["Packera anonyma", "Packera dubia"],
                "determiner_tier": ["Tier_1_Gold", "Tier_1_Gold"],
                "is_primary_duplicate": [True, True],
                "duplicate_count": [1, 1],
            }).to_csv(tables_dir / "curated_vouchers.csv", index=False)

            # Canonical SQLite vouchers include Packera obovata (from a taxon-restricted harvest)
            conn = get_db_connection(db_path=root_data / "packera_pipeline.sqlite", init_schema=True)
            try:
                upsert_vouchers_records(conn, [
                    {"catalogNumber": "A1", "scientificName": "Packera anonyma", "determiner_tier": "Tier_1_Gold"},
                    {"catalogNumber": "D1", "scientificName": "Packera dubia", "determiner_tier": "Tier_1_Gold"},
                    {"catalogNumber": "O1", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold",
                     "species_standardized": "Packera obovata", "is_primary_duplicate": 1, "duplicate_count": 1},
                    {"catalogNumber": "P1", "scientificName": "Packera paupercula", "determiner_tier": "Tier_1_Gold",
                     "species_standardized": "Packera paupercula", "is_primary_duplicate": 1, "duplicate_count": 1},
                ])
            finally:
                conn.close()

            res = _run_preflight(run_dir, root_data, "Packera anonyma,Packera dubia,Packera obovata")
            self.assertEqual(res.returncode, 0, res.stderr)
            self.assertIn("[Preflight: SQLite Voucher Sync] Restored 1 vouchers", res.stdout)
            out_vouchers = pd.read_csv(tables_dir / "curated_vouchers.csv")
            self.assertEqual(set(out_vouchers["catalogNumber"]), {"A1", "D1", "O1"})

            # A requested taxon with zero vouchers anywhere must fail fast rather than be silently dropped
            res_missing = _run_preflight(run_dir, root_data, "Packera anonyma,Packera aurea")
            self.assertNotEqual(res_missing.returncode, 0)
            self.assertIn("Packera aurea", res_missing.stderr)
            self.assertIn("zero vouchers", res_missing.stderr)

    def test_run_pipeline_sh_exclude_misidentified_and_triage_source_filtering(self):
        """Verify run_pipeline.sh --exclude-misidentified (with default ALL and explicit --triage-source) filters out Misidentification_Reassign while retaining Confirmed_Species, Hybrid_Introgressant, Subspecies_Ecotype, Ecophenotypic_Plasticity, and Putative_Cryptic_Species across curated_vouchers.csv, Step 02 manifests, and cropped_patches."""
        import subprocess
        import pandas as pd

        sh_path = PROJECT_ROOT / "run_pipeline.sh"
        sh_src = sh_path.read_text(encoding="utf-8")
        self.assertIn("--exclude-misidentified)", sh_src)
        self.assertIn("--triage-source)", sh_src)
        self.assertIn("--triage-source=*)", sh_src)
        self.assertIn("Misidentification_Reassign", sh_src)

        fn_start = sh_src.index("run_preflight_subset_filter() {")
        fn_end = sh_src.index("\nrun_preflight_subset_filter\n", fn_start)
        fn_body = sh_src[fn_start:fn_end]

        all_cats = ["CAT_CONF", "CAT_HYB", "CAT_SUB", "CAT_ECO", "CAT_CRYPT", "CAT_MISID_1", "CAT_MISID_2"]
        expected_keep_all = {"CAT_CONF", "CAT_HYB", "CAT_SUB", "CAT_ECO", "CAT_CRYPT"}
        expected_keep_run_a = {"CAT_CONF", "CAT_HYB", "CAT_SUB", "CAT_ECO", "CAT_CRYPT", "CAT_MISID_2"}

        for mode, triage_src, expected_keep in [
            ("all_default", "ALL", expected_keep_all),
            ("single_run", "run_A", expected_keep_run_a),
        ]:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_path = Path(tmpdir)
                root_data = tmp_path / "root_data"
                runs_root = tmp_path / "runs"
                run_dir = runs_root / f"test_{mode}"
                tables_dir = run_dir / "data" / "tables"
                tables_dir.mkdir(parents=True, exist_ok=True)

                # Create two historical runs with final_curation_triage.csv
                run_a_rep = runs_root / "run_A" / "outputs" / "reports"
                run_b_rep = runs_root / "run_B" / "outputs" / "reports"
                run_a_rep.mkdir(parents=True, exist_ok=True)
                run_b_rep.mkdir(parents=True, exist_ok=True)

                pd.DataFrame({
                    "catalogNumber": ["CAT_CONF", "CAT_HYB", "CAT_SUB", "CAT_ECO", "CAT_CRYPT", "CAT_MISID_1"],
                    "taxonomic_treatment_subtype": [
                        "Confirmed_Species",
                        "Hybrid_Introgressant",
                        "Subspecies_Ecotype",
                        "Ecophenotypic_Plasticity",
                        "Putative_Cryptic_Species",
                        "Misidentification_Reassign",
                    ],
                }).to_csv(run_a_rep / "final_curation_triage.csv", index=False)

                pd.DataFrame({
                    "catalogNumber": ["CAT_CONF", "CAT_MISID_2"],
                    "taxonomic_treatment_subtype": [
                        "Confirmed_Species",
                        "Misidentification_Reassign",
                    ],
                }).to_csv(run_b_rep / "final_curation_triage.csv", index=False)

                for sub in ("capitula",):
                    pdir = root_data / "cropped_patches" / sub
                    pdir.mkdir(parents=True, exist_ok=True)
                    for c in all_cats:
                        (pdir / f"{c}_cap1.jpg").write_bytes(b"jpg")

                pd.DataFrame({
                    "catalogNumber": all_cats,
                    "scientificName": ["Packera dubia"] * len(all_cats),
                    "determiner_tier": ["Tier_1_Gold"] * len(all_cats),
                    "authority_weight": [3.0] * len(all_cats),
                }).to_csv(tables_dir / "curated_vouchers.csv", index=False)

                for tbl_name in (
                    "extracted_leaf_manifest.csv",
                    "extracted_leaves_manifest.csv",
                    "failed_qc_vouchers.csv",
                    "voucher_reproductive_metrics.csv",
                    "voucher_phenological_states.csv",
                ):
                    pd.DataFrame({
                        "catalogNumber": all_cats,
                        "plant_individual_id": [0] * len(all_cats),
                    }).to_csv(tables_dir / tbl_name, index=False)

                test_bash = f"""
set -euo pipefail
PROJECT_ROOT="{PROJECT_ROOT}"
RUN_DIR="{run_dir}"
PYTHON_MAIN="{sys.executable}"
FOCAL_TAXON=""
EXCLUDE_TAXON_CSV=""
INTRASPECIFIC_MODE="false"
GOLD_ANCHORS=30
GOLD_ONLY=false
DETERMINER_TIER=""
EXCLUDE_MISIDENTIFIED=true
TRIAGE_SOURCE="{triage_src}"
REUSE_CANONICAL_DATA="true"
{fn_body}
"""
                test_bash = test_bash.replace('"${PROJECT_ROOT}/data"', f'"{root_data}"')
                test_bash += "\nrun_preflight_subset_filter\n"
                res = subprocess.run(["bash", "-c", test_bash], cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=True)
                self.assertIn("[Preflight: Exclude Misidentified]", res.stdout)

                out_vouchers = pd.read_csv(tables_dir / "curated_vouchers.csv")
                self.assertEqual(set(out_vouchers["catalogNumber"]), expected_keep)

                for tbl_name in (
                    "extracted_leaf_manifest.csv",
                    "extracted_leaves_manifest.csv",
                    "failed_qc_vouchers.csv",
                    "voucher_reproductive_metrics.csv",
                    "voucher_phenological_states.csv",
                ):
                    tdf = pd.read_csv(tables_dir / tbl_name)
                    self.assertEqual(set(tdf["catalogNumber"]), expected_keep, f"Failed filter on {tbl_name} in {mode}")

                for sub in ("capitula",):
                    link_p = run_dir / "data" / "cropped_patches" / sub
                    self.assertTrue(link_p.is_symlink(), f"Expected directory-level symlink for {sub}")
                    self.assertEqual(link_p.resolve(), (root_data / "cropped_patches" / sub).resolve())

                import sqlite3
                db_file = root_data / "packera_pipeline.sqlite"
                self.assertTrue(db_file.is_file())
                with sqlite3.connect(db_file) as conn:
                    cohort_cats = {
                        r[0]
                        for r in conn.execute(
                            "SELECT catalogNumber FROM run_active_cohort WHERE run_id = ?;",
                            (run_dir.name,),
                        ).fetchall()
                    }
                self.assertEqual(cohort_cats, expected_keep)

        # Also verify pure-SQLite --exclude-misidentified (with zero CSV triage reports on disk)
        with tempfile.TemporaryDirectory() as tmpdir_sql:
            tmp_sql = Path(tmpdir_sql)
            root_data_sql = tmp_sql / "root_data"
            root_data_sql.mkdir(parents=True, exist_ok=True)
            run_dir_sql = tmp_sql / "runs" / "test_sql_only"
            tables_dir_sql = run_dir_sql / "data" / "tables"
            tables_dir_sql.mkdir(parents=True, exist_ok=True)

            from scripts.core.db import ensure_pipeline_run, get_db_connection, upsert_rows_generic, upsert_vouchers_records

            conn_seed = get_db_connection(db_path=root_data_sql / "packera_pipeline.sqlite", init_schema=True)
            try:
                upsert_vouchers_records(
                    conn_seed,
                    [{"catalogNumber": c, "scientificName": "Packera dubia", "determiner_tier": "Tier_1_Gold"} for c in all_cats],
                )
                ensure_pipeline_run(conn_seed, run_id="hist_sql_run")
                upsert_rows_generic(
                    conn_seed,
                    table_name="curation_triage",
                    rows=[
                        {"run_id": "hist_sql_run", "catalogNumber": "CAT_CONF", "plant_individual_id": 0, "taxonomic_treatment_subtype": "Confirmed_Species"},
                        {"run_id": "hist_sql_run", "catalogNumber": "CAT_MISID_1", "plant_individual_id": 0, "taxonomic_treatment_subtype": "Misidentification_Reassign"},
                    ],
                    pk_cols=("run_id", "catalogNumber", "plant_individual_id"),
                )
            finally:
                conn_seed.close()

            pd.DataFrame({
                "catalogNumber": all_cats,
                "scientificName": ["Packera dubia"] * len(all_cats),
                "determiner_tier": ["Tier_1_Gold"] * len(all_cats),
                "authority_weight": [3.0] * len(all_cats),
            }).to_csv(tables_dir_sql / "curated_vouchers.csv", index=False)

            test_bash_sql = f"""
set -euo pipefail
PROJECT_ROOT="{PROJECT_ROOT}"
RUN_DIR="{run_dir_sql}"
PYTHON_MAIN="{sys.executable}"
FOCAL_TAXON=""
EXCLUDE_TAXON_CSV=""
INTRASPECIFIC_MODE="false"
GOLD_ANCHORS=30
GOLD_ONLY=false
DETERMINER_TIER=""
EXCLUDE_MISIDENTIFIED=true
TRIAGE_SOURCE="hist_sql_run"
REUSE_CANONICAL_DATA="false"
{fn_body}
"""
            test_bash_sql = test_bash_sql.replace('"${PROJECT_ROOT}/data"', f'"{root_data_sql}"')
            test_bash_sql += "\nrun_preflight_subset_filter\n"
            res_sql = subprocess.run(["bash", "-c", test_bash_sql], cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=True)
            self.assertIn("SQLite runs=1", res_sql.stdout)
            out_v_sql = pd.read_csv(tables_dir_sql / "curated_vouchers.csv")
            self.assertEqual(set(out_v_sql["catalogNumber"]), expected_keep_run_a)

    def test_cli_threads_and_gpu_batch_size_parsing_and_env_export(self):
        """Verify -j/--threads and --gpu-batch-size parsing, config resolution, and subprocess forwarding."""
        import importlib
        import os
        from scripts.core.config import resolve_cpu_threads, resolve_gpu_batch_size

        parser = main.build_parser()

        # Global and subcommand -j/--threads and --gpu-batch-size (plus hidden --batch-size alias)
        args_seg = parser.parse_args(["segment", "-j", "12", "--gpu-batch-size", "24", "--enable-geometric-qc"])
        self.assertEqual(args_seg.threads, 12)
        self.assertEqual(args_seg.gpu_batch_size, 24)
        self.assertTrue(args_seg.enable_geometric_qc)

        args_alias = parser.parse_args(["synthesis", "--threads", "6", "--batch-size", "16", "--ssurgo-batch-size", "50"])
        self.assertEqual(args_alias.threads, 6)
        self.assertEqual(args_alias.gpu_batch_size, 16)
        self.assertEqual(args_alias.ssurgo_batch_size, 50)

        # Verify _apply_explicit_db_arg exports PACKERA_MAX_WORKERS and PACKERA_GPU_BATCH_SIZE
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("PACKERA_MAX_WORKERS", None)
            os.environ.pop("PACKERA_GPU_BATCH_SIZE", None)
            main._apply_explicit_db_arg(args_seg)
            self.assertEqual(os.environ.get("PACKERA_MAX_WORKERS"), "12")
            self.assertEqual(os.environ.get("PACKERA_GPU_BATCH_SIZE"), "24")
            self.assertEqual(resolve_cpu_threads(None, default=4), 12)
            self.assertEqual(resolve_gpu_batch_size(None, default=16), 24)
            self.assertEqual(resolve_cpu_threads(2, default=4), 2)
            self.assertEqual(resolve_gpu_batch_size(8, default=16), 8)

        # Verify 01_voucher_harvester.py and 02_segment_and_extract.py CLI parsers
        fake_mods = {
            k: MagicMock()
            for k in (
                "pandas", "numpy", "cv2", "PIL", "PIL.Image", "scipy",
                "scipy.ndimage", "scipy.spatial", "scipy.spatial.distance",
                "sklearn", "sklearn.cluster", "sklearn.linear_model",
                "sklearn.calibration", "sklearn.model_selection", "sklearn.metrics",
                "sklearn.preprocessing",
                "matplotlib", "matplotlib.pyplot", "matplotlib.patches", "matplotlib.gridspec",
                "pygbif", "aiohttp", "tqdm",
                "torch", "torch.nn", "torch.nn.functional", "torch.utils", "torch.utils.data",
                "torchvision", "torchvision.transforms",
            )
            if k not in sys.modules
        }
        with patch.dict(sys.modules, fake_mods):
            mod_01 = importlib.import_module("scripts.data_prep.01_voucher_harvester")
            p01 = mod_01.parse_args(["--taxon", "Packera dubia", "--geographic-scope", "north_america", "--max-records", "250"])
            self.assertEqual(p01.taxa, ["Packera dubia"])
            self.assertEqual(p01.geographic_scope, "north_america")
            self.assertEqual(p01.max_records, 250)
            p01_legacy = mod_01.parse_args(["--taxa", "Packera dubia", "--include-western", "--max-records-per-taxon", "120", "--overwrite"])
            self.assertEqual(p01_legacy.taxa, ["Packera dubia"])
            self.assertTrue(p01_legacy.include_western)
            self.assertEqual(p01_legacy.max_records, 120)
            self.assertTrue(p01_legacy.force)

            mod_02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")
            with patch.object(sys, "argv", ["02_segment_and_extract.py", "-j", "10", "--gpu-batch-size", "20", "--enable-geometric-qc"]):
                p02 = mod_02.parse_args()
                self.assertEqual(p02.threads, 10)
                self.assertEqual(p02.gpu_batch_size, 20)
                self.assertTrue(p02.enable_geometric_qc)

    def test_subcommand_help_hides_suppressed_internal_and_legacy_flags(self):
        """Verify --help output exposes consolidated flags and hides internal/retired aliases."""
        parser = main.build_parser()
        subparsers_action = next(
            a for a in parser._actions if isinstance(a, argparse._SubParsersAction)
        )
        synth_help = subparsers_action.choices["synthesis"].format_help()
        self.assertIn("--threads", synth_help)
        self.assertIn("--gpu-batch-size", synth_help)
        self.assertIn("--ssurgo-batch-size", synth_help)
        for hidden_flag in ("--ssurgo-table", "--contingency-table", "--synthesis-report"):
            self.assertNotIn(hidden_flag, synth_help)

        run_all_help = subparsers_action.choices["run-all"].format_help()
        self.assertIn("--taxon", run_all_help)
        self.assertIn("--threads", run_all_help)
        self.assertIn("--gpu-batch-size", run_all_help)
        for hidden_flag in ("--overwrite", "--include-western", "--allometry-report"):
            self.assertNotIn(hidden_flag, run_all_help)

    def test_run_pipeline_sh_help_and_compute_flags(self):
        """Verify ./run_pipeline.sh --help exits 0 and documents -j/--threads and --gpu-batch-size."""
        import subprocess
        res = subprocess.run(
            ["bash", str(PROJECT_ROOT / "run_pipeline.sh"), "--help"],
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertEqual(res.returncode, 0)
        self.assertIn("-j, --threads", res.stdout)
        self.assertIn("--gpu-batch-size", res.stdout)
        self.assertIn("--taxon", res.stdout)
        for active_flag in ("--skip-harvest", "--resegment", "--purge-segments", "--purge-old-data", "--gold-only", "--auto-neighbors"):
            self.assertIn(active_flag, res.stdout)

    def test_purge_segments_cli_and_purge_segments_and_runs_execution(self):
        """Verify --purge-segments / --purge-old-data CLI parsing and purge_segments_and_runs execution."""
        import importlib
        import sqlite3
        from scripts.core.db import get_db_connection, purge_segments_and_runs

        parser = main.build_parser()

        # Top-level flag without subcommand
        args_top = parser.parse_args(["--purge-segments"])
        self.assertTrue(args_top.purge_segments)
        self.assertIsNone(args_top.subcommand)

        args_top_alias = parser.parse_args(["--purge-old-data", "--no-purge-runs"])
        self.assertTrue(args_top_alias.purge_segments)
        self.assertFalse(args_top_alias.purge_runs)

        # Standalone subcommand and alias
        args_sub = parser.parse_args(["purge-segments"])
        self.assertEqual(args_sub.subcommand, "purge-segments")
        args_sub_alias = parser.parse_args(["purge-old-data", "--no-purge-runs"])
        self.assertEqual(args_sub_alias.subcommand, "purge-old-data")
        self.assertFalse(args_sub_alias.purge_runs)

        # Subcommand flag on segment and run-all
        args_seg_purge = parser.parse_args(["segment", "--purge-segments"])
        self.assertTrue(args_seg_purge.purge_segments)
        args_run_all_purge = parser.parse_args(["run-all", "--purge-old-data"])
        self.assertTrue(args_run_all_purge.purge_segments)

        # 02_segment_and_extract.py flag parsing
        fake_mods = {
            k: MagicMock()
            for k in (
                "pandas", "numpy", "cv2", "PIL", "PIL.Image", "scipy",
                "scipy.ndimage", "scipy.spatial", "scipy.spatial.distance",
                "sklearn", "sklearn.cluster", "sklearn.linear_model",
                "sklearn.calibration", "sklearn.model_selection", "sklearn.metrics",
                "sklearn.preprocessing",
                "matplotlib", "matplotlib.pyplot", "matplotlib.patches", "matplotlib.gridspec",
                "pygbif", "aiohttp", "tqdm",
                "torch", "torch.nn", "torch.nn.functional", "torch.utils", "torch.utils.data",
                "torchvision", "torchvision.transforms",
            )
            if k not in sys.modules
        }
        with patch.dict(sys.modules, fake_mods):
            mod_02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")
            with patch.object(sys, "argv", ["02_segment_and_extract.py", "--purge-segments"]):
                p02 = mod_02.parse_args()
                self.assertTrue(p02.purge_segments)

        # End-to-end execution of purge_segments_and_runs in isolated temp workspace
        with tempfile.TemporaryDirectory() as tmpdir:
            ws = Path(tmpdir)
            (ws / "data" / "raw_vouchers").mkdir(parents=True)
            (ws / "data" / "raw_vouchers" / "NCU0001.webp").write_bytes(b"WEBP")
            (ws / "data" / "tables").mkdir(parents=True)
            (ws / "data" / "tables" / "curated_vouchers.csv").write_text("catalogNumber\nNCU0001\n", encoding="utf-8")
            (ws / "data" / "tables" / "extracted_leaf_traits.csv").write_text("catalogNumber\nNCU0001\n", encoding="utf-8")
            (ws / "data" / "contours").mkdir(parents=True)
            (ws / "data" / "contours" / "NCU0001_leaf_1.csv").write_text("x,y\n0,0\n", encoding="utf-8")
            (ws / "data" / "masks" / "tier1_intact").mkdir(parents=True)
            (ws / "data" / "masks" / "tier1_intact" / "NCU0001_leaf_1.png").write_bytes(b"PNG")
            (ws / "data" / "ect").mkdir(parents=True)
            (ws / "data" / "ect" / "NCU0001_leaf_1_ect.npy").write_bytes(b"NPY")
            (ws / "data" / "cropped_patches" / "tier1_intact").mkdir(parents=True)
            (ws / "data" / "cropped_patches" / "tier1_intact" / "NCU0001_leaf_1.webp").write_bytes(b"WEBP")
            (ws / "runs" / "20260101_000000").mkdir(parents=True)
            (ws / "runs" / "20260101_000000" / "run.log").write_text("log", encoding="utf-8")

            db_file = ws / "data" / "packera_pipeline.sqlite"
            conn = get_db_connection(db_file, init_schema=True)
            conn.execute(
                "INSERT INTO vouchers (catalogNumber, scientificName, qc_status) VALUES ('NCU0001', 'Packera aurea', 'PASS')"
            )
            conn.execute(
                "INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('NCU0001', 1)"
            )
            conn.execute(
                "INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, assigned_tier, contour_xy_json) VALUES ('NCU0001', 1, 'L1', 'Tier_1_Direct', '[[0,0]]')"
            )
            conn.commit()
            conn.close()

            summary = purge_segments_and_runs(workspace_root=ws, db_path=db_file, purge_runs=True, vacuum=True)
            self.assertGreaterEqual(summary["removed_files"], 5)
            self.assertEqual(summary["purged_db_rows"], 2)

            # Step 01 preserved
            self.assertTrue((ws / "data" / "raw_vouchers" / "NCU0001.webp").exists())
            self.assertTrue((ws / "data" / "tables" / "curated_vouchers.csv").exists())
            # Step 02-07 files removed and directories recreated
            self.assertFalse((ws / "data" / "tables" / "extracted_leaf_traits.csv").exists())
            self.assertEqual(list((ws / "data" / "contours").glob("*.csv")), [])
            self.assertEqual(list((ws / "data" / "ect").glob("*.npy")), [])
            self.assertTrue((ws / "data" / "masks" / "tier1_intact").is_dir())
            self.assertEqual(list((ws / "data" / "masks" / "tier1_intact").glob("*")), [])
            self.assertEqual(list((ws / "runs").iterdir()), [])

            # SQLite check: vouchers preserved with qc_status=NULL, Step 02 tables empty
            conn_check = sqlite3.connect(str(db_file))
            v_row = conn_check.execute("SELECT catalogNumber, qc_status FROM vouchers").fetchone()
            self.assertEqual(v_row, ("NCU0001", None))
            self.assertEqual(conn_check.execute("SELECT COUNT(*) FROM plant_individuals").fetchone()[0], 0)
            self.assertEqual(conn_check.execute("SELECT COUNT(*) FROM extracted_leaves").fetchone()[0], 0)
            conn_check.close()


if __name__ == "__main__":
    unittest.main()

