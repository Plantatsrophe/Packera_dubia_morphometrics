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

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            dummy_weights = tmp_path / "lm2_packera_pcd_finetuned.pth"
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

            for sub in ("rosettes_dense", "basal_leaves_raw", "capitula"):
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

            # Create cropped patches for Gold, Silver, and Bronze vouchers (including underscored catalogNumber)
            (root_data / "cropped_patches" / "rosettes_dense" / "NCU001_p0_rosette.jpg").write_bytes(b"jpg")
            (root_data / "cropped_patches" / "rosettes_dense" / "UNCC_49308_p0_rosette.jpg").write_bytes(b"jpg")
            (root_data / "cropped_patches" / "rosettes_dense" / "US002_p0_rosette.jpg").write_bytes(b"jpg")
            (root_data / "cropped_patches" / "basal_leaves_raw" / "NCU001_leaf1.jpg").write_bytes(b"jpg")
            (root_data / "cropped_patches" / "basal_leaves_raw" / "UNCC_49308_leaf2.jpg").write_bytes(b"jpg")
            (root_data / "cropped_patches" / "basal_leaves_raw" / "NY003_leaf1.jpg").write_bytes(b"jpg")
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

            ros_files = {p.name for p in (run_dir / "data" / "cropped_patches" / "rosettes_dense").iterdir()}
            leaf_files = {p.name for p in (run_dir / "data" / "cropped_patches" / "basal_leaves_raw").iterdir()}
            cap_files = {p.name for p in (run_dir / "data" / "cropped_patches" / "capitula").iterdir()}
            self.assertEqual(ros_files, {"NCU001_p0_rosette.jpg", "UNCC_49308_p0_rosette.jpg"})
            self.assertEqual(leaf_files, {"NCU001_leaf1.jpg", "UNCC_49308_leaf2.jpg"})
            self.assertEqual(cap_files, {"NCU001_cap1.jpg"})

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

                for sub in ("rosettes_dense", "basal_leaves_raw", "capitula"):
                    pdir = root_data / "cropped_patches" / sub
                    pdir.mkdir(parents=True, exist_ok=True)
                    for c in all_cats:
                        suffix = "p0_rosette.jpg" if sub == "rosettes_dense" else ("leaf1.jpg" if sub == "basal_leaves_raw" else "cap1.jpg")
                        (pdir / f"{c}_{suffix}").write_bytes(b"jpg")

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

                ros_files = {p.name for p in (run_dir / "data" / "cropped_patches" / "rosettes_dense").iterdir()}
                leaf_files = {p.name for p in (run_dir / "data" / "cropped_patches" / "basal_leaves_raw").iterdir()}
                cap_files = {p.name for p in (run_dir / "data" / "cropped_patches" / "capitula").iterdir()}
                self.assertEqual(ros_files, {f"{c}_p0_rosette.jpg" for c in expected_keep})
                self.assertEqual(leaf_files, {f"{c}_leaf1.jpg" for c in expected_keep})
                self.assertEqual(cap_files, {f"{c}_cap1.jpg" for c in expected_keep})


if __name__ == "__main__":
    unittest.main()

