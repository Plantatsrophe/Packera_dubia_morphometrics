"""
Unit tests for centralized configuration and root CLI runner (main.py).
"""

import sys
from pathlib import Path
import pytest

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import main
from config import load_config, DEFAULT_CONFIG_PATH


def test_load_config_structure():
    """Verify config.yaml loads with all required pipeline sections."""
    cfg = load_config(DEFAULT_CONFIG_PATH)
    assert "paths" in cfg
    assert "taxonomy" in cfg
    assert "harvesting" in cfg
    assert "segmentation" in cfg
    assert "morphometrics" in cfg

    # Verify resolved paths
    assert "resolved_paths" in cfg
    assert cfg["resolved_paths"]["raw_vouchers_dir"].is_absolute()
    assert cfg["resolved_paths"]["curated_vouchers_csv"].is_absolute()

    # Verify key hyperparameters
    assert cfg["harvesting"]["max_records_per_taxon"] == 5000
    assert cfg["segmentation"]["score_thresh"] == 0.40
    assert cfg["morphometrics"]["harmonics"] == 12


def test_build_parser_subcommands():
    """Verify CLI parser parses all subcommands correctly."""
    parser = main.build_parser()

    # harvest
    args_harvest = parser.parse_args(["harvest", "--max-records", "100", "--download-images"])
    assert args_harvest.subcommand == "harvest"
    assert args_harvest.max_records == 100
    assert args_harvest.download_images is True

    # segment
    args_seg = parser.parse_args(["segment", "--device", "cpu", "--limit", "10"])
    assert args_seg.subcommand == "segment"
    assert args_seg.device == "cpu"
    assert args_seg.limit == 10

    # morphometrics
    args_morph = parser.parse_args(["morphometrics", "--harmonics", "15", "--num-pcs", "6"])
    assert args_morph.subcommand == "morphometrics"
    assert args_morph.harmonics == 15
    assert args_morph.num_pcs == 6

    # run-all
    args_all = parser.parse_args(["run-all", "--download-images"])
    assert args_all.subcommand == "run-all"
    assert args_all.download_images is True


def test_gpu_availability_check():
    """Verify GPU availability check runs and returns boolean without unhandled exception."""
    result = main.check_gpu_availability(warn_only=True)
    assert isinstance(result, bool)


def test_verify_file_exists_fails_cleanly():
    """Verify preflight check raises SystemExit on missing upstream file."""
    with pytest.raises(SystemExit):
        main.verify_file_exists(
            Path("/non/existent/path/voucher_table.csv"),
            description="test table",
            hint_cmd="python main.py harvest"
        )


def test_verify_dir_has_files_fails_cleanly(tmp_path):
    """Verify preflight check raises SystemExit on empty directory."""
    empty_dir = tmp_path / "empty_dir"
    empty_dir.mkdir()
    with pytest.raises(SystemExit):
        main.verify_dir_has_files(
            empty_dir,
            pattern="*.csv",
            description="test empty dir",
            hint_cmd="python main.py segment"
        )
