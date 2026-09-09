#!/usr/bin/env python3
"""
===============================================================================
Script: main.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unified command-line interface (CLI) runner orchestrating the end-to-end
    botanical morphometrics pipeline:
      - harvest: Queries GBIF occurrences, scores determiner tiers, downloads imagery.
      - segment: PointRend deep segmentation, bilateral midrib reflection, contour export.
      - morphometrics: Rscript EFA (Momocs) and GMM/CDA clustering (MorphoTools2).
      - synthesis: Canonical Phases 5-7 (Vision XAI, Spatial RF, & Triage Dashboard).
      - run-all: Executes the complete production workflow sequentially.
===============================================================================
"""

from __future__ import annotations

import argparse
import importlib
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.artifact_manager import ensure_model_weights
from scripts.core.config import PipelineConfig
from scripts.core.logger import setup_logging

logger = setup_logging(name="PackeraPipeline")



# =============================================================================
# Defensive Preflight Sanity Checks
# =============================================================================

def check_gpu_availability(warn_only: bool = True) -> bool:
    """Checks whether PyTorch detects an active CUDA GPU device."""
    try:
        import torch
        if torch.cuda.is_available():
            device_name = torch.cuda.get_device_name(0)
            logger.info(f"GPU Preflight: CUDA device detected: {device_name}")
            return True
        else:
            msg = "GPU Preflight: No CUDA-capable device detected by PyTorch."
            if warn_only:
                logger.warning(f"{msg} Operations will run on CPU, which may be slower.")
            else:
                logger.error(msg)
            return False
    except ImportError:
        logger.warning("GPU Preflight: PyTorch not found in current environment.")
        return False


def verify_file_exists(file_path: Path, description: str, hint_cmd: Optional[str] = None) -> None:
    """Verifies that a required input file exists, failing early with clear instructions."""
    if not file_path.exists():
        err = f"Missing required {description} at: {file_path}"
        if hint_cmd:
            err += f"\n  -> Upstream phase skipped or incomplete! Please run: `{hint_cmd}`"
        logger.error(err)
        sys.exit(1)


def verify_dir_has_files(dir_path: Path, pattern: str, description: str, hint_cmd: Optional[str] = None) -> None:
    """Verifies that a required directory exists and is non-empty."""
    if not dir_path.exists() or not any(dir_path.glob(pattern)):
        err = f"Missing or empty {description} directory at: {dir_path} (pattern: {pattern})"
        if hint_cmd:
            err += f"\n  -> Upstream phase skipped or incomplete! Please run: `{hint_cmd}`"
        logger.error(err)
        sys.exit(1)


def check_environment(args: Optional[argparse.Namespace] = None) -> bool:
    """Performs comprehensive diagnostic checks across Python, CUDA, R,

    configuration files, filesystem directories, and model checkpoints.

    Returns:
        bool: True if all critical checks pass, False if any critical check fails.
    """
    strict = getattr(args, "strict", False) if args else False
    custom_config = getattr(args, "config", None) if args else None
    critical_failure = False
    warning_count = 0

    print("=" * 79)
    print("Packera dubia Morphometrics Pipeline: Environment & Dependency Diagnostic")
    print("=" * 79)

    # -------------------------------------------------------------------------
    # 1. Python & CUDA Environment
    # -------------------------------------------------------------------------
    print("\n[1/4] Python & CUDA Environment:")

    # Python version check (>= 3.10)
    py_ver_str = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    if sys.version_info >= (3, 10):
        print(f"  [PASS] Python version: {py_ver_str} (>= 3.10 required)")
    else:
        print(f"  [FAIL] Python version: {py_ver_str} is unsupported (>= 3.10 required).")
        critical_failure = True

    # PyTorch & CUDA check
    try:
        import torch
        if torch.cuda.is_available():
            dev_name = torch.cuda.get_device_name(0)
            props = torch.cuda.get_device_properties(0)
            vram_gb = props.total_memory / (1024 ** 3)
            print(f"  [PASS] PyTorch CUDA acceleration: {dev_name} ({vram_gb:.2f} GB VRAM)")
        else:
            print("  [WARN] CUDA acceleration unavailable; PyTorch operations will run on CPU.")
            warning_count += 1
            if strict:
                critical_failure = True
    except ImportError:
        print("  [WARN] PyTorch is not installed in current environment; cannot probe CUDA.")
        warning_count += 1
        if strict:
            critical_failure = True

    # Core package imports
    core_packages = [
        ("cv2", "OpenCV"),
        ("torch", "PyTorch"),
        ("cleanlab", "Cleanlab"),
        ("pygbif", "pygbif"),
        ("yaml", "PyYAML"),
    ]
    for pkg_module, pkg_display in core_packages:
        try:
            mod = importlib.import_module(pkg_module)
            ver = getattr(mod, "__version__", "installed")
            print(f"  [PASS] Python package '{pkg_module}' ({pkg_display}): version {ver}")
        except ImportError:
            print(f"  [FAIL] Python package '{pkg_module}' ({pkg_display}) is NOT installed.")
            print(f"         Remediation: pip install {pkg_module}")
            critical_failure = True

    # -------------------------------------------------------------------------
    # 2. R Runtime & Statistical Packages
    # -------------------------------------------------------------------------
    print("\n[2/4] R Runtime & Statistical Packages:")
    rscript_bin = shutil.which("Rscript")
    required_r_pkgs = ['Momocs', 'mclust', 'MorphoTools2', 'spatialRF', 'terra', 'tidyverse', 'optparse']
    r_pkg_str = ", ".join(f"'{p}'" for p in required_r_pkgs)

    if not rscript_bin:
        print("  [FAIL] 'Rscript' binary not found in system $PATH.")
        print("         R (>= 4.3) is required for Momocs, MorphoTools2, and spatialRF.")
        print("         Remediation:")
        print("           1. Install R (>= 4.3) via system package manager (e.g., apt-get install r-base r-base-dev).")
        print(f"           2. Install required packages:")
        print(f"              Rscript -e \"install.packages(c({r_pkg_str}), repos='https://cloud.r-project.org')\"")
        critical_failure = True
    else:
        try:
            r_ver_proc = subprocess.run(
                [rscript_bin, "-e", "cat(as.character(getRversion()))"],
                capture_output=True, text=True, timeout=5,
            )
            r_ver_text = r_ver_proc.stdout.strip()
            if r_ver_proc.returncode == 0 and r_ver_text and r_ver_text[0].isdigit():
                print(f"  [PASS] Rscript executable: {rscript_bin} (R {r_ver_text})")
            else:
                print(f"  [PASS] Rscript executable: {rscript_bin}")
        except Exception:
            print(f"  [PASS] Rscript executable: {rscript_bin}")

        r_probe_cmd = (
            "pkgs <- c('Momocs', 'mclust', 'MorphoTools2', 'spatialRF', 'terra', 'tidyverse', 'optparse'); "
            "missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; "
            "if(length(missing)>0) { cat('MISSING:', paste(missing, collapse=','), fill=TRUE); quit(status=1) } "
            "else { cat('ALL_INSTALLED', fill=TRUE); quit(status=0) }"
        )
        try:
            probe_proc = subprocess.run(
                [rscript_bin, "-e", r_probe_cmd],
                capture_output=True, text=True, timeout=10,
            )
            if probe_proc.returncode == 0 and "ALL_INSTALLED" in probe_proc.stdout:
                print(f"  [PASS] Required R packages verified: {', '.join(required_r_pkgs)}")
            else:
                combined_out = probe_proc.stdout + " " + probe_proc.stderr
                missing_pkgs: List[str] = []
                for line in combined_out.splitlines():
                    if "MISSING:" in line:
                        pkgs_part = line.split("MISSING:", 1)[1].strip()
                        missing_pkgs = [p.strip() for p in pkgs_part.split(",") if p.strip()]
                        break
                if not missing_pkgs:
                    missing_pkgs = required_r_pkgs

                print(f"  [FAIL] Missing R package(s): {', '.join(missing_pkgs)}")
                missing_quoted = ", ".join(f"'{p}'" for p in missing_pkgs)
                print("         Remediation:")
                print(f"           Rscript -e \"install.packages(c({missing_quoted}), repos='https://cloud.r-project.org')\"")
                critical_failure = True
        except subprocess.TimeoutExpired:
            print("  [FAIL] R package validation probe timed out (>10s).")
            critical_failure = True
        except Exception as exc:
            print(f"  [FAIL] Error executing R package probe: {exc}")
            critical_failure = True

    # -------------------------------------------------------------------------
    # 3. Directory & Configuration Verification
    # -------------------------------------------------------------------------
    print("\n[3/4] Directory & Configuration Verification:")
    cfg_file = Path(custom_config or (PROJECT_ROOT / "config" / "config.yaml"))
    parsed_config: Optional[Dict[str, Any]] = None

    if not cfg_file.exists():
        print(f"  [FAIL] Pipeline configuration file not found at: {cfg_file}")
        critical_failure = True
    else:
        try:
            with open(cfg_file, "r", encoding="utf-8") as f:
                parsed_config = yaml.safe_load(f)
            if not isinstance(parsed_config, dict) or "paths" not in parsed_config:
                print(f"  [FAIL] Configuration at {cfg_file} is invalid or missing 'paths' section.")
                critical_failure = True
            else:
                rel_cfg = cfg_file.relative_to(PROJECT_ROOT) if cfg_file.is_relative_to(PROJECT_ROOT) else cfg_file
                print(f"  [PASS] Configuration file parsed & valid: {rel_cfg}")
        except Exception as exc:
            print(f"  [FAIL] Error parsing YAML configuration at {cfg_file}: {exc}")
            critical_failure = True

    target_dirs = [
        "data/raw_vouchers",
        "data/tables",
        "data/contours",
        "outputs/figures",
    ]
    for dir_rel in target_dirs:
        dir_path = PROJECT_ROOT / dir_rel
        try:
            dir_path.mkdir(parents=True, exist_ok=True)
            probe_file = dir_path / ".perm_check_probe"
            probe_file.touch()
            probe_file.unlink()
            print(f"  [PASS] Directory exists and writable: {dir_rel}/")
        except Exception as exc:
            print(f"  [FAIL] Directory missing or write-permission denied: {dir_rel}/ ({exc})")
            critical_failure = True

    # -------------------------------------------------------------------------
    # 4. Model Checkpoint & Sub-Environments
    # -------------------------------------------------------------------------
    print("\n[4/4] Model Checkpoint & Sub-Environments:")
    venv_lm2 = PROJECT_ROOT / ".venv_LM2"
    venv_lm2_python = venv_lm2 / "bin" / "python"
    if venv_lm2.is_dir() and (venv_lm2_python.exists() or (venv_lm2 / "Scripts" / "python.exe").exists()):
        print(f"  [PASS] LeafMachine2 virtual environment found: {venv_lm2.name}/")
    else:
        print(f"  [WARN] LeafMachine2 environment (.venv_LM2) not found at {venv_lm2}")
        warning_count += 1
        if strict:
            critical_failure = True

    model_rel = (
        parsed_config.get("paths", {}).get("model_weights", "models/lm2_packera_pcd_finetuned.pth")
        if parsed_config
        else "models/lm2_packera_pcd_finetuned.pth"
    )
    model_path = PROJECT_ROOT / model_rel
    if not (model_path.exists() and model_path.is_file() and model_path.stat().st_size > 0):
        print(f"  [INFO] Model checkpoint missing at {model_path}. Invoking automated artifact manager...")
        try:
            cfg_obj = PipelineConfig.from_yaml(cfg_file) if cfg_file.exists() else (parsed_config or {})
            ensure_model_weights(cfg_obj)
        except Exception as exc:
            print(f"  [WARN] Automatic weights acquisition failed: {exc}")

    if model_path.exists() and model_path.is_file() and model_path.stat().st_size > 0:
        size_mb = model_path.stat().st_size / (1024 * 1024)
        print(f"  [PASS] Fine-tuned model checkpoint verified: {model_rel} ({size_mb:.2f} MB)")
    else:
        print(f"  [FAIL] Fine-tuned model checkpoint not found at: {model_path}")
        print("         Remediation: Run `python main.py download-weights` or place trained weights manually.")
        critical_failure = True

    # -------------------------------------------------------------------------
    # Readiness Verdict
    # -------------------------------------------------------------------------
    print("=" * 79)
    if critical_failure:
        print("Readiness Verdict: [FAIL] Environment is NOT ready for production.")
        print("Please resolve the [FAIL] dependencies listed above before starting runs.")
        print("=" * 79)
        return False
    else:
        if warning_count > 0:
            print(f"Readiness Verdict: [PASS] Ready with {warning_count} warning(s).")
        else:
            print("Readiness Verdict: [PASS] Complete pipeline environment verified & ready!")
        print("=" * 79)
        return True


def run_check_env(args: argparse.Namespace, cfg: Optional[PipelineConfig] = None) -> None:
    """CLI handler executing environment diagnostics and terminating with status code."""
    success = check_environment(args)
    sys.exit(0 if success else 1)


def run_download_weights(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Downloads and verifies fine-tuned model checkpoint weights."""
    logger.info("=== Checking & Downloading Model Weights Checkpoint ===")
    try:
        saved_path = ensure_model_weights(cfg, force=getattr(args, "force", False))
        logger.info(f"Model weights checkpoint ready at: {saved_path}")
    except Exception as exc:
        logger.error(f"Failed to acquire model weights: {exc}")
        sys.exit(1)


# =============================================================================
# Phase Execution Handlers
# =============================================================================

def run_harvest(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Executes Phase 1: Voucher harvesting and authority tier stratification."""
    logger.info("=== Starting Phase 1: Voucher Ingestion & Authority Stratification ===")
    
    script_path = PROJECT_ROOT / "scripts" / "data_prep" / "01_voucher_harvester.py"
    if not script_path.exists():
        script_path = PROJECT_ROOT / "scripts" / "pipeline" / "01_voucher_harvester.py"
    verify_file_exists(script_path, "harvesting script")

    out_dir = Path(args.out_dir or cfg["paths"]["raw_vouchers_dir"])
    output_csv = Path(args.output_csv or cfg["paths"]["curated_vouchers_csv"])
    out_dir.mkdir(parents=True, exist_ok=True)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, str(script_path),
        "--out-dir", str(out_dir),
        "--output-csv", str(output_csv),
        "--max-records", str(args.max_records or cfg["harvesting"]["max_records_per_taxon"]),
        "--max-uncertainty", str(args.max_uncertainty or cfg["harvesting"]["max_uncertainty_meters"]),
        "--min-megapixels", str(args.min_megapixels or cfg["harvesting"]["min_megapixels"]),
        "--min-file-size-kb", str(args.min_file_size_kb or cfg["harvesting"]["min_file_size_kb"]),
        "--concurrency", str(args.concurrency or cfg["harvesting"]["download_concurrency"]),
    ]
    if args.taxa:
        cmd.extend(["--taxa", *args.taxa])
    elif "taxonomy" in cfg and "target_taxa" in cfg["taxonomy"]:
        cmd.extend(["--taxa", *cfg["taxonomy"]["target_taxa"]])

    if args.download_images:
        cmd.append("--download-images")
    if getattr(args, "force", False):
        cmd.append("--force")
    if args.verbose:
        cmd.append("--verbose")

    logger.info(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    verify_file_exists(output_csv, "harvested vouchers CSV output")
    logger.info("Phase 1 completed successfully.")


def get_lm2_python_executable() -> Path:
    """Resolves the Python interpreter executable for LeafMachine2 execution.

    Checks for:
      - Unix/Ubuntu: .venv_LM2/bin/python
      - Windows: .venv_LM2/Scripts/python.exe

    Returns:
        Path: Path to the LM2 virtual environment Python interpreter if found,
        otherwise falls back to sys.executable and logs a warning advising setup
        via setup_leafmachine2.sh.
    """
    venv_lm2 = PROJECT_ROOT / ".venv_LM2"
    candidates = [
        venv_lm2 / "bin" / "python",
        venv_lm2 / "Scripts" / "python.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    logger.warning(
        "LeafMachine2 dedicated virtual environment (.venv_LM2) not found at %s. "
        "Falling back to current active interpreter (%s). "
        "Please run 'bash setup_leafmachine2.sh' to configure the dedicated LM2 environment.",
        venv_lm2,
        sys.executable,
    )
    return Path(sys.executable)


def run_segment(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Executes Phase 2: PointRend segmentation, midrib reflection, and contour extraction."""
    logger.info("=== Starting Phase 2: Segmentation & Geometric Leaf Extraction ===")
    
    script_path = PROJECT_ROOT / "scripts" / "pipeline" / "02_segment_and_extract.py"
    verify_file_exists(script_path, "segmentation script")

    vouchers_csv = Path(args.vouchers or cfg["paths"]["curated_vouchers_csv"])
    verify_file_exists(vouchers_csv, "curated vouchers metadata table", "python main.py harvest")

    if args.weights:
        model_weights = Path(args.weights)
    else:
        try:
            model_weights = ensure_model_weights(cfg)
        except Exception as exc:
            logger.warning(f"Automated weights acquisition could not complete: {exc}")
            model_weights = Path(cfg["paths"]["model_weights"])
    verify_file_exists(model_weights, "PointRend model weights checkpoint")

    device = args.device or cfg["segmentation"].get("device", "cuda")
    if device == "cuda":
        check_gpu_availability(warn_only=True)

    lm2_python = get_lm2_python_executable()
    try:
        display_path = lm2_python.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        display_path = str(lm2_python)
    logger.info(f"Executing LeafMachine2 segmentation via: {display_path}")

    out_dir = Path(args.output_dir or (PROJECT_ROOT / "data"))
    cmd = [
        str(lm2_python), str(script_path),
        "--vouchers", str(vouchers_csv),
        "--model-weights", str(model_weights),
        "--output-dir", str(out_dir),
        "--device", device,
        "--min-solidity", str(args.min_solidity or cfg["segmentation"]["min_solidity"]),
        "--min-ucs", str(args.min_ucs or cfg["segmentation"]["min_ucs"]),
        "--score-thresh", str(args.score_thresh or cfg["segmentation"]["score_thresh"]),
    ]
    if args.limit:
        cmd.extend(["--limit", str(args.limit)])
    if getattr(args, "force", False):
        cmd.append("--force")

    logger.info(f"Running command: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        logger.error(f"LeafMachine2 segmentation failed with exit code {exc.returncode}")
        sys.exit(exc.returncode)

    contours_dir = Path(cfg["paths"]["contours_dir"])
    verify_dir_has_files(contours_dir, "*.csv", "extracted leaf contours")
    logger.info("Phase 2 completed successfully.")


def run_morphometrics(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Executes Phase 3: EFA, GMM, and passive sample CDA via R scripts."""
    logger.info("=== Starting Phase 3: Morphometrics, EFA & Discriminant Analysis ===")

    rscript_bin = shutil.which("Rscript")
    if not rscript_bin:
        logger.error(
            "Rscript binary not found in system PATH.\n"
            "Phase 3 requires R (>= 4.3) with Momocs, MorphoTools2, and mclust installed."
        )
        sys.exit(1)

    contours_dir = Path(getattr(args, "input", None) or args.contours_dir or cfg["paths"]["contours_dir"])
    verify_dir_has_files(contours_dir, "*.csv", "leaf contours", "python main.py segment")

    vouchers_csv = Path(args.vouchers or cfg["paths"]["curated_vouchers_csv"])
    verify_file_exists(vouchers_csv, "curated vouchers table", "python main.py harvest")

    manifest_csv = Path(getattr(args, "manifest", None) or "data/tables/extracted_leaf_manifest.csv")

    efa_script = PROJECT_ROOT / "scripts" / "morphometrics" / "03_fourier_extractor.R"
    verify_file_exists(efa_script, "EFA extractor R script")

    efa_out = Path(args.harmonics_out or cfg["paths"]["leaf_efa_harmonics_csv"])
    efa_out.parent.mkdir(parents=True, exist_ok=True)

    harmonics = args.harmonics or cfg["morphometrics"]["harmonics"]
    num_pcs = args.num_pcs or cfg["morphometrics"]["num_pcs"]

    # Step 3A: Fourier EFA
    cmd_efa = [
        rscript_bin, str(efa_script),
        "--input", str(contours_dir),
        "--manifest", str(manifest_csv),
        "--vouchers", str(vouchers_csv),
        "--output", str(efa_out),
        "--harmonics", str(harmonics),
        "--num-pcs", str(num_pcs),
    ]
    logger.info(f"Running EFA: {' '.join(cmd_efa)}")
    subprocess.run(cmd_efa, check=True)
    verify_file_exists(efa_out, "EFA harmonics table")

    # Step 3B: GMM & CDA
    gmm_script = PROJECT_ROOT / "scripts" / "morphometrics" / "04_gmm_morphotools.R"
    verify_file_exists(gmm_script, "GMM / MorphoTools2 R script")

    flags_out = Path(args.flags_out or cfg["paths"]["morphometric_flags_csv"])
    plot_out = Path(args.plot_out or cfg["paths"]["cda_biplot_pdf"])
    report_out = Path(args.report_out or cfg["paths"]["gmm_report_csv"])
    flags_out.parent.mkdir(parents=True, exist_ok=True)
    plot_out.parent.mkdir(parents=True, exist_ok=True)
    report_out.parent.mkdir(parents=True, exist_ok=True)

    max_k = args.max_k or cfg["morphometrics"]["max_k"]

    cmd_gmm = [
        rscript_bin, str(gmm_script),
        "--vouchers", str(vouchers_csv),
        "--harmonics", str(efa_out),
        "--output-flags", str(flags_out),
        "--output-plot", str(plot_out),
        "--output-report", str(report_out),
        "--max-k", str(max_k),
        "--num-pcs", str(num_pcs),
    ]
    logger.info(f"Running GMM/CDA: {' '.join(cmd_gmm)}")
    subprocess.run(cmd_gmm, check=True)
    verify_file_exists(flags_out, "morphometric flags CSV")
    logger.info("Phase 3 completed successfully.")


def run_synthesis(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Executes Phases 5, 6, and 7: Vision XAI, Spatial Macroecology, and Synthesis Triage."""
    logger.info("=== Starting Synthesis Workflow (Phases 5, 6, 7) ===")

    rscript_bin = shutil.which("Rscript")
    if not rscript_bin:
        logger.error(
            "Rscript binary not found in system PATH.\n"
            "Spatial macroecology and synthesis dashboard require R (>= 4.3) with spatialRF, ENMTools, and terra."
        )
        sys.exit(1)

    vouchers_csv = Path(getattr(args, "vouchers", None) or cfg["paths"]["curated_vouchers_csv"])
    verify_file_exists(vouchers_csv, "curated vouchers table", "python main.py harvest")

    # -------------------------------------------------------------------------
    # Stage 5: DINOv2 Deep Vision Feature Extraction & Cleanlab XAI Audit
    # -------------------------------------------------------------------------
    xai_script = PROJECT_ROOT / "scripts" / "analysis" / "05_cleanlab_vision_xai.py"
    verify_file_exists(xai_script, "Cleanlab Vision XAI script")

    rosette_dir = Path(getattr(args, "rosette_dir", None) or "data/cropped_patches")
    audit_csv = Path(getattr(args, "vision_audit", None) or "data/tables/label_noise_audit.csv")
    xai_fig = Path(getattr(args, "xai_plot", None) or "outputs/figures/GradCAM_audit_panel.png")
    cleanlab_thresh = getattr(args, "cleanlab_threshold", 0.85)

    cmd_xai = [
        sys.executable, str(xai_script),
        "--rosette-dir", str(rosette_dir),
        "--vouchers-csv", str(vouchers_csv),
        "--output-csv", str(audit_csv),
        "--output-figure", str(xai_fig),
        "--cleanlab-threshold", str(cleanlab_thresh),
    ]
    if getattr(args, "export_figures", True):
        cmd_xai.append("--export-figures")
    else:
        cmd_xai.append("--no-export-figures")

    logger.info(f"Running Phase 5 Cleanlab Vision XAI: {' '.join(cmd_xai)}")
    subprocess.run(cmd_xai, check=True)
    verify_file_exists(audit_csv, "Cleanlab label noise audit CSV")

    # -------------------------------------------------------------------------
    # Stage 6: Multimodal Spatial Random Forests & Warren's Niche Identity Tests
    # -------------------------------------------------------------------------
    spatial_rf_script = PROJECT_ROOT / "scripts" / "analysis" / "06_multimodal_spatial_rf.R"
    verify_file_exists(spatial_rf_script, "Spatial RF R script")

    morph_flags = Path(getattr(args, "morphometrics", None) or cfg["paths"]["morphometric_flags_csv"])
    env_dir = Path(getattr(args, "env_dir", None) or "data/environmental")
    conflict_flags = Path(getattr(args, "conflict_flags", None) or "data/tables/multimodal_conflict_flags.csv")
    spatial_plot = Path(getattr(args, "spatial_plot", None) or "outputs/figures/spatial_rf_niche_importance.pdf")
    spatial_summary = Path(getattr(args, "spatial_summary", None) or "outputs/reports/multimodal_spatial_rf_summary.csv")

    cmd_spatial = [
        rscript_bin, str(spatial_rf_script),
        "--vouchers", str(vouchers_csv),
        "--morphometrics", str(morph_flags),
        "--vision-audit", str(audit_csv),
        "--env-dir", str(env_dir),
        "--output-flags", str(conflict_flags),
        "--output-plot", str(spatial_plot),
        "--output-summary", str(spatial_summary),
    ]
    if getattr(args, "permutations", None) is not None:
        cmd_spatial.extend(["--permutations", str(args.permutations)])
    if getattr(args, "n_trees", None) is not None:
        cmd_spatial.extend(["--n-trees", str(args.n_trees)])

    logger.info(f"Running Phase 6 Multimodal Spatial RF: {' '.join(cmd_spatial)}")
    subprocess.run(cmd_spatial, check=True)
    verify_file_exists(conflict_flags, "multimodal conflict flags CSV")

    # -------------------------------------------------------------------------
    # Stage 7: Multi-Evidence Taxonomic Decision Matrix & Triage Dashboard
    # -------------------------------------------------------------------------
    triage_script = PROJECT_ROOT / "scripts" / "analysis" / "07_triage_dashboard_synthesis.R"
    verify_file_exists(triage_script, "Triage Dashboard R script")

    gmm_summary = Path(getattr(args, "gmm_summary", None) or cfg["paths"]["gmm_report_csv"])
    output_queue = Path(getattr(args, "output_queue", None) or "data/tables/triage_queue.csv")
    synthesis_plot = Path(getattr(args, "synthesis_plot", None) or "outputs/figures/Figure_Integrative_Packera_dubia_Revision.pdf")
    synthesis_report = Path(getattr(args, "synthesis_report", None) or "outputs/reports/Packera_dubia_Taxonomic_Revision_Summary.md")

    cmd_triage = [
        rscript_bin, str(triage_script),
        "--vouchers", str(vouchers_csv),
        "--morphometrics", str(morph_flags),
        "--vision-audit", str(audit_csv),
        "--multimodal-flags", str(conflict_flags),
        "--gmm-summary", str(gmm_summary),
        "--niche-summary", str(spatial_summary),
        "--output-queue", str(output_queue),
        "--output-plot", str(synthesis_plot),
        "--output-report", str(synthesis_report),
    ]
    logger.info(f"Running Phase 7 Triage Dashboard Synthesis: {' '.join(cmd_triage)}")
    subprocess.run(cmd_triage, check=True)
    verify_file_exists(output_queue, "triage queue CSV")
    verify_file_exists(synthesis_report, "taxonomic revision summary markdown report")

    logger.info("=== Synthesis Workflow (Phases 5, 6, 7) Completed Successfully ===")


def run_all(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Executes the end-to-end production workflow sequentially."""
    logger.info("==================================================================")
    logger.info("Executing End-to-End Production Pipeline")
    logger.info("  (harvest -> segment -> morphometrics -> synthesis)")
    logger.info("==================================================================")
    run_harvest(args, cfg)
    run_segment(args, cfg)
    run_morphometrics(args, cfg)
    run_synthesis(args, cfg)
    logger.info("==================================================================")
    logger.info("End-to-End Pipeline Execution Completed Successfully!")
    logger.info("==================================================================")


def run_calibrate_geometry(args: argparse.Namespace, cfg: PipelineConfig) -> None:
    """Executes empirical geometry threshold calibration from annotated training leaves."""
    logger.info("=== Starting Geometric Threshold Calibration ===")
    calib_script = PROJECT_ROOT / "scripts" / "vision" / "tune_geometry_parameters.py"
    verify_file_exists(calib_script, "Geometry tuning script")

    ann_path = Path(getattr(args, "annotations", None) or (PROJECT_ROOT / "data" / "annotations" / "packera_train_coco.json"))
    verify_file_exists(ann_path, "COCO annotations file")

    output_plot = Path(
        getattr(args, "output_plot", None) or
        (PROJECT_ROOT / cfg["paths"].get("figures_dir", "outputs/figures") / "geometry_parameter_distributions.pdf")
    )
    output_plot.parent.mkdir(parents=True, exist_ok=True)

    config_path = Path(getattr(args, "config_path", None) or getattr(args, "config", None) or (PROJECT_ROOT / "config" / "config.yaml"))

    # Resolve Python interpreter with required dependencies
    py_candidates = [
        PROJECT_ROOT / ".venv" / "bin" / "python",
        PROJECT_ROOT / ".venv_LM2" / "bin" / "python",
        get_lm2_python_executable(),
    ]
    python_bin = sys.executable
    for cand in py_candidates:
        if cand.exists():
            python_bin = str(cand)
            break

    cmd = [
        python_bin, str(calib_script),
        "--annotations", str(ann_path),
        "--output-plot", str(output_plot),
        "--config-path", str(config_path),
    ]
    if getattr(args, "update_config", False):
        cmd.append("--update-config")
    if getattr(args, "min_area", None):
        cmd.extend(["--min-area", str(args.min_area)])
    if getattr(args, "category_id", None) is not None:
        cmd.extend(["--category-id", str(args.category_id)])
    if getattr(args, "category_name", None):
        cmd.extend(["--category-name", str(args.category_name)])

    logger.info(f"Running calibration command: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        logger.error(f"Geometry calibration failed with exit code {exc.returncode}")
        sys.exit(exc.returncode)

    verify_file_exists(output_plot, "Diagnostic distribution PDF")
    logger.info("Geometry calibration completed successfully.")


# =============================================================================
# CLI Parser Setup
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="main.py",
        description="Packera dubia Species Delimitation & Morphometrics Production Pipeline.",
    )
    parser.add_argument(
        "-c", "--config",
        type=Path,
        default=None,
        help="Path to YAML configuration file (default: config/config.yaml).",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Enable verbose debugging output.",
    )

    subparsers = parser.add_subparsers(dest="subcommand", required=True, help="Pipeline phase to execute")

    # Subcommand: check-env
    p_check = subparsers.add_parser(
        "check-env",
        help="Run comprehensive preflight environment & dependency diagnostics",
    )
    p_check.add_argument(
        "--strict",
        action="store_true",
        default=False,
        help="Treat warnings (such as missing CUDA) as critical failures.",
    )

    # Subcommand: download-weights
    p_dl = subparsers.add_parser(
        "download-weights",
        help="Download and verify fine-tuned PointRend model weights from remote repository",
    )
    p_dl.add_argument(
        "--force",
        "--overwrite",
        dest="force",
        action="store_true",
        default=False,
        help="Force re-download even if model weights file exists locally.",
    )

    # Subcommand: harvest
    p_harvest = subparsers.add_parser("harvest", help="Phase 1: Ingest GBIF occurrences & download images")
    p_harvest.add_argument("--taxa", nargs="+", default=None, help="Target taxonomic binomials")
    p_harvest.add_argument("--out-dir", type=Path, default=None, help="Output directory for raw vouchers")
    p_harvest.add_argument("--output-csv", type=Path, default=None, help="Curated vouchers output CSV")
    p_harvest.add_argument("--max-records", type=int, default=None, help="Max records per taxon")
    p_harvest.add_argument("--max-uncertainty", type=float, default=None, help="Max georeference uncertainty (m)")
    p_harvest.add_argument("--min-megapixels", type=float, default=None, help="Min resolution in Megapixels")
    p_harvest.add_argument("--min-file-size-kb", type=float, default=None, help="Min file size in KB")
    p_harvest.add_argument("--concurrency", type=int, default=None, help="Download concurrency")
    p_harvest.add_argument("--download-images", action="store_true", default=False, help="Download specimen images")
    p_harvest.add_argument(
        "--force",
        "--overwrite",
        dest="force",
        action="store_true",
        default=False,
        help="Force re-download and overwrite existing cached voucher images.",
    )

    # Subcommand: segment
    p_segment = subparsers.add_parser("segment", help="Phase 2: Segment leaves & export contours")
    p_segment.add_argument("--vouchers", type=Path, default=None, help="Input curated vouchers CSV")
    p_segment.add_argument("--weights", type=Path, default=None, help="Model weights checkpoint (.pth)")
    p_segment.add_argument("--output-dir", type=Path, default=None, help="Root output directory")
    p_segment.add_argument("--device", type=str, default=None, choices=["cuda", "cpu"], help="Inference device")
    p_segment.add_argument("--min-solidity", type=float, default=None, help="Tier 1 min solidity")
    p_segment.add_argument("--min-ucs", type=float, default=None, help="Tier 1 min UCS score")
    p_segment.add_argument("--score-thresh", type=float, default=None, help="PCD detection score threshold")
    p_segment.add_argument("--limit", type=int, default=None, help="Limit number of vouchers to process")
    p_segment.add_argument(
        "--force",
        "--overwrite",
        dest="force",
        action="store_true",
        default=False,
        help="Force re-segmentation and overwrite existing contour CSVs and masks.",
    )

    # Subcommand: morphometrics
    p_morph = subparsers.add_parser("morphometrics", help="Phase 3: Fourier EFA, GMM, & passive CDA in R")
    p_morph.add_argument("--input", type=Path, default=None, help="Input contours directory")
    p_morph.add_argument("--manifest", type=Path, default=None, help="Extracted leaf manifest CSV")
    p_morph.add_argument("--contours-dir", type=Path, default=None, help="Input contours directory")
    p_morph.add_argument("--vouchers", type=Path, default=None, help="Curated vouchers CSV")
    p_morph.add_argument("--harmonics-out", type=Path, default=None, help="Output harmonics CSV")
    p_morph.add_argument("--flags-out", type=Path, default=None, help="Output morphometric flags CSV")
    p_morph.add_argument("--plot-out", type=Path, default=None, help="Output CDA biplot PDF")
    p_morph.add_argument("--report-out", type=Path, default=None, help="Output BIC summary report")
    p_morph.add_argument("--harmonics", type=int, default=None, help="Number of Fourier harmonics (k)")
    p_morph.add_argument("--num-pcs", type=int, default=None, help="Number of PCA dimensions")
    p_morph.add_argument("--max-k", type=int, default=None, help="Max GMM mixture components")

    # Subcommand: synthesis
    p_synth = subparsers.add_parser("synthesis", help="Phases 5-7: Vision XAI, Spatial RF, & Triage Dashboard")
    p_synth.add_argument("--vouchers", type=Path, default=None, help="Curated vouchers CSV")
    p_synth.add_argument("--rosette-dir", type=Path, default=None, help="Directory with cropped rosette patches")
    p_synth.add_argument("--morphometrics", type=Path, default=None, help="Morphometrics flags CSV")
    p_synth.add_argument("--vision-audit", type=Path, default=None, help="Output vision audit CSV")
    p_synth.add_argument("--cleanlab-threshold", type=float, default=0.85, help="Cleanlab noise cutoff threshold")
    p_synth.add_argument("--export-figures", action="store_true", default=True, help="Export diagnostic figures")
    p_synth.add_argument("--no-export-figures", dest="export_figures", action="store_false", help="Disable figure export")
    p_synth.add_argument("--env-dir", type=Path, default=None, help="Directory containing environmental rasters")
    p_synth.add_argument("--conflict-flags", type=Path, default=None, help="Output multimodal conflict flags CSV")
    p_synth.add_argument("--spatial-plot", type=Path, default=None, help="Output spatial RF PDF plot")
    p_synth.add_argument("--spatial-summary", type=Path, default=None, help="Output Warren's identity summary CSV")
    p_synth.add_argument("--permutations", type=int, default=100, help="Warren's identity test permutations")
    p_synth.add_argument("--n-trees", type=int, default=500, help="Random forest trees")
    p_synth.add_argument("--gmm-summary", type=Path, default=None, help="GMM Bayes factor report CSV")
    p_synth.add_argument("--output-queue", type=Path, default=None, help="Output expert triage queue CSV")
    p_synth.add_argument("--synthesis-plot", type=Path, default=None, help="Output 6-panel synthesis plate PDF")
    p_synth.add_argument("--synthesis-report", type=Path, default=None, help="Output revision summary Markdown")
    p_synth.add_argument("--xai-plot", type=Path, default=None, help="Output Grad-CAM panel PNG")

    # Subcommand: run-all
    p_all = subparsers.add_parser("run-all", help="Execute complete pipeline (harvest -> segment -> morphometrics -> synthesis)")
    # Merge key arguments from previous stages
    p_all.add_argument("--taxa", nargs="+", default=None, help="Target taxonomic binomials")
    p_all.add_argument("--out-dir", type=Path, default=None, help="Output directory for raw vouchers")
    p_all.add_argument("--output-csv", type=Path, default=None, help="Curated vouchers output CSV")
    p_all.add_argument("--max-records", type=int, default=None, help="Max records per taxon")
    p_all.add_argument("--max-uncertainty", type=float, default=None, help="Max georeference uncertainty (m)")
    p_all.add_argument("--min-megapixels", type=float, default=None, help="Min resolution in Megapixels")
    p_all.add_argument("--min-file-size-kb", type=float, default=None, help="Min file size in KB")
    p_all.add_argument("--concurrency", type=int, default=None, help="Download concurrency")
    p_all.add_argument("--download-images", action="store_true", default=False, help="Download specimen images")
    p_all.add_argument(
        "--force",
        "--overwrite",
        dest="force",
        action="store_true",
        default=False,
        help="Force re-harvesting and re-segmentation, overwriting existing artifacts.",
    )
    p_all.add_argument("--vouchers", type=Path, default=None, help="Input curated vouchers CSV")
    p_all.add_argument("--weights", type=Path, default=None, help="Model weights checkpoint (.pth)")
    p_all.add_argument("--output-dir", type=Path, default=None, help="Root output directory")
    p_all.add_argument("--device", type=str, default=None, choices=["cuda", "cpu"], help="Inference device")
    p_all.add_argument("--min-solidity", type=float, default=None, help="Tier 1 min solidity")
    p_all.add_argument("--min-ucs", type=float, default=None, help="Tier 1 min UCS score")
    p_all.add_argument("--score-thresh", type=float, default=None, help="PCD detection score threshold")
    p_all.add_argument("--limit", type=int, default=None, help="Limit number of vouchers to process")
    p_all.add_argument("--input", type=Path, default=None, help="Input contours directory")
    p_all.add_argument("--manifest", type=Path, default=None, help="Extracted leaf manifest CSV")
    p_all.add_argument("--contours-dir", type=Path, default=None, help="Input contours directory")
    p_all.add_argument("--harmonics-out", type=Path, default=None, help="Output harmonics CSV")
    p_all.add_argument("--flags-out", type=Path, default=None, help="Output morphometric flags CSV")
    p_all.add_argument("--plot-out", type=Path, default=None, help="Output CDA biplot PDF")
    p_all.add_argument("--report-out", type=Path, default=None, help="Output BIC summary report")
    p_all.add_argument("--harmonics", type=int, default=None, help="Number of Fourier harmonics (k)")
    p_all.add_argument("--num-pcs", type=int, default=None, help="Number of PCA dimensions")
    p_all.add_argument("--max-k", type=int, default=None, help="Max GMM mixture components")
    p_all.add_argument("--rosette-dir", type=Path, default=None, help="Directory with cropped rosette patches")
    p_all.add_argument("--cleanlab-threshold", type=float, default=0.85, help="Cleanlab noise cutoff threshold")
    p_all.add_argument("--export-figures", action="store_true", default=True, help="Export diagnostic figures")
    p_all.add_argument("--no-export-figures", dest="export_figures", action="store_false", help="Disable figure export")

    # Subcommand: calibrate-geometry
    p_calib = subparsers.add_parser(
        "calibrate-geometry",
        help="Derive empirical geometric thresholds for fold detection and lyrate dissection from SAM 2 COCO annotations",
    )
    p_calib.add_argument(
        "--annotations",
        type=Path,
        default=None,
        help="Path to COCO JSON annotations (default: data/annotations/packera_train_coco.json)",
    )
    p_calib.add_argument(
        "--output-plot",
        type=Path,
        default=None,
        help="Destination path for 4-panel diagnostic distribution PDF (default: outputs/figures/geometry_parameter_distributions.pdf)",
    )
    p_calib.add_argument(
        "--config-path",
        type=Path,
        default=None,
        help="Target configuration YAML path to update (default: config/config.yaml)",
    )
    p_calib.add_argument(
        "--min-area",
        type=float,
        default=150.0,
        help="Minimum contour area threshold in pixels to filter noise (default: 150.0)",
    )
    p_calib.add_argument(
        "--category-id",
        type=int,
        default=None,
        help="Specific COCO category ID to target (default: auto-detect Class 0 / ideal_leaf)",
    )
    p_calib.add_argument(
        "--category-name",
        type=str,
        default=None,
        help="Specific COCO category name to target (e.g. basal_leaf_blade, ideal_leaf)",
    )
    p_calib.add_argument(
        "--update-config",
        action="store_true",
        default=False,
        help="Update thresholds.fold_detection and thresholds.solidity in config.yaml",
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    if args.subcommand == "check-env":
        run_check_env(args)
        return

    cfg = PipelineConfig.from_yaml(args.config)

    dispatch = {
        "check-env": run_check_env,
        "download-weights": run_download_weights,
        "harvest": run_harvest,
        "segment": run_segment,
        "morphometrics": run_morphometrics,
        "synthesis": run_synthesis,
        "run-all": run_all,
        "calibrate-geometry": run_calibrate_geometry,
    }

    handler = dispatch.get(args.subcommand)
    if handler:
        handler(args, cfg)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
