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
      - run-all: Executes the complete production workflow sequentially.
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import load_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("PackeraPipeline")


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


# =============================================================================
# Phase Execution Handlers
# =============================================================================

def run_harvest(args: argparse.Namespace, cfg: Dict[str, Any]) -> None:
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
    if args.verbose:
        cmd.append("--verbose")

    logger.info(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    verify_file_exists(output_csv, "harvested vouchers CSV output")
    logger.info("Phase 1 completed successfully.")


def run_segment(args: argparse.Namespace, cfg: Dict[str, Any]) -> None:
    """Executes Phase 2: PointRend segmentation, midrib reflection, and contour extraction."""
    logger.info("=== Starting Phase 2: Segmentation & Geometric Leaf Extraction ===")
    
    script_path = PROJECT_ROOT / "scripts" / "pipeline" / "02_segment_and_extract.py"
    verify_file_exists(script_path, "segmentation script")

    vouchers_csv = Path(args.vouchers or cfg["paths"]["curated_vouchers_csv"])
    verify_file_exists(vouchers_csv, "curated vouchers metadata table", "python main.py harvest")

    model_weights = Path(args.weights or cfg["paths"]["model_weights"])
    verify_file_exists(model_weights, "PointRend model weights checkpoint")

    device = args.device or cfg["segmentation"].get("device", "cuda")
    if device == "cuda":
        check_gpu_availability(warn_only=True)

    out_dir = Path(args.output_dir or cfg["paths"]["workspace_root"])
    cmd = [
        sys.executable, str(script_path),
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

    logger.info(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    contours_dir = Path(cfg["paths"]["contours_dir"])
    verify_dir_has_files(contours_dir, "*.csv", "extracted leaf contours")
    logger.info("Phase 2 completed successfully.")


def run_morphometrics(args: argparse.Namespace, cfg: Dict[str, Any]) -> None:
    """Executes Phase 3: EFA, GMM, and passive sample CDA via R scripts."""
    logger.info("=== Starting Phase 3: Morphometrics, EFA & Discriminant Analysis ===")

    rscript_bin = shutil.which("Rscript")
    if not rscript_bin:
        logger.error(
            "Rscript binary not found in system PATH.\n"
            "Phase 3 requires R (>= 4.3) with Momocs, MorphoTools2, and mclust installed."
        )
        sys.exit(1)

    contours_dir = Path(args.contours_dir or cfg["paths"]["contours_dir"])
    verify_dir_has_files(contours_dir, "*.csv", "leaf contours", "python main.py segment")

    vouchers_csv = Path(args.vouchers or cfg["paths"]["curated_vouchers_csv"])
    verify_file_exists(vouchers_csv, "curated vouchers table", "python main.py harvest")

    efa_script = PROJECT_ROOT / "scripts" / "morphometrics" / "03_fourier_extractor.R"
    verify_file_exists(efa_script, "EFA extractor R script")

    efa_out = Path(args.harmonics_out or cfg["paths"]["leaf_efa_harmonics_csv"])
    efa_out.parent.mkdir(parents=True, exist_ok=True)

    harmonics = args.harmonics or cfg["morphometrics"]["harmonics"]
    num_pcs = args.num_pcs or cfg["morphometrics"]["num_pcs"]

    # Step 3A: Fourier EFA
    cmd_efa = [
        rscript_bin, str(efa_script),
        "--contours-dir", str(contours_dir),
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


def run_all(args: argparse.Namespace, cfg: Dict[str, Any]) -> None:
    """Executes the end-to-end production workflow sequentially."""
    logger.info("==================================================================")
    logger.info("Executing End-to-End Production Pipeline (harvest -> segment -> morphometrics)")
    logger.info("==================================================================")
    run_harvest(args, cfg)
    run_segment(args, cfg)
    run_morphometrics(args, cfg)
    logger.info("==================================================================")
    logger.info("End-to-End Pipeline Execution Completed Successfully!")
    logger.info("==================================================================")


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

    # Subcommand: morphometrics
    p_morph = subparsers.add_parser("morphometrics", help="Phase 3: Fourier EFA, GMM, & passive CDA in R")
    p_morph.add_argument("--contours-dir", type=Path, default=None, help="Input contours directory")
    p_morph.add_argument("--vouchers", type=Path, default=None, help="Curated vouchers CSV")
    p_morph.add_argument("--harmonics-out", type=Path, default=None, help="Output harmonics CSV")
    p_morph.add_argument("--flags-out", type=Path, default=None, help="Output morphometric flags CSV")
    p_morph.add_argument("--plot-out", type=Path, default=None, help="Output CDA biplot PDF")
    p_morph.add_argument("--report-out", type=Path, default=None, help="Output BIC summary report")
    p_morph.add_argument("--harmonics", type=int, default=None, help="Number of Fourier harmonics (k)")
    p_morph.add_argument("--num-pcs", type=int, default=None, help="Number of PCA dimensions")
    p_morph.add_argument("--max-k", type=int, default=None, help="Max GMM mixture components")

    # Subcommand: run-all
    p_all = subparsers.add_parser("run-all", help="Execute complete pipeline (harvest -> segment -> morphometrics)")
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
    p_all.add_argument("--vouchers", type=Path, default=None, help="Input curated vouchers CSV")
    p_all.add_argument("--weights", type=Path, default=None, help="Model weights checkpoint (.pth)")
    p_all.add_argument("--output-dir", type=Path, default=None, help="Root output directory")
    p_all.add_argument("--device", type=str, default=None, choices=["cuda", "cpu"], help="Inference device")
    p_all.add_argument("--min-solidity", type=float, default=None, help="Tier 1 min solidity")
    p_all.add_argument("--min-ucs", type=float, default=None, help="Tier 1 min UCS score")
    p_all.add_argument("--score-thresh", type=float, default=None, help="PCD detection score threshold")
    p_all.add_argument("--limit", type=int, default=None, help="Limit number of vouchers to process")
    p_all.add_argument("--contours-dir", type=Path, default=None, help="Input contours directory")
    p_all.add_argument("--harmonics-out", type=Path, default=None, help="Output harmonics CSV")
    p_all.add_argument("--flags-out", type=Path, default=None, help="Output morphometric flags CSV")
    p_all.add_argument("--plot-out", type=Path, default=None, help="Output CDA biplot PDF")
    p_all.add_argument("--report-out", type=Path, default=None, help="Output BIC summary report")
    p_all.add_argument("--harmonics", type=int, default=None, help="Number of Fourier harmonics (k)")
    p_all.add_argument("--num-pcs", type=int, default=None, help="Number of PCA dimensions")
    p_all.add_argument("--max-k", type=int, default=None, help="Max GMM mixture components")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    cfg = load_config(args.config)

    dispatch = {
        "harvest": run_harvest,
        "segment": run_segment,
        "morphometrics": run_morphometrics,
        "run-all": run_all,
    }

    handler = dispatch.get(args.subcommand)
    if handler:
        handler(args, cfg)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
