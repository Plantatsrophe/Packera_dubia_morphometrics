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

from scripts.core.config import PipelineConfig

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
    if args.verbose:
        cmd.append("--verbose")

    logger.info(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    verify_file_exists(output_csv, "harvested vouchers CSV output")
    logger.info("Phase 1 completed successfully.")


def run_segment(args: argparse.Namespace, cfg: PipelineConfig) -> None:
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
    p_all.add_argument("--rosette-dir", type=Path, default=None, help="Directory with cropped rosette patches")
    p_all.add_argument("--cleanlab-threshold", type=float, default=0.85, help="Cleanlab noise cutoff threshold")
    p_all.add_argument("--export-figures", action="store_true", default=True, help="Export diagnostic figures")
    p_all.add_argument("--no-export-figures", dest="export_figures", action="store_false", help="Disable figure export")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    cfg = PipelineConfig.from_yaml(args.config)

    dispatch = {
        "harvest": run_harvest,
        "segment": run_segment,
        "morphometrics": run_morphometrics,
        "synthesis": run_synthesis,
        "run-all": run_all,
    }

    handler = dispatch.get(args.subcommand)
    if handler:
        handler(args, cfg)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
