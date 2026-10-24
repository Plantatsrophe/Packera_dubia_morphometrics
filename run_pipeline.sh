#!/usr/bin/env bash
# ==============================================================================
# Script: run_pipeline.sh
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   End-to-end run-scoped orchestrator for the 7-phase Packera dubia pipeline.
#   Every execution writes all generated artifacts (masks, contours, cropped
#   patches, CSV tables, reports, logs, and publication figures) into a
#   self-contained run directory:
#       runs/<RUN_ID>/   (default: runs/run_<YYYYMMDD_HHMMSS>/)
#   while preserving shared read-only caches (data/raw_vouchers/, data/raw_images/,
#   and data/environmental/) in the repository root.
# ==============================================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PROJECT_ROOT}"

ORIG_ARGS=("$@")
CUSTOM_RUN_ID=""
SKIP_HARVEST=false
FORCE_FLAG=""
LIMIT_ARG=()
SAM2_FLAG=()
DEMO_FLAG=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-id)
      CUSTOM_RUN_ID="$2"
      shift 2
      ;;
    --run-id=*)
      CUSTOM_RUN_ID="${1#*=}"
      shift
      ;;
    --skip-harvest)
      SKIP_HARVEST=true
      shift
      ;;
    --force|--overwrite)
      FORCE_FLAG="--force"
      shift
      ;;
    --limit)
      LIMIT_ARG=("--limit" "$2")
      shift 2
      ;;
    --sam2-refine)
      SAM2_FLAG=("--sam2-refine")
      shift
      ;;
    --demo)
      DEMO_FLAG=("--demo")
      shift
      ;;
    *)
      shift
      ;;
  esac
done

if [[ -n "${CUSTOM_RUN_ID}" ]]; then
  RUN_ID="${CUSTOM_RUN_ID}"
else
  RUN_ID="${RUN_ID:-run_$(date +%Y%m%d_%H%M%S)}"
fi

RUN_DIR="runs/${RUN_ID}"
export RUN_ID
export PACKERA_RUN_DIR="${PROJECT_ROOT}/${RUN_DIR}"

# 1. Create standardized run-scoped directory hierarchy
mkdir -p \
  "${RUN_DIR}/data/masks/tier1_intact" \
  "${RUN_DIR}/data/masks/tier2_reflected" \
  "${RUN_DIR}/data/contours" \
  "${RUN_DIR}/data/cropped_patches/basal_leaves_raw" \
  "${RUN_DIR}/data/cropped_patches/capitula" \
  "${RUN_DIR}/data/cropped_patches/rosettes_dense" \
  "${RUN_DIR}/data/tables" \
  "${RUN_DIR}/outputs/figures" \
  "${RUN_DIR}/outputs/reports" \
  "${RUN_DIR}/logs"

# Maintain runs/latest symlink pointing to the most recent RUN_DIR
ln -sfn "${RUN_ID}" "runs/latest"

LOG_FILE="${RUN_DIR}/logs/pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

START_TIME="$(date -Iseconds)"
GIT_SHA="$(git rev-parse HEAD 2>/dev/null || echo 'unknown')"

echo "=============================================================================="
echo "Packera dubia Morphometrics Pipeline — Run-Scoped Execution"
echo "  Run ID          : ${RUN_ID}"
echo "  Run Directory   : ${RUN_DIR}"
echo "  Start Timestamp : ${START_TIME}"
echo "  Git Commit SHA  : ${GIT_SHA}"
echo "=============================================================================="

# Resolve Python interpreters (.venv for primary orchestration/DINOv2, .venv_LM2 for LeafMachine2)
if [[ -x "${PROJECT_ROOT}/.venv/bin/python" ]]; then
  PYTHON_MAIN="${PROJECT_ROOT}/.venv/bin/python"
else
  PYTHON_MAIN="python3"
fi

if [[ -x "${PROJECT_ROOT}/.venv_LM2/bin/python" ]]; then
  PYTHON_LM2="${PROJECT_ROOT}/.venv_LM2/bin/python"
else
  PYTHON_LM2="${PYTHON_MAIN}"
fi

# Initialize run_metadata.json at start of run
"${PYTHON_MAIN}" -c "
import sys
from pathlib import Path
from scripts.core.config import write_run_metadata
write_run_metadata(
    run_dir=Path('${RUN_DIR}').resolve(),
    run_id='${RUN_ID}',
    start_time='${START_TIME}',
    cli_args=sys.argv[1:],
)
" -- "${ORIG_ARGS[@]}"

# ------------------------------------------------------------------------------
# Step 01: Voucher Harvesting & Authority Stratification
# ------------------------------------------------------------------------------
if [[ "${SKIP_HARVEST}" == "false" && ! -s "data/tables/curated_vouchers.csv" && ! -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ]]; then
  echo "[Step 01/07] Running Voucher Harvester..."
  "${PYTHON_MAIN}" scripts/data_prep/01_voucher_harvester.py \
    --out-dir "data/raw_vouchers" \
    --output-csv "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    ${FORCE_FLAG}
fi

# Ensure ${RUN_DIR}/data/tables/curated_vouchers.csv exists prior to Step 02
if [[ ! -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ]]; then
  if [[ -s "data/tables/curated_vouchers.csv" ]]; then
    echo "[Step 02 Preflight] Copying shared data/tables/curated_vouchers.csv -> ${RUN_DIR}/data/tables/curated_vouchers.csv"
    cp -f "data/tables/curated_vouchers.csv" "${RUN_DIR}/data/tables/curated_vouchers.csv.tmp"
    mv -f "${RUN_DIR}/data/tables/curated_vouchers.csv.tmp" "${RUN_DIR}/data/tables/curated_vouchers.csv"
  else
    echo "[ERROR] Missing curated_vouchers.csv in both ${RUN_DIR}/data/tables/ and data/tables/." >&2
    exit 1
  fi
fi

# ------------------------------------------------------------------------------
# Step 02: PointRend Segmentation, Geometric Gatekeeper & 2-Tier Leaf Extraction
# ------------------------------------------------------------------------------
echo "[Step 02/07] Running Segmentation & Organ Extraction..."
"${PYTHON_LM2}" scripts/pipeline/02_segment_and_extract.py \
  --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
  --output-dir "${RUN_DIR}/data" \
  ${FORCE_FLAG} \
  "${LIMIT_ARG[@]}" \
  "${SAM2_FLAG[@]}"

# ------------------------------------------------------------------------------
# Step 03: Elliptic Fourier Analysis (Momocs) & Bilateral Symmetry Decomposition
# ------------------------------------------------------------------------------
echo "[Step 03/07] Running Elliptic Fourier Extractor..."
Rscript scripts/morphometrics/03_fourier_extractor.R \
  --input "${RUN_DIR}/data/contours" \
  --contours-dir "${RUN_DIR}/data/contours" \
  --manifest "${RUN_DIR}/data/tables/extracted_leaf_manifest.csv" \
  --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
  --output "${RUN_DIR}/data/tables/leaf_efa_harmonics.csv" \
  --output-individual "${RUN_DIR}/data/tables/leaf_efa_harmonics_individual.csv" \
  --report-out "${RUN_DIR}/outputs/reports/tier_symmetry_validation.csv" \
  --plot-out "${RUN_DIR}/outputs/figures/tier1_vs_tier2_density_overlay.pdf" \
  --figures-dir "${RUN_DIR}/outputs/figures"

# ------------------------------------------------------------------------------
# Step 04: Gaussian Mixture Modeling (mclust), Passive CDA & Allometry
# ------------------------------------------------------------------------------
echo "[Step 04/07] Running GMM Clustering & Passive CDA..."
Rscript scripts/morphometrics/04_gmm_morphotools.R \
  --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
  --harmonics "${RUN_DIR}/data/tables/leaf_efa_harmonics.csv" \
  --reproductive "${RUN_DIR}/data/tables/voucher_reproductive_metrics.csv" \
  --allometry-report "${RUN_DIR}/outputs/reports/allometry_significance_audit.csv" \
  --allometry-pairwise "${RUN_DIR}/outputs/reports/allometry_pairwise_trajectories.csv" \
  --allometry-plot "${RUN_DIR}/outputs/figures/shape_allometry_regression.pdf" \
  --output-batch-audit "${RUN_DIR}/outputs/reports/institutional_batch_effect_audit.csv" \
  --output-flags "${RUN_DIR}/data/tables/morphometric_flags.csv" \
  --output-plot "${RUN_DIR}/outputs/figures/cda_passive_projection.pdf" \
  --output-report "${RUN_DIR}/outputs/reports/gmm_bayes_factors_summary.csv"

# ------------------------------------------------------------------------------
# Step 05: DINOv2 Deep Vision Feature Extraction & Cleanlab Confident Learning
# ------------------------------------------------------------------------------
echo "[Step 05/07] Running DINOv2 Vision XAI & Cleanlab Audit..."
"${PYTHON_MAIN}" scripts/analysis/05_cleanlab_vision_xai.py \
  --rosette-dir "${RUN_DIR}/data/cropped_patches" \
  --vouchers-csv "${RUN_DIR}/data/tables/curated_vouchers.csv" \
  --output-csv "${RUN_DIR}/data/tables/label_noise_audit.csv" \
  --output-figure "${RUN_DIR}/outputs/figures/GradCAM_audit_panel.png" \
  --export-figures \
  "${DEMO_FLAG[@]}"

# ------------------------------------------------------------------------------
# Step 06: Multimodal Spatial Random Forests, SoilGrids/SSURGO & Hopkins Phenology
# ------------------------------------------------------------------------------
echo "[Step 06/07] Running Multimodal Spatial RF & Niche Identity Tests..."
Rscript scripts/analysis/06_multimodal_spatial_rf.R \
  --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
  --morphometrics "${RUN_DIR}/data/tables/morphometric_flags.csv" \
  --vision-audit "${RUN_DIR}/data/tables/label_noise_audit.csv" \
  --pheno-states "${RUN_DIR}/data/tables/voucher_phenological_states.csv" \
  --env-dir "data/environmental" \
  --ssurgo-table "${RUN_DIR}/data/tables/ssurgo_edaphic_validation.csv" \
  --contingency-table "${RUN_DIR}/outputs/reports/micro_edaphic_outcrop_contingency.csv" \
  --pheno-anomalies "${RUN_DIR}/data/tables/phenological_anomalies.csv" \
  --pheno-summary "${RUN_DIR}/outputs/reports/phenological_anomaly_summary.csv" \
  --pheno-plot "${RUN_DIR}/outputs/figures/phenological_latitudinal_anomaly.pdf" \
  --output-flags "${RUN_DIR}/data/tables/multimodal_conflict_flags.csv" \
  --output-summary "${RUN_DIR}/outputs/reports/multimodal_spatial_rf_summary.csv" \
  --output-plot "${RUN_DIR}/outputs/figures/spatial_rf_niche_importance.pdf" \
  "${DEMO_FLAG[@]}"

# ------------------------------------------------------------------------------
# Step 07: Multi-Evidence Taxonomic Decision Matrix & Triage Dashboard Synthesis
# ------------------------------------------------------------------------------
echo "[Step 07/07] Running Triage Dashboard Synthesis..."
Rscript scripts/analysis/07_triage_dashboard_synthesis.R \
  --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
  --morphometrics "${RUN_DIR}/data/tables/morphometric_flags.csv" \
  --vision-audit "${RUN_DIR}/data/tables/label_noise_audit.csv" \
  --multimodal-flags "${RUN_DIR}/data/tables/multimodal_conflict_flags.csv" \
  --gmm-summary "${RUN_DIR}/outputs/reports/gmm_bayes_factors_summary.csv" \
  --niche-summary "${RUN_DIR}/outputs/reports/multimodal_spatial_rf_summary.csv" \
  --output-queue "${RUN_DIR}/data/tables/triage_queue.csv" \
  --output-triage "${RUN_DIR}/outputs/reports/final_curation_triage.csv" \
  --output-figure "${RUN_DIR}/outputs/figures/Figure_Integrative_Packera_dubia_Revision.pdf" \
  --output-report "${RUN_DIR}/outputs/reports/Packera_dubia_Taxonomic_Revision_Summary.md" \
  --output-html "${RUN_DIR}/outputs/reports/Triage_Dashboard_Interactive.html"

# Finalize run_metadata.json with end timestamp and table row counts
END_TIME="$(date -Iseconds)"
"${PYTHON_MAIN}" -c "
import sys
from pathlib import Path
from scripts.core.config import write_run_metadata
write_run_metadata(
    run_dir=Path('${RUN_DIR}').resolve(),
    run_id='${RUN_ID}',
    start_time='${START_TIME}',
    end_time='${END_TIME}',
    cli_args=sys.argv[1:],
)
" -- "${ORIG_ARGS[@]}"

echo "=============================================================================="
echo "Pipeline Run Complete: ${RUN_DIR}"
echo "  Metadata JSON : ${RUN_DIR}/outputs/reports/run_metadata.json"
echo "  Figures Dir   : ${RUN_DIR}/outputs/figures/"
echo "  Tables Dir    : ${RUN_DIR}/data/tables/"
echo "  Latest Link   : runs/latest -> ${RUN_ID}"
echo "=============================================================================="
