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
START_STEP=1
FOCAL_TAXON=""
INTRASPECIFIC_MODE="false"
GOLD_ANCHORS=30
LIMIT_ARG=()
SAM2_FLAG=()
DEMO_FLAG=()
FOCAL_ARGS=()
INTRASPECIFIC_ARGS=()

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
    --start-step)
      START_STEP="$2"
      shift 2
      ;;
    --start-step=*)
      START_STEP="${1#*=}"
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
    --focal-taxon)
      FOCAL_TAXON="$2"
      shift 2
      ;;
    --focal-taxon=*)
      FOCAL_TAXON="${1#*=}"
      shift
      ;;
    --intraspecific-mode)
      INTRASPECIFIC_MODE="true"
      shift
      ;;
    --gold-anchors)
      GOLD_ANCHORS="$2"
      shift 2
      ;;
    --gold-anchors=*)
      GOLD_ANCHORS="${1#*=}"
      shift
      ;;
    *)
      shift
      ;;
  esac
done

if [[ -n "${FOCAL_TAXON}" ]]; then
  FOCAL_ARGS=("--focal-taxon" "${FOCAL_TAXON}")
fi
if [[ "${INTRASPECIFIC_MODE}" == "true" ]]; then
  INTRASPECIFIC_ARGS=("--intraspecific-mode")
fi

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
if [[ "${START_STEP}" -le 1 ]]; then
  if [[ "${SKIP_HARVEST}" == "false" && ( -n "${FORCE_FLAG}" || ! -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ) ]]; then
    echo "[Step 01/07] Running Voucher Harvester (with high-res image downloads & quality audit)..."
    HARVEST_TAXON_ARGS=()
    if [[ "${INTRASPECIFIC_MODE}" == "true" && -n "${FOCAL_TAXON}" ]]; then
      if [[ "${FOCAL_TAXON}" == "Packera dubia" ]]; then
        HARVEST_TAXON_ARGS=(--taxon "Packera dubia" "Senecio tomentosus" "Packera tomentosa")
      else
        HARVEST_TAXON_ARGS=(--taxon "${FOCAL_TAXON}")
      fi
    fi
    "${PYTHON_MAIN}" scripts/data_prep/01_voucher_harvester.py \
      --out-dir "data/raw_vouchers" \
      --output-csv "${RUN_DIR}/data/tables/curated_vouchers.csv" \
      "${HARVEST_TAXON_ARGS[@]}" \
      --download-images \
      ${FORCE_FLAG}
  fi
fi

# Ensure ${RUN_DIR}/data/tables/curated_vouchers.csv exists prior to Step 02+
if [[ ! -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ]]; then
  if [[ -s "data/tables/curated_vouchers.csv" ]]; then
    echo "[Preflight] Copying shared data/tables/curated_vouchers.csv -> ${RUN_DIR}/data/tables/curated_vouchers.csv"
    cp -f "data/tables/curated_vouchers.csv" "${RUN_DIR}/data/tables/curated_vouchers.csv.tmp"
    mv -f "${RUN_DIR}/data/tables/curated_vouchers.csv.tmp" "${RUN_DIR}/data/tables/curated_vouchers.csv"
  else
    echo "[ERROR] Missing curated_vouchers.csv in both ${RUN_DIR}/data/tables/ and data/tables/." >&2
    exit 1
  fi
fi

# Focal-Taxon Voucher Subset Builder (Strategy 1: Focal + Gold Anchors; Strategy 2: Pure Intraspecific)
if [[ -n "${FOCAL_TAXON}" ]]; then
  echo "[Preflight] Filtering ${RUN_DIR}/data/tables/curated_vouchers.csv for focal taxon '${FOCAL_TAXON}' (intraspecific_mode=${INTRASPECIFIC_MODE}, gold_anchors=${GOLD_ANCHORS})..."
  "${PYTHON_MAIN}" -c "
import importlib
import os
import sys
from pathlib import Path
import pandas as pd

csv_path = Path(sys.argv[1])
focal_raw = sys.argv[2]
intraspecific_mode = sys.argv[3].lower() == 'true'
gold_anchors = int(sys.argv[4])

mod05 = importlib.import_module('scripts.analysis.05_cleanlab_vision_xai')
standardize_packera_taxon = mod05.standardize_packera_taxon
focal_taxon = standardize_packera_taxon(focal_raw)

df = pd.read_csv(csv_path)
sp_col = 'species_raw' if 'species_raw' in df.columns else ('scientificName' if 'scientificName' in df.columns else 'species')
df['species_standardized'] = df[sp_col].apply(standardize_packera_taxon)

focal_df = df[df['species_standardized'] == focal_taxon].copy()
if intraspecific_mode:
    subset_df = focal_df
    print(f'[Strategy 2: Pure Intraspecific] Retained {len(subset_df)} {focal_taxon} vouchers (0 outgroup anchors).')
else:
    outgroup_dfs = []
    outgroup_counts = {}
    outgroup_taxa = [t for t in mod05.TARGET_TAXA if t != focal_taxon]
    for tx in outgroup_taxa:
        tx_df = df[(df['species_standardized'] == tx) & (df.get('determiner_tier', '') == 'Tier_1_Gold')].copy()
        if not tx_df.empty:
            def _has_local(row):
                img_p = str(row.get('image_path', '')).strip()
                cat_p = f\"data/raw_vouchers/{str(row.get('catalogNumber', '')).strip()}.jpg\"
                return (bool(img_p) and Path(img_p).exists()) or Path(cat_p).exists()
            def _has_url_or_path(row):
                for c in ('image_path', 'image_url', 'identifier', 'good_quality_AccessURI'):
                    v = row.get(c, None)
                    if pd.notna(v) and str(v).strip():
                        return True
                return False
            tx_df['_has_local'] = tx_df.apply(_has_local, axis=1)
            tx_df['_has_valid'] = tx_df['_has_local'] | tx_df.apply(_has_url_or_path, axis=1)
            tx_df = tx_df[tx_df['_has_valid']].sort_values(by='_has_local', ascending=False, kind='mergesort')
            tx_sel = tx_df.head(gold_anchors).drop(columns=['_has_local', '_has_valid'])
            outgroup_dfs.append(tx_sel)
            outgroup_counts[tx] = len(tx_sel)
        else:
            outgroup_counts[tx] = 0
    subset_df = pd.concat([focal_df] + outgroup_dfs, ignore_index=True) if outgroup_dfs else focal_df
    print(f'[Strategy 1: Focal + Outgroup Gold Anchors] Focal {focal_taxon}: {len(focal_df)} | Outgroup Tier_1_Gold anchors: {outgroup_counts} (Total: {len(subset_df)})')

tmp_path = csv_path.with_suffix(csv_path.suffix + '.tmp')
with open(tmp_path, 'w', encoding='utf-8', newline='') as f:
    subset_df.to_csv(f, index=False)
    f.flush()
    os.fsync(f.fileno())
tmp_path.replace(csv_path)
" "${RUN_DIR}/data/tables/curated_vouchers.csv" "${FOCAL_TAXON}" "${INTRASPECIFIC_MODE}" "${GOLD_ANCHORS}"
fi

# ------------------------------------------------------------------------------
# Step 02: PointRend Segmentation, Geometric Gatekeeper & 2-Tier Leaf Extraction
# ------------------------------------------------------------------------------
if [[ "${START_STEP}" -le 2 ]]; then
  echo "[Step 02/07] Running Segmentation & Organ Extraction..."
  "${PYTHON_LM2}" scripts/pipeline/02_segment_and_extract.py \
    --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --output-dir "${RUN_DIR}/data" \
    ${FORCE_FLAG} \
    "${LIMIT_ARG[@]}" \
    "${SAM2_FLAG[@]}"
fi

# ------------------------------------------------------------------------------
# Step 03: Elliptic Fourier Analysis (Momocs) & Bilateral Symmetry Decomposition
# ------------------------------------------------------------------------------
if [[ "${START_STEP}" -le 3 ]]; then
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
fi

# ------------------------------------------------------------------------------
# Step 04: Gaussian Mixture Modeling (mclust), Passive CDA & Allometry
# ------------------------------------------------------------------------------
if [[ "${START_STEP}" -le 4 ]]; then
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
    --output-report "${RUN_DIR}/outputs/reports/gmm_bayes_factors_summary.csv" \
    "${FOCAL_ARGS[@]}" \
    "${INTRASPECIFIC_ARGS[@]}"
fi

# ------------------------------------------------------------------------------
# Step 05: DINOv2 Deep Vision Feature Extraction & Cleanlab Confident Learning
# ------------------------------------------------------------------------------
if [[ "${START_STEP}" -le 5 ]]; then
  echo "[Step 05/07] Running DINOv2 Vision XAI & Cleanlab Audit..."
  "${PYTHON_MAIN}" scripts/analysis/05_cleanlab_vision_xai.py \
    --rosette-dir "${RUN_DIR}/data/cropped_patches" \
    --vouchers-csv "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --morphometrics-csv "${RUN_DIR}/data/tables/morphometric_flags.csv" \
    --output-csv "${RUN_DIR}/data/tables/label_noise_audit.csv" \
    --output-figure "${RUN_DIR}/outputs/figures/GradCAM_audit_panel.png" \
    --export-figures \
    "${FOCAL_ARGS[@]}" \
    "${INTRASPECIFIC_ARGS[@]}" \
    "${DEMO_FLAG[@]}"
fi

# ------------------------------------------------------------------------------
# Step 06: Multimodal Spatial Random Forests, SoilGrids/SSURGO & Hopkins Phenology
# ------------------------------------------------------------------------------
if [[ "${START_STEP}" -le 6 ]]; then
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
    "${FOCAL_ARGS[@]}" \
    "${INTRASPECIFIC_ARGS[@]}" \
    "${DEMO_FLAG[@]}"
fi

# ------------------------------------------------------------------------------
# Step 07: Multi-Evidence Taxonomic Decision Matrix & Triage Dashboard Synthesis
# ------------------------------------------------------------------------------
if [[ "${START_STEP}" -le 7 ]]; then
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
    --output-html "${RUN_DIR}/outputs/reports/Triage_Dashboard_Interactive.html" \
    "${FOCAL_ARGS[@]}" \
    "${INTRASPECIFIC_ARGS[@]}"
fi

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
