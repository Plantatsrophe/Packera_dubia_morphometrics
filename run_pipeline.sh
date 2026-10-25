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
RESEGMENT=false
START_STEP=1
START_STEP_EXPLICIT=false
FOCAL_TAXON=""
EXCLUDE_TAXA_RAW=()
INTRASPECIFIC_MODE="false"
GOLD_ANCHORS=30
GOLD_ONLY=false
DETERMINER_TIER=""
EXCLUDE_MISIDENTIFIED=false
TRIAGE_SOURCE="ALL"
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
      START_STEP_EXPLICIT=true
      shift 2
      ;;
    --start-step=*)
      START_STEP="${1#*=}"
      START_STEP_EXPLICIT=true
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
    --resegment|--force-segment)
      RESEGMENT=true
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
    --exclude-taxon)
      EXCLUDE_TAXA_RAW+=("$2")
      shift 2
      ;;
    --exclude-taxon=*)
      EXCLUDE_TAXA_RAW+=("${1#*=}")
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
    --gold-only)
      GOLD_ONLY=true
      DETERMINER_TIER="Tier_1_Gold"
      shift
      ;;
    --determiner-tier)
      DETERMINER_TIER="$2"
      shift 2
      ;;
    --determiner-tier=*)
      DETERMINER_TIER="${1#*=}"
      shift
      ;;
    --exclude-misidentified)
      EXCLUDE_MISIDENTIFIED=true
      shift
      ;;
    --exclude-misidentified=*)
      EXCLUDE_MISIDENTIFIED="${1#*=}"
      shift
      ;;
    --triage-source)
      if [[ "${TRIAGE_SOURCE}" == "ALL" || -z "${TRIAGE_SOURCE}" ]]; then
        TRIAGE_SOURCE="$2"
      else
        TRIAGE_SOURCE="${TRIAGE_SOURCE},$2"
      fi
      shift 2
      ;;
    --triage-source=*)
      triage_val="${1#*=}"
      if [[ "${TRIAGE_SOURCE}" == "ALL" || -z "${TRIAGE_SOURCE}" ]]; then
        TRIAGE_SOURCE="${triage_val}"
      else
        TRIAGE_SOURCE="${TRIAGE_SOURCE},${triage_val}"
      fi
      shift
      ;;
    *)
      shift
      ;;
  esac
done

EXCLUDE_TAXON_CSV="$(IFS=,; echo "${EXCLUDE_TAXA_RAW[*]:-}")"

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

# Determine whether to reuse canonical root data/ store for Steps 01-02
REUSE_CANONICAL_DATA=false
if [[ -s "${PROJECT_ROOT}/data/tables/curated_vouchers.csv" && -s "${PROJECT_ROOT}/data/tables/extracted_leaf_manifest.csv" ]]; then
  if [[ "${RESEGMENT}" == "false" && ${#LIMIT_ARG[@]} -eq 0 ]]; then
    if [[ -z "${FORCE_FLAG}" ]] || [[ "${START_STEP}" -ge 3 ]]; then
      REUSE_CANONICAL_DATA=true
    fi
  fi
fi

if [[ "${REUSE_CANONICAL_DATA}" == "true" ]]; then
  echo "[Preflight] Seeding Step 01–02 tables, masks, and contours from canonical root data/ store..."
  for tbl in curated_vouchers.csv extracted_leaf_manifest.csv extracted_leaves_manifest.csv failed_qc_vouchers.csv voucher_reproductive_metrics.csv voucher_phenological_states.csv; do
    if [[ -s "${PROJECT_ROOT}/data/tables/${tbl}" ]]; then
      cp -f "${PROJECT_ROOT}/data/tables/${tbl}" "${RUN_DIR}/data/tables/${tbl}.tmp"
      mv -f "${RUN_DIR}/data/tables/${tbl}.tmp" "${RUN_DIR}/data/tables/${tbl}"
    fi
  done
  rm -rf "${RUN_DIR}/data/masks" "${RUN_DIR}/data/contours"
  ln -sfn "${PROJECT_ROOT}/data/masks" "${RUN_DIR}/data/masks"
  ln -sfn "${PROJECT_ROOT}/data/contours" "${RUN_DIR}/data/contours"
else
  mkdir -p \
    "${RUN_DIR}/data/masks/tier1_intact" \
    "${RUN_DIR}/data/masks/tier2_reflected" \
    "${RUN_DIR}/data/contours" \
    "${RUN_DIR}/data/cropped_patches/basal_leaves_raw" \
    "${RUN_DIR}/data/cropped_patches/capitula" \
    "${RUN_DIR}/data/cropped_patches/rosettes_dense"
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
if [[ "${START_STEP}" -le 1 && "${REUSE_CANONICAL_DATA}" == "false" ]]; then
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
    if [[ -z "${FOCAL_TAXON}" && -z "${EXCLUDE_TAXON_CSV}" && -z "${DETERMINER_TIER}" && "${GOLD_ONLY}" != "true" && "${EXCLUDE_MISIDENTIFIED}" != "true" && -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ]]; then
      mkdir -p "${PROJECT_ROOT}/data/tables"
      cp -f "${RUN_DIR}/data/tables/curated_vouchers.csv" "${PROJECT_ROOT}/data/tables/curated_vouchers.csv.tmp"
      mv -f "${PROJECT_ROOT}/data/tables/curated_vouchers.csv.tmp" "${PROJECT_ROOT}/data/tables/curated_vouchers.csv"
    fi
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

# ------------------------------------------------------------------------------
# Preflight Taxon, Determiner-Tier & Misidentification Subset Filter & Cropped Patch Symlinking
# (--focal-taxon / --exclude-taxon / --gold-only / --determiner-tier / --exclude-misidentified / --triage-source)
# Executed BEFORE running Steps 02-07
# ------------------------------------------------------------------------------
run_preflight_subset_filter() {
  if [[ -n "${FOCAL_TAXON:-}" || -n "${EXCLUDE_TAXON_CSV:-}" || -n "${DETERMINER_TIER:-}" || "${GOLD_ONLY:-false}" == "true" || "${EXCLUDE_MISIDENTIFIED:-false}" == "true" ]]; then
    echo "[Preflight] Filtering ${RUN_DIR}/data/tables/ and linking cropped_patches for focal='${FOCAL_TAXON:-}' exclude='${EXCLUDE_TAXON_CSV:-}' determiner_tier='${DETERMINER_TIER:-}' exclude_misidentified='${EXCLUDE_MISIDENTIFIED:-false}' triage_source='${TRIAGE_SOURCE:-ALL}' (intraspecific_mode=${INTRASPECIFIC_MODE:-false}, gold_anchors=${GOLD_ANCHORS:-30})..."
    "${PYTHON_MAIN}" -c "
import importlib
import os
import re
import shutil
import sys
from pathlib import Path
import pandas as pd

run_dir = Path(sys.argv[1]).resolve()
root_data = Path(sys.argv[2]).resolve()
focal_raw = sys.argv[3].strip()
exclude_csv = sys.argv[4].strip()
intraspecific_mode = sys.argv[5].lower() == 'true'
gold_anchors = int(sys.argv[6])
determiner_tier_raw = sys.argv[7].strip()
exclude_misidentified = sys.argv[8].lower() == 'true' if len(sys.argv) > 8 else False
triage_source_raw = sys.argv[9].strip() if len(sys.argv) > 9 else 'ALL'

mod05 = importlib.import_module('scripts.analysis.05_cleanlab_vision_xai')
standardize_packera_taxon = mod05.standardize_packera_taxon

ALLOWED_TREATMENT_SUBTYPES = {
    'Confirmed_Species',
    'Hybrid_Introgressant',
    'Subspecies_Ecotype',
    'Ecophenotypic_Plasticity',
    'Putative_Cryptic_Species',
}

tier_map = {
    'tier_1_gold': 'Tier_1_Gold', 'gold': 'Tier_1_Gold', 'tier_1': 'Tier_1_Gold', 'tier1': 'Tier_1_Gold', '1': 'Tier_1_Gold',
    'tier_2_silver': 'Tier_2_Silver', 'silver': 'Tier_2_Silver', 'tier_2': 'Tier_2_Silver', 'tier2': 'Tier_2_Silver', '2': 'Tier_2_Silver',
    'tier_3_bronze': 'Tier_3_Bronze', 'bronze': 'Tier_3_Bronze', 'tier_3': 'Tier_3_Bronze', 'tier3': 'Tier_3_Bronze', '3': 'Tier_3_Bronze',
}
target_tier = tier_map.get(determiner_tier_raw.lower(), determiner_tier_raw) if determiner_tier_raw else ''

exclude_taxa = set()
if exclude_csv:
    for item in exclude_csv.split(','):
        item = item.strip()
        if item:
            exclude_taxa.add(standardize_packera_taxon(item))

tables_dir = run_dir / 'data' / 'tables'
csv_path = tables_dir / 'curated_vouchers.csv'
df = pd.read_csv(csv_path)
sp_col = 'species_raw' if 'species_raw' in df.columns else ('scientificName' if 'scientificName' in df.columns else 'species')
df['species_standardized'] = df[sp_col].apply(standardize_packera_taxon)

misidentified_cats = set()
if exclude_misidentified:
    runs_candidates = [run_dir.parent, root_data.parent / 'runs', Path.cwd() / 'runs']
    runs_base = next(
        (rc for rc in runs_candidates if rc.is_dir() and any(rc.glob('*/outputs/reports/final_curation_triage.csv'))),
        run_dir.parent,
    )
    triage_files = []
    seen_resolved = set()
    self_triage = (run_dir / 'outputs' / 'reports' / 'final_curation_triage.csv').resolve()

    tokens = [t.strip() for t in (triage_source_raw or 'ALL').split(',') if t.strip()]
    if not tokens:
        tokens = ['ALL']

    for tok in tokens:
        if tok.upper() in ('ALL', '*'):
            for p in sorted(runs_base.glob('*/outputs/reports/final_curation_triage.csv')):
                rp = p.resolve()
                if rp != self_triage and rp not in seen_resolved and p.is_file() and p.stat().st_size > 0:
                    seen_resolved.add(rp)
                    triage_files.append(p)
        else:
            cand_paths = [
                Path(tok),
                Path(tok) / 'outputs' / 'reports' / 'final_curation_triage.csv',
                runs_base / tok / 'outputs' / 'reports' / 'final_curation_triage.csv',
                root_data.parent / 'runs' / tok / 'outputs' / 'reports' / 'final_curation_triage.csv',
                Path.cwd() / 'runs' / tok / 'outputs' / 'reports' / 'final_curation_triage.csv',
            ]
            matched = next((cp for cp in cand_paths if cp.is_file() and cp.stat().st_size > 0), None)
            if matched is None:
                print(
                    f'[ERROR] --triage-source {tok!r} could not be resolved to an existing final_curation_triage.csv file.',
                    file=sys.stderr,
                )
                sys.exit(1)
            rp = matched.resolve()
            if rp not in seen_resolved:
                seen_resolved.add(rp)
                triage_files.append(matched)

    if not triage_files:
        print(
            f'[ERROR] --exclude-misidentified requested (triage_source={triage_source_raw!r}), but no runs/*/outputs/reports/final_curation_triage.csv files were found.',
            file=sys.stderr,
        )
        sys.exit(1)

    for tf in triage_files:
        tdf_tri = pd.read_csv(tf, low_memory=False)
        sub_col = (
            'taxonomic_treatment_subtype'
            if 'taxonomic_treatment_subtype' in tdf_tri.columns
            else ('taxonomic_treatment_sub' if 'taxonomic_treatment_sub' in tdf_tri.columns else None)
        )
        if sub_col is None or 'catalogNumber' not in tdf_tri.columns:
            print(
                f'[ERROR] Triage report {tf} is missing required columns (catalogNumber, taxonomic_treatment_subtype).',
                file=sys.stderr,
            )
            sys.exit(1)
        sub_series = tdf_tri[sub_col].fillna('').astype(str).str.strip()
        cat_series = tdf_tri['catalogNumber'].fillna('').astype(str).str.strip()
        is_misid = (sub_series == 'Misidentification_Reassign') | (
            ~sub_series.isin(ALLOWED_TREATMENT_SUBTYPES) & (sub_series != '')
        )
        for c_val in cat_series[is_misid]:
            if c_val:
                misidentified_cats.add(c_val)

    before_misid_n = len(df)
    df = df[~df['catalogNumber'].astype(str).str.strip().isin(misidentified_cats)].copy()
    print(
        f'[Preflight: Exclude Misidentified] Source={triage_source_raw!r} ({len(triage_files)} report(s)): '
        f'excluded {before_misid_n - len(df)} vouchers with taxonomic_treatment_subtype == \'Misidentification_Reassign\' '
        f'({len(misidentified_cats)} unique flagged catalogNumbers across reports; retaining only {sorted(ALLOWED_TREATMENT_SUBTYPES)}): '
        f'retained {len(df)} / {before_misid_n} vouchers.'
    )

if target_tier and 'determiner_tier' in df.columns:
    before_tier_n = len(df)
    df = df[df['determiner_tier'].astype(str).str.strip() == target_tier].copy()
    print(f'[Preflight: Determiner Tier] Filtered strictly to determiner_tier == {target_tier!r}: retained {len(df)} / {before_tier_n} vouchers.')

if exclude_taxa:
    before_n = len(df)
    df = df[~df['species_standardized'].isin(exclude_taxa)].copy()
    retained_taxa = sorted(df['species_standardized'].unique().tolist())
    print(f'[Preflight: Exclude Taxa] Excluded {sorted(exclude_taxa)}: retained {len(df)} / {before_n} vouchers across {retained_taxa}.')

if focal_raw:
    focal_taxon = standardize_packera_taxon(focal_raw)
    focal_df = df[df['species_standardized'] == focal_taxon].copy()
    if intraspecific_mode:
        subset_df = focal_df
        print(f'[Strategy 2: Pure Intraspecific] Retained {len(subset_df)} {focal_taxon} vouchers (0 outgroup anchors).')
    else:
        outgroup_dfs = []
        outgroup_counts = {}
        outgroup_taxa = [t for t in mod05.TARGET_TAXA if t != focal_taxon and t not in exclude_taxa]
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
else:
    subset_df = df

if target_tier and 'determiner_tier' in subset_df.columns:
    subset_df = subset_df[subset_df['determiner_tier'].astype(str).str.strip() == target_tier].copy()
if exclude_misidentified and misidentified_cats:
    subset_df = subset_df[~subset_df['catalogNumber'].astype(str).str.strip().isin(misidentified_cats)].copy()

def _write_atomic(df_obj, target_path):
    tmp_p = target_path.with_suffix(target_path.suffix + '.tmp')
    with open(tmp_p, 'w', encoding='utf-8', newline='') as f:
        df_obj.to_csv(f, index=False)
        f.flush()
        os.fsync(f.fileno())
    tmp_p.replace(target_path)

_write_atomic(subset_df, csv_path)
keep_cats = set(subset_df['catalogNumber'].astype(str).str.strip())

# Filter all Step 02 tables in ${RUN_DIR}/data/tables/ to retained catalogNumbers
for tbl_name in ('extracted_leaf_manifest.csv', 'extracted_leaves_manifest.csv', 'failed_qc_vouchers.csv', 'voucher_reproductive_metrics.csv', 'voucher_phenological_states.csv'):
    tp = tables_dir / tbl_name
    if tp.exists() and tp.stat().st_size > 0:
        tdf = pd.read_csv(tp)
        if 'catalogNumber' in tdf.columns:
            tdf_sub = tdf[tdf['catalogNumber'].astype(str).str.strip().isin(keep_cats)].copy()
            _write_atomic(tdf_sub, tp)
            print(f'  Filtered {tbl_name}: {len(tdf_sub)} / {len(tdf)} rows')

# Symlink matching .jpg patches from canonical root data/cropped_patches/
run_patches_root = run_dir / 'data' / 'cropped_patches'
if run_patches_root.is_symlink():
    run_patches_root.unlink()
run_patches_root.mkdir(parents=True, exist_ok=True)

for sub in ('rosettes_dense', 'basal_leaves_raw', 'capitula'):
    src_sub = root_data / 'cropped_patches' / sub
    dst_sub = run_patches_root / sub
    if dst_sub.is_symlink():
        dst_sub.unlink()
    elif dst_sub.exists():
        shutil.rmtree(dst_sub)
    dst_sub.mkdir(parents=True, exist_ok=True)
    if not src_sub.exists():
        continue
    linked = 0
    with os.scandir(src_sub) as it:
        for entry in it:
            if entry.name.endswith('.jpg'):
                stem_cat = re.sub(r'_(?:p\d+_rosette|rosette|leaf\d+|cap\d+)\.jpg$', '', entry.name)
                cat = stem_cat if stem_cat in keep_cats else entry.name.split('_')[0]
                if cat in keep_cats:
                    os.symlink(entry.path, dst_sub / entry.name)
                    linked += 1
    print(f'  Linked {linked} patches in cropped_patches/{sub}')
" "${RUN_DIR}" "${PROJECT_ROOT}/data" "${FOCAL_TAXON:-}" "${EXCLUDE_TAXON_CSV:-}" "${INTRASPECIFIC_MODE:-false}" "${GOLD_ANCHORS:-30}" "${DETERMINER_TIER:-}" "${EXCLUDE_MISIDENTIFIED:-false}" "${TRIAGE_SOURCE:-ALL}"
  elif [[ "${REUSE_CANONICAL_DATA}" == "true" ]]; then
    mkdir -p "${RUN_DIR}/data/cropped_patches"
    for sub in rosettes_dense basal_leaves_raw capitula; do
      rm -rf "${RUN_DIR}/data/cropped_patches/${sub}"
      ln -sfn "${PROJECT_ROOT}/data/cropped_patches/${sub}" "${RUN_DIR}/data/cropped_patches/${sub}"
    done
  fi
}

run_preflight_subset_filter

# ------------------------------------------------------------------------------
# Step 02: PointRend Segmentation, Geometric Gatekeeper & 2-Tier Leaf Extraction
# ------------------------------------------------------------------------------
if [[ "${START_STEP}" -le 2 && "${REUSE_CANONICAL_DATA}" == "false" ]]; then
  echo "[Step 02/07] Running Segmentation & Organ Extraction..."
  find "${RUN_DIR}/data/cropped_patches" -type l -delete 2>/dev/null || true
  "${PYTHON_LM2}" scripts/pipeline/02_segment_and_extract.py \
    --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --output-dir "${RUN_DIR}/data" \
    ${FORCE_FLAG} \
    "${LIMIT_ARG[@]}" \
    "${SAM2_FLAG[@]}"

  # Sync newly extracted masks, contours, cropped patches, and Step 02 manifests back to canonical root data/
  echo "[Step 02/07] Syncing newly extracted Step 02 assets back to canonical root data/ store..."
  mkdir -p \
    "${PROJECT_ROOT}/data/masks/tier1_intact" \
    "${PROJECT_ROOT}/data/masks/tier2_reflected" \
    "${PROJECT_ROOT}/data/contours" \
    "${PROJECT_ROOT}/data/cropped_patches/basal_leaves_raw" \
    "${PROJECT_ROOT}/data/cropped_patches/capitula" \
    "${PROJECT_ROOT}/data/cropped_patches/rosettes_dense" \
    "${PROJECT_ROOT}/data/tables"
  cp -alf "${RUN_DIR}/data/masks/." "${PROJECT_ROOT}/data/masks/" 2>/dev/null || cp -rf "${RUN_DIR}/data/masks/." "${PROJECT_ROOT}/data/masks/"
  cp -alf "${RUN_DIR}/data/contours/." "${PROJECT_ROOT}/data/contours/" 2>/dev/null || cp -rf "${RUN_DIR}/data/contours/." "${PROJECT_ROOT}/data/contours/"
  cp -alf "${RUN_DIR}/data/cropped_patches/." "${PROJECT_ROOT}/data/cropped_patches/" 2>/dev/null || cp -rf "${RUN_DIR}/data/cropped_patches/." "${PROJECT_ROOT}/data/cropped_patches/"
  if [[ -z "${FOCAL_TAXON}" && -z "${EXCLUDE_TAXON_CSV}" && -z "${DETERMINER_TIER}" && "${GOLD_ONLY}" != "true" && "${EXCLUDE_MISIDENTIFIED}" != "true" && ${#LIMIT_ARG[@]} -eq 0 ]]; then
    "${PYTHON_MAIN}" -c "
import os
from pathlib import Path
run_tables = Path('${RUN_DIR}/data/tables')
root_tables = Path('${PROJECT_ROOT}/data/tables')
for tbl in ('extracted_leaf_manifest.csv', 'extracted_leaves_manifest.csv', 'failed_qc_vouchers.csv', 'voucher_reproductive_metrics.csv', 'voucher_phenological_states.csv'):
    src = run_tables / tbl
    if not src.exists():
        continue
    dst = root_tables / tbl
    txt = src.read_text(encoding='utf-8')
    if 'manifest' in tbl:
        txt = txt.replace('/${RUN_DIR}/data/', '/data/')
    tmp = dst.with_suffix(dst.suffix + '.tmp')
    with open(tmp, 'w', encoding='utf-8', newline='') as f:
        f.write(txt)
        f.flush()
        os.fsync(f.fileno())
    tmp.replace(dst)
"
  fi
  run_preflight_subset_filter
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
