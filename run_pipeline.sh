#!/usr/bin/env bash
# ==============================================================================
# Script: run_pipeline.sh
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   End-to-end run-scoped orchestrator for the 7-phase Packera morphometrics pipeline.
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
CUSTOM_DB_PATH=""
SKIP_HARVEST=false
FORCE_FLAG=""
RESEGMENT=false
START_STEP=1
START_STEP_EXPLICIT=false
FOCAL_TAXON=""
CUSTOM_TAXA_RAW=()
EXCLUDE_TAXA_RAW=()
INTRASPECIFIC_MODE="false"
GOLD_ANCHORS=30
GOLD_ONLY=false
DETERMINER_TIER=""
EXCLUDE_MISIDENTIFIED=false
TRIAGE_SOURCE="ALL"
GEOGRAPHIC_SCOPE=""
INCLUDE_WESTERN=false
COUNTRIES_RAW=()
EXPAND_SYNONYMS_FLAG="--expand-synonyms"
DISCOVER_COMPLEXES=false
AUTO_NEIGHBORS=false
MIN_AFFINITY="0.03"
SPATIAL_THIN_KM="0"
MAX_PER_STATE=""
MAX_PER_REGION=""
NO_SPATIAL_WEIGHTS=false
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
    --db)
      CUSTOM_DB_PATH="$2"
      shift 2
      ;;
    --db=*)
      CUSTOM_DB_PATH="${1#*=}"
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
    --taxon|--taxa)
      shift
      while [[ $# -gt 0 && "$1" != -* ]]; do
        CUSTOM_TAXA_RAW+=("$1")
        shift
      done
      ;;
    --taxon=*|--taxa=*)
      CUSTOM_TAXA_RAW+=("${1#*=}")
      shift
      ;;
    --geographic-scope)
      GEOGRAPHIC_SCOPE="$2"
      shift 2
      ;;
    --geographic-scope=*)
      GEOGRAPHIC_SCOPE="${1#*=}"
      shift
      ;;
    --include-western|--no-exclude-western)
      INCLUDE_WESTERN=true
      shift
      ;;
    --countries)
      shift
      while [[ $# -gt 0 && "$1" != -* ]]; do
        COUNTRIES_RAW+=("$1")
        shift
      done
      ;;
    --countries=*)
      COUNTRIES_RAW+=("${1#*=}")
      shift
      ;;
    --expand-synonyms)
      EXPAND_SYNONYMS_FLAG="--expand-synonyms"
      shift
      ;;
    --no-expand-synonyms)
      EXPAND_SYNONYMS_FLAG="--no-expand-synonyms"
      shift
      ;;
    --discover-complexes)
      DISCOVER_COMPLEXES=true
      shift
      ;;
    --auto-neighbors)
      AUTO_NEIGHBORS=true
      DISCOVER_COMPLEXES=true
      shift
      ;;
    --min-affinity)
      MIN_AFFINITY="$2"
      shift 2
      ;;
    --min-affinity=*)
      MIN_AFFINITY="${1#*=}"
      shift
      ;;
    --spatial-thin-km)
      SPATIAL_THIN_KM="$2"
      shift 2
      ;;
    --spatial-thin-km=*)
      SPATIAL_THIN_KM="${1#*=}"
      shift
      ;;
    --max-per-state)
      MAX_PER_STATE="$2"
      shift 2
      ;;
    --max-per-state=*)
      MAX_PER_STATE="${1#*=}"
      shift
      ;;
    --max-per-region)
      MAX_PER_REGION="$2"
      shift 2
      ;;
    --max-per-region=*)
      MAX_PER_REGION="${1#*=}"
      shift
      ;;
    --no-spatial-weights)
      NO_SPATIAL_WEIGHTS=true
      shift
      ;;
    *)
      shift
      ;;
  esac
done

CUSTOM_TAXA_CSV="$(IFS=,; echo "${CUSTOM_TAXA_RAW[*]:-}")"
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
if [[ ${#DEMO_FLAG[@]} -gt 0 ]]; then
  export PACKERA_DEMO_MODE=1
  export PACKERA_DB_PATH="${CUSTOM_DB_PATH:-${PROJECT_ROOT}/${RUN_DIR}/data/demo_ephemeral.sqlite}"
else
  export PACKERA_DB_PATH="${CUSTOM_DB_PATH:-${PROJECT_ROOT}/data/packera_pipeline.sqlite}"
fi

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
echo "Packera Morphometrics Pipeline — Run-Scoped Execution"
echo "  Run ID          : ${RUN_ID}"
echo "  Run Directory   : ${RUN_DIR}"
echo "  SQLite Database : ${PACKERA_DB_PATH}"
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

# Determine whether a non-default geographic scope or custom taxon harvest was requested
CUSTOM_GEO_SCOPE=false
if [[ ( -n "${GEOGRAPHIC_SCOPE}" && "${GEOGRAPHIC_SCOPE}" != "eastern_na" ) || "${INCLUDE_WESTERN}" == "true" || ${#COUNTRIES_RAW[@]} -gt 0 || ${#CUSTOM_TAXA_RAW[@]} -gt 0 ]]; then
  CUSTOM_GEO_SCOPE=true
fi
if [[ "${CUSTOM_GEO_SCOPE}" == "false" && -n "${FOCAL_TAXON}" && "${SKIP_HARVEST}" == "false" && "${START_STEP}" -le 1 && -s "${PROJECT_ROOT}/data/tables/curated_vouchers.csv" ]]; then
  if ! "${PYTHON_MAIN}" -c "
import sys
import pandas as pd
from scripts.core.harvester import standardize_packera_taxon
csv_p, focal_raw = sys.argv[1], sys.argv[2]
target = standardize_packera_taxon(focal_raw)
if not target:
    sys.exit(1)
df = pd.read_csv(csv_p, usecols=lambda c: c in ('scientificName', 'species_raw'))
present = set()
for col in ('scientificName', 'species_raw'):
    if col in df.columns:
        present.update(df[col].dropna().map(standardize_packera_taxon).dropna())
sys.exit(0 if target in present else 1)
" "${PROJECT_ROOT}/data/tables/curated_vouchers.csv" "${FOCAL_TAXON}" 2>/dev/null; then
    CUSTOM_GEO_SCOPE=true
  fi
fi

# Determine whether to reuse canonical root data/ store for Steps 01-02
REUSE_CANONICAL_DATA=false
if [[ -s "${PROJECT_ROOT}/data/tables/curated_vouchers.csv" && -s "${PROJECT_ROOT}/data/tables/extracted_leaf_manifest.csv" ]]; then
  if [[ "${RESEGMENT}" == "false" && ${#LIMIT_ARG[@]} -eq 0 ]]; then
    if [[ "${START_STEP}" -ge 3 ]] || [[ -z "${FORCE_FLAG}" && ( "${CUSTOM_GEO_SCOPE}" == "false" || "${SKIP_HARVEST}" == "true" ) ]]; then
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
  if [[ ${#DEMO_FLAG[@]} -eq 0 ]]; then
    "${PYTHON_MAIN}" -c "
from pathlib import Path
from scripts.core.db import get_db_connection, ingest_existing_csvs_to_db
ws = Path('${PROJECT_ROOT}').resolve()
conn = get_db_connection(init_schema=True, workspace_root=ws)
try:
    v_n = conn.execute('SELECT COUNT(*) FROM vouchers;').fetchone()[0]
    el_n = conn.execute('SELECT COUNT(*) FROM extracted_leaves;').fetchone()[0]
finally:
    conn.close()
if v_n == 0 or el_n == 0:
    counts = ingest_existing_csvs_to_db(ws)
    print(f'[Preflight: SQLite Bootstrap] Synchronized canonical assets into SQLite: {counts}')
"
  fi
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
    if [[ ${#CUSTOM_TAXA_RAW[@]} -gt 0 ]]; then
      HARVEST_TAXON_ARGS=(--taxon "${CUSTOM_TAXA_RAW[@]}" "${EXPAND_SYNONYMS_FLAG}")
    elif [[ -n "${FOCAL_TAXON}" && "${INTRASPECIFIC_MODE}" == "true" ]]; then
      HARVEST_TAXON_ARGS=(--taxon "${FOCAL_TAXON}" "${EXPAND_SYNONYMS_FLAG}")
    else
      HARVEST_TAXON_ARGS=("${EXPAND_SYNONYMS_FLAG}")
    fi
    HARVEST_GEO_ARGS=()
    if [[ -n "${GEOGRAPHIC_SCOPE}" ]]; then
      HARVEST_GEO_ARGS+=(--geographic-scope "${GEOGRAPHIC_SCOPE}")
    fi
    if [[ "${INCLUDE_WESTERN}" == "true" ]]; then
      HARVEST_GEO_ARGS+=(--include-western)
    fi
    if [[ ${#COUNTRIES_RAW[@]} -gt 0 ]]; then
      HARVEST_GEO_ARGS+=(--countries "${COUNTRIES_RAW[@]}")
    fi
    HARVEST_SPATIAL_ARGS=()
    if [[ -n "${SPATIAL_THIN_KM}" && "${SPATIAL_THIN_KM}" != "0" && "${SPATIAL_THIN_KM}" != "0.0" ]]; then
      HARVEST_SPATIAL_ARGS+=(--spatial-thin-km "${SPATIAL_THIN_KM}")
    fi
    if [[ -n "${MAX_PER_STATE}" ]]; then
      HARVEST_SPATIAL_ARGS+=(--max-per-state "${MAX_PER_STATE}")
    fi
    if [[ -n "${MAX_PER_REGION}" ]]; then
      HARVEST_SPATIAL_ARGS+=(--max-per-region "${MAX_PER_REGION}")
    fi
    "${PYTHON_MAIN}" scripts/data_prep/01_voucher_harvester.py \
      --out-dir "data/raw_vouchers" \
      --output-csv "${RUN_DIR}/data/tables/curated_vouchers.csv" \
      "${HARVEST_TAXON_ARGS[@]}" \
      "${HARVEST_GEO_ARGS[@]}" \
      "${HARVEST_SPATIAL_ARGS[@]}" \
      --download-images \
      ${FORCE_FLAG}
    if [[ "${CUSTOM_GEO_SCOPE}" != "true" && -z "${FOCAL_TAXON}" && -z "${EXCLUDE_TAXON_CSV}" && -z "${DETERMINER_TIER}" && "${GOLD_ONLY}" != "true" && "${EXCLUDE_MISIDENTIFIED}" != "true" && ( -z "${SPATIAL_THIN_KM}" || "${SPATIAL_THIN_KM}" == "0" || "${SPATIAL_THIN_KM}" == "0.0" ) && -z "${MAX_PER_STATE}" && -z "${MAX_PER_REGION}" && -s "${RUN_DIR}/data/tables/curated_vouchers.csv" ]]; then
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
  echo "[Preflight] Balancing/filtering ${RUN_DIR}/data/tables/ for taxa='${CUSTOM_TAXA_CSV:-}' focal='${FOCAL_TAXON:-}' exclude='${EXCLUDE_TAXON_CSV:-}' determiner_tier='${DETERMINER_TIER:-}' exclude_misidentified='${EXCLUDE_MISIDENTIFIED:-false}' spatial_thin_km='${SPATIAL_THIN_KM:-0}' max_per_state='${MAX_PER_STATE:-}' max_per_region='${MAX_PER_REGION:-}' no_spatial_weights='${NO_SPATIAL_WEIGHTS:-false}'..."
  "${PYTHON_MAIN}" -c "
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
auto_neighbors_flag = sys.argv[10].lower() == 'true' if len(sys.argv) > 10 else False
min_affinity_val = float(sys.argv[11]) if len(sys.argv) > 11 and sys.argv[11] else 0.03
spatial_thin_km = float(sys.argv[12]) if len(sys.argv) > 12 and sys.argv[12] else 0.0
max_per_state = int(sys.argv[13]) if len(sys.argv) > 13 and sys.argv[13].strip() else None
max_per_region = int(sys.argv[14]) if len(sys.argv) > 14 and sys.argv[14].strip() else None
no_spatial_weights = sys.argv[15].lower() == 'true' if len(sys.argv) > 15 else False
reuse_canonical = sys.argv[16].lower() == 'true' if len(sys.argv) > 16 else False
custom_taxa_csv = sys.argv[17].strip() if len(sys.argv) > 17 else ''

from scripts.core.harvester import (
    select_cascading_outgroup_anchors,
    spatially_balance_and_thin_vouchers,
    standardize_packera_taxon,
)

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

custom_taxa = set()
if custom_taxa_csv:
    for item in custom_taxa_csv.split(','):
        item = item.strip()
        if item:
            std_item = standardize_packera_taxon(item)
            if std_item:
                custom_taxa.add(std_item)

exclude_taxa = set()
if exclude_csv:
    for item in exclude_csv.split(','):
        item = item.strip()
        if item:
            exclude_taxa.add(standardize_packera_taxon(item))

tables_dir = run_dir / 'data' / 'tables'
csv_path = tables_dir / 'curated_vouchers.csv'
df = pd.read_csv(csv_path)
initial_voucher_count = len(df)
if 'scientificName' in df.columns and 'species_raw' in df.columns:
    df['species_standardized'] = [
        standardize_packera_taxon(sn) or standardize_packera_taxon(sr)
        for sn, sr in zip(df['scientificName'], df['species_raw'])
    ]
else:
    sp_col = 'scientificName' if 'scientificName' in df.columns else ('species_raw' if 'species_raw' in df.columns else 'species')
    df['species_standardized'] = df[sp_col].apply(standardize_packera_taxon)

if custom_taxa:
    before_custom_n = len(df)
    df = df[df['species_standardized'].isin(custom_taxa)].copy()
    print(f'[Preflight: Custom Taxa] Restricted strictly to {sorted(custom_taxa)}: retained {len(df)} / {before_custom_n} vouchers.')

from scripts.core.db import (
    ensure_pipeline_run,
    get_db_connection,
    populate_run_active_cohort,
    query_misidentified_vouchers_from_db,
)
target_db_path = os.environ.get('PACKERA_DB_PATH') or str(root_data / 'packera_pipeline.sqlite')

misidentified_cats = set()
if exclude_misidentified:
    tokens = [t.strip() for t in (triage_source_raw or 'ALL').split(',') if t.strip()]
    if not tokens:
        tokens = ['ALL']

    db_run_count = 0
    if os.environ.get('PACKERA_DEMO_MODE') != '1':
        try:
            conn_mis = get_db_connection(db_path=target_db_path, init_schema=True, workspace_root=root_data.parent)
            try:
                db_flagged, db_run_count = query_misidentified_vouchers_from_db(
                    conn_mis,
                    triage_sources=tokens,
                    exclude_run_id=os.environ.get('RUN_ID') or run_dir.name,
                    allowed_subtypes=ALLOWED_TREATMENT_SUBTYPES,
                )
                misidentified_cats.update(db_flagged)
            finally:
                conn_mis.close()
        except Exception:
            db_run_count = 0

    runs_candidates = [run_dir.parent, root_data.parent / 'runs', Path.cwd() / 'runs']
    runs_base = next(
        (rc for rc in runs_candidates if rc.is_dir() and any(rc.glob('*/outputs/reports/final_curation_triage.csv'))),
        run_dir.parent,
    )
    triage_files = []
    seen_resolved = set()
    self_triage = (run_dir / 'outputs' / 'reports' / 'final_curation_triage.csv').resolve()

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
            if matched is None and db_run_count == 0:
                print(
                    f'[ERROR] --triage-source {tok!r} could not be resolved to a SQLite curation_triage run or final_curation_triage.csv file.',
                    file=sys.stderr,
                )
                sys.exit(1)
            if matched is not None:
                rp = matched.resolve()
                if rp not in seen_resolved:
                    seen_resolved.add(rp)
                    triage_files.append(matched)

    if db_run_count == 0 and not triage_files:
        print(
            f'[ERROR] --exclude-misidentified requested (triage_source={triage_source_raw!r}), but no SQLite curation_triage rows or runs/*/outputs/reports/final_curation_triage.csv files were found.',
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
        f'[Preflight: Exclude Misidentified] Source={triage_source_raw!r} (SQLite runs={db_run_count}, CSV reports={len(triage_files)}): '
        f'excluded {before_misid_n - len(df)} vouchers with taxonomic_treatment_subtype == \'Misidentification_Reassign\' '
        f'({len(misidentified_cats)} unique flagged catalogNumbers; retaining only {sorted(ALLOWED_TREATMENT_SUBTYPES)}): '
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
        present_taxa = sorted([t for t in df['species_standardized'].dropna().unique().tolist() if str(t).startswith('Packera ')])
        base_outgroups = present_taxa
        outgroup_taxa = [t for t in base_outgroups if t != focal_taxon and t not in exclude_taxa]

        if auto_neighbors_flag:
            discovered_neighbors = []
            neigh_candidates = [
                run_dir / 'outputs' / 'reports' / 'focal_auto_neighbors.txt',
                run_dir.parent / 'latest' / 'outputs' / 'reports' / 'focal_auto_neighbors.txt',
            ] + sorted(run_dir.parent.glob('*/outputs/reports/focal_auto_neighbors.txt'), reverse=True)
            comp_candidates = [
                run_dir / 'outputs' / 'reports' / 'empirical_taxon_complexes.csv',
                run_dir.parent / 'latest' / 'outputs' / 'reports' / 'empirical_taxon_complexes.csv',
            ] + sorted(run_dir.parent.glob('*/outputs/reports/empirical_taxon_complexes.csv'), reverse=True)
            aff_candidates = [
                run_dir / 'outputs' / 'reports' / 'pairwise_taxon_affinity_matrix.csv',
                run_dir.parent / 'latest' / 'outputs' / 'reports' / 'pairwise_taxon_affinity_matrix.csv',
            ] + sorted(run_dir.parent.glob('*/outputs/reports/pairwise_taxon_affinity_matrix.csv'), reverse=True)

            comp_file = next((p for p in comp_candidates if p.is_file() and p.stat().st_size > 0), None)
            aff_file = next((p for p in aff_candidates if p.is_file() and p.stat().st_size > 0), None)
            if comp_file is not None:
                cdf = pd.read_csv(comp_file)
                if 'taxon' in cdf.columns and 'consensus_complex_id' in cdf.columns:
                    cdf['taxon_std'] = cdf['taxon'].apply(standardize_packera_taxon)
                    focal_rows = cdf[cdf['taxon_std'] == focal_taxon]
                    if not focal_rows.empty:
                        f_cid = focal_rows['consensus_complex_id'].iloc[0]
                        same_c = cdf[(cdf['consensus_complex_id'] == f_cid) & (cdf['taxon_std'] != focal_taxon)]['taxon_std'].tolist()
                        discovered_neighbors.extend(same_c)
            if aff_file is not None:
                adf = pd.read_csv(aff_file)
                if 'taxon' in adf.columns:
                    adf['taxon_std'] = adf['taxon'].apply(standardize_packera_taxon)
                    frow = adf[adf['taxon_std'] == focal_taxon]
                    if not frow.empty:
                        for col in adf.columns:
                            if col in ('taxon', 'taxon_std'):
                                continue
                            col_std = standardize_packera_taxon(col)
                            if col_std != focal_taxon:
                                try:
                                    val = float(frow[col].iloc[0])
                                    if val >= min_affinity_val:
                                        discovered_neighbors.append(col_std)
                                except Exception:
                                    pass
            if not discovered_neighbors:
                nf = next((p for p in neigh_candidates if p.is_file() and p.stat().st_size > 0), None)
                if nf is not None:
                    discovered_neighbors = [
                        standardize_packera_taxon(line.strip())
                        for line in nf.read_text(encoding='utf-8').splitlines()
                        if line.strip() and standardize_packera_taxon(line.strip()) != focal_taxon
                    ]
            if discovered_neighbors:
                disc_set = set(discovered_neighbors) - exclude_taxa - {focal_taxon}
                outgroup_taxa = [t for t in outgroup_taxa if t in disc_set] or sorted(disc_set)
                print(f'[Auto-Neighbors] Empirical complex & affinity neighbors for {focal_taxon}: {outgroup_taxa}')

        outgroup_dfs, outgroup_counts, outgroup_tiers = select_cascading_outgroup_anchors(
            df,
            outgroup_taxa,
            max_anchors=gold_anchors,
            workspace_dir=root_data.parent,
        )
        subset_df = pd.concat([focal_df] + outgroup_dfs, ignore_index=True) if outgroup_dfs else focal_df
        print(
            f'[Strategy 1: Focal + Outgroup Cascading Anchors] Focal {focal_taxon}: {len(focal_df)} | '
            f'Outgroup anchors (Gold->Silver->Bronze): {outgroup_counts} | Tier breakdown: {outgroup_tiers} '
            f'(Total: {len(subset_df)})'
        )
else:
    subset_df = df

if target_tier and 'determiner_tier' in subset_df.columns:
    subset_df = subset_df[subset_df['determiner_tier'].astype(str).str.strip() == target_tier].copy()
if exclude_misidentified and misidentified_cats:
    subset_df = subset_df[~subset_df['catalogNumber'].astype(str).str.strip().isin(misidentified_cats)].copy()

manifest_path = tables_dir / 'extracted_leaf_manifest.csv'
manifest_df = pd.read_csv(manifest_path) if manifest_path.exists() and manifest_path.stat().st_size > 0 else None

pre_balance_n = len(subset_df)
subset_df = spatially_balance_and_thin_vouchers(
    subset_df,
    spatial_thin_km=spatial_thin_km,
    max_per_state=max_per_state,
    max_per_region=max_per_region,
    primary_duplicates_only=True,
    preserve_all_gold=True,
    manifest_df=manifest_df,
    workspace_dir=root_data.parent,
)
if no_spatial_weights:
    subset_df['spatial_balance_weight'] = 1.0

has_spatial_thin_or_quota = (
    spatial_thin_km > 0.0
    or (max_per_state is not None and max_per_state > 0)
    or (max_per_region is not None and max_per_region > 0)
)
if has_spatial_thin_or_quota or len(subset_df) < pre_balance_n:
    print(
        f'[Preflight: Spatial Balancing & Deduplication] Retained {len(subset_df)} / {pre_balance_n} vouchers '
        f'(spatial_thin_km={spatial_thin_km}, max_per_state={max_per_state}, max_per_region={max_per_region}).'
    )
if has_spatial_thin_or_quota and not subset_df.empty:
    if 'stateProvince' in subset_df.columns:
        st_summary = subset_df['stateProvince'].fillna('UNKNOWN').astype(str).value_counts().to_dict()
        print(f'  Retained counts per stateProvince: {st_summary}')
    if 'regional_group' in subset_df.columns:
        reg_summary = subset_df['regional_group'].fillna('NA').astype(str).value_counts().to_dict()
        print(f'  Retained counts per regional_group: {reg_summary}')

def _write_atomic(df_obj, target_path):
    tmp_p = target_path.with_suffix(target_path.suffix + '.tmp')
    with open(tmp_p, 'w', encoding='utf-8', newline='') as f:
        df_obj.to_csv(f, index=False)
        f.flush()
        os.fsync(f.fileno())
    tmp_p.replace(target_path)

_write_atomic(subset_df, csv_path)
keep_cats = set(subset_df['catalogNumber'].astype(str).str.strip())

has_subset_filter = bool(
    custom_taxa_csv or focal_raw or exclude_csv or target_tier or exclude_misidentified or has_spatial_thin_or_quota or len(subset_df) < initial_voucher_count
)

if has_subset_filter:
    # Filter all Step 02 tables in ${RUN_DIR}/data/tables/ to retained catalogNumbers
    for tbl_name in ('extracted_leaf_manifest.csv', 'extracted_leaves_manifest.csv', 'failed_qc_vouchers.csv', 'voucher_reproductive_metrics.csv', 'voucher_phenological_states.csv'):
        tp = tables_dir / tbl_name
        if tp.exists() and tp.stat().st_size > 0:
            tdf = pd.read_csv(tp)
            if 'catalogNumber' in tdf.columns:
                tdf_sub = tdf[tdf['catalogNumber'].astype(str).str.strip().isin(keep_cats)].copy()
                _write_atomic(tdf_sub, tp)
                print(f'  Filtered {tbl_name}: {len(tdf_sub)} / {len(tdf)} rows')

# Register pipeline_runs and populate run_active_cohort in SQLite (skipping canonical DB under --demo)
if os.environ.get('PACKERA_DEMO_MODE') != '1':
    run_id_val = os.environ.get('RUN_ID') or run_dir.name
    conn_coh = get_db_connection(db_path=target_db_path, init_schema=True, workspace_root=root_data.parent)
    try:
        ensure_pipeline_run(
            conn_coh,
            run_id=run_id_val,
            run_dir=run_dir,
            focal_taxon=standardize_packera_taxon(focal_raw) if focal_raw else None,
            intraspecific_mode=intraspecific_mode,
        )
        n_cohort = populate_run_active_cohort(
            conn_coh,
            run_id=run_id_val,
            voucher_rows=subset_df.to_dict(orient='records'),
            clear_existing=True,
        )
        print(f'[Preflight: SQLite Cohort] Registered run_id={run_id_val!r} with {n_cohort} active (catalogNumber, plant_individual_id) rows.')
    finally:
        conn_coh.close()

# Use directory-level symlinks to canonical data/cropped_patches/{rosettes_dense,basal_leaves_raw,capitula}
# instead of creating thousands of per-voucher .jpg symlinks (SQLite run_active_cohort & CSV snapshots filter the cohort)
if reuse_canonical:
    run_patches_root = run_dir / 'data' / 'cropped_patches'
    if run_patches_root.is_symlink():
        run_patches_root.unlink()
    run_patches_root.mkdir(parents=True, exist_ok=True)
    for sub in ('rosettes_dense', 'basal_leaves_raw', 'capitula'):
        src_sub = root_data / 'cropped_patches' / sub
        dst_sub = run_patches_root / sub
        if dst_sub.is_symlink() or dst_sub.is_file():
            dst_sub.unlink()
        elif dst_sub.exists():
            shutil.rmtree(dst_sub)
        if src_sub.exists():
            os.symlink(src_sub, dst_sub)
" "${RUN_DIR}" "${PROJECT_ROOT}/data" "${FOCAL_TAXON:-}" "${EXCLUDE_TAXON_CSV:-}" "${INTRASPECIFIC_MODE:-false}" "${GOLD_ANCHORS:-30}" "${DETERMINER_TIER:-}" "${EXCLUDE_MISIDENTIFIED:-false}" "${TRIAGE_SOURCE:-ALL}" "${AUTO_NEIGHBORS:-false}" "${MIN_AFFINITY:-0.03}" "${SPATIAL_THIN_KM:-0}" "${MAX_PER_STATE:-}" "${MAX_PER_REGION:-}" "${NO_SPATIAL_WEIGHTS:-false}" "${REUSE_CANONICAL_DATA:-false}" "${CUSTOM_TAXA_CSV:-}"
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
  if [[ "${CUSTOM_GEO_SCOPE}" != "true" && -z "${FOCAL_TAXON}" && -z "${EXCLUDE_TAXON_CSV}" && -z "${DETERMINER_TIER}" && "${GOLD_ONLY}" != "true" && "${EXCLUDE_MISIDENTIFIED}" != "true" && ${#LIMIT_ARG[@]} -eq 0 ]]; then
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

  echo "[Step 04/07] Regenerating post-QC Step 02 publication figures..."
  "${PYTHON_MAIN}" scripts/pipeline/02_segment_and_extract.py \
    --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --output-dir "${RUN_DIR}/data" \
    --regenerate-figures-only
fi

# ------------------------------------------------------------------------------
# Step 05: DINOv2 Deep Vision Feature Extraction & Cleanlab Confident Learning
# ------------------------------------------------------------------------------
if [[ "${START_STEP}" -le 5 ]]; then
  echo "[Step 05/07] Running DINOv2 Vision XAI & Cleanlab Audit..."
  "${PYTHON_MAIN}" scripts/analysis/05_cleanlab_vision_xai.py \
    --rosette-dir "${RUN_DIR}/data/cropped_patches" \
    --vouchers-csv "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --manifest-csv "${RUN_DIR}/data/tables/extracted_leaf_manifest.csv" \
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
# Step 08 (Optional): Unsupervised Quantitative Taxon-Complex Discovery
# (--discover-complexes / --auto-neighbors)
# ------------------------------------------------------------------------------
if [[ "${DISCOVER_COMPLEXES}" == "true" || "${AUTO_NEIGHBORS}" == "true" ]]; then
  echo "[Step 08] Running Quantitative Taxon-Complex Discovery Engine..."
  Rscript scripts/analysis/08_discover_taxon_complexes.R \
    --morphometrics "${RUN_DIR}/data/tables/morphometric_flags.csv" \
    --harmonics "${RUN_DIR}/data/tables/leaf_efa_harmonics.csv" \
    --vision-audit "${RUN_DIR}/data/tables/label_noise_audit.csv" \
    --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --multimodal "${RUN_DIR}/data/tables/multimodal_conflict_flags.csv" \
    --min-affinity "${MIN_AFFINITY}" \
    --output-complexes "${RUN_DIR}/outputs/reports/empirical_taxon_complexes.csv" \
    --output-affinity "${RUN_DIR}/outputs/reports/pairwise_taxon_affinity_matrix.csv" \
    --output-neighbors "${RUN_DIR}/outputs/reports/focal_auto_neighbors.txt" \
    --output-figure "${RUN_DIR}/outputs/figures/fig_empirical_taxon_complex_network.pdf" \
    "${FOCAL_ARGS[@]}"
fi

# ------------------------------------------------------------------------------
# Step 07: Multi-Evidence Taxonomic Decision Matrix & Triage Dashboard Synthesis
# ------------------------------------------------------------------------------
if [[ "${START_STEP}" -le 7 ]]; then
  echo "[Step 07/07] Running Triage Dashboard Synthesis..."
  if [[ -n "${FOCAL_TAXON}" ]]; then
    FOCAL_TAXON_SLUG="$(printf '%s' "${FOCAL_TAXON}" | sed -E 's/[^A-Za-z0-9]+/_/g; s/^_+|_+$//g')"
    STEP07_FIGURE="${RUN_DIR}/outputs/figures/Figure_Integrative_${FOCAL_TAXON_SLUG}_Revision.pdf"
    STEP07_REPORT="${RUN_DIR}/outputs/reports/${FOCAL_TAXON_SLUG}_Taxonomic_Revision_Summary.md"
  else
    STEP07_FIGURE="${RUN_DIR}/outputs/figures/Figure_Integrative_Packera_Revision.pdf"
    STEP07_REPORT="${RUN_DIR}/outputs/reports/Packera_Taxonomic_Revision_Summary.md"
  fi
  Rscript scripts/analysis/07_triage_dashboard_synthesis.R \
    --vouchers "${RUN_DIR}/data/tables/curated_vouchers.csv" \
    --morphometrics "${RUN_DIR}/data/tables/morphometric_flags.csv" \
    --vision-audit "${RUN_DIR}/data/tables/label_noise_audit.csv" \
    --multimodal-flags "${RUN_DIR}/data/tables/multimodal_conflict_flags.csv" \
    --gmm-summary "${RUN_DIR}/outputs/reports/gmm_bayes_factors_summary.csv" \
    --niche-summary "${RUN_DIR}/outputs/reports/multimodal_spatial_rf_summary.csv" \
    --output-queue "${RUN_DIR}/data/tables/triage_queue.csv" \
    --output-triage "${RUN_DIR}/outputs/reports/final_curation_triage.csv" \
    --output-figure "${STEP07_FIGURE}" \
    --output-report "${STEP07_REPORT}" \
    --output-html "${RUN_DIR}/outputs/reports/Triage_Dashboard_Interactive.html" \
    "${FOCAL_ARGS[@]}" \
    "${INTRASPECIFIC_ARGS[@]}"
fi

# Finalize run_metadata.json and SQLite pipeline_runs with end timestamp and table row counts
END_TIME="$(date -Iseconds)"
"${PYTHON_MAIN}" -c "
import os
import sys
from pathlib import Path
from scripts.core.config import write_run_metadata
from scripts.core.db import ensure_pipeline_run, get_db_connection
rdir = Path('${RUN_DIR}').resolve()
write_run_metadata(
    run_dir=rdir,
    run_id='${RUN_ID}',
    start_time='${START_TIME}',
    end_time='${END_TIME}',
    cli_args=sys.argv[1:],
)
if os.environ.get('PACKERA_DEMO_MODE') != '1':
    conn = get_db_connection(init_schema=True, workspace_root=Path('${PROJECT_ROOT}').resolve())
    try:
        ensure_pipeline_run(
            conn,
            run_id='${RUN_ID}',
            run_dir=rdir,
            start_timestamp='${START_TIME}',
            end_timestamp='${END_TIME}',
            git_commit_sha='${GIT_SHA}',
            cli_arguments=sys.argv[1:],
        )
    finally:
        conn.close()
" -- "${ORIG_ARGS[@]}"

echo "=============================================================================="
echo "Pipeline Run Complete: ${RUN_DIR}"
echo "  Metadata JSON : ${RUN_DIR}/outputs/reports/run_metadata.json"
echo "  Figures Dir   : ${RUN_DIR}/outputs/figures/"
echo "  Tables Dir    : ${RUN_DIR}/data/tables/"
echo "  Latest Link   : runs/latest -> ${RUN_ID}"
echo "=============================================================================="
