---
name: pipeline-contract-audit
description: Audits cross-stage Python and R CSV schema contracts, checks for forbidden _archive or synthetic fallback leaks, and runs the test suite. Use before committing code changes or when debugging pipeline stage handshakes.
---

# Pipeline Contract & Schema Audit (`UNC-BOT-SOP-2026-04-REV4`)

Use this skill before committing code changes, after modifying any pipeline stage (`01` through `07`), or when diagnosing missing-column / join-key errors across Python and R scripts.

---

## Step 1: Scan Active `.py` and `.R` Files for Forbidden `_archive` & Synthetic Leaks

Active code must never import from, read from, or fall back to `scripts/_archive/` or `data/_archive/`, and must never silently synthesize random embeddings, Gaussian heatmaps, or trigonometric environmental variables outside an explicit `--demo` flag.

Run the following verification script from the repository root:

```bash
python3 - <<'EOF'
import re
import sys
from pathlib import Path

root = Path(".").resolve()
active_files = [
    p for p in list(root.rglob("*.py")) + list(root.rglob("*.R"))
    if "_archive" not in p.parts
    and ".venv" not in p.parts
    and ".venv_LM2" not in p.parts
    and "LeafMachine2" not in p.parts
]

archive_violations = []
for f in sorted(active_files):
    rel = f.relative_to(root)
    # Skip test files that explicitly assert _archive isolation
    if rel.parts[:2] == ("scripts", "tests"):
        continue
    for lineno, line in enumerate(f.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if "_archive" in line:
            archive_violations.append(f"{rel}:{lineno}: {stripped}")

if archive_violations:
    print("[FAIL] Forbidden _archive references found in active code:")
    for v in archive_violations:
        print("  -", v)
    sys.exit(1)
else:
    print(f"[PASS] Scanned {len(active_files)} active .py/.R files: 0 _archive leaks detected.")
EOF
```

---

## Step 2: Audit Cross-Stage Python & R CSV Schema Contracts

Verify that all six core production stages (`02` through `07`) agree on primary keys, column names, join semantics, and dynamic reporting:

### Canonical Stage-to-Stage Handshake Matrix

| Stage Script | Primary Inputs | Primary Outputs | Required Keys & Shared Columns | Critical Handshake Invariants |
| :--- | :--- | :--- | :--- | :--- |
| **`scripts/pipeline/02_segment_and_extract.py`** | `data/tables/curated_vouchers.csv` | `data/tables/extracted_leaf_manifest.csv`, `data/tables/voucher_reproductive_metrics.csv`, `data/contours/*.csv` | `catalogNumber`, `plant_individual_id`, `leaf_id`, `assigned_tier`, `pixels_per_mm`, `is_calibrated_mm`, `contour_path`, `mask_path` | Stage 1 YOLOv5 processes both `Class 0` and `Class 1` candidate boxes with high-DPI upper cauline/tiny fragment pre-filtering; post-gatekeeper basal filtering (`cy_rel >= 0.42`, Tier 1 area `>= 12000 * s^2`, Tier 2 area `>= 18000 * s^2` + `cy_rel >= 0.48`, boxy clump rejection, sheet-relative floor `>= 0.28 * max_area`, top-6 cap with disk cleanup); must preserve multi-plant `plant_individual_id` clusters via DBSCAN; `assigned_tier` exports `"Tier_1_Direct"` / `"Tier_2_Reflected"`. |
| **`scripts/morphometrics/03_fourier_extractor.R`** | `data/tables/extracted_leaf_manifest.csv`, `data/contours/*.csv`, `data/tables/curated_vouchers.csv` | `data/tables/fourier_pcs.csv`, `data/tables/leaf_efa_harmonics.csv`, `data/tables/leaf_efa_harmonics_individual.csv` | `(catalogNumber, plant_individual_id)`, `assigned_tier`, `reconstruction_tier`, `pixels_per_mm`, `is_calibrated_mm`, `centroid_size`, `log_centsize`, `PC1`–`PC5`, `Sym_PC1`–`Sym_PC5`, `sym_var`, `asym_var`, `pct_asymmetry`, `harmonic_roughness`, `is_contour_outlier`, `leaf_quality_weight`, `A1`–`D12` | Computes continuous bilateral-symmetric harmonic matrix PCA (`Sym_PC1`–`Sym_PC5` via `Momocs::rm_asym(ef)`) and PERMANOVA across `assigned_tier`; segment-level quality defenses: bilateral asymmetry ratio (`pct_asymmetry`), harmonic roughness (`harmonic_roughness`, $n \ge 6$ vs. $n=1..3$), outlier gating (`is_contour_outlier`: `pct_asymmetry > 0.22 | harmonic_roughness > 0.15`) excluded before fitting PCA axes; exports `leaf_quality_weight`; matches tier case-insensitively via `grepl("^tier_?[12]", assigned_tier, ignore.case = TRUE)`; aggregates specimen medians by `(catalogNumber, plant_individual_id)`. |
| **`scripts/morphometrics/04_gmm_morphotools.R`** | `data/tables/leaf_efa_harmonics.csv` (or `data/tables/fourier_pcs.csv`), `data/tables/voucher_reproductive_metrics.csv` | `data/tables/morphometric_flags.csv` (and `morphometrics_misidentification_flags.csv`), `outputs/reports/gmm_bayes_factors_summary.csv` | `(catalogNumber, plant_individual_id)`, `determiner_tier`, `species_standardized`, `cda_predicted_taxon`, `cda_posterior_prob`, `can1`, `can2`, `gmm_cluster`, `gmm_uncertainty`, `morphometric_outlier_flag`, `loocv_accuracy` | Within-plant quality-weighted (`leaf_quality_weight`) and medoid-filtered aggregation to `(catalogNumber, plant_individual_id)`; robust Mahalanobis (`MASS::cov.rob`) core fitting; parsimonious covariance family restriction (`EII, VII, EEI, VEI, EVI, VVI` when $N/P < 5$) prior to `mclust::Mclust()`; honest Leave-One-Out Cross-Validation accuracy (`loocv_accuracy`) in CDA; exports `morphometric_outlier_flag` (aliased to `misidentification_flag` for downstream consumers). |
| **`scripts/analysis/05_cleanlab_vision_xai.py`** | `data/tables/curated_vouchers.csv`, `data/cropped_patches/` | `data/tables/label_noise_audit.csv` | `catalogNumber`, `plant_individual_id`, `predicted_label` (`vision_predicted_label`), `confidence_predicted_class`, `c_error`, `is_label_corrupted`, `surface_orientation`, `tomentum_taxon_discordance` | L2 cosine normalization of 768-d DINOv2 embeddings prior to fold-internal `StandardScaler`, balanced class weights (`class_weight="balanced"`), diagnostic `tomentum_taxon_discordance` flagging; must neutralize mounting paper (`RGB(128,128,128)`) and control for `surface_orientation` (`"abaxial"` vs. `"adaxial"`). Never synthesize random embeddings unless `--demo` is passed. |
| **`scripts/analysis/06_multimodal_spatial_rf.R`** | `data/tables/curated_vouchers.csv`, `data/tables/morphometric_flags.csv`, `data/tables/label_noise_audit.csv`, `data/tables/voucher_phenological_states.csv` (optional) | `data/tables/multimodal_conflict_flags.csv`, `outputs/reports/multimodal_spatial_rf_summary.csv`, `data/tables/ssurgo_edaphic_validation.csv` | `catalogNumber`, `plant_individual_id`, `misidentification_flag` (resolved from `morphometric_outlier_flag`), `c_error`, `edaphic_best_fit_taxon`, `multimodal_concordance`, `triage_category`, `spatial_block_cv_acc` | Resilient batched `soilDB::SDA_spatialQuery()` chunking with per-batch schema normalization (`dplyr::bind_rows`) and SoilGrids 250m fallback so HTML/HTTP error frames never crash; spatially blocked cross-validation ($0.5^\circ \times 0.5^\circ$ grid blocks, `spatial_block_cv_acc`), centered Hopkins bioclimatic baseline; must use **`dplyr::left_join`** (never `inner_join`) when merging optional phenology tables. Never cache synthetic SSURGO/SoilGrids values to disk. |
| **`scripts/analysis/07_triage_dashboard_synthesis.R`** | `data/tables/curated_vouchers.csv`, `data/tables/morphometric_flags.csv`, `data/tables/label_noise_audit.csv`, `data/tables/multimodal_conflict_flags.csv`, `outputs/reports/gmm_bayes_factors_summary.csv` | `data/tables/triage_queue.csv`, `outputs/figures/Figure_Integrative_Packera_dubia_Revision.pdf`, `outputs/reports/Packera_dubia_Taxonomic_Revision_Summary.md` | `(catalogNumber, plant_individual_id)`, `morphometric_outlier_flag` / `misidentification_flag`, `taxonomic_status_call`, `recommended_determination`, `synthesis_triage_priority` | Enforces `validate_upstream_contracts()` schema verification; must resolve `morphometric_outlier_flag` <-> `misidentification_flag` seamlessly, dynamically compute CDA axis variance labels and GMM $2\Delta\text{BIC}$ values from input tables/summaries (never hardcode static percentages or $2\Delta\text{BIC}$ constants), and write outputs atomically. |

### Automated Contract Verification Script

Run this Python audit to verify column-name agreement and join-key rules across `02`–`07`:

```bash
python3 - <<'EOF'
import sys
from pathlib import Path

scripts_to_check = [
    "scripts/pipeline/02_segment_and_extract.py",
    "scripts/morphometrics/03_fourier_extractor.R",
    "scripts/morphometrics/04_gmm_morphotools.R",
    "scripts/analysis/05_cleanlab_vision_xai.py",
    "scripts/analysis/06_multimodal_spatial_rf.R",
    "scripts/analysis/07_triage_dashboard_synthesis.R",
]

missing = [s for s in scripts_to_check if not Path(s).exists()]
if missing:
    print(f"[FAIL] Missing pipeline scripts: {missing}")
    sys.exit(1)

s02 = Path("scripts/pipeline/02_segment_and_extract.py").read_text()
s03 = Path("scripts/morphometrics/03_fourier_extractor.R").read_text()
s04 = Path("scripts/morphometrics/04_gmm_morphotools.R").read_text()
s05 = Path("scripts/analysis/05_cleanlab_vision_xai.py").read_text()
s06 = Path("scripts/analysis/06_multimodal_spatial_rf.R").read_text()
s07 = Path("scripts/analysis/07_triage_dashboard_synthesis.R").read_text()

checks = [
    ("02 exports assigned_tier & pixels_per_mm", "assigned_tier" in s02 and "pixels_per_mm" in s02),
    ("02 supports Class 0 and Class 1 candidate boxes", "s1_cls not in (0, 1)" in s02 and "LM2_LEAF_CLASS_IDS" in s02),
    ("03 groups by (catalogNumber, plant_individual_id)", "plant_individual_id" in s03 and "catalogNumber" in s03),
    ("03 matches tier via case-insensitive regex", "tier_?1" in s03.lower() or "tier_?[12]" in s03.lower()),
    ("03 exports harmonic_roughness, is_contour_outlier & leaf_quality_weight", "harmonic_roughness" in s03 and "is_contour_outlier" in s03 and "leaf_quality_weight" in s03),
    ("04 aggregates with leaf_quality_weight & cov.rob", "leaf_quality_weight" in s04 and "cov.rob" in s04),
    ("04 exports morphometric_outlier_flag & loocv_accuracy", "morphometric_outlier_flag" in s04 and "loocv_accuracy" in s04),
    ("05 exports c_error, is_label_corrupted & tomentum_taxon_discordance", "c_error" in s05 and "is_label_corrupted" in s05 and "tomentum_taxon_discordance" in s05),
    ("06 uses left_join & spatial_block_cv_acc", "left_join" in s06 and "spatial_block_cv_acc" in s06),
    ("07 validates contracts & ingests CDA/GMM & Vision flags", "validate_upstream_contracts" in s07 and "cda_predicted_taxon" in s07 and "c_error" in s07),
]

failed = False
for label, ok in checks:
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {label}")
    if not ok:
        failed = True

if failed:
    sys.exit(1)
EOF
```

---

## Step 2B: Verify Run-Scoped Directory Layout & Publication Figure Suite Contracts

When executing via `./run_pipeline.sh` or `python main.py --run-id <name>` (or `PACKERA_RUN_DIR`), verify that all generated tables, masks, contours, logs, and figures are isolated inside `${RUN_DIR}` (`runs/<run_id>/`):

* **Run-Scoped Tables (`${RUN_DIR}/data/tables/`):** `curated_vouchers.csv` (snapshotted automatically from `data/tables/curated_vouchers.csv` before Step 02 if not yet present), `extracted_leaf_manifest.csv`, `failed_qc_vouchers.csv`, `voucher_reproductive_metrics.csv`, `voucher_phenological_states.csv`, `fourier_pcs.csv`, `symmetric_fourier_pcs.csv`, `morphometric_flags.csv`, `label_noise_audit.csv`, and `multimodal_conflict_flags.csv`.
* **Run-Scoped Reports & Metadata (`${RUN_DIR}/outputs/reports/`):** `final_curation_triage.csv`, `multimodal_spatial_rf_summary.csv`, and `run_metadata.json` (recording `run_id`, start/end timestamps, `git_commit_sha`, CLI arguments, and table row counts).
* **Publication Figure Suite (`${RUN_DIR}/outputs/figures/`):** Dual vector `.pdf` and 300-DPI `.png` figures with Okabe-Ito palettes (`#0072B2`, `#D55E00`, `#009E73`, `#CC79A7`) and italicized *Packera* taxon labels:
  * Step 02: `fig02a_segmentation_qc_distributions.{pdf,png}`, `fig02b_silhouette_montage.{pdf,png}`
  * Step 03: `fig03a_efa_morphospace_pca.{pdf,png}`, `fig03b_harmonic_power_and_asymmetry.{pdf,png}`
  * Step 04: `fig04a_gmm_bic_and_uncertainty.{pdf,png}`, `fig04b_cda_biplot_and_allometry.{pdf,png}`
  * Step 05: `GradCAM_audit_panel.{png,pdf}`, `fig05b_cleanlab_quality_and_tomentum.{pdf,png}`
  * Step 06: `spatial_rf_niche_importance.{pdf,png}`, `fig06b_geographic_map_and_hopkins_phenology.{pdf,png}`
  * Step 07: `fig07_multimodal_discordance_summary.{pdf,png}`, `Figure_Integrative_Packera_dubia_Revision.{pdf,png}`
* **Shared Read-Only Caches:** Raw voucher imagery (`data/raw_vouchers/` / `data/raw_images/`) and environmental rasters (`data/environmental/`) remain in the repository root and are never duplicated per run.

---

## Step 3: Execute the Central Test Suite (`run_audit.py`)

Always run the centralized QA test suite to verify 100% pass rate across all unit and integration tests:

```bash
python scripts/tests/run_audit.py
```

* Inspect `/tmp/lm2_audit_results.txt` and stdout.
* If any test fails or errors, fix the underlying regression before committing.

