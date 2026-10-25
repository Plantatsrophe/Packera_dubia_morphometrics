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
| **`scripts/pipeline/02_segment_and_extract.py`** | `data/tables/curated_vouchers.csv` | `data/tables/extracted_leaf_manifest.csv`, `data/tables/voucher_reproductive_metrics.csv`, `data/contours/*.csv`, `data/masks/{tier1_intact,tier2_reflected}/*.png`, `data/cropped_patches/{basal_leaves_raw,capitula,rosettes_dense}/` | `catalogNumber`, `plant_individual_id`, `leaf_id`, `assigned_tier`, `pixels_per_mm`, `is_calibrated_mm`, `contour_path`, `mask_path` | Stage 1 YOLOv5 processes both `Class 0` and `Class 1` candidate boxes with high-DPI upper cauline/tiny fragment pre-filtering; resolution-adaptive `w_span = max(min(max(h, w), 350), int(0.65 * max(h, w)))` and resolution-scaled petiole width thresholds; spatial `rosette_centroids_by_pid` wired via `rosette_centroid`; `petiole_attachment_pt=p_base` in all `resample_normalized_contour()` calls; post-gatekeeper basal filtering (`cy_rel >= 0.42`, Tier 1 area `>= 12000 * s^2`, Tier 2 area `>= 18000 * s^2` + `cy_rel >= 0.48`, boxy clump rejection, sheet-relative floor `>= 0.28 * max_area`, top-6 cap with disk cleanup); strengthened `is_cauline_or_stem_contaminant()` rejecting upper-stem sessile/clasping cauline bracts (`aspect_ratio > 3.8` & `solidity < 0.78`, or displaced `> 2.2 * med_rosette_radius` above basal rosette zone) even when `capitula_count == 0`; preserves multi-plant `plant_individual_id` clusters via DBSCAN; `assigned_tier` exports `"Tier_1_Direct"` / `"Tier_2_Reflected"`. |
| **`scripts/morphometrics/03_fourier_extractor.R`** | `data/tables/extracted_leaf_manifest.csv`, `data/contours/*.csv`, `data/tables/curated_vouchers.csv` | `data/tables/fourier_pcs.csv`, `data/tables/leaf_efa_harmonics.csv`, `data/tables/leaf_efa_harmonics_individual.csv` | `(catalogNumber, plant_individual_id)`, `assigned_tier`, `reconstruction_tier`, `pixels_per_mm`, `is_calibrated_mm`, `centroid_size`, `log_centsize`, `PC1`–`PC5`, `Sym_PC1`–`Sym_PC5`, `sym_var`, `asym_var`, `pct_asymmetry`, `harmonic_roughness`, `is_contour_outlier`, `leaf_quality_weight`, `n_tier1_leaves`, `n_tier2_leaves`, `morph_tier_source`, `A1`–`D12` | Merges voucher metadata prior to PCA; fits weighted covariance eigendecomposition (`stats::cov.wt(..., wt = tier_wt * leaf_quality_weight)`) using determiner authority weights (`Tier_1_Gold = 3.0`, `Tier_2_Silver = 1.5`, `Tier_3_Bronze = 0.5`) **strictly on intact `Tier_1_Direct` contours** (`t1_pca_mask`), projecting `Tier_2_Reflected` passively; segment-level quality defenses (`pct_asymmetry`, `harmonic_roughness`, `is_contour_outlier`: `pct_asymmetry > 0.22 | harmonic_roughness > 0.15`) excluded before fitting PCA axes; **Tier-1 Pure Priority aggregation**: whenever `n_tier1_leaves >= 1`, specimen medians in `leaf_efa_harmonics.csv` are computed **100% from `Tier_1_Direct` leaves** (`morph_tier_source = "Tier_1_Pure"`), using `Tier_2_Reflected` only when `n_tier1_leaves == 0` (`morph_tier_source = "Tier_2_Fallback"`). |
| **`scripts/morphometrics/04_gmm_morphotools.R`** | `data/tables/leaf_efa_harmonics.csv` (or `data/tables/fourier_pcs.csv`), `data/tables/voucher_reproductive_metrics.csv` | `data/tables/morphometric_flags.csv` (and `morphometrics_misidentification_flags.csv`), `outputs/reports/gmm_bayes_factors_summary.csv` | `(catalogNumber, plant_individual_id)`, `determiner_tier`, `authority_weight`, `species_standardized`, `morph_tier_source`, `cda_predicted_taxon`, `cda_posterior_prob`, `cda_second_posterior_prob`, `cda_posterior_margin`, `can1`, `can2`, `can3`, `gmm_cluster`, `gmm_uncertainty`, `morphometric_outlier_flag`, `loocv_accuracy` | Within-plant quality-weighted (`leaf_quality_weight`) medoid aggregation to `(catalogNumber, plant_individual_id)`; restricts pooled within-species allometric regression (`stats::lm(..., weights = authority_weight)`), label-blind GMM (`clustvarsel` / `mclust::Mclust` on Gold+Silver `!mahal_outliers` with parsimonious covariance restriction `EII..VVI` when $N/P < 5$), and supervised CDA training anchors strictly to **`Tier_1_Pure` vouchers (`n_tier1_leaves >= 1`, `morph_tier_source == "Tier_1_Pure"`)**, projecting `Tier_2_Fallback` and `Tier_3_Bronze` passively via `predict()`; passive CDA trains on Gold ($3.0$) and Silver ($1.5$) `Tier_1_Pure` anchors with 95% $\chi^2$ Mahalanobis outlier trimming (`qchisq(0.95, df = p)`), uniform priors (`rep(1/g, g)`), `loocv_accuracy`, and 3D canonical axes (`can1`, `can2`, `can3` with orthogonalized `PC3`/`PC2` fallback when $g - 1 < 3$); `audit_institutional_batch_effects` partials out `species_standardized` before testing `institutionCode`; exports `morphometric_outlier_flag` (aliased to `misidentification_flag`). |
| **`scripts/analysis/05_cleanlab_vision_xai.py`** | `data/tables/curated_vouchers.csv`, `data/tables/extracted_leaf_manifest.csv`, `data/cropped_patches/basal_leaves_raw/`, `data/masks/tier1_intact/` | `data/tables/label_noise_audit.csv`, `outputs/figures/GradCAM_audit_panel.{png,pdf}`, `outputs/figures/fig05b_cleanlab_quality_and_tomentum.{pdf,png}` | `catalogNumber`, `plant_individual_id`, `predicted_label` (`vision_predicted_label`), `confidence_predicted_class`, `c_error`, `is_label_corrupted`, `surface_orientation`, `tomentum_taxon_discordance` | Never loads whole-clump `rosettes_dense/` stem crops when `extracted_leaf_manifest.csv` is present; filters strictly to verified `Tier_1_Direct` basal leaf crops (`basal_leaves_raw/{catalogNumber}_leaf{id}.jpg`), dilates each leaf's binary mask (`masks/tier1_intact/{catalogNumber}_leaf{id}.png`) by 5 px (`apply_tier1_leaf_mask_to_rgb`), and replaces exterior stems, cauline bracts, and labels inside the bounding box with neutral herbarium paper RGB (`[242, 239, 230]`); L2 cosine normalization of 768-d DINOv2 embeddings prior to fold-internal `StandardScaler`, authority weights (`3.0 / 1.5 / 0.5`) combined with `class_weight="balanced"` in `LogisticRegression`, dynamically remaps class indices to contiguous `0..K-1` (`active_taxa`) when `--exclude-taxon` is active, and exports `tomentum_taxon_discordance` alongside `is_label_corrupted`. Never synthesize random embeddings unless `--demo` is passed. |
| **`scripts/analysis/06_multimodal_spatial_rf.R`** | `data/tables/curated_vouchers.csv`, `data/tables/morphometric_flags.csv`, `data/tables/label_noise_audit.csv`, `data/environmental/fine_grain_edaphic_cache.csv`, `data/tables/voucher_phenological_states.csv` (optional) | `data/tables/multimodal_conflict_flags.csv`, `outputs/reports/multimodal_spatial_rf_summary.csv`, `data/tables/ssurgo_edaphic_validation.csv` | `catalogNumber`, `plant_individual_id`, `misidentification_flag` (resolved from `morphometric_outlier_flag`), `cda_predicted_taxon`, `cda_posterior_prob`, `cda_second_posterior_prob`, `cda_posterior_margin`, `vision_predicted_label`, `c_error`, `edaphic_best_fit_taxon`, `multimodal_concordance`, `triage_category`, `spatial_block_cv_acc`, `ssurgo_rock_outcrop_pct`, `ssurgo_bedrock_depth_cm`, `fine_soil_ph`, `fine_soil_sand`, `fine_bedrock_depth_cm`, `topo_tpi_30m`, `topo_twi_30m` | Integrates the 10 m – 30 m rock-outcrop and edaphic predictor stack cached persistently in `data/environmental/fine_grain_edaphic_cache.csv` (`scripts/analysis/edaphic_fine_grain_utils.R`): 50 m-buffered USDA-NRCS SSURGO queries (`soilDB::SDA_query`: `ssurgo_rock_outcrop_pct`, `ssurgo_bedrock_depth_cm`, `ssurgo_aws_0_50cm`, component-weighted surface pH, % sand, % clay, OM, CEC), 30 m POLARIS (`polaris_30m/`) & 100 m SOLUS100 (`solus_100m/`: `resdept`, `fragvol`) fused into `fine_soil_ph`, `fine_soil_sand`, `fine_soil_clay`, `fine_bedrock_depth_cm` with 250 m SoilGrids gap-free fallback, and 10 m – 30 m USGS 3DEP micro-topography (`elevatr` + `terra::terrain`: `topo_tpi_30m`, `topo_roughness_30m`, `topo_slope_deg`, `topo_twi_30m`). Uses authority-weighted circular phenology & edaphic envelopes, centered Hopkins' Law baseline (`pred_at_min_lat + 4.0 * (decimalLatitude - min_lat)`), `case.weights = case_wts` (`3.0 / 1.5 / 0.5`) in `ranger::ranger()`, $0.5^\circ \times 0.5^\circ$ `spatial_block_cv_acc`, and **`dplyr::left_join`**. |
| **`scripts/analysis/07_triage_dashboard_synthesis.R`** | `data/tables/curated_vouchers.csv`, `data/tables/morphometric_flags.csv`, `data/tables/label_noise_audit.csv`, `data/tables/multimodal_conflict_flags.csv`, `outputs/reports/gmm_bayes_factors_summary.csv` | `data/tables/triage_queue.csv`, `outputs/reports/final_curation_triage.csv`, `outputs/figures/Figure_Integrative_Packera_dubia_Revision.pdf`, `outputs/reports/Packera_dubia_Taxonomic_Revision_Summary.md` | `(catalogNumber, plant_individual_id)`, `morphometric_outlier_flag` / `misidentification_flag`, `morph_tier_source`, `can1`, `can2`, `can3`, `cda_predicted_taxon`, `cda_posterior_prob`, `cda_second_posterior_prob`, `cda_posterior_margin`, `vision_predicted_label`, `taxonomic_status_call`, `recommended_determination`, `taxonomic_treatment_broad`, `taxonomic_treatment_subtype`, `is_completed_voucher`, `synthesis_triage_priority` | Enforces `validate_upstream_contracts()`; joins `multimodal_conflict_flags.csv` and coalesces non-georeferenced vouchers directly from `morphometric_flags.csv` (Step 04) and `label_noise_audit.csv` (Step 05) via `unique(c("catalogNumber", ...))` achieving 98.4% voucher completion (`5,117 / 5,201` vouchers); evaluates cross-modal known-species reassignments (`both_alt_agree`, `morph_edaph_alt_agree`, `vis_edaph_alt_agree`, `is_severe_misid`) **BEFORE** `Hybrid_Intergrade_Swarm` ($\rightarrow$ `Misidentification_Reassign`); reserves `Putative_Cryptic_Species` (`is_crypt`) strictly for `morph_tier_source == "Tier_1_Pure"` (`n_tier1_leaves >= 1`), `morphometric_outlier_flag == TRUE`, absence of high-confidence known-taxon match, and placement in peripheral GMM clusters (`gmm_cluster > 4`, $|\Delta\text{DOY}| \ge 14.0\text{ d}$) or extreme allochrony ($|\Delta\text{DOY}| \ge 21.0\text{ d}$, `cda_posterior_prob < 0.55`). |

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
    ("02 exports assigned_tier, pixels_per_mm & rosette_centroids_by_pid", "assigned_tier" in s02 and "pixels_per_mm" in s02 and "rosette_centroids_by_pid" in s02),
    ("02 supports Class 0 and Class 1 candidate boxes & is_cauline_or_stem_contaminant", "s1_cls not in (0, 1)" in s02 and "LM2_LEAF_CLASS_IDS" in s02 and "is_cauline_or_stem_contaminant" in s02),
    ("03 groups by (catalogNumber, plant_individual_id), uses cov.wt & exports morph_tier_source", "plant_individual_id" in s03 and "catalogNumber" in s03 and "cov.wt" in s03 and "morph_tier_source" in s03),
    ("03 matches tier via case-insensitive regex", "tier_?1" in s03.lower() or "tier_?[12]" in s03.lower()),
    ("03 exports harmonic_roughness, is_contour_outlier & leaf_quality_weight", "harmonic_roughness" in s03 and "is_contour_outlier" in s03 and "leaf_quality_weight" in s03),
    ("04 aggregates with leaf_quality_weight, authority_weight, Tier_1_Pure & cov.rob", "leaf_quality_weight" in s04 and "authority_weight" in s04 and "Tier_1_Pure" in s04 and "cov.rob" in s04),
    ("04 exports morphometric_outlier_flag, loocv_accuracy, cda_posterior_margin & can3", "morphometric_outlier_flag" in s04 and "loocv_accuracy" in s04 and "cda_posterior_margin" in s04 and "can3" in s04),
    ("05 masks Tier_1_Direct basal leaves & exports c_error, is_label_corrupted & tomentum_taxon_discordance", "apply_tier1_leaf_mask_to_rgb" in s05 and "c_error" in s05 and "is_label_corrupted" in s05 and "tomentum_taxon_discordance" in s05),
    ("06 uses left_join, case.weights, spatial_block_cv_acc & fine_grain_edaphic_cache", "left_join" in s06 and "case.weights" in s06 and "spatial_block_cv_acc" in s06 and "fine_grain_edaphic_cache" in s06),
    ("07 validates contracts, coalesces non-georeferenced vouchers & enforces Tier_1_Pure cryptic gate", "validate_upstream_contracts" in s07 and "taxonomic_treatment_broad" in s07 and "taxonomic_treatment_subtype" in s07 and "is_completed_voucher" in s07 and "Tier_1_Pure" in s07),
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

## Step 2B: Verify Canonical Two-Tier Storage Architecture (`data/` vs. `runs/<RUN_ID>/`) & Figure Suite Contracts

When executing via `./run_pipeline.sh` or `python main.py --run-id <name>` (or `PACKERA_RUN_DIR`), verify adherence to the **Canonical Two-Tier Storage Architecture**:

* **Tier A — Canonical Shared Store (`data/` at repository root):** Houses all invariant Step 01 and Step 02 specimen, environmental, and organ extraction assets:
  * `data/raw_vouchers/` & `data/environmental/` (including `data/environmental/fine_grain_edaphic_cache.csv`, `polaris_30m/`, `solus_100m/`, `dem_3dep_30m/`, `ssurgo_cache/`)
  * `data/masks/tier1_intact/` & `data/masks/tier2_reflected/`
  * `data/contours/`
  * `data/cropped_patches/{rosettes_dense,basal_leaves_raw,capitula}/`
  * The 5 Step 01–02 CSV manifests in `data/tables/`: `curated_vouchers.csv`, `extracted_leaf_manifest.csv`, `failed_qc_vouchers.csv`, `voucher_reproductive_metrics.csv`, and `voucher_phenological_states.csv`.
* **Tier B — Run-Scoped Analysis Directories (`runs/<RUN_ID>/` with `runs/latest` symlink):**
  * Automatically seeds Step 01–02 CSV tables (applying `--exclude-taxon`, `--focal-taxon`, or `--limit` filters when specified) and symlinks `data/masks`, `data/contours`, and `data/cropped_patches/` directly from canonical root `data/` without referencing previous `runs/` directories.
  * Isolates all Step 03–07 statistical tables (`leaf_efa_harmonics.csv`, `fourier_pcs.csv`, `symmetric_fourier_pcs.csv`, `morphometric_flags.csv`, `label_noise_audit.csv`, `multimodal_conflict_flags.csv`), reports (`final_curation_triage.csv`, `multimodal_spatial_rf_summary.csv`, `gmm_bayes_factors_summary.csv`, `Packera_dubia_Taxonomic_Revision_Summary.md`, `run_metadata.json`), logs (`logs/pipeline.log`), and figures (`outputs/figures/`) inside `runs/<RUN_ID>/`.
* **Publication Figure Suite, 3D Ordinations & Universal Unfilled Marker Standard (`${RUN_DIR}/outputs/figures/`):** Dual vector `.pdf` and 300-DPI `.png` figures with Okabe-Ito palettes (`#0072B2`, `#D55E00`, `#009E73`, `#CC79A7`), italicized *Packera* taxon labels, 3D PCA/CDA ordinations (`PC1 x PC2 x PC3`, `Sym_PC1 x Sym_PC2 x Sym_PC3`, `Can1 x Can2 x Can3` via `plot_3d_scatter_box` in `pub_figures_utils.R`), and unfilled (hollow) scatterplot markers (`pch = 1, 2, 5, 0` / `bg = NA` in base R; `shape = 21..25, fill = NA` in `ggplot2`; `facecolors = 'none'` in `matplotlib`):
  * Step 02: `fig02a_segmentation_qc_distributions.{pdf,png}`, `fig02b_silhouette_montage.{pdf,png}`
  * Step 03: `fig03a_efa_morphospace_pca.{pdf,png}`, `fig03b_harmonic_power_and_asymmetry.{pdf,png}`
  * Step 04: `fig04a_gmm_bic_and_uncertainty.{pdf,png}`, `fig04b_cda_biplot_and_allometry.{pdf,png}`
  * Step 05: `GradCAM_audit_panel.{png,pdf}`, `fig05b_cleanlab_quality_and_tomentum.{pdf,png}`
  * Step 06: `spatial_rf_niche_importance.{pdf,png}`, `fig06b_geographic_map_and_hopkins_phenology.{pdf,png}`
  * Step 07: `fig07_multimodal_discordance_summary.{pdf,png}`, `Figure_Integrative_Packera_dubia_Revision.{pdf,png}`

---

## Step 3: Execute the Central Test Suite (`run_audit.py`)

Always run the centralized QA test suite to verify 100% pass rate across all unit and integration tests:

```bash
python scripts/tests/run_audit.py
```

* Inspect `/tmp/lm2_audit_results.txt` and stdout.
* If any test fails or errors, fix the underlying regression before committing.

