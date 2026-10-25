# Workspace Rules & Non-Negotiable Invariants (`Packera_dubia_morphometrics`)

**Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV4`  
**Affiliation:** University of North Carolina at Chapel Hill Herbarium (NCU)

This file defines the always-on architectural, mathematical, and cross-language data contracts for the `Packera_dubia_morphometrics` repository. Every AI agent and contributor modifying Python (`.py`) or R (`.R`) code in this workspace must enforce the five non-negotiable invariants below without exception.

---

## 1. Strict Archive Isolation

* **Zero `_archive` Coupling:** Never import from, source, read from, or fall back to `scripts/_archive/` or `data/_archive/` (or `outputs/_archive/`) in any active `.py` or `.R` file.
* **No Legacy Fallbacks:** Decommissioned scripts (such as `scripts/_archive/03_fourier_extractor.py`) and archived tables exist solely for historical provenance. Active production modules in `scripts/core/`, `scripts/data_prep/`, `scripts/pipeline/`, `scripts/vision/`, `scripts/morphometrics/`, `scripts/analysis/`, and `scripts/annotation_and_training/` must fail fast with a clear diagnostic error if an expected active input is missing rather than falling back to `_archive` paths.

---

## 2. Zero Silent Synthetic Data

* **No Unflagged Mock Generation:** Never generate random or mock embeddings (`np.random.randn`), static Gaussian attention/Grad-CAM heatmaps, or deterministic trigonometric environmental/pedological (SoilGrids / WorldClim / USDA SSURGO) values during normal pipeline execution.
* **Explicit `--demo` Requirement:** Any synthetic or simulated data generation path must be strictly gated behind an explicit `--demo` CLI flag. If real rasters, SDA network queries, or neural checkpoints are unavailable and `--demo` is not explicitly passed, raise a clear runtime error or mark records as `NA` / uncalibrated rather than fabricating values silently.
* **Never Cache Synthetic Records to Disk:** Even when running under `--demo` or test harnesses, never persist synthetic/fallback records into production cache directories or canonical tables via `saveRDS()`, `write.csv()`, `readr::write_csv()`, or `pd.DataFrame.to_csv()`.

---

## 3. Dual Virtualenv & Submodule Boundaries

* **Immutable `LeafMachine2/` Submodule:** Never modify files inside the git-ignored `LeafMachine2/` directory directly. All custom botanical adapters, configuration builders, PointRend wrappers, and geometric gatekeepers must reside in `scripts/vision/` or `scripts/pipeline/`.
* **Strict Dual-Environment Separation:**
  * **`.venv` (Primary Orchestration & Modern Vision):** Houses Python orchestration (`main.py`), Segment Anything Model 2 (`SAM 2`), `DINOv2` (`ViT-B/14`), `cleanlab` confident learning, and general scientific utilities.
  * **`.venv_LM2` (LeafMachine2 & Detectron2 PointRend Isolation):** Houses LeafMachine2, YOLOv5 (`LeafPriority.pt` and the 38-class Ruler Detector), and compiled Detectron2 PointRend (`leaf_seg_rcnn_v2_1`) with strict ABI pins: `numpy<2.0` and `opencv-contrib-python<=4.10.0.84`. Never upgrade `numpy>=2.0` or `opencv>4.10.0.84` inside `.venv_LM2`.

---

## 4. Canonical Cross-Language CSV Contracts

All Python (`scripts/pipeline/`, `scripts/vision/`, `scripts/analysis/`) and R (`scripts/morphometrics/`, `scripts/analysis/`) stages must adhere to the exact primary keys and shared column schemas below:

| Contract Element | Canonical Column / Pattern | Cross-Language Rule & Semantics |
| :--- | :--- | :--- |
| **Composite Grouping Key** | `(catalogNumber, plant_individual_id)` | Every leaf-level and specimen-level aggregation across Python and R must group by `(catalogNumber, plant_individual_id)`. **Never** collapse multi-plant herbarium sheets to `plant_individual_id = 0` or group by `catalogNumber` alone when aggregating organ traits. |
| **Routing Tier Column** | `assigned_tier` | Canonical closed-silhouette tiers are `"Tier_1_Direct"` and `"Tier_2_Reflected"` (with shorthand `"tier1"` / `"tier2"` in intermediate routing objects). In R, always match case-insensitively via `grepl("^tier_?[12]", assigned_tier, ignore.case = TRUE)`. |
| **Scale Calibration Columns** | `pixels_per_mm`, `is_calibrated_mm` | `pixels_per_mm` stores the empirical ruler conversion factor (`float` or `NA`), and `is_calibrated_mm` (`bool` / `logical`) explicitly flags whether absolute millimeter traits ($\text{CS}_{\text{mm}}$, involucre height/width in mm) are calibrated vs. uncalibrated pixel units. |
| **Authority-Tier Weighting** | `determiner_tier`, `authority_weight` | Canonical authority weights (`Tier_1_Gold = 3.0`, `Tier_2_Silver = 1.5`, `Tier_3_Bronze = 0.5`) propagate across Steps 03–06: Step 03 (`stats::cov.wt(..., wt = tier_wt * leaf_quality_weight)` after prior voucher metadata merge), Step 04 (`stats::lm(..., weights = authority_weight)`, `Tier_1_Gold` + `Tier_2_Silver` anchor fitting in `clustvarsel`/`mclust::Mclust` with passive `predict()` for `Tier_3_Bronze`, and Gold/Silver CDA with 95% $\chi^2$ Mahalanobis trimming and uniform priors `rep(1/g, g)`), Step 05 (`3.0 / 1.5 / 0.5` + `class_weight="balanced"` in `LogisticRegression`), and Step 06 (`case.weights = case_wts` in `ranger::ranger()` and weighted phenology/edaphic envelopes). |
| **Fourier Quality, Symmetry & Tier-1 Pure Priority** | `harmonic_roughness`, `is_contour_outlier`, `leaf_quality_weight`, `pct_asymmetry`, `sym_var`, `asym_var`, `Sym_PC1`–`Sym_PC5`, `PC1`–`PC5`, `n_tier1_leaves`, `n_tier2_leaves`, `morph_tier_source` | Exported by `03_fourier_extractor.R` (`fourier_pcs.csv`, `leaf_efa_harmonics.csv`). Voucher metadata merge occurs prior to PCA; both `Sym_PC1`–`Sym_PC5` (`Momocs::rm_asym(ef)`) and primary `PC1`–`PC5` use weighted covariance eigendecomposition (`stats::cov.wt(..., wt = tier_wt * leaf_quality_weight)`) fitted **strictly on intact `Tier_1_Direct` contours** (`t1_pca_mask`), projecting `Tier_2_Reflected` contours passively so synthetic mirrored outlines never warp principal axes. Outlier threshold `pct_asymmetry > 0.22 | harmonic_roughness > 0.15` gates PCA fitting (`clean_pca_mask`). **Tier-1 Pure Priority Aggregation:** Whenever a voucher has `n_tier1_leaves >= 1`, specimen-level median harmonics and PCs in `leaf_efa_harmonics.csv` are computed **100% from its `Tier_1_Direct` leaves** (`morph_tier_source = "Tier_1_Pure"`), completely ignoring any `Tier_2_Reflected` leaves on that sheet. `Tier_2_Reflected` contours are used only as a passive fallback when `n_tier1_leaves == 0` (`morph_tier_source = "Tier_2_Fallback"`). |
| **Morphometric Outlier Flag, 3D CDA Posteriors & Tier-1 Pure Anchor Fitting** | `can1`, `can2`, `can3`, `cda_predicted_taxon`, `cda_posterior_prob`, `cda_second_posterior_prob`, `cda_posterior_margin`, `morphometric_outlier_flag`, `loocv_accuracy`, `morph_tier_source` | Exported by `04_gmm_morphotools.R` (`data/tables/morphometric_flags.csv` / `morphometrics_misidentification_flags.csv`) with `morphometric_outlier_flag` aliased to `misidentification_flag` in Steps 06–07. Pooled within-species allometric regression (`test_and_correct_allometry`), label-blind GMM clustering (`mclust::Mclust`), and supervised CDA training anchors (`run_canonical_discriminant_analysis`) are restricted strictly to **`Tier_1_Pure` vouchers (`n_tier1_leaves >= 1`, `morph_tier_source == "Tier_1_Pure"`)**, while `Tier_2_Fallback` vouchers are projected passively via `predict()`. GMM fits on `Tier_1_Gold` + `Tier_2_Silver` `Tier_1_Pure` anchors (`!mahal_outliers` via `MASS::cov.rob`) with parsimonious covariance restriction (`EII`..`VVI` when $N/P < 5$) and projects `Tier_3_Bronze` and `Tier_2_Fallback` via `predict(gmm_fit, newdata = X_pcs)`. Passive CDA trains on Gold ($3.0$) and Silver ($1.5$) `Tier_1_Pure` anchors with per-taxon 95% $\chi^2$ Mahalanobis outlier trimming (`qchisq(0.95, df = p)`), uniform class priors (`rep(1/g, g)`), honest Leave-One-Out Cross-Validation (`loocv_accuracy`), and 3D canonical axes (`can1`, `can2`, `can3`, with orthogonalized `PC3`/`PC2` fallback when $g - 1 < 3$). Institutional herbarium batch-effect auditing (`audit_institutional_batch_effects`) partials out taxonomic composition (`species_standardized`) before testing `institutionCode`. |
| **Masked Basal-Leaf Vision XAI & Tomentum Flag** | `predicted_label` (`vision_predicted_label`), `c_error`, `is_label_corrupted`, `surface_orientation`, `tomentum_taxon_discordance` | Exported by `05_cleanlab_vision_xai.py` (`data/tables/label_noise_audit.csv`). Never loads whole-clump `rosettes_dense/` stem crops when `extracted_leaf_manifest.csv` is present; filters strictly to verified `Tier_1_Direct` basal leaf crops (`basal_leaves_raw/{catalogNumber}_leaf{id}.jpg`), dilates each leaf's binary mask (`masks/tier1_intact/{catalogNumber}_leaf{id}.png`) by 5 px (`apply_tier1_leaf_mask_to_rgb`), and replaces exterior stems, cauline bracts, and labels inside the bounding box with neutral herbarium paper RGB (`[242, 239, 230]`) prior to computing CIE Lab luminance/saturation (`surface_orientation`), L2-normalized 768-d DINOv2 embeddings, authority-weighted (`3.0 / 1.5 / 0.5`) `class_weight="balanced"` `LogisticRegression`, Cleanlab `ConfidentLearning`, and `GradCAM_audit_panel.png`. Dynamically remaps class indices to contiguous `0..K-1` (`active_taxa`) when `--exclude-taxon` is active. |
| **Fine-Grain (10 m – 30 m) Edaphic/Outcrop ENM Stack & Hopkins Baseline** | `spatial_block_cv_acc`, `edaphic_best_fit_taxon`, `multimodal_concordance`, `triage_category`, `ssurgo_rock_outcrop_pct`, `ssurgo_bedrock_depth_cm`, `fine_soil_ph`, `fine_soil_sand`, `fine_bedrock_depth_cm`, `topo_tpi_30m`, `topo_twi_30m` | Exported by `06_multimodal_spatial_rf.R` (`multimodal_conflict_flags.csv`, `multimodal_spatial_rf_summary.csv`) and cached persistently in `data/environmental/fine_grain_edaphic_cache.csv` (`scripts/analysis/edaphic_fine_grain_utils.R`). Integrates (1) 50 m-buffered USDA-NRCS SSURGO queries (`soilDB::SDA_query` with WKT `ST_Buffer`: `ssurgo_rock_outcrop_pct` matching *Rock outcrop*, *Wake*, *Louisburg*, *Ashlar*, *Rion*, *Wateree*, *Wilkes*, *Iredell*, *Mecklenburg*, etc., `ssurgo_bedrock_depth_cm`, `ssurgo_aws_0_50cm`, and component-weighted surface pH, % sand, % clay, OM, and CEC), (2) 30 m POLARIS (`data/environmental/polaris_30m/`) & 100 m SOLUS100 (`data/environmental/solus_100m/`: `resdept`, `fragvol`) fused into `fine_soil_ph`, `fine_soil_sand`, `fine_soil_clay`, and `fine_bedrock_depth_cm` with 250 m SoilGrids gap-free fallback, and (3) 10 m – 30 m USGS 3DEP micro-topography (`elevatr` + `terra::terrain`: `topo_tpi_30m`, `topo_roughness_30m`, `topo_slope_deg`, `topo_twi_30m`). Uses authority-weighted circular phenology & edaphic envelopes, centered Hopkins' Bioclimatic Law baseline (`pred_at_min_lat + 4.0 * (decimalLatitude - min_lat)`), `case.weights = case_wts` (`3.0 / 1.5 / 0.5`) in `ranger::ranger()`, and $0.5^\circ \times 0.5^\circ$ `spatial_block_cv_acc`. |
| **Triage Coalescing & Disambiguation (`06`–`07`)** | `recommended_determination`, `taxonomic_treatment_broad`, `taxonomic_treatment_subtype`, `is_completed_voucher`, `can1`, `can2`, `can3` | Exported by `07_triage_dashboard_synthesis.R` (`triage_queue.csv`, `final_curation_triage.csv`). Step 07 joins `multimodal_conflict_flags.csv` (georeferenced subset) and coalesces non-georeferenced vouchers directly from `morphometric_flags.csv` (Step 04) and `label_noise_audit.csv` (Step 05) using `unique(c("catalogNumber", ...))`, achieving **98.4% voucher completion (`5,117 / 5,201` vouchers)**. Cross-modal known-species reassignments (`both_alt_agree`, `morph_edaph_alt_agree`, `vis_edaph_alt_agree`, `is_severe_misid`) route to `"Severe_Triple_Stream_Conflict"` $\rightarrow$ `Misidentification_Severe` $\rightarrow$ `Misidentification_Reassign` **before** evaluating `Hybrid_Intergrade_Swarm` so sympatric misidentifications are never misclassified as hybrids. `Putative_Cryptic_Species` (`is_crypt`) requires **`morph_tier_source == "Tier_1_Pure"` (`n_tier1_leaves >= 1`)**, `morphometric_outlier_flag == TRUE`, absence of high-confidence known-taxon assignment (`cda_posterior_prob < 0.70` and `vision_c_error < 0.85`), and placement in peripheral GMM clusters (`gmm_cluster > 4`, $|\Delta\text{DOY}| \ge 14.0\text{ d}$) or extreme allochrony ($|\Delta\text{DOY}| \ge 21.0\text{ d}$, `cda_posterior_prob < 0.55`). |

---

## 5. Atomic Persistence & Schema Verification

* **Write-Then-Replace Pattern:** Every tabular export (`.csv`), JSON manifest (`.json`), and model checkpoint write must write to a temporary `.tmp` file in the target directory, flush buffers (`f.flush()` and `os.fsync(f.fileno())` in Python; atomic `file.rename()` from `.tmp` in R), and atomically replace the destination via `Path.replace()` (or `os.replace()`).
* **Zero Partial Artifacts:** Never write directly to final output paths in-place where an interrupted process or OOM signal could leave truncated CSVs or corrupted checkpoints.
* **Upstream Contract Verification:** Synthesis stages (such as `07_triage_dashboard_synthesis.R`) must invoke `validate_upstream_contracts()` at initialization to ensure all upstream tables comply with required primary keys and columns prior to computing downstream triage calls.

---

## 6. Phase 2 Native-DPI Basal Filtering, Cauline Rejection & Geometry Invariants

* **Dual-Class Stage 1 Box Processing:** Stage 1 YOLOv5 (`LeafPriority.pt`) must process both `Class 0` (`leaf_whole`) and `Class 1` (`leaf_partial`) candidate boxes to capture rosette-attached blades on native-DPI scans (`>1000px`), with cauline pre-filtering (`cy < 0.40 * sheet_h`) and fragment rejection (`max(bh, bw) < int(200 * scale_f)`).
* **Post-Gatekeeper Basal Filtering & Strengthened Cauline/Stem Rejection (`is_cauline_or_stem_contaminant`):** High-DPI basal candidates must strictly enforce `cy_rel >= 0.42`, minimum area floors (`>= 12000 * scale_f^2` for Tier 1, `>= 18000 * scale_f^2` and `cy_rel >= 0.48` for Tier 2), boxy clump rejection (`min(w,h)/max(w,h) > 0.78` when `solidity < 0.84`), sheet-relative area floor (`>= 0.28 * max_sheet_area`), and top-6 cap per voucher ranked by $Q = \text{solidity} \cdot \text{ucs} \cdot \sqrt{\text{area}}$ with automatic cleanup of dropped mask and contour files from disk. Furthermore, `is_cauline_or_stem_contaminant()` rejects upper-stem sessile/clasping cauline bracts (`aspect_ratio > 3.8` and `solidity < 0.78`, or displaced `> 2.2 * med_rosette_radius` above the basal rosette zone) even on herbarium sheets where capitula were not detected (`capitula_count == 0`).
* **Botanical Geometry Conventions:**
  * **Resolution-Adaptive Cross-Sectional Width Profiling:** `_sample_cross_sectional_widths()` and `trim_petiole_tail()` enforce `w_span = max(min(max(h, w), 350), int(0.65 * max(h, w)))` so native-DPI basal blades ($>700\text{ px}$ wide) are never truncated.
  * **Resolution-Scaled Petiole Width Thresholds:** `_disambiguate_base_apex_endpoints()` scales `eff_max_petiole_w` when `axis_len > 450.0`, and `trim_petiole_tail()` scales `max_petiole_width = max(16.0, 0.105 * blade_width_est)` and `cordate_guard_w = max(..., 0.135 * blade_width_est)` when `blade_width_est > 250.0 and min(h, w) > 500`, preserving exact backward compatibility on $400 \times 400$ synthetic unit-test fixtures.
  * **Spatial Rosette Centroid Polarity & Base Homologization:** `02_segment_and_extract.py` computes `rosette_centroids_by_pid` and passes `rosette_centroid` into `classify_and_route_leaf()` and `compute_geometric_metrics()`, while anchoring all contours via `petiole_attachment_pt=p_base` in every `resample_normalized_contour()` call across Tier 1, folded, and Tier 2 reflected leaves. Diagonal-crease fold detection requires `has_diagonal_crease and solidity < 0.72`.

---

## 7. Canonical Two-Tier Storage Architecture (`data/` vs. `runs/<RUN_ID>/`) & Orchestrator CLI

The pipeline implements a strict **Two-Tier Storage Architecture** separating shared invariant Step 01–02 specimen/organ assets at the repository root (`data/`) from run-scoped Step 03–07 statistical analysis directories (`runs/<RUN_ID>/` with a `runs/latest` symlink):

* **Tier A — Canonical Shared Store (`data/` at repository root):** Houses all invariant Step 01 and Step 02 specimen, environmental, and organ extraction assets:
  * `data/raw_vouchers/` & `data/raw_images/` (herbarium sheet image cache)
  * `data/environmental/` (`fine_grain_edaphic_cache.csv`, `polaris_30m/`, `solus_100m/`, `dem_3dep_30m/`, `ssurgo_cache/`, WorldClim v2.1, SoilGrids 250m)
  * `data/masks/tier1_intact/` & `data/masks/tier2_reflected/`
  * `data/contours/` (standardized 120-vertex contours)
  * `data/cropped_patches/{rosettes_dense,basal_leaves_raw,capitula}/`
  * The 5 canonical Step 01–02 CSV manifests in `data/tables/`: `curated_vouchers.csv`, `extracted_leaf_manifest.csv`, `failed_qc_vouchers.csv`, `voucher_reproductive_metrics.csv`, and `voucher_phenological_states.csv`.
* **Tier B — Run-Scoped Analysis Directories (`runs/<RUN_ID>/` with `runs/latest` symlink):** Automatically seeds Step 01–02 CSV tables and symlinks `data/masks`, `data/contours`, and `data/cropped_patches/` directly from canonical root `data/` without referencing previous `runs/` directories, while isolating all Step 03–07 statistical tables, reports, logs (`logs/pipeline.log`), and figures (`outputs/figures/`) inside `runs/<RUN_ID>/`:

```text
Packera_dubia_morphometrics/
├── data/                                    # Tier A: Canonical Shared Store (Steps 01–02 & Environmental)
│   ├── raw_vouchers/                        # Shared read-only herbarium sheet cache (never duplicated)
│   ├── environmental/                       # WorldClim, SoilGrids, 30m POLARIS, 100m SOLUS100, 3DEP & fine_grain_edaphic_cache.csv
│   ├── masks/{tier1_intact,tier2_reflected}/# Canonical binary leaf silhouettes
│   ├── contours/                            # Canonical 120-vertex (x, y) CSV contours
│   ├── cropped_patches/{basal_leaves_raw,capitula,rosettes_dense}/
│   └── tables/                              # Canonical Step 01–02 CSV manifests
│       ├── curated_vouchers.csv
│       ├── extracted_leaf_manifest.csv
│       ├── failed_qc_vouchers.csv
│       ├── voucher_reproductive_metrics.csv
│       └── voucher_phenological_states.csv
└── runs/                                    # Tier B: Run-Scoped Analysis Directories (Steps 03–07)
    ├── latest -> run_20261024_143000        # Symlink pointing to the most recent run directory
    └── run_20261024_143000/                 # Self-contained run directory (${RUN_DIR})
        ├── data/
        │   ├── masks -> ../../../data/masks # Symlinked from canonical root data/
        │   ├── contours -> ../../../data/contours
        │   ├── cropped_patches/             # Symlinked or taxon-filtered from canonical root data/
        │   └── tables/
        │       ├── curated_vouchers.csv     # Seeded (and optionally taxon-filtered) from root data/tables/
        │       ├── extracted_leaf_manifest.csv
        │       ├── failed_qc_vouchers.csv
        │       ├── voucher_reproductive_metrics.csv
        │       ├── voucher_phenological_states.csv
        │       ├── leaf_efa_harmonics.csv
        │       ├── fourier_pcs.csv
        │       ├── symmetric_fourier_pcs.csv
        │       ├── morphometric_flags.csv
        │       ├── label_noise_audit.csv
        │       └── multimodal_conflict_flags.csv
        ├── outputs/
        │   ├── figures/                     # Dual vector .pdf and 300-DPI .png publication figures
        │   └── reports/
        │       ├── final_curation_triage.csv
        │       ├── multimodal_spatial_rf_summary.csv
        │       ├── gmm_bayes_factors_summary.csv
        │       ├── Packera_dubia_Taxonomic_Revision_Summary.md
        │       └── run_metadata.json        # Run ID, timestamps, git SHA, CLI args, & table row counts
        └── logs/
            └── pipeline.log                 # Complete execution log for the run
```

* **`run_pipeline.sh` CLI Flags & Execution Modes:**
  * `--run-id <NAME>`: Custom run directory under `runs/<NAME>/`
  * `--skip-harvest`: Reuse cached `data/tables/curated_vouchers.csv` and local images in `data/raw_vouchers/`
  * `--start-step <1..7>`: Resume execution from any pipeline step (e.g., `--start-step 3` or `--start-step 6`) seeded automatically from canonical root `data/`
  * `--exclude-taxon "<TAXON>"`: Exclude one or more taxa (e.g., `--exclude-taxon "Packera plattensis"`), filtering `curated_vouchers.csv`, Step 02 manifests, and `cropped_patches/` while dynamically remapping Step 05 class indices to contiguous `0..K-1` (`active_taxa`)
  * `--focal-taxon "<TAXON>"` & `--gold-anchors <INT>` (default `30`): **Strategy 1 (Focal Taxon + Outgroup Gold Reference Anchors)**—processes all focal taxon vouchers plus `N` `Tier_1_Gold` anchors per outgroup species for Steps 04–06 calibration, then filters Step 07's final triage queue and figures strictly to the focal taxon
  * `--focal-taxon "<TAXON>" --intraspecific-mode`: **Strategy 2 (Pure Intraspecific Delimitation Mode)**—runs strictly on a single species (`n_unique_taxa == 1`) with zero outgroup specimens, switching Steps 04–06 from `species_standardized` to label-blind GMM morphotypes (`Morphotype_K1..K_n` or geographic/edaphic ecotypes)
  * `--gold-only` & `--determiner-tier <TIER>`: Filters `curated_vouchers.csv`, all Step 02 manifests, and `cropped_patches/` in `${RUN_DIR}/data/` strictly to rows where `determiner_tier == 'Tier_1_Gold'` (or `<TIER>`) before running Steps 02–07
  * `--exclude-misidentified` & `--triage-source <RUN_ID_OR_ALL>` (default `ALL` across `runs/*/outputs/reports/final_curation_triage.csv`): Filters out any `catalogNumber` with `taxonomic_treatment_subtype == 'Misidentification_Reassign'` (retaining only `Confirmed_Species`, `Hybrid_Introgressant`, `Subspecies_Ecotype`, `Ecophenotypic_Plasticity`, and `Putative_Cryptic_Species`) across `curated_vouchers.csv`, Step 02 manifests, and `cropped_patches/` during preflight
  * `--force` / `--overwrite`, `--resegment` / `--force-segment`, `--limit <N>`, and `--sam2-refine`
* **Backward Compatibility Guard:** When any phase script is invoked directly by unit tests without `PACKERA_RUN_DIR` or `--run-id`, it honors explicit CLI path arguments (or default `data/` / `outputs/` paths) so all unit tests pass unchanged.

---

## 8. Publication Figure Catalog, 3D Ordinations & Universal Unfilled Marker Standard

Every pipeline run exports publication-ready figures in both vector `.pdf` and 300-DPI `.png` formats inside `${RUN_DIR}/outputs/figures/` using the colorblind-safe Okabe-Ito palette (`#0072B2` *Packera anonyma*, `#D55E00` *P. dubia*, `#009E73` *P. paupercula*, `#CC79A7` *P. plattensis*) and italicized taxon labels, enforcing two universal plotting standards:

* **3D PCA & CDA Ordination Standard (`PC1 x PC2 x PC3`, `Sym_PC1 x Sym_PC2 x Sym_PC3`, `Can1 x Can2 x Can3`):** Rendered via `plot_3d_scatter_box()` in `scripts/morphometrics/pub_figures_utils.R`, projecting 3D specimen coordinates with vertical drop-lines to the floor grid and exporting `can1`, `can2`, and `can3` in `morphometric_flags.csv` and `triage_queue.csv` (with orthogonalized `PC3`/`PC2` fallback via `orth_against` when $g - 1 < 3$).
* **Universal Unfilled (Hollow) Scatterplot Marker Standard:** All scatterplots across R and Python use unfilled (hollow) open markers with colored strokes (`pch = 1, 2, 5, 0` or `bg = NA` in base R; `shape = 21..25, fill = NA` in `ggplot2`; `facecolors = 'none'` in `matplotlib`) to prevent overplotting occlusion in dense morphospace and ecological ordinations.

| Figure Stem (`.pdf` & `.png`) | Producing Script | Biological / Statistical Hypothesis Visualized |
| :--- | :--- | :--- |
| `fig02a_segmentation_qc_distributions` | `scripts/pipeline/02_segment_and_extract.py` | (A) Lamina solidity vs. Unoccluded Completeness Score (UCS) with unfilled markers (`facecolors='none'`) by routing tier (`Tier_1_Direct`, `Tier_2_Reflected`) against gatekeeper thresholds (`0.72`, `0.85`); (B) Basal leaf yield per voucher & QC failure breakdown. |
| `fig02b_silhouette_montage` | `scripts/pipeline/02_segment_and_extract.py` | 16-tile montage of representative extracted binary basal leaf silhouettes across `Tier_1_Direct` and `Tier_2_Reflected` with homologized petiole base anchor (`p_base`) and midrib axis overlaid. |
| `fig03a_efa_morphospace_pca` | `scripts/morphometrics/03_fourier_extractor.R` | (A) 3D Full 12-harmonic Elliptic Fourier PCA (`PC1 x PC2 x PC3`) and (B) 3D bilateral-symmetric (`Momocs::rm_asym`) morphospace (`Sym_PC1 x Sym_PC2 x Sym_PC3`) fitted on `Tier_1_Direct` leaves (`plot_3d_scatter_box`, unfilled markers). |
| `fig03b_harmonic_power_and_asymmetry` | `scripts/morphometrics/03_fourier_extractor.R` | (A) Cumulative Fourier harmonic power spectrum (`n = 1..12`) and (B) bilateral asymmetry (`pct_asymmetry`) vs. high-frequency perimeter roughness (`harmonic_roughness`) gating contour outliers (`pch = 1, 2`). |
| `fig04a_gmm_bic_and_uncertainty` | `scripts/morphometrics/04_gmm_morphotools.R` | (A) Unsupervised `mclust` BIC model selection curves across $G = 1..9$ on `Tier_1_Pure` Gold+Silver anchors; (B) 3D morphospace (`Sym_PC1 x Sym_PC2 x Sym_PC3`) colored by posterior cluster with point size scaled by taxonomic uncertainty $U_i = 1 - \max(P_{ik})$. |
| `fig04b_cda_biplot_and_allometry` | `scripts/morphometrics/04_gmm_morphotools.R` | (A) 3D Canonical Discriminant Analysis (`Can1 x Can2 x Can3`) ordination with honest Leave-One-Out Cross-Validation (`loocv_accuracy`) and (B) vegetative-to-reproductive allometric scaling ($\log_{10}$ leaf size vs. $\log_{10}$ involucre scale). |
| `GradCAM_audit_panel` | `scripts/analysis/05_cleanlab_vision_xai.py` | Side-by-side background-masked `Tier_1_Direct` basal leaf crops (`[242, 239, 230]` neutral herbarium paper outside 5 px dilated `tier1_intact` mask) and DINOv2 ViT-B/14 $16 \times 16$ patch-token attribution heatmaps for Cleanlab-flagged label noise candidates. |
| `fig05b_cleanlab_quality_and_tomentum` | `scripts/analysis/05_cleanlab_vision_xai.py` | (A) Cleanlab label quality ($1 - c_{\text{error}}$) distributions across taxa and determiner authority tiers (`Gold`, `Silver`, `Bronze`); (B) CIE $L^*$ luminance vs. HSV saturation (`facecolors='none'`) on masked `Tier_1_Direct` basal leaves separating adaxial vs. abaxial arachnoid tomentum and `tomentum_taxon_discordance`. |
| `spatial_rf_niche_importance` | `scripts/analysis/06_multimodal_spatial_rf.R` | Spatial Random Forest variable importance across 10m–30m SSURGO/POLARIS/SOLUS100/3DEP edaphic & rock-outcrop predictors, WorldClim bioclimatic variables, and phenological anomaly residuals under $0.5^\circ \times 0.5^\circ$ spatial block CV. |
| `fig06b_geographic_map_and_hopkins_phenology` | `scripts/analysis/06_multimodal_spatial_rf.R` | (A) Eastern/central North America geographic voucher distribution (`fill = NA` hollow markers) colored by taxon and shaped by `consensus_conflict_flag`; (B) Collection DOY vs. Latitude against the centered Hopkins Bioclimatic Law baseline ($+4.0\text{ d/}^\circ\text{N}$). |
| `fig07_multimodal_discordance_summary` | `scripts/analysis/07_triage_dashboard_synthesis.R` | (A) Modality-specific conflict rates across Morphology (CDA/GMM), DINOv2 Vision (Cleanlab), and Ecological Niche (Spatial RF); (B) Consensus curatorial triage priority breakdown across the active *Packera* taxa. |
| `Figure_Integrative_Packera_dubia_Revision` | `scripts/analysis/07_triage_dashboard_synthesis.R` | 4-panel synthesis dashboard combining (A) 3D CDA morphospace (`Can1 x Can2 x Can3`), (B) DINOv2 Cleanlab confidence vs. morphological distance (`pch = 1, 4`), (C) Latitudinal vs. flowering DOY allochrony (`pch = 1, 2`), and (D) Taxonomic treatment summary. |

---

## Workspace Skills Reference

Activate the repository's specialized skills in `.agents/skills/` whenever working on the corresponding domains:

* **`/morphometric-cv-guardrails`** (`.agents/skills/morphometric-cv-guardrails/SKILL.md`): Mandatory whenever modifying scripts in `scripts/vision/`, `scripts/pipeline/`, `scripts/annotation_and_training/`, or `scripts/morphometrics/`.
* **`/pipeline-contract-audit`** (`.agents/skills/pipeline-contract-audit/SKILL.md`): Mandatory before committing code changes or when debugging cross-stage Python/R pipeline handshakes.

