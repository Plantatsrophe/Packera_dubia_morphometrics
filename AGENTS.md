# Workspace Rules & Non-Negotiable Invariants (`Packera_dubia_morphometrics`)

**Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV5`  
**Affiliation:** University of North Carolina at Chapel Hill Herbarium (NCU)

This file defines the always-on architectural, mathematical, and cross-language data contracts for the `Packera_dubia_morphometrics` repository. Every AI agent and contributor modifying Python (`.py`) or R (`.R`) code in this workspace must enforce the non-negotiable invariants below without exception.

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
  * **`.venv` (Primary Orchestration & Modern Vision):** Houses Python orchestration (`main.py`), `DINOv2` (`dinov2_vitb14_reg`), `cleanlab` confident learning, and general scientific utilities.
  * **`.venv_LM2` (LeafMachine2 & Detectron2 PointRend Isolation):** Houses stock LeafMachine2 v2.2+, YOLOv5 (`PLANT_LeafPriority` / `LeafPriority.pt` and the 38-class Ruler Detector), and compiled Detectron2 PointRend (`leaf_seg_rcnn_v2_1`) with strict ABI pins: `numpy<2.0` and `opencv-contrib-python<=4.10.0.84`. Never upgrade `numpy>=2.0` or `opencv>4.10.0.84` inside `.venv_LM2`.

---

## 4. Canonical Cross-Language CSV Contracts

All Python (`scripts/pipeline/`, `scripts/vision/`, `scripts/analysis/`) and R (`scripts/morphometrics/`, `scripts/analysis/`) stages must adhere to the exact primary keys and shared column schemas below:

| Contract Element | Canonical Column / Pattern | Cross-Language Rule & Semantics |
| :--- | :--- | :--- |
| **Ingestion Schema & EPA Level III Ecoregions** | `catalogNumber`, `institutionCode`, `scientificName`, `species_raw`, `recordedBy`, `recordNumber`, `identifiedBy`, `determiner_raw`, `determiner_tier`, `type_status`, `county`, `stateProvince`, `decimalLatitude`, `latitude`, `decimalLongitude`, `longitude`, `coordinateUncertainty`, `year`, `month`, `day`, `eventDate`, `regional_group`, `image_path`, `exsiccatae_key`, `is_primary_duplicate`, `duplicate_count` | Exported by `01_voucher_harvester.py` (`curated_vouchers.csv`). Standardizes dynamic taxon binomials/basionyms across all 24+ Eastern North American *Packera* taxa. Evaluates coordinates against official US EPA Level III Ecoregions (`data/environmental/epa_ecoregions_l3/`, Omernik & Griffith 2014); for ungeoreferenced sheets, assigns ecoregion via county polygon containment ($\ge 95\%$ of county land area in a single ecoregion) or marks `NA`. |
| **Composite Grouping Key** | `(catalogNumber, plant_individual_id)` | Every leaf-level and specimen-level aggregation across Python and R must group by `(catalogNumber, plant_individual_id)`. **Never** collapse multi-plant herbarium sheets to `plant_individual_id = 0` or group by `catalogNumber` alone when aggregating organ traits. |
| **Routing Tier Column** | `assigned_tier` | Canonical closed-silhouette EFA tier is `"Tier_1_Direct"` (with `"Tier_2_Reflected"` retained strictly as an auxiliary non-EFA linear/count diagnostic tier). In R, always match case-insensitively via `grepl("^tier_?[12]", assigned_tier, ignore.case = TRUE)`. |
| **Scale Calibration & Botanical Geometry Invariants** | `pixels_per_mm`, `is_calibrated_mm`, `petiole_length_px`, `petiole_length_mm`, `plant_height_px`, `plant_height_mm`, `is_bent_stem` | `pixels_per_mm` stores empirical ruler scale (`float` or `NA`) and `is_calibrated_mm` (`bool`) flags absolute metric calibration. Phase 2 computes per-leaf petiole length $L_{\text{petiole}} = \max(d_{\text{geodesic}}(p_{\text{base}}, p_{\text{crown}}), \|p_{\text{base}} - p_{\text{crown}}\|_2)$ from crown point $p_{\text{crown}}$, aggregating outer-rosette leaves at the 75th percentile ($P_{75}$). Scape height $H_{\text{plant}} = \|p_{\text{bend}} - p_{\text{crown}}\|_2 + \|p_{\text{corymb}} - p_{\text{bend}}\|_2$ accounts for bent stem posture. Post-QC montage is regenerated via `--regenerate-figures-only`. |
| **Strict Active/Passive Authority Split & Dawid-Skene Weighting** | `determiner_tier`, `authority_weight` | Enforces a strict active/passive authority split across Steps 04–06 (Šlenker et al. 2022; Dawid & Skene 1979): `Tier_1_Gold` (monographic/type-verified) and `Tier_2_Silver` (regional specialist) vouchers form the active reference anchor set (`authority_weight = 1.0` or empirical Dawid-Skene precision $w_k = 1 - \text{error\_rate}_k$ after per-taxon $\chi^2_{0.95}$ Mahalanobis outlier pruning), whereas `Tier_3_Bronze` vouchers receive `authority_weight = 0.0` during model fitting and are projected passively (`predict()` / `predict_proba()`) into the fitted morphospace, CDA canonical space, Cleanlab classifier, and Spatial Random Forests. |
| **Strict `Tier_1_Direct` EFA, Inertial Pre-Alignment, Symmetrization & OTU Medians** | `harmonic_roughness`, `is_contour_outlier`, `leaf_quality_weight`, `pct_asymmetry`, `sym_var`, `asym_var`, `Sym_PC1`–`Sym_PC5`, `PC1`–`PC5`, `n_tier1_leaves`, `n_tier2_leaves`, `morph_tier_source` | Exported by `03_fourier_extractor.R` (`fourier_pcs.csv`, `leaf_efa_harmonics.csv`). Only intact `Tier_1_Direct` basal leaf contours enter `Momocs::efourier` (never pooling mirrored hemi-blades or open-curve approximations into EFA). Contours undergo inertial axis pre-alignment (`Momocs::coo_alignxax` with circularity $> 0.75$ suborbicular guard), empirical $\ge 99\%$ cumulative harmonic power calibration at $H = 12$ (`calibrate_harmonicpower_efourier`), uniform bilateral symmetrization (`Momocs::rm_asym(ef)`, Claude 2008), and Minimum Covariance Determinant (MCD) Mahalanobis outlier filtering (`MASS::cov.rob`, $p < 0.001$ combined with `pct_asymmetry > 0.22 | harmonic_roughness > 0.15`), followed by **unweighted specimen-level (OTU) median** aggregation (`morph_tier_source = "Tier_1_Pure"`). |
| **Klingenberg Allometry, Dual Clustering Architecture & MANOVA** | `can1`, `can2`, `can3`, `cda_predicted_taxon`, `cda_posterior_prob`, `cda_second_posterior_prob`, `cda_posterior_margin`, `morphometric_outlier_flag`, `loocv_accuracy`, `morph_tier_source` | Exported by `04_gmm_morphotools.R` (`morphometric_flags.csv`, `allometry_homogeneity_test.csv`, `gmm_bayes_factors_summary.csv`, `institutional_batch_effect_audit.csv`, `ecoregion_supervised_cda_summary.csv`, `manova_group_significance.csv`, `manova_pairwise_posthoc.csv`). Klingenberg allometry tests slope homogeneity (`PC ~ log_centsize * species`) on calibrated anchors. Retains Horn / Broken-Stick PCs capped by Williams & Titus (1988) $P \le \lfloor N_{\min} / 3 \rfloor$. Executes Dual Clustering Architecture: Track 1 Unsupervised GMM Morphotype CDA (`fig04_gmm_cda_clustering.{pdf,png}`) and Track 2 Supervised EPA Level III Ecoregion CDA (`fig04c_supervised_ecoregion_cda.{pdf,png}`). Computes parametric omnibus and post-hoc pairwise MANOVA with BH FDR. Institutional batch audit confirms partial $R^2 < 0.05$. |
| **Indumentum-Blind Vision XAI, Foreground Spatial Token Pooling & OOD Scoring** | `predicted_label` (`vision_predicted_label`), `c_error`, `is_label_corrupted`, `surface_orientation`, `tomentum_taxon_discordance`, `vision_ood_score`, `vision_is_ood`, `foreground_attention_ratio`, `is_background_shortcut` | Exported by `05_cleanlab_vision_xai.py` (`data/tables/label_noise_audit.csv`). Applies a $\sigma = 1.0\text{ px}$ Gaussian low-pass filter at $224 \times 224$ and composites 5-px morphologically dilated `Tier_1_Direct` leaf masks onto neutral herbarium paper `[242, 239, 230]`. Extracts `dinov2_vitb14_reg` features via foreground-masked $16 \times 16$ spatial patch token average pooling (`compute_foreground_patch_mask`), trains out-of-fold L2-penalized `LogisticRegression(C=0.1, class_weight="balanced")` inside a fold-internal `StandardScaler` + `PCA(0.95)` pipeline strictly on Gold/Silver anchors using `institutionCode`-grouped `StratifiedGroupKFold` with Platt/temperature calibration, runs class-conditional Cleanlab `prune_by_noise_rate`, computes feature-space Mahalanobis + cosine $k$-NN OOD scores (`vision_ood_score`), and enforces Pointing Game `foreground_attention_ratio >= 0.65`. |
| **$0\text{–}15\text{ cm}$ Rhizosphere/Outcrop ENM Stack, Collinearity Pruning & Niche Equivalency** | `spatial_block_cv_acc`, `edaphic_best_fit_taxon`, `multimodal_concordance`, `triage_category`, `ssurgo_rock_outcrop_pct`, `ssurgo_bedrock_depth_cm`, `fine_soil_ph`, `fine_soil_sand`, `fine_bedrock_depth_cm`, `topo_tpi_30m`, `topo_twi_30m`, `sgmc_lithology_class` | Exported by `06_multimodal_spatial_rf.R` (`multimodal_conflict_flags.csv`, `multimodal_spatial_rf_summary.csv`, `edaphic_niche_equivalency.csv`) and cached in `data/environmental/fine_grain_edaphic_cache.csv`. Integrates $0\text{–}15\text{ cm}$ trapezoidal rhizosphere depth-weighting (`(5 * v_0_5 + 10 * v_5_15) / 15`) with strict non-blended coalesce hierarchy (100% SSURGO / 30m POLARIS / 100m SOLUS100 in CONUS; 250m SoilGrids strictly non-CONUS fallback when both NA; raw `soil_*` excluded from active `env_vars`), depth-to-bedrock/CFVO/AWC/clr-texture (`clr_sand_silt`, `clr_clay_geom`), 50m-buffered SSURGO outcrop queries, 30m 3DEP micro-topography, 6-class USGS SGMC bedrock lithology, coordinate-uncertainty buffered focal extraction (`100–1000 m`), $|r| < 0.70$ and $\text{VIF} < 5.0$ collinearity pruning, partitioned Climatic ($D_{\text{clim}}, I_{\text{clim}}$) vs. Edaphic-Geologic ($D_{\text{edaph}}, I_{\text{edaph}}$) Niche Equivalency & Background Similarity tests, and $0.5^\circ \times 0.5^\circ$ Spatial Block CV (`spatial_block_cv_acc`) trained on active Gold/Silver anchors. |
| **Reliability-Weighted Consensus, 5-Stage Decision Matrix & Complex Discovery (`07`–`08`)** | `recommended_determination`, `taxonomic_treatment_broad`, `taxonomic_treatment_subtype`, `is_completed_voucher`, `can1`, `can2`, `can3` | Exported by `07_triage_dashboard_synthesis.R` (`triage_queue.csv`, `final_curation_triage.csv`) and `08_discover_taxon_complexes.R` (`discovered_taxon_complexes.csv`, located at `scripts/analysis/08_discover_taxon_complexes.R`). Step 07 weights each modality by out-of-fold Cohen's $\kappa_m$ ($w_m = \max(0.05, \kappa_m)$) dynamically renormalized over observed modalities $\mathcal{O}_i$ (`98.4%` voucher completion), computes Hopkins-adjusted circular anthesis phenology (`atan2` residuals controlling for latitude, elevation, and MAT gated on active anthesis), runs taxon-adjusted `institutionCode` ANCOVA with BH FDR, and applies the 5-stage hierarchical decision matrix (`Misidentification_Reassign` evaluated before `Hybrid_Introgressant`; multi-modal Integration by Congruence for `Putative_Cryptic_Species`). Step 08 performs unsupervised Leiden/Louvain graph modularity and Ward.D2 clustering to discover Eastern North American *Packera* species complexes, resolving Taxon x GMM Morphotype subclusters when $N_{\text{taxa}} \le 2$. |

---

## 5. Atomic Persistence & Schema Verification

* **Write-Then-Replace Pattern:** Every tabular export (`.csv`), JSON manifest (`.json`), and model checkpoint write must write to a temporary `.tmp` file in the target directory, flush buffers (`f.flush()` and `os.fsync(f.fileno())` in Python; atomic `file.rename()` from `.tmp` in R), and atomically replace the destination via `Path.replace()` (or `os.replace()`).
* **Zero Partial Artifacts:** Never write directly to final output paths in-place where an interrupted process or OOM signal could leave truncated CSVs or corrupted checkpoints.
* **Upstream Contract Verification:** Synthesis stages (such as `07_triage_dashboard_synthesis.R`) must invoke `validate_upstream_contracts()` at initialization to ensure all upstream tables comply with required primary keys and columns prior to computing downstream triage calls.

---

## 6. Phase 2 Stock LeafMachine2 v2.2+, Classical Gatekeeping & Geometry Invariants

* **Stock LeafMachine2 v2.2+ Configuration (`configure_leafmachine2.py`):** Uses stock LeafMachine2 v2.2+ `PLANT_LeafPriority` (`minimum_confidence_threshold: 0.65`) and `leaf_seg_rcnn_v2_1` weights out of the box without requiring custom PCD retraining or manual SAM 2 annotation, enforcing `treat_leaflet_as_leaf_whole: False`, `use_efds_for_png_masks: False`, and `do_apply_landmarks: False`.
* **Otsu Boundary-Continuity Lobe Rejection & Calibrated Area Floor:** `is_truncated_lobe_or_leaflet()` rejects unattached terminal lobes and severed fragments where non-paper tissue occupies $> 15\%$ of the bounding-box perimeter (`max_tissue_border_ratio: 0.15`), and calibrated sheets enforce a physical blade area floor of `min_blade_area_cm2: 0.75`.
* **Classical Dimensionless Descriptors & Dissection Gates:** `classify_and_route_leaf()` evaluates classical dimensionless shape descriptors (`Solidity >= 0.72`, `Circularity >= 0.25`, `Extent >= 0.35`, `UCS >= 0.85`), routing pinnatifid blades ($0.48 \le \text{Solidity} < 0.72$) and bipinnatisect blades ($0.18 \le \text{Solidity} < 0.48$ with $\ge 4$ deep sinuses) via `is_botanical_dissection()`.
* **Post-Gatekeeper Basal Filtering & Strengthened Cauline/Stem Rejection (`is_cauline_or_stem_contaminant`):** High-DPI basal candidates enforce `cy_rel >= 0.42`, minimum area floors, boxy clump rejection, sheet-relative area floor (`>= 0.28 * max_sheet_area`), and top-6 cap per voucher with automatic disk cleanup, while `is_cauline_or_stem_contaminant()` rejects upper-stem sessile/clasping cauline bracts (`aspect_ratio > 3.8` and `solidity < 0.78`, or displaced `> 2.2 * med_rosette_radius` above the basal rosette zone).
* **Botanical Geometry Conventions:**
  * **Resolution-Adaptive Cross-Sectional Width Profiling:** `_sample_cross_sectional_widths()` and `trim_petiole_tail()` enforce `w_span = max(min(max(h, w), 350), int(0.65 * max(h, w)))` so native-DPI basal blades ($>700\text{ px}$ wide) are never truncated.
  * **Resolution-Scaled Petiole Width Thresholds:** `_disambiguate_base_apex_endpoints()` scales `eff_max_petiole_w` when `axis_len > 450.0`, and `trim_petiole_tail()` scales `max_petiole_width = max(16.0, 0.105 * blade_width_est)` and `cordate_guard_w = max(..., 0.135 * blade_width_est)` when `blade_width_est > 250.0 and min(h, w) > 500`, preserving exact backward compatibility on $400 \times 400$ synthetic unit-test fixtures.
  * **Spatial Rosette Centroid Polarity & Base Homologization:** `02_segment_and_extract.py` computes `rosette_centroids_by_pid` and passes `rosette_centroid` into `classify_and_route_leaf()` and `compute_geometric_metrics()`, while anchoring all contours via `petiole_attachment_pt=p_base` in every `resample_normalized_contour()` call.

---

## 7. Canonical Two-Tier Storage Architecture (`data/` vs. `runs/<RUN_ID>/`) & Orchestrator CLI

The pipeline implements a strict **Two-Tier Storage Architecture** separating shared invariant Step 01–02 specimen/organ assets at the repository root (`data/`) from run-scoped Step 03–08 statistical analysis directories (`runs/<RUN_ID>/` with a `runs/latest` symlink):

* **Tier A — Canonical Shared Store (`data/` at repository root):** Houses all invariant Step 01 and Step 02 specimen, environmental, and organ extraction assets:
  * `data/raw_vouchers/` & `data/raw_images/` (herbarium sheet image cache)
  * `data/environmental/` (`epa_ecoregions_l3/`, `fine_grain_edaphic_cache.csv`, `polaris_30m/`, `solus_100m/`, `dem_3dep_30m/`, `ssurgo_cache/`, WorldClim v2.1, SoilGrids 250m)
  * `data/masks/tier1_intact/` & `data/masks/tier2_reflected/`
  * `data/contours/` (standardized 120-vertex contours)
  * `data/cropped_patches/{rosettes_dense,basal_leaves_raw,capitula}/`
  * The 5 canonical Step 01–02 CSV manifests in `data/tables/`: `curated_vouchers.csv`, `extracted_leaf_manifest.csv`, `failed_qc_vouchers.csv`, `voucher_reproductive_metrics.csv`, and `voucher_phenological_states.csv`.
* **Tier B — Run-Scoped Analysis Directories (`runs/<RUN_ID>/` with `runs/latest` symlink):** Automatically seeds Step 01–02 CSV tables and symlinks `data/masks`, `data/contours`, and `data/cropped_patches/` directly from canonical root `data/` without referencing previous `runs/` directories, while isolating all Step 03–08 statistical tables, reports, logs (`logs/pipeline.log`), and figures (`outputs/figures/`) inside `runs/<RUN_ID>/`:

```text
Packera_dubia_morphometrics/
├── data/                                    # Tier A: Canonical Shared Store (Steps 01–02 & Environmental)
│   ├── raw_vouchers/                        # Shared read-only herbarium sheet cache (never duplicated)
│   ├── environmental/                       # EPA Level III, WorldClim, SoilGrids, 30m POLARIS, 100m SOLUS100, 3DEP & fine_grain_edaphic_cache.csv
│   ├── masks/{tier1_intact,tier2_reflected}/# Canonical binary leaf silhouettes
│   ├── contours/                            # Canonical 120-vertex (x, y) CSV contours
│   ├── cropped_patches/{basal_leaves_raw,capitula,rosettes_dense}/
│   └── tables/                              # Canonical Step 01–02 CSV manifests
│       ├── curated_vouchers.csv
│       ├── extracted_leaf_manifest.csv
│       ├── failed_qc_vouchers.csv
│       ├── voucher_reproductive_metrics.csv
│       └── voucher_phenological_states.csv
└── runs/                                    # Tier B: Run-Scoped Analysis Directories (Steps 03–08)
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
        │       ├── edaphic_niche_equivalency.csv
        │       ├── gmm_bayes_factors_summary.csv
        │       ├── allometry_homogeneity_test.csv
        │       ├── institutional_batch_effect_audit.csv
        │       ├── ecoregion_supervised_cda_summary.csv
        │       ├── manova_group_significance.csv
        │       ├── manova_pairwise_posthoc.csv
        │       ├── discovered_taxon_complexes.csv
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
| `fig03_efa_morphospace_reconstruction` | `scripts/morphometrics/03_fourier_extractor.R` | Full 12-harmonic empirical EFA contour reconstructions and representative basal silhouettes across taxa. |
| `fig03a_efa_morphospace_pca` | `scripts/morphometrics/03_fourier_extractor.R` | (A) 3D Full 12-harmonic Elliptic Fourier PCA (`PC1 x PC2 x PC3`) and (B) 3D bilateral-symmetric (`Momocs::rm_asym`) morphospace (`Sym_PC1 x Sym_PC2 x Sym_PC3`) fitted strictly on `Tier_1_Direct` leaves (`plot_3d_scatter_box`, unfilled markers). |
| `fig03b_harmonic_power_and_asymmetry` | `scripts/morphometrics/03_fourier_extractor.R` | (A) Cumulative Fourier harmonic power spectrum (`n = 1..12` achieving $\ge 99\%$ cumulative power) and (B) bilateral asymmetry (`pct_asymmetry`) vs. high-frequency perimeter roughness (`harmonic_roughness`) gating contour outliers (`pch = 1, 2`). |
| `fig04_gmm_cda_clustering` | `scripts/morphometrics/04_gmm_morphotools.R` | Track 1 Unsupervised GMM Morphotype CDA ordination biplot (`Can1 x Can2 x Can3`) with 95% confidence ellipses and LOOCV accuracy. |
| `fig04a_gmm_bic_and_uncertainty` | `scripts/morphometrics/04_gmm_morphotools.R` | (A) Unsupervised `mclust` BIC model selection curves ($2\Delta\text{BIC}$ Kass-Raftery thresholds) across $G = 1..9$ on `Tier_1_Pure` Gold+Silver anchors; (B) 3D morphospace (`Sym_PC1 x Sym_PC2 x Sym_PC3`) colored by posterior cluster with point size scaled by taxonomic uncertainty $U_i = 1 - \max(P_{ik})$. |
| `fig04b_cda_biplot_and_allometry` | `scripts/morphometrics/04_gmm_morphotools.R` | (A) 3D Canonical Discriminant Analysis (`Can1 x Can2 x Can3`) ordination with honest Leave-One-Out Cross-Validation (`loocv_accuracy`) and (B) Klingenberg (2016) allometric scaling ($\log_{10}$ leaf size vs. shape PCs / involucre scale). |
| `fig04c_supervised_ecoregion_cda` | `scripts/morphometrics/04_gmm_morphotools.R` | Track 2 Supervised EPA Level III Ecoregion CDA biplot projecting ecotypes across physiographic provinces. |
| `shape_allometry_regression` | `scripts/morphometrics/04_gmm_morphotools.R` | Multivariate shape allometry regressions ($\log_{10}\text{CS}$ vs. PC1/PC2) testing Klingenberg slope homogeneity. |
| `GradCAM_audit_panel` | `scripts/analysis/05_cleanlab_vision_xai.py` | Side-by-side indumentum-blind `Tier_1_Direct` crops and `dinov2_vitb14_reg` $16 \times 16$ foreground patch-token attribution heatmaps (`foreground_attention_ratio >= 0.65`) for Cleanlab-flagged label noise candidates. |
| `fig05b_cleanlab_quality_and_tomentum` | `scripts/analysis/05_cleanlab_vision_xai.py` | (A) Cleanlab label quality ($1 - c_{\text{error}}$) distributions across taxa and determiner authority tiers (`Gold`, `Silver`, `Bronze`); (B) CIE $L^*$ luminance vs. HSV saturation (`facecolors='none'`) separating adaxial vs. abaxial arachnoid tomentum and `tomentum_taxon_discordance`. |
| `spatial_rf_niche_importance` | `scripts/analysis/06_multimodal_spatial_rf.R` | Spatial Random Forest variable importance across $0\text{–}15\text{ cm}$ rhizosphere edaphic predictors, USGS SGMC bedrock lithology, rock-outcrop variables, WorldClim bioclimatic variables, and partitioned Climatic vs. Edaphic niche equivalency under $0.5^\circ \times 0.5^\circ$ spatial block CV. |
| `fig06_ecological_niche_phenology` | `scripts/analysis/06_multimodal_spatial_rf.R` | Integrated 4-panel macroecological ordination combining bioclimatic PCA, edaphic-geologic PCA, and anthesis phenology. |
| `fig06b_geographic_map_and_hopkins_phenology` | `scripts/analysis/06_multimodal_spatial_rf.R` | (A) Eastern/central North America geographic voucher distribution (`fill = NA` hollow markers) colored by taxon and shaped by `consensus_conflict_flag`; (B) Collection DOY vs. Latitude against the Hopkins Bioclimatic Law baseline. |
| `fig07_multimodal_discordance_summary` | `scripts/analysis/07_triage_dashboard_synthesis.R` | (A) Modality-specific conflict rates weighted by out-of-fold Cohen's $\kappa$ across Morphology (CDA/GMM), DINOv2 Vision (Cleanlab), and Ecological Niche (Spatial RF); (B) Consensus curatorial triage priority breakdown across the active *Packera* taxa. |
| `Figure_Integrative_Packera_dubia_Revision` | `scripts/analysis/07_triage_dashboard_synthesis.R` | 4-panel synthesis dashboard combining (A) 3D CDA morphospace (`Can1 x Can2 x Can3`), (B) DINOv2 Cleanlab confidence vs. morphological distance (`pch = 1, 4`), (C) Hopkins-adjusted circular flowering DOY allochrony (`pch = 1, 2`), and (D) 5-stage taxonomic treatment summary. |
| `fig_empirical_taxon_complex_network` | `scripts/analysis/08_discover_taxon_complexes.R` | Unsupervised Leiden/Louvain community network graph resolving Eastern North American species complexes and Taxon x GMM Morphotype subclusters. |

---

## Workspace Skills Reference

Activate the repository's specialized skills in `.agents/skills/` whenever working on the corresponding domains:

* **`/morphometric-cv-guardrails`** (`.agents/skills/morphometric-cv-guardrails/SKILL.md`): Mandatory whenever modifying scripts in `scripts/vision/`, `scripts/pipeline/`, `scripts/annotation_and_training/`, or `scripts/morphometrics/`.
* **`/pipeline-contract-audit`** (`.agents/skills/pipeline-contract-audit/SKILL.md`): Mandatory before committing code changes or when debugging cross-stage Python/R pipeline handshakes.
