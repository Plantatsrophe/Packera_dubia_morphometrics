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
| **Fourier Quality & Symmetry** | `harmonic_roughness`, `is_contour_outlier`, `leaf_quality_weight`, `pct_asymmetry`, `sym_var`, `asym_var`, `Sym_PC1`–`Sym_PC5`, `PC1`–`PC5` | Exported by `03_fourier_extractor.R` (`fourier_pcs.csv`, `leaf_efa_harmonics.csv`). Voucher metadata merge occurs prior to PCA; both `Sym_PC1`–`Sym_PC5` (`Momocs::rm_asym(ef)`) and primary `PC1`–`PC5` use weighted covariance eigendecomposition (`stats::cov.wt(..., wt = tier_wt * leaf_quality_weight)`). Outlier threshold `pct_asymmetry > 0.22 | harmonic_roughness > 0.15` gates PCA fitting (`clean_pca_mask`); `leaf_quality_weight` weights specimen medoid aggregation. |
| **Morphometric Outlier Flag, CDA Posteriors & Validation** | `cda_predicted_taxon`, `cda_posterior_prob`, `cda_second_posterior_prob`, `cda_posterior_margin`, `morphometric_outlier_flag`, `loocv_accuracy` | Exported by `04_gmm_morphotools.R` (`data/tables/morphometric_flags.csv` / `morphometrics_misidentification_flags.csv`) with `morphometric_outlier_flag` aliased to `misidentification_flag` in Steps 06–07. Label-blind GMM fits on `Tier_1_Gold` + `Tier_2_Silver` anchors (`!mahal_outliers` via `MASS::cov.rob`) with parsimonious covariance restriction (`EII`..`VVI` when $N/P < 5$) and projects `Tier_3_Bronze` via `predict(gmm_fit, newdata = X_pcs)`. Passive CDA trains on Gold ($3.0$) and Silver ($1.5$) anchors with per-taxon 95% $\chi^2$ Mahalanobis outlier trimming (`qchisq(0.95, df = p)`), uniform class priors (`rep(1/g, g)`), and honest Leave-One-Out Cross-Validation (`loocv_accuracy`). |
| **Vision Noise & Tomentum Flag** | `predicted_label` (`vision_predicted_label`), `c_error`, `is_label_corrupted`, `tomentum_taxon_discordance` | Exported by `05_cleanlab_vision_xai.py` (`data/tables/label_noise_audit.csv`). Requires L2 cosine normalization of 768-d DINOv2 embeddings prior to fold-internal `StandardScaler`, authority weights (`3.0 / 1.5 / 0.5`) combined with `class_weight="balanced"` in `LogisticRegression`, and `tomentum_taxon_discordance` flagging where vision predictions diverge due to abaxial/adaxial vestiture and arachnoid tomentum variation. |
| **Spatial Model Validation & Hopkins Baseline** | `spatial_block_cv_acc`, `edaphic_best_fit_taxon`, `multimodal_concordance`, `triage_category` | Exported by `06_multimodal_spatial_rf.R` (`multimodal_conflict_flags.csv`, `multimodal_spatial_rf_summary.csv`). Uses authority-weighted circular phenology & edaphic envelopes, centers Hopkins' Bioclimatic Law baseline at `pred_at_min_lat + 4.0 * (decimalLatitude - min_lat)`, passes `case.weights = case_wts` (`3.0 / 1.5 / 0.5`) into `ranger::ranger()`, and evaluates $0.5^\circ \times 0.5^\circ$ `spatial_block_cv_acc` with resilient batched `soilDB::SDA_spatialQuery()` and SoilGrids 250m fallback. |
| **Triage Coalescing & Disambiguation (`06`–`07`)** | `recommended_determination`, `taxonomic_treatment_broad`, `taxonomic_treatment_subtype`, `is_completed_voucher` | Exported by `07_triage_dashboard_synthesis.R` (`triage_queue.csv`, `final_curation_triage.csv`). Step 07 joins `multimodal_conflict_flags.csv` (georeferenced subset) and coalesces non-georeferenced vouchers directly from `morphometric_flags.csv` (Step 04) and `label_noise_audit.csv` (Step 05) using `unique(c("catalogNumber", ...))`, achieving **98.4% voucher completion (`5,117 / 5,201` vouchers)**. Dual-stream alternate agreement (`both_alt_agree`) or high-confidence conflict (`is_severe_misid`) routes to `"Severe_Triple_Stream_Conflict"` $\rightarrow$ `Misidentification_Severe` $\rightarrow$ `Misidentification_Reassign` **before** evaluating `Hybrid_Intergrade_Swarm`. `Putative_Cryptic_Species` (`is_crypt`) is reserved strictly for true morphometric outliers (`morphometric_outlier_flag == TRUE`) lacking high-confidence known-taxon assignment (`cda_posterior_prob < 0.70` and `vision_c_error < 0.85`) in peripheral GMM clusters (`gmm_cluster > 4`, $|\Delta\text{DOY}| \ge 14.0\text{ d}$) or extreme allochrony ($|\Delta\text{DOY}| \ge 21.0\text{ d}$, `cda_posterior_prob < 0.55`). |

---

## 5. Atomic Persistence & Schema Verification

* **Write-Then-Replace Pattern:** Every tabular export (`.csv`), JSON manifest (`.json`), and model checkpoint write must write to a temporary `.tmp` file in the target directory, flush buffers (`f.flush()` and `os.fsync(f.fileno())` in Python; atomic `file.rename()` from `.tmp` in R), and atomically replace the destination via `Path.replace()` (or `os.replace()`).
* **Zero Partial Artifacts:** Never write directly to final output paths in-place where an interrupted process or OOM signal could leave truncated CSVs or corrupted checkpoints.
* **Upstream Contract Verification:** Synthesis stages (such as `07_triage_dashboard_synthesis.R`) must invoke `validate_upstream_contracts()` at initialization to ensure all upstream tables comply with required primary keys and columns prior to computing downstream triage calls.

---

## 6. Phase 2 Native-DPI Basal Filtering & Geometry Invariants

* **Dual-Class Stage 1 Box Processing:** Stage 1 YOLOv5 (`LeafPriority.pt`) must process both `Class 0` (`leaf_whole`) and `Class 1` (`leaf_partial`) candidate boxes to capture rosette-attached blades on native-DPI scans (`>1000px`), with cauline pre-filtering (`cy < 0.40 * sheet_h`) and fragment rejection (`max(bh, bw) < int(200 * scale_f)`).
* **Post-Gatekeeper Basal Filtering:** High-DPI basal candidates must strictly enforce `cy_rel >= 0.42`, minimum area floors (`>= 12000 * scale_f^2` for Tier 1, `>= 18000 * scale_f^2` and `cy_rel >= 0.48` for Tier 2), boxy clump rejection (`min(w,h)/max(w,h) > 0.78` when `solidity < 0.84`), sheet-relative area floor (`>= 0.28 * max_sheet_area`), and top-6 cap per voucher ranked by $Q = \text{solidity} \cdot \text{ucs} \cdot \sqrt{\text{area}}$ with automatic cleanup of dropped mask and contour files from disk.
* **Botanical Geometry Conventions:**
  * **Resolution-Adaptive Cross-Sectional Width Profiling:** `_sample_cross_sectional_widths()` and `trim_petiole_tail()` enforce `w_span = max(min(max(h, w), 350), int(0.65 * max(h, w)))` so native-DPI basal blades ($>700\text{ px}$ wide) are never truncated.
  * **Resolution-Scaled Petiole Width Thresholds:** `_disambiguate_base_apex_endpoints()` scales `eff_max_petiole_w` when `axis_len > 450.0`, and `trim_petiole_tail()` scales `max_petiole_width = max(16.0, 0.105 * blade_width_est)` and `cordate_guard_w = max(..., 0.135 * blade_width_est)` when `blade_width_est > 250.0 and min(h, w) > 500`, preserving exact backward compatibility on $400 \times 400$ synthetic unit-test fixtures.
  * **Spatial Rosette Centroid Polarity & Base Homologization:** `02_segment_and_extract.py` computes `rosette_centroids_by_pid` and passes `rosette_centroid` into `classify_and_route_leaf()` and `compute_geometric_metrics()`, while anchoring all contours via `petiole_attachment_pt=p_base` in every `resample_normalized_contour()` call across Tier 1, folded, and Tier 2 reflected leaves. Diagonal-crease fold detection requires `has_diagonal_crease and solidity < 0.72`.

---

## 7. Run-Scoped Directory Structure (`runs/run_<YYYYMMDD_HHMMSS>/`)

Every pipeline execution launched via `./run_pipeline.sh` or `python main.py` (with `--run-id <name>` or `PACKERA_RUN_DIR`) writes all generated artifacts into an isolated, self-contained run directory under `runs/${RUN_ID}` (`RUN_ID="${RUN_ID:-run_$(date +%Y%m%d_%H%M%S)}"`), while preserving shared read-only caches (`data/raw_vouchers/`, `data/raw_images/`, and `data/environmental/`) at the repository root:

```text
Packera_dubia_morphometrics/
├── data/
│   ├── raw_vouchers/                        # Shared read-only herbarium sheet cache (never duplicated)
│   ├── raw_images/                          # Shared read-only raw image cache
│   ├── environmental/                       # Shared read-only WorldClim v2.1 & SoilGrids 250m rasters
│   └── tables/curated_vouchers.csv          # Shared master voucher table (snapshotted into each run)
└── runs/
    ├── latest -> run_20261024_143000        # Symlink pointing to the most recent run directory
    └── run_20261024_143000/                 # Self-contained run directory (${RUN_DIR})
        ├── data/
        │   ├── masks/
        │   │   ├── tier1_intact/            # Direct intact & dissected binary leaf silhouettes
        │   │   └── tier2_reflected/         # Bilaterally reflected half-blade silhouettes
        │   ├── contours/                    # Standardized 120-vertex clockwise (x, y) CSV contours
        │   ├── cropped_patches/
        │   │   ├── basal_leaves_raw/        # Native-DPI RGB crops of extracted basal leaves
        │   │   ├── capitula/                # Cropped reproductive capitula for phenology scoring
        │   │   └── rosettes_dense/          # Tier 4 dense basal rosette crops per plant_individual_id
        │   └── tables/
        │       ├── curated_vouchers.csv
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

* **CLI Execution, Resumption & Dual-Strategy Single-Taxon Modes:**
  * New timestamped run: `./run_pipeline.sh` or `python main.py run-all`
  * Named or resumed run: `./run_pipeline.sh --run-id revision_v1 [--skip-harvest] [--start-step <1..7>] [--force|--overwrite] [--limit <N>] [--sam2-refine]` or `python main.py segment --run-id revision_v1`
  * **Strategy 1 (Focal Taxon + Outgroup Gold Reference Anchors):** `./run_pipeline.sh --run-id dubia_triage --skip-harvest --focal-taxon "Packera dubia" --gold-anchors 30` (processes all focal taxon vouchers plus `N` `Tier_1_Gold` anchors per outgroup species for Steps 04–06 calibration, then filters Step 07's final triage queue and figures strictly to `--focal-taxon`).
  * **Strategy 2 (Pure Intraspecific Delimitation Mode):** `./run_pipeline.sh --run-id dubia_cryptic --skip-harvest --focal-taxon "Packera dubia" --intraspecific-mode` (runs strictly on a single species `n_unique_taxa == 1` with zero outgroup specimens, switching Steps 04–06 from `species_standardized` to label-blind GMM morphotypes `Morphotype_K1..K_n` or geographic/edaphic ecotypes to test for cryptic lineages and granite-outcrop ecotypes).
* **Backward Compatibility Guard:** When any phase script is invoked directly by unit tests without `PACKERA_RUN_DIR` or `--run-id`, it honors explicit CLI path arguments (or default `data/` / `outputs/` paths) so all 253 unit tests pass unchanged.

---

## 8. Publication Figure Catalog

Every pipeline run exports publication-ready figures in both vector `.pdf` and 300-DPI `.png` formats inside `${RUN_DIR}/outputs/figures/` using the colorblind-safe Okabe-Ito palette (`#0072B2` *Packera anonyma*, `#D55E00` *P. dubia*, `#009E73` *P. paupercula*, `#CC79A7` *P. plattensis*) and italicized taxon labels:

| Figure Stem (`.pdf` & `.png`) | Producing Script | Biological / Statistical Hypothesis Visualized |
| :--- | :--- | :--- |
| `fig02a_segmentation_qc_distributions` | `scripts/pipeline/02_segment_and_extract.py` | (A) Lamina solidity vs. Unoccluded Completeness Score (UCS) by routing tier (`Tier_1_Direct`, `Tier_2_Reflected`) against gatekeeper thresholds (`0.72`, `0.85`); (B) Basal leaf yield per voucher & QC failure breakdown. |
| `fig02b_silhouette_montage` | `scripts/pipeline/02_segment_and_extract.py` | 16-tile montage of representative extracted binary basal leaf silhouettes across `Tier_1_Direct` and `Tier_2_Reflected` with homologized petiole base anchor (`p_base`) and midrib axis overlaid. |
| `fig03a_efa_morphospace_pca` | `scripts/morphometrics/03_fourier_extractor.R` | (A) Full 12-harmonic Elliptic Fourier PCA (`PC1` vs. `PC2`) and (B) bilateral-symmetric (`Momocs::rm_asym`) morphospace (`Sym_PC1` vs. `Sym_PC2`) with 95% taxon confidence ellipses. |
| `fig03b_harmonic_power_and_asymmetry` | `scripts/morphometrics/03_fourier_extractor.R` | (A) Cumulative Fourier harmonic power spectrum (`n = 1..12`) and (B) bilateral asymmetry (`pct_asymmetry`) vs. high-frequency perimeter roughness (`harmonic_roughness`) gating contour outliers. |
| `fig04a_gmm_bic_and_uncertainty` | `scripts/morphometrics/04_gmm_morphotools.R` | (A) Unsupervised `mclust` BIC model selection curves across $G = 1..9$ and parsimonious covariance families; (B) Morphospace colored by posterior cluster with point size scaled by taxonomic uncertainty $U_i = 1 - \max(P_{ik})$. |
| `fig04b_cda_biplot_and_allometry` | `scripts/morphometrics/04_gmm_morphotools.R` | (A) Canonical Discriminant Analysis (`CV1` vs. `CV2`) biplot with honest Leave-One-Out Cross-Validation (`loocv_accuracy`) and (B) vegetative-to-reproductive allometric scaling ($\log_{10}$ leaf size vs. $\log_{10}$ involucre scale). |
| `GradCAM_audit_panel` | `scripts/analysis/05_cleanlab_vision_xai.py` | Side-by-side neutral-background foliar patches and DINOv2 ViT-B/14 $16 \times 16$ patch-token attribution heatmaps for Cleanlab-flagged label noise candidates. |
| `fig05b_cleanlab_quality_and_tomentum` | `scripts/analysis/05_cleanlab_vision_xai.py` | (A) Cleanlab label quality ($1 - c_{\text{error}}$) distributions across taxa and determiner authority tiers (`Gold`, `Silver`, `Bronze`); (B) CIE $L^*$ luminance vs. HSV saturation separating adaxial vs. abaxial arachnoid tomentum and `tomentum_taxon_discordance`. |
| `spatial_rf_niche_importance` | `scripts/analysis/06_multimodal_spatial_rf.R` | Spatial Random Forest variable importance across SoilGrids 250m edaphic predictors, WorldClim bioclimatic variables, and phenological anomaly residuals under $0.5^\circ \times 0.5^\circ$ spatial block CV. |
| `fig06b_geographic_map_and_hopkins_phenology` | `scripts/analysis/06_multimodal_spatial_rf.R` | (A) Eastern/central North America geographic voucher distribution colored by taxon and shaped by `consensus_conflict_flag`; (B) Collection DOY vs. Latitude against the centered Hopkins Bioclimatic Law baseline ($+4.0\text{ d/}^\circ\text{N}$). |
| `fig07_multimodal_discordance_summary` | `scripts/analysis/07_triage_dashboard_synthesis.R` | (A) Modality-specific conflict rates across Morphology (CDA/GMM), DINOv2 Vision (Cleanlab), and Ecological Niche (Spatial RF); (B) Consensus curatorial triage priority breakdown across the four *Packera* taxa. |

---

## Workspace Skills Reference

Activate the repository's specialized skills in `.agents/skills/` whenever working on the corresponding domains:

* **`/morphometric-cv-guardrails`** (`.agents/skills/morphometric-cv-guardrails/SKILL.md`): Mandatory whenever modifying scripts in `scripts/vision/`, `scripts/pipeline/`, `scripts/annotation_and_training/`, or `scripts/morphometrics/`.
* **`/pipeline-contract-audit`** (`.agents/skills/pipeline-contract-audit/SKILL.md`): Mandatory before committing code changes or when debugging cross-stage Python/R pipeline handshakes.

