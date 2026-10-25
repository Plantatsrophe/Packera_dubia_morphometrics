---
name: morphometric-cv-guardrails
description: Enforces scientific, mathematical, and computer vision protocols for LeafMachine2, Detectron2 PointRend, OpenCV geometry, and R Momocs/MorphoTools2 morphometrics. Use whenever modifying scripts in scripts/vision/, scripts/pipeline/, scripts/annotation_and_training/, scripts/morphometrics/, or scripts/analysis/.
---

# Morphometric & Computer Vision Guardrails (`UNC-BOT-SOP-2026-04-REV4`)

Activate and follow these scientific, mathematical, and computer vision guardrails whenever creating, refactoring, or debugging code in:
* `scripts/vision/` (`lm2_geometry_utils.py`, `configure_leafmachine2.py`, `capitulum_phenology_classifier.py`, `tune_geometry_parameters.py`)
* `scripts/pipeline/` (`02_segment_and_extract.py`)
* `scripts/annotation_and_training/` (`annotate_with_sam2.py`, `sam2_annotator_utils.py`, `finetune_lm2_pcd.py`)
* `scripts/morphometrics/` (`03_fourier_extractor.R`, `04_gmm_morphotools.R`)
* `scripts/analysis/` (`05_cleanlab_vision_xai.py`, `06_multimodal_spatial_rf.R`, `07_triage_dashboard_synthesis.R`)

---

## 1. Two-Stage LeafMachine2 Handshake & Basal Leaf Candidate Filtering

Never run Detectron2 PointRend directly on full-resolution herbarium sheets, as PointRend is trained on cropped single-organ bounding boxes and fails when presented with whole-sheet scale variance and background clutter.

1. **Stage 1 (Full-Sheet Component & Ruler Detection):**
   * Run LeafMachine2's stock YOLOv5 Plant Component Detector (`LeafPriority.pt`) and the **LM2 38-class Ruler Detector** on full herbarium sheets to localize organ bounding boxes (`bbox`) and extract empirical scale bars (`pixels_per_mm`, `is_calibrated_mm`).
   * **Dual-Class Candidate Box Processing (`Class 0` & `Class 1`):** Stage 1 processes both `Class 0` (`leaf_whole`) and `Class 1` (`leaf_partial`) candidate boxes (`if s1_cls not in (0, 1): continue`) so rosette-attached basal blades on native-DPI (`3753x5634`) scans are captured despite partial overlap or occlusions.
   * **Upper-Sheet Cauline & Tiny Fragment Pre-Filtering:** On high-DPI images (`is_high_dpi`: `max(sheet_h, sheet_w) > 1000`), pre-filter cauline leaf candidates in the upper inflorescence zone (`cy < 0.40 * sheet_h`), tiny fragments (`max(bh, bw) < int(200 * scale_f)`), and square partial fragments above mid-sheet (`s1_cls == 1 and bh == bw and cy < 0.60 * sheet_h`).
   * Keep `detect_ruler_scale_hough()` only as an asynchronous secondary fallback if the 38-class Ruler Detector finds no ruler block.
2. **Stage 2 (Single-Leaf Native-DPI PointRend Segmentation):**
   * For each candidate leaf bounding box from Stage 1, extract a **12%-padded native-DPI single-leaf crop** from the original high-resolution voucher image (or 6% padding for high-DPI scans).
   * Run the fine-tuned Detectron2 PointRend checkpoint (`leaf_seg_rcnn_v2_1` / `models/lm2_packera_pcd_finetuned.pth`) **strictly on these padded native-DPI single-leaf crops** (never on full sheets).
   * Project the segmented binary crop mask back into full-sheet coordinates for spatial clustering (`cluster_plant_individuals`) via DBSCAN (`plant_individual_id`), while preserving native-DPI crop resolution for contour extraction.

---

## 2. Label-Blind Gatekeeping & Post-Gatekeeper Basal Filtering

To prevent taxonomic circularity under the Unified Species Concept (de Queiroz 2007), all geometric quality gatekeeping and tier routing in `classify_and_route_leaf()` must be **100% blind to the specimen's species determination (`taxon` / `scientificName` / `species_raw`)**.

* **Pristine Tier 1 (`Tier_1_Direct`):** Requires $\text{UCS} \ge 0.85$ and $\text{Solidity} \ge 0.72$.
* **Dissected / Lyrate-Pinnatifid Tier 1 (`Tier 1 (Dissected)`):**
  * Low-solidity leaves in the range $0.50 \le \text{Solidity} < 0.72$ must be routed into `Tier 1 (Dissected)` **solely** on the basis of intrinsic geometric bilateral sinus topology via `is_botanical_dissection()` (evaluating paired convexity defects, sinus depth ratios $\ge 0.08$, and bilateral lobe symmetry across the midrib).
  * **Forbidden Pattern:** Never use `taxon_has_lyrate_tendency` or any species-name lookup (`Packera paupercula`, `Packera plattensis`, etc.) to bypass or override `is_botanical_dissection()`.
* **Post-Gatekeeper High-DPI Basal-Leaf Filtering (`02_segment_and_extract.py`):**
  * **Vertical Position Floor:** Enforce `cy_rel >= 0.42` relative to sheet height, ensuring upper stem foliage is rejected.
  * **Minimum Pixel Area Thresholds:** Require `m_area >= 12000 * scale_f^2` for Tier 1 pristine blades; require `m_area >= 18000 * scale_f^2` and `cy_rel >= 0.48` for Tier 2 reflected candidates.
  * **Boxy Clump Rejection:** Reject non-foliar clumps exhibiting aspect ratios `min(w, h) / max(w, h) > 0.78` when `solidity < 0.84`.
  * **Sheet-Relative Area Floor:** Discard minor fragments with area `_area < 0.28 * max_sheet_area`.
  * **Top-6 Basal Leaves Cap & Disk Cleanup:** Retain at most the top 6 basal leaves per voucher, ranked by composite quality $Q = \text{solidity} \cdot \text{ucs} \cdot \sqrt{\text{area}}$. Unconditionally clean up dropped binary masks, 2D contour CSVs, and RGB patch files from disk to prevent orphaned artifacts.

---

## 3. OpenCV Affine & Contour Math

Every geometric operation in `scripts/vision/lm2_geometry_utils.py` and `scripts/pipeline/02_segment_and_extract.py` must obey these exact coordinate conventions:

1. **Crop with Padding Before Rotating:**
   * Before calling `cv2.warpAffine()`, always crop the binary mask to its tight foreground bounding box plus sufficient zero-border padding (e.g., `pad = int(math.ceil(0.5 * math.hypot(h, w)))` via `cv2.copyMakeBorder`) so that affine rotation never clips apices, basal lobes, or reflected lamina margins against image canvas boundaries.
2. **Sign of Rotation Angle (`+angle` in OpenCV):**
   * Given the directed midrib vector $\overrightarrow{p_{\text{base}} p_{\text{apex}}} = (\Delta x, \Delta y) = (x_{\text{apex}} - x_{\text{base}},\, y_{\text{apex}} - y_{\text{base}})$ with $\text{angle} = \text{degrees}(\text{atan2}(\Delta y, \Delta x))$ in image coordinates ($y$ pointing downward), OpenCV's `cv2.getRotationMatrix2D(center, theta, 1.0)` rotates **counter-clockwise** in screen space by `theta`.
   * Therefore, always pass **`+angle`** (`cv2.getRotationMatrix2D((cx, cy), +angle, 1.0)`) to rotate the directed vector $\overrightarrow{p_{\text{base}} p_{\text{apex}}}$ onto the positive horizontal $+x$ axis ($0^\circ$). Passing `-angle` doubles the inclination and breaks horizontal alignment.
3. **Disambiguate `p_base` vs. `p_apex` via Resolution-Adaptive Width Tapering & Rosette Centroids:**
   * Principal axis projection (`cv2.fitLine`) identifies the two longitudinal extrema along the midrib but is sign-ambiguous ($180^\circ$ arbitrary).
   * **Resolution-Adaptive Cross-Sectional Width Profiling (`w_span`):** In `_sample_cross_sectional_widths()` and `trim_petiole_tail()`, compute `w_span = max(min(max(h, w), 350), int(0.65 * max(h, w)))` so native-DPI basal blades ($>700\text{ px}$ wide) are never truncated when profiling orthogonal widths across the terminal $5\%\text{–}15\%$ zones of both extrema.
   * **Resolution-Scaled Petiole Width Thresholds:** In `_disambiguate_base_apex_endpoints()`, scale `eff_max_petiole_w = max(float(max_petiole_width), 0.085 * axis_len)` when `axis_len > 450.0`. In `trim_petiole_tail()`, scale `max_petiole_width = max(16.0, 0.105 * blade_width_est)` and `cordate_guard_w = max(..., 0.135 * blade_width_est)` when `blade_width_est > 250.0 and min(h, w) > 500`, preserving exact backward compatibility on $400 \times 400$ synthetic unit-test fixtures.
   * **Spatial Rosette Centroid Wiring (`rosette_centroids_by_pid` $\rightarrow$ `rosette_centroid`):** In `02_segment_and_extract.py`, compute per-individual DBSCAN rosette centroids (`rosette_centroids_by_pid`) and pass `rosette_centroid` into `classify_and_route_leaf()` and `compute_geometric_metrics()`. The endpoint closer to `rosette_centroid` is favored as `p_base`, resolving polarity on rounded or truncated blades.
4. **Anatomical Base Homologization (`petiole_attachment_pt=p_base`) & Native-DPI Petiole Junction Trimming (`trim_petiole_tail`):**
   * Pass `petiole_attachment_pt=p_base` into **all** `resample_normalized_contour()` calls across Tier 1 intact, folded, and Tier 2 reflected leaves in `02_segment_and_extract.py` so Index 0 is always anchored at the anatomical petiole insertion (`p_base`).
   * Resolution-adaptive cross-sectional width sampling identifies the petiole-lamina junction via width gradient inflection, trimming fibrous non-blade petiole tails prior to computing circularity, solidity, and harmonic decomposition.
5. **Diagonal-Crease Fold Detection Logic (`detect_folded_leaf`):**
   * Diagonal crease detections must only trigger reflection routing when accompanied by compromised blade solidity (`has_diagonal_crease and solidity < 0.72`), preventing minor surface venation lines from falsely routing intact Tier 1 leaves to Tier 2 reflection.
6. **Isotropic Contour Normalization (`max(span_x, span_y)`):**
   * When exporting normalized contour coordinates (`x_norm`, `y_norm`) in `export_standardized_contour()`, always scale **both** axes isotropically by `max_span = max(span_x, span_y, 1e-6)`. Never scale axes independently, which collapses diagnostic length-to-width aspect ratios.

---

## 4. R Elliptic Fourier Analysis (EFA), Authority Weighting (`3.0 / 1.5 / 0.5`), Allometry & GMM Morphometrics

When modifying `scripts/morphometrics/03_fourier_extractor.R` and `scripts/morphometrics/04_gmm_morphotools.R`, enforce the following mathematical protocols and the canonical determiner authority weighting hierarchy (`authority_weight`: `Tier_1_Gold = 3.0`, `Tier_2_Silver = 1.5`, `Tier_3_Bronze = 0.5`):

1. **Pre-Aligned EFA without First-Harmonic $180^\circ$ Flips (`norm = FALSE, start = FALSE`):**
   * Contours exported by Phase 2 are already homologized to start at `p_base`, aligned along $\overrightarrow{p_{\text{base}} p_{\text{apex}}}$, centered, and scaled isotropically. Run:
     ```r
     Momocs::efourier(out_obj, nb.h = opts$harmonics, norm = FALSE, start = FALSE)
     ```
2. **Bilateral-Symmetric Harmonic Decomposition & Authority-Weighted PCA (`03_fourier_extractor.R`):**
   * Merge voucher metadata (`curated_vouchers.csv`) **prior to PCA** so determiner authority weights (`tier_wt`: `3.0 / 1.5 / 0.5`) are available during eigendecomposition.
   * Decompose the Fourier matrix into continuous bilateral-symmetric harmonics ($A_n, D_n$) via `Momocs::rm_asym(ef)`, and compute both `Sym_PC1`..`Sym_PC5` and primary `PC1`..`PC5` using weighted covariance eigendecomposition (`stats::cov.wt(..., wt = tier_wt * leaf_quality_weight)`).
   * Validate Tier 1 pristine vs. Tier 2 reflected parity across `assigned_tier` via randomized residual permutation PERMANOVA (`RRPP::lm.rrpp()` or `vegan::adonis2()`).
3. **Segment-Level Quality Defenses in EFA (`03_fourier_extractor.R`):**
   * **Bilateral Asymmetry Ratio (`pct_asymmetry`):** Ratio of asymmetric harmonic power ($B_n, C_n$) to total power ($A_n, B_n, C_n, D_n$).
   * **High-Frequency Harmonic Roughness (`harmonic_roughness`):** Ratio of high-frequency power ($n \ge 6$) to low-frequency power ($n = 1..3$).
   * **Contour Outlier Gating (`is_contour_outlier`):** Flags outlines with `pct_asymmetry > 0.22 | harmonic_roughness > 0.15`. Outliers are strictly excluded before fitting PCA axes (`clean_pca_mask`).
   * **Leaf Quality Weight (`leaf_quality_weight`):** Continuous weighting $W = \max(0.05, \text{solidity} \cdot \text{ucs} \cdot (1 - \text{pct\_asymmetry})^2 \cdot (0.10 \text{ if outlier else } 1.0))$ exported in `fourier_pcs.csv`.
4. **Within-Plant Quality-Weighted Medoid Aggregation (`04_gmm_morphotools.R`):**
   * Multi-leaf vouchers are aggregated to `(catalogNumber, plant_individual_id)` using `leaf_quality_weight`-weighted medoid-filtered averages, preserving true individual centroids.
5. **Anchor-Fitted GMM & Passive CDA with Mahalanobis Outlier Trimming (`04_gmm_morphotools.R`):**
   * **Label-Blind GMM (`run_gmm_cluster_analysis`):** Fit `clustvarsel::clustvarsel()` and `mclust::Mclust()` exclusively on high-reliability `Tier_1_Gold` + `Tier_2_Silver` anchors (`!mahal_outliers` from `MASS::cov.rob`) with parsimonious covariance restriction (`EII`, `VII`, `EEI`, `VEI`, `EVI`, `VVI` when `n_samples / n_pcs < 5`), then project `Tier_3_Bronze` vouchers onto `gmm_fit` via `predict(gmm_fit, newdata = X_pcs)` so noisy sheets do not bloat cluster covariances or inflate `gmm_uncertainty`.
   * **Passive CDA (`run_canonical_discriminant_analysis`):** Train on `Tier_1_Gold` ($w = 3.0$) and `Tier_2_Silver` ($w = 1.5$) anchors with per-taxon 95% chi-square Mahalanobis outlier trimming (`qchisq(0.95, df = p)`), uniform class priors (`rep(1/g, g)`), and honest Leave-One-Out Cross-Validation (`loocv_accuracy`), exporting `cda_predicted_taxon`, `cda_posterior_prob`, `cda_second_posterior_prob`, and `cda_posterior_margin`.
   * **Single-Taxon Intraspecific Mode (`--intraspecific-mode`):** When `n_unique_taxa == 1`, switch Steps 04–06 from `species_standardized` to label-blind GMM morphotypes (`Morphotype_K1..K_n` or geographic/edaphic ecotypes).
6. **Authority-Weighted Pooled Within-Species Allometric Regression (`test_and_correct_allometry`):**
   * Fit allometry on calibrated millimeter Centroid Size ($\text{CS}_{\text{mm}}$) via `stats::lm(Y ~ log_centsize + scientificName, weights = authority_weight)` and subtract size effects using common slope $\hat{\beta}_{\text{within}}$, preserving species intercepts $\alpha_{\text{species}}$.

---

## 5. Deep Vision, Spatial Macroecology & Synthesis Guardrails

1. **DINOv2 Embeddings, Authority Weights & Confident Learning (`05_cleanlab_vision_xai.py`):**
   * L2 cosine normalization of 768-d DINOv2 embeddings prior to fold-internal `StandardScaler`.
   * Apply determiner authority weights (`Tier_1_Gold = 3.0`, `Tier_2_Silver = 1.5`, `Tier_3_Bronze = 0.5`) combined with `class_weight="balanced"` in out-of-fold `LogisticRegression`.
   * Export `tomentum_taxon_discordance` alongside `is_label_corrupted` when abaxial leaf surfaces on glabrate taxa (*P. anonyma*, *P. paupercula*) are predicted as tomentose taxa (*P. dubia*, *P. plattensis*) due to mounting orientation variants.
2. **Authority-Weighted Spatial RF, Centered Hopkins Baseline & Resilient SSURGO (`06_multimodal_spatial_rf.R`):**
   * Use authority-weighted circular phenology and edaphic envelopes in `evaluate_multimodal_concordance()`, center Hopkins' Bioclimatic Law baseline at `pred_at_min_lat + 4.0 * (decimalLatitude - min_lat)`, and pass `case.weights = case_wts` (`3.0 / 1.5 / 0.5`) into `ranger::ranger()`.
   * Resilient batched `soilDB::SDA_spatialQuery()` chunking with per-batch schema normalization (`dplyr::bind_rows`) and SoilGrids 250m fallback to withstand USDA SDA network timeouts and HTML error frames without crashing.
   * Spatially blocked cross-validation ($0.5^\circ \times 0.5^\circ$ grid blocks, `blockCV::cv_spatial()`, exporting `spatial_block_cv_acc`) to mitigate spatial autocorrelation inflation.
3. **Non-Georeferenced Voucher Coalescing (`98.4%` Completion) & Misidentification vs. Hybrid vs. Cryptic Disambiguation (`06`–`07`):**
   * `07_triage_dashboard_synthesis.R` invokes `validate_upstream_contracts()`, joins `multimodal_conflict_flags.csv` (georeferenced subset), and coalesces non-georeferenced vouchers directly from `morphometric_flags.csv` (Step 04) and `label_noise_audit.csv` (Step 05) using `unique(c("catalogNumber", ...))`, achieving **98.4% voucher completion (`5,117 / 5,201` vouchers in full production)** rather than dropping non-georeferenced herbarium sheets.
   * **Pure Misidentification vs. Hybrid Swarm:** When both Morphology (`cda_predicted_taxon`) and Vision (`vision_predicted_label`) agree on the same alternate species (`both_alt_agree`), or when either stream shows high-confidence conflict (`is_severe_misid`), Step 06 classifies the voucher as `"Severe_Triple_Stream_Conflict"` and Step 07 routes it to `Misidentification_Severe` $\rightarrow$ `Misidentification_Reassign` **before** evaluating `Hybrid_Intergrade_Swarm`.
   * **Reassign vs. Putative Cryptic Species (`is_crypt`):** Within `Misidentification_Severe`, ordinary misidentifications where both streams agree on a known alternate species or where either stream assigns high confidence to a known congener (`cda_posterior_prob >= 0.70` or `vision_c_error >= 0.85`) are classified as `Misidentification_Reassign`. `Putative_Cryptic_Species` is reserved strictly for true morphometric outliers (`morphometric_outlier_flag == TRUE`) lacking high-confidence known-taxon assignment that occupy peripheral GMM clusters (`gmm_cluster > 4` with $|\Delta\text{DOY}| \ge 14.0\text{ d}$) or exhibit extreme allochronic flowering divergence ($|\Delta\text{DOY}| \ge 21.0\text{ d}$ with `cda_posterior_prob < 0.55`).
4. **Run-Scoped Execution & Publication Figure Suite (`02`–`07` & `pub_figures_utils.R`):**
   * When `--run-id <name>` or `PACKERA_RUN_DIR` is active, write all masks (`tier1_intact/`, `tier2_reflected/`), contours, cropped patches, CSV manifests (`fourier_pcs.csv`, `symmetric_fourier_pcs.csv`, etc.), reports (`final_curation_triage.csv`, `run_metadata.json`), logs (`pipeline.log`), and figures into `runs/<run_id>/` while preserving backward-compatible CLI paths when invoked directly by unit tests.
   * Every phase (`02`–`07`) must export publication-ready figures in both vector `.pdf` and 300-DPI `.png` formats (`fig02a`, `fig02b`, `fig03a`, `fig03b`, `fig04a`, `fig04b`, `GradCAM_audit_panel`, `fig05b`, `spatial_rf_niche_importance`, `fig06b`, `fig07`) using colorblind-safe Okabe-Ito palettes (`#0072B2`, `#D55E00`, `#009E73`, `#CC79A7`), italicized *Packera* taxon labels, and atomic write-then-replace persistence.

---

## Verification Checklist Before Completing Changes

1. Confirm no changes were made inside `LeafMachine2/`.
2. Run unit tests covering geometry, orientation, Fourier decomposition, and contracts:
   ```bash
   .venv/bin/python3 -m unittest discover -s scripts/tests -p "test_*.py" -v
   ```
3. Verify that `test_code_hygiene.py` confirms 0 `_archive` references across active files.
