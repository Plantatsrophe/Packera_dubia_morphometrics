---
name: morphometric-cv-guardrails
description: Enforces scientific, mathematical, and computer vision protocols for LeafMachine2, Detectron2 PointRend, OpenCV geometry, and R Momocs/MorphoTools2 morphometrics. Use whenever modifying scripts in scripts/vision/, scripts/pipeline/, scripts/annotation_and_training/, scripts/morphometrics/, or scripts/analysis/.
---

# Morphometric & Computer Vision Guardrails (`UNC-BOT-SOP-2026-04-REV5`)

Activate and follow these scientific, mathematical, and computer vision guardrails whenever creating, refactoring, or debugging code in:
* `scripts/vision/` (`lm2_geometry_utils.py`, `configure_leafmachine2.py`, `capitulum_phenology_classifier.py`, `tune_geometry_parameters.py`)
* `scripts/pipeline/` (`02_segment_and_extract.py`)
* `scripts/annotation_and_training/` (`annotate_with_sam2.py`, `sam2_annotator_utils.py`, `finetune_lm2_pcd.py` — optional standalone utilities)
* `scripts/morphometrics/` (`03_fourier_extractor.R`, `04_gmm_morphotools.R`, `08_discover_taxon_complexes.R`)
* `scripts/analysis/` (`05_cleanlab_vision_xai.py`, `06_multimodal_spatial_rf.R`, `07_triage_dashboard_synthesis.R`)

---

## 1. Two-Stage Stock LeafMachine2 v2.2+ Handshake & Basal Leaf Candidate Filtering

Never run Detectron2 PointRend directly on full-resolution herbarium sheets, as PointRend is trained on cropped single-organ bounding boxes and fails when presented with whole-sheet scale variance and background clutter.

1. **Stage 1 (Full-Sheet Stock LeafMachine2 v2.2+ Component & Ruler Detection):**
   * Run LeafMachine2's stock YOLOv5 Plant Component Detector (`PLANT_LeafPriority` / `LeafPriority.pt`, `minimum_confidence_threshold: 0.65`) and the **LM2 38-class Ruler Detector** on full herbarium sheets to localize organ bounding boxes (`bbox`) and extract empirical scale bars (`pixels_per_mm`, `is_calibrated_mm`).
   * **Stock LeafMachine2 v2.2+ Configuration (`configure_leafmachine2.py`):** Enforce `treat_leaflet_as_leaf_whole: False` (preventing individual lobes of pinnatifid *Packera* blades from being segmented as standalone leaves), `use_efds_for_png_masks: False` (exporting raw pixel-accurate PointRend binary masks without premature LM2 Fourier smoothing), and `do_apply_landmarks: False` (disabling angiosperm venation landmarking). No custom PCD retraining or manual annotation is required by default.
   * **Dual-Class Candidate Box Processing (`Class 0` & `Class 1`):** Stage 1 processes both `Class 0` (`leaf_whole`) and `Class 1` (`leaf_partial`) candidate boxes (`if s1_cls not in (0, 1): continue`) so rosette-attached basal blades on native-DPI (`3753x5634`) scans are captured despite partial overlap or occlusions.
   * **Upper-Sheet Cauline & Tiny Fragment Pre-Filtering:** On high-DPI images (`is_high_dpi`: `max(sheet_h, sheet_w) > 1000`), pre-filter cauline leaf candidates in the upper inflorescence zone (`cy < 0.40 * sheet_h`), tiny fragments (`max(bh, bw) < int(200 * scale_f)`), and square partial fragments above mid-sheet (`s1_cls == 1 and bh == bw and cy < 0.60 * sheet_h`).
   * Keep `detect_ruler_scale_hough()` only as an asynchronous secondary fallback if the 38-class Ruler Detector finds no ruler block.
2. **Stage 2 (Single-Leaf Native-DPI PointRend Segmentation):**
   * For each candidate leaf bounding box from Stage 1, extract a **12%-padded native-DPI single-leaf crop** from the original high-resolution voucher image (or 6% padding for high-DPI scans).
   * Run the stock LeafMachine2 v2.2+ Detectron2 PointRend checkpoint (`leaf_seg_rcnn_v2_1`) **strictly on these padded native-DPI single-leaf crops** (never on full sheets).
   * Project the segmented binary crop mask back into full-sheet coordinates for spatial clustering (`cluster_plant_individuals`) via DBSCAN (`plant_individual_id`), while preserving native-DPI crop resolution for contour extraction.

---

## 2. Label-Blind Classical Geometric Gatekeeping & Boundary-Continuity Lobe Rejection

To prevent taxonomic circularity under the Unified Species Concept (de Queiroz 2007), all geometric quality gatekeeping and tier routing in `classify_and_route_leaf()` must be **100% blind to the specimen's species determination (`taxon` / `scientificName` / `species_raw`)**.

* **Otsu Boundary-Continuity Lobe & Fragment Rejection (`is_truncated_lobe_or_leaflet`, `max_tissue_border_ratio: 0.15`):**
  * Inspects a 2-pixel perimeter band around the candidate bounding box using Otsu's threshold in CIE $L^*a^*b^*$ luminance space. If non-paper plant tissue occupies $> 15\%$ of the bounding-box perimeter (excluding a single narrow petiole exit $\le 12\%$ on one edge), the candidate is rejected as an unattached terminal lobe or severed fragment of a larger pinnatifid leaf.
* **Ruler-Calibrated Physical Area Floor (`min_blade_area_cm2: 0.75`):**
  * When `is_calibrated_mm == True` and `pixels_per_mm > 0`, compute physical blade area $\text{Area}_{\text{cm}^2} = \text{area\_px} / (100 \cdot \text{pixels\_per\_mm}^2)$ and reject fragments below $0.75\text{ cm}^2$.
* **Classical Dimensionless Shape Descriptors (`Tier_1_Direct`):**
  * Intact closed-contour EFA routing (`Tier_1_Direct`) evaluates classical dimensionless morphological descriptors:
    * **Pristine / Lyrate-Entire Tier 1 (`Tier_1_Direct`):** Requires $\text{Solidity} \ge 0.72$, isoperimetric $\text{Circularity} = \frac{4\pi A}{P^2} \ge 0.25$, bounding-box $\text{Extent} = \frac{A}{W \cdot H} \ge 0.35$, and intact boundary continuity ($\text{UCS} \ge 0.85$).
    * **Pinnatifid / Lyrate-Pinnatifid Tier 1 (`Tier 1 (Dissected)`):** Low-solidity leaves in $0.48 \le \text{Solidity} < 0.72$ are routed into `Tier_1_Direct` solely on the basis of intrinsic geometric bilateral sinus topology via `is_botanical_dissection()` (evaluating paired convexity defects, sinus depth ratios $\ge 0.08$, and bilateral lobe symmetry across the midrib).
    * **Bipinnatisect / Deeply Dissected Tier 1 (`Tier 1 (Bipinnatisect)`):** Deeply dissected leaves in $0.18 \le \text{Solidity} < 0.48$ (e.g., *Packera millefolium*) are routed into `Tier_1_Direct` via `is_botanical_dissection(..., bipinnatisect_min_solidity=0.18)` when $\ge 4$ deep convexity-hull sinuses confirm true botanical dissection rather than insect herbivory or tearing.
  * **Forbidden Pattern:** Never use `taxon_has_lyrate_tendency` or any species-name lookup (`Packera paupercula`, `Packera plattensis`, etc.) to bypass or override `is_botanical_dissection()`.
* **Strict Tier 1 Closed-Contour EFA Scope (`Tier_2_Reflected` Excluded from EFA):**
  * Only intact `Tier_1_Direct` contours enter Elliptic Fourier Analysis (`Momocs::efourier`) in Step 03. Folded or torn leaves (`Tier_2_Reflected`) and linear scalar traits are retained only for auxiliary linear measurements and count diagnostics; mirrored hemi-blades and open-curve approximations are **never** pooled into closed-contour EFA harmonic matrices.
* **Post-Gatekeeper High-DPI Basal-Leaf Filtering & Strengthened Cauline/Stem Rejection (`02_segment_and_extract.py`):**
  * **Vertical Position Floor:** Enforce `cy_rel >= 0.42` relative to sheet height, ensuring upper stem foliage is rejected.
  * **Minimum Pixel Area Thresholds:** Require `m_area >= 12000 * scale_f^2` for Tier 1 intact blades; require `m_area >= 18000 * scale_f^2` and `cy_rel >= 0.48` for Tier 2 auxiliary candidates.
  * **Boxy Clump Rejection:** Reject non-foliar clumps exhibiting aspect ratios `min(w, h) / max(w, h) > 0.78` when `solidity < 0.84`.
  * **Strengthened Basal-vs.-Cauline Anatomical Filter (`is_cauline_or_stem_contaminant`):** Rejects upper-stem sessile/clasping cauline bracts (`aspect_ratio > 3.8` and `solidity < 0.78`, or displaced `> 2.2 * med_rosette_radius` above the basal rosette zone) even on herbarium sheets where capitula were not detected (`capitula_count == 0`).
  * **Sheet-Relative Area Floor & Top-6 Cap with Disk Cleanup:** Discard minor fragments with area `_area < 0.28 * max_sheet_area`, retain at most the top 6 basal leaves per voucher ranked by $Q = \text{solidity} \cdot \text{ucs} \cdot \sqrt{\text{area}}$, and unconditionally clean up dropped binary masks, 2D contour CSVs, and RGB patch files from disk.

---

## 3. OpenCV Affine & Contour Math

Every geometric operation in `scripts/vision/lm2_geometry_utils.py` and `scripts/pipeline/02_segment_and_extract.py` must obey these exact coordinate conventions:

1. **Crop with Padding Before Rotating:**
   * Before calling `cv2.warpAffine()`, always crop the binary mask to its tight foreground bounding box plus sufficient zero-border padding (`pad = int(math.ceil(0.5 * math.hypot(h, w)))` via `cv2.copyMakeBorder`) so that affine rotation never clips apices or basal lobes against image canvas boundaries.
2. **Sign of Rotation Angle (`+angle` in OpenCV):**
   * Given the directed midrib vector $\overrightarrow{p_{\text{base}} p_{\text{apex}}} = (\Delta x, \Delta y) = (x_{\text{apex}} - x_{\text{base}},\, y_{\text{apex}} - y_{\text{base}})$ with $\text{angle} = \text{degrees}(\text{atan2}(\Delta y, \Delta x))$ in image coordinates ($y$ pointing downward), OpenCV's `cv2.getRotationMatrix2D(center, theta, 1.0)` rotates **counter-clockwise** in screen space by `theta`.
   * Therefore, always pass **`+angle`** (`cv2.getRotationMatrix2D((cx, cy), +angle, 1.0)`) to rotate the directed vector $\overrightarrow{p_{\text{base}} p_{\text{apex}}}$ onto the positive horizontal $+x$ axis ($0^\circ$).
3. **Disambiguate `p_base` vs. `p_apex` via Resolution-Adaptive Width Tapering & Rosette Centroids:**
   * **Resolution-Adaptive Cross-Sectional Width Profiling (`w_span`):** In `_sample_cross_sectional_widths()` and `trim_petiole_tail()`, compute `w_span = max(min(max(h, w), 350), int(0.65 * max(h, w)))` so native-DPI basal blades ($>700\text{ px}$ wide) are never truncated when profiling orthogonal widths across the terminal $5\%\text{–}15\%$ zones of both extrema.
   * **Resolution-Scaled Petiole Width Thresholds:** In `_disambiguate_base_apex_endpoints()`, scale `eff_max_petiole_w = max(float(max_petiole_width), 0.085 * axis_len)` when `axis_len > 450.0`. In `trim_petiole_tail()`, scale `max_petiole_width = max(16.0, 0.105 * blade_width_est)` and `cordate_guard_w = max(..., 0.135 * blade_width_est)` when `blade_width_est > 250.0 and min(h, w) > 500`, preserving exact backward compatibility on $400 \times 400$ synthetic unit-test fixtures.
   * **Spatial Rosette Centroid Wiring (`rosette_centroids_by_pid` $\rightarrow$ `rosette_centroid`):** In `02_segment_and_extract.py`, compute per-individual DBSCAN rosette centroids (`rosette_centroids_by_pid`) and pass `rosette_centroid` into `classify_and_route_leaf()` and `compute_geometric_metrics()`. The endpoint closer to `rosette_centroid` is favored as `p_base`.
4. **Anatomical Base Homologization (`petiole_attachment_pt=p_base`) & Native-DPI Petiole Junction Trimming (`trim_petiole_tail`):**
   * Pass `petiole_attachment_pt=p_base` into **all** `resample_normalized_contour()` calls in `02_segment_and_extract.py` so Index 0 is always anchored at the anatomical petiole insertion (`p_base`).
   * Resolution-adaptive cross-sectional width sampling identifies the petiole-lamina junction via width gradient inflection, trimming fibrous non-blade petiole tails prior to computing circularity, solidity, and harmonic decomposition.
5. **Isotropic Contour Normalization (`max(span_x, span_y)`):**
   * When exporting normalized contour coordinates (`x_norm`, `y_norm`) in `export_standardized_contour()`, always scale **both** axes isotropically by `max_span = max(span_x, span_y, 1e-6)`. Never scale axes independently, which collapses diagnostic length-to-width aspect ratios.

---

## 4. R Elliptic Fourier Analysis (EFA), Strict Active/Passive Authority Split, Allometry & GMM Morphometrics

When modifying `scripts/morphometrics/03_fourier_extractor.R`, `scripts/morphometrics/04_gmm_morphotools.R`, `scripts/morphometrics/08_discover_taxon_complexes.R`, and `scripts/morphometrics/pub_figures_utils.R`, enforce the following mathematical protocols:

1. **Strict `Tier_1_Direct` Closed-Contour EFA, Inertial Pre-Alignment & $\ge 99\%$ Cumulative Harmonic Power (`03_fourier_extractor.R`):**
   * **Strict Tier 1 Filtering:** Only verified `Tier_1_Direct` intact basal leaf contours enter `Momocs::efourier`. Never pool reflected hemi-blades (`Tier_2_Reflected`) or open-curve approximations into EFA.
   * **Inertial Axis Pre-Alignment (`Momocs::coo_alignxax`) with Suborbicular Guard:** Contours are centered (`Momocs::coo_center`) and aligned along their major inertial axis via `Momocs::coo_alignxax()` whenever blade circularity $\le 0.75$ (suppressing $90^\circ$ inertial axis flipping on suborbicular *P. obovata* / *P. aurea* blades while preserving $\overrightarrow{p_{\text{base}} p_{\text{apex}}}$ orientation via `coo_slidedirection("left")`).
   * **Empirical Harmonic Power Calibration (`calibrate_harmonicpower_efourier`):** Verify that $H = 12$ harmonics capture $\ge 99\%$ cumulative harmonic power across the corpus (`Momocs::calibrate_harmonicpower_efourier`), then run `Momocs::efourier(out_obj, nb.h = 12, norm = FALSE, start = FALSE)`.
   * **Uniform Bilateral Symmetrization (Claude 2008):** Apply `Momocs::rm_asym(ef)` uniformly across all contours so primary morphospace ordinates symmetric lamina architecture ($A_n, D_n$) without taphonomic pressing curvature ($B_n, C_n$ retained as scalar `pct_asymmetry` and `asym_var`).
   * **Minimum Covariance Determinant (MCD) Mahalanobis Outlier Filtering & Unweighted OTU Median Aggregation:** Flag extreme contour extraction anomalies in symmetric harmonic space using robust MCD Mahalanobis distances (`MASS::cov.rob(..., method = "mcd")`, $\chi^2$ cutoff $p < 0.001$ combined with `pct_asymmetry > 0.22 | harmonic_roughness > 0.15`), fit unweighted covariance PCA (`stats::cov.wt` with uniform clean-contour weights) strictly on non-outlier `Tier_1_Direct` blades, and aggregate multi-leaf vouchers to the Specimen / Operational Taxonomic Unit (`(catalogNumber, plant_individual_id)`) via **unweighted coordinate-wise median** harmonics and PCs (`morph_tier_source = "Tier_1_Pure"`).
2. **Strict Active/Passive Determiner Authority Split & Empirical Dawid-Skene Weighting (`04_gmm_morphotools.R`, `06_multimodal_spatial_rf.R`):**
   * **Active Reference Anchor Set (`Tier_1_Gold = 1.0`, `Tier_2_Silver = 1.0` or $w_k = 1 - \text{error\_rate}_k$):** Supervised training in CDA, allometric regression, Cleanlab cross-validation, and Spatial Random Forests is restricted strictly to `Tier_1_Gold` (monographic/type-verified) and `Tier_2_Silver` (regional specialist) `Tier_1_Pure` vouchers after per-taxon $\chi^2_{0.95}$ Mahalanobis outlier pruning (Šlenker et al. 2022). When empirical per-tier leave-one-out discordance rates are available, assign Dawid-Skene precision weights $w_k = 1 - \text{error\_rate}_k$ (`compute_empirical_tier_weights`, Dawid & Skene 1979); otherwise assign uniform active anchor weights ($1.0$).
   * **Passive Projection of Unverified Vouchers (`Tier_3_Bronze = 0.0`):** `Tier_3_Bronze` vouchers receive authority weight `0.0` during model fitting and are projected passively (`predict()`) into the fitted morphospace, CDA canonical space, and GMM clusters so noisy crowd/legacy labels never distort canonical axes or cluster covariance matrices.
3. **Klingenberg (2016) Allometric Slope Homogeneity Testing (`test_and_correct_allometry`):**
   * Restrict allometric calibration to ruler-calibrated specimens (`is_calibrated_mm == TRUE`) in the active anchor set.
   * For each shape PC, first test homogeneity of allometric slopes across taxa via the interaction model `stats::lm(PC ~ log_centsize * species_standardized)` ($F$-test on `log_centsize:species_standardized` at $\alpha = 0.01$, Klingenberg 2016).
   * Apply pooled within-group allometric subtraction (`stats::lm(PC ~ log_centsize + species_standardized)`) **only** when slopes are parallel ($p \ge 0.01$); when slopes diverge significantly ($p < 0.01$), retain raw shape PCs and record `slope_homogeneity_p` and `allometry_correction_mode = "Skipped_Heterogeneous_Slopes"` in `allometry_homogeneity_test.csv`.
4. **Horn's Parallel Analysis, Williams & Titus (1988) Sample-to-Variable Cap & Kass-Raftery GMM (`04_gmm_morphotools.R`):**
   * **Subspace Dimensionality Selection (`select_informative_pcs`):** Select the number of informative PCs $P$ via Horn's Parallel Analysis / Broken-Stick eigenvalue threshold, capped for supervised CDA by the Williams & Titus (1988) rule $P_{\text{cda}} \le \min(P_{\text{horn}}, \lfloor N_{\min} / 3 \rfloor)$ where $N_{\min}$ is the smallest active reference taxon sample size.
   * **Label-Blind GMM (`run_gmm_cluster_analysis`):** Fit `clustvarsel::clustvarsel()` and `mclust::Mclust()` on active `Tier_1_Pure` anchors (`!mahal_outliers` from `MASS::cov.rob`), evaluating evidence via Kass & Raftery (1995) twice-log-Bayes-Factor thresholds ($2\Delta\text{BIC} = \text{BIC}_{\text{alt}} - \text{BIC}_{\text{null}}$, requiring $2\Delta\text{BIC} > 6$ for strong evidence and $> 10$ for very strong evidence of discrete multimodal clusters).
   * **3D Passive CDA (`run_canonical_discriminant_analysis`):** Train on Gold and Silver `Tier_1_Pure` anchors with per-taxon 95% $\chi^2$ Mahalanobis outlier trimming (`qchisq(0.95, df = p)`), uniform class priors (`rep(1/g, g)`), and honest Leave-One-Out Cross-Validation (`loocv_accuracy`), exporting 3D canonical scores (`can1`, `can2`, `can3`, with orthogonalized `PC3`/`PC2` fallback via `orth_against` when $g - 1 < 3$).
   * **Taxon-Adjusted Institutional Batch-Effect Audit (`audit_institutional_batch_effects`):** Fits marginal ANCOVA models `lm(PC1 ~ species_standardized + institutionCode)` and `lm(PC2 ~ species_standardized + institutionCode)` when $\ge 2$ taxa are present, computing partial $R^2$ (`SS_inst / (SS_inst + SS_res)`) for `institutionCode` so biogeographic species turnover across regional herbaria is never conflated with mounting bias.
5. **Step 08 Unsupervised Eastern North American Taxon-Complex Discovery (`08_discover_taxon_complexes.R`):**
   * Computes pairwise morphological ($D_{\text{Mahal}}$), environmental ($1 - D_{\text{Schoener}}$), phenological ($d_{\text{circ}}$), and vision centroid cosine distances across all regional *Packera* species, fuses them into a multi-modal distance matrix $D_{\text{total}}$, and partitions species into candidate complexes via $k$-NN graph modularity (`igraph::cluster_leiden` / `cluster_louvain`) and Ward's hierarchical clustering (`hclust(..., method = "ward.D2")`) with silhouette width validation (`discovered_taxon_complexes.csv`).

---

## 5. Deep Vision, Rhizosphere/Outcrop ENM Stack & Reliability-Weighted Synthesis Guardrails

1. **Taphonomic Indumentum-Blind Macro-Architecture, Foreground Spatial Token Pooling & Cleanlab (`05_cleanlab_vision_xai.py`):**
   * **Indumentum-Blind Gaussian Low-Pass Preprocessing ($\sigma = 1.0\text{ px}$ at $224 \times 224$):** Apply a $\sigma = 1.0\text{ px}$ Gaussian low-pass filter (`apply_indumentum_blind_lowpass`) after resizing crop patches to $224 \times 224$ so high-frequency trichome/floccose hair textures (which abrade taphonomically on older herbarium sheets or vary across adaxial/abaxial mounting) are attenuated while macro-architectural venation, sinus geometry, and petiole winging are preserved.
   * **Foreground-Masked $16 \times 16$ Spatial Patch Pooling (`dinov2_vitb14_reg`):** Load `dinov2_vitb14_reg` (4 register tokens preventing background high-norm artifacts; Darcet et al. 2024), downsample the binary foreground mask to the $16 \times 16$ patch grid (`compute_foreground_patch_mask`, patch threshold $\ge 25\%$ foreground), and average-pool strictly over foreground spatial patch tokens (`forward_features()["x_norm_patchtokens"]`) rather than relying on an unmasked global `[CLS]` token.
   * **Strict Gold/Silver Anchor Training, `StratifiedGroupKFold(groups=institutionCode)` & Temperature/Platt Calibration:** Fit the out-of-fold L2-penalized classifier (`LogisticRegression(C=0.1, class_weight="balanced")`) inside a `Pipeline` (fold-internal `StandardScaler` + `PCA(0.95)`) strictly on `Tier_1_Gold` and `Tier_2_Silver` anchors using `StratifiedGroupKFold` grouped by `institutionCode` (holding out entire herbaria per fold to prevent institutional mounting/background shortcuts), calibrate probabilities via Platt/temperature scaling (`calibrate_oof_probabilities`), and project `Tier_3_Bronze` vouchers passively (`predict_proba`).
   * **Class-Conditional Cleanlab Noise Pruning & Feature-Space OOD Scoring (`vision_ood_score`):** Run `cleanlab.filter.find_label_issues(..., filter_by="prune_by_noise_rate")` so minority species are not disproportionately flagged, and compute Mahalanobis + cosine $k$-NN feature-space distance to the nearest reference taxon centroid (`vision_ood_score`, `vision_is_ood`) so novel or unmodeled morphotypes are flagged as out-of-distribution rather than forced into a high-confidence known class.
   * **Pointing Game Attention Attribution (`foreground_attention_ratio >= 0.65`):** Compute inside-mask vs. total attribution energy (`foreground_attention_ratio`) and flag predictions with `< 65%` foreground plant attribution (`is_background_shortcut`) so background shortcuts cannot drive taxonomic reassignment.
2. **$0\text{–}15\text{ cm}$ Rhizosphere Depth Integration, Uncertainty-Buffered Extraction & Partitioned Niche Equivalency (`06_multimodal_spatial_rf.R` + `edaphic_fine_grain_utils.R`):**
   * **$0\text{–}15\text{ cm}$ Trapezoidal Rhizosphere Depth-Weighting:** Combine $0\text{–}5\text{ cm}$ and $5\text{–}15\text{ cm}$ soil horizons via trapezoidal depth weighting $X_{0\text{–}15\text{ cm}} = \frac{5 X_{0\text{–}5} + 10 X_{5\text{–}15}}{15}$ across SoilGrids 250m and 30m POLARIS/gSSURGO predictors, paired with isometric log-ratio / centered log-ratio (`clr_sand_silt`, `clr_clay_geom`) soil texture transformations (Compositional Data Analysis), coarse fragment volume (`cfvo_0_15cm`), available water capacity (`awc_0_15cm`), depth to bedrock (`fine_bedrock_depth_cm`), 50m-buffered SSURGO rock outcrop percentage (`ssurgo_rock_outcrop_pct`), 30m 3DEP micro-topography (`topo_tpi_30m`, `topo_twi_30m`), and 6-class USGS SGMC bedrock lithology (`Felsic_Granitic_Gneissic`, `Mafic_Diabase_Gabbro`, `Ultramafic_Serpentine`, `Carbonate_Limestone_Dolostone`, `Siliciclastic_Sandstone_Shale`, `Unconsolidated_Coastal_Alluvial`).
   * **Coordinate-Uncertainty Buffered Extraction (`100–1000 m`) & Collinearity Pruning:** Extract focal median raster values within `r_buf = clamp(coordinateUncertaintyInMeters, 100, 1000)` m (downweighting sheets with $> 5000\text{ m}$ uncertainty and excluding $> 10000\text{ m}$ from anchor training), and prune collinear predictors at $|r| < 0.70$ and $\text{VIF} < 5.0$ (`prune_collinear_predictors`) before fitting `ranger::ranger` under $0.5^\circ \times 0.5^\circ$ Spatial Block Cross-Validation (`spatial_block_cv_acc`).
   * **Partitioned Climatic vs. Edaphic-Geologic Niche Equivalency & Background Similarity (`ecospat` / Warren et al. 2008):** Compute Schoener's $D$ and Hellinger's $I$ separately in Climatic PCA space ($D_{\text{clim}}, I_{\text{clim}}$) and Edaphic-Geologic PCA space ($D_{\text{edaph}}, I_{\text{edaph}}$) alongside 100-replicate permutation $p$-values (`edaphic_niche_equivalency.csv`).
3. **Reliability-Weighted Consensus, Hopkins-Adjusted Circular Anthesis Phenology & 5-Stage Hierarchical Decision Matrix (`07_triage_dashboard_synthesis.R`):**
   * **Out-of-Fold Cohen's $\kappa$ Reliability Weighting:** Weight each diagnostic modality $m \in \{\text{Morph}, \text{Vis}, \text{Edaph}, \text{Pheno}\}$ by its empirical out-of-fold validation Cohen's $\kappa_m$ ($w_m = \max(0.05, \kappa_m)$), dynamically renormalizing weights over the subset of observed modalities $\mathcal{O}_i$ on each specimen ($S_i = \frac{\sum_{m \in \mathcal{O}_i} w_m d_{i,m}}{\sum_{m \in \mathcal{O}_i} w_m}$) so non-georeferenced vouchers are evaluated without missing-data penalties (`98.4%` voucher completion).
   * **Hopkins-Adjusted Circular Anthesis Phenology:** Restrict phenological divergence testing strictly to specimens verified in active anthesis (`in_anthesis == TRUE`, Stage 2 flowering capitula), regressing circular flowering angle $(\sin\theta_i, \cos\theta_i)$ on latitude, elevation, and MAT (`atan2` circular residuals) before computing angular deviation $|\Delta\text{DOY}_{\text{Hopkins}}|$.
   * **5-Stage Hierarchical Decision Matrix:**
     1. **Stage 1 (`Misidentification_Reassign`):** Cross-modal alternative-species agreement (`both_alt_agree`, `morph_edaph_alt_agree`, `vis_edaph_alt_agree`, `is_severe_misid`) is evaluated **FIRST**, before hybrid or cryptic states, so sympatric misidentifications are routed directly to `Misidentification_Reassign`.
     2. **Stage 2 (`Hybrid_Introgressant`):** Intermediate CDA posteriors ($0.35 \le P(\text{taxon}_1) \le 0.65$ and $P(\text{taxon}_2) \ge 0.25$) in verified sympatry where no single alternative taxon achieves cross-modal consensus.
     3. **Stage 3 (`Putative_Cryptic_Species` — Integration by Congruence; Padial et al. 2010):** Requires **`morph_tier_source == "Tier_1_Pure"` (`n_tier1_leaves >= 1`)**, `morphometric_outlier_flag == TRUE`, high OOD / low known-taxon posterior across both morphology (`cda_posterior_prob < 0.70`) and vision (`vision_c_error < 0.85`), placement in a discrete GMM cluster supported by $2\Delta\text{BIC} > 6$, AND corroboration by at least one independent non-morphological axis (significant edaphic niche divergence $p < 0.05$, Hopkins-adjusted allochrony $|\Delta\text{DOY}| \ge 14\text{ d}$, or distinct geographic cluster).
     4. **Stage 4 (`Subspecies_Ecotype` vs. `Ecophenotypic_Plasticity`):** Intraspecific morphs on specialized substrates (e.g., granite flatrocks or serpentine barrens) with $2\Delta\text{BIC} \le 6$ route to `Subspecies_Ecotype`; continuous size/allometric variation with overlapping niches routes to `Ecophenotypic_Plasticity`.
     5. **Stage 5 (`Confirmed_Species`):** Concordant baseline vouchers (`Retain_Original_Label`).

---

## Verification Checklist Before Completing Changes

1. Confirm no changes were made inside `LeafMachine2/`.
2. Run unit tests covering geometry, orientation, Fourier decomposition, and contracts:
   ```bash
   python3 -m unittest discover -s scripts/tests -p "test_*.py" -v
   ```
3. Verify that `test_code_hygiene.py` confirms 0 `_archive` references across active files.
