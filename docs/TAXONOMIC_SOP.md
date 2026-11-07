# Standard Operating Procedure (SOP)

## Automated High-Throughput Vegetative Phenotyping and Herbarium Triage Pipeline for North American *Packera*
### Operationalizing Morphological Diagnosability under the Unified Species Concept (*Packera dubia* Complex Baseline)

- **Protocol Identifier:** `UNC-BOT-SOP-2026-04-REV5`
- **Principal Investigator:** J. Brandon Fuller (PhD Candidate, Department of Biology, UNC-CH)
- **Faculty Advisor:** Dr. Alan S. Weakley (Director, NCU Herbarium; UNC Biology)
- **Target Taxa:** ***Packera dubia*** (Spreng.) Trock & Mabb. (Asteraceae: Senecioneae) as primary baseline, evaluated dynamically against all regional and continent-wide North American *Packera* Á. Löve & D. Löve congeners without hardcoded species groupings.

---

## 1. Executive Summary & Theoretical Framework

This Standard Operating Procedure establishes the protocol for an **Automated High-Throughput Vegetative Phenotyping and Herbarium Triage Pipeline** applied to North American ***Packera* Á. Löve & D. Löve** (Asteraceae: Senecioneae), utilizing the taxonomically recalcitrant ***Packera dubia* complex** as primary dissertation baseline.

The pipeline operationalizes morphological diagnosability within the **Unified Species Concept (USC; de Queiroz 2007)** via **Integration by Congruence (Padial et al. 2010)**. Under the USC, species are separately evolving metapopulation lineages; morphological diagnosability, reproductive isolation, ecological divergence, and reciprocal monophyly serve as secondary operational criteria. Automated vegetative morphometrics serves as a rigorous operational screen evaluated alongside five independent lines of evolutionary evidence:

1. **Reproductive Macro-Morphology & Anthesis:** Quantitative capitulum aspect ratios, phyllary counts, floret dimensions, and 3-tier optical anthesis staging.
2. **Cytotaxonomic Cytotypes:** Chromosome numbers and polyploid series ($2n = 44, 46, 88$; Kowal 1975).
3. **Reduced-Representation nextRAD Phylogenomics:** High-density SNP datasets resolving reticulate evolution, hybrid swarms, and lineage boundaries.
4. **$0\text{–}15\text{ cm}$ Rhizosphere & USGS SGMC Bedrock Stack:** Cached in SQLite `environmental_cache` + `environmental_rtree`—integrating $0\text{–}15\text{ cm}$ depth-weighting (`(5 * v_0_5 + 10 * v_5_15) / 15`), SSURGO map units (`ssurgo_rock_outcrop_pct`, `ssurgo_bedrock_depth_cm`), 30 m POLARIS, 100 m SOLUS100, 6-class USGS SGMC lithology (`sgmc_lithology_class`), 10 m – 30 m USGS 3DEP topography, and Climatic ($D_{\text{clim}}$) vs. Edaphic ($D_{\text{edaph}}$) niche tests (Warren et al. 2008).
5. **Central Relational WAL Architecture (`data/packera_pipeline.sqlite`):** Transacts against a centralized SQLite WAL database housing Tier A shared tables (`taxon_synonyms`, `vouchers`, `plant_individuals`, `extracted_leaves`, `voucher_reproductive_phenology`, `dinov2_feature_cache`, `environmental_cache` + `rtree`) and Tier B run-scoped tables (`pipeline_runs`, `run_active_cohort`, `leaf_efa_harmonics`, `otu_morphometrics`, `vision_xai_audit`, `spatial_rf_predictions`, `curation_triage`, and master view `v_integrative_triage_master`).

---

## 2. Phase 1 Protocol: Ingestion, Standardization & US EPA Level III Ecoregions (`01_voucher_harvester.py`)

* **Darwin Core Schema Contract & Relational UPSERT (`vouchers`):**
  Harvests occurrences into SQLite `vouchers` (PK: `catalogNumber`) via atomic UPSERT while resolving synonyms through `taxon_synonyms`, exporting synchronized CSV `data/tables/curated_vouchers.csv` (26 DwC columns including `catalogNumber`, `scientificName`, `determiner_tier`, `type_status`, coordinates, ecoregion, duplicate/conflict flags).
* **Dynamic Taxon Standardization & GBIF Synonym Expansion:**
  Standardizes binomials and basionyms dynamically via `expand_taxon_synonyms_gbif()`. Running with no taxon flags (`target_species: []`) queries all *Packera* species across the genus (`discover_genus_packera_species_gbif()`). Passing `--taxon` or `--focal-taxon` restricts analysis strictly to specified species with zero hardcoded defaults or injected outgroups.
* **Configurable Geographic Scope:**
  Supports `--geographic-scope {eastern_na,north_america,global}` (default: `eastern_na`, enforcing `decimalLongitude` bounds `[-106.65, -50.0]` during GBIF queries with `typeStatus="*"` rescue query for historical types; `north_america` and `global` include western occurrences). When `exclude_western` is active (default), excludes Western US states and Western/Northern Canadian provinces (`BC, AB, SK, YT, NT, NU`) and requires $\ge 3$ georeferenced records east of $-106.65^\circ\text{W}$. Supported by `--countries` (default `["US"]`).
* **Pre-Download Deduplication & WebP Storage:**
  `VoucherHarvester.harvest()` deduplicates exsiccatae (`primary_duplicates_only=True`) *prior* to image download, retaining primary sheets (`is_primary_duplicate == TRUE`) while unconditionally preserving verified nomenclatural types. Full-sheet imagery in `data/raw_vouchers/` and crops in `data/cropped_patches/` default to quality-85 `.webp` compression (`method=4` via PIL), providing $\sim 70\text{–}80\%$ disk reduction relative to JPEG/PNG with zero edge degradation.
* **Official US EPA Level III Ecoregion Point-in-Polygon Delineation:**
  Assigns `regional_group` strictly via spatial intersection against official **US EPA Level III Ecoregions** (Omernik & Griffith 2014; `data/environmental/epa_ecoregions_l3/`):
  $$\text{regional\_group}_i = \text{Ecoregion}_k \iff (x_i, y_i) \in \mathcal{P}_k$$
* **Ungeoreferenced Specimen Fallback Policy:**
  When coordinates are missing ($y_i = \text{NA}$), evaluate US Census county boundary polygon $\mathcal{C}_i$:
  $$\text{regional\_group}_i = \begin{cases} \text{Ecoregion}_k & \text{if } \frac{\text{Area}(\mathcal{C}_i \cap \mathcal{P}_k)}{\text{Area}(\mathcal{C}_i)} \ge 0.95 \\ \text{NA} & \text{otherwise} \end{cases}$$
* **Nomenclatural Type Hierarchy & Priority Engine:**
  - `PRIMARY_TYPE_STATUSES = {"HOLOTYPE", "LECTOTYPE", "NEOTYPE", "EPITYPE", "ISOTYPE", "ISOLECTOTYPE", "ISONEOTYPE"}` (`Tier_0_Primary_Type`, rank 2).
  - `SECONDARY_TYPE_STATUSES = {"SYNTYPE", "ISOSYNTYPE", "PARATYPE", "ISOPARATYPE", "TYPE", "COTYPE"}` (`Tier_0_Secondary_Type`, rank 1).
  - `classify_nomenclatural_type_subrank()` resolves fine-grained ICN priority (`2.2` primary, `2.0` isotype, `1.5` syntype, `1.0` paratype; Turland et al. 2018). Verified types are exempt from uncertainty/date filters and unconditionally retained during spatial thinning and outgroup cascades.
* **3-Layer Dynamic Authority Engine (`parse_determiner_tier`):**
  - **Layer 1 (Monographer Heuristics & Determiner Isolation):** `SPECIALIST_PATTERNS` matching is restricted strictly to `identifiedBy` and `identificationHistory` (only falling back to `recordedBy` when `identifiedBy` is unrecorded), preventing collector-as-determiner false Gold promotions.
  - **Layer 2 (Protologue Author Extraction):** `extract_protologue_authors()` parses basionym and combination authors directly from Darwin Core `scientificName` strings; sheets annotated by a nomenclatural author are classified as `Tier_1_Gold`.
  - **Layer 3 (Corpus Determiners & Regional Herbaria):** Determiners identifying $\ge 5$ specimens across $\ge 3$ herbaria become `Tier_1_Gold` (weighting type annotations $2\times$ and requiring $\le 25\%$ duplicate conflict). Curated herbaria with $\ge 5$ sheets holding $\ge 1$ type are classified as `Tier_2_Silver`. Unverified collections receive `Tier_3_Bronze`.

---

## 3. Phase 2 Protocol: Pure LeafMachine3 ONNX Pipeline & Botanical Gatekeeping (`02_segment_and_extract.py`)

Candidate organs are detected, segmented, and extracted using **LeafMachine3** (`lm3`) ONNX Runtime models coupled with specialized *Packera* botanical gatekeepers:

### 3.1 LeafMachine3 18-Module Pipeline Components
* **Scale Calibration (Modules 1, 6 & 7):** Module 1 (`lm3_mp_conversion_factor__sqrt_fit`, $R^2 = 0.905$, MAE 2.3%) scale prior; Module 6 (`yolo26x_cls_224` + `yolo26n_cls_224` + `dinov2_frozen_mlp`, 99.0% accuracy) ruler detection; Module 7 (Tick-Lattice Solver & FieldPrism markers) physical scale calibration (`pixels_per_mm`, `is_calibrated_mm`).
* **Organ Detection & Whole-Sheet Specimen Segmentation (Modules 2, 3 & 4):** Modules 2 & 3 (`lm3_plant_detector__yolo26x_det_1280` Archival & Plant Detectors; `Leaflet` merged into `Leaf_WHOLE`); Module 4 (`lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain` [UNet++ v9]) with `LLR_Hyst_KeepBlobs: true` color likelihood removing paper pockets in deep sinuses while preserving disconnected pinnae.
* **Single-Leaf Instance Segmentation & 31-Keypoint Pose Skeleton (Modules 8, 10, 12, 13 & 14):** Module 8 (`lm3_leaf_segmenter__yolo26x_seg_1024`) segments `Leaf` lamina (0.933 mAP50), `Petiole` (0.818 mAP50), and `Hole` perforations (0.429 mAP50); Module 10 (`lm3_landmark_detector__yolo26x_pose_640`) infers 31-keypoint `mid15_pet5` pose (15 points along midvein, 6 along petiole); Module 12 rotates leaf to standard vertical tip-up orientation; Modules 13 & 14 extract petiole width and midvein bilateral symmetry.

### 3.2 *Packera* Botanical Gatekeeping (`leaf_gatekeeper.py` & `organ_crown_metrics.py`)
* **DBSCAN Plant Decomposition (`cluster_voucher_plants_dbscan`):** Groups organs into individual plants (`plant_individual_id = 0, 1, ...`), estimating rosette centroids and extents.
* **Rosette Caudex Crown Estimation (`estimate_rosette_crown_point`):** Identifies root-shoot transition zone ($p_{\text{crown}}$). Petiole length is evaluated as outer-rosette geodesic distance:
  $$L_{\text{petiole}} = \max\left(d_{\text{geodesic}}(p_{\text{base}}, p_{\text{crown}}),\, \|p_{\text{base}} - p_{\text{crown}}\|_2\right)$$
  aggregated at $P_{75}$ across outer basal leaves. Scape height ($H_{\text{plant}}$) accounts for bent stems:
  $$H_{\text{plant}} = \|p_{\text{bend}} - p_{\text{crown}}\|_2 + \|p_{\text{corymb}} - p_{\text{bend}}\|_2$$
* **Cauline / Stem Rejection (`is_cauline_or_stem_contaminant`):** Rejects cauline bracts and stems via aspect ratio $> 3.8$, solidity $< 0.78$, or height $> 2.2 \times \text{median rosette radius}$ above crown.
* **Boundary-Continuity Lobe Rejection (`check_boundary_tissue_continuity`):** Inspects a 2-px perimeter in CIE $L^*a^*b^*$. Candidates with $> 15\%$ non-paper tissue crossing borders are rejected.
* **Dissected-Leaf Rachis Bridging:** In deeply pinnatifid leaves with narrow segmentation breaks, UNet++ v9 `specimen_paper_cleaned_mask` acts as an anatomical constraint (`bridge_dissected_rachis_components`), reconnecting dissected pinnae without bridging across mounting paper.
* **Tier 2 Hemi-Blade Reflection (`reflect_folded_mask`):** For folded blades, reflects intact half across 15-point midvein into `Tier_2_Reflected` for auxiliary linear traits.
* **Classical Gatekeeper & Dissection Gates (`Tier_1_Direct`):** Entire blades require $\text{Solidity} \ge 0.72$, isoperimetric $\text{Circularity} = \frac{4\pi A}{P^2} \ge 0.25$, $\text{Extent} = \frac{A}{W \cdot H} \ge 0.35$, and $\text{UCS} \ge 0.85$. Pinnatifid ($0.48 \le \text{Solidity} < 0.72$) and bipinnatisect ($0.18 \le \text{Solidity} < 0.48$, $\ge 4$ deep sinuses) blades are validated via `is_botanical_dissection()` and routed into `Tier_1_Direct`.
* **Cohort-Adaptive Physical Blade Area Floor:**
  $$\text{effective\_min\_area\_cm2} = \min(0.75, \max(0.15, 0.25 \times \text{median\_sheet\_leaf\_area\_cm2}))$$
  preventing dwarf or alpine specimens from being arbitrarily discarded.
* **Post-QC Basal Masking (`crop_dense_rosette_tier4`):** Replaces non-basal-leaf pixels outside 5-px dilated union of verified basal leaf masks with neutral paper `[242, 239, 230]`.

---

## 4. Phase 3 & 4 Protocol: Morphometrics, EFA, Allometry & Ordination (`03` & `04`)

### 4.1 Shape Quantification & Strict `Tier_1_Direct` EFA (`03_fourier_extractor.R`)
* **Momocs Outline Export & Euler Characteristic Transform (ECT):**
  LeafMachine3 exports native tip-up, base-anchored clockwise outlines (`momocs_outlines.json`) paired with 360-direction **Euler Characteristic Transform (ECT)** coordinates (`Coordinates/*.h5`, Ayub et al. 2026; Weaver 2026).
* **Strict `Tier_1_Direct` Closed-Contour Restriction:**
  Only intact `Tier_1_Direct` contours enter `Momocs::efourier`. Folded/reflected `Tier_2_Reflected` leaves are never pooled into EFA. Ingests contours via batch SQL from `extracted_leaves.contour_xy_json`.
* **Inertial Axis Pre-Alignment & Harmonic Power Calibration:**
  Contours are aligned along major inertial axis via `Momocs::coo_alignxax()` with circularity $> 0.75$ suborbicular guard. Harmonic power:
  $$\text{Power}_n = \frac{A_n^2 + B_n^2 + C_n^2 + D_n^2}{2},\quad \frac{\sum_{n=1}^{H_{\text{opt}}} \text{Power}_n}{\sum_{n=1}^{H_{\max}} \text{Power}_n} \ge 0.99$$
  dynamically determines minimum harmonic count $H_{\text{opt}} \in [6, 24]$ capturing $\ge 99\%$ cumulative power, executing `Momocs::efourier(out_obj, nb.h = H_opt, norm = FALSE, start = FALSE)`.
* **Uniform Bilateral Symmetrization (Claude 2008):**
  $$A_n^{\text{sym}} = A_n,\quad D_n^{\text{sym}} = D_n,\quad B_n^{\text{sym}} = 0,\quad C_n^{\text{sym}} = 0,\quad \text{pct\_asymmetry} = \frac{\sum_{n=1}^{H_{\text{opt}}} (B_n^2 + C_n^2)}{\sum_{n=1}^{H_{\text{opt}}} (A_n^2 + B_n^2 + C_n^2 + D_n^2)}$$
* **Type Outlier Immunity & OTU Median Aggregation:**
  Intact type contours are exempt from within-taxon Mahalanobis and MCD outlier screening (`mcd_out[is_type_specimen] <- FALSE`). Active type anchors receive a $1.5\times$ precision boost in `stats::cov.wt()`. Exports PCs accounting for $\ge 95\%$ variance ($k_{\text{export}} = \max(5, k_{95})$), followed by unweighted coordinate-wise OTU median aggregation to `(catalogNumber, plant_individual_id)`.

### 4.2 Determiner Authority Weighting & Adaptive Reference Anchors (`authority_weighting_utils.R`)
* **Adaptive Reference Anchors & Type-Medoid Shrinkage:**
  Reference anchors consist of `Tier_1_Gold` and `Tier_2_Silver` vouchers pruned once at $\chi^2_{0.95}$ FastMCD. Intact types receive retention immunity (`anchor_tier_source = "Tier_0_Type_Anchor"`), while extreme types ($> \chi^2_{0.99}$) are isolated from FastMCD scatter estimation. Taxon reference centers shrink toward the type medoid:
  $$\boldsymbol{\mu}_{\text{ref}, k} \leftarrow (1 - \gamma_k) \boldsymbol{\mu}_{\text{ref}, k} + \gamma_k \boldsymbol{\mu}_{\text{type}, k},\quad \gamma_k = \min\!\left(0.35,\, \frac{2.5\,n_{\text{eff, type}, k}}{N_{\text{gs}, k} + 2.5\,n_{\text{eff, type}, k}}\right)$$
  Taxa with $< 8$ Gold+Silver anchors promote central `Tier_3_Bronze` vouchers within $\chi^2_{0.80}$ (`is_promoted_anchor = TRUE`).
* **2-Step Dawid-Skene / Beta-Binomial Empirical Bayes Authority Weighting:**
  Estimates observer error rates $\hat{\epsilon}_i$ exogenously (`target_modality = ...`; Dawid & Skene 1979) and pools priors ($\alpha_0 = 5.0$, prior reliabilities $\pi_{\text{PrimaryType}} = 0.999$, $\pi_{\text{SecType}} = 0.985$, $\pi_{\text{Gold}} = 0.95$, $\pi_{\text{Silver}} = 0.85$, $\pi_{\text{PromotedBronze}} = 0.65$, $\pi_{\text{Bronze}} = 0.35$). Types are exempt from downward shrinkage. Unpromoted Bronze receive $0.0$ in `active_passive` mode.

### 4.3 Klingenberg Allometry, Dual-Track CDA & MANOVA (`04_gmm_morphotools.R`)
* **Two-Sample Prediction-Variance Corrected Type Distance:**
  $$D^2_{\text{type}, i} = \left(1 + \frac{1}{n_{\text{eff, type}, k}}\right)^{-1} (\mathbf{x}_i - \boldsymbol{\mu}_{\text{type}, k})^\top \mathbf{S}_{\text{reg}}^{-1} (\mathbf{x}_i - \boldsymbol{\mu}_{\text{type}, k})$$
  persisted as `type_morph_mahal_d2`, alongside `nearest_type_morph_taxon` and `is_typical_morphotype` (GMM cluster $c^*$).
* **Klingenberg (2016) Allometric Slope Homogeneity:**
  Evaluates slope homogeneity across ruler-calibrated anchors (`is_calibrated_mm == TRUE`) via multivariate interaction models in `RRPP::lm.rrpp(..., weights = eff_weight)`:
  $$\mathbf{y}_{\text{shape}} = \beta_0 + \beta_1 \log_{10}(\text{CS}) + \boldsymbol{\gamma}_{\text{species}} + \boldsymbol{\delta} (\log_{10}(\text{CS}) \times \boldsymbol{\gamma}_{\text{species}}) + \boldsymbol{\epsilon}$$
  If $p \ge 0.01$, common-slope correction is applied: $\mathbf{y}_{\text{corrected}} = \mathbf{y}_{\text{shape}} - \hat{\beta}_1 (\log_{10}(\text{CS}) - \overline{\log_{10}(\text{CS})})$.
* **Subspace Selection & Dual-Track CDA:**
  Horn's Parallel Analysis capped by Williams & Titus (1988) $P_{\max} = \lfloor (N_{\min,\text{anchor}} - 1) / 3 \rfloor$. Evaluates honest LOOCV on weighted Ledoit-Wolf shrinkage CDA. Track 1: Unsupervised GMM Morphotypes (`Morphotype_K1..Kn`) via Kass & Raftery (1995) $2\Delta\text{BIC}$ thresholds ($> 6$ strong, $> 10$ very strong). Track 2: Supervised EPA Level III Ecoregion CDA biplot and MANOVA with Benjamini-Hochberg FDR.

---

## 5. Phase 5 Protocol: Taphonomic Indumentum-Blind Vision XAI & Cleanlab (`05_cleanlab_vision_xai.py`)

* **Taphonomic Indumentum-Blind Preprocessing:**
  Applies $\sigma = 1.0\text{ px}$ Gaussian low-pass filter on $224 \times 224$ crops (`apply_indumentum_blind_lowpass`) with non-leaf pixels masked to neutral paper `[242, 239, 230]`, attenuating abraded hairs while preserving venation and lobing.
* **Foreground-Masked $16 \times 16$ Patch Token Pooling (`dinov2_vitb14_reg`):**
  Uses `dinov2_vitb14_reg` (4 register tokens) and pools strictly over foreground spatial patch tokens ($\ge 25\%$ foreground coverage). Flags `foreground_attention_ratio < 0.65` as `is_artifact_driven = TRUE`.
* **768-d DINOv2 Type Cosine Control (`compute_type_visual_control`):**
  Computes L2-normalized class visual prototypes $\mathbf{v}_{\text{type}, c}$ as weighted spherical mean across primary ($2.0\times$) and secondary ($1.0\times$) types, blending single-specimen type vectors with Gold class medoid:
  $$\text{vision\_type\_cosine\_sim}_i = \frac{\mathbf{v}_i \cdot \mathbf{v}_{\text{type}, c}}{\|\mathbf{v}_i\|_2 \|\mathbf{v}_{\text{type}, c}\|_2},\quad \text{vision\_type\_margin}_i = \cos(\mathbf{v}_i, \mathbf{v}_{\text{type}, \text{nearest\_other}}) - \cos(\mathbf{v}_i, \mathbf{v}_{\text{type}, \text{given}})$$
* **Class-Conditional Cleanlab Noise Auditing & OOD Scoring:**
  Trains fold-internal `LogisticRegressionCV` propagating `best_c` into `CalibratedClassifierCV` under `StratifiedGroupKFold(groups=institutionCode)` with weights bounded to `[0.35, 3.0]` and types assigned weight 4.5. Primary types are immune to corruption flags (`is_label_corrupted = False`). Computes feature-space Mahalanobis + cosine $k$-NN distance (`vision_ood_score`, `vision_is_ood = TRUE` when $\ge 0.65$). Caches embeddings in `dinov2_feature_cache`.

---

## 6. Phase 6 Protocol: Multimodal Spatial RF, *Locus Classicus* Ecological Control & SDM Ensembles (`06_multimodal_spatial_rf.R`)

* **Type Locality (*Locus Classicus*) Ecological Baseline:**
  Direct point-raster extraction on types is gated by coordinate uncertainty ($\le 2000\text{ m}$ edaphic/topographic, $\le 10000\text{ m}$ climatic; Daru et al. 2018). Coarse type sheets derive *locus classicus* edaphic profiles from proximate Gold vouchers within $150\text{ km}$ or sharing the EPA Level III ecoregion. Physiological divergence is evaluated via weighted RMS distance:
  $$z_{\text{type}, i} = \sqrt{\frac{\sum_{j \in \text{Primary}} w_j \left(\frac{x_{i, j} - \mu_{\text{type}, k, j}}{\max(\text{MAD}_{\text{type}, k, j},\, s_{\text{floor}, j},\, 10^{-4})}\right)^2}{\sum_{j \in \text{Primary}} w_j}}$$
  exported as `type_edaphic_z_score`, alongside `nearest_type_niche_taxon` and boolean `matches_type_lithology`.
* **Full 19-Bioclimatic & 0–15 cm Rhizosphere Stack:**
  Extracts 19 WorldClim v2.1 bioclimatic variables (`bio1`..`bio19`) and 0–15 cm rhizosphere depth-weighting (`(5 * v_0_5 + 10 * v_5_15) / 15`) across SSURGO, 30m POLARIS, 100m SOLUS100, 3DEP topography, and 6-class USGS SGMC lithology.
* **3-Stage Dynamic Primary Predictor Gating:**
  (1) Prunes collinear predictors ($|r| < 0.70$, $\text{VIF} < 5.0$) across Climatic, Edaphic, and Topographic blocks; (2) Trains `ranger::ranger` under $0.5^\circ \times 0.5^\circ$ Spatial Block CV with continuous Empirical Bayes weights (primary types pinned to 3.0 ceiling) and winsorized case weights `[0.25, 4.0]`, isolating variables accounting for top 85% cumulative importance; (3) Restricts downstream niche distances strictly to discoverable primary predictors.
* **SDM Ensembles (MaxEnt + Spatial RF):**
  Fits authority-weighted MaxEnt (`maxnet`, `cloglog` link) and Spatial RF ensembles with target-group background selection, exporting spatial-block $AUC_{\text{block}}$ and Continuous Boyce Index (CBI) to reports and GeoTIFFs.

---

## 7. Phase 7 & 8 Protocol: Hierarchical Nomenclatural Type Decision Matrix & Complex Discovery (`07` & `08`)

* **Master View Synthesis & Reliability-Weighted Consensus:**
  Consensus synthesis queries `v_integrative_triage_master`, left-joining all modalities onto `plant_individuals`. Weights modalities by out-of-fold Cohen's $\kappa_m$ and class-conditional precision $C^{(m)}_{k,j}$:
  $$w_{i,m} = w_m \cdot (0.5 + 0.5 \cdot \text{Precision}_m(\hat{y}_{i,m}))$$
  achieving 98.4% voucher completion across georeferenced and non-georeferenced sheets.
* **Hopkins Circular Phenology:**
  Regresses circular anthesis DOY ($\theta_i = \frac{2\pi \cdot \text{DOY}_i}{365.25}$) on latitude, elevation, and temperature, scaled by $\alpha = \min(1.0, \max(0.0, (N_{\text{anthesis}} - 5)/15))$. Flowering dispersion scales thresholds ($\tau_{\text{cryptic}} = \min(35, \max(14, 2.5\,\sigma_{\text{DOY}}))$, $\tau_{\text{ecotype}} = \min(28, \max(10, 2.0\,\sigma_{\text{DOY}}))$).
* **Taxon-Agnostic 5-Stage Hierarchical Decision Matrix:**
  - **Rule A (Nomenclatural Type Override):** Type specimens are never reassigned (`recommended_determination = given_sp`), nor routed to `Misidentification_Reassign`, `Hybrid_Introgressant`, or `Putative_Cryptic_Species`. Concordant types receive `taxonomic_status_call = "Confirmed_Type_Control"` (`RESOLVED`). Cross-modal divergence triggers `taxonomic_status_call = "Type_Synonymy_Candidate"` (`CRITICAL`).
  - **Rule B (Type-Control Reassignment):** Non-types deviating from type baselines ($D^2_{\text{type}} > \chi^2_{0.95}$, $z_{\text{type}} > 2.0$) with $\ge 2$ type streams agreeing on $T_{\text{alt}} \ne \text{given\_sp}$ without classifier contradiction are reassigned to $T_{\text{alt}}$. Vision reassignments require `!is_artifact_driven`.
  - **Rule C (*Sensu Stricto* Typical Morphotype Locking):** Specimens in typical morphotype GMM cluster ($c^*$) with $z_{\text{type}} \le 2.0$ and `!vision_is_ood` are locked as `Confirmed_Species` (*sensu stricto*).
  - **Stage 2 (`Hybrid_Introgressant`):** Intermediate CDA posteriors ($0.35 \le P_1 \le 0.65$, $P_2 \ge 0.25$) in sympatry, resolving parents via `cda_second_taxon`.
  - **Stage 3 (`Putative_Cryptic_Species`):** Requires `morph_tier_source == "Tier_1_Pure"`, GMM $2\Delta\text{BIC} > 6$, deep visual divergence (`vision_is_ood`), and independent edaphic/phenological divergence.
  - **Stage 4 (`Subspecies_Ecotype` vs. `Ecophenotypic_Plasticity`):** Evaluates cohesive geographic/edaphic divergence for ecotypes versus overlapping envelopes for plasticity.
* **Step 08 Taxon-Complex Discovery (`08_discover_taxon_complexes.R`):**
  Fuses pairwise morphological, fine-scale edaphic, climatic, phenological, and vision distances, shrinking centroids toward type baselines with weight $\gamma_k$. Discovers complexes via unsupervised graph modularity (Leiden/Louvain) and Ward.D2 clustering, resolving Taxon x GMM Morphotype subclusters when $N_{\text{taxa}} \le 2$.

### Integrative Taxonomic Decision Matrix

| Taxonomic Status Call | Vegetative Morphometrics | Deep Vision XAI | Cytology ($2n$) | nextRAD Phylogenomics | Rhizosphere Edaphics & Niche | Actionable Outcome |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`Confirmed_Type_Control`** | Type centroid ($\boldsymbol{\mu}_{\text{type}}$); `is_typical_morphotype` | High cosine to type ($\mathbf{v}_{\text{type}}$) | Reference cytotype | Type/paratype genetic profile | High fidelity to *locus classicus* | Immutable nomenclatural anchor; retain determination (`RESOLVED`) |
| **`Type_Synonymy_Candidate`** | Projects into congener morphospace ($P \ge 0.85$) | High visual cosine to congener | Shared/intermediate cytotype | Clusters with congener lineages | Overlaps congener ecological envelope | Critical priority: evaluate formal taxonomic synonymy (`CRITICAL`) |
| **`Putative_Cryptic_Species`** | Discrete GMM ($2\Delta\text{BIC} > 6$); CDA $P < 0.70$ | High OOD (`vision_is_ood`); low confidence | Distinct/polyploid ($2n = 44, 46, 88$) | Diagnostic monophyletic cluster | Significant edaphic divergence or allochrony | Formal species-level recognition / resurrection under USC |
| **`Misidentification_Reassign`** | High posterior for known congener; Rule B | High confidence for alternate; clean | Matches reassigned congener | Clusters within reassigned congener | Matches edaphic/climatic envelope of congener | Annotate & reassign determination label |
| **`Hybrid_Introgressant`** | Intermediate CDA ($0.35\text{–}0.65$, $P_2 \ge 0.25$) | Intermediate vision probabilities | Aneuploid, triploid, or mixed | Admixed ancestry; reticulation | Sympatric contact zones or ecotones | Designation as nothospecies / hybrid swarm |
| **`Subspecies_Ecotype`** | Partial shift ($2\Delta\text{BIC} \le 6$); conserved outline | Conserved macro-architecture | Shared cytotype with nominal | Shallow regional structure | Specialized substrate ($D_{\text{edaph}} < D_{\text{clim}}$) | Infraspecific / varietal recognition |
| **`Ecophenotypic_Plasticity`** | Unimodal GMM ($K=1$); parallel slope | Concordant with given taxon | Identical cytotype | Panmictic with parental populations | Overlapping climatic/edaphic envelopes | Retain under polymorphic species |

---

## 8. Literature Cited

1. Ayub, M. S., et al. 2026. Topological shape descriptors and ECT in plant phenomics. *bioRxiv* https://doi.org/10.64898/2026.04.13.718293
2. Barkley, T. M. 1988. Variation among Senecioneae in North America. *Brittonia* 40: 211–221.
3. Broennimann, O., et al. 2012. Measuring ecological niche overlap. *Glob Ecol Biogeogr* 21: 481–497.
4. Claude, J. 2008. *Morphometrics with R*. Springer, New York.
5. Collyer, M. L., & D. C. Adams. 2018. RRPP: Linear models via residual randomization. *Methods Ecol Evol* 9: 1772–1779.
6. Darcet, T., et al. 2024. Vision transformers need registers. *ICLR*.
7. Daru, B. H., et al. 2018. Widespread sampling biases in herbaria. *New Phytol* 217: 939–955.
8. Dawid, A. P., & A. M. Skene. 1979. Maximum likelihood estimation of observer error-rates. *Appl Stat* 28: 20–28.
9. de Queiroz, K. 2007. Species concepts and species delimitation. *Syst Biol* 56: 879–886.
10. Hirzel, A. H., et al. 2006. Evaluating ability of habitat suitability models to predict presence. *Ecol Modell* 199: 142–152.
11. Kass, R. E., & A. E. Raftery. 1995. Bayes Factors. *JASA* 90: 773–795.
12. Klingenberg, C. P. 2016. Size, shape, and form: concepts of allometry in geometric morphometrics. *Dev Genes Evol* 226: 113–137.
13. Kowal, R. R. 1975. Systematics of *Senecio aureus* and allied species on Gaspé Peninsula. *Mem Torrey Bot Club* 23: 1–113.
14. Kuhl, F. P., & C. R. Giardina. 1982. Elliptic Fourier features of a closed contour. *Comput Graph Image Process* 18: 236–258.
15. Northcutt, C. G., et al. 2021. Confident Learning: Estimating uncertainty in dataset labels. *JAIR* 70: 1373–1411.
16. Omernik, J. M., & G. E. Griffith. 2014. Ecoregions of conterminous United States. *Environ Manage* 54: 1249–1266.
17. Padial, J. M., et al. 2010. The integrative future of taxonomy. *Front Zool* 7: 16.
18. Phillips, S. J., et al. 2009. Sample selection bias in Maxent. *Ecography* 32: 181–197.
19. Phillips, S. J., et al. 2017. Opening the black box: open-source release of Maxent. *Ecography* 40: 887–893.
20. Scrucca, L., & A. E. Raftery. 2018. clustvarsel: Variable selection for model-based clustering. *J Stat Softw* 84: 1–28.
21. Šlenker, M., et al. 2022. MorphoTools2: multivariate morphometric analysis. *Bioinformatics* 38: 2954–2955.
22. Trock, D. K. 2006. *Packera*. In *Flora of North America*, Vol. 20, 570–602. Oxford University Press, New York.
23. Turland, N. J., et al. 2018. *International Code of Nomenclature (Shenzhen Code)*. Regnum Veg 159. Koeltz, Glashütten.
24. Valavi, R., et al. 2019. blockCV: Spatially separated folds for k-fold cross-validation. *Methods Ecol Evol* 10: 225–232.
25. Warren, D. L., et al. 2008. Environmental niche equivalency versus conservatism. *Evolution* 62: 2868–2883.
26. Weakley, A. S. 2026. *Flora of the Southeastern United States*. UNC Herbarium (NCU), Chapel Hill.
27. Weaver, W. N. 2026. LeafMachine3: Automated machine vision for herbarium trait phenotyping. https://github.com/Gene-Weaver/LeafMachine3
28. Williams, B. K., & K. Titus. 1988. Sampling stability in discriminant analysis. *Ecology* 69: 1275–1285.
29. Wright, M. N., & A. Ziegler. 2017. ranger: Fast implementation of random forests. *J Stat Softw* 77: 1–17.
