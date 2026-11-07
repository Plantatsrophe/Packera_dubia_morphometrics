# Standard Operating Procedure (SOP)

## Automated High-Throughput Vegetative Phenotyping and Herbarium Triage Pipeline for North American *Packera*
### Quantitative Morphological Diagnosability and Species Boundaries under the Unified Species Concept (*Packera dubia* Complex Baseline)

- **Protocol Identifier:** `UNC-BOT-SOP-2026-04-REV5`
- **Principal Investigator:** J. Brandon Fuller (PhD Candidate, Department of Biology, UNC-CH)
- **Faculty Advisor:** Dr. Alan S. Weakley (Director, NCU Herbarium; UNC Biology)
- **Target Taxa:** ***Packera dubia*** (Spreng.) Trock & Mabb. (Asteraceae: Senecioneae) as primary baseline, evaluated dynamically against all regional and continent-wide North American *Packera* Á. Löve & D. Löve congeners.

---

## 1. Executive Summary and Theoretical Framework

I established this protocol to conduct automated high-throughput vegetative phenotyping and herbarium triage across North American ***Packera* Á. Löve & D. Löve** (Asteraceae: Senecioneae), using the taxonomically recalcitrant ***Packera dubia* complex** as my primary dissertation baseline. The lineage presents acute identification challenges on pressed sheets, as developmental heterophylly produces radical shape shifts between juvenile rosettes and mature cauline leaves, whereas manual curatorial annotations frequently confuse sympatric congeners.

I frame morphological diagnosability within the **Unified Species Concept (USC; de Queiroz 2007)** through **Integration by Congruence (Padial et al. 2010)**. Under the USC, species represent separately evolving metapopulation lineages. Lineage divergence generates morphological diagnosability, reproductive isolation, ecological divergence, and reciprocal monophyly as contingent properties along the speciation continuum. I evaluate automated vegetative shape data alongside four independent lines of evolutionary evidence:

1. **Reproductive Macro-Morphology and Anthesis:** Quantitative capitulum aspect ratios, phyllary counts, floret dimensions, and 3-tier optical anthesis staging.
2. **Cytotaxonomic Cytotypes:** Chromosome numbers and polyploid series ($2n = 44, 46, 88$; Kowal 1975) that separate diploid and tetraploid lineages.
3. **Reduced-Representation nextRAD Phylogenomics:** High-density SNP datasets resolving reticulate evolution, hybrid swarms, and lineage boundaries.
4. **0 to 15 cm Rhizosphere and USGS SGMC Bedrock Stack:** Edaphic data cached in SQLite `environmental_cache` and `environmental_rtree`, integrating 0 to 15 cm depth-weighting (`(5 * v_0_5 + 10 * v_5_15) / 15`), SSURGO map units (`ssurgo_rock_outcrop_pct`, `ssurgo_bedrock_depth_cm`), 30 m POLARIS, 100 m SOLUS100, 6-class USGS SGMC lithology (`sgmc_lithology_class`), 10 m to 30 m USGS 3DEP topography, and Climatic ($D_{\text{clim}}$) versus Edaphic ($D_{\text{edaph}}$) niche tests (Warren et al. 2008).
5. **Central Relational WAL Architecture (`data/packera_pipeline.sqlite`):** All pipeline stages transact against a centralized SQLite WAL database housing Tier A shared tables (`taxon_synonyms`, `vouchers`, `plant_individuals`, `extracted_leaves`, `voucher_reproductive_phenology`, `environmental_cache` alongside `rtree`) and Tier B run-scoped tables (`pipeline_runs`, `run_active_cohort`, `leaf_efa_harmonics`, `otu_morphometrics`, `spatial_rf_predictions`, `curation_triage`, and master view `v_integrative_triage_master`).

This multi-tiered relational design separates immutable raw observations from run-specific statistical cohorts, so I can test alternative taxonomic hypotheses without risking data contamination.

---

## 2. Phase 1 Protocol: Ingestion, Standardization and US EPA Level III Ecoregions (`01_voucher_harvester.py`)

I designed Phase 1 to solve the foundational challenge of herbarium data intake: raw specimen databases contain corrupted taxonomic names, collector duplicate conflicts, and geographic clustering that distort statistical analyses if ingested uncritically.

* **Darwin Core Schema Contract and Relational UPSERT (`vouchers`):**
  I harvest occurrences into SQLite `vouchers` (PK: `catalogNumber`) via atomic UPSERT while resolving synonyms through `taxon_synonyms`, exporting synchronized CSV `data/tables/curated_vouchers.csv` with 26 Darwin Core columns: `catalogNumber`, `institutionCode`, `scientificName`, `species_raw`, `recordedBy`, `recordNumber`, `identifiedBy`, `determiner_raw`, `determiner_tier`, `type_status`, `county`, `stateProvince`, `decimalLatitude`, `latitude`, `decimalLongitude`, `longitude`, `coordinateUncertainty`, `year`, `month`, `day`, `eventDate`, `regional_group`, `image_path`, `exsiccatae_key`, `is_primary_duplicate`, `duplicate_count`, `exsiccatae_taxonomic_conflict`, and `exsiccatae_gold_taxon`.
* **Dynamic Taxon Standardization and GBIF Synonym Expansion:**
  I standardize binomials and basionyms dynamically through `expand_taxon_synonyms_gbif()`. Omitting taxon flags queries all *Packera* species across the genus through `discover_genus_packera_species_gbif()`. When I pass `--taxon` or `--focal-taxon`, the harvester restricts analysis strictly to those specified taxa.
* **Configurable Geographic Scope:**
  I support `--geographic-scope {eastern_na,north_america,global}` (default: `eastern_na`, enforcing `decimalLongitude` bounds `[-106.65, -50.0]` during GBIF queries with a `typeStatus="*"` rescue query for historical types; `north_america` and `global` include western occurrences). When `exclude_western` is active, I exclude Western US states and Western or Northern Canadian provinces (`BC, AB, SK, YT, NT, NU`) and require $\ge 3$ georeferenced records east of $-106.65^\circ\text{W}$, supported by `--countries` (default `["US"]`).
* **Pre-Download Deduplication and WebP Storage:**
  I deduplicate exsiccatae (`primary_duplicates_only=True`) before downloading images, retaining primary sheets (`is_primary_duplicate == TRUE`) while unconditionally preserving verified nomenclatural types. Full-sheet imagery in `data/raw_vouchers/` and crops in `data/cropped_patches/` default to quality-85 `.webp` compression (`method=4` via PIL), which yields a 70% to 80% disk reduction relative to JPEG or PNG with zero edge degradation.
* **Official US EPA Level III Ecoregion Point-in-Polygon Delineation:**
  I assign `regional_group` strictly through spatial intersection against official **US EPA Level III Ecoregions** (Omernik & Griffith 2014; `data/environmental/epa_ecoregions_l3/`):
  $$\text{regional\_group}_i = \text{Ecoregion}_k \iff (x_i, y_i) \in \mathcal{P}_k$$
* **Ungeoreferenced Specimen Fallback Policy:**
  When coordinate data are missing ($y_i = \text{NA}$), I evaluate the US Census county boundary polygon $\mathcal{C}_i$:
  $$\text{regional\_group}_i = \begin{cases} \text{Ecoregion}_k & \text{if } \frac{\text{Area}(\mathcal{C}_i \cap \mathcal{P}_k)}{\text{Area}(\mathcal{C}_i)} \ge 0.95 \\ \text{NA} & \text{otherwise} \end{cases}$$
* **Nomenclatural Type Hierarchy and Priority Engine:**
  - `PRIMARY_TYPE_STATUSES = {"HOLOTYPE", "LECTOTYPE", "NEOTYPE", "EPITYPE", "ISOTYPE", "ISOLECTOTYPE", "ISONEOTYPE"}` (`Tier_0_Primary_Type`, rank 2).
  - `SECONDARY_TYPE_STATUSES = {"SYNTYPE", "ISOSYNTYPE", "PARATYPE", "ISOPARATYPE", "TYPE", "COTYPE"}` (`Tier_0_Secondary_Type`, rank 1).
  - I resolve fine-grained ICN priority through `classify_nomenclatural_type_subrank()` (`2.2` primary, `2.0` isotype, `1.5` syntype, `1.0` paratype; Turland et al. 2018). I exempt verified types from uncertainty or partial-date filters and retain them unconditionally during spatial thinning and outgroup selection.
* **3-Layer Dynamic Authority Engine (`parse_determiner_tier`):**
  - **Layer 1 (Monographer Heuristics and Determiner Isolation):** I restrict `SPECIALIST_PATTERNS` matching strictly to `identifiedBy` and `identificationHistory` (falling back to `recordedBy` only when `identifiedBy` is absent), preventing collector-as-determiner false Gold promotions.
  - **Layer 2 (Protologue Author Extraction):** I parse basionym and combination authors directly from Darwin Core `scientificName` strings with `extract_protologue_authors()`; sheets annotated by a nomenclatural author become `Tier_1_Gold`.
  - **Layer 3 (Corpus Determiners and Regional Herbaria):** Determiners identifying $\ge 5$ specimens across $\ge 3$ herbaria become `Tier_1_Gold` (weighting type annotations $2\times$ and requiring $\le 25\%$ duplicate conflict). Curated herbaria with $\ge 5$ sheets holding $\ge 1$ type become `Tier_2_Silver`, whereas unverified collections receive `Tier_3_Bronze`.

Now that Phase 1 establishes an authoritative, georeferenced specimen cohort, I feed these vetted vouchers directly into the vision pipeline.

---

## 3. Phase 2 Protocol: Pure LeafMachine3 ONNX Pipeline and Botanical Gatekeeping (`02_segment_and_extract.py`)

Herbarium sheets present challenging visual artifacts: white mounting tape, fragment packets, discolored paper, and overlapping plants frequently cause standard bounding-box detectors to fail. I address this challenge by deploying **LeafMachine3** (`lm3`) ONNX Runtime models alongside botanical gatekeeper algorithms.

### 3.1 LeafMachine3 18-Module Pipeline Components
* **Scale Calibration (Modules 1, 6 and 7):** Module 1 (`lm3_mp_conversion_factor__sqrt_fit`, $R^2 = 0.905$, MAE 2.3%) predicts scale priors; Module 6 ensemble detects rulers; Module 7 (Tick-Lattice Solver and FieldPrism markers) establishes physical scale (`pixels_per_mm`, `is_calibrated_mm`).
* **Organ Detection and Whole-Sheet Specimen Segmentation (Modules 2, 3 and 4):** Modules 2 and 3 (`lm3_plant_detector__yolo26x_det_1280` Archival and Plant Detectors; `Leaflet` merged into `Leaf_WHOLE`) identify structures; Module 4 (`lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain` [UNet++ v9]) with `LLR_Hyst_KeepBlobs: true` color likelihood removes paper pockets in deep sinuses while preserving disconnected pinnae.
* **Single-Leaf Instance Segmentation and 31-Keypoint Pose Skeleton (Modules 8, 10, 12, 13 and 14):** Module 8 (`lm3_leaf_segmenter__yolo26x_seg_1024`) segments `Leaf` lamina (0.933 mAP50), `Petiole` (0.818 mAP50), and `Hole` perforations (0.429 mAP50); Module 10 (`lm3_landmark_detector__yolo26x_pose_640`) infers 31-keypoint `mid15_pet5` pose (15 points along midvein, 6 along petiole); Module 12 rotates leaves to standard vertical tip-up orientation; Modules 13 and 14 extract petiole width and midvein bilateral symmetry.

### 3.2 *Packera* Botanical Gatekeeping (`leaf_gatekeeper.py` and `organ_crown_metrics.py`)
* **DBSCAN Plant Decomposition (`cluster_voucher_plants_dbscan`):** I cluster organs into individual plants (`plant_individual_id = 0, 1, ...`), estimating rosette centroids and extents.
* **Rosette Caudex Crown Estimation (`estimate_rosette_crown_point`):** I identify the root-shoot transition zone ($p_{\text{crown}}$). I measure petiole length as outer-rosette geodesic distance:
  $$L_{\text{petiole}} = \max\left(d_{\text{geodesic}}(p_{\text{base}}, p_{\text{crown}}),\, \|p_{\text{base}} - p_{\text{crown}}\|_2\right)$$
  aggregated at $P_{75}$ across outer basal leaves. I evaluate scape height ($H_{\text{plant}}$) accounting for bent stems:
  $$H_{\text{plant}} = \|p_{\text{bend}} - p_{\text{crown}}\|_2 + \|p_{\text{corymb}} - p_{\text{bend}}\|_2$$
* **Cauline and Stem Rejection (`is_cauline_or_stem_contaminant`):** I reject cauline bracts and stem fragments through aspect ratio $> 3.8$, solidity $< 0.78$, or height $> 2.2 \times \text{median rosette radius}$ above the caudex crown.
* **Boundary-Continuity Lobe Rejection (`check_boundary_tissue_continuity`):** I inspect a 2-pixel perimeter in CIE $L^*a^*b^*$ and reject candidates when $> 15\%$ non-paper tissue crosses borders.
* **Dissected-Leaf Rachis Bridging:** In deeply pinnatifid leaves with narrow segmentation breaks, UNet++ v9 `specimen_paper_cleaned_mask` acts as an anatomical constraint (`bridge_dissected_rachis_components`), reconnecting dissected pinnae without bridging across mounting paper.
* **Tier 2 Hemi-Blade Reflection (`reflect_folded_mask`):** For folded blades, I reflect the intact half across the 15-point midvein into `Tier_2_Reflected` for auxiliary linear traits.
* **Classical Gatekeeper and Dissection Gates (`Tier_1_Direct`):** Entire blades require $\text{Solidity} \ge 0.72$, isoperimetric $\text{Circularity} = \frac{4\pi A}{P^2} \ge 0.25$, $\text{Extent} = \frac{A}{W \cdot H} \ge 0.35$, and $\text{UCS} \ge 0.85$. Pinnatifid ($0.48 \le \text{Solidity} < 0.72$) and bipinnatisect ($0.18 \le \text{Solidity} < 0.48$, $\ge 4$ deep sinuses) blades are validated through `is_botanical_dissection()` and routed into `Tier_1_Direct`.
* **Cohort-Adaptive Physical Blade Area Floor:**
  $$\text{effective\_min\_area\_cm2} = \min(0.75, \max(0.15, 0.25 \times \text{median\_sheet\_leaf\_area\_cm2}))$$
  I enforce this floor to prevent dwarf or alpine specimens from being discarded.
* **Post-QC Basal Masking (`crop_dense_rosette_tier4`):** I replace non-basal pixels outside the 5-pixel dilated union of verified basal masks with neutral paper `[242, 239, 230]`.

By filtering out cauline leaves and mounting tape artifacts at this stage, I guarantee that only verified, biologically homologous basal leaves enter geometric shape analysis.

---

## 4. Phase 3 and 4 Protocol: Morphometrics, EFA, Allometry and Ordination (`03` and `04`)

In Phase 3 and Phase 4, I transform discrete leaf outlines into continuous mathematical functions. At first glance, elliptic Fourier analysis and landmark morphometrics seem fundamentally different, but that difference is superficial because both approaches quantify boundary variance relative to biological homology.

### 4.1 Shape Quantification and Strict `Tier_1_Direct` EFA (`03_fourier_extractor.R`)
* **Momocs Outline Export and Euler Characteristic Transform (ECT):**
  LeafMachine3 exports native tip-up, base-anchored clockwise outlines (`momocs_outlines.json`) paired with 360-direction **Euler Characteristic Transform (ECT)** coordinates (`Coordinates/*.h5`, Ayub et al. 2026; Weaver 2026).
* **Strict `Tier_1_Direct` Closed-Contour Restriction:**
  I restrict `Momocs::efourier` strictly to intact `Tier_1_Direct` contours. I never pool folded or reflected `Tier_2_Reflected` leaves into EFA. I ingest contours through batch SQL from `extracted_leaves.contour_xy_json`.
* **Inertial Axis Pre-Alignment and Harmonic Power Calibration:**
  I align contours along their major inertial axis via `Momocs::coo_alignxax()` with a circularity $> 0.75$ suborbicular guard. Harmonic power:
  $$\text{Power}_n = \frac{A_n^2 + B_n^2 + C_n^2 + D_n^2}{2},\quad \frac{\sum_{n=1}^{H_{\text{opt}}} \text{Power}_n}{\sum_{n=1}^{H_{\max}} \text{Power}_n} \ge 0.99$$
  I dynamically determine the harmonic count $H_{\text{opt}} \in [6, 24]$ capturing $\ge 99\%$ cumulative power, executing `Momocs::efourier(out_obj, nb.h = H_opt, norm = FALSE, start = FALSE)`.
* **Uniform Bilateral Symmetrization (Claude 2008):**
  $$A_n^{\text{sym}} = A_n,\quad D_n^{\text{sym}} = D_n,\quad B_n^{\text{sym}} = 0,\quad C_n^{\text{sym}} = 0,\quad \text{pct\_asymmetry} = \frac{\sum_{n=1}^{H_{\text{opt}}} (B_n^2 + C_n^2)}{\sum_{n=1}^{H_{\text{opt}}} (A_n^2 + B_n^2 + C_n^2 + D_n^2)}$$
* **Type Outlier Immunity and OTU Median Aggregation:**
  I exempt intact type contours from within-taxon Mahalanobis and MCD outlier screening (`mcd_out[is_type_specimen] <- FALSE`). Active type anchors receive a $1.5\times$ precision boost in `stats::cov.wt()`. I export PCs accounting for $\ge 95\%$ variance ($k_{\text{export}} = \max(5, k_{95})$), followed by unweighted coordinate-wise OTU median aggregation to `(catalogNumber, plant_individual_id)`.

### 4.2 Determiner Authority Weighting and Adaptive Reference Anchors (`authority_weighting_utils.R`)
* **Adaptive Reference Anchors and Type-Medoid Shrinkage:**
  Reference anchors consist of `Tier_1_Gold` and `Tier_2_Silver` vouchers pruned once at $\chi^2_{0.95}$ FastMCD. Intact types receive retention immunity (`anchor_tier_source = "Tier_0_Type_Anchor"`), whereas extreme types ($> \chi^2_{0.99}$) are isolated from FastMCD scatter estimation. Taxon reference centers shrink toward the type medoid:
  $$\boldsymbol{\mu}_{\text{ref}, k} \leftarrow (1 - \gamma_k) \boldsymbol{\mu}_{\text{ref}, k} + \gamma_k \boldsymbol{\mu}_{\text{type}, k},\quad \gamma_k = \min\!\left(0.35,\, \frac{2.5\,n_{\text{eff, type}, k}}{N_{\text{gs}, k} + 2.5\,n_{\text{eff, type}, k}}\right)$$
  Taxa with $< 8$ Gold and Silver anchors promote central `Tier_3_Bronze` vouchers within $\chi^2_{0.80}$ (`is_promoted_anchor = TRUE`).
* **2-Step Dawid-Skene and Beta-Binomial Empirical Bayes Authority Weighting:**
  I estimate observer error rates $\hat{\epsilon}_i$ exogenously (`target_modality = ...`; Dawid & Skene 1979) and pool priors ($\alpha_0 = 5.0$, prior reliabilities $\pi_{\text{PrimaryType}} = 0.999$, $\pi_{\text{SecType}} = 0.985$, $\pi_{\text{Gold}} = 0.95$, $\pi_{\text{Silver}} = 0.85$, $\pi_{\text{PromotedBronze}} = 0.65$, $\pi_{\text{Bronze}} = 0.35$). Types are exempt from downward shrinkage. Unpromoted Bronze receive $0.0$ in `active_passive` mode.

### 4.3 Klingenberg Allometry, Dual-Track CDA and MANOVA (`04_gmm_morphotools.R`)
* **Two-Sample Prediction-Variance Corrected Type Distance:**
  $$D^2_{\text{type}, i} = \left(1 + \frac{1}{n_{\text{eff, type}, k}}\right)^{-1} (\mathbf{x}_i - \boldsymbol{\mu}_{\text{type}, k})^\top \mathbf{S}_{\text{reg}}^{-1} (\mathbf{x}_i - \boldsymbol{\mu}_{\text{type}, k})$$
  persisted as `type_morph_mahal_d2`, alongside `nearest_type_morph_taxon` and `is_typical_morphotype` (GMM cluster $c^*$).
* **Klingenberg (2016) Allometric Slope Homogeneity:**
  I evaluate slope homogeneity across ruler-calibrated anchors (`is_calibrated_mm == TRUE`) through multivariate interaction models in `RRPP::lm.rrpp(..., weights = eff_weight)`:
  $$\mathbf{y}_{\text{shape}} = \beta_0 + \beta_1 \log_{10}(\text{CS}) + \boldsymbol{\gamma}_{\text{species}} + \boldsymbol{\delta} (\log_{10}(\text{CS}) \times \boldsymbol{\gamma}_{\text{species}}) + \boldsymbol{\epsilon}$$
  If $p \ge 0.01$, I apply common-slope correction: $\mathbf{y}_{\text{corrected}} = \mathbf{y}_{\text{shape}} - \hat{\beta}_1 (\log_{10}(\text{CS}) - \overline{\log_{10}(\text{CS})})$.
* **Subspace Selection and Dual-Track CDA:**
  I apply Horn's Parallel Analysis capped by Williams & Titus (1988) $P_{\max} = \lfloor (N_{\min,\text{anchor}} - 1) / 3 \rfloor$. I evaluate honest LOOCV on weighted Ledoit-Wolf shrinkage CDA. Track 1 evaluates Unsupervised GMM Morphotypes (`Morphotype_K1..Kn`) via Kass & Raftery (1995) $2\Delta\text{BIC}$ thresholds ($> 6$ strong, $> 10$ very strong). Track 2 evaluates Supervised EPA Level III Ecoregion CDA biplots and MANOVA with Benjamini-Hochberg FDR.

---

## 5. Phase 5 Protocol: Multimodal Spatial RF, *Locus Classicus* Ecological Control and SDM Ensembles (`05_multimodal_spatial_rf.R`)

Historical collections frequently lack fine-scale GPS coordinates, as eighteenth and nineteenth-century collectors recorded only general county or river valley descriptions. If an ecological model extracts raster values directly at an uncertain county centroid, the extracted soil variables become distorted.

* **Type Locality (*Locus Classicus*) Ecological Baseline:**
  I gate direct point-raster extraction on types by coordinate uncertainty ($\le 2000\text{ m}$ edaphic and topographic, $\le 10000\text{ m}$ climatic; Daru et al. 2018). Coarse type sheets derive their *locus classicus* edaphic profiles from proximate Gold vouchers within 150 km or sharing the EPA Level III ecoregion. I calculate physiological divergence through weighted RMS distance:
  $$z_{\text{type}, i} = \sqrt{\frac{\sum_{j \in \text{Primary}} w_j \left(\frac{x_{i, j} - \mu_{\text{type}, k, j}}{\max(\text{MAD}_{\text{type}, k, j},\, s_{\text{floor}, j},\, 10^{-4})}\right)^2}{\sum_{j \in \text{Primary}} w_j}}$$
  exported as `type_edaphic_z_score`, alongside `nearest_type_niche_taxon` and boolean `matches_type_lithology`.
* **Full 19-Bioclimatic and 0 to 15 cm Rhizosphere Stack:**
  I extract 19 WorldClim v2.1 bioclimatic variables (`bio1` to `bio19`) and 0 to 15 cm rhizosphere depth-weighting (`(5 * v_0_5 + 10 * v_5_15) / 15`) across SSURGO, 30 m POLARIS, 100 m SOLUS100, 3DEP topography, and 6-class USGS SGMC lithology (`edaphic_fine_grain_utils.R`).
* **3-Stage Dynamic Primary Predictor Gating:**
  First, I prune collinear predictors ($|r| < 0.70$, $\text{VIF} < 5.0$) across Climatic, Edaphic, and Topographic blocks. Second, I train `ranger::ranger` under $0.5^\circ \times 0.5^\circ$ Spatial Block CV with continuous Empirical Bayes weights (primary types pinned to 3.0 ceiling) and winsorized case weights `[0.25, 4.0]`, isolating variables accounting for the top 85% cumulative importance. Third, I restrict downstream niche distances strictly to discoverable primary predictors.
* **SDM Ensembles (MaxEnt and Spatial RF):**
  I fit authority-weighted MaxEnt (`maxnet`, `cloglog` link; `step05_sdm_maxent.R`) and Spatial RF ensembles with target-group background selection, exporting spatial-block $AUC_{\text{block}}$ and Continuous Boyce Index (CBI) to reports and GeoTIFFs (`fig05`, `fig05b`).
* **Hopkins Circular Phenology (`step05_phenology_ecotypes.R`):**
  I regress circular anthesis DOY ($\theta_i = \frac{2\pi \cdot \text{DOY}_i}{365.25}$) on latitude, elevation, and temperature, scaled by $\alpha = \min(1.0, \max(0.0, (N_{\text{anthesis}} - 5)/15))$. Flowering dispersion scales thresholds ($\tau_{\text{cryptic}} = \min(35, \max(14, 2.5\,\sigma_{\text{DOY}}))$, $\tau_{\text{ecotype}} = \min(28, \max(10, 2.0\,\sigma_{\text{DOY}}))$).

---

## 6. Phase 6 Protocol: Reliability-Weighted Consensus and Hierarchical Nomenclatural Type Decision Matrix (`06_triage_dashboard_synthesis.R`)

In Phase 6, I synthesize all lines of evidence. If for say, a specimen displays an anomalous leaf shape but matches its type locality in soils and phenology, a single-modality model would misclassify it as a new species. By synthesizing across independent data streams, I resolve genuine taxonomic status.

* **Master View Synthesis and Reliability-Weighted Consensus:**
  I synthesize consensus by querying `v_integrative_triage_master`, which joins all modalities onto `plant_individuals`. I weight modalities by out-of-fold Cohen's $\kappa_m$ and class-conditional precision $C^{(m)}_{k,j}$:
  $$w_{i,m} = w_m \cdot (0.5 + 0.5 \cdot \text{Precision}_m(\hat{y}_{i,m}))$$
  achieving high voucher completion across georeferenced and non-georeferenced sheets.
* **Taxon-Agnostic 5-Stage Hierarchical Decision Matrix:**
  - **Rule A (Nomenclatural Type Override):** I never reassign type specimens (`recommended_determination = given_sp`), nor route them to `Misidentification_Reassign`, `Hybrid_Introgressant`, or `Putative_Cryptic_Species`. Concordant types receive `taxonomic_status_call = "Confirmed_Type_Control"` (`RESOLVED`). Cross-modal divergence triggers `taxonomic_status_call = "Type_Synonymy_Candidate"` (`CRITICAL`).
  - **Rule B (Type-Control Reassignment):** When non-types deviate from type baselines ($D^2_{\text{type}} > \chi^2_{0.95}$, $z_{\text{type}} > 2.0$) with $\ge 2$ type streams agreeing on $T_{\text{alt}} \ne \text{given\_sp}$ without classifier contradiction, I reassign them to $T_{\text{alt}}$.
  - **Rule C (*Sensu Stricto* Typical Morphotype Locking):** Specimens in typical morphotype GMM cluster ($c^*$) with $z_{\text{type}} \le 2.0$ are locked as `Confirmed_Species` (*sensu stricto*).
  - **Stage 2 (`Hybrid_Introgressant`):** Intermediate CDA posteriors ($0.35 \le P_1 \le 0.65$, $P_2 \ge 0.25$) in sympatry, resolving parents via `cda_second_taxon`.
  - **Stage 3 (`Putative_Cryptic_Species`):** Requires `morph_tier_source == "Tier_1_Pure"`, GMM $2\Delta\text{BIC} > 6$, alongside independent edaphic or phenological divergence.
  - **Stage 4 (`Subspecies_Ecotype` vs. `Ecophenotypic_Plasticity`):** Evaluates cohesive geographic or edaphic divergence for ecotypes versus overlapping envelopes for plasticity.
* **Synthesis Figures and Curatorial Reports:**
  I generate integrative multi-panel plates (`pub_synthesis_plates.R`; `fig06`) and produce actionable curation summaries (`Packera_Taxonomic_Revision_Summary.md`, `curation_triage`).

### Integrative Taxonomic Decision Matrix

| Taxonomic Status Call | Vegetative Morphometrics | Reproductive Anthesis & Scape | Cytology ($2n$) | nextRAD Phylogenomics | Rhizosphere Edaphics & Niche | Actionable Outcome |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`Confirmed_Type_Control`** | Type centroid ($\boldsymbol{\mu}_{\text{type}}$); `is_typical_morphotype` | Concordant with type phenology | Reference cytotype | Type and paratype genetic profile | High fidelity to *locus classicus* | Immutable nomenclatural anchor; retain determination (`RESOLVED`) |
| **`Type_Synonymy_Candidate`** | Projects into congener morphospace ($P \ge 0.85$) | Overlaps congener flowering schedule | Shared or intermediate cytotype | Clusters with congener lineages | Overlaps congener ecological envelope | Critical priority: evaluate formal taxonomic synonymy (`CRITICAL`) |
| **`Putative_Cryptic_Species`** | Discrete GMM ($2\Delta\text{BIC} > 6$); CDA $P < 0.70$ | Divergent flowering window or scape height | Distinct or polyploid ($2n = 44, 46, 88$) | Diagnostic monophyletic cluster | Significant edaphic divergence or allochrony | Formal species-level recognition or resurrection under USC |
| **`Misidentification_Reassign`** | High posterior for known congener; Rule B | Concordant with reassigned congener | Matches reassigned congener | Clusters within reassigned congener | Matches edaphic and climatic envelope of congener | Annotate and reassign determination label |
| **`Hybrid_Introgressant`** | Intermediate CDA ($0.35\text{ to }0.65$, $P_2 \ge 0.25$) | Intermediate flowering timing | Aneuploid, triploid, or mixed | Admixed ancestry; reticulation | Sympatric contact zones or ecotones | Designation as nothospecies or hybrid swarm |
| **`Subspecies_Ecotype`** | Partial shift ($2\Delta\text{BIC} \le 6$); conserved outline | Conserved flowering window | Shared cytotype with nominal | Shallow regional structure | Specialized substrate ($D_{\text{edaph}} < D_{\text{clim}}$) | Infraspecific or varietal recognition |
| **`Ecophenotypic_Plasticity`** | Unimodal GMM ($K=1$); parallel slope | Concordant with given taxon | Identical cytotype | Panmictic with parental populations | Overlapping climatic and edaphic envelopes | Retain under polymorphic species |

---

## 7. Phase 7 Protocol: Unsupervised Quantitative Taxon-Complex Discovery (`07_discover_taxon_complexes.R`)

In Phase 7, I scale beyond single-species comparisons to delimit species complexes across the genus. At first glance, local ecotypes and distinct species complexes seem to require different analytical approaches, but that difference is superficial because both represent continuous divergence surfaces across morphospace and ecological landscapes.

* **Multimodal Pairwise Affinity Calculation:**
  I calculate pairwise morphological, fine-scale edaphic, climatic, and phenological distances across all active taxa, shrinking taxon centroids toward nomenclatural type baselines with weight $\gamma_k$.
* **Unsupervised Graph Modularity and Hierarchical Delimitation:**
  I discover species complexes through unsupervised graph modularity (Leiden and Louvain algorithms) and Ward.D2 hierarchical clustering, resolving Taxon by GMM Morphotype subclusters when $N_{\text{taxa}} \le 2$.
* **Complex Discovery Plates (`pub_synthesis_plates.R`; `fig07`):**
  I export affinity networks and hierarchical dendrograms alongside pairwise distance summaries (`empirical_taxon_complexes.csv`, `taxon_complex_pairwise_affinity.csv`), providing quantitative evidence for species complex boundaries.

---

## 8. Literature Cited

1. Ayub, M. S., et al. 2026. Topological shape descriptors and ECT in plant phenomics. *bioRxiv* https://doi.org/10.64898/2026.04.13.718293
2. Barkley, T. M. 1988. Variation among Senecioneae in North America. *Brittonia* 40: 211 to 221.
3. Broennimann, O., et al. 2012. Measuring ecological niche overlap. *Glob Ecol Biogeogr* 21: 481 to 497.
4. Claude, J. 2008. *Morphometrics with R*. Springer, New York.
5. Collyer, M. L., & D. C. Adams. 2018. RRPP: Linear models via residual randomization. *Methods Ecol Evol* 9: 1772 to 1779.
6. Daru, B. H., et al. 2018. Widespread sampling biases in herbaria. *New Phytol* 217: 939 to 955.
7. Dawid, A. P., & A. M. Skene. 1979. Maximum likelihood estimation of observer error-rates. *Appl Stat* 28: 20 to 28.
8. de Queiroz, K. 2007. Species concepts and species delimitation. *Syst Biol* 56: 879 to 886.
9. Hirzel, A. H., et al. 2006. Evaluating ability of habitat suitability models to predict presence. *Ecol Modell* 199: 142 to 152.
10. Kass, R. E., & A. E. Raftery. 1995. Bayes Factors. *JASA* 90: 773 to 795.
11. Klingenberg, C. P. 2016. Size, shape, and form: concepts of allometry in geometric morphometrics. *Dev Genes Evol* 226: 113 to 137.
12. Kowal, R. R. 1975. Systematics of *Senecio aureus* and allied species on Gaspé Peninsula. *Mem Torrey Bot Club* 23: 1 to 113.
13. Kuhl, F. P., & C. R. Giardina. 1982. Elliptic Fourier features of a closed contour. *Comput Graph Image Process* 18: 236 to 258.
14. Omernik, J. M., & G. E. Griffith. 2014. Ecoregions of conterminous United States. *Environ Manage* 54: 1249 to 1266.
15. Padial, J. M., et al. 2010. The integrative future of taxonomy. *Front Zool* 7: 16.
16. Phillips, S. J., et al. 2009. Sample selection bias in Maxent. *Ecography* 32: 181 to 197.
17. Phillips, S. J., et al. 2017. Opening the black box: open-source release of Maxent. *Ecography* 40: 887 to 893.
18. Scrucca, L., & A. E. Raftery. 2018. clustvarsel: Variable selection for model-based clustering. *J Stat Softw* 84: 1 to 28.
19. Šlenker, M., et al. 2022. MorphoTools2: multivariate morphometric analysis. *Bioinformatics* 38: 2954 to 2955.
20. Trock, D. K. 2006. *Packera*. In *Flora of North America*, Vol. 20, 570 to 602. Oxford University Press, New York.
21. Turland, N. J., et al. 2018. *International Code of Nomenclature (Shenzhen Code)*. Regnum Veg 159. Koeltz, Glashütten.
22. Valavi, R., et al. 2019. blockCV: Spatially separated folds for k-fold cross-validation. *Methods Ecol Evol* 10: 225 to 232.
23. Warren, D. L., et al. 2008. Environmental niche equivalency versus conservatism. *Evolution* 62: 2868 to 2883.
24. Weakley, A. S. 2026. *Flora of the Southeastern United States*. UNC Herbarium (NCU), Chapel Hill.
25. Weaver, W. N. 2026. LeafMachine3: Automated machine vision for herbarium trait phenotyping. https://github.com/Gene-Weaver/LeafMachine3
26. Williams, B. K., & K. Titus. 1988. Sampling stability in discriminant analysis. *Ecology* 69: 1275 to 1285.
27. Wright, M. N., & A. Ziegler. 2017. ranger: Fast implementation of random forests. *J Stat Softw* 77: 1 to 17.
