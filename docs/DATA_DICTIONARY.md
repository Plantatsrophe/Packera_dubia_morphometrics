# Relational Schema and Output Data Dictionary (`PLANTMATICS_SCHEMA_USER_VERSION = 20261127`)

**Standard Operating Procedure:** [`docs/TAXONOMIC_SOP.md`](TAXONOMIC_SOP.md) (`UNC-BOT-SOP-2026-04-REV5`)
**Authoritative Schema Definitions:** `scripts/core/db_schema.py` and `scripts/core/db_utils.R`

## 1. Architectural Conventions and Key Invariants

Herbarium sheets frequently carry multiple pressed plants (`ramets`) mounted under a single barcode, whereas repeated exploratory runs over different species subsets or authority filters must never overwrite raw computer vision extractions. I therefore split pipeline storage across two transactional tiers inside `data/plantmatics.sqlite` (Write-Ahead Logging mode, `PRAGMA user_version = 20261127`, supporting automatic zero-loss atomic rename migration from legacy `data/packera_pipeline.sqlite` with `-wal` and `-shm` sidecars) paired with atomic CSV snapshots on disk. Environment overrides use `PLANTMATICS_DB_PATH`, `PLANTMATICS_RUN_DIR`, `PLANTMATICS_GENOMICS_DIR`, `PLANTMATICS_CYTOLOGY_DIR`, and `PLANTMATICS_TAXONOMY_DIR` while preserving `PACKERA_*` fallbacks.

Tier A (`data/`) stores invariant specimen metadata, full-sheet WebP scans, segmented leaf masks, 120-vertex contours, $32 \times 360$ Euler Characteristic Transform (ECT) matrices, spatially cached environmental plus USGS `DS-801` geochemical and USFWS NWI wetland extractions, flat `data/taxonomy/` CoL WCVP `2232` alongside BHL API v3 and NatureServe Explorer taxonomic backbones (`taxon_backbone_profiles` and `natureserve_subnational_ranges`), staged genomic sample metadata paired with `HDplot`/`Compositae-ParaLoss-1272` locus weights and inferred cytotype records (`genomic_samples`), plus curated public and voucher-linked cytological registries (`cytological_records` and `taxon_cytology_summary`). Tier B (`runs/<RUN_ID>/` and `runs/latest`) isolates run-specific cohort memberships, Elliptic Fourier harmonics, canonical Gaussian Mixture Models, Ledoit-Wolf shrinkage discriminant scores, Spatial Random Forest posteriors, Step 06 phylogenomic and coalescent delimitation metrics (`genetic_lineage_assignments`), and five-stage curatorial triage calls.

Three cross-language invariants govern every column documented below:

* **Composite Individual Key `(run_id, catalogNumber, plant_individual_id)`:** A single barcode (`catalogNumber`) can yield multiple DBSCAN-clustered individuals (`plant_individual_id = 1, 2, ...`, or `0` when no whole-plant cluster is resolved). Every Python and R join across Steps 02 through 07 groups strictly on `(catalogNumber, plant_individual_id)` so multi-plant sheets never collapse to a single row, whereas `v_integrative_triage_master` deduplicates replicate genomic samples on `(run_id, catalogNumber)` via `ROW_NUMBER()` so replicate genetic libraries never inflate specimen row counts.
* **Native Pixel vs. Calibrated Millimeter Coordinates:** All `_px` columns represent coordinates on the original unscaled full-sheet scan (`coord_work / work_scale`). Columns ending in `_mm`, `_cm`, or `_mm2` contain physical measurements only when `is_calibrated_mm = 1` (`scale_source` in `{'ruler_calibrated', 'fieldprism_calibrated'}`), whereas uncalibrated sheets (`mp_estimated` or `uncalibrated`) populate dimensionless ratios (`aspect_ratio`, `solidity`, `petiole_to_blade_ratio`, `AR_inv`) and unit-centroid-size EFA/ECT descriptors while leaving physical linear columns honestly `NULL` (`NA`) in strict scale-gated models.
* **Honest Missingness Outside Demo Mode:** When `--demo` is inactive, missing environmental rasters, high-uncertainty coordinates exceeding layer gates, or absent floral heads record `NULL` in SQLite and `NA` in exported CSVs rather than imputed placeholders. Running with `--demo` (`PLANTMATICS_DEMO_MODE=1`, with `PACKERA_DEMO_MODE=1` compatibility) redirects all database writes to an ephemeral SQLite file (`demo_ephemeral.sqlite`) so canonical Tier A tables remain untouched.

## 2. Tier A Shared Canonical SQLite Tables (`data/plantmatics.sqlite`)

### 2.1 `taxon_synonyms` (Step 01 Dynamic Basionym Cache)

**Primary Key:** `raw_name_lower` | **Populated By:** `scripts/core/harvester.py`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `raw_name_lower` | `TEXT` | Lowercase binomial/trinomial | Normalized raw taxon string queried from GBIF, iDigBio, or historical *Senecio* and *Cineraria* labels. |
| `canonical_packera` | `TEXT` | *Packera* taxon string | Accepted *Packera* binomial or infraspecific taxon (`subsp.` / `var.`) resolved via GBIF Backbone and Weakley (2026). |
| `protologue_authors_json` | `TEXT` | JSON string array | Optional legacy metadata array retained for schema compatibility; Step 01 never promotes non-type sheets from author-string matches. |
| `updated_at` | `TEXT` | ISO-8601 UTC timestamp | Timestamp of the most recent taxonomic backbone verification. |

### 2.2 `usda_county_ranges` (Step 01 County Distribution Cache)

**Primary Key:** `(canonical_packera, state_norm, county_norm)` | **Populated By:** `scripts/core/usda_plants_ranges.py`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `canonical_packera` | `TEXT` | *Packera* taxon string | Standardized *Packera* species or variety name. |
| `usda_symbol` | `TEXT` | USDA PLANTS code | Accepted USDA PLANTS symbol (for example `PADU` for *Packera dubia*). |
| `state_norm` | `TEXT` | 2-letter postal code | Standardized US state or Canadian province abbreviation. |
| `county_norm` | `TEXT` | Normalized county name | Lowercase county or parish name stripped of administrative suffixes (`""` when `is_state_only = 1`). |
| `fips_code` | `TEXT` | 5-digit FIPS string | US county FIPS identifier when available. |
| `native_status` | `TEXT` | Status code | USDA PLANTS nativity string (`N`, `I`, or `W`). |
| `is_state_only` | `INTEGER` | `{0, 1}` | Equals `1` when USDA PLANTS records presence only at the state level without county breakdown. |
| `updated_at` | `TEXT` | ISO-8601 UTC timestamp | Cache refresh timestamp. |

### 2.3 `vouchers` (Step 01 Canonical Specimen Registry)

**Primary Key:** `catalogNumber` | **Populated By:** `scripts/data_prep/01_voucher_harvester.py` and Step 02 calibration updates
**Companion CSV:** `curated_vouchers.csv` exports the 25 core Darwin Core, authority, ecoregion, and duplicate-audit columns from this table for R consumption.

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `catalogNumber` | `TEXT` | Herbarium barcode / UUID | Unique specimen identifier (`institutionCode` + barcode or occurrence UUID fallback). |
| `institutionCode` | `TEXT` | Index Herbariorum code | Source herbarium code (for example `NCU`, `US`, `NY`, `MO`, `GH`, `DUKE`). |
| `scientificName` | `TEXT` | Darwin Core string | Raw verbatim scientific name with authorship as served by GBIF or iDigBio. |
| `species_raw` | `TEXT` | Binomial/trinomial | Parsed raw taxon name prior to *Senecio*/*Cineraria* synonym mapping. |
| `species_standardized` | `TEXT` | *Packera* taxon string | Accepted *Packera* binomial or infraspecific rank (`subsp.` / `var.`) used as the label hypothesis $y_i$. |
| `infraspecific_taxon` | `TEXT` | Infraspecific epithet | Preserved `subsp.`, `var.`, or `f.` rank token when present on the label. |
| `recordedBy` | `TEXT` | Collector string | Primary collector name(s); isolated strictly from `identifiedBy` unless label explicitly states `det. by collector`. |
| `recordNumber` | `TEXT` | Collector number | Alphanumeric gathering number paired with `recordedBy` and `eventDate` to form `exsiccatae_key`. |
| `identifiedBy` | `TEXT` | Determiner string | Taxonomist who annotated the determination slip. |
| `determiner_tier` | `TEXT` | Authority tier enum | Stratified reliability tier: `Tier_0_Type_Anchor`, `Tier_1_Gold`, `Tier_2_Silver`, `Tier_3_Bronze`, or `Tier_4_Unverified`. |
| `type_status` | `TEXT` | Nomenclatural type enum | Standardized ICN type rank (`Holotype`, `Lectotype`, `Neotype`, `Isotype`, `Syntype`, `Paratype`, `Epitype`, or `""`). |
| `county` | `TEXT` | County / parish | Administrative county parsed from Darwin Core metadata. |
| `stateProvince` | `TEXT` | State / province | Standardized state or province name. |
| `locality` | `TEXT` | Free text | Verbatim collection locality description used for outcrop and habitat parsing. |
| `decimalLatitude` | `REAL` | Decimal degrees (`WGS84`) | Latitude in `[-90, 90]`; `NULL` when un-georeferenced or scrubbed by centroid QC. |
| `decimalLongitude` | `REAL` | Decimal degrees (`WGS84`) | Longitude in `[-180, 180]`; `NULL` when un-georeferenced or scrubbed by centroid QC. |
| `coordinateUncertainty` | `REAL` | Meters (`m`) | Reported horizontal georeferencing uncertainty (`coordinateUncertaintyInMeters`). |
| `georef_qc_status` | `TEXT` | Georeference QC enum | Coordinate audit tag (`PASS`, `UN_GEOREFERENCED_RESCUED_TYPE`, `COUNTY_CENTROID_SCRUBBED`, `HIGH_UNCERTAINTY`). |
| `regional_group` | `TEXT` | US EPA Level III name | Assigned EPA Level III Ecoregion from point-in-polygon lookup (or $\ge 95\%$ single-ecoregion county fallback). |
| `spatial_balance_weight` | `REAL` | Dimensionless $(0, 3]$ | Set to neutral `1.0` during lossless Step 01 Tier A ingestion, then recomputed downstream across non-quarantined active anchors in Steps 03 to 05 and Tier B preflight. |
| `year`, `month`, `day` | `INTEGER` | Calendar integers | Parsed collection year, month (`1` to `12`), and day (`1` to `31`). |
| `eventDate` | `TEXT` | ISO-8601 date | Standardized collection date (`YYYY-MM-DD`). |
| `doy` | `REAL` | Day of Year (`1.0` to `366.0`) | Ordinal collection day; set to `NULL` when month/day are absent or represent `YYYY-01-01` placeholder artifacts. |
| `image_path` | `TEXT` | Relative file path | Local path to the lossless/high-Q full-sheet WebP scan under `data/raw_vouchers/`. |
| `exsiccatae_key` | `TEXT` | Hash / normalized key | Gathering deduplication key combining normalized primary collector surname, `recordNumber`, and `eventDate`. |
| `is_primary_duplicate` | `INTEGER` | `{0, 1}` | Equals `1` for the highest-resolution / highest-authority sheet retained from a duplicate gathering set. |
| `duplicate_count` | `INTEGER` | Count $\ge 1$ | Total duplicate sheets detected across herbaria for `exsiccatae_key`. |
| `exsiccatae_taxonomic_conflict` | `INTEGER` | `{0, 1}` | Equals `1` when duplicate sheets of the exact same gathering carry conflicting *Packera* determinations across herbaria. |
| `exsiccatae_gold_taxon` | `TEXT` | *Packera* taxon string | Non-mutating diagnostic `Tier_1_Gold` or `Tier_0_Type_Anchor` determination recorded alongside `exsiccatae_taxonomic_conflict` without altering conflicting sheet labels. |
| `usda_symbol` | `TEXT` | USDA PLANTS code | Matched USDA PLANTS symbol for `species_standardized`. |
| `usda_county_status` | `TEXT` | Range audit enum | `in_range`, `buffered_50km`, `out_of_range`, `state_only_match`, or `unverified` (`50 km` county border tolerance; audit-only by default). |
| `pixels_per_mm` | `REAL` | $\text{px}\cdot\text{mm}^{-1}$ | Sheet optical resolution updated by Step 02 ruler tick-lattice or FieldPrism calibration. |
| `is_calibrated_mm` | `INTEGER` | `{0, 1}` | Equals `1` when `pixels_per_mm` derives from physical ruler or FieldPrism ticks rather than megapixel estimation. |
| `qc_status` | `TEXT` | `CHECK(qc_status IN ('PASS', 'FAILED_QC', 'PENDING'))` | Step 02 specimen processing status (`PENDING` after Step 01 harvest; `PASS` or `FAILED_QC` after Step 02). Sole SQLite `CHECK` constraint across the schema; Step 02 Tier A and Tier B tables enforce all domain checks explicitly in Python and R. |
| `failure_reason` | `TEXT` | Audit code | Concise categorical failure code when `qc_status = 'FAILED_QC'` (for example `NO_VALID_LEAVES` or `CORRUPT_IMAGE`). |
| `failure_details` | `TEXT` | Free text | Diagnostic details accompanying `failure_reason`. |

### 2.4 `plant_individuals` (Step 02 Multi-Plant Ramet Architecture)

**Primary Key:** `(catalogNumber, plant_individual_id)` | **Populated By:** `scripts/vision/organ_crown_metrics.py`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `catalogNumber` | `TEXT` | Herbarium barcode | Foreign key to `vouchers.catalogNumber`. |
| `plant_individual_id` | `INTEGER` | Integer $\ge 0$ | DBSCAN individual index on the sheet (`1..m` for clustered plants; `0` fallback). |
| `rosette_centroid_x`, `rosette_centroid_y` | `REAL` | Unscaled pixels (`px`) | Full-sheet coordinates of the basal rosette centroid. |
| `crown_x`, `crown_y` | `REAL` | Unscaled pixels (`px`) | Estimated root-crown / caudex attachment coordinates anchoring basal vs. cauline leaf zoning. |
| `plant_height_px` | `REAL` | Unscaled pixels (`px`) | Geodesic stem/scape length from caudex crown to terminal inflorescence apex. |
| `plant_height_mm`, `plant_height_cm` | `REAL` | `mm` / `cm` | Calibrated plant height (`NULL` when `is_calibrated_mm = 0`). |
| `is_bent_stem` | `INTEGER` | `{0, 1}` | Equals `1` when the mounted stem is folded (`V` or `Z` bend) and measured piecewise along skeleton segments. |
| `petiole_length_px`, `petiole_length_mm` | `REAL` | `px` / `mm` | Individual-level $P_{75}$ upper-quartile basal petiole length resilient to partial overlapping. |
| `petiole_width_px`, `petiole_width_mm` | `REAL` | `px` / `mm` | Median basal petiole width measured orthogonal to the petiole skeleton. |
| `blade_length_px`, `blade_length_mm` | `REAL` | `px` / `mm` | Median basal lamina length from `lamina_base` to `lamina_tip`. |
| `petiole_to_blade_ratio` | `REAL` | Dimensionless $\ge 0$ | Scale-invariant ratio `petiole_length_px / blade_length_px` diagnostic across *Packera dubia* and North American congeners. |
| `scape_to_basal_leaf_ratio` | `REAL` | Dimensionless $\ge 0$ | Scale-invariant ratio of scape height to total basal leaf length (`petiole + blade`). |
| `rosette_patch_path`, `rosette_mask_path` | `TEXT` | Relative file path | Optional cropped rosette RGB image and binary mask paths. |

### 2.5 `extracted_leaves` (Step 02 Read-Only LeafMachine3 Organ Table)

**Primary Key:** `(catalogNumber, plant_individual_id, leaf_id)` | **Populated By:** `scripts/pipeline/02_segment_and_extract.py`
**Companion CSV:** `extracted_leaf_manifest.csv`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `catalogNumber` | `TEXT` | Herbarium barcode | Foreign key to `vouchers.catalogNumber`. |
| `plant_individual_id` | `INTEGER` | Integer $\ge 0$ | Parent plant individual on the sheet. |
| `leaf_id` | `INTEGER` | Integer $\ge 1$ | Sequential leaf index within `(catalogNumber, plant_individual_id)`. |
| `shape_id` | `TEXT` | Unique organ ID | Composite string `<catalogNumber>_ind<pid>_leaf<lid>` linking SQLite rows to contour and ECT files. |
| `solidity` | `REAL` | Dimensionless $(0, 1]$ | Ratio of hole-subtracted lamina area to its convex hull area; inversely tracks pinnatifid sinus depth. |
| `circularity` | `REAL` | Dimensionless $(0, 1]$ | Isoperimetric quotient $4\pi A / P^2$ reflecting margin serration and lobing complexity. |
| `extent` | `REAL` | Dimensionless $(0, 1]$ | Ratio of lamina area to its oriented bounding-box area. |
| `aspect_ratio` | `REAL` | Dimensionless $> 0$ | Anatomical lamina length-to-width ratio ($L / W$); values $< 1.0$ denote reniform or transversely obovate blades. |
| `area_px`, `area_mm2` | `REAL` | $\text{px}^2$ / $\text{mm}^2$ | Hole-subtracted lamina surface area (`area_mm2` populated when calibrated). |
| `midrib_angle_deg` | `REAL` | Degrees $[-180, 180]$ | Native midrib orientation angle on the herbarium sheet prior to tip-up rotation. |
| `petiole_length_px`, `petiole_length_mm` | `REAL` | `px` / `mm` |Segmented Class-1 petiole arc length attached to this leaf (`0` for sessile cauline leaves). |
| `petiole_width_px`, `petiole_width_mm` | `REAL` | `px` / `mm` | Median width of the attached Class-1 petiole mask. |
| `blade_length_px`, `blade_length_mm` | `REAL` | `px` / `mm` | Keypoint-guided geodesic midrib length from `lamina_base` to `lamina_tip`. |
| `blade_width_px`, `blade_width_mm` | `REAL` | `px` / `mm` | Maximum lamina width orthogonal to the `lamina_base`-to-`lamina_tip` axis. |
| `petiole_to_blade_ratio` | `REAL` | Dimensionless $\ge 0$ | Organ-level petiole-to-lamina length ratio. |
| `pixels_per_mm`, `is_calibrated_mm` | `REAL` / `INTEGER` | $\text{px}\cdot\text{mm}^{-1}$ / `{0, 1}` | Sheet calibration inherited by the organ row. |
| `scale_source` | `TEXT` | Calibration source enum | `ruler_calibrated`, `fieldprism_calibrated`, `tick_lattice_unboxed`, `mp_estimated`, or `uncalibrated` (`tick_lattice_unboxed`, `mp_estimated`, and `uncalibrated` force `is_calibrated_mm = 0` and gate `*_mm` and `centroid_size` outputs to `NULL` in Step 03). |
| `mask_path`, `crop_path` | `TEXT` | Relative file path | Paths to the binary lamina mask (`data/masks/`) and cropped RGB patch. |
| `mask_is_upright` | `INTEGER` | `{0, 1}` | Equals `1` when the saved mask and contour are rotated apex-upward along the midrib vector. |
| `contour_path` | `TEXT` | Relative file path | Path to the 120-vertex clockwise `lamina_base`-anchored contour CSV (`data/contours/`). |
| `ect_path` | `TEXT` | Relative file path | Path to the $32 \times 360$ polar Euler Characteristic Transform `.npy` array (`data/ect/*_ect.npy`). |
| `contour_n_pts` | `INTEGER` | `120` | Vertex count of `contour_xy_json` (`120` arc-length-equidistant vertices). |
| `contour_xy_json` | `TEXT` | JSON `[[x,y], ...]` | Inline 120-vertex tip-up clockwise contour consumed directly by Step 03 in R without filesystem bottlenecks. |
| `leaf_zone` | `TEXT` | `{'basal', 'cauline'}` | Architectural position along the shoot relative to the caudex crown and petiole attachment. |
| `leaf_zone_confidence` | `REAL` | Probability $[0, 1]$ | Confidence score of the `basal` vs. `cauline` classifier. |
| `dist_from_crown_norm` | `REAL` | Dimensionless $[0, 1]$ | Euclidean distance from `lamina_base` to `(crown_x, crown_y)` normalized by `plant_height_px`. |
| `cy_rel_plant_bbox` | `REAL` | Dimensionless $[0, 1]$ | Vertical relative position within the DBSCAN plant bounding box (`0` = base, `1` = apex). |
| `archetype_score` | `REAL` | Score $[0, 1]$ | Native `LeafMachine3` organ completeness score used when selecting visual plate medoids. |
| `bilateral_dice` | `REAL` | Score $[0, 1]$ | Dice overlap between left and right lamina halves reflected across the midrib. |
| `si_a` | `REAL` | Dimensionless $\ge 0$ | Native `LeafMachine3` areal asymmetry index (`0` = perfectly symmetric). |
| `lm3_gates_pass` | `INTEGER` | `{0, 1}` | Native `LeafMachine3` completeness flag audited alongside R-side clump/boundary checks in Step 03. |
| `n_lamina_components` | `INTEGER` | Integer $\ge 1$ | Connected foreground components on the focal `Leaf` mask (`1` = single contiguous outline; `> 1` withheld from EFA and retained for multi-component ECT). |
| `detector_class_id` | `INTEGER` | `{0, 1}` | Stage 1 plant detector organ class (`0` = `Leaf_WHOLE`, `1` = `Leaf_PARTIAL` rescue). |
| `is_crop_clipped` | `INTEGER` | `{0, 1}` | Equals `1` when lamina pixels are truncated flat along a Stage 1 crop box or sheet edge. |
| `hole_area_frac` | `REAL` | Fraction $[0, 1]$ | Fraction of focal blade area occupied by Stage 2 Class-2 insect or herbivory holes (`> 0.035` flags `SEVERE_HERBIVORY_HOLES`). |
| `enclosed_void_frac` | `REAL` | Fraction $[0, 1]$ | Fraction of lamina area enclosed inside unopened inter-lobe overlap loops (`> 0.020` flags `ENCLOSED_INTERLOBE_VOID`, withholding single-loop EFA while retaining ECT). |
| `bilateral_lobe_symmetry` | `INTEGER` | `{0, 1}` | Equals `1` when both left and right lamina halves show balanced multi-lobe inflections and `|si_a| < 0.28`, gating dissected `bilateral_dice` tolerance. |
| `efa_qc_pass` | `INTEGER` | `{0, 1}` | Equals `1` when the leaf passes all Step 02 damage, fold, clipping, symmetry, and single-Jordan-loop EFA QC gates. |
| `qc_flag_reason` | `TEXT` | Semicolon codes | `'PASS'` or semicolon-delimited diagnostic flags (`LEAF_PARTIAL_FRAGMENT`, `CROP_BORDER_CLIPPED`, `SEVERE_HERBIVORY_HOLES`, `LOW_LANDMARK_CONFIDENCE`, `FOLDED_MIDVEIN_TRACE`, `APICAL_TIP_FOLD`, `LONGITUDINAL_MARGIN_CREASE`, `UNILATERAL_ASYMMETRY_OR_TEAR`, `UNRESOLVED_CLUMP`, `ENCLOSED_INTERLOBE_VOID`, `MULTI_COMPONENT_PINNAE`). |
| `blade_cie_l_norm`, `blade_cie_l_sd`, `blade_cie_l_sheet_delta`, `blade_cie_a_delta`, `blade_cie_b_delta` | `REAL` | Standardized CIE $L^*a^*b^*$ units | Eroded-lamina paper-normalized CIE $L^*$ reflectance (`blade_cie_l_norm`, active trait in Step 04), within-blade luminance SD (`blade_cie_l_sd`, held out from GMM/CDA covariance estimation), and paper-relative color/browning deltas (`blade_cie_l_sheet_delta`, `blade_cie_a_delta`, `blade_cie_b_delta`). |

### 2.6 `voucher_reproductive_phenology` (Step 02 Optical Capitulum Staging)

**Primary Key:** `(catalogNumber, plant_individual_id)` | **Populated By:** `scripts/vision/capitulum_phenology_classifier.py`
**Companion CSVs:** `voucher_reproductive_metrics.csv` and `voucher_phenological_states.csv`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `catalogNumber`, `plant_individual_id` | `TEXT` / `INTEGER` | Composite PK | Specimen and individual identifiers. |
| `capitula_count`, `total_capitula` | `INTEGER` | Count $\ge 0$ | Total detected reproductive heads (`flower` + `bud` + `fruit` boxes) assigned to the individual. |
| `involucre_height_px`, `involucre_height_mm` | `REAL` | `px` / `mm` | Median phyllary cup height from receptacle base to phyllary apex. |
| `involucre_width_px`, `involucre_width_mm` | `REAL` | `px` / `mm` | Median phyllary cup width (`NULL` on overly flattened heads failing compression screening). |
| `AR_inv` | `REAL` | Dimensionless $> 0$ | Scale-invariant involucre aspect ratio (`involucre_height_px / involucre_width_px`). |
| `n_anthesis`, `n_bud`, `n_fruiting`, `n_senescent` | `INTEGER` | Counts $\ge 0$ | Head counts per optical developmental stage (`ray`/`disc` ligule chroma vs. expanded pappus bristles). |
| `voucher_phenological_state` | `TEXT` | `{'anthesis','fruit','bud','sterile'}` | Individual-level phenophase resolved under determinate corymbose cyme rules. |
| `is_verified_anthesis` | `INTEGER` | `{0, 1}` | Equals `1` strictly when active ligulate/tubular anthesis is verified; gates Step 05 Hopkins circular phenology. |
| `reproductive_stage_confidence` | `REAL` | Probability $[0, 1]$ | Optical phenophase classification confidence. |
| `dominant_pappus_ratio` | `REAL` | Dimensionless $[0, 1]$ | Fraction of capitulum crown area occupied by bright white capillary pappus bristles. |

### 2.7 `environmental_cache` and `environmental_rtree` (Tier A Spatial Extraction Cache)

**Primary Key:** `lat_lon_key` (`%.5f_%.5f` rounded WGS84 coordinate pair, $\sim 1.1\text{ m}$ precision) paired with virtual R*Tree index `environmental_rtree(id, min_lon, max_lon, min_lat, max_lat)` maintained automatically by `INSERT`, `UPDATE`, and `DELETE` SQLite triggers.
**Populated By:** `scripts/analysis/edaphic_fine_grain_utils.R` alongside `scripts/data_prep/download_geochem_wetland_layers.py`, `scripts/analysis/edaphic_geochem_wetland_utils.R`, and `scripts/analysis/05_multimodal_spatial_rf.R`

| Column Group | Column Names | SQLite Type | Units / Domain | Description |
| :--- | :--- | :---: | :--- | :--- |
| **Keys & Coordinates** | `lat_lon_key`, `catalogNumber`, `latitude`, `longitude` | `TEXT` / `REAL` | `WGS84` deg | Rounded coordinate cache key, specimen catalog identifier, exact decimal latitude, and exact decimal longitude indexed by `environmental_rtree`. |
| **Uncertainty & Eligibility Gates** | `buffer_radius_m`, `high_coord_uncertainty_flag`, `is_training_eligible_coord`, `is_clim_eligible_coord`, `is_edaphic_eligible_coord`, `edaphic_coord_weight` | `REAL` / `INTEGER` | `m` / `{0, 1}` / $[0, 1]$ | General extraction radius (`250 m` default; clamped to `[100, 1000] m` across `<= 5000 m`, whereas `min_bedrock_depth_buffer_cm` caps at `[100, 300] m` and `30 m` microtopography plus `wetland_seepage_index` gate at `<= 250 m`), high-uncertainty flag ($> 10000\text{ m}$ training gate), separate climatic ($\le 10000\text{ m}$) vs. point edaphic ($\le 2000\text{ m}$) eligibility gates, and inverse coordinate-uncertainty weight. |
| **SSURGO Map Units (`10 to 30 m`)** | `musym`, `muname`, `taxsuborder`, `taxgreatgroup`, `is_rock_outcrop`, `outcrop_class`, `ssurgo_rock_outcrop_pct`, `ssurgo_bedrock_depth_cm`, `ssurgo_sand_pct`, `ssurgo_clay_pct`, `ssurgo_ph`, `ssurgo_om_pct`, `ssurgo_cec`, `ssurgo_aws_0_50cm` | `TEXT` / `INTEGER` / `REAL` | `%`, `cm`, `pH`, $\text{cmol}_+\cdot\text{kg}^{-1}$, `cm` | USDA-NRCS SSURGO map unit taxonomy, `outcrop_class` (`Granitic flatrock`, `Limestone glade`, `Sandstone/Siliceous`, `Serpentine`, `Diabase/Mafic glade`, `Non-outcrop matrix`), lithic contact depth, and $0\text{ to }15\text{ cm}$ horizon properties. |
| **Retired `POLARIS` Schema Columns** | `polaris_ph_30m`, `polaris_sand_30m`, `polaris_clay_30m`, `polaris_om_30m`, `polaris_bd_30m`, `polaris_ksat_30m` | `REAL` | Legacy `NULL` (`NA_real_`) | Legacy $30\text{ m}$ columns retained solely for SQLite schema compatibility and initialized to `NULL` (`NA_real_`) via `sanitize_disabled_polaris_columns()` after retiring partial `polaris_30m` overrides. |
| **SOLUS100 (`100 m`)** | `solus_resdept_cm`, `solus_fragvol_pct`, `solus_sand_pct`, `solus_clay_pct`, `solus_silt_pct`, `solus_cec` | `REAL` | `cm`, `%`, $\text{cmol}_+\cdot\text{kg}^{-1}$ | USDA-NRCS $100\text{ m}$ depth to restriction/lithic contact, coarse fragments, particle-size fractions, and CEC over $0\text{ to }15\text{ cm}$. |
| **Harmonized Rhizosphere (`0 to 15 cm`)** | `fine_soil_ph`, `fine_soil_sand`, `fine_soil_clay`, `soil_silt_pct`, `ilr_texture_1`, `ilr_texture_2`, `fine_soil_cec`, `fine_cfvo_pct`, `fine_bedrock_depth_cm` | `REAL` | `pH`, `%`, ILR, $\text{cmol}_+\cdot\text{kg}^{-1}$, `cm` | Priority-harmonized (`SOLUS100 -> measured SSURGO -> SoilGrids -> NULL`) $0\text{ to }15\text{ cm}$ rhizosphere predictors and two orthogonal Isometric Log-Ratio (`ILR`) texture coordinates freeing sand/silt/clay from the $100\%$ unit-sum constraint. |
| **Uncertainty-Buffered Extremes** | `ph_q10_buffer`, `ph_q90_buffer`, `min_bedrock_depth_buffer_cm`, `max_cfvo_buffer_pct`, `azonal_lithology_presence_buffer` | `REAL` / `INTEGER` | `pH`, `cm`, `%`, `{0, 1}` | Sub-grid micro-habitat quantiles (`[100, 1000] m` radius for `ph_q10_buffer`, `ph_q90_buffer`, and `max_cfvo_buffer_pct`, whereas `min_bedrock_depth_buffer_cm` caps radius at `[100, 300] m` via `clamp_bedrock_buffer_radius_m()` and aggregates finite raster cells via `min()` capped to `[0, 200]` cm). |
| **3DEP Micro-Topography (`30 m`)** | `elevation_m`, `topo_slope_deg`, `topo_tpi_30m`, `topo_hli_30m`, `topo_twi_30m` | `REAL` | `m`, `deg`, `m`, $[0, 1]$, $\ln(a/\tan\beta)$ | Elevation and slope retained across valid coordinates, paired with `coordinateUncertaintyInMeters <= 250 m`-gated Topographic Position Index (`TPI`), McCune-Keon Heat Load Index (`HLI`), and D8 Topographic Wetness Index (`TWI`). |
| **USGS SGMC Bedrock** | `sgmc_lithology_class`, `is_azonal_lithology` | `TEXT` / `INTEGER` | 6-class enum / `{0, 1}` | Standardized bedrock class (`Ultramafic_Serpentine`, `Mafic_Igneous_Metamorphic`, `Carbonate_Calcareous`, `Felsic_Granitic_Gneiss`, `Siliciclastic_Consolidated`, `Unconsolidated_Sedimentary`) and azonal indicator (`1` for ultramafic, mafic, or carbonate). |
| **WorldClim v2.1 (`1 km`)** | `bio1_temp_mean` through `bio19_precip_coldest_qtr` | `REAL` | $^\circ\text{C}$, `%`, `mm` | All 19 standard bioclimatic temperature and precipitation rasters (`bio1` to `bio19`) gated at `<= 10000 m`. |
| **USGS `DS-801` Soil Geochemistry** | `soil_ca_mg_molar_ratio`, `soil_ultramafic_ni_cr_index`, `soil_carbonate_c_pct` | `REAL` | Molar ratio / Index / `%` | Topsoil A-horizon IDW geochemistry capped strictly at `10.0 km` (`max_dist_km = 10.0`) for Step 07 hydro-edaphic auditing (audit-only alongside `ultramafic_calcareous_gradient`, never entering `X_env` or `PCA-env`). |
| **USFWS NWI Wetlands** | `is_nwi_wetland`, `nwi_wetland_system`, `nwi_wetland_dist_m` | `INTEGER` / `TEXT` / `REAL` | `{0, 1}` / System / `m` | National Wetlands Inventory intersection indicator, Cowardin wetland system class (`Palustrine`, `Riverine`, `Upland`), and Euclidean distance (`m`) to nearest mapped wetland polygon. |
| **Coarse Fallbacks** | `soil_ph`, `soil_cec`, `soil_sand`, `soil_clay`, `soil_silt`, `soil_cfvo`, `soil_bulk_density`, `soil_ph_soilgrids` through `soil_bulk_density_soilgrids` | `REAL` | Standard pedological units | Coarse SoilGrids and legacy compatibility columns retained for non-CONUS occurrences outside SSURGO and SOLUS100 coverage. |

### 2.8 `genomic_samples` (Tier A Staged Sequencing Metadata & Polyploid QC Registry)

**Primary Key:** `genetic_sample_id` | **Populated By:** `scripts/data_prep/06_ingest_genetic_data.py` (`preserve_existing_on_null=True` default, backed by `step06_synonym_utils.py`, `step06_sample_metadata_utils.py`, `step06_voucher_linker_utils.py`, `step06_vcf_dosage_utils.py`, and `step06_ploidy_hdplot_utils.py`) alongside `scripts/data_prep/skim_offtarget_organellar.py` and `scripts/data_prep/download_public_packera_sequences.py`
**Companion Exports:** `runs/<RUN_ID>/outputs/genomics_parsed/parsed_genomic_samples.csv` and `allele_dosage_matrix.tsv`

| Column Group | Column Names | SQLite Type | Units / Domain | Description |
| :--- | :--- | :---: | :--- | :--- |
| **Identifiers & Linkage** | `genetic_sample_id`, `catalogNumber`, `project_id`, `population_id`, `replicate_pair_id`, `linkage_tier` | `TEXT` | Sample ID / Tier enum | Cleaned sequencing library ID, linked specimen barcode (`NULL` when unvouchered), project/collection label, locality population token, basal replicate pair key (`M`/`P` maternal vs. population sibling), and `G1_Direct_Voucher` through `G4_Ungeoreferenced_Tip` linkage tier (`HybSeq_TargetCapture` promotes vouchered matches to `G1_Direct_Voucher`, single-locus `Sanger_Barcode` matches stay quarantined under `G3_Public_Sanger_NonMutating` with `NA` predicted taxon and probability, and `G3` downgrades to `G4` when `(0, 0)`, US/CA sign-flipped, or herbarium-rooftop coordinates are nullified). |
| **Taxonomy & Collection Context** | `taxon_raw`, `species_standardized`, `is_putative_novel_or_hybrid_query`, `ploidy_level`, `marker_platform`, `decimalLatitude`, `decimalLongitude`, `stateProvince`, `county`, `eventDate` | `TEXT` / `INTEGER` / `REAL` | Taxon / Ploidy / `WGS84` | Verbatim field code or metadata taxon, resolved *Packera* binomial, hybrid or novel query indicator (`{0, 1}`), user or metadata ploidy integer, assay platform (`nextRAD` / `ddRAD` / `HiFi` / `HybSeq_TargetCapture` / `Sanger_Barcode`), and sanitized field locality coordinates. |
| **Read QC & Paralog Filter** | `hq_reads`, `mapped_reads`, `pct_repetitive`, `n_called_loci`, `missing_genotype_rate`, `observed_het`, `geno_quality_pass`, `hdplot_retained_loci` | `INTEGER` / `REAL` | Counts / $[0, 1]$ / `{0, 1}` | Total high-quality reads, reference-mapped reads, repetitive read fraction, callable SNP count, missingness rate, individual heterozygosity, binary QC pass indicator, and singleton SNP loci retained under McKinney `HDplot` ($H \le 0.60$, $|D| \le 4.0$, using aligned companion VCF depths when dosage tables are staged) without MAF pruning. |
| **Cytotype & Organellar Skim** | `inferred_cytotype`, `nquire_delta_logl`, `plastid_clade`, `its_ribotype`, `cytonuclear_discordant` | `TEXT` / `REAL` / `INTEGER` | `'2x'|'3x'|'4x'|'6x'|'8x'|'mixed/uncertain'` / Signed $\Delta\text{BIC}$ / `{0, 1}` | Weiß et al. (2018) EM Gaussian plus uniform-noise allele-balance cytotype (`M_2x`..`M_6x`, `8x` priors), signed `BIC` support relative to `2x` (`<= -4.0` for empirical `2x`, `0.0` for prior-rescued calls, `>= +4.0` for empirical polyploids), skimmed plastid cpDNA clade, nrITS ribotype label, and empirical-nuclear reconciled cytonuclear discordance indicator. |
| **Public Taxon Placement** | `public_plastome_match`, `public_its_match`, `waxy_lineage_match` | `TEXT` | Taxon match string | Best 11-mer containment and Jaccard taxon match scored against `data/genomics/public_organellar_and_loci/` complete plastomes and nrITS/ETS alignments alongside nuclear `waxy` references (`skimmed_waxy_consensus.fasta`). |

### 2.9 `cytological_records` (Tier A Public & Specimen-Linked Chromosome/C-Value Registry)

**Primary Key:** `cytology_id` | **Populated By:** `scripts/data_prep/download_public_cytology.py`
**Companion Exports:** `data/cytology/public_packera_cytology_records.csv` and `data/cytology/public_packera_cytology_manifest.json`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `cytology_id` | `TEXT` | Primary Key | Deterministic hash identifier built from source database and record ID with standardized taxon and raw count. |
| `source_db`, `source_record_id` | `TEXT` | Registry / ID | Source archive (`CCDB`, `NBDC_ICNA`, `KEW_CVALUES`, `LOCAL_CSV`, `VOUCHER_LABEL`) alongside upstream accession ID. |
| `original_taxon_name`, `basionym_name`, `species_standardized` | `TEXT` | Botanical name | Verbatim source name paired with historical *Senecio* or *Cineraria* basionym when applicable, resolved into an accepted North American *Packera* binomial. |
| `catalogNumber`, `exsiccatae_key`, `cytology_linkage_tier` | `TEXT` | Voucher / Tier | Direct barcode match paired with normalized collector-number exsiccatae key and linkage tier (`C1_Direct_Cytovoucher` through `C4_Taxon_MultiCytotype`). |
| `raw_count_str`, `gametophytic_n`, `sporophytic_2n` | `TEXT` / `INTEGER` | Count string / Chromosomes | Original count string (`22II`, `2n = 44`, `ca. 92`) paired with meiotic haploid count $n$ and resolved mitotic count $2n$. |
| `base_x`, `ploidy_level`, `cytotype_label` | `INTEGER` / `TEXT` | Base $x$ / Ploidy | Inferred paleopolyploid base number ($x \in \{22, 23\}$ or local series) paired with somatic ploidy multiplier ($2..16$) and canonical label (`2x`, `3x`, `4x`, `6x`, `8x`). |
| `b_chromosomes`, `is_aneuploid_or_odd_ploid`, `is_approximate` | `INTEGER` | Count / `{0, 1}` | Supernumerary B-chromosome tally paired with odd-ploid/aneuploid indicator plus `ca.` or range approximation flag. |
| `c_value_1c_pg`, `c_value_2c_pg`, `breeding_system` | `REAL` / `TEXT` | Picograms (`pg`) / Enum | Holoploid genome size in `1C pg` and `2C pg` alongside reported apomixis or sexual breeding system tokens. |
| `genbank_accession`, `collector_surname`, `record_number_clean` | `TEXT` | Accession / Exsiccatae | Linked GenBank sequence accession paired with uppercase collector surname and stripped collection number used for `C2` matching. |
| `institution_code`, `voucher_citation`, `locality_text`, `reference_citation` | `TEXT` | Herbarium / Provenance | Holding herbarium code alongside verbatim cytovoucher string with collection locality and primary literature citation. |

### 2.10 `taxon_cytology_summary` (Tier A Taxon Modal Ploidy & Cytotype Series Registry)

**Primary Key:** `species_standardized` | **Populated By:** `scripts/data_prep/download_public_cytology.py`
**Companion Export:** `data/cytology/public_packera_taxon_cytology_summary.csv`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `species_standardized` | `TEXT` | Primary Key | Accepted *Packera* binomial aggregating all curated chromosome and genome-size records. |
| `n_cytology_records`, `n_direct_cytovouchers` | `INTEGER` | Counts $\ge 0$ | Total curated cytological records paired with directly vouchered `C1`/`C2` counts for the taxon. |
| `modal_sporophytic_2n`, `modal_base_x`, `modal_ploidy_level` | `INTEGER` | Modal counts | Weighted modal $2n$ chromosome number paired with modal base $x$ and modal ploidy level (`2..16`) used as the Step 06 Empirical Bayes monomorphic prior. |
| `reported_cytotypes`, `is_multi_cytotype_taxon` | `TEXT` / `INTEGER` | Set / `{0, 1}` | Sorted comma-separated cytotype series (such as `2x, 4x`) and flag equal to `1` when two or more ploidy levels or dysploid base numbers occur. |
| `has_dysploid_series`, `has_odd_polyploid_or_apomixis` | `INTEGER` | `{0, 1}` | Flags set to `1` when multiple base $x$ numbers (`x = 22` vs `x = 23`) or triploid/odd-ploid/apomictic records appear within the taxon. |
| `mean_2c_pg`, `breeding_system_summary` | `REAL` / `TEXT` | `pg` / Summary | Arithmetic mean `2C` genome size in picograms and concatenated breeding-system annotations. |

### 2.11 `taxon_backbone_profiles` (Tier A Taxonomic Backbone, Protologue, and NWPL Wetland Registry)

**Primary Key:** `species_standardized` | **Populated By:** `scripts/data_prep/download_taxonomic_backbone.py` and `scripts/data_prep/download_geochem_wetland_layers.py`
**Companion Exports:** `data/taxonomy/public_taxon_backbone_profiles.csv` and `data/taxonomy/public_taxonomic_backbone_manifest.json`

| Column Group | Column Names | SQLite Type | Units / Domain | Description |
| :--- | :--- | :---: | :--- | :--- |
| **Identity & Nomenclature** | `species_standardized`, `target_genus`, `accepted_scientific_name`, `accepted_authorship` | `TEXT` | Botanical binomial | Standardized binomial primary key paired with configured `TaxaConfig.target_genus`, accepted CoL WCVP `2232` name, and full authorship string. |
| **Basionym & Synonymy** | `basionym_name`, `basionym_authorship`, `homotypic_synonyms_json`, `heterotypic_synonyms_json` | `TEXT` | Binomial / JSON array | Protologue basionym plus authorship paired with JSON lists partitioning homotypic (nomenclatural) versus heterotypic (taxonomic) synonyms. |
| **Hybrid Formula** | `is_hybrid_taxon`, `hybrid_formula` | `INTEGER` / `TEXT` | `{0, 1}` / Formula | Named nothotaxon flag alongside parsed parent-cross hybrid formula when documented in CoL WCVP `2232`. |
| **Protologue & Type** | `publication_citation`, `publication_year`, `protologue_url`, `locus_classicus_text` | `TEXT` / `INTEGER` | Citation / Year / URL | Protologue bibliography, priority publication year, Biodiversity Heritage Library API v3 deep-link URL, and transcribed *locus classicus* string. |
| **Distribution & Status** | `tdwg_native_codes_json`, `tdwg_introduced_codes_json`, `natureserve_element_id`, `natureserve_global_rank` | `TEXT` | JSON / Element / `G`-rank | TDWG Level-3 native alongside introduced country/state codes plus NatureServe Explorer global conservation rank (`G1` through `G5`). |
| **Wetland Indicator & Audit** | `nwpl_indicator_code`, `nwpl_wetland_score`, `backbone_source`, `updated_at` | `TEXT` / `REAL` | NWPL code / `1..5` / ISO | USACE National Wetland Plant List indicator (`OBL`, `FACW`, `FAC`, `FACU`, `UPL`) paired with numeric wetland fidelity score, provenance string, and UTC sync timestamp. |

### 2.12 `natureserve_subnational_ranges` (Tier A NatureServe State & Province Range Table)

**Primary Key:** `(species_standardized, country_code, subnation_code)` | **Populated By:** `scripts/data_prep/download_taxonomic_backbone.py`
**Companion Export:** `data/taxonomy/public_natureserve_subnational_ranges.csv`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `species_standardized` | `TEXT` | Composite PK | Accepted standardized species binomial matching `taxon_backbone_profiles`. |
| `country_code`, `subnation_code` | `TEXT` | ISO-2 / State-Prov | Two-letter country (`US` / `CA`) paired with state or Canadian province postal code (`NC`, `VA`, `ON`, `QC`). |
| `s_rank` | `TEXT` | NatureServe `S`-rank | Subnational conservation rank (`S1` through `S5`, `SNA`, `SH`, `SX`, or `SNR`). |
| `is_native`, `is_exotic`, `is_historical_or_extirpated` | `INTEGER` | `{0, 1}` | Subnational native residency flag paired with introduced/waif flag (`natureserve_exotic_or_waif`) and historical/extirpated flag queried by `filter_vouchers_by_usda_county_range()`. |

## 3. Tier B Run-Scoped SQLite Tables (`runs/<RUN_ID>/` & SQLite)

### 3.1 `pipeline_runs` and `run_active_cohort` (Preflight Cohort Isolation)

* **`pipeline_runs` (PK `run_id`):** Tracks the 8 run provenance columns (`run_id`, `run_dir`, `start_timestamp`, `end_timestamp`, `git_commit_sha`, `cli_arguments_json`, `focal_taxon`, and `intraspecific_mode`).
* **`run_active_cohort` (PK `(run_id, catalogNumber, plant_individual_id)`):** Snapshots every individual included in `<RUN_ID>` along with its preflight-refreshed `species_standardized`, `determiner_tier`, `type_status`, `is_type_specimen`, `is_promoted_bronze`, `authority_weight` (`3.0` Type, `2.0` Gold, `0.75` Silver, `0.50` `Promoted_Bronze`, `0.0` unverified or unpromoted Bronze), and `spatial_balance_weight`.

### 3.2 `leaf_efa_harmonics` (Step 03 Individual-Leaf EFA & ECT Table)

**Primary Key:** `(run_id, catalogNumber, plant_individual_id, leaf_id)` | **Populated By:** `scripts/morphometrics/03_fourier_extractor.R`
**Companion CSV:** `runs/<RUN_ID>/data/tables/leaf_efa_harmonics_individual.csv` (which also exports CSV-joined `shape_id`, `species_standardized`, and `determiner_tier` alongside SQLite DDL columns)

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `run_id`, `catalogNumber`, `plant_individual_id`, `leaf_id` | `TEXT` / `INTEGER` | Composite PK | Run, specimen, individual, and leaf keys (`shape_id` joins via companion CSV). |
| `leaf_zone`, `leaf_zone_confidence`, `dist_from_crown_norm`, `cy_rel_plant_bbox`, `mask_is_upright` | `TEXT` / `REAL` / `INTEGER` | Zone / $[0, 1]$ / `{0, 1}` | Shoot zone (`{'basal', 'cauline'}`), continuous Step 02 zone probability, normalized distance from caudex crown, relative plant position, and apex-up orientation flag. |
| `type_status`, `is_type_specimen`, `is_promoted_bronze`, `anchor_tier_source`, `authority_weight`, `is_passive`, `spatial_balance_weight` | `TEXT` / `INTEGER` / `REAL` | Authority metadata | Taxon authority provenance (`Tier_0_Type_Anchor`, `Tier_1_Gold`, `Tier_2_Silver`, `Promoted_Bronze`, `Tier_3_Bronze`), active weight, passive projection flag, and spatial balance weight synchronized from `run_active_cohort`. |
| `solidity`, `dissection_index`, `aspect_ratio`, `area_px`, `area_mm2`, `is_dissected`, `blade_cie_l_norm`, `blade_cie_l_sd`, `blade_cie_l_sheet_delta`, `blade_cie_a_delta`, `blade_cie_b_delta`, `pixels_per_mm`, `is_calibrated_mm`, `scale_source`, `petiole_length_mm`, `petiole_width_mm`, `blade_length_mm`, `blade_width_mm`, `petiole_to_blade_ratio` | `REAL` / `INTEGER` / `TEXT` | Morphometrics & `mm` | Leaf-level shape ratios, calibrated linear dimensions (`NULL` when `is_calibrated_mm = 0`), eroded-lamina CIE $L^*$ reflectance (`blade_cie_l_norm`), within-blade SD (`blade_cie_l_sd`), sheet-relative color deltas (`blade_cie_l_sheet_delta`, `blade_cie_a_delta`, `blade_cie_b_delta`), and scale provenance (`scale_source`). |
| `A1..A24`, `B1..B24`, `C1..C24`, `D1..D24` | `REAL` | Dimensionless EFA | 96 Elliptic Fourier coefficients (`Momocs::efourier(..., norm = FALSE, start = FALSE)` up to $H_{\max} = 24$; harmonics above $H_{\text{opt}}$ are `0` or `NULL`). |
| `PC1..PC15` | `REAL` | Dimensionless EFA PCs | Authority-weighted (`stats::cov.wt`) symmetric EFA principal component scores ($A_n, D_n$ lateral sub-matrix). |
| `sym_var`, `asym_var`, `pct_asymmetry` | `REAL` | Power / `%` $[0, 100]$ | Symmetric ($A_n^2 + D_n^2$) vs. asymmetric ($B_n^2 + C_n^2$) harmonic power and bilateral asymmetry percentage (`> 30%` flags pressing/taphonomic distortion). |
| `harmonic_roughness` | `REAL` | Dimensionless $[0, 1]$ | High-frequency harmonic energy fraction ($n = 9..\min(H_{\text{opt}}, 24)$) screening jagged segmentation tears. |
| `is_contour_outlier` | `INTEGER` | `{0, 1}` | Equals `1` when leaf fails Step 03 R-side clump/boundary QC, high-frequency roughness, or within-taxon FastMCD ($p < 0.001$; `0` for protected types). |
| `leaf_quality_weight` | `REAL` | Weight $[0, 1]$ | Continuous Tukey bisquare quality weight used in weighted median OTU aggregation. |
| `centroid_size`, `log_centsize` | `REAL` | `mm` (or `px`) / $\log_{10}$ | Contour centroid size and base-10 log centroid size used in Step 04 `RRPP` allometry. |
| `efd_harmonic_opt`, `efd_k95_pcs` | `INTEGER` | $[6, 24]$ / $[5, 15]$ | Cohort-calibrated cumulative $\ge 99\%$ Fourier harmonic cutoff ($H_{\text{opt}}$) and $\ge 95\%$ PCA subspace rank ($k_{95}$). |
| `ECT_PC1`, `ECT_PC2`, `ECT_PC3` | `REAL` | Dimensionless SVD scores | Top 3 authority-weighted SVD axes of the $32 \times 360$ polar Euler Characteristic Transform capturing pinnatifid and lyrate sinus topology without perimeter parameterization drift. |
| `ect_lobe_complexity` | `REAL` | Dimensionless $\ge 1$ | Mean positive topological component multiplicity across radial ECT slices (`1.0` for entire blades; elevated across deeply pinnatisect *Packera* taxa). |
| `ect_bilateral_asym` | `REAL` | Dimensionless $\ge 0$ | Polar reflectional L2 discrepancy across the $32 \times 360$ ECT field. |

### 3.3 `otu_morphometrics` (Step 03 OTU Medians & Step 04 GMM, CDA, Allometry & OSR)

**Primary Key:** `(run_id, catalogNumber, plant_individual_id)` | **Populated By:** `03_fourier_extractor.R`, `04_gmm_morphotools.R`, and `step04_allometry_merger_utils.R` (restoring raw Step 03 baselines from SQLite `leaf_efa_harmonics` via `.restore_raw_step03_from_leaf_harmonics_db()` during SQLite-only re-runs)
**Companion CSVs:** `leaf_efa_harmonics.csv` (Step 03 OTU output), `morphometric_flags.csv` (Step 04 enriched OTU output), `gmm_bayes_factors_summary.csv` (exporting `ICL`, `Delta_ICL_vs_Null`, `Boot_LRT_Stat_K1_vs_K2`, `Boot_LRT_P_K1_vs_K2`, `Boot_LRT_Stat_prev`, and `Boot_LRT_P_prev`), and `batch_variance_audit.csv` (marginal drop-1 ANCOVA controlling for `species_standardized` plus non-aliased `regional_group`)

| Column Group | Column Names | SQLite Type | Units / Domain | Description |
| :--- | :--- | :---: | :--- | :--- |
| **OTU Provenance & Authority** | `species_standardized`, `determiner_tier`, `morph_tier_source`, `leaf_count`, `basal_leaf_count`, `cauline_leaf_count`, `is_promoted_bronze`, `is_passive`, `authority_weight`, `spatial_balance_weight`, `scale_source`, `is_calibrated_mm`, `pixels_per_mm` | `TEXT` / `INTEGER` / `REAL` | Zone / Counts | Records whether primary EFA/ECT medians come from high-confidence `Basal_Rosette` (`>= 0.65` confidence), `Transitional_Passive` (`< 0.65` confidence), or `Cauline_Only_Fallback` (`0` basal leaves), alongside active (`is_passive = 0`) vs. passive (`1`) projection status (`population_deme_id` exports in `morphometric_flags.csv`); `is_step04_calibrated_mm()` masks `tick_lattice_unboxed`, `mp_estimated`, `predicted_from_megapixels`, and `uncalibrated` rows out of physical-mm trait merges and allometry models. |
| **Pure Basal Shape, Size & Reflectance** | `A1..A24`, `D1..D24`, `PC1..PC15`, `ECT_PC1..ECT_PC3`, `ect_lobe_complexity`, `solidity`, `q25_solidity`, `dissection_index`, `dissected_leaf_fraction`, `is_dissected`, `aspect_ratio`, `area_px`, `area_mm2`, `centroid_size`, `log_centsize`, `blade_length_mm`, `blade_width_mm`, `petiole_length_mm`, `petiole_width_mm`, `petiole_to_blade_ratio`, `blade_cie_l_norm`, `blade_cie_l_sd`, `blade_cie_l_sheet_delta`, `blade_cie_a_delta`, `blade_cie_b_delta`, `sym_var`, `asym_var`, `pct_asymmetry`, `harmonic_roughness`, `is_contour_outlier`, `leaf_quality_weight`, `foliar_variance`, `within_rosette_efa_dispersion`, `efd_harmonic_opt`, `efd_k95_pcs` | `REAL` / `INTEGER` | Dimensionless & `mm` | Unweighted basal rosette medians (`A_n, D_n` symmetric harmonics only; `B_n, C_n` excluded from OTU aggregation as taphonomic bending terms), lower-quartile solidity `q25_solidity`, perimeter dissection index $P / (2\sqrt{\pi A})$, paper-normalized reflectance `blade_cie_l_norm` alongside browning/yellowing deltas (`blade_cie_l_sheet_delta`, `blade_cie_a_delta`/`blade_cie_a_browning_delta`, `blade_cie_b_delta`/`blade_cie_b_yellowing_delta`, with only `^blade_cie_` columns centered by batch correction and `blade_cie_l_sd` withheld from covariance fitting), `NULL` physical `*_mm` and `centroid_size` outputs when `is_calibrated_mm = 0` (`petiole_to_blade_ratio` retained across all leaves), and within-rosette dispersion (`plant_height_cm` and `scape_to_basal_leaf_ratio` join from Tier A `plant_individuals` in `v_integrative_triage_master`, whereas `circularity`, `extent`, and `ect_bilateral_asym` reside at leaf level). |
| **Cauline & Heteroblasty Track** | `cauline_solidity_med`, `cauline_dissection_index_med`, `cauline_aspect_ratio_med`, `cauline_blade_length_mm`, `cauline_ect_lobe_complexity`, `basal_vs_cauline_length_ratio`, `basal_vs_cauline_solidity_delta`, `heteroblasty_efa_distance`, `heteroblasty_ect_distance`, `Cauline_PC1..Cauline_PC3`, `Cauline_ECT_PC1..Cauline_ECT_PC2` | `REAL` | Dimensionless & `mm` | Orthogonal stem-leaf morphology and within-plant basal-to-cauline developmental divergence (`NULL` unless both basal and cauline leaves occur on the individual). |
| **Reproductive Geometry** | `AR_inv`, `AR_inv_std`, `log1p_capitulum_count`, `log1p_capitulum_count_std`, `reproductive_trait_status` | `REAL` / `TEXT` | Ratio / Z-score | Involucre aspect ratio, standardized scores, log capitulum count computed from `pmax(capitula_count, corymb_head_count_est)` (`NA` for sterile `0`-head sheets), and availability status (`Observed` vs. `Vegetative_Missing`). |
| **Canonical GMM & Shrinkage CDA** | `selected_morph_features`, `gmm_cluster`, `gmm_uncertainty`, `cda_predicted_taxon`, `cda_posterior_prob`, `cda_oof_posterior_prob`, `cl_target_threshold`, `self_morph_mahal_d2`, `self_morph_df`, `cda_second_taxon`, `cda_second_posterior_prob`, `cda_posterior_margin`, `loocv_accuracy`, `can1`, `can2`, `can3`, `triage_priority`, `discordance_reason` | `TEXT` / `INTEGER` / `REAL` | Semicolon list / `TEXT` ID / $[0, 1]$ / $D^2$ | Semicolon-delimited active predictors, canonical GMM cluster (`TEXT` cluster ID `1..G` after sequential covariance-locked bootstrap gating), in-sample and 10 km deme-blocked OOF shrinkage CDA posteriors (computed with fold-internal median imputation, 75%-core Mahalanobis training outlier trimming, and subspace re-selection), Confident Learning self-confidence threshold, LOO within-taxon Mahalanobis distance and `df`, canonical variates `can1..can3`, and preliminary Step 04 triage metadata. |
| **EPA Ecoregion CDA** | `ecoregion_cda_predicted`, `ecoregion_cda_posterior_prob`, `ecoregion_can1`, `ecoregion_can2` | `TEXT` / `REAL` | Ecoregion / $[0, 1]$ / Canonical SD | Supervised biogeographic discriminant axis scores testing regional morphological structuring across US EPA Level III Ecoregions; `is_centroid_or_inst_coord` sheets mask `lat`/`lon` before polygon lookup and project passively without entering active CDA or MANOVA training pools. |
| **Type Control, OSR & Label Noise** | `type_status`, `is_type_specimen`, `has_empirical_type_anchor`, `type_morph_mahal_d2`, `nearest_type_morph_taxon`, `nearest_type_basionym`, `nearest_type_basionym_d2`, `type_morph_concordance`, `is_typical_morphotype`, `min_congener_mahal_d2`, `min_congener_taxon`, `open_set_novelty_p`, `open_set_df`, `morphometric_outlier_flag`, `misidentification_flag`, `is_confident_label_noise` | `INTEGER` / `REAL` / `TEXT` | `{0, 1}`, $D^2$, $\chi^2$ tail $p \in [0, 1]$ | Squared Mahalanobis distance to nomenclatural type morphotype ($c^*_k$) and basionym type centroids, Open-Set Recognition (`OSR`) minimum squared distance and exact feature dimensionality `open_set_df` across congeneric centroids, and Northcutt Confident Learning label-noise flag. |

### 3.4 `spatial_rf_predictions` (Step 05 Independent Niche Regimes, Spatial RF & Hopkins Phenology)

**Primary Key:** `(run_id, catalogNumber, plant_individual_id)` | **Populated By:** `05_multimodal_spatial_rf.R`, `step05_config_geo_utils.R`, `step05_phenology_ecotypes.R`, and `step05_sdm_maxent.R` *(100% independent of Step 04 morphometrics)*
**Companion CSVs:** `runs/<RUN_ID>/data/tables/multimodal_conflict_flags.csv`, `runs/<RUN_ID>/outputs/reports/phenological_anomaly_summary.csv`, `runs/<RUN_ID>/outputs/reports/hopkins_circular_model_selection.csv`, `runs/<RUN_ID>/outputs/reports/multimodal_spatial_rf_summary.csv`, `runs/<RUN_ID>/outputs/reports/micro_edaphic_outcrop_contingency.csv`, `runs/<RUN_ID>/outputs/reports/environmental_outlier_quarantine_audit.csv`, and `runs/<RUN_ID>/outputs/reports/sdm_model_performance.csv`

| Column Group | Column Names | SQLite Type | Units / Domain | Description |
| :--- | :--- | :---: | :--- | :--- |
| **Keys & Independent Regime / Habitat Indices** | `lat_lon_key`, `ecotype_regime`, `lithic_outcrop_index`, `summer_moisture_deficit_idx`, `maritime_continentality_idx` | `TEXT` / `REAL` | Coordinate key / Regime / Indices | Joined spatial cache key (`%.5f_%.5f`), strictly morphometric-independent Step 05 `resolve_analysis_group()` regime (`Ecotype_1_...` through `Ecotype_4_...`, or canonical species name when $K_{\text{eco}} = 1$, excluding quarantined rooftop/centroid and non-type label-noise rows from discovery), `[0, 1]` lithic outcrop confinement integrating point and buffered `min_bedrock_depth_buffer_cm`, warm-season moisture deficit, and coastal vs. interior continentality index wiring `bio2_mean_diurnal_range` directly alongside `bio6` and `bio4`. |
| **Spatial Block-CV RF** | `rf_predicted_taxon`, `rf_predicted_prob`, `rf_given_taxon_prob`, `rf_prediction_source`, `edaphic_best_fit_taxon`, `spatial_block_cv_acc`, `spatial_autocorr_range_km`, `spatial_block_id` | `TEXT` / `REAL` | Taxon / $[0, 1]$ / `km` / Block ID | Out-of-fold Spatial Block-CV Random Forest top taxon (or `ecotype_regime` in $K = 1$ discovery), max probability, probability assigned to recorded label, prediction provenance (`Spatial_Block_OOF_Anchor`, `Fallback_OOB_Anchor`, `Passive_Projection_Bronze`, `Passive_Projection_High_Uncertainty`, or `Passive_Ungeoreferenced`), empirical variogram block range (`km`), and fold block ID (`n_spatial_rf_predictors` exports in companion CSV); `Spatial_Block_CV_Accuracy` scores strictly held-out fold anchors (`compute_fold_reclustered_ecotype_ref()` on post-selection `Ecotype_*` groups) with `X_env` fit independently of `MEM1..MEM5`, while spherical Haversine `km` matrices (`compute_haversine_dist_km_matrix()`) drive `MEM` eigenfunctions, Moran's *I* weights, and `125 km` accessible `M_i` envelopes. |
| **Selected Predictor Manifests & Climate** | `primary_env_predictors`, `primary_clim_predictors`, `primary_edaph_predictors`, `bio1_temp_mean` through `bio19_precip_coldest_qtr` | `TEXT` / `REAL` | Semicolon lists / WorldClim units | Non-collinear (`VIF <= 5.0`) predictor subsets across environmental, macro-climatic, and $0\text{ to }15\text{ cm}$ edaphic/topographic blocks read dynamically by Step 07 Complex Discovery (excluding `ultramafic_calcareous_gradient` and `polaris_*_30m`), plus cached WorldClim `bio1..bio19` values (`bio2_mean_diurnal_range` wired into `maritime_continentality_idx`). |
| **Circular Hopkins Phenology** | `doy`, `expected_doy`, `hopkins_expected`, `delta_doy`, `hopkins_adjusted_doy`, `pheno_timing_category`, `anthesis_filter_status` | `REAL` / `TEXT` | Days (`d`) / Category | Observed DOY hydrated from `voucher_phenological_states.csv` or SQLite `voucher_reproductive_phenology` (`fetch_step02_phenology_from_sqlite()` integrating `corymb_head_count_est` and `anthesis_head_count`, paired with rooftop baseline quarantine and `load_step05_macroecology_config()` at `min_flowering_doy = 30`), empirical Hopkins baseline DOY, signed isolation residual $\Delta\text{DOY}$ (`[-182.625, +182.625]` days via $\text{atan2}$ circular wrap), Hopkins-adjusted flowering day, and optical anthesis filter tag. |
| **Locus Classicus & Substrate** | `edaphic_z_score`, `climatic_z_score`, `type_edaphic_z_score`, `nearest_type_niche_taxon`, `matches_type_lithology`, `substrate_specialist_flag` | `REAL` / `TEXT` / `INTEGER` | Robust RMS $Z$ / `{0, 1}` | Multivariate rhizosphere and climatic deviations from taxon reference anchors and *locus classicus* type localities, bedrock match indicator, and outcrop/azonal specialist flag (`edaphic_schoener_d` and Warren/contingency summaries export to companion CSVs with post-selection `p`-values suppressed to `NA` on `Ecotype_*` partitions). |
| **Multimodal Concordance** | `n_observed_modalities`, `multimodal_concordance`, `multimodal_conflict_flag`, `triage_category`, `conflict_discordance_notes` | `INTEGER` / `TEXT` / `INTEGER` | `0..3` / Concordance enum / `{0, 1}` | Count of observed Step 05 streams (`environmental`, `phenology`), concordance status string, binary ecological/phenological conflict flag, preliminary Step 05 triage label, and explanatory conflict notes. |

*(Fine-grain `environmental_cache` SSURGO, SOLUS100, 3DEP topography, and USGS SGMC lithology columns join onto Step 05 metrics through `lat_lon_key` inside `v_integrative_triage_master`, while legacy `polaris_*_30m` columns stay `NULL` and `load_and_update_fine_grain_cache()` skips redundant raster queries on `was_in_cache` hits.)*

#### Companion Step 05 Macroecology & Phenology Report CSV Schemas

| Companion CSV | Producer Module | Key Columns & Statistical Role |
| :--- | :--- | :--- |
| **`phenological_anomaly_summary.csv`** | `step05_phenology_ecotypes.R` & `step05_config_geo_utils.R` | Exports `Analysis_Section` rows (`Baseline_Regression`, `Circular_Phase_Test` via weighted Watson-Williams `F` and permutation `p`, `Taxon_Distribution`, and BH-adjusted pairwise permutation `Circular_Pairwise_PostHoc` via `summarize_weighted_circular_phenology()`) paired with `Selected_Hopkins_Model`, `Hopkins_AICc_Weight`, `Rayleigh_R_Bar`, `Rayleigh_P`, `Circ_Watson_F`, and `Circ_Perm_P`. |
| **`hopkins_circular_model_selection.csv`** | `step05_circular_hopkins_models.R` | Compares circular von Mises baseline models `M0_Constant_Phase`, `M1_Hopkins_1918_Fixed`, `M2_Empirical_Lat`, and `M3_Empirical_TopoClimate_Year` by `AICc` and Akaike weight. |
| **`multimodal_spatial_rf_summary.csv`** | `05_multimodal_spatial_rf.R` | Stores partitioned `PCA-env` Warren niche identity/background metrics (`D_clim`, `D_edaph`, `D_joint`, `Schoener_D`, `p_value_D`, `p_value_I`, `p_bg_D_*`, and `p_equiv_D_*`) paired with `CV_Mean_Fold_Predictors` and `CV_Fold_Predictor_Jaccard`; all `p`-values suppress to `NA` on post-selection `Ecotype_*` partitions. |
| **`micro_edaphic_outcrop_contingency.csv`** | `edaphic_fine_grain_utils.R` | Reports substrate/outcrop odds ratios alongside BH-adjusted probabilities (`P_Value`, `P_Value_BH`, `Outcrop_Fisher_P`, and `Top_Substrate_P_Value_BH`), and suppresses `SSURGO_Outcrop` alongside binary `Outcrop_Fisher_P` `p`-values to `NA` on post-selection `Ecotype_*` rows. |
| **`sdm_model_performance.csv`** | `step05_sdm_maxent.R` | Records `Selected_Feature_Class`, `Selected_Reg_Multiplier`, `MaxEnt_AICc`, and `MaxEnt_NonZero_Betas` alongside held-out-only `AUC_Block_*` and `CBI_*` (`ok_oof_mx`, `ok_oof_rf`, and `ok_oof_ens`) across `MaxEnt`, `SpatialRF`, and `Ensemble`. |

### 3.5 `genetic_lineage_assignments` (Step 06 Independent Phylogenomics & Coalescent Delimitation)

**Primary Key:** `(run_id, genetic_sample_id)` | **Populated By:** `scripts/analysis/06_genomic_phylogenomics.R`, `scripts/analysis/step06_dynamic_tier_engine.R`, `scripts/analysis/step06_phylogenomics_utils.R`, `scripts/analysis/step06_model_selection.R`, `scripts/analysis/step06_public_backbone.R`, `scripts/analysis/step06_reticulation_coancestry.R`, `scripts/analysis/step06_fbranch_quartets.R`, `scripts/analysis/step06_coalescent_models.py`, `scripts/analysis/step06_coalescent_alignment_utils.py`, `scripts/analysis/step06_coalescent_delimitation.py`, and `scripts/analysis/step06_spatial_landscape_genomics.R` *(100% independent of Step 04 morphometrics and Step 05 ecological niche modeling prior to post-hoc MMRR)*
**Companion CSVs & Manifest:** `runs/<RUN_ID>/outputs/reports/genomic_tier_manifest.json` (mirrored to `outputs/genomics/`), `runs/<RUN_ID>/outputs/genomics/genomic_assignments.csv`, `admixture_k_model_selection.csv`, `genomic_model_metrics.csv`, `genomic_clade_cohesion_summary.csv`, `coancestry_matrix_finerad.csv`, `hybrid_index_q12_triangle.csv`, `introgressive_fbranch_summary.csv`, `quartet_concordance_factors_snps2cf.csv`, `svdquartets_support.csv`, `substitution_model_selection.csv`, `coalescent_model_adequacy.csv`, `coalescent_gdi_delimitation.csv`, `mmrr_model_selection.csv`, `landscape_mmrr_feems_summary.csv`, `public_backbone_concordance.csv`, and `public_organellar_placement.csv`

| Column Group | Column Names | SQLite Type | Units / Domain | Description |
| :--- | :--- | :---: | :--- | :--- |
| **Keys & DAPC Classification** | `run_id`, `genetic_sample_id`, `catalogNumber`, `population_id`, `geno_predicted_taxon`, `geno_posterior_prob`, `geno_margin`, `genetic_cluster` | `TEXT` / `REAL` | Keys / Taxon / $[0, 1]$ | Composite key, optional voucher foreign link, population token, Leave-One-Replicate-Out DAPC top taxon calibrated via fold-internal analytical Ledoit-Wolf shrinkage (`[0.01, 0.95]`, logged in `genomic_model_metrics.csv`), posterior probability, top-two margin, and unsupervised PAM silhouette cluster (`Geno_Cluster_k`). |
| **Admixture & Phylogeny** | `q_max`, `q_admixture_json`, `is_geno_admixed`, `geno_parent_1`, `geno_parent_2`, `nearest_tree_sister_taxon`, `patristic_dist_to_dubia`, `replicate_pair_concordance` | `REAL` / `TEXT` / `INTEGER` | Ancestry / Path / `{0, 1}` | Dominant `sNMF` simplex ancestry coefficient under 5-fold 1-SE masked cross-entropy `K` selection (`admixture_k_model_selection.csv`), JSON ancestry dictionary, admixture indicator, parent candidates, nearest Newick neighbor taxon, patristic distance to focal *Packera* clade, and `M`/`P` replicate allele-sharing concordance. |
| **Hybrid Triangle & Coancestry** | `hybrid_index_h`, `interspecific_ancestry_q12`, `hybrid_generation_class`, `coancestry_cluster` | `REAL` / `TEXT` | $[0, 1]$ / Class enum | Pairwise $(h, Q_{12})$ genome-wide ancestry and interspecific heterozygosity scored across 7 multinomial models (`Pure_Parent_1`, `Pure_Parent_2`, `F1_Hybrid`, `F2_Hybrid`, `BC1_Parent_1`, `BC1_Parent_2`, `Advanced_Introgressant`, or `Divergent_NonHybrid_Lineage`) by `AICc` and BH-LRT (`best_model_aicc`, `akaike_weight`, `delta_aicc_vs_pure`, `lrt_q_bh` in `hybrid_index_q12_triangle.csv`), alongside `fineRADstructure` shared-haplotype Ward community. |
| **Reticulation & Coalescent `gdi`** | `f_branch_max`, `f_branch_donor_taxon`, `snaq_gamma_minor`, `coalescent_gdi`, `bpp_introgression_bf` | `REAL` / `TEXT` | Ratio $[0, 1]$ / BF | Maximum `Dsuite`-$f_b$ branch excess allele sharing (`polarization_source` in `introgressive_fbranch_summary.csv`), leading donor taxon, 16-quartet `SNPs2CF` minor inheritance $\gamma$ gated on `delta_aicc_reticulation >= 2.0` and `g_test_p < 0.05`, locus-bootstrapped invariant-site `122 bp` model-corrected `gdi` (`substitution_model` plus `bf_method` in `coalescent_gdi_delimitation.csv`), and coalescent introgression Bayes factor (`bpp_stepping_stone` or `bic_coalescent_mixture`). |
| **Public Nuclear Backbone** | `public_backbone_clade`, `public_phytop_ils_index`, `public_phytop_ih_index` | `TEXT` / `REAL` | Clade / $[0, 1]$ | Nuclear species-tree clade assignment across `data/genomics/public_nuclear_backbone/public_packera_*_Astral_lpp.tre` backbones (Moore-Pollard et al. 2026) paired with `PhyTop` incomplete lineage sorting (`public_phytop_ils_index`) and introgression (`public_phytop_ih_index`) conflict indices (`public_backbone_concordance.csv` and `public_organellar_placement.csv` alongside `genomic_clade_cohesion_summary.csv`). |
| **Post-Hoc Landscape `MMRR`** | `mmrr_beta_ibd`, `mmrr_beta_ibe`, `mmrr_p_ibe` | `REAL` | Partial $\beta$ / Permutation $p$ | Read-only post-hoc Wang (2013) `MMRR` standardized partial coefficients for log-geographic IBD vs. environmental IBE alongside 9,999-matrix Mantel significance, paired with nested `MMRR` comparisons (`mmrr_model_selection.csv`) plus leave-one-node-out `feems_cv_lambda` and `small_n_guardrail_active` inside `landscape_mmrr_feems_summary.csv`. |

#### Companion Step 06 Model Selection & Hypothesis CSV Schemas

| Companion CSV | Producer Module | Key Columns & Statistical Role |
| :--- | :--- | :--- |
| **`admixture_k_model_selection.csv`** | `step06_model_selection.R` | Evaluates `K = 1..k_max` (including $K = 1$ panmixia) via 5-fold masked cross-entropy and selects `best_k` under the 1-SE parsimony rule. |
| **`genomic_model_metrics.csv`** | `06_genomic_phylogenomics.R` / `step06_dynamic_tier_engine.R` | Stores the 13 baseline NextRAD metrics (`n_samples_total` through `mean_replicate_concordance`, including `selected_shrinkage_lambda`, `loro_cv_logloss`, and `selected_k_cross_entropy`) plus the 6 dynamic tier columns (`active_genomic_tier`, `focal_taxon_scoped`, `n_focal_samples`, `barcode_gap_ratio`, `focal_monophyly_status`, `sampling_gap_flag`). |
| **`genomic_clade_cohesion_summary.csv`** | `step06_model_selection.R` | Reports ML and ME substitution adequacy plus clade bootstrap support relative to the public ASTRAL backbone (`Robinson-Foulds` distance). |
| **`hybrid_index_q12_triangle.csv`** | `step06_reticulation_coancestry.R` | Extends `hybrid_index_h` and `interspecific_ancestry_q12` with 7-model multinomial `best_model_aicc`, `akaike_weight`, `delta_aicc_vs_pure`, and `lrt_q_bh`. |
| **`quartet_concordance_factors_snps2cf.csv`** | `step06_reticulation_coancestry.R` | Reports 16-quartet mean concordance (`cf_12_34`, `cf_13_24`, `cf_14_23`) paired with `delta_aicc_reticulation` and `g_test_p` plus hypothesis-gated `snaq_gamma_minor`. |
| **`substitution_model_selection.csv`** | `step06_coalescent_models.py` | Ranks `JC69`, `K80`, `F81`, `HKY85`, and `TN93` with `+G4` and conditional Lewis (2001) `+ASC` (`prop_invariant < 0.05` and `n_sites >= 30`); columns store `best_model`, `base_model`, `has_gamma`, `asc_applied`, `bic_weight`, `ti_tv_kappa`, and `gamma_alpha`. |
| **`coalescent_model_adequacy.csv`** | `step06_coalescent_models.py` | Compares observed block divergence against expected model-corrected block divergence (`compute_model_corrected_divergence_blocks()`). |
| **`svdquartets_support.csv`** | `step06_coalescent_delimitation.py` | Adds 200-replicate quartet block-bootstrap `quartet_boot_support` and `quartet_p_value` onto quartet SVD splits. |
| **`coalescent_gdi_delimitation.csv`** | `step06_coalescent_delimitation.py` | Stores selected `substitution_model` (`HKY` or `JC69` synced to `bpp_msci.ctl`), `bf_method` (`bpp_stepping_stone` or `bic_coalescent_mixture`), `bf_10`, `bic_m0`, `bic_m1`, `w_introgr`, `tau_introgr`, `tau_sp`, $\widehat{gdi}$, and `gdi_source` (`nextrad_msc`, `hybseq_supermatrix_msc`, or `sanger_multilocus_approx`). |
| **`mmrr_model_selection.csv`** | `step06_spatial_landscape_genomics.R` | Compares nested `MMRR` models (`IBD_only` vs. `IBD_plus_IBE`) alongside leave-one-node-out `feems_cv_lambda` in `landscape_mmrr_feems_summary.csv`. |

### 3.6 `curation_triage` (Step 07 Five-Stage Integrative Synthesis)

**Primary Key:** `(run_id, catalogNumber, plant_individual_id)` | **Populated By:** `scripts/analysis/07_triage_dashboard_synthesis.R` (`step07_ancova_merge_utils.R` paired with `step07_hopkins_consensus_utils.R`, `step07_biogeo_hybrid_utils.R` alongside `step07_cluster_cohesion_utils.R`, `step07_decision_stages_utils.R` paired with `step07_postprocess_guards_utils.R`, `step07_cytology_integration.R`, `step07_genomic_propagation.R`, `step07_deme_propagation_utils.R`, `step07_organellar_triage_utils.R`), `genomic_triage_utils.R`, plus `edaphic_geochem_wetland_utils.R`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `run_id`, `catalogNumber`, `plant_individual_id` | `TEXT` / `INTEGER` | Composite PK | Run, specimen, and DBSCAN individual identifiers. |
| `w_morph_kappa`, `w_edaphic_kappa`, `w_pheno_kappa`, `w_geno_kappa` | `REAL` | Reliability $\in [0.05, 0.95]$ | Empirical out-of-fold Cohen's $\kappa$ reliability weights computed on active `Tier_0_Type_Anchor`, `Tier_1_Gold`, and `Tier_2_Silver` anchors (plus `G1`/`G2` sequenced samples for `w_geno_kappa`). |
| `n_available_modalities` | `INTEGER` | `1..4` | Number of non-missing evidence streams (Morphology, Rhizosphere/Climate, Verified-Anthesis Phenology, and Step 06 Genomics) observed for this individual; weights renormalize strictly over observed streams. |
| `probabilistic_consensus_conflict`, `probabilistic_consensus_score` | `REAL` | Probability $[0, 1]$ | Precision-modulated cross-modal conflict probability $P_{\text{conflict},i}$ and concordance score $1 - P_{\text{conflict},i}$. |
| `hopkins_adjusted_doy`, `delta_doy_resid` | `REAL` | DOY / Days (`d`) | Circular Hopkins-standardized flowering day and signed phenological residual (days from conspecific expectation). |
| `taxonomic_status_call` | `TEXT` | Status call string | Diagnostic status call summarizing the active hierarchical decision stage. |
| `recommended_determination` | `TEXT` | *Packera* taxon / novel call | Recommended herbarium annotation label (reassigned *Packera* species, hybrid formula, `sp. nov. ined.`, or confirmed taxon). |
| `taxonomic_treatment_broad` | `TEXT` | Broad treatment enum | High-level curatorial routing bucket (`Confirmed_Taxon`, `Reassign_Determination`, `Interspecific_Hybrid`, `Novel_Lineage_Candidate`, `Infraspecific_or_Plastic`). |
| `taxonomic_treatment_subtype` | `TEXT` | Canonical 8-state enum | Five-stage biological determination: `Confirmed_Species`, `Confirmed_Type_Control`, `Type_Synonymy_Candidate`, `Misidentification_Reassign`, `Hybrid_Introgressant`, `Putative_Cryptic_Species`, `Subspecies_Ecotype`, or `Ecophenotypic_Plasticity`. |
| `cryptic_cluster_cohesion_pass` | `INTEGER` | `{0, 1}` | Equals `1` only when a non-core morphocluster ($g_{\text{cl}} > 1$) satisfies exsiccatae-deduplicated gathering cohesion ($\ge 5$ gatherings across $\ge 3$ 10 km spatial demes and $\ge 2$ collectors or herbaria), exact `open_set_df` congener rejection, multimodal barrier congruence, and `cryptic_non_clinal_pass == 1` (gated by `coalescent_gdi >= 0.20` or `cluster_propagated_gdi >= 0.20` when sequenced or cluster-governed). |
| `cryptic_non_clinal_pass` | `INTEGER` | `{0, 1}` | Equals `1` when a secondary morphocluster is maintained across parapatric contact ($\le 50\text{ km}$), Hopkins reproductive allochrony, or azonal lithology rather than a smooth allopatric geographic cline ($R^2 < 0.60$). |
| `resurrection_basionym` | `TEXT` | Basionym / `NULL` | Unstandardized heterotypic nomenclatural basionym matched within $\chi^2_{0.95, p}$ for Stage 3 `Taxonomic_Resurrection_Candidate` lineages (`NULL` for unnamed *de novo* cryptic candidates, cleared automatically to `NA`/`NULL` by `evaluate_step07_hydro_edaphic_and_basionym_guards()` whenever `basionym_synonymy_type == "Homotypic"`). |
| `hybrid_zone_cohesion_pass` | `INTEGER` | `{0, 1}` | Equals `1` only when a multi-species specimen satisfies all five Conservative Contact-Zone Hybrid Swarm Gates (local two-parent syntopy within `15 km`, co-flowering within `21 d`, pairwise Mahalanobis chord intermediacy $t_{AB} \in [0.22, 0.78]$ with bounded $d^2_{\perp}$, multi-organ corroboration, and $\ge 2$ admixed vouchers in the local contact zone; always `0` in single-taxon `--intraspecific-mode` unless directly supported by F1 genomic $(h, Q_{12})$). |
| `hybrid_parent_1`, `hybrid_parent_2` | `TEXT` | *Packera* binomial / `NULL` | Recorded maternal/baseline parent (`given_sp`) and verified sympatric congener parent for confirmed `Hybrid_Introgressant` swarm vouchers (`NULL` for all non-hybrid calls and in single-taxon `--intraspecific-mode`). |
| `synthesis_triage_priority` | `TEXT` | `{'CRITICAL','HIGH','MEDIUM','LOW','RESOLVED'}` | Physical herbarium loan and annotation priority (`CRITICAL` for cohesive cryptic candidates and type synonymy conflicts; `HIGH` for confident misidentifications and hybrids; `RESOLVED` for concordant anchors). |
| `synthesis_triage_action` | `TEXT` | Curatorial action code | Concrete curatorial instruction (for example `Annotate & Re-folder`, `Inspect for Hybridity / Cytotype`, `Prioritize Sequenced Loan`, or `Retain Current Determination`). |
| `synthesis_rationale` | `TEXT` | Human-readable sentence | Quantitative biological explanation citing the individual's CDA posterior, OSR $D^2_{\min}$, rhizosphere $Z$-score, lithology, $\Delta\text{DOY}$ residual, Step 06 genomic `gdi` or $(h, Q_{12})$ evidence, Step 07 genomic propagation or `[Tier_0 Abstention: ...]` rationales, quality-gated cytonuclear capture provenance (`cpDNA=...`), confirmed Stage 2 hybrid maternal cpDNA polarity (`[Maternal cpDNA Parent: ...]`), and non-mutating out-of-panel `[Out-of-Panel Nuclear Locus Match: ITS+waxy=...]` notes. |
| `genomic_propagation_tier` | `TEXT` | Tier enum or `NULL` | Unsequenced genomic propagation provenance (`G1_Exsiccatae_Isoduplicate`, `G2_Conditional_Deme_Hull`, `Tier0_Projected_Calibrated`, or `Cluster_GDI_Governed`) that remains `NULL` when unpropagated or when genomic sample inputs are absent; `G2_Conditional_Deme_Hull` gates strictly on observed cytotype conflicts and verified micro-sympatry anchors. |
| `tier0_calibrated_taxon` | `TEXT` | Taxon, abstention tag, or `NULL` | Leave-one-population-out Ledoit-Wolf shrinkage discriminant prediction trained on `Tier_0` DNA-verified anchors (`weight = 4.0`), or `Sensu_Lato_DNA_Required` when any abstention condition fires. |
| `tier0_calibrated_prob` | `REAL` | Probability $[0, 1]$ | Temperature-scaled softmax posterior probability supporting `tier0_calibrated_taxon` under the `Tier_0`-calibrated shrinkage discriminant model (`NULL` when uncalibrated). |
| `tier0_abstention_flag` | `INTEGER` | `{0, 1}` or `NULL` | Equals `1` when phenotypic overlap, low margin, high `LOPO-CV` error, or out-of-envelope distance triggers `Sensu_Lato_DNA_Required` abstention, whereas `0` records confident `Tier_0` projection and `NULL` marks inactive calibration. |
| `deme_mahal_d2` | `REAL` | Squared distance $D^2 \ge 0$ | Squared Mahalanobis morphometric distance from an unsequenced `Basal_Rosette` specimen to its nearest candidate sequenced `2.0 km` deme anchor taxon centroid. |
| `geno_deme_weight_scale` | `REAL` | $[0.45, 0.75]$ or `NULL` | Soft reliability multiplier assigned to `G2_Conditional_Deme_Hull` recipients from baseline `0.75` after multiplicative attenuation for `1000 to 2000 m` precision or `NA` precision beyond `1.0 km` (`0.85x`) alongside unobserved multi-cytotype ploidy (`0.85x`) plus outer-buffer soil or anthesis offsets (`0.88x` each), clamped to `[0.45, 0.75]`. |
| `cluster_propagated_gdi` | `REAL` | $[0, 1]$ or `NULL` | Genealogical divergence index propagated from sequenced `G1` or `G2` cluster anchors across non-abstained, non-admixed cluster members unless observed ploidies conflict within `coancestry_cluster` or `2.0 km` or against the recipient; blocks Stage 3 oversplits below `0.20` and corroborates cohesive novel lineages above `0.70`. |
| `cytology_linkage_tier` | `TEXT` | Tier enum or `NULL` | Voucher-to-cytology provenance tier (`C1_Direct_Cytovoucher`, `C2_Exsiccatae_Cytoduplicate`, `C3_Taxon_Monomorphic`, or `C4_Taxon_MultiCytotype`). |
| `voucher_sporophytic_2n` | `INTEGER` | Chromosomes or `NULL` | Unanimously corroborated somatic chromosome number $2n$ resolved from direct `C1` or duplicate `C2` cytovoucher records. |
| `voucher_ploidy_level` | `INTEGER` | Integer ploidy or `NULL` | Ploidy level (`2..16`) assigned strictly from `C1`/`C2` cytovoucher counts or high-support Step 06 genomic cytotype calls (`|\Delta\ell| \ge 4.0`). |
| `taxon_reported_cytotypes` | `TEXT` | Series string or `NULL` | Comma-separated cytotype series documented for the sheet's standardized taxon across curated cytological summaries and typed cohort vouchers. |
| `is_multi_cytotype_taxon` | `INTEGER` | `{0, 1}` or `NULL` | Equals `1` when the taxon exhibits multiple ploidy levels or dysploid base numbers so Step 07 gates cross-cytotype deme propagation. |
| `cytotype_reproductive_barrier_pass` | `INTEGER` | `{0, 1}` or `NULL` | Equals `1` when a secondary GMM cluster carries at least two independent typed anchors that hold disjoint ploidy or dysploid base $x$ relative to the primary type-anchored cluster. |
| `hydro_edaphic_niche_regime` | `TEXT` | Regime enum or `NULL` | Step 05/07 hydro-edaphic habitat classification derived from USGS `DS-801` Ca:Mg/Ni-Cr/carbonate C chemistry and USFWS NWI wetland state (e.g., `Ultramafic_Serpentine_Barren`, `Calcareous_Glade`, `Palustrine_Wetland`). |
| `wetland_fidelity_concordance` | `REAL` | Concordance $[0, 1]$ | Quantitative agreement between site NWI wetland occupancy (`is_nwi_wetland` / distance) and taxon USACE NWPL wetland fidelity score (`nwpl_wetland_score`). |
| `basionym_synonymy_type` | `TEXT` | Synonymy class or `NULL` | Compares candidate `resurrection_basionym` against `taxon_backbone_profiles` homotypic/heterotypic JSON arrays (`Homotypic` vs `Heterotypic`); clears homotypic combinations out of `resurrection_basionym`. |
| `is_completed_voucher` | `INTEGER` | `{0, 1}` | Equals `1` when the individual completed integrative evaluation without upstream QC rejection. |

### 3.7 `taxon_complexes` (Step 07 Unsupervised Species-Complex Delimitation)

**Primary Key:** `(run_id, species_standardized)` | **Populated By:** `scripts/analysis/taxon_complex_discovery_utils.R` paired with `scripts/analysis/step07_complex_affinity_utils.R`
**Companion CSV:** `runs/<RUN_ID>/outputs/reports/empirical_taxon_complexes.csv` (exact 23-column twin-free schema excluding `run_id`)

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `run_id` | `TEXT` | Run ID *(SQLite only)* | Foreign key to `pipeline_runs.run_id`. |
| `species_standardized` | `TEXT` | Taxon or `[K#]` subcluster | Evaluated *Packera* taxon (or `P. <epithet> [K#]` Taxon $\times$ GMM Morphotype subcluster when $N_{\text{taxa}} \le 2$). |
| `n_vouchers` | `INTEGER` | Count $\ge 1$ | Total evaluated individuals assigned to this taxon node. |
| `n_gold_silver_anchors` | `INTEGER` | Count $\ge 0$ | Count of active reference anchors (`Tier_0_Type_Anchor`, `Tier_1_Gold`, `Tier_2_Silver`, or promoted fallback anchors) defining the taxon centroid. |
| `n_type_specimens` | `INTEGER` | Count $\ge 0$ | Count of verified ICN nomenclatural type specimens shrinking the taxon centroid with weight $\gamma_k$. |
| `consensus_complex_id` | `INTEGER` | Sequential integer $\ge 1$ | Consensus species-complex cluster ID selected between Leiden modularity and metric Ward.D2 silhouette optimization. |
| `complex_label` | `TEXT` | `Complex_<ID>_<Ep1>_<Ep2>` | Deterministic complex name constructed from the top two highest-degree constituent epithets. |
| `leiden_complex_id` | `INTEGER` | Integer $\ge 1$ | Community assignment from multi-resolution Leiden modularity (`igraph::cluster_leiden`, with Louvain fallback). |
| `ward_complex_id` | `INTEGER` | Integer $\ge 1$ | Cluster assignment from silhouette-optimal Ward.D2 hierarchical clustering on $D_{ij} = \sqrt{1 - W_{ij}}$. |
| `silhouette_width` | `REAL` | $[-1, 1]$ | Taxon silhouette width on the metric multimodal dissimilarity matrix $D_{ij}$. |
| `within_complex_mean_affinity` | `REAL` | Affinity $[0, 1]$ | Mean multimodal affinity $W_{ij}$ between this taxon and all other members of its `consensus_complex_id` (`0` for singleton complexes). |
| `nearest_neighbor_taxon`, `nearest_neighbor_affinity` | `TEXT` / `REAL` | Taxon / $[0, 1]$ | Highest-affinity congeneric neighbor in the cohort and its composite multimodal affinity $W_{i,(1)}$. |
| `second_neighbor_taxon`, `second_neighbor_affinity` | `TEXT` / `REAL` | Taxon / $[0, 1]$ | Second-closest congeneric neighbor and its composite multimodal affinity $W_{i,(2)}$. |
| `hopkins_pheno_overlap` | `REAL` | Overlap $[0, 1]$ | Circular Hopkins flowering overlap $S^{\text{pheno}}_{i,\text{nn1}}$ with `nearest_neighbor_taxon` (`NA` when verified-anthesis phenology is unavailable). |
| `eco_schoener_overlap` | `REAL` | Overlap $[0, 1]$ | Standardized environmental / rhizosphere overlap $S^{\text{env}}_{i,\text{nn1}}$ with `nearest_neighbor_taxon` (`NA` when environmental data are absent). |
| `integrative_congruence` | `INTEGER` | `{0, 1}` (`logical` in CSV) | Equals `TRUE` (`1`) when $\ge 2$ modalities show non-trivial affinity ($S^{\text{morph}} \ge \tau_{\min}$, $S^{\text{pheno}} \ge 0.25$, $S^{\text{env}} \ge 0.25$) with the nearest neighbor. |
| `w_morph_kappa`, `w_edaphic_kappa`, `w_pheno_kappa` | `REAL` | Weights summing to $1.0$ | Cohort-calibrated OOF Cohen's $\kappa$ weights combining $S^{\text{morph}}_{ij}$, $S^{\text{env}}_{ij}$, and $S^{\text{pheno}}_{ij}$ into $W_{ij}$. |
| `n_misidentification_reassign`, `n_hybrid_introgressant`, `n_putative_cryptic` | `INTEGER` | Counts $\ge 0$ | Specimen-level anomaly counts tallied within this taxon during complex evaluation. |

## 4. Relational Master View and Curatorial Triage Exports

### 4.1 `v_integrative_triage_master` (Composite Seven-Stage SQL View)

Defined identically in `scripts/core/db_schema.py` and `scripts/core/db_utils.R`, `v_integrative_triage_master` selects from `run_active_cohort` left-joined onto `vouchers`, `plant_individuals`, `voucher_reproductive_phenology` (applying `plant_individual_id = 0` sheet-level fallback `vrp0` only when an individual row is absent), `otu_morphometrics`, `spatial_rf_predictions`, `environmental_cache` (via `lat_lon_key`), `taxon_backbone_profiles` (exposing `backbone_basionym_name`, `protologue_url`, `natureserve_global_rank`, `nwpl_indicator_code`, and `nwpl_wetland_score`), a `ROW_NUMBER()`-deduplicated join of `genetic_lineage_assignments` with `genomic_samples` (ranked by `G1_Direct_Voucher` > `G2_Population_Voucher_Proxy` > QC pass > HQ reads > posterior), and `curation_triage` on `(run_id, catalogNumber, plant_individual_id)`, and the view exposes `is_promoted_bronze`, `blade_cie_l_sheet_delta`, `blade_cie_a_delta`, `blade_cie_b_delta`, `ecotype_regime`, `lithic_outcrop_index`, `summer_moisture_deficit_idx`, `maritime_continentality_idx`, `soil_ca_mg_molar_ratio`, `soil_ultramafic_ni_cr_index`, `soil_carbonate_c_pct`, `is_nwi_wetland`, `nwi_wetland_system`, `nwi_wetland_dist_m`, `n_genomic_samples_linked`, `inferred_cytotype`, `hybrid_generation_class`, `coalescent_gdi`, `public_plastome_match`, `public_its_match`, `waxy_lineage_match`, `public_backbone_clade`, `public_phytop_ils_index`, `public_phytop_ih_index`, `genomic_propagation_tier`, `tier0_calibrated_taxon`, `tier0_calibrated_prob`, `tier0_abstention_flag`, `deme_mahal_d2`, `cluster_propagated_gdi`, `cytology_linkage_tier`, `voucher_sporophytic_2n`, `voucher_ploidy_level`, `taxon_reported_cytotypes`, `is_multi_cytotype_taxon`, `cytotype_reproductive_barrier_pass`, `hydro_edaphic_niche_regime`, `wetland_fidelity_concordance`, `basionym_synonymy_type`, `morph_discordance_reason`, and `conflict_discordance_notes`.

Three explicit precedence rules prevent stale values during partial re-runs (`--start-step 3` or `--start-step 5`):
1. `run_active_cohort` overrides raw `vouchers` columns for `species_standardized`, `determiner_tier`, `type_status`, `is_promoted_bronze`, and `spatial_balance_weight` via `COALESCE(rac.<col>, ...)` so dynamic authority or synonym updates apply without re-harvesting images.
2. Step 05 `spatial_rf_predictions` BioClim variables take precedence over `environmental_cache` via `COALESCE(srf.bio*, ec.bio*)` while `environmental_cache` supplies fine-grain soil, 3DEP, and SGMC columns.
3. Step 05 `spatial_rf_predictions.hopkins_adjusted_doy` takes precedence over `curation_triage.hopkins_adjusted_doy` via `COALESCE(srf.hopkins_adjusted_doy, ct.hopkins_adjusted_doy)` so re-running Step 05 immediately propagates updated Hopkins phenology into Step 07.

### 4.2 `triage_queue.csv`, `genomic_field_triage.csv`, and Interactive Dashboard Filtered CSV

* **`runs/<RUN_ID>/data/tables/triage_queue.csv`:** Full wide-table export written atomically by `build_triage_queue()` in `scripts/analysis/07_triage_dashboard_synthesis.R`, sorted by `synthesis_triage_priority` (`CRITICAL -> HIGH -> MEDIUM -> LOW -> RESOLVED`) and descending `coordinateUncertainty`. `build_triage_queue()` also hydrates 11-mer organellar placement metrics (`plastid_kmer_containment`, `plastid_shared_11mers`, `its_kmer_containment`, `its_shared_11mers`, `waxy_kmer_containment`, `waxy_shared_11mers`) from `public_organellar_placement.csv` via `step07_organellar_triage_utils.R`. When `--focal-taxon` is supplied, SQLite `curation_triage` and `taxon_complexes` persist the complete multi-taxon reference cohort first, whereas `triage_queue.csv` filters to specimens originally filed under or reassigned into `--focal-taxon`.
* **`runs/<RUN_ID>/data/tables/genomic_field_triage.csv`:** Stage 0 DNA-plate and unvouchered field-sample triage table exported by Step 07 (`scripts/analysis/genomic_triage_utils.R`) whenever Step 06 genomic predictions are staged, summarizing G1 to G4 voucher linkage, cytotype, hybrid class, $\widehat{gdi}$, and prioritized herbarium loan requests for `G3`/`G4` unvouchered genetic tips. Four-tier unsequenced genomic propagation counts (`n_g1_exsiccatae`, `n_g2_deme_hull`, `n_tier0_anchors`, `n_tier0_evicted`, `n_tier0_calibrated`, `n_tier0_abstained`, `n_cluster_gdi`, and `lopo_cv_kappa`) attach via `attr(df, "genomic_propagation_summary")` (`scripts/analysis/step07_genomic_propagation.R`) and render into `Packera_Taxonomic_Revision_Summary.md` via `build_genomic_markdown_summary()`. When propagation is active, `build_triage_queue()` also writes the per-taxon propagation summary table to `runs/<RUN_ID>/outputs/reports/unsequenced_genomic_propagation_summary.csv`. Whenever multi-cytotype specimens exist across Step 07 cohorts, `run_step07_cytology_integration()` (`scripts/analysis/step07_cytology_integration.R`) writes the post-hoc Gigas trait audit across `area_mm2` alongside `petiole_width_mm` and `log1p_capitulum_count` to `runs/<RUN_ID>/outputs/reports/cytotype_concordance_summary.csv` without circular feedback into morphometrics or consensus training weights.
* **`triage_queue_filtered.csv` (Browser Export from `Triage_Dashboard_Interactive.html`):** Exports the active curatorial bench columns (`catalogNumber`, `plant_individual_id`, `institutionCode`, `recordedBy`, `determiner_tier`, `species_standardized`, `recommended_determination`, `synthesis_triage_priority`, `taxonomic_treatment_subtype`, `cda_predicted_taxon`, `cda_posterior_prob`, `edaphic_best_fit_taxon`, `edaphic_z_score`, `delta_doy_resid`, `probabilistic_consensus_score`, `synthesis_triage_action`, `synthesis_rationale`) matching the active search, priority, treatment, and taxon filters in the offline HTML viewer.
