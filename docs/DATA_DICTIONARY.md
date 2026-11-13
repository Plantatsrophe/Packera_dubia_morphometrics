# Relational Schema and Output Data Dictionary (`PACKERA_SCHEMA_USER_VERSION = 20261121`)

**Standard Operating Procedure:** [`docs/TAXONOMIC_SOP.md`](TAXONOMIC_SOP.md) (`UNC-BOT-SOP-2026-04-REV5`)
**Authoritative Schema Definitions:** `scripts/core/db_schema.py` and `scripts/core/db_utils.R`

## 1. Architectural Conventions and Key Invariants

Herbarium sheets frequently carry multiple pressed plants (`ramets`) mounted under a single barcode, whereas repeated exploratory runs over different species subsets or authority filters must never overwrite raw computer vision extractions. I therefore split pipeline storage across two transactional tiers inside `data/packera_pipeline.sqlite` (Write-Ahead Logging mode, `PRAGMA user_version = 20261121`) paired with atomic CSV snapshots on disk.

Tier A (`data/`) stores invariant specimen metadata, full-sheet WebP scans, segmented leaf masks, 120-vertex contours, $32 \times 360$ Euler Characteristic Transform (ECT) matrices, and spatially cached environmental extractions. Tier B (`runs/<RUN_ID>/` and `runs/latest`) isolates run-specific cohort memberships, Elliptic Fourier harmonics, canonical Gaussian Mixture Models, Ledoit-Wolf shrinkage discriminant scores, Spatial Random Forest posteriors, and five-stage curatorial triage calls.

Three cross-language invariants govern every column documented below:

* **Composite Individual Key `(run_id, catalogNumber, plant_individual_id)`:** A single barcode (`catalogNumber`) can yield multiple DBSCAN-clustered individuals (`plant_individual_id = 1, 2, ...`, or `0` when no whole-plant cluster is resolved). Every Python and R join across Steps 02 through 06 groups strictly on `(catalogNumber, plant_individual_id)` so multi-plant sheets never collapse to a single row.
* **Native Pixel vs. Calibrated Millimeter Coordinates:** All `_px` columns represent coordinates on the original unscaled full-sheet scan (`coord_work / work_scale`). Columns ending in `_mm`, `_cm`, or `_mm2` contain physical measurements only when `is_calibrated_mm = 1` (`scale_source` in `{'ruler_calibrated', 'fieldprism_calibrated'}`), whereas uncalibrated sheets (`mp_estimated`) populate dimensionless ratios (`aspect_ratio`, `solidity`, `petiole_to_blade_ratio`, `AR_inv`) and unit-centroid-size EFA/ECT descriptors while leaving physical linear columns honestly `NULL` (`NA`) in strict scale-gated models.
* **Honest Missingness Outside Demo Mode:** When `--demo` is inactive, missing environmental rasters, high-uncertainty coordinates exceeding layer gates, or absent floral heads record `NULL` in SQLite and `NA` in exported CSVs rather than imputed placeholders. Running with `--demo` (`PACKERA_DEMO_MODE=1`) redirects all database writes to an ephemeral SQLite file (`demo_ephemeral.sqlite`) so canonical Tier A tables remain untouched.

## 2. Tier A Shared Canonical SQLite Tables (`data/packera_pipeline.sqlite`)

### 2.1 `taxon_synonyms` (Step 01 Dynamic Basionym Cache)

**Primary Key:** `raw_name_lower` | **Populated By:** `scripts/core/harvester.py`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `raw_name_lower` | `TEXT` | Lowercase binomial/trinomial | Normalized raw taxon string queried from GBIF, iDigBio, or historical *Senecio* and *Cineraria* labels. |
| `canonical_packera` | `TEXT` | *Packera* taxon string | Accepted *Packera* binomial or infraspecific taxon (`subsp.` / `var.`) resolved via GBIF Backbone and Weakley (2026). |
| `protologue_authors_json` | `TEXT` | JSON string array | Parsed basionym and combination author surnames used by the monographer loan guard in Step 01. |
| `updated_at` | `TEXT` | ISO-8601 UTC timestamp | Timestamp of the most recent taxonomic backbone verification. |

### 2.2 `usda_county_ranges` (Step 01 County Distribution Cache)

**Primary Key:** `(canonical_packera, state_norm, county_norm)` | **Populated By:** `scripts/core/usda_plants_ranges.py`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `canonical_packera` | `TEXT` | *Packera* taxon string | Standardized *Packera* species or variety name. |
| `usda_symbol` | `TEXT` | USDA PLANTS code | Accepted USDA PLANTS symbol (for example `PADU` or `PAAU3`). |
| `state_norm` | `TEXT` | 2-letter postal code | Standardized US state or Canadian province abbreviation. |
| `county_norm` | `TEXT` | Normalized county name | Lowercase county or parish name stripped of administrative suffixes (`""` when `is_state_only = 1`). |
| `fips_code` | `TEXT` | 5-digit FIPS string | US county FIPS identifier when available. |
| `native_status` | `TEXT` | Status code | USDA PLANTS nativity string (`N`, `I`, or `W`). |
| `is_state_only` | `INTEGER` | `{0, 1}` | Equals `1` when USDA PLANTS records presence only at the state level without county breakdown. |
| `updated_at` | `TEXT` | ISO-8601 UTC timestamp | Cache refresh timestamp. |

### 2.3 `vouchers` (Step 01 Canonical Specimen Registry)

**Primary Key:** `catalogNumber` | **Populated By:** `scripts/data_prep/01_voucher_harvester.py` and Step 02 calibration updates
**Companion CSV:** `curated_vouchers.csv` exports the 26 core Darwin Core, authority, ecoregion, and duplicate-audit columns from this table for R consumption.

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `catalogNumber` | `TEXT` | Herbarium barcode / UUID | Unique specimen identifier (`institutionCode` + barcode or occurrence UUID fallback). |
| `institutionCode` | `TEXT` | Index Herbariorum code | Source herbarium code (for example `NCU`, `US`, `NY`, `MO`, `GH`, `DUKE`). |
| `scientificName` | `TEXT` | Darwin Core string | Rawverbatim scientific name with authorship as served by GBIF or iDigBio. |
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
| `spatial_balance_weight` | `REAL` | Dimensionless $(0, 1]$ | Inverse spatial density weight moderating over-collected counties in downstream models. |
| `year`, `month`, `day` | `INTEGER` | Calendar integers | Parsed collection year, month (`1` to `12`), and day (`1` to `31`). |
| `eventDate` | `TEXT` | ISO-8601 date | Standardized collection date (`YYYY-MM-DD`). |
| `doy` | `INTEGER` | Day of Year (`1` to `366`) | Ordinal collection day; set to `NULL` when month/day are absent or represent `YYYY-01-01` placeholder artifacts. |
| `image_path` | `TEXT` | Relative file path | Local path to the lossless/high-Q full-sheet WebP scan under `data/raw_vouchers/`. |
| `exsiccatae_key` | `TEXT` | Hash / normalized key | Gathering deduplication key combining normalized primary collector surname, `recordNumber`, and `eventDate`. |
| `is_primary_duplicate` | `INTEGER` | `{0, 1}` | Equals `1` for the highest-resolution / highest-authority sheet retained from a duplicate gathering set. |
| `duplicate_count` | `INTEGER` | Count $\ge 1$ | Total duplicate sheets detected across herbaria for `exsiccatae_key`. |
| `exsiccatae_taxonomic_conflict` | `INTEGER` | `{0, 1}` | Equals `1` when duplicate sheets of the exact same gathering carry conflicting *Packera* determinations across herbaria. |
| `exsiccatae_gold_taxon` | `TEXT` | *Packera* taxon string | Consensus `Tier_1_Gold` or `Tier_0_Type_Anchor` determination propagated across a split duplicate gathering when present. |
| `usda_symbol` | `TEXT` | USDA PLANTS code | Matched USDA PLANTS symbol for `species_standardized`. |
| `usda_county_status` | `TEXT` | Range audit enum | `in_range`, `buffered_50km`, `out_of_range`, `state_only_match`, or `unverified` (`50 km` county border tolerance; audit-only by default). |
| `pixels_per_mm` | `REAL` | $\text{px}\cdot\text{mm}^{-1}$ | Sheet optical resolution updated by Step 02 ruler tick-lattice or FieldPrism calibration. |
| `is_calibrated_mm` | `INTEGER` | `{0, 1}` | Equals `1` when `pixels_per_mm` derives from physical ruler or FieldPrism ticks rather than megapixel estimation. |
| `qc_status` | `TEXT` | `{'PASS','FAILED_QC','PENDING'}` | Step 02 specimen processing status (`PENDING` after Step 01 harvest; `PASS` or `FAILED_QC` after Step 02). |
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
| `petiole_to_blade_ratio` | `REAL` | Dimensionless $\ge 0$ | Scale-invariant ratio `petiole_length_px / blade_length_px` diagnostic in the *P. aurea* and *P. paupercula* complexes. |
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
| `scale_source` | `TEXT` | Calibration source enum | `ruler_calibrated`, `fieldprism_calibrated`, or `mp_estimated`. |
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
**Populated By:** `scripts/analysis/edaphic_fine_grain_utils.R` and `scripts/analysis/05_multimodal_spatial_rf.R`

| Column Group | Column Names | SQLite Type | Units / Domain | Description |
| :--- | :--- | :---: | :--- | :--- |
| **Keys & Coordinates** | `lat_lon_key`, `decimalLatitude`, `decimalLongitude`, `updated_at` | `TEXT` / `REAL` | `WGS84` deg | Rounded coordinate cache key, exact decimal coordinates, and extraction timestamp. |
| **Uncertainty & Eligibility Gates** | `buffer_radius_m`, `high_coord_uncertainty_flag`, `is_training_eligible_coord`, `is_clim_eligible_coord`, `is_edaphic_eligible_coord`, `edaphic_coord_weight` | `REAL` / `INTEGER` | `m` / `{0, 1}` / $[0, 1]$ | Extraction radius (`50 m` default), high-uncertainty flag ($> 2000\text{ m}$), separate climatic ($\le 10000\text{ m}$) vs. edaphic ($\le 2000\text{ m}$) eligibility gates, and inverse coordinate-uncertainty weight. |
| **SSURGO Map Units (`10 to 30 m`)** | `musym`, `muname`, `taxsuborder`, `taxgreatgroup`, `is_rock_outcrop`, `outcrop_class`, `ssurgo_rock_outcrop_pct`, `ssurgo_bedrock_depth_cm`, `ssurgo_sand_pct`, `ssurgo_clay_pct`, `ssurgo_ph`, `ssurgo_om_pct`, `ssurgo_cec`, `ssurgo_aws_0_50cm` | `TEXT` / `INTEGER` / `REAL` | `%`, `cm`, `pH`, $\text{cmol}_+\cdot\text{kg}^{-1}$, `cm` | USDA-NRCS SSURGO map unit taxonomy, `outcrop_class` (`Granitic flatrock`, `Limestone glade`, `Sandstone/Siliceous`, `Serpentine`, `Diabase/Mafic glade`, `Non-outcrop matrix`), lithic contact depth, and $0\text{ to }15\text{ cm}$ horizon properties. |
| **POLARIS (`30 m`)** | `polaris_ph_30m`, `polaris_sand_30m`, `polaris_clay_30m`, `polaris_om_30m`, `polaris_bd_30m`, `polaris_ksat_30m` | `REAL` | `pH`, `%`, $\text{g}\cdot\text{cm}^{-3}$, $\text{cm}\cdot\text{hr}^{-1}$ | Optional $30\text{ m}$ probabilistic rhizosphere properties integrated over $0\text{ to }15\text{ cm}$ when local POLARIS tiles are staged. |
| **SOLUS100 (`100 m`)** | `solus_resdept_cm`, `solus_fragvol_pct`, `solus_sand_pct`, `solus_clay_pct`, `solus_silt_pct`, `solus_cec` | `REAL` | `cm`, `%`, $\text{cmol}_+\cdot\text{kg}^{-1}$ | USDA-NRCS $100\text{ m}$ depth to restriction/lithic contact, coarse fragments, particle-size fractions, and CEC over $0\text{ to }15\text{ cm}$. |
| **Harmonized Rhizosphere (`0 to 15 cm`)** | `fine_soil_ph`, `fine_soil_sand`, `fine_soil_clay`, `soil_silt_pct`, `ilr_texture_1`, `ilr_texture_2`, `fine_soil_cec`, `fine_cfvo_pct`, `fine_bedrock_depth_cm` | `REAL` | `pH`, `%`, ILR, $\text{cmol}_+\cdot\text{kg}^{-1}$, `cm` | Priority-harmonized (`POLARIS -> SOLUS100 -> SSURGO -> SoilGrids`) $0\text{ to }15\text{ cm}$ rhizosphere predictors and two orthogonal Isometric Log-Ratio (`ILR`) texture coordinates freeing sand/silt/clay from the $100\%$ unit-sum constraint. |
| **50 m Buffer Quantiles** | `ph_q10_buffer`, `ph_q90_buffer`, `min_bedrock_depth_buffer_cm`, `max_cfvo_buffer_pct`, `azonal_lithology_presence_buffer` | `REAL` / `INTEGER` | `pH`, `cm`, `%`, `{0, 1}` | Sub-grid micro-habitat quantiles within `50 m` capturing narrow rock outcrops embedded inside forest soil matrices. |
| **3DEP Micro-Topography (`30 m`)** | `elevation_m`, `topo_slope_deg`, `topo_tpi_30m`, `topo_hli_30m`, `topo_twi_30m` | `REAL` | `m`, `deg`, `m`, $[0, 1]$, $\ln(a/\tan\beta)$ | Elevation, slope, Topographic Position Index (`TPI`), McCune-Keon Heat Load Index (`HLI`), and D8 Topographic Wetness Index (`TWI`). |
| **USGS SGMC Bedrock** | `sgmc_lithology_class`, `is_azonal_lithology` | `TEXT` / `INTEGER` | 6-class enum / `{0, 1}` | Standardized bedrock class (`Ultramafic_Serpentine`, `Mafic_Igneous_Metamorphic`, `Carbonate_Calcareous`, `Felsic_Granitic_Gneiss`, `Siliciclastic_Consolidated`, `Unconsolidated_Sedimentary`) and azonal indicator (`1` for ultramafic, mafic, or carbonate). |
| **WorldClim v2.1 (`1 km`)** | `bio1_temp_mean` through `bio19_precip_coldest_qtr` | `REAL` | $^\circ\text{C}$, `%`, `mm` | All 19 standard bioclimatic temperature and precipitation rasters (`bio1` to `bio19`) gated at `<= 10000 m`. |
| **Coarse Fallbacks** | `soil_ph`, `soil_clay`, `soil_sand`, `soil_nitrogen`, `soil_soc`, `soil_cec`, `soil_cfvo`, `soil_bdod`, `*_soilgrids` | `REAL` | Standard pedological units | Coarse SoilGrids / legacy compatibility columns retained for non-CONUS occurrences outside SSURGO and SOLUS100 coverage. |

## 3. Tier B Run-Scoped SQLite Tables (`runs/<RUN_ID>/` & SQLite)

### 3.1 `pipeline_runs` and `run_active_cohort` (Preflight Cohort Isolation)

* **`pipeline_runs` (PK `run_id`):** Tracks `created_at`, `completed_at`, `status` (`RUNNING`, `COMPLETED`, `FAILED`), `git_commit`, `cli_args_json`, `config_hash`, `focal_taxon`, `run_mode`, and `run_dir` (`runs/<RUN_ID>`).
* **`run_active_cohort` (PK `(run_id, catalogNumber, plant_individual_id)`):** Snapshots every individual included in `<RUN_ID>` along with its preflight-refreshed `species_standardized`, `determiner_tier`, `type_status`, `is_type_specimen`, `authority_weight` (`3.0` Type, `2.0` Gold, `0.75` Silver, `0.0` unverified/unpromoted Bronze), and `spatial_balance_weight`.

### 3.2 `leaf_efa_harmonics` (Step 03 Individual-Leaf EFA & ECT Table)

**Primary Key:** `(run_id, catalogNumber, plant_individual_id, leaf_id)` | **Populated By:** `scripts/morphometrics/03_fourier_extractor.R`
**Companion CSV:** `runs/<RUN_ID>/data/tables/leaf_efa_harmonics_individual.csv`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `run_id`, `catalogNumber`, `plant_individual_id`, `leaf_id` | `TEXT` / `INTEGER` | Composite PK | Run, specimen, individual, and leaf keys. |
| `shape_id`, `leaf_zone` | `TEXT` | Organ ID / `{'basal','cauline'}` | Unique organ key and architectural shoot zone. |
| `species_standardized`, `determiner_tier`, `type_status`, `is_type_specimen`, `authority_weight` | `TEXT` / `INTEGER` / `REAL` | Authority metadata | Taxon hypothesis and authority metadata synchronized from `run_active_cohort`. |
| `A1..A24`, `B1..B24`, `C1..C24`, `D1..D24` | `REAL` | Dimensionless EFA | 96 Elliptic Fourier coefficients (`Momocs::efourier(..., norm = FALSE, start = FALSE)` up to $H_{\max} = 24$; harmonics above $H_{\text{opt}}$ are `0` or `NULL`). |
| `PC1..PC15` | `REAL` | Dimensionless EFA PCs | Authority-weighted (`stats::cov.wt`) symmetric EFA principal component scores ($A_n, D_n$ lateral sub-matrix). |
| `sym_var`, `asym_var`, `pct_asymmetry` | `REAL` | Power / `%` $[0, 100]$ | Symmetric ($A_n^2 + D_n^2$) vs. asymmetric ($B_n^2 + C_n^2$) harmonic power and bilateral asymmetry percentage (`> 30%` flags pressing/taphonomic distortion). |
| `harmonic_roughness` | `REAL` | Dimensionless $[0, 1]$ | High-frequency harmonic energy fraction ($n = 9..\min(H_{\text{opt}}, 24)$) screening jagged segmentation tears. |
| `is_contour_outlier` | `INTEGER` | `{0, 1}` | Equals `1` when leaf fails Step 03 R-side clump/boundary QC, high-frequency roughness, or within-taxon FastMCD ($p < 0.001$; `0` for protected types). |
| `leaf_quality_weight` | `REAL` | Weight $[0, 1]$ | Continuous Tukey bisquare quality weight used in weighted median OTU aggregation. |
| `centroid_size`, `log_centsize` | `REAL` | `mm` (or `px`) / $\log_{10}$ | Contour centroid size and base-10 log centroid size used in Step 04 `RRPP` allometry. |
| `efd_harmonic_opt`, `efd_k95_pcs` | `INTEGER` | $[6, 24]$ / $[5, 15]$ | Cohort-calibrated cumulative $\ge 99\%$ Fourier harmonic cutoff ($H_{\text{opt}}$) and $\ge 95\%$ PCA subspace rank ($k_{95}$). |
| `ECT_PC1`, `ECT_PC2`, `ECT_PC3` | `REAL` | Dimensionless SVD scores | Top 3 authority-weighted SVD axes of the $32 \times 360$ polar Euler Characteristic Transform capturing pinnatifid and lyrate sinus topology without perimeter parameterization drift. |
| `ect_lobe_complexity` | `REAL` | Dimensionless $\ge 1$ | Mean positive topological component multiplicity across radial ECT slices (`1.0` for entire blades; elevated in pinnatifid taxa such as *P. millefolium* and *P. glabellus*). |
| `ect_bilateral_asym` | `REAL` | Dimensionless $\ge 0$ | Polar reflectional L2 discrepancy across the $32 \times 360$ ECT field. |

### 3.3 `otu_morphometrics` (Step 03 OTU Medians & Step 04 GMM, CDA, Allometry & OSR)

**Primary Key:** `(run_id, catalogNumber, plant_individual_id)` | **Populated By:** `03_fourier_extractor.R` and `04_gmm_morphotools.R`
**Companion CSVs:** `leaf_efa_harmonics.csv` (Step 03 OTU output) and `morphometric_flags.csv` (Step 04 enriched OTU output)

| Column Group | Column Names | SQLite Type | Units / Domain | Description |
| :--- | :--- | :---: | :--- | :--- |
| **OTU Provenance & Demes** | `morph_tier_source`, `population_deme_id`, `leaf_count`, `basal_leaf_count`, `cauline_leaf_count`, `is_passive`, `authority_weight`, `scale_source`, `is_calibrated_mm`, `pixels_per_mm` | `TEXT` / `INTEGER` / `REAL` | Zone / Deme ID / Counts | Indicates whether primary EFA/ECT medians come from `Basal_Rosette` (`basal_leaf_count >= 1`) or `Cauline_Only` fallback, $10\text{ km}$ conspecific spatial deme ID, leaf counts, and active (`is_passive = 0`) vs. passive (`1`) projection status. |
| **Pure Basal Shape & Size** | `A1..A24`, `D1..D24`, `PC1..PC15`, `ECT_PC1..ECT_PC3`, `ect_lobe_complexity`, `ect_bilateral_asym`, `solidity`, `q25_solidity`, `circularity`, `extent`, `aspect_ratio`, `dissection_index`, `dissected_leaf_fraction`, `centroid_size`, `log_centsize`, `blade_length_mm`, `blade_width_mm`, `petiole_length_mm`, `petiole_to_blade_ratio`, `plant_height_cm`, `scape_to_basal_leaf_ratio` | `REAL` | Dimensionless & `mm`/`cm` | Quality-weighted basal rosette medians (`A_n, D_n` symmetric harmonics only; `B_n, C_n` excluded from OTU aggregation as taphonomic bending terms), lower-quartile solidity `q25_solidity`, perimeter dissection index $P / (2\sqrt{\pi A})$, and linear/ratio architecture. |
| **Cauline & Heteroblasty Track** | `cauline_solidity_med`, `cauline_dissection_index_med`, `cauline_aspect_ratio_med`, `cauline_blade_length_mm`, `cauline_ect_lobe_complexity`, `basal_vs_cauline_length_ratio`, `basal_vs_cauline_solidity_delta`, `heteroblasty_efa_distance`, `heteroblasty_ect_distance`, `Cauline_PC1..Cauline_PC3`, `Cauline_ECT_PC1..Cauline_ECT_PC2` | `REAL` | Dimensionless & `mm` | Orthogonal stem-leaf morphology and within-plant basal-to-cauline developmental divergence (`NULL` unless both basal and cauline leaves occur on the individual). |
| **Reproductive Geometry** | `AR_inv`, `AR_inv_std`, `reproductive_trait_status` | `REAL` / `TEXT` | Ratio / Z-score | Involucre aspect ratio, standardized score, and availability status (`Observed` vs. `Vegetative_Missing`). |
| **Canonical GMM & Shrinkage CDA** | `selected_morph_features`, `gmm_cluster`, `gmm_uncertainty`, `cda_predicted_taxon`, `cda_posterior_prob`, `cda_oof_posterior_prob`, `cl_target_threshold`, `self_morph_mahal_d2`, `self_morph_df`, `cda_second_taxon`, `cda_second_posterior_prob`, `cda_posterior_margin`, `loocv_accuracy`, `can1`, `can2`, `can3` | `TEXT` / `INTEGER` / `REAL` | Semicolon list / Cluster ID / $[0, 1]$ / $D^2$ | Semicolon-delimited active morphometric predictors, canonical GMM cluster (`1..G`), in-sample and 10 km deme-blocked out-of-fold shrinkage CDA posteriors, Northcutt Confident Learning per-taxon self-confidence threshold, LOO within-taxon Mahalanobis tolerance distance and degrees of freedom, and canonical variates `can1..can3`. |
| **EPA Ecoregion CDA** | `ecoregion_cda_predicted`, `ecoregion_cda_posterior_prob`, `ecoregion_can1`, `ecoregion_can2` | `TEXT` / `REAL` | Ecoregion / $[0, 1]$ / Canonical SD | Supervised biogeographic discriminant axis scores testing regional morphological structuring across US EPA Level III Ecoregions. |
| **Type Control, OSR & Label Noise** | `has_empirical_type_anchor`, `type_morph_mahal_d2`, `nearest_type_morph_taxon`, `nearest_type_basionym`, `nearest_type_basionym_d2`, `type_morph_concordance`, `is_typical_morphotype`, `min_congener_mahal_d2`, `min_congener_taxon`, `open_set_novelty_p`, `open_set_df`, `morphometric_outlier_flag`, `misidentification_flag`, `is_confident_label_noise` | `INTEGER` / `REAL` / `TEXT` | `{0, 1}`, $D^2$, $\chi^2$ tail $p \in [0, 1]$ | Squared Mahalanobis distance to nomenclatural type morphotype ($c^*_k$) and basionym type centroids, Open-Set Recognition (`OSR`) minimum squared distance and exact feature dimensionality `open_set_df` across congeneric centroids, and Northcutt Confident Learning label-noise flag. |

### 3.4 `spatial_rf_predictions` (Step 05 Rhizosphere RF, Hopkins Phenology & Niche Control)

**Primary Key:** `(run_id, catalogNumber, plant_individual_id)` | **Populated By:** `05_multimodal_spatial_rf.R` and `step05_phenology_ecotypes.R`
**Companion CSV:** `runs/<RUN_ID>/data/tables/multimodal_conflict_flags.csv`

| Column Group | Column Names | SQLite Type | Units / Domain | Description |
| :--- | :--- | :---: | :--- | :--- |
| **Spatial Block-CV RF** | `rf_predicted_taxon`, `rf_predicted_prob`, `rf_given_taxon_prob`, `rf_prediction_source`, `edaphic_best_fit_taxon`, `spatial_block_cv_acc`, `spatial_autocorr_range_km`, `spatial_block_id`, `n_spatial_rf_predictors` | `TEXT` / `REAL` / `INTEGER` | Taxon / $[0, 1]$ / `km` / Block ID | Out-of-fold Spatial Block-CV Random Forest top taxon, max probability, probability assigned to recorded label, prediction provenance (`Spatial_Block_OOF_Anchor` vs. `Passive_Projection_*`), empirical variogram block range (`km`), and predictor count after 3-block VIF pruning. |
| **Selected Predictor Manifests** | `primary_env_predictors`, `primary_clim_predictors`, `primary_edaph_predictors` | `TEXT` | Semicolon-delimited names | Non-collinear (`VIF <= 5.0`) predictor subsets retained across all environmental, macro-climatic, and $0\text{ to }15\text{ cm}$ edaphic/topographic blocks; read dynamically by Step 06 Complex Discovery. |
| **Circular Hopkins Phenology** | `doy`, `expected_doy`, `hopkins_expected`, `delta_doy`, `hopkins_adjusted_doy`, `pheno_timing_category`, `anthesis_filter_status` | `REAL` / `TEXT` | Days (`d`) / Category | Observed DOY, empirical Hopkins latitude/elevation/BioClim expected DOY, signed isolation residual $\Delta\text{DOY}$ (`[-182.5, +182.5]` days via $\text{atan2}$ circular wrap), Hopkins-adjusted flowering day, and optical anthesis filter tag. |
| **Locus Classicus & Substrate** | `edaphic_z_score`, `climatic_z_score`, `type_edaphic_z_score`, `nearest_type_niche_taxon`, `matches_type_lithology`, `substrate_specialist_flag`, `edaphic_schoener_d` | `REAL` / `TEXT` / `INTEGER` | Robust RMS $Z$ / `{0, 1}` / $[0, 1]$ | Multivariate rhizosphere and climatic deviations from taxon reference anchors and *locus classicus* type localities, bedrock match indicator, outcrop/azonal specialist flag, and `PCA-env` Schoener's $D$ niche overlap. |
| **Multimodal Concordance** | `n_observed_modalities`, `multimodal_concordance`, `multimodal_conflict_flag`, `triage_category` | `INTEGER` / `REAL` / `TEXT` | `0..3` / $[0, 1]$ / `{0, 1}` | Count of observed evidence streams in Step 05, cross-stream agreement score, binary conflict flag, and preliminary Step 05 triage label. |

*(For convenience in standalone SQL queries, `spatial_rf_predictions` also mirrors the individual's extracted `bio1..bio19`, `fine_soil_*`, `ilr_texture_*`, `ssurgo_*`, `polaris_*`, `solus_*`, `topo_*`, and `sgmc_lithology_class` values alongside `environmental_cache`.)*

### 3.5 `curation_triage` (Step 06 Five-Stage Integrative Synthesis)

**Primary Key:** `(run_id, catalogNumber, plant_individual_id)` | **Populated By:** `scripts/analysis/06_triage_dashboard_synthesis.R`

| Column | SQLite Type | Units / Domain | Description |
| :--- | :---: | :--- | :--- |
| `run_id`, `catalogNumber`, `plant_individual_id` | `TEXT` / `INTEGER` | Composite PK | Run, specimen, and DBSCAN individual identifiers. |
| `w_morph_kappa`, `w_edaphic_kappa`, `w_pheno_kappa` | `REAL` | Reliability $\in [0.05, 0.95]$ | Empirical out-of-fold Cohen's $\kappa$ reliability weights computed on active `Tier_0_Type_Anchor`, `Tier_1_Gold`, and `Tier_2_Silver` anchors. |
| `n_available_modalities` | `INTEGER` | `1..3` | Number of non-missing evidence streams (Morphology, Rhizosphere/Climate, Verified-Anthesis Phenology) observed for this individual; weights renormalize strictly over observed streams. |
| `probabilistic_consensus_conflict`, `probabilistic_consensus_score` | `REAL` | Probability $[0, 1]$ | Precision-modulated cross-modal conflict probability $P_{\text{conflict},i}$ and concordance score $1 - P_{\text{conflict},i}$. |
| `hopkins_adjusted_doy`, `delta_doy_resid` | `REAL` | DOY / Days (`d`) | Circular Hopkins-standardized flowering day and signed phenological residual (days from conspecific expectation). |
| `taxonomic_status_call` | `TEXT` | Status call string | Diagnostic status call summarizing the active hierarchical decision stage. |
| `recommended_determination` | `TEXT` | *Packera* taxon / novel call | Recommended herbarium annotation label (reassigned *Packera* species, hybrid formula, `sp. nov. ined.`, or confirmed taxon). |
| `taxonomic_treatment_broad` | `TEXT` | Broad treatment enum | High-level curatorial routing bucket (`Confirmed_Taxon`, `Reassign_Determination`, `Interspecific_Hybrid`, `Novel_Lineage_Candidate`, `Infraspecific_or_Plastic`). |
| `taxonomic_treatment_subtype` | `TEXT` | Canonical 8-state enum | Five-stage biological determination: `Confirmed_Species`, `Confirmed_Type_Control`, `Type_Synonymy_Candidate`, `Misidentification_Reassign`, `Hybrid_Introgressant`, `Putative_Cryptic_Species`, `Subspecies_Ecotype`, or `Ecophenotypic_Plasticity`. |
| `cryptic_cluster_cohesion_pass` | `INTEGER` | `{0, 1}` | Equals `1` only when a non-core morphocluster ($g_{\text{cl}} > 1$) satisfies exsiccatae-deduplicated gathering cohesion ($\ge 5$ gatherings across $\ge 3$ 10 km spatial demes and $\ge 2$ collectors or herbaria), exact `open_set_df` congener rejection, multimodal barrier congruence, and `cryptic_non_clinal_pass == 1`. |
| `cryptic_non_clinal_pass` | `INTEGER` | `{0, 1}` | Equals `1` when a secondary morphocluster is maintained across parapatric contact ($\le 50\text{ km}$), Hopkins reproductive allochrony, or azonal lithology rather than a smooth allopatric geographic cline ($R^2 < 0.60$). |
| `resurrection_basionym` | `TEXT` | Basionym / `NULL` | Unstandardized heterotypic nomenclatural basionym matched within $\chi^2_{0.95, p}$ for Stage 3 `Taxonomic_Resurrection_Candidate` lineages (`NULL` for unnamed *de novo* cryptic candidates and all other stages). |
| `hybrid_zone_cohesion_pass` | `INTEGER` | `{0, 1}` | Equals `1` only when a multi-species specimen satisfies all five Conservative Contact-Zone Hybrid Swarm Gates (local two-parent syntopy within `15 km`, co-flowering within `21 d`, pairwise Mahalanobis chord intermediacy $t_{AB} \in [0.22, 0.78]$ with bounded $d^2_{\perp}$, multi-organ corroboration, and $\ge 2$ admixed vouchers in the local contact zone; always `0` in single-taxon `--intraspecific-mode`). |
| `hybrid_parent_1`, `hybrid_parent_2` | `TEXT` | *Packera* binomial / `NULL` | Recorded maternal/baseline parent (`given_sp`) and verified sympatric congener parent for confirmed `Hybrid_Introgressant` swarm vouchers (`NULL` for all non-hybrid calls and in single-taxon `--intraspecific-mode`). |
| `synthesis_triage_priority` | `TEXT` | `{'CRITICAL','HIGH','MEDIUM','LOW','RESOLVED'}` | Physical herbarium loan and annotation priority (`CRITICAL` for cohesive cryptic candidates and type synonymy conflicts; `HIGH` for confident misidentifications and hybrids; `RESOLVED` for concordant anchors). |
| `synthesis_triage_action` | `TEXT` | Curatorial action code | Concrete curatorial instruction (for example `Annotate & Re-folder`, `Inspect for Hybridity / Cytotype`, `Prioritize Sequenced Loan`, or `Retain Current Determination`). |
| `synthesis_rationale` | `TEXT` | Human-readable sentence | Quantitative biological explanation citing the individual's CDA posterior, OSR $D^2_{\min}$, rhizosphere $Z$-score, lithology, and $\Delta\text{DOY}$ residual. |
| `is_completed_voucher` | `INTEGER` | `{0, 1}` | Equals `1` when the individual completed integrative evaluation without upstream QC rejection. |

### 3.6 `taxon_complexes` (Step 06 Unsupervised Species-Complex Delimitation)

**Primary Key:** `(run_id, species_standardized)` | **Populated By:** `scripts/analysis/taxon_complex_discovery_utils.R`
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

### 4.1 `v_integrative_triage_master` (Composite Six-Stage SQL View)

Defined identically in `scripts/core/db_schema.py` and `scripts/core/db_utils.R`, `v_integrative_triage_master` joins `pipeline_runs`, `run_active_cohort`, `vouchers`, `plant_individuals`, leaf counts from `extracted_leaves`, `voucher_reproductive_phenology`, `environmental_cache` (via rounded `lat_lon_key`), `otu_morphometrics`, `spatial_rf_predictions`, `curation_triage`, and `taxon_complexes` on `(run_id, catalogNumber, plant_individual_id)`.

Two explicit precedence rules prevent stale values during partial re-runs (`--start-step 3` or `--start-step 5`):
1. `run_active_cohort` overrides raw `vouchers` columns for `species_standardized`, `determiner_tier`, `type_status`, and `spatial_balance_weight` via `COALESCE(rac.<col>, v.<col>)` so dynamic authority or synonym updates apply without re-harvesting images.
2. Step 05 `spatial_rf_predictions.hopkins_adjusted_doy` takes precedence over `curation_triage.hopkins_adjusted_doy` via `COALESCE(s.hopkins_adjusted_doy, t.hopkins_adjusted_doy)` so re-running Step 05 immediately propagates updated Hopkins phenology into Step 06.

### 4.2 `triage_queue.csv` and Interactive Dashboard Filtered CSV

* **`runs/<RUN_ID>/data/tables/triage_queue.csv`:** Full wide-table export written atomically by `build_triage_queue()` in `scripts/analysis/06_triage_dashboard_synthesis.R`, sorted by `synthesis_triage_priority` (`CRITICAL -> HIGH -> MEDIUM -> LOW -> RESOLVED`) and descending `coordinateUncertainty`. When `--focal-taxon` is supplied, SQLite `curation_triage` and `taxon_complexes` persist the complete multi-taxon reference cohort first, whereas `triage_queue.csv` filters to specimens originally filed under or reassigned into `--focal-taxon`.
* **`triage_queue_filtered.csv` (Browser Export from `Triage_Dashboard_Interactive.html`):** Exports the 17 curatorial bench columns (`catalogNumber`, `plant_individual_id`, `institutionCode`, `recordedBy`, `determiner_tier`, `species_standardized`, `recommended_determination`, `synthesis_triage_priority`, `taxonomic_treatment_subtype`, `cda_predicted_taxon`, `cda_posterior_prob`, `edaphic_best_fit_taxon`, `edaphic_z_score`, `delta_doy_resid`, `probabilistic_consensus_score`, `synthesis_triage_action`, `synthesis_rationale`) matching the active search, priority, treatment, and taxon filters in the offline HTML viewer.
