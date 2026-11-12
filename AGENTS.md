# Workspace Rules and Non-Negotiable Invariants (`Packera_dubia_morphometrics`)

**Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV5` | **Affiliation:** UNC-Chapel Hill Herbarium (NCU)

I define here the non-negotiable architectural and mathematical invariants alongside cross-language data contracts for this repository. Every contributor or agent that modifies Python (`.py`) or R (`.R`) code must enforce these rules without exception.

## 1. Non-Negotiable Code and Documentation Hygiene

* **Zero `_archive` Coupling:** Never import or source files from `scripts/_archive/`, `data/_archive/`, or `outputs/_archive/`. Active code must fail fast when required inputs are missing.
* **Zero Obsolete Framework Remnants:** The pipeline operates purely on LeafMachine3 (`lm3`) ONNX Runtime and modern Python and R statistical stacks without modifying files inside `LeafMachine3/`. Active files must contain zero obsolete segmentation framework remnants and zero retired helper tokens.
* **300-DPI PNG-Only Figure Architecture:** All 12 canonical publication plates (`fig02a` through `fig06b`) and native `LeafMachine3` visual caches (`lm3_assets/`, `ect_overlays/`, and `lm3_sheet_overlays/`) export strictly as 300-DPI `.png` rasters (`0` `.pdf` figure outputs) with unfilled hollow scatterplot symbols, combinatorial 60-taxon palette-and-shape scaling, and ICN rank-aware italicization.
* **Dynamic Genus-Wide Discovery:** The harvester queries genus-wide occurrences dynamically through `discover_genus_packera_species_gbif()`. Explicit `taxon` or `focal_taxon` runs restrict analysis strictly to the specified taxa and reconcile SQLite `vouchers` with `curated_vouchers.csv` via `reconcile_curated_vouchers_csv_with_db()`.
* **Positive Present-State Documentation:** Maintain all documentation strictly in terms of active capabilities and positive requirements alongside current implementations.
* **Zero Silent Synthetic Data:** Never generate unflagged random embeddings (`np.random.randn`) or synthetic raster values when `demo` mode is inactive. Synthetic data generation requires explicit `demo` mode (`PACKERA_DEMO_MODE=1`), which isolates writes to ephemeral SQLite (`demo_ephemeral.sqlite`) and leaves canonical stores untouched.
* **Authentic Voice and Documentation Standard (`@write-in-my-voice`):** All prose and SOP descriptions alongside `README.md`, `AGENTS.md`, `CONTRIBUTING.md`, `CHANGELOG.md`, and `docs/*.md` must follow Brandon Fuller's authentic voice (Writing DNA v5.2 in [`.agents/skills/write-in-my-voice/SKILL.md`](.agents/skills/write-in-my-voice/SKILL.md)), with zero em dashes, zero en dashes, zero tricolons, zero trailing participial clauses, and active first-person empirical framing.
* **Dual Virtualenv Boundaries:**
  * `.venv`: Primary orchestration (`main.py`, `run_pipeline.sh`), environmental layer downloads, and DuckDB/SQLite management.
  * `.venv_LM3` (`LeafMachine3/.venv`): LeafMachine3 ONNX Runtime isolation (YOLO26x, UNet++ v9, Ruler solver, ECT, and visual product rendering), relying strictly on standard-library `sqlite3`.
* **Verification Gate:** Verify code hygiene and cross-stage contracts via `.venv/bin/python scripts/tests/run_audit.py`.

## 2. Refactored Module Ownership Map

Parallel work must maintain file-disjoint ownership across the active directory structure:

| Directory | Active Source Files | Core Responsibility and Ownership Scope |
| :--- | :--- | :--- |
| **`scripts/core/`** | `config.py`, `db.py`, `db_schema.py`, `db_utils.R`, `harvester.py`, `taxonomy_authority.py`, `epa_ecoregions.py`, `usda_plants_ranges.py`, `artifact_manager.py`, `logger.py` | Configuration schemas, SQLite WAL engine (`20261121`), dynamic migrations, atomic CSV bridge, GBIF/iDigBio harvester, EPA Level III ecoregions, USDA PLANTS county ranges, and authority engine. |
| **`scripts/data_prep/`** | `01_voucher_harvester.py`, `download_environmental_layers.py` | Step 01 CLI parser and occurrence harvester; automated pre-Step-05 SOLUS100, USGS SGMC `.gpkg`, and cohort 3DEP 30 m DEM downloader. |
| **`scripts/vision/`** | `configure_leafmachine3.py`, `lm3_compat.py`, `lm3_data_loader.py`, `lm3_inference_engines.py`, `leaf_geometry_utils.py`, `organ_crown_metrics.py`, `capitulum_phenology_classifier.py`, `step02_figures.py` | LeafMachine3 config, ONNX Runtime engines, DBSCAN ramet clustering, anatomical `leaf_zone` classification, petiole/scape metrics, optical phenology, and Step 02 300-DPI PNG plates + native LM3 asset caches. |
| **`scripts/pipeline/`** | `01_voucher_harvester.py`, `02_segment_and_extract.py` | Step 01 pipeline alias; Step 02 read-only LeafMachine3 ingestion adapter, 3-stage wave engine, $32 \times 360$ ECT export, and per-voucher SQLite transactions. |
| **`scripts/morphometrics/`** | `03_fourier_extractor.R`, `04_gmm_morphotools.R`, `authority_weighting_utils.R`, `dissected_leaf_utils.R`, `fourier_parallel_utils.R`, `plasticity_mitigation_utils.R`, `pub_figures_utils.R`, `pub_ecoregion_manova_utils.R`, `pub_synthesis_plates.R` | Step 03 LM3 EFA QC (`<= 800` lines), homologized EFA, $32 \times 360$ ECT SVD, basal/cauline OTU medians; Step 04 (`<= 650` lines) deme plasticity mitigation, `RRPP` allometry, Stratified Authority Pool Cap, canonical GMM, shrinkage CDA, OSR, and EPA Ecoregion MANOVA. |
| **`scripts/figures/`** | `theme_packera_pub.R`, `step03_figures.R`, `step04_figures.R`, `step05_figures.R`, `step06_figures.R` | Modular `ggplot2` + `patchwork` + `sf` 300-DPI PNG figure suite (`fig03a` to `fig06b`), 60-taxon hollow symbol/palette scaling, ICN rank-aware formatting, and `EPSG:5070` projected mapping. |
| **`scripts/analysis/`** | `05_multimodal_spatial_rf.R`, `step05_phenology_ecotypes.R`, `step05_sdm_maxent.R`, `edaphic_fine_grain_utils.R`, `06_triage_dashboard_synthesis.R`, `taxon_complex_discovery_utils.R` | Step 05 19-BioClim and 0 to 15 cm rhizosphere Spatial RF, circular Hopkins phenology, and MaxEnt SDMs; Step 06 5-stage consensus triage, offline HTML dashboard, and Leiden/Ward.D2 complex discovery (`<= 800` lines). |
| **`scripts/tests/`** | `run_audit.py` + 15 `test_*.py` modules | 250-test unit, SQLite handshake, R statistical integration, figure contract, and static code-hygiene suite. |

## 3. Two-Tier Storage Architecture and Relational Contracts

The central SQLite WAL database (`data/packera_pipeline.sqlite`) pairs with synchronized atomic CSV snapshots:
* **Tier A: Canonical Shared Store (`data/`):** Houses invariant Step 01 and Step 02 specimen records, WebP scans, binary masks, 120-vertex tip-up contours, $32 \times 360$ ECT `.npy` matrices, organ crops, canonical manifests, and the environmental cache.
* **Tier B: Run-Scoped Stores (`runs/<RUN_ID>/` with `runs/latest`):** Isolates active cohort filtering and Step 03 to 06 statistical tables alongside reports, rasters, and 300-DPI PNG figures.
* **Segment and Run Reset (`purge-segments` / `purge-old-data`):** Invoking `purge_segments_and_runs()` in `scripts/core/db.py` clears Step 02 to 06 masks, contours, ECT matrices, organ crops, downstream CSVs, `runs/`, and SQLite Step 02 to 06 rows (`VACUUM`), whereas Step 01 `vouchers`, `taxon_synonyms`, `environmental_cache`, `data/raw_vouchers/`, and `data/tables/curated_vouchers.csv` remain intact.

### Database Tables and Master View (`data/packera_pipeline.sqlite`)

| Tier & Table / View | Primary Key | Key Invariants |
| :--- | :--- | :--- |
| **A: `taxon_synonyms`** | `raw_name_lower` | Maps raw names and *Senecio*/*Cineraria* basionyms to accepted *Packera* taxa; protologue authors. |
| **A: `vouchers`** | `catalogNumber` | Darwin Core metadata, coordinates, determiner tiers, USDA county status, duplicate conflicts, and scale calibration. |
| **A: `plant_individuals`** | `(catalogNumber, plant_individual_id)` | DBSCAN rosette centroids, caudex crowns, $P_{75}$ petiole length, and straight/bent scape height. |
| **A: `extracted_leaves`** | `(catalogNumber, plant_individual_id, leaf_id)` | Read-only LM3 leaf descriptors, `lm3_gates_pass`, `leaf_zone`, `contour_xy_json`, and `ect_path`. |
| **A: `voucher_reproductive_phenology`** | `(catalogNumber, plant_individual_id)` | Capitulum counts, involucre dimensions, and optical bud/anthesis/fruit state. |
| **A: `environmental_cache` / `rtree`** | `lat_lon_key` / `id` (R*Tree) | 0 to 15 cm rhizosphere SSURGO, POLARIS, SOLUS100, 3DEP, SGMC lithology, and WorldClim v2.1. |
| **B: `pipeline_runs`** | `run_id` | Execution timestamps, git commit SHA, CLI arguments, focal taxon, and run mode. |
| **B: `run_active_cohort`** | `(run_id, catalogNumber, plant_individual_id)` | Preflight active cohort: refreshed tiers, type status, and authority/spatial weights. |
| **B: `leaf_efa_harmonics`** | `(run_id, catalogNumber, plant_individual_id, leaf_id)` | Symmetrized harmonics ($A_1..D_{24}$), roughness, asymmetry, EFA PCs, and `ECT_PC1..PC3` / `ect_lobe_complexity`. |
| **B: `otu_morphometrics`** | `(run_id, catalogNumber, plant_individual_id)` | Pure basal OTU medians, 10 km spatial demes (`population_deme_id`), `RRPP` allometry, `scale_source`, `is_passive`, `AR_inv_std`, `reproductive_trait_status`, canonical GMM clusters, shrinkage CDA + Ecoregion CDA (`ecoregion_can1..2`), OSR metrics (`min_congener_mahal_d2`, `open_set_novelty_p`), and `is_confident_label_noise`. |
| **B: `spatial_rf_predictions`** | `(run_id, catalogNumber, plant_individual_id)` | Spatial Block-CV RF OOF predictions, `hopkins_adjusted_doy`, `delta_doy`, `type_edaphic_z_score`, and passive non-Tier-1 / high-uncertainty projections. |
| **B: `curation_triage` & `taxon_complexes`** | `(run_id, catalogNumber, plant_individual_id)` / `(run_id, species_standardized)` | OOF $\kappa$-weighted 5-stage treatments, `cryptic_cluster_cohesion_pass`, and Leiden/Ward.D2 empirical species complexes. |
| **B: `v_integrative_triage_master`** | SQL View | Pre-joined master view left-joining all modalities onto `plant_individuals` (`PACKERA_SCHEMA_USER_VERSION = 20261121`), giving upstream Step 05 `hopkins_adjusted_doy` precedence over stale triage rows. |

### Mandatory Cross-Language Schema and Handoff Invariants

* **Composite Grouping Key `(catalogNumber, plant_individual_id)`:** All joins across Python and R must group by this composite key and never collapse multi-plant sheets to `plant_individual_id = 0`.
* **Schema Synchronization (`PACKERA_SCHEMA_USER_VERSION = 20261121`):** Python (`db_schema.py`) and R (`db_utils.R`) must maintain identical table schemas and `v_integrative_triage_master` definitions.
* **Read-Only LeafMachine3 Step 02 Adapter and Step 03 R EFA QC:** Step 02 persists 100% of `LeafMachine3` segmented leaves, native symmetry/QC metrics, 120-vertex tip-up `lamina_base`-anchored contours, and $32 \times 360$ ECT matrices to SQLite without destructive Python-side leaf dropping; Step 03 (`evaluate_step03_lm3_efa_qc`) filters unresolved clumps in R before running `coo_center() %>% coo_scale()` and `Momocs::efourier(..., norm = FALSE, start = FALSE)` (preserving reniform/obovate $L/W < 1.0$ homology without `norm = TRUE` ellipse flipping).
* **Stratified Authority Pool Cap and Type Control:** Active anchors enforce `Tier_0_Type = 3.0`, `Tier_1_Gold = 2.0`, `Tier_2_Silver = 0.75` (capped at $\le 40\%$ of active pool weight when $N_{\text{type+gold}} \ge 5$), `Promoted_Bronze = 0.50`, and unpromoted `Tier_3_Bronze = 0.0`. Nomenclatural types retain MCD immunity (`Tier_0_Type_Anchor`), anchor per-species typical morphotypes ($c^*_k$), and follow Rules A, B, and C alongside the Stage 3 Population-Cohesion Lineage Gate (`cryptic_cluster_cohesion_pass`).

## 4. Documentation Pointers

* **Operational Manual and CLI Reference:** [`README.md`](README.md)
* **Mathematical SOP and Botanical Protocol:** [`docs/TAXONOMIC_SOP.md`](docs/TAXONOMIC_SOP.md)
* **Relational Schema and CSV Data Dictionary:** [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md)
* **Environmental and Pedological Data Provenance:** [`docs/ENVIRONMENTAL_DATA_SOURCES.md`](docs/ENVIRONMENTAL_DATA_SOURCES.md)
* **Developer and QA Contribution Guide:** [`CONTRIBUTING.md`](CONTRIBUTING.md)
* **Version and Schema Release History:** [`CHANGELOG.md`](CHANGELOG.md)
* **Citation Metadata and Software License:** [`CITATION.cff`](CITATION.cff) and [`LICENSE`](LICENSE)
* **Automated Contract and Hygiene Audit Skill:** [`.agents/skills/pipeline-contract-audit/SKILL.md`](.agents/skills/pipeline-contract-audit/SKILL.md)
* **Writing Voice Skill:** [`.agents/skills/write-in-my-voice/SKILL.md`](.agents/skills/write-in-my-voice/SKILL.md)
