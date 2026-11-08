# Workspace Rules and Non-Negotiable Invariants (`Packera_dubia_morphometrics`)

**Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV5` | **Affiliation:** UNC-Chapel Hill Herbarium (NCU)

Defines non-negotiable architectural, mathematical, and cross-language data contracts. Every agent modifying Python (`.py`) or R (`.R`) code must enforce these rules without exception.

---

## 1. Non-Negotiable Code and Documentation Hygiene

* **Zero `_archive` Coupling:** Never import, source, read, or fall back to `scripts/_archive/`, `data/_archive/`, or `outputs/_archive/`. Active code must fail fast if inputs are missing.
* **Zero Obsolete Framework Remnants:** Operates purely on LeafMachine3 (`lm3`) ONNX Runtime and modern PyTorch/R stacks. Active files must contain zero obsolete ML or segmentation framework remnants and zero obsolete helper tokens.
* **Dynamic Genus-Wide Discovery:** Queries genus-wide occurrences dynamically through `discover_genus_packera_species_gbif()`. CLI `--taxon` or `--focal-taxon` runs restrict analysis strictly to specified taxa.
* **Positive Present-State Documentation:** Maintain all documentation strictly in terms of active capabilities, positive requirements, and current implementations. Never define system behavior by referencing removed subsystems or documenting negative omissions.
* **Zero Silent Synthetic Data:** Never generate unflagged random embeddings (`np.random.randn`), static Gaussian heatmaps, or trigonometric environmental values. Synthetic generation requires explicit `--demo` (`PACKERA_DEMO_MODE=1`), which isolates writes to ephemeral SQLite (`demo_ephemeral.sqlite`) and never touches canonical stores.
* **Authentic Voice and Documentation Standard (`@write-in-my-voice`):** All prose, documentation, SOP descriptions, and README sections must follow Brandon Fuller's authentic voice (Writing DNA v5.2 in `.agents/skills/write-in-my-voice/SKILL.md`). Enforce zero em dashes, zero en dashes, zero banned AI buzzwords, and active first-person empirical framing.
* **Dual Virtualenv Boundaries:**
  * `.venv`: Primary orchestration (`main.py`, `run_pipeline.sh`), DuckDB, general scientific tools.
  * `.venv_LM3`: LeafMachine3 ONNX Runtime isolation (YOLO26x, UNet++ v9, Ruler solver, ECT), relying strictly on standard-library `sqlite3` without compiled C/C++ database drivers.
* **Verification Gate:** Verify code hygiene via `pytest scripts/tests/test_code_hygiene.py -q`.

---

## 2. Refactored Module Ownership Map

Subagents dispatched for parallel work must maintain file-disjoint ownership across the refactored directory structure:

| Directory | Active Source Files | Core Responsibility and Ownership Scope |
| :--- | :--- | :--- |
| **`scripts/core/`** | `config.py`, `db.py`, `db_schema.py`, `harvester.py`, `epa_ecoregions.py`, `taxonomy_authority.py`, `artifact_manager.py`, `logger.py`, `db_utils.R` | Configuration schemas, SQLite WAL engine, dynamic migrations, atomic CSV bridge, GBIF/iDigBio harvester, EPA L3 ecoregions, authority engine. |
| **`scripts/vision/`** | `configure_leafmachine3.py`, `lm3_data_loader.py`, `lm3_inference_engines.py`, `leaf_geometry_utils.py`, `leaf_gatekeeper.py`, `organ_crown_metrics.py`, `capitulum_phenology_classifier.py`, `step02_figures.py` | LeafMachine3 config, ONNX Runtime engines, DBSCAN plant clustering, classical gatekeeper, petiole/scape metrics, optical phenology. |
| **`scripts/pipeline/`** | `01_voucher_harvester.py`, `02_segment_and_extract.py` | Step 01 occurrence harvesting/curation; Step 02 LeafMachine3 segmentation, boundary filtering, and per-voucher atomic SQLite transactions. |
| **`scripts/morphometrics/`** | `03_fourier_extractor.R`, `04_gmm_morphotools.R`, `authority_weighting_utils.R`, `dissected_leaf_utils.R`, `fourier_parallel_utils.R`, `pub_figures_utils.R`, `pub_ecoregion_manova_utils.R`, `pub_synthesis_plates.R` | Step 03 strict Tier 1 EFA, symmetrization, OTU medians; Step 04 allometry, Dawid-Skene weighting, GMM morphotypes, Track 1/2 CDA + MANOVA. |
| **`scripts/analysis/`** | `05_multimodal_spatial_rf.R`, `step05_phenology_ecotypes.R`, `step05_sdm_maxent.R`, `06_triage_dashboard_synthesis.R`, `07_discover_taxon_complexes.R`, `edaphic_fine_grain_utils.R` | Step 05 19-BioClim/0 to 15 cm Spatial RF and SDMs; Step 06 consensus triage; Step 07 unsupervised complex discovery. |
| **`scripts/tests/`** | `test_code_hygiene.py`, `test_sql_pipeline_handshake.py`, `test_voucher_harvester.py`, `test_segment_and_extract.py`, `test_fourier_extractor.py`, `test_gmm_morphotools.py`, `test_multimodal_spatial_rf.py`, `test_triage_dashboard_synthesis.py`, `test_discover_taxon_complexes.py`, `run_audit.py` | Unit, integration, and cross-stage contract test suites (`pytest scripts/tests/ -q`). |

---

## 3. Two-Tier Storage Architecture and Relational Contracts

Central SQLite WAL database (`data/packera_pipeline.sqlite`) paired with synchronized atomic CSV snapshots:
* **Tier A: Canonical Shared Store (`data/`):** Invariant Step 01 to 02 assets, masks, 120-vertex contours, crops, environmental cache, and canonical tables.
* **Tier B: Run-Scoped Stores (`runs/<RUN_ID>/` with `runs/latest`):** Active cohort filtering, Step 03 to 07 statistics, reports, figures, and run logs.

### Database Tables and Master View (`data/packera_pipeline.sqlite`)
| Tier & Table / View | Primary Key | Key Invariants |
| :--- | :--- | :--- |
| **A: `taxon_synonyms`** | `raw_name_lower` | Standardizes raw names to accepted *Packera* taxa; protologue authors. |
| **A: `vouchers`** | `catalogNumber` | Darwin Core metadata, coordinates, determiner tiers, scale calibration. |
| **A: `plant_individuals`** | `(catalogNumber, plant_individual_id)` | Rosette centroids, caudex crowns, petiole/height metrics, bent stems. |
| **A: `extracted_leaves`** | `(catalogNumber, plant_individual_id, leaf_id)` | Morphological descriptors, `assigned_tier`, `contour_xy_json`. |
| **A: `voucher_reproductive_phenology`** | `(catalogNumber, plant_individual_id)` | Capitulum counts, involucre dimensions, optical anthesis state. |
| **A: `environmental_cache` / `rtree`** | `lat_lon_key` / `id` (R*Tree) | 0 to 15 cm rhizosphere SSURGO, POLARIS, SOLUS100, 3DEP, SGMC lithology, WorldClim. |
| **B: `pipeline_runs`** | `run_id` | Execution timestamps, git commit SHA, CLI arguments, focal taxon, mode. |
| **B: `run_active_cohort`** | `(run_id, catalogNumber, plant_individual_id)` | Preflight active cohort: tiers, type status, authority and spatial weights. |
| **B: `leaf_efa_harmonics`** | `(run_id, catalogNumber, plant_individual_id, leaf_id)` | Symmetrized harmonics ($A_1..D_{12}$), roughness, asymmetry variance, PCs. |
| **B: `otu_morphometrics`** | `(run_id, catalogNumber, plant_individual_id)` | Specimen OTU medians, Klingenberg allometry, GMM clusters, CDA scores, $D^2_{\text{type}}$. |
| **B: `spatial_rf_predictions`** | `(run_id, catalogNumber, plant_individual_id)` | Spatial RF OOF predictions, block CV accuracy, Hopkins residuals, $z_{\text{type}}$. |
| **B: `curation_triage`** | `(run_id, catalogNumber, plant_individual_id)` | Consensus triage calls, OOF kappa weights, 5-stage treatments, priority. |
| **B: `v_integrative_triage_master`** | SQL View | Pre-joined master relational view left-joining all modalities onto `plant_individuals`. |

### Mandatory Cross-Language Schema and Handoff Invariants
* **Composite Grouping Key `(catalogNumber, plant_individual_id)`:** All joins across Python and R must group by this composite key. Never collapse multi-plant sheets to `plant_individual_id = 0` or group by `catalogNumber` alone.
* **Routing Tier `assigned_tier`:** `'Tier_1_Direct'` (intact silhouettes for EFA) vs `'Tier_2_Reflected'` (folded blades for linear traits). In R, match case-insensitively: `grepl("^tier_?[12]", assigned_tier, ignore.case = TRUE)`.
* **Scale and Geometry:** Empirical `pixels_per_mm` and boolean `is_calibrated_mm`. Native LeafMachine3 stride-32 letterboxing (`auto=True, stride=32`), `overlap = 0.95` Stage 1 nested-box deduplication, tight `crop_pad = 0.0` Stage 2 leaf crops, and $32 \times 32$ ECT matrices. Per-plant (`plant_individual_id`) adaptive blade floor $\min(0.75, \max(0.15, 0.25 \times \text{median\_area}))$ and relative area floor ($0.15 \times \max\text{area}$) without a sheet-wide leaf cap. Petiole $L_{\text{petiole}}$ from crown $p_{\text{crown}}$ at $P_{75}$. Scape height $H_{\text{plant}}$ accounts for bent stems. Rosette background masked to `[242, 239, 230]`.
* **Nomenclatural Types and Authority:** Ranks 2 and 1 exempt from uncertainty/date drop filters, retained in spatial thinning, prioritized in outgroups. Supervised anchors (`Tier_1_Gold`, `Tier_2_Silver`, promoted Bronze within $\chi^2_{0.80}$) fit models; unpromoted Bronze receive $0.0$. Type retention immunity (`Tier_0_Type_Anchor`); extreme types ($> \chi^2_{0.99}$) isolated from FastMCD. Centers shrink toward type medoids ($\gamma_k$).
* **Strict Tier-1 EFA (Step 03):** Batch SQL from `contour_xy_json`. `coo_alignxax`, dynamic $H_{\text{opt}} \in [6, 24]$ ($\ge 99\%$ power), `rm_asym`, $k_{95}$ variance PCs. Type contours exempt from MCD filtering. Unweighted OTU medians.
* **Morphometrics and Allometry (Step 04):** Klingenberg slope homogeneity (`lm.rrpp`). Horn/Williams-Titus subspace. Shrinkage CDA with honest LOOCV. Exports $D^2_{\text{type}, i}$, `nearest_type_morph_taxon`, `cda_second_taxon`, and cluster $c^*$.
* **Spatial RF and *Locus Classicus* (Step 05):** 19 BioClim + 0 to 15 cm edaphics. Uncertainty-gated type extraction ($\le 2000\text{ m}$ edaphic, $\le 10000\text{ m}$ climatic) with proximate Gold fallback ($\le 150\text{ km}$). Predictor gating (VIF $< 5.0$, $|r| < 0.70$, RF importance top 85%). Spatial Block CV ($0.5^\circ$) with weights `[0.25, 4.0]`. MaxEnt + Spatial RF.
* **Consensus Triage and Complexes (Steps 06 to 07):** Reliability weighting via out-of-fold $\kappa_m$. Circular anthesis DOY windows. Rule A (type immunity / `Confirmed_Type_Control`), Rule B (type-control reassignment on $D^2_{\text{type}} > \chi^2_{0.95}$, $z_{\text{type}} > 2.0$), Rule C (*sensu stricto* locking for $c^*$). Unsupervised graph modularity (Leiden/Louvain) and Ward.D2 clustering across North American *Packera*.

---

## 4. Documentation Pointers and External References

* **CLI Execution and Operational Manual:** See [`README.md`](README.md) for installation, WSL 2 setup, environment configuration, quick-start commands, and CLI flag tables.
* **Mathematical SOP and Scientific Methodology:** See [`docs/TAXONOMIC_SOP.md`](docs/TAXONOMIC_SOP.md) for complete mathematical formulas, biological justifications, and literature citations.
* **Contract and Quality Verification:** See [`.agents/skills/pipeline-contract-audit/SKILL.md`](.agents/skills/pipeline-contract-audit/SKILL.md) for automated cross-stage contract audits and LeafMachine3 computer vision guardrails.
* **Writing Voice Skill:** See [`.agents/skills/write-in-my-voice/SKILL.md`](.agents/skills/write-in-my-voice/SKILL.md) for Writing DNA v5.2 guidelines and calibration rules.
