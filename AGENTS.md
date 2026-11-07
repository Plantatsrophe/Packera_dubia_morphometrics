# Workspace Rules & Non-Negotiable Invariants (`Packera_dubia_morphometrics`)

**Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV5` | **Affiliation:** UNC-Chapel Hill Herbarium (NCU)

Defines non-negotiable architectural, mathematical, and cross-language data contracts. Every agent modifying Python (`.py`) or R (`.R`) code must enforce these rules without exception.

---

## 1. Non-Negotiable Code Hygiene

* **Zero `_archive` Coupling:** Never import, source, read, or fall back to `scripts/_archive/`, `data/_archive/`, or `outputs/_archive/`. Active code must fail fast if inputs are missing.
* **Zero Legacy Framework Remnants:** Operates purely on LeafMachine3 (`lm3`) ONNX Runtime and modern PyTorch/R stacks. Active files must contain zero legacy ML/segmentation framework remnants or obsolete helper tokens.
* **Zero Hardcoded Species Lists:** Dynamic genus-wide discovery (`discover_genus_packera_species_gbif()`, `target_species: []`). CLI `--taxon` or `--focal-taxon` runs restrict analysis strictly to specified taxa with zero hardcoded defaults or injected outgroups.
* **Zero Silent Synthetic Data:** Never generate unflagged random embeddings (`np.random.randn`), static Gaussian heatmaps, or trigonometric environmental values. Synthetic generation requires explicit `--demo` (`PACKERA_DEMO_MODE=1`), which isolates writes to ephemeral SQLite (`demo_ephemeral.sqlite`) and never touches canonical stores.
* **Dual Virtualenv Boundaries:**
  * `.venv`: Primary orchestration (`main.py`, `run_pipeline.sh`), DINOv2 (`dinov2_vitb14_reg`), Cleanlab, DuckDB, general scientific tools.
  * `.venv_LM3`: LeafMachine3 ONNX Runtime isolation (YOLO26x, UNet++ v9, Ruler solver, ECT), relying strictly on standard-library `sqlite3` without compiled C/C++ database drivers.
* **Verification Gate:** Verify code hygiene via `pytest scripts/tests/test_code_hygiene.py -q`.

---

## 2. Refactored Module Ownership Map

Subagents dispatched for parallel work must maintain file-disjoint ownership across the refactored directory structure:

| Directory | Active Source Files | Core Responsibility & Ownership Scope |
| :--- | :--- | :--- |
| **`scripts/core/`** | `config.py`, `db.py`, `db_schema.py`, `harvester.py`, `epa_ecoregions.py`, `taxonomy_authority.py`, `artifact_manager.py`, `logger.py`, `db_utils.R` | Configuration schemas, SQLite WAL engine, dynamic migrations, atomic CSV bridge, GBIF/iDigBio harvester, EPA L3 ecoregions, authority engine. |
| **`scripts/vision/`** | `configure_leafmachine3.py`, `lm3_data_loader.py`, `lm3_inference_engines.py`, `leaf_geometry_utils.py`, `leaf_gatekeeper.py`, `organ_crown_metrics.py`, `capitulum_phenology_classifier.py`, `step02_figures.py` | LeafMachine3 config, ONNX Runtime engines, DBSCAN plant clustering, classical gatekeeper, petiole/scape metrics, optical phenology. |
| **`scripts/pipeline/`** | `01_voucher_harvester.py`, `02_segment_and_extract.py` | Step 01 occurrence harvesting/curation; Step 02 LeafMachine3 segmentation, boundary filtering, and per-voucher atomic SQLite transactions. |
| **`scripts/morphometrics/`** | `03_fourier_extractor.R`, `04_gmm_morphotools.R`, `authority_weighting_utils.R`, `dissected_leaf_utils.R`, `fourier_parallel_utils.R`, `pub_figures_utils.R`, `pub_ecoregion_manova_utils.R`, `pub_synthesis_plates.R` | Step 03 strict Tier 1 EFA, symmetrization, OTU medians; Step 04 allometry, Dawid-Skene weighting, GMM morphotypes, Track 1/2 CDA + MANOVA. |
| **`scripts/analysis/`** | `05_cleanlab_vision_xai.py`, `06_multimodal_spatial_rf.R`, `step06_phenology_ecotypes.R`, `step06_sdm_maxent.R`, `07_triage_dashboard_synthesis.R`, `08_discover_taxon_complexes.R`, `edaphic_fine_grain_utils.R`, `vision_dissection_utils.py` | Step 05 DINOv2 patch XAI & Cleanlab; Step 06 19-BioClim/0–15cm Spatial RF & SDMs; Step 07 consensus triage; Step 08 unsupervised complex discovery. |
| **`scripts/tests/`** | `test_code_hygiene.py`, `test_sql_pipeline_handshake.py`, `test_voucher_harvester.py`, `test_segment_and_extract.py`, `test_fourier_extractor.py`, `test_gmm_morphotools.py`, `test_cleanlab_vision_xai.py`, `test_multimodal_spatial_rf.py`, `test_triage_dashboard_synthesis.py`, `test_discover_taxon_complexes.py`, `run_audit.py` | Unit, integration, and cross-stage contract test suites (`pytest scripts/tests/ -q`). |

---

## 3. Two-Tier Storage Architecture & Relational Contracts

Central SQLite WAL database (`data/packera_pipeline.sqlite`) paired with synchronized atomic CSV snapshots:
* **Tier A — Canonical Shared Store (`data/`):** Invariant Step 01–02 assets, masks, 120-vertex contours, crops, environmental cache, and canonical tables.
* **Tier B — Run-Scoped Stores (`runs/<RUN_ID>/` with `runs/latest`):** Active cohort filtering, Step 03–08 statistics, reports, figures, and run logs.

### Database Tables & Master View (`data/packera_pipeline.sqlite`)
| Tier & Table / View | Primary Key | Key Invariants |
| :--- | :--- | :--- |
| **A: `taxon_synonyms`** | `raw_name_lower` | Standardizes raw names to accepted *Packera* taxa; protologue authors. |
| **A: `vouchers`** | `catalogNumber` | Darwin Core metadata, coordinates, determiner tiers, scale calibration. |
| **A: `plant_individuals`** | `(catalogNumber, plant_individual_id)` | Rosette centroids, caudex crowns, petiole/height metrics, bent stems. |
| **A: `extracted_leaves`** | `(catalogNumber, plant_individual_id, leaf_id)` | Morphological descriptors, `assigned_tier`, `contour_xy_json`. |
| **A: `voucher_reproductive_phenology`** | `(catalogNumber, plant_individual_id)` | Capitulum counts, involucre dimensions, optical anthesis state. |
| **A: `dinov2_feature_cache`** | `(catalogNumber, plant_individual_id, model_name, sigma)` | Caches 768-d embeddings (`embedding_json`) and foreground attention ratio. |
| **A: `environmental_cache` / `rtree`** | `lat_lon_key` / `id` (R*Tree) | 0–15cm rhizosphere SSURGO, POLARIS, SOLUS100, 3DEP, SGMC lithology, WorldClim. |
| **B: `pipeline_runs`** | `run_id` | Execution timestamps, git commit SHA, CLI arguments, focal taxon, mode. |
| **B: `run_active_cohort`** | `(run_id, catalogNumber, plant_individual_id)` | Preflight active cohort: tiers, type status, authority and spatial weights. |
| **B: `leaf_efa_harmonics`** | `(run_id, catalogNumber, plant_individual_id, leaf_id)` | Symmetrized harmonics ($A_1..D_{12}$), roughness, asymmetry variance, PCs. |
| **B: `otu_morphometrics`** | `(run_id, catalogNumber, plant_individual_id)` | Specimen OTU medians, Klingenberg allometry, GMM clusters, CDA scores, $D^2_{\text{type}}$. |
| **B: `vision_xai_audit`** | `(run_id, catalogNumber, plant_individual_id)` | Cleanlab OOF predictions, confidence, margins, OOD score, type cosines. |
| **B: `spatial_rf_predictions`** | `(run_id, catalogNumber, plant_individual_id)` | Spatial RF OOF predictions, block CV accuracy, Hopkins residuals, $z_{\text{type}}$. |
| **B: `curation_triage`** | `(run_id, catalogNumber, plant_individual_id)` | Consensus triage calls, OOF kappa weights, 5-stage treatments, priority. |
| **B: `v_integrative_triage_master`** | SQL View | Pre-joined master relational view left-joining all modalities onto `plant_individuals`. |

### Mandatory Cross-Language Schema & Handoff Invariants
* **Composite Grouping Key `(catalogNumber, plant_individual_id)`:** All joins across Python and R must group by this composite key. Never collapse multi-plant sheets to `plant_individual_id = 0` or group by `catalogNumber` alone.
* **Routing Tier `assigned_tier`:** `'Tier_1_Direct'` (intact silhouettes for EFA) vs `'Tier_2_Reflected'` (folded blades for linear traits). In R, match case-insensitively: `grepl("^tier_?[12]", assigned_tier, ignore.case = TRUE)`.
* **Scale & Geometry:** Empirical `pixels_per_mm` and boolean `is_calibrated_mm`. Adaptive blade floor $\min(0.75, \max(0.15, 0.25 \times \text{median\_area}))$. Petiole $L_{\text{petiole}}$ from crown $p_{\text{crown}}$ at $P_{75}$. Scape height $H_{\text{plant}}$ accounts for bent stems. Rosette background masked to `[242, 239, 230]`.
* **Nomenclatural Types & Authority:** Ranks 2 & 1 exempt from uncertainty/date drop filters, retained in spatial thinning, prioritized in outgroups. Supervised anchors (`Tier_1_Gold`, `Tier_2_Silver`, promoted Bronze within $\chi^2_{0.80}$) fit models; unpromoted Bronze receive $0.0$. Type retention immunity (`Tier_0_Type_Anchor`); extreme types ($> \chi^2_{0.99}$) isolated from FastMCD. Centers shrink toward type medoids ($\gamma_k$).
* **Strict Tier-1 EFA (Step 03):** Batch SQL from `contour_xy_json`. `coo_alignxax`, dynamic $H_{\text{opt}} \in [6, 24]$ ($\ge 99\%$ power), `rm_asym`, $k_{95}$ variance PCs. Type contours exempt from MCD filtering. Unweighted OTU medians.
* **Morphometrics & Allometry (Step 04):** Klingenberg slope homogeneity (`lm.rrpp`). Horn/Williams-Titus subspace. Shrinkage CDA with honest LOOCV. Exports $D^2_{\text{type}, i}$, `nearest_type_morph_taxon`, `cda_second_taxon`, and cluster $c^*$.
* **Vision XAI & Cleanlab (Step 05):** Indumentum-blind filter ($\sigma = 1.0\text{ px}$) on paper `[242, 239, 230]`. Foreground $16 \times 16$ patch pooling (`dinov2_vitb14_reg`). Types receive weight 4.5 and spherical prototype $\mathbf{v}_{\text{type}, c}$. Cleanlab bounded weights `[0.35, 3.0]`; primary types immune to corruption flags. Flag `foreground_attention_ratio < 0.65` as artifact.
* **Spatial RF & *Locus Classicus* (Step 06):** 19 BioClim + 0–15cm edaphics. Uncertainty-gated type extraction ($\le 2000\text{ m}$ edaphic, $\le 10000\text{ m}$ climatic) with proximate Gold fallback ($\le 150\text{ km}$). Predictor gating (VIF $< 5.0$, $|r| < 0.70$, RF importance top 85%). Spatial Block CV ($0.5^\circ$) with weights `[0.25, 4.0]`. MaxEnt + Spatial RF.
* **Consensus Triage & Complexes (Steps 07–08):** Reliability weighting via out-of-fold $\kappa_m$. Circular anthesis DOY windows. Rule A (type immunity / `Confirmed_Type_Control`), Rule B (type-control reassignment on $D^2_{\text{type}} > \chi^2_{0.95}$, $z_{\text{type}} > 2.0$), Rule C (*sensu stricto* locking for $c^*$). Unsupervised graph modularity (Leiden/Louvain) and Ward.D2 clustering across North American *Packera*.

---

## 4. Documentation Pointers & External References

* **CLI Execution & Operational Manual:** See [`README.md`](README.md) for installation, WSL 2 setup, environment configuration, quick-start commands, and CLI flag tables.
* **Mathematical SOP & Scientific Methodology:** See [`docs/TAXONOMIC_SOP.md`](docs/TAXONOMIC_SOP.md) for complete mathematical formulas, biological justifications, and literature citations.
* **Contract & Quality Verification:** See [`.agents/skills/pipeline-contract-audit/SKILL.md`](.agents/skills/pipeline-contract-audit/SKILL.md) for automated cross-stage contract audits and LeafMachine3 computer vision guardrails.
