# Automated High-Throughput Vegetative Phenotyping and Herbarium Triage Pipeline for North American *Packera*

### Operationalizing Morphological Diagnosability under the Unified Species Concept

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/) [![R 4.3+](https://img.shields.io/badge/R-4.3+-blue.svg)](https://www.r-project.org/) [![LeafMachine3](https://img.shields.io/badge/LeafMachine3-ONNX_Runtime-brightgreen.svg)](https://github.com/Gene-Weaver/LeafMachine3)

- **Principal Investigator:** J. Brandon Fuller (PhD Candidate, Department of Biology, UNC-CH)
- **Faculty Advisor:** Dr. Alan S. Weakley (Director, NCU Herbarium; UNC Biology)
- **Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV5` (Detailed theory, equations, and proofs in [`docs/TAXONOMIC_SOP.md`](docs/TAXONOMIC_SOP.md))
- **Focal Dissertation Taxon:** ***Packera dubia*** (Spreng.) Trock & Mabb. (Asteraceae: Senecioneae) evaluated dynamically against all regional and continent-wide North American congeners.

---

## 1. Project Overview & 8-Step Architecture

This repository contains the production computational pipeline for automated vegetative phenotyping, geometric morphometrics, deep vision XAI, and integrative taxonomic triage of North American *Packera*. The system operates under the **Unified Species Concept (USC; de Queiroz 2007)** via **Integration by Congruence (Padial et al. 2010)**, evaluating vegetative morphometrics alongside capitulum phenology, nextRAD phylogenomics, and multi-scale $0\text{–}15\text{ cm}$ rhizosphere edaphics.

Data transactions synchronize across a centralized transactional SQLite database (`data/packera_pipeline.sqlite`) in WAL mode and atomic CSV snapshots under a Canonical Two-Tier Storage Architecture. For complete mathematical formulations, statistical proofs, and biological justifications, refer to [`docs/TAXONOMIC_SOP.md`](docs/TAXONOMIC_SOP.md).

### 8-Step Pipeline Architecture

| Step / Phase | Entrypoint Script | Helper Modules | Key Inputs | Key Outputs | Core Purpose |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **01: Ingestion** | `scripts/pipeline/01_voucher_harvester.py` | `scripts/core/` (`harvester.py`, `epa_ecoregions.py`, `taxonomy_authority.py`, `db.py`) | GBIF/iDigBio, DwC queries, EPA L3 ecoregions | `vouchers`, `taxon_synonyms`, `curated_vouchers.csv` | Ingests occurrences; expands synonymy; enforces type hierarchy; 3-layer authority; EPA L3 ecoregions. |
| **02: Extraction** | `scripts/pipeline/02_segment_and_extract.py` | `scripts/vision/` (`configure_leafmachine3.py`, `lm3_data_loader.py`, `lm3_inference_engines.py`, `leaf_gatekeeper.py`, `organ_crown_metrics.py`) | `vouchers`, images, LeafMachine3 ONNX models | `plant_individuals`, `extracted_leaves`, masks, crops | LeafMachine3 ONNX (YOLO26x, UNet++ v9, Ruler, 31-kp pose); botanical gatekeeper; petiole/height metrics. |
| **03: Fourier EFA** | `scripts/morphometrics/03_fourier_extractor.R` | `scripts/morphometrics/` (`authority_weighting_utils.R`, `fourier_parallel_utils.R`), `scripts/core/db_utils.R` | `extracted_leaves.contour_xy_json`, `vouchers` | `leaf_efa_harmonics`, `leaf_efa_harmonics.csv` | Batch SQL contour ingestion; inertial alignment; harmonic power calibration ($H_{\text{opt}} \ge 99\%$); symmetrization; OTU medians. |
| **04: Morphometrics** | `scripts/morphometrics/04_gmm_morphotools.R` | `scripts/morphometrics/` (`authority_weighting_utils.R`, `pub_figures_utils.R`), `scripts/core/db_utils.R` | `leaf_efa_harmonics`, `vouchers`, phenology | `otu_morphometrics`, `morphometric_flags.csv` | Klingenberg allometry; Horn/Williams-Titus subspace; Dawid-Skene authority; Track 1 GMM morphotypes; Track 2 EPA CDA. |
| **05: Vision XAI** | `scripts/analysis/05_cleanlab_vision_xai.py` | `scripts/analysis/vision_dissection_utils.py`, `scripts/core/db.py` | `vouchers`, `plant_individuals`, `cropped_patches/` | `dinov2_feature_cache`, `vision_xai_audit` | Indumentum-blind filter; foreground $16 \times 16$ patch pooling (`dinov2_vitb14_reg`); spherical prototypes; Cleanlab OOD. |
| **06: Spatial RF** | `scripts/analysis/06_multimodal_spatial_rf.R` | `scripts/analysis/` (`step06_phenology_ecotypes.R`, `step06_sdm_maxent.R`, `edaphic_fine_grain_utils.R`) | `otu_morphometrics`, `vision_xai_audit`, WorldClim, SSURGO | `spatial_rf_predictions`, `multimodal_conflict_flags.csv` | 19 BioClim + 0–15cm edaphics; *locus classicus* baseline ($z_{\text{type}}$); predictor gating; spatial block CV; SDMs. |
| **07: Triage** | `scripts/analysis/07_triage_dashboard_synthesis.R` | `scripts/morphometrics/pub_figures_utils.R`, `scripts/core/db_utils.R` | `v_integrative_triage_master` (master relational view) | `curation_triage`, `final_curation_triage.csv` | Out-of-fold kappa weighting; Hopkins anthesis envelopes; 5-stage decision matrix; Rules A, B, and C type controls. |
| **08: Complexes** | `scripts/analysis/08_discover_taxon_complexes.R` | `scripts/morphometrics/authority_weighting_utils.R`, `scripts/core/db_utils.R` | `v_integrative_triage_master`, multimodal flags | `discovered_taxon_complexes.csv`, network graph | Multimodal affinity network; type-shrunk centroids; unsupervised Leiden/Louvain modularity and Ward.D2 clustering. |

---

## 2. Refactored & LeafMachine3 Repository Directory Tree

```text
Packera_dubia_morphometrics/
├── setup_leafmachine3.sh              # LeafMachine3 isolated environment setup script (uv 0.12.23)
├── main.py                            # Unified Python CLI entrypoint (check-env, harvest, segment, etc.)
├── run_pipeline.sh                    # Master run-scoped shell orchestrator CLI
├── config/
│   └── config.yaml                    # Central paths, taxon synonyms, authority patterns, thresholds
├── docs/
│   └── TAXONOMIC_SOP.md               # Scientific Standard Operating Procedure (equations, citations, theory)
├── data/                              # Tier A: Canonical Shared Store (invariant data & central DB)
│   ├── packera_pipeline.sqlite        # Central transactional SQLite WAL relational store
│   ├── raw_vouchers/                  # High-resolution voucher images (*.webp)
│   ├── masks/                         # Canonical binary leaf silhouettes (tier1_intact/, tier2_reflected/)
│   ├── contours/                      # Canonical 120-vertex CSV contour coordinates
│   ├── cropped_patches/               # Organ crops (basal_leaves_raw/, capitula/, rosettes_dense/)
│   ├── environmental/                 # EPA L3 ecoregions, SSURGO, POLARIS, SOLUS100, 3DEP, SGMC stack
│   └── tables/                        # Step 01–02 canonical CSV manifests (curated_vouchers.csv, etc.)
├── runs/                              # Tier B: Run-Scoped Analysis Directories (Steps 03–08)
│   ├── latest -> run_<TIMESTAMP>      # Symlink to most recent run directory
│   └── <RUN_ID>/                      # Isolated run workspace (data/tables/, outputs/figures/, reports/, logs/)
├── scripts/
│   ├── core/                          # config.py, db.py, db_schema.py, harvester.py, epa_ecoregions.py,
│   │                                  # taxonomy_authority.py, artifact_manager.py, logger.py, db_utils.R
│   ├── vision/                        # configure_leafmachine3.py, lm3_data_loader.py, lm3_inference_engines.py,
│   │                                  # leaf_geometry_utils.py, leaf_gatekeeper.py, organ_crown_metrics.py,
│   │                                  # capitulum_phenology_classifier.py, step02_figures.py
│   ├── pipeline/                      # 01_voucher_harvester.py, 02_segment_and_extract.py
│   ├── morphometrics/                 # 03_fourier_extractor.R, 04_gmm_morphotools.R, authority_weighting_utils.R,
│   │                                  # dissected_leaf_utils.R, fourier_parallel_utils.R, pub_figures_utils.R,
│   │                                  # pub_ecoregion_manova_utils.R, pub_synthesis_plates.R
│   ├── analysis/                      # 05_cleanlab_vision_xai.py, 06_multimodal_spatial_rf.R, step06_phenology_ecotypes.R,
│   │                                  # step06_sdm_maxent.R, 07_triage_dashboard_synthesis.R, 08_discover_taxon_complexes.R,
│   │                                  # edaphic_fine_grain_utils.R, vision_dissection_utils.py
│   └── tests/                         # Full automated test suite (pytest scripts/tests/ -q)
├── requirements.txt                   # Primary Python dependencies
└── README.md                          # Single operational, setup, and CLI manual
```

---

## 3. Cross-Platform Environment & Setup (Linux, macOS & Windows WSL 2)

Supported across native Linux, macOS (Apple Silicon / Intel), and Windows via WSL 2.

### 3.1 Windows WSL 2 Setup & The Critical ext4 Storage Rule

For Windows users, WSL 2 with Ubuntu 22.04 LTS (or newer) is the officially supported environment:
1. **Install WSL 2 & Ubuntu:** `wsl --install -d Ubuntu-22.04`
2. **CRITICAL ext4 Storage Rule:**
   > [!IMPORTANT]
   > You **MUST** clone and run this repository inside the native Linux `ext4` filesystem (`~/Packera_dubia_morphometrics`), **NEVER on a Windows mounted drive** (`/mnt/c/...`).
   > Running on `/mnt/c/` causes severe 9P bridging bottlenecks, breaks SQLite WAL concurrent file locks, and corrupts raster and image extraction.
3. **NVIDIA GPU Acceleration in WSL 2:** Install NVIDIA driver on Windows host ($\ge 528.33$ Windows / $\ge 525.60$ Linux). No separate Linux driver or CUDA toolkit needed inside WSL 2; precompiled `onnxruntime-gpu` wheels automatically interface with the virtualized GPU.

### 3.2 System Dependencies

```bash
# Ubuntu / Debian / WSL 2:
sudo apt update && sudo apt install -y build-essential git curl libgl1 libglib2.0-0 r-base r-base-dev libgdal-dev libgeos-dev libproj-dev libudunits2-dev libssl-dev libxml2-dev

# macOS (Homebrew):
brew install r gdal geos proj udunits
```

### 3.3 Dual Virtual Environment Setup

#### 1. Primary Python Environment (`.venv`)
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install --upgrade pip && pip install -r requirements.txt
```

#### 2. LeafMachine3 ONNX Runtime Environment (`.venv_LM3`)
LeafMachine3 is provisioned automatically via `setup_leafmachine3.sh`, which installs `uv 0.12.23`, creates `.venv_LM3`, installs model checkpoints, and runs diagnostics:

```bash
bash setup_leafmachine3.sh --extra gpu    # Linux / Windows WSL 2 (NVIDIA GPU)
bash setup_leafmachine3.sh --extra cpu    # CPU-only systems
bash setup_leafmachine3.sh --extra macos  # Apple Silicon macOS
```
Performs: (1) Installs `uv 0.12.23`; (2) `uv sync --frozen --extra <gpu|cpu|macos>`; (3) Downloads models via `uv run --frozen --no-sync lm3 models install --yes`; (4) Verifies hardware via `uv run --frozen --no-sync lm3 doctor`.

#### 3. R Statistical Dependencies
```r
install.packages(c("DBI", "RSQLite", "tidyverse", "remotes", "optparse", "mclust", "clustvarsel", "RRPP", "vegan", "terra", "sf", "geodata", "soilDB", "spatialRF", "ranger", "adespatial", "blockCV", "ecospat", "ENMTools", "igraph", "patchwork", "viridis"))
remotes::install_github("MomX/Momocs")
remotes::install_github("MarekSlenker/MorphoTools2")
```

---

## 4. Quick-Start Terminal Cheat-Sheet & Complete CLI Flag Reference

### Quick-Start Execution Patterns

```bash
source .venv/bin/activate

# 1. Environment preflight check & SQLite bootstrap:
python main.py check-env
python3 -c "from scripts.core.db import ingest_existing_csvs_to_db; print(ingest_existing_csvs_to_db())"

# 2. Pipeline runs: default genus-wide, custom taxon subset, or resume from step:
./run_pipeline.sh
./run_pipeline.sh --run-id packera_subset --taxon "Packera dubia" "Packera anonyma" "Packera aurea"
./run_pipeline.sh --run-id packera_revision_v1 --start-step 3

# 3. Focal species, cryptic delimitation, or complex discovery:
./run_pipeline.sh --run-id dubia_triage --start-step 2 --focal-taxon "Packera dubia" --gold-anchors 30
./run_pipeline.sh --run-id dubia_cryptic --start-step 2 --focal-taxon "Packera dubia" --intraspecific-mode
./run_pipeline.sh --run-id complex_discovery --start-step 2 --discover-complexes
```

---

### Tier 1: Master Run-Scoped Shell Orchestrator (`./run_pipeline.sh`)

| Category | Flag | Arguments | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| **Run & Compute** | `--run-id` | `<ID>` | `run_<TIMESTAMP>` | Run workspace directory under `runs/<ID>/`; updates `runs/latest`. |
| | `--db` | `<PATH>` | `data/packera_pipeline.sqlite` | SQLite database path (`PACKERA_DB_PATH`). |
| | `--demo` | *(none)* | `false` | Ephemeral test mode (`PACKERA_DEMO_MODE=1`, writes to `demo_ephemeral.sqlite`). |
| | `-j`, `--threads` | `<INT>` | `auto` | CPU worker thread limit across Python and R stages (`PACKERA_MAX_WORKERS`). |
| | `--gpu-batch-size` | `<INT>` | `16` / `64` | GPU batch size for Step 02 LeafMachine3 and Step 05 DINOv2. |
| | `--start-step` | `<1..8>` | `1` | Start/resume step; skips Step 01 when $\ge 2$. |
| | `--force` | *(none)* | `false` | Overwrite cached artifacts across active phases. |
| | `--limit` | `<N>` | `0` *(unlimited)* | Max vouchers processed in Step 02 segmentation. |
| | `--regenerate-figures-only` | *(none)* | `false` | Re-render Step 02 figures from validated silhouettes without re-segmenting. |
| **Taxonomy & Scope** | `--focal-taxon` | `"<Binomial>"` | *(none)* | Target focal species for targeted CDA, niche modeling, and synthesis. |
| | `--taxon` | `"<Sp1>" "<Sp2>" ...` | *(all / genus-wide)* | Explicitly specify target binomials (omit for genus-wide discovery). |
| | `--exclude-taxon` | `"<Binomial>"` | *(none)* | Exclude one or more taxa from active run (repeatable). |
| | `--intraspecific-mode` | *(none)* | `false` | Run single-taxon intraspecific morphotype/ecotype delimitation. |
| | `--expand-synonyms` | *(none)* | `true` | Toggle expansion to historical *Senecio* / *Cineraria* basionyms. |
| | `--discover-complexes` | *(none)* | `false` | Trigger Step 08 unsupervised Leiden/Louvain and Ward.D2 complex discovery. |
| | `--min-affinity` | `<FLOAT>` | `0.03` | Minimum multimodal affinity weight $W_{ij}$ for complex graph edges. |
| **Authority & Sampling** | `--gold-anchors` | `<N>` | `30` | Tier 1 Gold reference anchors per outgroup taxon in focal runs. |
| | `--determiner-tier` | `<TIER>` | *(all)* | Filter vouchers to an authority tier (`Tier_1_Gold`, `Tier_2_Silver`, `Tier_3_Bronze`). |
| | `--exclude-misidentified` | *(none)* | `false` | Prune vouchers flagged as misidentified in prior runs. |
| | `--geographic-scope` | `{eastern_na,north_america,global}` | `eastern_na` | Longitudinal envelope (`eastern_na` enforces $\ge -106.65^\circ\text{W}$). |
| | `--countries` | `"<ISO2 ...>"` | `US` | Space- or comma-separated ISO-2 country codes for GBIF (e.g., `"US CA MX"`). |
| | `--spatial-thin-km` | `<FLOAT>` | `0.0` | Spatial grid thinning distance in km per species. |
| | `--max-per-state` | `<INT>` | *(unlimited)* | Cap vouchers retained per `(species, stateProvince)`. |
| | `--max-per-region` | `<INT>` | *(unlimited)* | Cap vouchers retained per EPA Level III Ecoregion. |
| | `--no-spatial-weights` | *(none)* | `false` | Disable inverse spatial-density weighting during model fitting. |

---

### Tier 2: Unified Python CLI (`python main.py <subcommand>`)

#### Global Flags
`-c/--config <PATH>` (default: `config/config.yaml`), `-v/--verbose`, `--run-id <ID>`, `--db <PATH>`, `--demo`, `-j/--threads <INT>`, `--gpu-batch-size <INT>`.

#### Subcommands Reference
| Subcommand | Syntax & Key Options | Operational Scope |
| :--- | :--- | :--- |
| **`check-env`** | `python main.py check-env [--strict]` | Verifies Python $\ge 3.10$, CUDA, R, core packages, and LeafMachine3. |
| **`download-weights`** | `python main.py download-weights [--force]` | Verifies or downloads LeafMachine3 ONNX model checkpoints with SHA-256 checks. |
| **`harvest`** | `python main.py harvest [--taxon ...] [--geographic-scope ...] [--download-images]` | Queries GBIF/iDigBio, standardizes synonyms, assigns EPA L3 ecoregions. |
| **`segment`** | `python main.py segment [--vouchers PATH] [--device cuda\|cpu] [--limit N]` | Runs LeafMachine3 ONNX detection/segmentation, gatekeeping, per-voucher commits. |
| **`morphometrics`** | `python main.py morphometrics [--harmonics K] [--num-pcs P] [--focal-taxon T]` | Step 03 EFA ($H_{\text{opt}} \ge 99\%$ power) and Step 04 allometry, GMM, dual-track CDA. |
| **`synthesis`** | `python main.py synthesis [--cleanlab-threshold FLOAT] [--export-figures]` | Step 05 DINOv2 patch XAI, Step 06 Spatial RF, Step 07 consensus triage. |
| **`run-all`** | `python main.py run-all [combined options]` | Full end-to-end execution from GBIF harvest through 5-stage revision synthesis. |

---

### Tier 3: Standalone Per-Step Script Flags

- **Step 01 (`scripts/pipeline/01_voucher_harvester.py`)**: `--taxon <NAMES>`, `--expand-synonyms`, `--geographic-scope <SCOPE>`, `--countries <CODES>`, `--max-records <N>`, `--max-uncertainty <M>`, `--download-images`, `--spatial-thin-km <KM>`, `--max-per-state <N>`, `--max-per-region <N>`, `--force`, `--db <PATH>`.
- **Step 02 (`scripts/pipeline/02_segment_and_extract.py`)**: `--vouchers <PATH>`, `--device <cuda|cpu>`, `-j/--threads <INT>`, `--gpu-batch-size <INT>`, `--enable-geometric-qc`, `--min-solidity <FLOAT>`, `--min-ucs <FLOAT>`, `--limit <N>`, `--regenerate-figures-only`, `--force`, `--db <PATH>`.
- **Step 03 (`scripts/morphometrics/03_fourier_extractor.R`)**: `-i/--input <DIR>`, `-f/--manifest <PATH>`, `-v/--vouchers <PATH>`, `-o/--output <PATH>`, `-k/--harmonics <N>`, `-p/--num-pcs <N>`, `--symmetrize-efa-harmonics`, `--db <PATH>`.
- **Step 04 (`scripts/morphometrics/04_gmm_morphotools.R`)**: `-v/--vouchers <PATH>`, `-e/--harmonics <PATH>`, `-m/--reproductive <PATH>`, `--correct-allometry <auto|true|false>`, `-k/--max-k <INT>`, `-d/--num-pcs <INT>`, `--focal-taxon <NAME>`, `--intraspecific-mode`, `--db <PATH>`.
- **Step 05 (`scripts/analysis/05_cleanlab_vision_xai.py`)**: `--rosette-dir <DIR>`, `--vouchers-csv <PATH>`, `--manifest-csv <PATH>`, `--output-csv <PATH>`, `-j/--threads <INT>`, `--gpu-batch-size <INT>`, `--cleanlab-threshold <FLOAT>`, `--focal-taxon <NAME>`, `--intraspecific-mode`, `--demo`, `--db <PATH>`.
- **Step 06 (`scripts/analysis/06_multimodal_spatial_rf.R`)**: `-v/--vouchers <PATH>`, `-m/--morphometrics <PATH>`, `-n/--vision-audit <PATH>`, `--pheno-states <PATH>`, `-e/--env-dir <DIR>`, `-k/--permutations <N>`, `-t/--n-trees <INT>`, `--focal-taxon <NAME>`, `--intraspecific-mode`, `--demo`, `--db <PATH>`.
- **Step 07 (`scripts/analysis/07_triage_dashboard_synthesis.R`)**: `-v/--vouchers <PATH>`, `-m/--morphometrics <PATH>`, `-n/--vision-audit <PATH>`, `-c/--multimodal-flags <PATH>`, `--focal-taxon <NAME>`, `--intraspecific-mode`, `--demo`, `--db <PATH>`.
- **Step 08 (`scripts/analysis/08_discover_taxon_complexes.R`)**: `-m/--morphometrics <PATH>`, `-e/--harmonics <PATH>`, `-n/--vision-audit <PATH>`, `-v/--vouchers <PATH>`, `-c/--multimodal <PATH>`, `--focal-taxon <NAME>`, `--min-affinity <FLOAT>`, `--db <PATH>`.

---

## 5. Publication Figure Suite Inventory

Every execution exports figures in vector `.pdf` and 300-DPI `.png` formats inside `${RUN_DIR}/outputs/figures/` using Okabe-Ito palettes, unfilled hollow markers, and 3D ordinations with drop-lines:

| Figure Stem (`.pdf` & `.png`) | Script Source | Biological / Statistical Hypothesis Visualized |
| :--- | :--- | :--- |
| `fig02a_segmentation_qc_distributions` | `scripts/pipeline/02_segment_and_extract.py` | Solidity vs. UCS by routing tier; basal leaf yield per voucher & QC breakdown. |
| `fig02b_silhouette_montage` | `scripts/pipeline/02_segment_and_extract.py` | 16-tile montage of basal leaf silhouettes (`Tier_1_Direct` and `Tier_2_Reflected`) with $p_{\text{base}}$ anchor. |
| `fig03_efa_morphospace_reconstruction` | `scripts/morphometrics/03_fourier_extractor.R` | Inverse Fourier contour reconstructions across $H_{\text{opt}} \in [6, 24]$ and eigenshape deformation vectors across PC1–PC3. |
| `fig03a_efa_morphospace_pca` | `scripts/morphometrics/03_fourier_extractor.R` | 3D Full Elliptic Fourier PCA and 3D bilateral-symmetric morphospace (`Sym_PC1 x Sym_PC2 x Sym_PC3`). |
| `fig03b_harmonic_power_and_asymmetry` | `scripts/morphometrics/03_fourier_extractor.R` | Cumulative harmonic power spectrum ($H_{\text{opt}} \ge 99\%$ power) and bilateral asymmetry vs. roughness. |
| `fig04a_gmm_bic_and_uncertainty` | `scripts/morphometrics/04_gmm_morphotools.R` | Unsupervised `mclust` BIC curves ($2\Delta\text{BIC}$) and 3D morphospace with uncertainty scaling. |
| `fig04b_cda_biplot_and_allometry` | `scripts/morphometrics/04_gmm_morphotools.R` | 3D CDA biplot across Horn/Williams-Titus subspace with honest LOOCV; Klingenberg allometric scaling. |
| `fig04_gmm_cda_clustering` | `scripts/morphometrics/04_gmm_morphotools.R` | Track 1 Unsupervised GMM Morphotype CDA biplot displaying discrete morphotypes with 95% ellipses. |
| `fig04c_supervised_ecoregion_cda` | `scripts/morphometrics/04_gmm_morphotools.R` | Track 2 Supervised EPA Level III Ecoregion CDA biplot displaying canonical ordination and MANOVA statistics. |
| `shape_allometry_regression` | `scripts/morphometrics/04_gmm_morphotools.R` | Bivariate shape PC vs. $\log_{10}(\text{CS})$ allometric regressions testing slope homogeneity across taxa. |
| `GradCAM_audit_panel` | `scripts/analysis/05_cleanlab_vision_xai.py` | 3-panel XAI audit: [1] Indumentum-blind crop; [2] `dinov2_vitb14_reg` $16 \times 16$ patch attention overlay; [3] Posteriors. |
| `fig05b_cleanlab_quality_and_tomentum` | `scripts/analysis/05_cleanlab_vision_xai.py` | Cleanlab label quality distributions and CIE $L^*$ luminance vs. HSV saturation separating tomentum discordance. |
| `spatial_rf_niche_importance` | `scripts/analysis/06_multimodal_spatial_rf.R` | Spatial RF out-of-bag permutation importance (`MeanDecreaseAccuracy`) identifying primary environmental predictors. |
| `fig06_ecological_niche_phenology` | `scripts/analysis/06_multimodal_spatial_rf.R` | 4-panel niche plate: climatic vs. edaphic overlap ($D_{\text{clim}}$ vs. $D_{\text{edaph}}$), SGMC lithology, anthesis phenology. |
| `fig06b_geographic_map_and_hopkins_phenology` | `scripts/analysis/06_multimodal_spatial_rf.R` | Geographic voucher map colored by taxon and collection DOY vs. Latitude against Hopkins Bioclimatic Law. |
| `fig06b_sdm_habitat_suitability_map` | `scripts/analysis/06_multimodal_spatial_rf.R` | 4-panel MaxEnt (`cloglog`) and Spatial RF ensemble habitat suitability plate with authority tier symbols and Type overlays. |
| `fig07_multimodal_discordance_summary` | `scripts/analysis/07_triage_dashboard_synthesis.R` | Modality-specific conflict rates weighted by out-of-fold Cohen's $\kappa$; consensus triage priority breakdown. |
| `Figure_Integrative_Packera_Revision` | `scripts/analysis/07_triage_dashboard_synthesis.R` | 4-panel synthesis dashboard: 3D CDA morphospace, DINOv2 confidence, Hopkins anthesis residuals, 5-stage treatment. |
| `fig_empirical_taxon_complex_network` | `scripts/analysis/08_discover_taxon_complexes.R` | Unsupervised multimodal affinity network graph depicting species complexes and Taxon x GMM Morphotype subclusters. |

---

## 6. License & Attribution

This project is licensed under the **MIT License**. Herbarium specimen imagery harvested through the pipeline remains subject to the data access agreements and institutional copyright policies of the contributing herbaria (NCU, NY, US, MO, GA, VDB, etc.).
