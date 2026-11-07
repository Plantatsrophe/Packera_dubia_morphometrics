# Automated Vegetative Phenotyping and Herbarium Triage Pipeline for North American *Packera*

### Quantitative Morphological Diagnosability and Species Boundaries under the Unified Species Concept

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/) [![R 4.3+](https://img.shields.io/badge/R-4.3+-blue.svg)](https://www.r-project.org/) [![LeafMachine3](https://img.shields.io/badge/LeafMachine3-ONNX_Runtime-brightgreen.svg)](https://github.com/Gene-Weaver/LeafMachine3)

- **Principal Investigator:** J. Brandon Fuller (PhD Candidate, Department of Biology, UNC-CH)
- **Faculty Advisor:** Dr. Alan S. Weakley (Director, NCU Herbarium; UNC Biology)
- **Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV5` (Detailed theory, equations, and proofs in [`docs/TAXONOMIC_SOP.md`](docs/TAXONOMIC_SOP.md))
- **Focal Dissertation Taxon:** ***Packera dubia*** (Spreng.) Trock & Mabb. (Asteraceae: Senecioneae) evaluated dynamically against all regional and continent-wide North American congeners.

---

## 1. Project Overview and 8-Step Architecture

I built this computational pipeline to conduct high-throughput vegetative phenotyping, geometric morphometrics, deep vision explainable AI, and integrative taxonomic triage across North American *Packera*. The pipeline evaluates morphological diagnosability within the theoretical framework of the **Unified Species Concept (USC; de Queiroz 2007)** through **Integration by Congruence (Padial et al. 2010)**. I analyze basal leaf morphometrics alongside capitulum phenology, nextRAD phylogenomic evidence, and multi-scale 0 to 15 cm rhizosphere edaphic profiles.

All pipeline stages synchronize through a central transactional SQLite database (`data/packera_pipeline.sqlite`) in WAL mode and atomic CSV snapshots under a Canonical Two-Tier Storage Architecture. For complete mathematical formulations, statistical proofs, and biological justifications, refer to [`docs/TAXONOMIC_SOP.md`](docs/TAXONOMIC_SOP.md).

### 8-Step Pipeline Architecture

| Step / Phase | Entrypoint Script | Helper Modules | Key Inputs | Key Outputs | Core Purpose |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **01: Ingestion** | `scripts/pipeline/01_voucher_harvester.py` | `scripts/core/` (`harvester.py`, `epa_ecoregions.py`, `taxonomy_authority.py`, `db.py`) | GBIF and iDigBio API, Darwin Core queries, EPA L3 ecoregions | `vouchers`, `taxon_synonyms`, `curated_vouchers.csv` | I harvest occurrences from GBIF and iDigBio, standardize synonyms, rank nomenclatural types, stratify determiner authority, and assign EPA Level III ecoregions. |
| **02: Extraction** | `scripts/pipeline/02_segment_and_extract.py` | `scripts/vision/` (`configure_leafmachine3.py`, `lm3_data_loader.py`, `lm3_inference_engines.py`, `leaf_gatekeeper.py`, `organ_crown_metrics.py`) | `vouchers`, herbarium images, LeafMachine3 ONNX models | `plant_individuals`, `extracted_leaves`, masks, crops | I segment individual plants and leaves with LeafMachine3 ONNX engines, apply botanical gatekeepers, extract petiole dimensions, and record phenology. |
| **03: Fourier EFA** | `scripts/morphometrics/03_fourier_extractor.R` | `scripts/morphometrics/` (`authority_weighting_utils.R`, `fourier_parallel_utils.R`), `scripts/core/db_utils.R` | `extracted_leaves.contour_xy_json`, `vouchers` | `leaf_efa_harmonics`, `leaf_efa_harmonics.csv` | I read contours from SQLite, align outlines along the major axis, calibrate harmonic count at 99% power, symmetrize harmonics, and calculate OTU medians. |
| **04: Morphometrics** | `scripts/morphometrics/04_gmm_morphotools.R` | `scripts/morphometrics/` (`authority_weighting_utils.R`, `pub_figures_utils.R`), `scripts/core/db_utils.R` | `leaf_efa_harmonics`, `vouchers`, phenology metrics | `otu_morphometrics`, `morphometric_flags.csv` | I evaluate allometric slope homogeneity, select subspaces with Horn's analysis, weight determiners via Dawid-Skene priors, and run dual-track CDA. |
| **05: Vision XAI** | `scripts/analysis/05_cleanlab_vision_xai.py` | `scripts/analysis/vision_dissection_utils.py`, `scripts/core/db.py` | `vouchers`, `plant_individuals`, `cropped_patches/` | `dinov2_feature_cache`, `vision_xai_audit` | I apply low-pass filtering to neutralize abraded hairs, pool foreground patch tokens with DINOv2, build type prototypes, and score label noise. |
| **06: Spatial RF** | `scripts/analysis/06_multimodal_spatial_rf.R` | `scripts/analysis/` (`step06_phenology_ecotypes.R`, `step06_sdm_maxent.R`, `edaphic_fine_grain_utils.R`) | `otu_morphometrics`, `vision_xai_audit`, WorldClim, SSURGO | `spatial_rf_predictions`, `multimodal_conflict_flags.csv` | I extract 19 BioClim variables and 0 to 15 cm rhizosphere soils, calculate type locality baselines, prune collinear predictors, and run spatial block CV. |
| **07: Triage** | `scripts/analysis/07_triage_dashboard_synthesis.R` | `scripts/morphometrics/pub_figures_utils.R`, `scripts/core/db_utils.R` | `v_integrative_triage_master` (master relational view) | `curation_triage`, `final_curation_triage.csv` | I weight modalities by out-of-fold kappa reliability, fit circular anthesis windows, apply a 5-stage decision matrix, and enforce type controls. |
| **08: Complexes** | `scripts/analysis/08_discover_taxon_complexes.R` | `scripts/morphometrics/authority_weighting_utils.R`, `scripts/core/db_utils.R` | `v_integrative_triage_master`, multimodal flags | `discovered_taxon_complexes.csv`, network graph | I build an affinity network across all modalities, pull centroids toward type specimens, and run unsupervised Leiden clustering to identify species complexes. |

---

## 2. Refactored and LeafMachine3 Repository Directory Tree

```text
Packera_dubia_morphometrics/
├── setup_leafmachine3.sh              # LeafMachine3 isolated environment setup script (uv 0.12.23)
├── main.py                            # Unified Python CLI entrypoint (check-env, harvest, segment, etc.)
├── run_pipeline.sh                    # Master run-scoped shell orchestrator CLI
├── config/
│   └── config.yaml                    # Central paths, taxon synonyms, authority patterns, thresholds
├── docs/
│   └── TAXONOMIC_SOP.md               # Scientific Standard Operating Procedure (equations, citations, theory)
├── data/                              # Tier A: Canonical Shared Store (invariant data and central DB)
│   ├── packera_pipeline.sqlite        # Central transactional SQLite WAL relational store
│   ├── raw_vouchers/                  # High-resolution voucher images (*.webp)
│   ├── masks/                         # Canonical binary leaf silhouettes (tier1_intact/, tier2_reflected/)
│   ├── contours/                      # Canonical 120-vertex CSV contour coordinates
│   ├── cropped_patches/               # Organ crops (basal_leaves_raw/, capitula/, rosettes_dense/)
│   ├── environmental/                 # EPA L3 ecoregions, SSURGO, POLARIS, SOLUS100, 3DEP, SGMC stack
│   └── tables/                        # Step 01 to 02 canonical CSV manifests (curated_vouchers.csv, etc.)
├── runs/                              # Tier B: Run-Scoped Analysis Directories (Steps 03 to 08)
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

## 3. Cross-Platform Environment and Setup (Linux, macOS and Windows WSL 2)

I tested and verified the pipeline on native Linux, macOS on Apple Silicon and Intel, and Windows through WSL 2.

### 3.1 Windows WSL 2 Setup and The Critical ext4 Storage Rule

For Windows users, WSL 2 with Ubuntu 22.04 LTS (or newer) provides the recommended execution environment:
1. **Install WSL 2 and Ubuntu:** Run `wsl --install -d Ubuntu-22.04` in an Administrator PowerShell prompt.
2. **CRITICAL ext4 Storage Rule:**
   > [!IMPORTANT]
   > You must clone and run this repository inside the native Linux `ext4` filesystem (`~/Packera_dubia_morphometrics`), and never on a Windows mounted drive (`/mnt/c/...`).
   > Running on `/mnt/c/` causes severe 9P bridging bottlenecks, breaks SQLite WAL concurrent file locks, and corrupts raster and image extraction routines.
3. **NVIDIA GPU Acceleration in WSL 2:** Install the latest NVIDIA driver on the Windows host ($\ge 528.33$ on Windows or $\ge 525.60$ on Linux). Do not install a separate Linux display driver or CUDA toolkit inside WSL 2, as precompiled `onnxruntime-gpu` wheels automatically interface with the virtualized GPU.

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
I provision LeafMachine3 automatically through `setup_leafmachine3.sh`, which installs `uv 0.12.23`, creates `.venv_LM3`, installs model checkpoints, and runs diagnostics:

```bash
bash setup_leafmachine3.sh --extra gpu    # Linux or Windows WSL 2 with NVIDIA GPU
bash setup_leafmachine3.sh --extra cpu    # CPU systems
bash setup_leafmachine3.sh --extra macos  # Apple Silicon macOS
```
This script completes four tasks: (1) installs standalone `uv 0.12.23`; (2) synchronizes dependencies with `uv sync --frozen --extra <gpu|cpu|macos>`; (3) downloads model checkpoints through `uv run --frozen --no-sync lm3 models install --yes`; (4) verifies hardware through `uv run --frozen --no-sync lm3 doctor`.

#### 3. R Statistical Dependencies
```r
install.packages(c("DBI", "RSQLite", "tidyverse", "remotes", "optparse", "mclust", "clustvarsel", "RRPP", "vegan", "terra", "sf", "geodata", "soilDB", "spatialRF", "ranger", "adespatial", "blockCV", "ecospat", "ENMTools", "igraph", "patchwork", "viridis"))
remotes::install_github("MomX/Momocs")
remotes::install_github("MarekSlenker/MorphoTools2")
```

---

## 4. Quick-Start Terminal Cheat-Sheet and Complete CLI Flag Reference

### Quick-Start Execution Patterns

```bash
source .venv/bin/activate

# 1. Environment preflight check and SQLite bootstrap:
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
| **Run and Compute** | `--run-id` | `<ID>` | `run_<TIMESTAMP>` | I assign this run workspace directory under `runs/<ID>/` and update `runs/latest`. |
| | `--db` | `<PATH>` | `data/packera_pipeline.sqlite` | I set the SQLite database file path through this argument or `PACKERA_DB_PATH`. |
| | `--demo` | *(none)* | `false` | I isolate writes into `demo_ephemeral.sqlite` when running in demo test mode. |
| | `-j`, `--threads` | `<INT>` | `auto` | I cap CPU worker thread counts across Python and R stages through this setting. |
| | `--gpu-batch-size` | `<INT>` | `16` / `64` | I set GPU inference batch sizes for Step 02 LeafMachine3 and Step 05 DINOv2. |
| | `--start-step` | `<1..8>` | `1` | I begin or resume execution at this step, automatically bypassing Step 01 when $\ge 2$. |
| | `--force` | *(none)* | `false` | I force full re-execution and overwrite existing cached artifacts across active phases. |
| | `--limit` | `<N>` | `0` | I limit the total count of vouchers processed during Step 02 segmentation. |
| | `--regenerate-figures-only` | *(none)* | `false` | I rebuild Step 02 diagnostic figures from stored silhouettes without running segmentation. |
| **Taxonomy and Scope** | `--focal-taxon` | `"<Binomial>"` | *(none)* | I target this single species for focused CDA ordination, niche models, and triage. |
| | `--taxon` | `"<Sp1>" "<Sp2>" ...` | *(all)* | I restrict analysis strictly to these named species instead of genus-wide discovery. |
| | `--exclude-taxon` | `"<Binomial>"` | *(none)* | I drop one or more named taxa from the active cohort (repeatable flag). |
| | `--intraspecific-mode` | *(none)* | `false` | I restrict Steps 04 to 07 to single-species morphotype and ecotype delimitation. |
| | `--expand-synonyms` | *(none)* | `true` | I toggle automatic query expansion to historical basionyms in *Senecio* and *Cineraria*. |
| | `--discover-complexes` | *(none)* | `false` | I trigger Step 08 unsupervised Leiden modularity and Ward.D2 complex discovery. |
| | `--min-affinity` | `<FLOAT>` | `0.03` | I require this minimum multimodal affinity weight $W_{ij}$ to form complex graph edges. |
| **Authority and Sampling** | `--gold-anchors` | `<N>` | `30` | I draw up to this many Tier 1 Gold anchors per outgroup species during focal runs. |
| | `--determiner-tier` | `<TIER>` | *(all)* | I restrict vouchers strictly to one authority tier (`Tier_1_Gold`, `Tier_2_Silver`, etc.). |
| | `--exclude-misidentified` | *(none)* | `false` | I remove vouchers that were flagged as misidentified during earlier pipeline runs. |
| | `--geographic-scope` | `{eastern_na,north_america,global}` | `eastern_na` | I enforce this longitudinal bounding box (`eastern_na` requires longitude $\ge -106.65^\circ\text{W}$). |
| | `--countries` | `"<ISO2 ...>"` | `US` | I pass these ISO-2 country codes to GBIF occurrence harvesting (such as `"US CA MX"`). |
| | `--spatial-thin-km` | `<FLOAT>` | `0.0` | I thin occurrences by this geographic grid distance in kilometers for each species. |
| | `--max-per-state` | `<INT>` | *(unlimited)* | I cap the maximum voucher count retained for each `(species, stateProvince)` combination. |
| | `--max-per-region` | `<INT>` | *(unlimited)* | I cap the maximum voucher count retained for each EPA Level III Ecoregion. |
| | `--no-spatial-weights` | *(none)* | `false` | I disable inverse spatial density weighting during statistical model fitting. |

---

### Tier 2: Unified Python CLI (`python main.py <subcommand>`)

#### Global Flags
I support `-c/--config <PATH>` (default: `config/config.yaml`), `-v/--verbose`, `--run-id <ID>`, `--db <PATH>`, `--demo`, `-j/--threads <INT>`, and `--gpu-batch-size <INT>`.

#### Subcommands Reference
| Subcommand | Syntax and Key Options | Operational Scope |
| :--- | :--- | :--- |
| **`check-env`** | `python main.py check-env [--strict]` | I verify Python versions, CUDA availability, R packages, and LeafMachine3 runtime health. |
| **`download-weights`** | `python main.py download-weights [--force]` | I verify and download LeafMachine3 ONNX model weights with SHA-256 integrity checks. |
| **`harvest`** | `python main.py harvest [--taxon ...] [--geographic-scope ...] [--download-images]` | I harvest records from GBIF and iDigBio, standardize synonyms, and assign EPA ecoregions. |
| **`segment`** | `python main.py segment [--vouchers PATH] [--device cuda\|cpu] [--limit N]` | I run LeafMachine3 ONNX detection and segmentation, apply gatekeepers, and save records. |
| **`morphometrics`** | `python main.py morphometrics [--harmonics K] [--num-pcs P] [--focal-taxon T]` | I run Step 03 EFA at 99% harmonic power followed by Step 04 allometry, GMM, and CDA. |
| **`synthesis`** | `python main.py synthesis [--cleanlab-threshold FLOAT] [--export-figures]` | I execute Step 05 DINOv2 vision XAI, Step 06 Spatial RF, and Step 07 consensus triage. |
| **`run-all`** | `python main.py run-all [combined options]` | I run the complete end-to-end pipeline from initial harvesting through taxonomic triage. |

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

| Figure Stem (`.pdf` and `.png`) | Script Source | Biological and Statistical Hypothesis Visualized |
| :--- | :--- | :--- |
| `fig02a_segmentation_qc_distributions` | `scripts/pipeline/02_segment_and_extract.py` | Solidity versus UCS distributions across routing tiers, basal leaf yields, and QC filter breakdown. |
| `fig02b_silhouette_montage` | `scripts/pipeline/02_segment_and_extract.py` | 16-tile montage displaying verified basal leaf silhouettes with basal anchor points and midvein traces. |
| `fig03_efa_morphospace_reconstruction` | `scripts/morphometrics/03_fourier_extractor.R` | Inverse Fourier contour reconstructions across 6 to 24 harmonics and eigenshape deformation vectors for PC1 to PC3. |
| `fig03a_efa_morphospace_pca` | `scripts/morphometrics/03_fourier_extractor.R` | 3D Full Elliptic Fourier PCA and 3D bilateral-symmetric morphospace (`Sym_PC1` by `Sym_PC2` by `Sym_PC3`). |
| `fig03b_harmonic_power_and_asymmetry` | `scripts/morphometrics/03_fourier_extractor.R` | Cumulative harmonic power spectrum reaching 99% power and bilateral asymmetry plotted against roughness. |
| `fig04a_gmm_bic_and_uncertainty` | `scripts/morphometrics/04_gmm_morphotools.R` | Unsupervised `mclust` BIC model selection curves and 3D morphospace with uncertainty point scaling. |
| `fig04b_cda_biplot_and_allometry` | `scripts/morphometrics/04_gmm_morphotools.R` | 3D CDA biplot across the Horn and Williams-Titus subspace alongside Klingenberg allometric regressions. |
| `fig04_gmm_cda_clustering` | `scripts/morphometrics/04_gmm_morphotools.R` | Track 1 Unsupervised GMM Morphotype CDA ordination biplot displaying discrete morphotypes with 95% ellipses. |
| `fig04c_supervised_ecoregion_cda` | `scripts/morphometrics/04_gmm_morphotools.R` | Track 2 Supervised EPA Level III Ecoregion CDA biplot displaying regional ordination and MANOVA tests. |
| `shape_allometry_regression` | `scripts/morphometrics/04_gmm_morphotools.R` | Bivariate shape PC against log-transformed centroid size testing slope homogeneity across taxa. |
| `GradCAM_audit_panel` | `scripts/analysis/05_cleanlab_vision_xai.py` | 3-panel XAI audit showing low-pass filtered crops, DINOv2 spatial patch attention, and calibrated posteriors. |
| `fig05b_cleanlab_quality_and_tomentum` | `scripts/analysis/05_cleanlab_vision_xai.py` | Cleanlab label quality distributions and CIE luminance plotted against HSV saturation to separate hair discordance. |
| `spatial_rf_niche_importance` | `scripts/analysis/06_multimodal_spatial_rf.R` | Spatial RF out-of-bag permutation importance identifying primary environmental predictors under block CV. |
| `fig06_ecological_niche_phenology` | `scripts/analysis/06_multimodal_spatial_rf.R` | 4-panel niche plate displaying climatic versus edaphic overlap, SGMC lithology, and anthesis timing. |
| `fig06b_geographic_map_and_hopkins_phenology` | `scripts/analysis/06_multimodal_spatial_rf.R` | Geographic voucher distribution map colored by species and anthesis day of year against Hopkins bioclimatic law. |
| `fig06b_sdm_habitat_suitability_map` | `scripts/analysis/06_multimodal_spatial_rf.R` | 4-panel MaxEnt and Spatial RF ensemble habitat suitability plate with authority tiers and type specimen overlays. |
| `fig07_multimodal_discordance_summary` | `scripts/analysis/07_triage_dashboard_synthesis.R` | Modality conflict rates weighted by out-of-fold Cohen's kappa and consensus curatorial triage breakdowns. |
| `Figure_Integrative_Packera_Revision` | `scripts/analysis/07_triage_dashboard_synthesis.R` | 4-panel synthesis dashboard uniting 3D CDA morphospace, DINOv2 confidence, phenology, and taxonomic treatments. |
| `fig_empirical_taxon_complex_network` | `scripts/analysis/08_discover_taxon_complexes.R` | Unsupervised multimodal affinity network graph depicting species complexes and morphotype subclusters. |

---

## 6. License and Attribution

This project is licensed under the **MIT License**. Herbarium specimen imagery harvested through the pipeline remains subject to the data access agreements and institutional copyright policies of the contributing herbaria (NCU, NY, US, MO, GA, VDB, and partner collections).
