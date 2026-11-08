# Automated High-Throughput Vegetative Phenotyping and Herbarium Triage Pipeline for North American *Packera*

### Integrative Morphometric and Edaphic Species Delimitation Anchored by Nomenclatural Types

[](https://opensource.org/licenses/MIT) [](https://www.python.org/downloads/) [](https://www.r-project.org/) [](https://github.com/Gene-Weaver/LeafMachine3)

## Project Overview

Species boundaries in North American *Packera* Á. Löve & D. Löve (Asteraceae: Senecioneae) are difficult to resolve from herbarium sheets alone. Across the ***Packera dubia* (Spreng.) Trock & Mabb. complex** and its eastern North American congeners, taxa exhibit pronounced heteroblasty between basal rosette leaves and upper cauline stem leaves alongside frequent polyploidy ($2n = 44, 46, 88$; Kowal 1975) and reticulate hybridization. Many lineages also occupy specialized edaphic habitats on granite flatrocks and diabase glades, alongside serpentine barrens and coastal plain sands. Because herbarium collections contain misidentifications and mixed duplicate gatherings (exsiccatae), pooling raw specimen measurements by folder labels obscures evolutionary boundaries.

I built this pipeline at the University of North Carolina at Chapel Hill and the UNC Herbarium (NCU) to evaluate morphological diagnosability across North American *Packera* under the Unified Species Concept (de Queiroz 2007) through Integration by Congruence (Padial et al. 2010). To account for variation in determiner reliability across herbarium holdings, I stratify specimens by taxonomic authority and pin verified nomenclatural type specimens as morphological and *locus classicus* ecological baselines. I isolate intact basal rosette leaves from cauline stem leaves before quantifying lamina shape, then test whether vegetative leaf geometry and reproductive macro-morphology align with flowering phenology, alongside fine-scale rhizosphere and bedrock niches ($0\text{ to }15\text{ cm}$ soil properties and $10\text{ to }30\text{ m}$ outcrop lithology).

Now that I have framed the taxonomic problem, let us examine how the data move through the system. I structure the pipeline around verified nomenclatural types, as historical primary and secondary types provide immutable morphological anchors. When I evaluate an unknown specimen, I measure its divergence from these type baselines to establish authentic morphological identity.

### Core Empirical Capabilities

* **Taxonomic Authority Stratification and Nomenclatural Type Pinning (Step 01):** I harvest Darwin Core archives from GBIF and iDigBio, resolve basionyms across the genus, deduplicate exsiccatae series prior to image download, store sheets as lossless-calibrated WebP images, intersect coordinates with official US EPA Level III Ecoregions (Omernik & Griffith 2014), and classify determiners across a 4-tier authority hierarchy (`Tier_0_Type`, `Tier_1_Gold` monographers and protologue authors, `Tier_2_Silver` curated regional herbaria, and `Tier_3_Bronze` general collections). The harvester retains historical nomenclatural type gatherings with county-level retro-georeferences or partial collection dates as permanent taxonomic anchors.
* **LeafMachine3 Computer Vision (`scripts/vision/`, `scripts/pipeline/02_segment_and_extract.py`):** I run [LeafMachine3](https://github.com/Gene-Weaver/LeafMachine3) (Weaver 2026) on ONNX Runtime to execute Stage 1 plant and organ detection (`lm3_plant_detector__yolo26x_det_1280`, `conf = 0.30`, nested-box deduplication at `overlap = 0.95`) with dynamic stride-32 letterboxing (`auto=True, stride=32`, `cv2.INTER_LINEAR_EXACT`). I perform whole-sheet substrate segmentation with `lm3_specimen_segmenter` (`UNet++ v9` with `LLR_Hyst_KeepBlobs` paper cleanup), and resolve millimeter scale calibration through `lm3_ruler_classifier` paired with `lm3_ruler_segmenter`. For individual leaves, I run tight Stage 2 per-leaf segmentation (`lm3_leaf_segmenter__yolo26x_seg_1024`, `crop_pad = 0.0`, native `process_mask_native` prototype decoding, unioning all valid `Leaf`, `Petiole`, and `Hole` instances `>= 64 px`), and infer 31-keypoint anatomical pose estimation (`lm3_landmark_detector__yolo26x_pose_640` with `10%` white border padding). I compute $32 \times 32$ radial-and-angular Euler Characteristic Transform (ECT) matrices (`data/ect/*.npy`), and enforce a per-plant (`plant_individual_id`) adaptive physical blade-area floor ($\min(0.75, \max(0.15, 0.25 \times \text{median\_area}))\text{ cm}^2$) without a sheet-wide leaf count cap.
* **Type-Anchored Elliptic Fourier Analysis and Multivariate Morphometrics (Steps 03 to 04):** I decompose closed `Tier_1_Direct` basal leaf outlines using Elliptic Fourier Analysis (`Momocs::efourier`, calibrated to capture $\ge 99\%$ cumulative harmonic power), symmetrize bilateral leaf outlines (`Momocs::rm_asym`; Claude 2008), and aggregate leaves to Operational Taxonomic Unit (OTU) medians per plant individual. In Step 04, I estimate determiner error rates via an exogenous Dawid-Skene and Beta-Binomial Empirical Bayes model, test Klingenberg (2016) allometric slope homogeneity, shrink taxon reference centroids toward nomenclatural type medoids in proportion to effective type sample size (ICN Art. 7.2), and run both unsupervised Gaussian Mixture Modeling (`Morphotype_K1..Kn`) and supervised Canonical Discriminant Analysis and Pillai/Wilks MANOVA across EPA Level III Ecoregions.
* **Fine-Grain Rhizosphere/Outcrop Niche Modeling and SDM Ensembles (Step 05):** I pair 19 WorldClim bioclimatic variables with a cached $0\text{ to }15\text{ cm}$ rhizosphere and bedrock stack combining 50 m-buffered USDA-NRCS SSURGO map units (`ssurgo_rock_outcrop_pct`, `ssurgo_bedrock_depth_cm`), 30 m POLARIS probabilistic soil grids, 100 m SOLUS100 bedrock depth and rock fragment volume, $10\text{ to }30\text{ m}$ USGS 3DEP micro-topography, and USGS SGMC bedrock lithology, alongside circular Hopkins flowering phenology anomalies. I fit spatially blocked Random Forest and authority-weighted MaxEnt (`maxnet`) habitat suitability models, and run partitioned climatic versus edaphic niche equivalency tests (Warren et al. 2008; Broennimann et al. 2012).
* **Reliability-Weighted Consensus and Quantitative Taxon-Complex Discovery (Steps 06 to 07):** In Step 06, I synthesize all lines of evidence through an out-of-fold Cohen's $\kappa$ reliability-weighted consensus query against `v_integrative_triage_master`. I execute a 5-stage hierarchical decision matrix that protects nomenclatural types under Rule A, reassigns misidentified specimens under Rule B, and locks typical morphotypes under Rule C. In Step 07, I discover species complexes through unsupervised graph modularity (Leiden and Louvain) and Ward.D2 clustering across multimodal affinity matrices to delimit species complexes and flag candidate cryptic lineages.

### Project Metadata

* **Principal Investigator:** J. Brandon Fuller (PhD Candidate, Department of Biology, University of North Carolina at Chapel Hill)
* **Faculty Advisor:** Dr. Alan S. Weakley (Director, UNC Herbarium [NCU])
* **Methodological Reference:** See [docs/TAXONOMIC_SOP.md](docs/TAXONOMIC_SOP.md) (`UNC-BOT-SOP-2026-04-REV5`) for complete mathematical equations, statistical proofs, and botanical criteria.

## Seven-Stage Pipeline Architecture

The pipeline executes under a **Two-Tier Canonical Storage Architecture**:
1. **Tier A: Canonical Shared Store (`data/` and `data/packera_pipeline.sqlite`):** Houses invariant Step 01 and Step 02 specimen records, WebP imagery (`data/raw_vouchers/`), organ masks and contours (`data/masks/`, `data/contours/`, `data/cropped_patches/`), and environmental caches (`environmental_cache`, `environmental_rtree`).
2. **Tier B: Run-Scoped Analysis Directories (`runs/<RUN_ID>/` and `runs/latest`):** Registers each execution in `pipeline_runs`, populates `run_active_cohort`, links shared Tier A image and mask assets via symbolic links, and isolates all Step 03 to 07 tables, reports, rasters, and publication figures inside `runs/<RUN_ID>/`.

| Step | Entrypoint & Helper Modules | Primary Inputs | Primary Outputs (SQLite & CSV/Figures) | Core Scientific Function |
| :--- | :--- | :--- | :--- | :--- |
| **01** | `scripts/pipeline/01_voucher_harvester.py`<br>`scripts/core/harvester.py`<br>`scripts/core/taxonomy_authority.py`<br>`scripts/core/epa_ecoregions.py` | GBIF & iDigBio APIs, `config/config.yaml`, US EPA Level III Ecoregion polygons | `vouchers`, `taxon_synonyms`, `data/tables/curated_vouchers.csv`, `data/raw_vouchers/*.webp` | Synonym expansion, pre-download exsiccatae deduplication, 3-layer determiner authority scoring, nomenclatural type retention, and EPA Level III ecoregion assignment |
| **02** | `scripts/pipeline/02_segment_and_extract.py`<br>`scripts/vision/lm3_inference_engines.py`<br>`scripts/vision/lm3_data_loader.py`<br>`scripts/vision/leaf_geometry_utils.py`<br>`scripts/vision/leaf_gatekeeper.py`<br>`scripts/vision/organ_crown_metrics.py`<br>`scripts/vision/step02_figures.py` | `curated_vouchers.csv`, `data/raw_vouchers/*.webp`, LeafMachine3 ONNX models | `plant_individuals`, `extracted_leaves`, `voucher_reproductive_phenology`, `extracted_leaf_manifest.csv`, `fig02a`, `fig02b` | LeafMachine3 ONNX detection (`yolo26x_det_1280`), UNet++ v9 paper cleanup, ruler lattice + MP scale calibration, 3-class `Leaf`/`Petiole`/`Hole` segmentation, 31-keypoint pose tracking, DBSCAN plant clustering, basal rosette selection, rachis bridging, and 120-vertex tip-up contour + 32x32 ECT export |
| **03** | `scripts/morphometrics/03_fourier_extractor.R`<br>`scripts/morphometrics/dissected_leaf_utils.R`<br>`scripts/morphometrics/fourier_parallel_utils.R` | `extracted_leaves.contour_xy_json`, `reports/Leaf_Momocs/momocs_outlines.json`, `Coordinates/*.h5` | `leaf_efa_harmonics`, `leaf_efa_harmonics.csv`, `fourier_pcs.csv`, `fig03` plates | Batch SQLite/Momocs tip-up contour ingestion, $\ge 99\%$ harmonic power calibration, bilateral symmetrization (`rm_asym`), type-anchored MCD outlier screening, and OTU median aggregation |
| **04** | `scripts/morphometrics/04_gmm_morphotools.R`<br>`scripts/morphometrics/authority_weighting_utils.R`<br>`scripts/morphometrics/pub_figures_utils.R`<br>`scripts/morphometrics/pub_ecoregion_manova_utils.R` | `leaf_efa_harmonics`, `curated_vouchers.csv`, `voucher_reproductive_phenology` | `otu_morphometrics`, `morphometric_audit.csv`, `gmm_bayes_factors_summary.csv`, `fig04a`, `fig04b`, `fig04c` | Dawid-Skene + Beta-Binomial Empirical Bayes authority weighting, Klingenberg allometric testing, Art. 7.2 type-medoid shrinkage, Track 1 unsupervised GMM clustering, and Track 2 supervised EPA Level III Ecoregion CDA & MANOVA |
| **05** | `scripts/analysis/05_multimodal_spatial_rf.R`<br>`scripts/analysis/step05_phenology_ecotypes.R`<br>`scripts/analysis/step05_sdm_maxent.R`<br>`scripts/analysis/edaphic_fine_grain_utils.R` | `otu_morphometrics`, WorldClim, SSURGO, POLARIS, SOLUS100, 3DEP, SGMC | `environmental_cache`, `spatial_rf_predictions`, `multimodal_conflict_flags.csv`, `sdm_*.tif`, `fig05`, `fig05b` | $0\text{ to }15\text{ cm}$ rhizosphere and outcrop extraction, *locus classicus* ecological baseline control, Hopkins circular phenology anomalies, Spatial Block-CV Random Forest, MaxEnt (`maxnet`) ensemble SDMs, and niche equivalency tests |
| **06** | `scripts/analysis/06_triage_dashboard_synthesis.R`<br>`scripts/morphometrics/pub_synthesis_plates.R` | `v_integrative_triage_master` (pre-joined SQLite view of Steps 01 to 05) | `curation_triage`, `final_curation_triage.csv`, `Packera_Taxonomic_Revision_Summary.md`, `fig06` | Out-of-fold Cohen's $\kappa$ reliability-weighted consensus, 5-stage hierarchical nomenclatural type decision matrix (Rules A, B, and C), *sensu stricto* type locking, and herbarium curation queue generation |
| **07** | `scripts/analysis/07_discover_taxon_complexes.R`<br>`scripts/morphometrics/pub_synthesis_plates.R` | `v_integrative_triage_master`, `otu_morphometrics`, `spatial_rf_predictions` | `empirical_taxon_complexes.csv`, `taxon_complex_pairwise_affinity.csv`, `fig07` | Unsupervised graph modularity (Leiden/Louvain) and Ward.D2 hierarchical clustering across morphometric, phenological, climatic, and edaphic affinity matrices to delimit species complexes |

```mermaid
flowchart TD
    DB[("Central Transactional SQLite WAL Store\n(data/packera_pipeline.sqlite)\nTier A: vouchers · plant_individuals · extracted_leaves · environmental_cache\nTier B: pipeline_runs · run_active_cohort · otu_morphometrics · spatial_rf_predictions · curation_triage\nView: v_integrative_triage_master")]

    B1["Step 01: GBIF/iDigBio Ingestion, Type Hierarchy and Authority Engine\n(01_voucher_harvester.py · taxonomy_authority.py · epa_ecoregions.py)"] -->|"Atomic UPSERT and WebP Store"| DB
    B1 --> S02["Step 02: LeafMachine3 Segmentation<br/>YOLO26x Det/Seg/Pose + UNet++ v9 + 32x32 ECT<br/>(scripts/pipeline/02_segment_and_extract.py)"]
    S02 --> B3["Step 02: Packera Basal-Rosette Gatekeeper and Dissected Rachis Bridging\n(leaf_gatekeeper.py · organ_crown_metrics.py · 120-Vertex Tip-Up Contours and 32x32 ECT)"]
    B3 -->|"Per-Voucher SQLite Transaction"| DB
    B3 -->|"Tip-Up Tier_1_Direct Outlines"| B4["Step 03: Closed-Contour EFA (>=99% Power), Symmetrization and OTU Medians\n(03_fourier_extractor.R · dissected_leaf_utils.R · Type MCD Anchoring)"]
    B4 -->|"leaf_efa_harmonics"| DB
    B4 --> B5["Step 04: Dawid-Skene Authority Weighting, Type-Medoid Shrinkage and Dual CDA/GMM\n(04_gmm_morphotools.R · authority_weighting_utils.R · pub_ecoregion_manova_utils.R)"]
    B5 -->|"otu_morphometrics"| DB
    B5 --> B6["Step 05: 0 to 15 cm Rhizosphere Stack, Spatial RF and MaxEnt SDMs\n(05_multimodal_spatial_rf.R · step05_phenology_ecotypes.R · step05_sdm_maxent.R)"]
    B6 -->|"environmental_rtree and spatial_rf_predictions"| DB
    B6 --> B7["Step 06: Reliability-Weighted Consensus and Hierarchical Type Decision Matrix\n(06_triage_dashboard_synthesis.R · pub_synthesis_plates.R)"]
    DB -->|"v_integrative_triage_master"| B7
    B7 -->|"curation_triage and final_curation_triage.csv"| DB
    B7 --> B8["Step 07: Unsupervised Quantitative Taxon-Complex Discovery\n(07_discover_taxon_complexes.R · pub_synthesis_plates.R)"]
    DB -->|"v_integrative_triage_master"| B8
    B8 -->|"empirical_taxon_complexes.csv"| DB
```

## Repository Structure

```text
Packera_dubia_morphometrics/
├── config/
│   └── config.yaml                                  # Central paths, LeafMachine3 settings and empirical thresholds
├── docs/
│   └── TAXONOMIC_SOP.md                             # Complete mathematical, statistical and botanical protocol (SOP)
├── .agents/
│   └── skills/
│       ├── pipeline-contract-audit/SKILL.md         # Cross-stage SQLite/CSV contract and CV guardrail verification
│       └── write-in-my-voice/SKILL.md               # Writing DNA v5.2 voice and punctuation enforcement rules
├── AGENTS.md                                        # Compact workspace invariants and module map
├── README.md                                        # Operational guide, setup instructions and CLI reference
├── setup_leafmachine3.sh                            # Automated uv + LeafMachine3 ONNX model installer and doctor check
├── run_pipeline.sh                                  # Master run-scoped shell orchestrator (Steps 01 to 07)
├── main.py                                          # Unified Python CLI runner (check-env, harvest, segment, morphometrics, synthesis, run-all)
├── requirements.txt                                 # Primary Python environment dependencies
├── data/                                            # Tier A: Canonical Shared Store (invariant across runs)
│   ├── packera_pipeline.sqlite                      # Transactional SQLite WAL database (Tier A + Tier B tables and views)
│   ├── raw_vouchers/                                # Full-sheet herbarium imagery (.webp)
│   ├── masks/                                       # Binary leaf silhouettes (tier1_intact/, tier2_reflected/)
│   ├── contours/                                    # Tip-up 120-vertex 2D coordinate contours (.csv, mirrored in SQLite)
│   ├── ect/                                         # 32x32 radial-and-angular ECT matrices (*_ect.npy)
│   ├── Segmented_Leaves/                            # Individual leaf PNG cutouts and overlays
│   ├── cropped_patches/                             # Organ RGB cutouts (basal_leaves_raw/, capitula/, rosettes_dense/)
│   ├── environmental/                               # EPA Level III ecoregions, SSURGO, POLARIS, SOLUS100, 3DEP and SGMC cache
│   └── tables/                                      # Canonical Step 01 to 02 CSV manifests
└── scripts/
    ├── core/
    │   ├── config.py                                # Typed configuration loader and hardware thread/batch resolvers
    │   ├── db.py                                    # Transactional SQLite WAL and DuckDB persistence operations
    │   ├── db_schema.py                             # SQLite/DuckDB DDL schema and dynamic column migration builders
    │   ├── db_utils.R                               # R SQLite/DuckDB batch reader and atomic writer utilities
    │   ├── harvester.py                             # GBIF/iDigBio occurrence harvester, synonym expansion and WebP migrator
    │   ├── taxonomy_authority.py                    # 3-layer determiner authority engine, type hierarchy and exsiccatae deduplication
    │   ├── epa_ecoregions.py                        # US EPA Level III Ecoregion point-in-polygon and county fallback index
    │   ├── artifact_manager.py                      # Checksum verification and remote model artifact acquisition
    │   └── logger.py                                # Standardized console and file logging
    ├── vision/
    │   ├── configure_leafmachine3.py                # Synthesizes LM3_settings.local.yaml for all 18 LeafMachine3 modules
    │   ├── lm3_data_loader.py                       # Working-frame unscaling (work_scale), ruler lattice loader and 31-keypoint reader
    │   ├── lm3_inference_engines.py                 # ONNX Runtime engines (YOLO26x Det/Seg/Pose and UNet++ v9 Specimen Segmenter)
    │   ├── leaf_geometry_utils.py                   # Tip-up contour resampling, Tier 1/Tier 2 extraction and geometry facade
    │   ├── leaf_gatekeeper.py                       # Boundary continuity, dissected rachis bridging and 4-tier leaf routing
    │   ├── organ_crown_metrics.py                   # DBSCAN plant clustering, rosette crown estimation and petiole/scape metrics
    │   ├── capitulum_phenology_classifier.py        # 3-tier reproductive anthesis staging and capitulum macro-metrics
    │   └── step02_figures.py                        # Step 02 QC distribution and taxon silhouette montage plates
    ├── pipeline/
    │   ├── 01_voucher_harvester.py                  # Step 01 CLI entrypoint
    │   └── 02_segment_and_extract.py                # Step 02 CLI entrypoint and wave-pipelined coordinator
    ├── morphometrics/
    │   ├── 03_fourier_extractor.R                   # Step 03 closed-contour EFA, symmetrization and OTU median aggregation
    │   ├── 04_gmm_morphotools.R                     # Step 04 allometry, type-medoid shrinkage, GMM clustering and CDA/MANOVA
    │   ├── authority_weighting_utils.R              # Dawid-Skene and Beta-Binomial Empirical Bayes determiner weighting
    │   ├── dissected_leaf_utils.R                   # Dissection index, harmonic power calibration and regime-stratified MCD
    │   ├── fourier_parallel_utils.R                 # Parallel EFA harmonic extraction workers
    │   ├── pub_figures_utils.R                      # Core ordination, taxon palettes and Step 03 to 04 figure generators
    │   ├── pub_ecoregion_manova_utils.R             # Track 2 supervised EPA Level III Ecoregion CDA and MANOVA audit suite
    │   └── pub_synthesis_plates.R                   # North American basemaps, pure-R GMM EM solver and Steps 05 to 07 plates
    ├── analysis/
    │   ├── 05_multimodal_spatial_rf.R               # Step 05 rhizosphere and outcrop extraction, Spatial Block RF and niche overlap
    │   ├── step05_phenology_ecotypes.R              # Hopkins circular phenology anomalies and locus classicus ecological controls
    │   ├── step05_sdm_maxent.R                      # Authority-weighted MaxEnt and Spatial RF ensemble SDMs and Boyce Index
    │   ├── edaphic_fine_grain_utils.R               # SSURGO SDA, POLARIS, SOLUS100, 3DEP and SGMC lithology extractors
    │   ├── 06_triage_dashboard_synthesis.R          # Step 06 reliability-weighted consensus and 5-stage type decision matrix
    │   └── 07_discover_taxon_complexes.R            # Step 07 Leiden/Louvain and Ward.D2 unsupervised complex discovery
    └── tests/                                       # Comprehensive Python and R unit and contract hygiene test suites
```

## Cross-Platform Environment and Setup (Linux, macOS, and Windows WSL 2)

### 1. Windows WSL 2 Storage and GPU Prerequisites (Windows Users Only)

On Windows 10 or Windows 11, run the pipeline inside **WSL 2 (Ubuntu 22.04 LTS or newer)** and store the repository on the native Linux `ext4` filesystem (`~/Packera_dubia_morphometrics`) for high-speed SQLite WAL concurrency and raster I/O.

Install the Windows NVIDIA driver ($\ge 528.33$ on Windows, or $\ge 525.60$ on native Linux) and verify that `nvidia-smi` prints your GPU table inside WSL 2. [LeafMachine3](https://github.com/Gene-Weaver/LeafMachine3) installs its own pinned CUDA runtime wheels directly through `uv.lock`, so only the host NVIDIA driver is required.

### 2. System Packages and Primary Python Environment (`.venv`)

```bash
# Ubuntu / WSL 2 system libraries for OpenCV, GDAL/PROJ/GEOS (sf/terra), and R compilation
sudo apt-get update && sudo apt-get install -y \
  build-essential git curl wget libgl1 libglib2.0-0 python3-venv python3-dev \
  libgdal-dev gdal-bin libgeos-dev libproj-dev libudunits2-dev libsqlite3-dev \
  r-base r-base-dev

# Create and populate the primary orchestration environment (.venv)
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

### 3. LeafMachine3 ONNX Environment and Model Setup (`setup_leafmachine3.sh`)

Run the automated setup script to install `uv` (`0.12.23`), clone [LeafMachine3](https://github.com/Gene-Weaver/LeafMachine3), sync the hardware-matched `uv.lock` environment (`--extra gpu`, `--extra macos`, or `--extra cpu`), download the locked Hugging Face ONNX models (~1.7 GB), and verify hardware binding with `lm3 doctor`:

```bash
bash setup_leafmachine3.sh
```

### 4. R Statistical Environment and Preflight Verification

```bash
# Install required CRAN packages for Steps 03 to 07
Rscript -e 'install.packages(c("Momocs", "RRPP", "vegan", "MASS", "cluster", "mclust", "ranger", "maxnet", "ecospat", "sf", "terra", "DBI", "RSQLite", "duckdb", "dplyr", "readr", "tibble", "ggplot2", "optparse", "jsonlite", "igraph"), repos="https://cloud.r-project.org")'

# Run the unified environment preflight check
python main.py check-env
```

## Quick-Start Terminal Cheat-Sheet

```bash
# 0. Run pipeline skipping Step 01 network harvest (seeds from existing SQLite / CSV vouchers)
./run_pipeline.sh --skip-harvest

# 1. Default Eastern North American genus-wide run
./run_pipeline.sh

# 2. Targeted taxon subset mode (harvests and analyzes the specified Packera species)
./run_pipeline.sh --taxon "Packera dubia" "Packera aurea" "Packera paupercula" "Packera glabella"

# 3. Custom named run with explicit CPU worker threads and GPU batch size
./run_pipeline.sh --run-id packera_eastern_rev1 --threads 12 --gpu-batch-size 16

# 4. Resume an existing run from a downstream step (seeds active cohort from data/)
./run_pipeline.sh --run-id packera_eastern_rev1 --start-step 3

# 5. Strategy 1: Focal Taxon + Cascading Outgroup Reference Anchors (default --gold-anchors 30)
./run_pipeline.sh --focal-taxon "Packera dubia" --gold-anchors 30

# 6. Strategy 2: Intraspecific Delimitation Mode (single-species morphotype and ecotype evaluation)
./run_pipeline.sh --focal-taxon "Packera dubia" --intraspecific-mode

# 7. Continent-wide North American expansion across the US, Canada, and Mexico
./run_pipeline.sh --geographic-scope north_america --countries "US,CA,MX"

# 9. Execute standalone Step 07 unsupervised complex discovery
Rscript scripts/analysis/07_discover_taxon_complexes.R \
  --morphometrics runs/latest/otu_morphometrics.csv \
  --harmonics runs/latest/leaf_efa_harmonics.csv \
  --vouchers runs/latest/curated_vouchers.csv \
  --multimodal runs/latest/multimodal_rf_predictions.csv \
  --output-complexes runs/latest/discovered_taxon_complexes.csv \
  --output-affinity runs/latest/pairwise_taxon_affinity.csv \
  --output-neighbors runs/latest/recommended_complex_neighbors.txt \
  --output-figure runs/latest/figures/fig07_taxon_complex_network.png \
  --focal-taxon "Packera dubia"
```

## CLI Feature Flags Reference (`./run_pipeline.sh` and `python main.py`)

### Table 4A: Master Shell Orchestrator Flags (`./run_pipeline.sh`)

| Flag | Arguments | Default | Description |
| :--- | :--- | :--- | :--- |
| `--run-id` | `<ID>` | `<YYYYMMDD_HHMMSS>_<scope>` | Assigns a custom run identifier under `runs/<ID>/`. |
| `--taxon` | `"<Name>"` | Genus-wide | Restricts pipeline execution to a single target taxon (triggers single-taxon mode unless combined with `--taxa` or `--auto-neighbors`). |
| `--focal-taxon` | `"<Name>"` | `"Packera dubia"` | Sets the focal reference taxon for supervised anchors, morphometric projections, and complex discovery. |
| `--taxa` | `"<T1,T2,...>"` | None | Restricts pipeline execution to a comma-separated list of target taxa. |
| `--complex` | `<preset>` | None | Expands `--taxa` to a predefined taxonomic complex (`aurea`, `paupercula`, `tomentosa`, `anonyma`, `obovata`, `schweinitziana`, or `eastern-all`). |
| `--auto-neighbors` | *(flag)* | `false` | Reads `runs/latest/recommended_complex_neighbors.txt` from a prior Step 07 run and appends high-affinity neighbor taxa ($\ge 0.20$) to the focal cohort. |
| `--skip-harvest` | *(flag)* | `false` | Skips Step 01 GBIF/iDigBio network harvesting and seeds the active run from existing SQLite `vouchers` and `data/curated_vouchers.csv`. |
| `--include-western` | *(flag)* | `false` | Retains western North American occurrences west of $-106.65^\circ\text{W}$ by forwarding `--no-exclude-western` to Step 01. |
| `--resegment` | *(flag)* | `false` | Forces Step 02 to clear cached segmentation entries for the active cohort and re-run LeafMachine3 ONNX inference (`--force-reprocess`). |
| `--gold-only` | *(flag)* | `false` | Restricts the active run cohort to Nomenclatural Types (`Tier_0_Type_Anchor`) and `Tier_1_Gold` authority vouchers. |
| `--max-per-taxon` | `<N>` | None | Caps the number of harvested vouchers per taxon in Step 01 (`0` disables the per-taxon cap). |
| `--workers` | `<N>` | Auto-detected | Sets the worker thread count for Step 01 image downloads, Step 02 ONNX inference, and R `future::plan(multisession)` parallelism. |
| `--batch-size` | `<N>` | Auto-detected | Overrides the Step 02 image batch size (`0` auto-tunes from available system RAM and GPU VRAM). |
| `--strict-env-gate` | *(flag)* | `false` | Enforces the 95% complete-coverage gate across all 6 environmental layers (`WorldClim`, `3DEP`, `SSURGO`, `POLARIS`, `SOLUS100`, `SGMC`) before Step 05. |
| `--skip-step07` | *(flag)* | `false` | Skips Step 07 unsupervised taxonomic complex discovery. |
| `--verify` | *(flag)* | `false` | Runs the contract verification suite (`pytest scripts/tests/ -q`) after pipeline completion. |
| `--dry-run` | *(flag)* | `false` | Prints planned execution stages and resolved cohort parameters without executing commands. |
| `--demo` | *(flag)* | `false` | Runs synthetic demonstration mode (`PACKERA_DEMO_MODE=1`) isolated inside `demo_ephemeral.sqlite` without touching canonical `data/` tables. |
| `--help`, `-h` | *(flag)* | `false` | Displays the `run_pipeline.sh` command-line usage reference. |

### Table 4B: Python Orchestrator (`python main.py`) Subcommands and Stage-Specific Flags

| Subcommand or Stage | Target Script | Key Flags | Default | Behavioral Scope |
| :--- | :--- | :--- | :--- | :--- |
| `check-env` | `main.py` | `[--strict]` | `false` | Audits Python and R virtual environments, system libraries, and ONNX Runtime execution providers. |
| `download-weights` | `scripts/core/artifact_manager.py` | `[--force]` | `false` | Downloads and verifies SHA-256 checksums for LeafMachine3 ONNX model weights from Hugging Face (`brandonfuller/lm3-weights`). |
| `harvest` | `scripts/pipeline/01_voucher_harvester.py` | `--taxon`, `--focal-taxon`<br>`--home-herbarium`<br>`--max-records`<br>`--max-per-taxon`<br>`--max-uncertainty`<br>`--exclude-western` / `--no-exclude-western`<br>`--min-megapixels`<br>`--min-file-size-kb`<br>`--check-sharpness`<br>`--min-sharpness`<br>`--concurrency`<br>`--download-images` | `NCU`<br>`5000`<br>None<br>`5000.0` m<br>`--exclude-western`<br>`2.0` MP<br>`150.0` KB<br>`false`<br>`25.0`<br>`15`<br>`true` | Queries GBIF and iDigBio Darwin Core APIs, standardizes taxon names, resolves nomenclatural type tiers, and downloads calibrated WebP voucher images. |
| `segment` | `scripts/pipeline/02_segment_and_extract.py` | `--enable-geometric-qc`<br><br>`--force-reprocess`<br>`--regenerate-figures-only`<br><br>`--lm3-dir <PATH>`<br>`--batch-size <N>`<br>`--workers <N>`<br>`--device {cuda,cpu}`<br>`--min-solidity <F>`<br>`--min-ucs <F>`<br>`--score-thresh <F>` | `false`<br><br>`false`<br>`false`<br><br>Auto<br>Auto<br>Auto<br>`cuda`<br>`0.72`<br>`0.85`<br>`0.30` | Executes LeafMachine3 ONNX segmentation, DBSCAN plant clustering, and 31-keypoint anatomical pose tracking. `--enable-geometric-qc` enables the *Packera* botanical gatekeeper, dissected-rachis bridging, and per-plant adaptive area floor (enabled by default in `run_pipeline.sh`, off by default in standalone `02_segment_and_extract.py` for unconstrained LeafMachine3 extraction). `--regenerate-figures-only` rebuilds Step 02 diagnostic plates from cached Step 02 to 04 tables without re-running ONNX segmentation. |
| `morphometrics` | `scripts/morphometrics/03_fourier_extractor.R`<br>`scripts/morphometrics/04_gmm_morphotools.R` | `--harmonics <N>`<br>`--num-pcs <N>`<br>`--max-k <N>`<br>`--workers <N>`<br>`--export-figures` / `--no-export-figures` | `12`<br>`5`<br>`9`<br>Auto<br>`--export-figures` | Runs closed-contour Elliptic Fourier Analysis on `Tier_1_Direct` basal leaves, Claude bilateral symmetrization, Klingenberg allometry testing, type-medoid shrinkage, and Track 1 GMM and Track 2 CDA/MANOVA. |
| `synthesis` | `scripts/analysis/05_multimodal_spatial_rf.R`<br>`scripts/analysis/06_triage_dashboard_synthesis.R` | `--ssurgo-batch-size <N>`<br>`--permutations <N>`<br>`--n-trees <N>`<br>`--workers <N>`<br>`--export-figures` / `--no-export-figures` | `50`<br>`100`<br>`500`<br>Auto<br>`--export-figures` | Extracts 0 to 15 cm rhizosphere and bedrock properties, fits Spatial Block Random Forest and MaxEnt SDMs, and executes Cohen's $\kappa$ reliability-weighted consensus triage across Rules A through C. |
| `run-all` | Pipeline Coordinator | Forwards flags to all pipeline stages | | Executes complete multi-stage pipeline from preflight checks to publication synthesis. |

## Publication Figure Suite Inventory

Every figure is exported atomically in both vector `.pdf` and 300-DPI `.png` formats inside `runs/<RUN_ID>/outputs/figures/`:

| Figure Stem | Pipeline Step | Visual Contents and Statistical Annotations |
| :--- | :---: | :--- |
| `fig02a_segmentation_qc_distributions` | Step 02 | 4-panel distribution of lamina solidity, UCS contour completeness, calibrated blade area ($\text{cm}^2$), and tier routing |
| `fig02b_silhouette_montage` | Step 02 / Post-Step 04 | Downstream-QC-validated `Tier_1_Direct` binary leaf silhouettes scaled to each taxon's empirical median aspect ratio |
| `tier1_vs_tier2_density_overlay` | Step 03 | 2D EFA morphospace density overlay and PERMANOVA test comparing `Tier_1_Direct` intact vs. `Tier_2_Reflected` leaves |
| `fig03_fourier_morphospace_3panel` | Step 03 | Symmetric EFA PCA ordination, mean $\pm 1\text{ SD}$ leaf outline reconstructions along PC1 to PC3, and harmonic power decay |
| `fig04a_allometry_diagnostic` | Step 04 | Log centroid size vs. shape regression testing Klingenberg (2016) allometric slope homogeneity via weighted `RRPP` |
| `fig04b_cda_4panel_plate` | Step 04 | 4-panel integrative plate: Track 1 shrunk CDA biplot (with double-ringed diamond Type anchors), GMM $\times$ Species contingency heatmap, diagnostic trait boxplots, and STRUCTURE-style posterior assignment barplot |
| `fig04c_supervised_ecoregion_cda` | Step 04 | Track 2 supervised US EPA Level III Ecoregion CDA biplot, Pillai's Trace MANOVA subtitle, GMM $\times$ Ecoregion heatmap, and ecoregion mean leaf silhouettes |
| `fig05_ecological_niche_phenology` | Step 05 | 4-panel ecological niche plate: Spatial Block RF permutation importance, `PCA-env` climate and edaphic biplot with Schoener's $D$, $0\text{ to }15\text{ cm}$ rhizosphere bedrock boxplots, and circular Hopkins flowering phenology |
| `fig05b_sdm_habitat_suitability_map` | Step 05 | Multi-panel North American MaxEnt (`maxnet`) and Spatial RF ensemble habitat suitability maps with $AUC_{\text{block}}$, Continuous Boyce Index, and Nomenclatural Type locality overlays |
| `fig06_integrative_synthesis_plate` | Step 06 | Cross-modal concordance matrix, Cohen's $\kappa$ reliability weights, and 5-stage taxonomic decision breakdown across active taxa |
| `fig07_complex_discovery_plate` | Step 07 | Leiden and Louvain morpho-ecological affinity network and Ward.D2 hierarchical dendrogram delimiting empirical species complexes |

## Automated Test Suite Verification

```bash
# Run the complete unit, SQLite handshake, R morphometrics, and static code-hygiene test suite
pytest scripts/tests/ -q
```

## License and Attribution

Released under the MIT License ([LeafMachine3](https://github.com/Gene-Weaver/LeafMachine3) is distributed under the GNU Affero General Public License v3.0). When citing this pipeline or its underlying methods, please cite [docs/TAXONOMIC_SOP.md](docs/TAXONOMIC_SOP.md), Weaver & Smith (2023), Weaver (2026), and Ayub et al. (2026).
