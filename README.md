# Automated High-Throughput Vegetative Phenotyping and Herbarium Triage Pipeline for the *Packera dubia* Complex

### Operationalizing Morphological Diagnosability under the Unified Species Concept

[![DOI]()](https://doi.org/10.5281/zenodo.xxxxxx) [![License: MIT]()](https://opensource.org/licenses/MIT) [![Python 3.10+]()](https://www.python.org/downloads/) [![R 4.3+]()](https://www.r-project.org/)

---

## 🌿 Project Overview

This repository houses the computational and statistical pipeline for the automated high-throughput vegetative phenotyping, morphological diagnosability assessment, and herbarium triage of the ***Packera dubia*** **(Spreng.) Trock & Mabb. complex** (Asteraceae: Senecioneae) across Eastern and Central North America.

Developed as part of doctoral research at the **University of North Carolina at Chapel Hill** in collaboration with the **UNC Herbarium (NCU)**, this project operationalizes morphological diagnosability under the **Unified Species Concept (USC; de Queiroz 2007)**. Under the USC, species are conceptualized as separately evolving metapopulation lineages. In taxonomically complex plant clades characterized by ecophenotypic plasticity and historical misidentification, vegetative morphometrics serves as a high-throughput operational screen evaluated in concert with multiple independent lines of evolutionary evidence:
- **Reproductive Macro-Morphology:** Capitulum aspect ratios, phyllary series counts, and ray/disc floret dimensions.
- **Cytotaxonomic Cytotypes:** Known chromosome numbers and ploidy races ($2n = 44, 46, 88$; Kowal 1975).
- **Reduced-Representation Phylogenomics:** High-throughput nextRAD single nucleotide polymorphism (SNP) datasets resolving reticulate evolution and lineage boundaries.
- **Micro-Edaphic Validation:** Regional 250m SoilGrids pedological rasters complemented by fine-scale 1:24,000 USDA NRCS SSURGO vector map units resolving localized rock outcrop endemics.

The computational infrastructure couples automated organ detection via **LeafMachine2 (LM2)**, deterministic geometric gatekeeping, bilateral symmetry reconstruction, label-blind Elliptic Fourier Analysis (EFA), self-supervised deep vision (DINOv2), and confident learning with a formal **Multi-Tiered Herbarium Misidentification Mitigation Architecture**.

- **Principal Investigator:** J. Brandon Fuller (PhD Candidate, Department of Biology, UNC-CH)  
- **Faculty Advisor:** Dr. Alan S. Weakley (Director, NCU Herbarium; UNC Biology)  
- **Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV4`  
- **Target Taxa:** *Packera dubia* (formerly *P. tomentosa* / *Senecio tomentosus*), *Packera anonyma*, *Packera plattensis*, *Packera paupercula* (including var. *paupercula* and var. *savannarum*), and allied southeastern lineages.

---

## ⚙️ Architecture & Two-Track Workflow

The pipeline decouples the **model training / assisted-annotation track** from the **high-throughput production inference track**:

1. **Assisted Annotation & LM2 Fine-Tuning (Offline Track):**  
   - Uses **Segment Anything Model 2 (SAM 2\)** to interactively annotate ground-truth basal leaf masks on a representative set of 50–100 challenging *Packera* vouchers (dense rosettes, overlapping blades, tomentum).  
   - Exports standardized COCO JSON annotations to fine-tune LeafMachine2’s Plant Component Detector (PCD / PointRend backbone), producing domain-adapted weights (`models/lm2_packera_pcd_finetuned.pth`).  
2. **Production Pipeline & Taxonomic Synthesis (Inference Track):**  
   - Orchestrated end-to-end via a unified CLI (`main.py`).  
   - Runs voucher harvesting, high-throughput LM2 extraction, midrib symmetry reconstruction, 12-harmonic normalized EFA in R (`Momocs`), unsupervised GMM clustering (`mclust`), passive sample CDA (`MorphoTools2`), and environmental niche modeling.

flowchart TD

    subgraph TrackA\["Track A: Assisted Annotation & Domain Adaptation (Run Once)"\]

        A1\["Raw Vouchers (50–100 Challenging Sheets)"\] \--\> A2\["Interactive Mask Prompting\\n(annotate\_with\_sam2.py)"\]

        A2 \--\> A3\["COCO Dataset Export\\n(sam2\_annotator\_utils.py)"\]

        A3 \--\> A4\["Fine-Tune LM2 Plant Component Detector\\n(finetune\_lm2\_pcd.py)"\]

        A4 \--\> A5\[("Domain-Adapted Weights:\\nlm2\_packera\_pcd\_finetuned.pth")\]

    end

    subgraph TrackB\["Track B: Unified Production Pipeline (main.py)"\]

        B1\["Phase 1: Ingestion\\n(python main.py harvest)"\] \--\> B2\["GBIF/iDigBio Ingestion & Determiner Authority Tiers\\n(Gold Specialists | Silver Herbaria | Bronze Unverified)"\]

        B2 \--\> B3\["Phase 2: Segmentation & Extraction\\n(python main.py segment)"\]

        A5 \-.-\> B3

        B3 \--\> B4{"Geometric Routing"}

        B4 \--\>|"Tier 1 (Pristine: UCS ≥ 0.85)"| B5\["Closed Silhouettes"\]

        B4 \--\>|"Tier 2 (Hemi-Blade)"| B6\["Midrib Bilateral Reflection"\]

        B6 \--\> B5

        B4 \--\>|"Failed / Clumped"| B7\["Excluded / Whole-Rosette Embeddings"\]

        B5 \--\> B8\["Phase 3: Morphometrics\\n(python main.py morphometrics)"\]

        B8 \--\> B9\["12-Harmonic Normalized EFA\\n(03\_fourier\_extractor.R · Momocs)"\]

        B9 \--\> B10\["GMM Cluster Testing & Passive CDA\\n(04\_gmm\_morphotools.R · Gold Anchors, Bronze Passive)"\]

        B10 \--\> B11\["Phase 4: Synthesis\\n(python main.py synthesis)"\]

        B11 \--\> B12\["Multi-Modal Niche Models & Digital Triage Queue\\n(SoilGrids 250m, WorldClim, triage\_queue.csv)"\]

    end

---

## 📁 Streamlined Repository Structure

Packera\_dubia\_morphometrics/

├── config/

│   └── config.yaml                    \# Central paths, taxon synonyms & model hyperparameters

├── main.py                            \# Unified CLI runner (harvest, segment, morphometrics, run-all)

├── LeafMachine2/                      \# LeafMachine2 engine submodule

│   ├── LeafMachine2.py

│   └── leafmachine2/

├── models/

│   └── lm2\_packera\_pcd\_finetuned.pth  \# Domain-adapted LM2 weights

├── data/

│   ├── raw\_vouchers/                  \# High-resolution specimen sheet imagery (.jpg)

│   ├── masks/                         \# Binary leaf silhouette masks

│   ├── contours/                      \# Extracted 2D coordinate contours (.csv)

│   ├── environmental/                 \# SoilGrids 250m and WorldClim 2.1 rasters

│   └── tables/

│       ├── curated\_vouchers.csv       \# Harvested vouchers with Determiner Authority Tiers

│       ├── leaf\_efa\_harmonics.csv     \# 12-harmonic EFA coefficients & PC scores

│       ├── morphometrics\_flags.csv    \# CDA predictions & GMM cluster assignments

│       └── triage\_queue.csv           \# Final ranked expert verification queue

├── scripts/

│   ├── core/                          \# Consolidated shared utilities

│   │   ├── harvester.py               \# Combined API queries, DwC parsing & media downloads

│   │   ├── config.py                  \# Python-side configuration schema

│   │   └── logger.py                  \# Structured console and file logger

│   ├── annotation\_and\_training/       \# Track A: SAM 2 labeling & LM2 fine-tuning

│   │   ├── annotate\_with\_sam2.py      \# Interactive point/box prompt segmentation tool

│   │   ├── sam2\_annotator\_utils.py    \# Geometry helpers & COCO export utilities

│   │   └── finetune\_lm2\_pcd.py        \# LM2 Plant Component Detector fine-tuning script

│   ├── pipeline/                      \# Track B: Production execution modules

│   │   ├── 01\_voucher\_harvester.py    \# Ingestion & 3-tier authority stratification

│   │   └── 02\_segment\_and\_extract.py  \# LM2 inference, midrib reflection & contour export

│   ├── morphometrics/                 \# Canonical R statistical morphometrics

│   │   ├── 03\_fourier\_extractor.R     \# 12-harmonic normalized EFA & PCA (Momocs)

│   │   └── 04\_gmm\_morphotools.R       \# mclust GMMs & passive sample CDA (MorphoTools2)

│   ├── analysis/                      \# Macroecology and synthesis

│   │   ├── 05\_cleanlab\_vision\_xai.py  \# DINOv2 embeddings & Grad-CAM XAI

│   │   ├── 06\_multimodal\_spatial\_rf.R \# Spatial Random Forests (SoilGrids, WorldClim)

│   │   └── 07\_triage\_dashboard.R      \# Decision matrix & triage queue generation

│   ├── tests/                         \# Unit and integration test suite

│   │   ├── test\_voucher\_harvester.py

│   │   ├── test\_segment\_and\_extract.py

│   │   └── test\_morphometrics.R

│   └── \_archive/                      \# Retired duplicate scripts (e.g., Python Fourier)

│       └── 03\_fourier\_extractor.py

├── outputs/

│   ├── figures/                       \# Publication plots, CDA biplots & QC panels

│   └── reports/                       \# Taxonomic summaries, BIC tables & manifests

├── requirements.txt                   \# Primary Python dependencies

└── README.md

---

## 🚀 Execution Guide

### Track A: Assisted Annotation & LM2 Fine-Tuning (Run Once)

To fine-tune LeafMachine2's PointRend Plant Component Detector (PCD) on dense, overlapping *Packera* basal rosettes:

```bash
# 1. Interactively annotate 50–100 vouchers with SAM 2
python scripts/annotation_and_training/annotate_with_sam2.py \
    --images-dir data/raw_vouchers/ \
    --output-coco data/annotations/packera_train_coco.json

# 2. Fine-tune LM2's Plant Component Detector (PCD) on generated COCO annotations
python scripts/annotation_and_training/finetune_lm2_pcd.py \
    --coco-annotations data/annotations/packera_train_coco.json \
    --epochs 30 \
    --output-weights models/lm2_packera_pcd_finetuned.pth
```

#### SAM 2 Botanical Annotator Hotkey Cheat-Sheet
| Key / Gesture | Function & Botanical Context |
| :--- | :--- |
| **Mouse Wheel** | Cursor-centered smooth zoom in / out ($1.0\times$–$16.0\times$). |
| **Middle Drag** / `Space`+Drag | Pan viewport smoothly across high-resolution herbarium scan. |
| **`f`** | Fit viewport to window ($1.0\times$). |
| **Left / Right Click** | Positive foreground / negative exclusion prompt points. |
| **Shift + Left-Drag** | Bounding box prompt constraint. |
| **`Tab`** | Cycle 3 SAM 2 candidate granularities (sub-lobe vs. blade vs. clump). |
| **`o`** | Toggle mask fill vs. 1-px boundary contour line (inspect margin crenations). |
| **Hold `v`** | Hold-to-peek: Temporarily hide all overlays to inspect bare pixels. |
| **`[` / `]`** | Adjust mask overlay alpha transparency ($0.10$ to $0.90$). |
| **`k`** | Two-click knife tool to sever petiole bases from caudex tissue. |
| **`+` / `-`** | 1-pixel binary dilation / erosion for tomentum margin tuning. |
| **`0`–`6`** | Instant botanical class commit (`0`: blade, `1`: petiole, `2`: cauline leaf, `3`: cauline stem, `4`: root, `5`: rosette clump, `6`: capitulum). |
| **`u` / `c`** | Undo last committed instance (`u`) / clear active prompts (`c`). |
| **`n` / `Enter`** | Advance/skip voucher (`n`) / save annotations and advance to next sheet (`Enter`). |

> Detailed protocols, class ontologies, and boundary best practices for dense rosettes are documented in [`docs/SAM2_Precision_Botanical_Annotation_Guide.txt`](file:///home/brandon/Packera_dubia_morphometrics/docs/SAM2_Precision_Botanical_Annotation_Guide.txt) and [`docs/WORKFLOW_GUIDE.md`](file:///home/brandon/Packera_dubia_morphometrics/docs/WORKFLOW_GUIDE.md).

---

### Track B: Unified Production Pipeline (`main.py`)

All production pipeline phases are orchestrated through the unified entry point `main.py`. Execution parameters and file paths default to the centralized settings in `config/config.yaml` and can be overridden via CLI arguments:

```bash
# Preflight Environment & Dependency Diagnostic:
python main.py check-env

# Pre-cache / Download Fine-Tuned Model Weights (Automatic on First Run):
python main.py download-weights

# Phase 1: Voucher Harvesting & Determiner Authority Stratification
python main.py harvest --max-records 5000 --download-images

# Phase 2: Segmentation, Midrib Symmetry Reflection & Contour Extraction (GPU-accelerated)
python main.py segment --weights models/lm2_packera_pcd_finetuned.pth

# Phase 3: Morphometrics (12-Harmonic Normalized EFA, GMM, Passive CDA in R)
python main.py morphometrics --harmonics 12

# Run the entire pipeline end-to-end:
python main.py run-all --download-images
```

#### Model Weights & Automated Artifact Acquisition
- **Automatic Retrieval:** Fine-tuned LeafMachine2 PointRend weights (`models/lm2_packera_pcd_finetuned.pth`, ~425 MB) are automatically fetched on the first run of `python main.py segment` (or `main.py run-all`) if not locally present.
- **Manual Pre-Caching:** Weights can be pre-cached or force-refreshed ahead of batch execution using the dedicated CLI command:
  ```bash
  python main.py download-weights [--force]
  ```
- **Persistent Archival:** Model weights are permanently deposited under the release / Zenodo archive (`https://github.com/Plantatsrophe/Packera_dubia_morphometrics/releases/download/v1.0-weights/lm2_packera_pcd_finetuned.pth` / Zenodo DOI: `10.5281/zenodo.xxxxxx`).
- **Integrity Verification:** Downloads stream atomically to temporary files and undergo automated SHA256 checksum verification (`config.models.pcd_weights_sha256`) to guarantee artifact integrity and prevent partial corrupted checkpoints.

#### Preflight Sanity Checks & Defensive Gating
`main.py` incorporates automated preflight sanity checks and a dedicated environment diagnostic subcommand:
- **`check-env` Subcommand**: Instantaneous (<5s) validation across Python version ($\ge 3.10$), active CUDA device & VRAM (without tensor allocation), core Python dependencies (`cv2`, `torch`, `cleanlab`, `pygbif`, `yaml`), `Rscript` availability and R package probe (`Momocs`, `mclust`, `MorphoTools2`, `spatialRF`, `terra`, `tidyverse`, `optparse`), directory write permissions, `.venv_LM2`, and fine-tuned model checkpoint.
- **GPU Availability**: Inspects CUDA availability before Phase 2 segmentation, warning if running on CPU or failing early if CUDA was explicitly required.
- **Input Gating**: Verifies existence and non-emptiness of upstream prerequisite files (`curated_vouchers.csv`, model checkpoints, `data/contours/`), failing early with informative diagnostic hints if an upstream step was skipped.
- **Automated Multi-Environment Orchestration**: `main.py` dynamically resolves the dedicated LeafMachine2 Python interpreter (`.venv_LM2/bin/python` on Linux/macOS or `.venv_LM2/Scripts/python.exe` on Windows) for Phase 2 segmentation tasks. Users run `main.py` from the primary pipeline environment without manual virtualenv activation/deactivation; `main.py` seamlessly invokes `02_segment_and_extract.py` under `.venv_LM2` via subprocess, providing real-time diagnostics (`[INFO] Executing LeafMachine2 segmentation via: .venv_LM2/bin/python`) and strict exit-code propagation.
- **Environment Tools**: Checks for system dependencies (such as `Rscript` for Phase 3 morphometrics) before starting processing.

---

## 🔬 Methodological Highlights

1. **Determiner Authority Stratification:** Partitions specimens into three confidence tiers:  
   - 🥇 **Tier 1 (Gold Standard):** Specialist determinations (Barkley, Trock, Kowal, Weakley, Bain, Mahoney, Fuller) and verified type specimens.  
   - 🥈 **Tier 2 (Silver Standard):** Vouchers determined at major herbaria (NCU, GA, US, NY, BRIT, MO, WIS, VDB, FLAS) with complete locality data.  
   - 🥉 **Tier 3 (Bronze Standard):** General floristic collections and unverified aggregator determinations.  
2. **Bilateral Symmetry Reconstruction & Symmetric EFA:** Asteraceae basal leaves frequently overlap in herbarium presses. When a leaf has an unobstructed half-blade from apex to base, the algorithm detects the primary midrib axis and symmetrically reflects the clean half across the midrib in OpenCV, synthesizing a complete, unoccluded bilateral silhouette for closed EFA. Downstream Fourier analysis isolates the **symmetric harmonic component** ($A_n, D_n$), eliminating fluctuating asymmetry artifacts and reflection bias between Tier 1 pristine and Tier 2 reflected leaves.  
3. **Decoupled Metric Scaling:** Closed Elliptic Fourier Analysis is inherently scale-invariant when normalized (`Momocs::efourier(..., norm = TRUE)`). Contour extraction and shape analysis proceed unconditionally; ruler scale detection is decoupled and used strictly for absolute scalar metrics (blade area, petiole length, centroid size).  
4. **Passive Sample Canonical Discriminant Analysis:** Canonical axes in `MorphoTools2` are trained strictly on verified Tier 1 Gold specimens. Tier 3 Bronze vouchers and conflicting specimens are projected passively, preventing aggregator label noise from distorting morphological taxon boundaries.  
5. **Micro-Edaphic Validation:** Regional 250m SoilGrids pedological rasters (pH, CEC, sand fraction, bulk density) are complemented by fine-scale 1:24,000 USDA NRCS SSURGO vector map units, resolving localized edaphic specialization for rock outcrop endemics (granite flatrocks, sandstone glades, and ultramafic barrens).  
6. **Batch-Effect Controls:** Rosette crops undergo mounting paper background neutralization to eliminate herbarium sheet aging and color artifacts in DINOv2 self-supervised embeddings. Institutional ANOVA audits across contributing herbaria verify that morphometric and latent vision clusters represent genuine biological lineages rather than digitization artifacts.  
7. **Allometry-Free Shape Analysis:** Testing multivariate shape coordinates (EFA harmonics) against log-transformed Centroid Size ($\log(CS)$) to prevent plant stature, environmental vigor, or developmental stage from confounding taxonomic clusters. When allometric scaling is detected ($R^2 \ge 0.10$), size-dependent variation is regressed out, ensuring that downstream GMM clustering and passive CDA reflect genuine lineage divergence rather than phenotypic plasticity or leaf size variation (Klingenberg 2016).  
8. **Latitude-Adjusted Flowering Anomalies:** Testing temporal reproductive isolation independent of continental latitudinal clines. Clinal spring progression advances northward at ~4.0 days per degree of latitude (Hopkins' Bioclimatic Law). Calculating phenological anomalies ($\Delta\text{DOY} = \text{DOY}_{\text{obs}} - \text{DOY}_{\text{expected}}$) eliminates latitudinal gradients, enabling unconfounded testing of allochronic speciation and prezygotic isolation across the species complex (Davis et al. 2015).

---

## 💾 Data & Model Availability

### Academic Citation & Replication Statement
For dissertation methods and publication reproduction:
> *"Model weights, extracted 2D contour coordinate matrices, and environmental rasters are persistently deposited on Zenodo under DOI: 10.5281/zenodo.xxxxxx."*

- **Fine-Tuned Model Weights (`models/`):** LeafMachine2 PointRend weights (`lm2_packera_pcd_finetuned.pth`) are deposited on Zenodo (DOI: `10.5281/zenodo.xxxxxx`) and mirrored via GitHub Releases (`v1.0-weights`).
- **Extracted 2D Contour Matrices (`data/contours/`):** Normalized 2D coordinate matrices and 12-harmonic Fourier coefficients (`data/tables/leaf_efa_harmonics.csv`) are persistently archived on Zenodo.
- **Environmental Rasters (`data/environmental/`):** Pre-computed WorldClim 2.1 (30s bioclimatic variables) and SoilGrids 250m GeoTIFF layers are packaged in the Zenodo archive for headless execution of spatial random forests and niche identity tests.

---

## 📦 System Requirements & Dependencies

- **Python $\ge$ 3.10:** `torch`, `torchvision`, `opencv-python`, `scikit-learn`, `scipy`, `pandas`, `numpy`, `pygbif`, `requests`, `pyyaml`.  
- **R $\ge$ 4.3:** `Momocs`, `MorphoTools2`, `mclust`, `spatialRF`, `terra`, `tidyverse`, `optparse`.

> [!TIP]
> **Windows Users**: The recommended and fully supported way to run this pipeline on Windows is via **WSL 2 (Ubuntu)**. See [`docs/WINDOWS_WSL2_SETUP.md`](docs/WINDOWS_WSL2_SETUP.md) for a 10-minute setup guide with full GPU and GUI support.

---

## 📚 Key Literature & Citations

1. Barkley, T. M. 1988. Variation among the Senecioneae (Asteraceae) in North America. *Brittonia* 40(2): 211–221. doi: 10.2307/2807005  
2. Davis, C. C., C. G. Willis, B. Connolly, C. Kelly, and A. M. Ellison. 2015. Herbarium records are a viable alternative to traditional phenological data for testing responses to climate change. *Journal of Ecology* 103(5): 1126–1134. doi: 10.1111/1365-2745.12457  
3. de Queiroz, K. 2007. Species concepts and species delimitation. *Systematic Biology* 56(6): 879–886. doi: 10.1080/10635150701701083  
4. Klingenberg, C. P. 2016. Size, shape, and form: concepts of allometry in geometric morphometrics. *Development Genes and Evolution* 226(3): 113–137. doi: 10.1007/s00427-016-0539-8  
5. Kowal, R. R. 1975. Systematics of *Senecio aureus* and allied species on the Gaspé Peninsula, Quebec. *Memoirs of the Torrey Botanical Club* 23(2): 1–113.  
6. Kuhl, F. P., and C. R. Giardina. 1982. Elliptic Fourier features of a closed contour. *Computer Graphics and Image Processing* 18(3): 236–258. doi: 10.1016/0146-664X(82)90034-X  
7. Mabberley, D. J., D. K. Trock, and A. S. Weakley. 2020. The nomenclature of *Packera dubia* (Asteraceae: Senecioneae). *Taxon* 69(6): 1334–1337. doi: 10.1002/tax.12351  
8. Northcutt, C. G., L. Jiang, and I. L. Chuang. 2021. Confident Learning: Estimating Uncertainty in Dataset Labels. *Journal of Artificial Intelligence Research* 70: 1373–1411. doi: 10.1613/jair.1.12125  
9. Šlenker, M., P. Koutecký, and P. Marhold. 2022. MorphoTools2: an R package for multivariate morphometric analysis. *Bioinformatics* 38(10): 2954–2955. doi: 10.1093/bioinformatics/btac173  
10. Trock, D. K. 2006. *Packera*. In Flora of North America Editorial Committee (eds.), *Flora of North America North of Mexico*, Vol. 20, 570–602. Oxford University Press, New York.  
11. Weakley, A. S. 2026. *Flora of the Southeastern United States*. University of North Carolina Herbarium (NCU), North Carolina Botanical Garden, Chapel Hill.  
12. Weaver, W. N., P. S. Ng, and R. LaFrance. 2024. LeafMachine2: Using machine learning to rapidly measure plant traits captured in herbarium specimens. *Applications in Plant Sciences* 12(1): e11545. doi: 10.1002/aps3.11545

---

## 📄 License & Attribution

This project is licensed under the **MIT License**. Herbarium specimen images harvested through the pipeline remain subject to the individual institutional data and copyright policies of the contributing herbaria.
  
