# Operational Workflow Guide: End-to-End Species Delimitation & Morphometrics Pipeline
### North American *Packera* Á. Löve & D. Löve (*P. dubia* Complex Baseline) (Asteraceae: Senecioneae)
**University of North Carolina at Chapel Hill Herbarium (NCU)**  
**Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV5`  
**Author:** J. Brandon Fuller (PhD Candidate, UNC-CH Department of Biology)  
**Advisor:** Dr. Alan S. Weakley (Director, NCU Herbarium)

---

## 📋 Table of Contents

1. [Environment & Setup](#1-environment--setup)
   - [Windows / WSL 2 Setup Guide (`docs/WINDOWS_WSL2_SETUP.md`)](WINDOWS_WSL2_SETUP.md)
   - [A. Primary Python Pipeline Environment (`.venv`)](#a-primary-python-pipeline-environment-venv)
   - [B. LeafMachine2 Dedicated Virtual Environment (`.venv_LM2`)](#b-leafmachine2-dedicated-virtual-environment-venv_lm2)
   - [C. R Statistical Computing Environment](#c-r-statistical-computing-environment)
   - [⚡ Quick-Start Copy-Paste Terminal Cheat-Sheet](#-quick-start-copy-paste-terminal-cheat-sheet)
   - [Recommended Multi-Stage Execution Sequence](#recommended-multi-stage-execution-sequence)
2. [Eight-Phase Automated Pipeline Overview](#2-eight-phase-automated-pipeline-overview)
3. [Phase 1: Voucher Ingestion & Strict Active/Passive Authority Stratification](#phase-1-voucher-ingestion--strict-activepassive-authority-stratification)
4. [Phase 2: Stock LeafMachine2 v2.2+ Segmentation, Boundary-Continuity Lobe Rejection & Classical Gatekeeping](#phase-2-stock-leafmachine2-v22-segmentation-boundary-continuity-lobe-rejection--classical-gatekeeping)
5. [Phase 3: Strict Tier-1 Closed-Contour EFA, Inertial Alignment, Symmetrization & OTU Medians](#phase-3-strict-tier-1-closed-contour-efa-inertial-alignment-symmetrization--otu-medians)
6. [Phase 4: Klingenberg Allometry, Dual CDA & Parametric MANOVA](#phase-4-klingenberg-allometry-dual-cda--parametric-manova)
7. [Phase 5: Taphonomic Indumentum-Blind Vision XAI, Foreground Patch Pooling & OOD Scoring](#phase-5-taphonomic-indumentum-blind-vision-xai-foreground-patch-pooling--ood-scoring)
8. [Phase 6: 0–15 cm Rhizosphere/Outcrop ENM Stack, Uncertainty Buffers & Partitioned Niche Equivalency](#phase-6-015-cm-rhizosphereoutcrop-enm-stack-uncertainty-buffers--partitioned-niche-equivalency)
9. [Phase 7 & 8: Reliability-Weighted Consensus, 5-Stage Decision Matrix & Complex Discovery](#phase-7--8-reliability-weighted-consensus-5-stage-decision-matrix--complex-discovery)
10. [Automated Test Suite Verification](#10-automated-test-suite-verification)
11. [Key Literature & Citations](#11-key-literature--citations)

---

## 1. Environment & Setup

> [!TIP]
> **Windows Users (WSL 2 Setup Guide):**  
> If running on Windows, the recommended and supported execution environment is **Windows Subsystem for Linux (WSL 2 with Ubuntu)**. Please follow our step-by-step setup guide in [`docs/WINDOWS_WSL2_SETUP.md`](WINDOWS_WSL2_SETUP.md).

The pipeline operates across two dedicated Python environments and an R statistical environment:

### A. Primary Python Pipeline Environment (`.venv`)
Used for orchestration (`main.py`), voucher harvesting, `dinov2_vitb14_reg` foreground patch embeddings, Cleanlab XAI, and quality assurance tests:
```bash
source .venv/bin/activate
pip install -r requirements.txt
```

### B. LeafMachine2 Dedicated Virtual Environment (`.venv_LM2`)
Houses stock LeafMachine2 v2.2+, YOLOv5 (`PLANT_LeafPriority` / `LeafPriority.pt` and the 38-class Ruler Detector), and compiled Detectron2 PointRend (`leaf_seg_rcnn_v2_1`) with strict ABI pins (`numpy<2.0`, `opencv-contrib-python<=4.10.0.84`):
```bash
bash setup_leafmachine2.sh
```

> [!NOTE]
> **Zero Custom Retraining or Manual Annotation Required:**  
> The automated production pipeline uses stock LeafMachine2 v2.2+ `PLANT_LeafPriority` and `leaf_seg_rcnn_v2_1` weights out of the box. The standalone utilities in `scripts/annotation_and_training/` and `docs/SAM2_Precision_Botanical_Annotation_Guide.txt` are archived/optional tools not required by the automated pipeline.

### C. R Statistical Computing Environment
Requires R $\ge 4.3$ with CRAN statistical and spatial macroecological packages (`remotes`, `tidyverse`, `optparse`, `mclust`, `clustvarsel`, `MorphoTools2`, `RRPP`, `vegan`, `terra`, `sf`, `geodata`, `soilDB`, `spatialRF`, `ranger`, `adespatial`, `blockCV`, `ecospat`, `ENMTools`, `igraph`, `patchwork`, `viridis`) and `Momocs` (`remotes::install_github("MomX/Momocs")`):
```bash
Rscript -e 'install.packages(c("remotes", "tidyverse", "optparse", "mclust", "clustvarsel", "MorphoTools2", "RRPP", "vegan", "terra", "sf", "geodata", "soilDB", "spatialRF", "ranger", "adespatial", "blockCV", "ecospat", "ENMTools", "igraph", "patchwork", "viridis"), repos="https://cloud.r-project.org"); if (!requireNamespace("Momocs", quietly=TRUE)) remotes::install_github("MomX/Momocs")'
```

---

### ⚡ Quick-Start Copy-Paste Terminal Cheat-Sheet

#### A. Environment Activation & Preflight Verification
```bash
source .venv/bin/activate
python main.py check-env
```

#### B. End-to-End Automated Pipeline (`run_pipeline.sh` CLI & Dual-Strategy Single-Taxon Modes)
Every pipeline execution via `./run_pipeline.sh` or `python main.py` operates under the **Canonical Two-Tier Storage Architecture**: invariant Step 01–02 specimen/organ assets reside in canonical root `data/` and are automatically seeded/symlinked into `runs/<RUN_ID>/` (updating `runs/latest`), while all Step 03–08 statistical tables, reports, logs, and figures (`.pdf` + 300-DPI `.png`) are isolated inside `runs/<RUN_ID>/`:
```bash
source .venv/bin/activate

# 1. Launch a new timestamped Eastern NA baseline run (creates runs/run_<YYYYMMDD_HHMMSS>/ and updates runs/latest)
./run_pipeline.sh

# 2. Launch a custom named run (--run-id <NAME>) reusing cached vouchers & local images (--skip-harvest)
./run_pipeline.sh --run-id packera_full_run --skip-harvest

# 3. Resume an existing run from any step (--start-step <1..7>, e.g., Step 03 EFA or Step 04 GMM/CDA onward)
./run_pipeline.sh --run-id packera_full_run --skip-harvest --start-step 3

# 4. Strategy 1: Focal Taxon + Outgroup Gold Reference Anchors (--focal-taxon "<TAXON>" --gold-anchors 30)
# (Cascades Gold -> Silver -> Bronze anchors when an outgroup taxon lacks Tier_1_Gold vouchers)
./run_pipeline.sh --run-id dubia_triage --skip-harvest --focal-taxon "Packera dubia" --gold-anchors 30

# 5. Strategy 2: Pure Intraspecific Delimitation Mode (--focal-taxon "<TAXON>" --intraspecific-mode)
./run_pipeline.sh --run-id dubia_cryptic --skip-harvest --focal-taxon "Packera dubia" --intraspecific-mode

# 6. Continent-wide North American expansion across US, Canada, and Mexico with dynamic synonym expansion
./run_pipeline.sh --run-id na_wide_run --geographic-scope north_america --include-western --countries "US,CA,MX" --expand-synonyms

# 7. Harvest and process all Eastern North American taxa with graph modularity complex discovery
./run_pipeline.sh --run-id eastern_complexes --all-eastern-taxa --discover-complexes --spatial-thin-km 5.0

# 8. Post-QC Step 02 montage regeneration (called automatically after Step 04)
python scripts/pipeline/02_segment_and_extract.py --regenerate-figures-only

# 9. Step 08: Standalone North American Taxon-Complex Discovery
Rscript scripts/analysis/08_discover_taxon_complexes.R \
    --harmonics data/tables/leaf_efa_harmonics.csv \
    --vouchers data/tables/curated_vouchers.csv \
    --multimodal data/tables/multimodal_conflict_flags.csv \
    --vision data/tables/label_noise_audit.csv \
    --output-dir outputs/reports
```

---

### Recommended Multi-Stage Execution Sequence

To execute a complete, rigorous systematic revision of *Packera dubia* and allied taxa while avoiding cold cache downloads, execute the following 5-stage protocol:

| Stage | Objective | Production Command |
| :--- | :--- | :--- |
| **Stage 1: Baseline Harmonization** | Full 4-taxon harvest, ruler calibration, LeafMachine2 segmentation, EFA, GMM/CDA, DINOv2 XAI, and multimodal random forests. | `./run_pipeline.sh --run-id stage1_baseline` |
| **Stage 2: Promote Manifests** | Verify and promote generated manifests into canonical shared store `data/tables/`. | `cp runs/stage1_baseline/data/tables/*.csv data/tables/` |
| **Stage 3: Pure Intraspecific Mode** | Delimit cryptic lineages and morphotypes within *P. dubia* strictly (`n_unique_taxa == 1`), gating Step 04–06 on GMM clusters and `PCA-env` ecotypes. | `./run_pipeline.sh --run-id stage3_dubia_cryptic --skip-harvest --focal-taxon "Packera dubia" --intraspecific-mode` |
| **Stage 4: Outgroup Gold Anchors** | Calibrate *P. dubia* against 30 reference anchors per outgroup species (cascading Gold -> Silver -> Bronze) for definitive curation triage. | `./run_pipeline.sh --run-id stage4_dubia_triage --skip-harvest --focal-taxon "Packera dubia" --gold-anchors 30` |
| **Stage 5: Full Tree Discovery** | Expand harvest across Eastern North America or continent-wide (`--geographic-scope north_america`) and perform unsupervised network community discovery. | `./run_pipeline.sh --run-id stage5_eastern_tree --all-eastern-taxa --discover-complexes` |

---

## 2. Eight-Phase Automated Pipeline Overview

```mermaid
flowchart TD
    subgraph P1["Phase 1: Ingestion, EPA Ecoregions & 3-Layer Authority Engine"]
        A["GBIF / iDigBio Harvesting (01_voucher_harvester.py)\n+ 26-col Darwin Core Schema · Dynamic Synonyms"] --> B["EPA Level III Ecoregions + Centroid Fallback (>= 95%)\n3-Layer Authority: Protologue + Specialists + Herbaria"]
    end

    subgraph P2["Phase 2: Stock LeafMachine2 v2.2+ & Classical Gatekeeping"]
        B --> C["Stock YOLOv5 PLANT_LeafPriority (Conf >= 0.65)\n+ 38-Class Ruler Calibration (lm2_data_loader.py)"]
        C --> D["12%-Padded Crops -> Stock PointRend (leaf_seg_rcnn_v2_1)\n+ Boundary Filter (<= 0.15) & Calibrated Area (>= 0.75 cm2)"]
        D --> D1["Intact Tier_1_Direct Contours (Solidity >= 0.72)\n+ 100% Label-Blind Dissection Gates & Petiole/Height"]
        D --> D2["Organ RGB Crops (rosettes_dense · basal_leaves_raw · capitula)\n+ Post-QC Montage (--regenerate-figures-only)"]
    end

    subgraph P3["Phase 3: Strict Tier-1 Closed-Contour EFA (Momocs)"]
        D1 --> E["Inertial Axis Alignment (coo_alignxax), >=99% Power (H=12),\nUniform Symmetrization (rm_asym) & MCD Outlier Gating (p < 0.001)"]
        E --> F["Unweighted OTU Median Harmonics & PCs\n(data/tables/leaf_efa_harmonics.csv · Tier_1_Pure)"]
    end

    subgraph P4["Phase 4: Klingenberg Allometry, Dual CDA & Parametric MANOVA"]
        F --> G["Slope Homogeneity Test (PC ~ log_centsize * taxon),\nAdaptive Bronze Promotion + Empirical Bayes Shrinkage (alpha0=5.0)"]
        G --> H["Track 1: Unsupervised Morphotype CDA (fig04)\nTrack 2: Supervised EPA Level III CDA (fig04c) + MANOVA Reports"]
    end

    subgraph P5["Phase 5: Indumentum-Blind Vision XAI & Cleanlab"]
        D2 --> I["224x224 Low-Pass (sigma=1.0) & 5-px Dilated Neutral Paper Masking\nForeground 16x16 Patch Pooling (dinov2_vitb14_reg), GroupKFold(institutionCode)"]
        I --> J["Platt/Temp Calibration, Class-Conditional Cleanlab & Pointing Game >= 0.65\nEmpirical CIE L* / HSV S Indumentum Discordance Audit"]
    end

    subgraph P6["Phase 6: 0-15cm Non-Blended Edaphic Stack & Partitioned Niche"]
        B & H & J --> K["0-15cm Trapezoidal Depth Weighting, Non-Blended Hierarchy\n(100% SSURGO / 30m POLARIS / 100m SOLUS100 vs SoilGrids Fallback)"]
        K --> L["30m 3DEP Microtopography, SGMC Bedrock, Uncertainty Buffers\nAll-Class SGMC/SSURGO Fisher BH-FDR + Organic PCA-env Ecotypes"]
    end

    subgraph P7["Phase 7 & 8: Consensus Triage & Complex Discovery"]
        L --> M["Out-of-Fold Kappa Consensus & Hopkins Circular Anthesis Residuals\nTaxon-Agnostic 5-Stage Decision Matrix & North American Leiden/Ward.D2"]
    end
```

---

## Phase 1: Voucher Ingestion & Strict Active/Passive Authority Stratification

- **Scripts:** [`scripts/data_prep/01_voucher_harvester.py`](../scripts/data_prep/01_voucher_harvester.py), [`scripts/core/harvester.py`](../scripts/core/harvester.py), [`scripts/core/config.py`](../scripts/core/config.py)
- **Empirical Protocol:**
  1. **Dynamic Taxon Standardization & GBIF Synonym Expansion (`standardize_packera_taxon`):** Harmonizes binomials and basionyms across all 24+ Eastern North American *Packera* taxa (e.g., *Senecio pauperculus*, *Senecio smallii*, *Senecio plattensis*, *Senecio anonymous*, *Senecio obovatus*, *Senecio aureus*, etc.) into valid accepted binomials while handling arbitrary dynamic North American, boreal, and Mexican *Packera* species via `expand_taxon_synonyms_gbif()`. Configurable geographic filtering via `--geographic-scope {eastern_na,north_america,global}` (default: `eastern_na`, excluding longitudes $< -106.65$; `north_america` and `global` disable the western longitude boundary), `--include-western` (retaining occurrences west of $-106.65$), and `--countries` (default `["US"]`, expandable to `["US", "CA", "MX"]`).
  2. **Canonical 26 Darwin Core Columns:** Exports `data/tables/curated_vouchers.csv` with exactly: `catalogNumber`, `institutionCode`, `scientificName`, `species_raw`, `recordedBy`, `recordNumber`, `identifiedBy`, `determiner_raw`, `determiner_tier`, `type_status`, `county`, `stateProvince`, `decimalLatitude`, `latitude`, `decimalLongitude`, `longitude`, `coordinateUncertainty`, `year`, `month`, `day`, `eventDate`, `regional_group`, `image_path`, `exsiccatae_key`, `is_primary_duplicate`, `duplicate_count`.
  3. **Official US EPA Level III Ecoregion Attribution (`regional_group`):** Points are evaluated against official Level III polygons (`data/environmental/epa_ecoregions_l3/`, Omernik & Griffith 2014). For ungeoreferenced sheets with valid county metadata, the county polygon is checked for containment: if $\ge 95\%$ of county area resides in a single Level III ecoregion, that ecoregion is assigned; otherwise `NA`.
  4. **3-Layer Dynamic Expert & Herbarium Authority Engine (`parse_determiner_tier`):**
     - **Layer 1 (Genus-Wide Monographer Heuristics):** Expanded `SPECIALIST_PATTERNS` covering North American *Packera* and *Senecioneae* monographers and flora authors (Barkley, Trock, Kowal, Weakley, Bain, Mahoney, Fuller, Freeman, Cronquist, Greenman, Weber, Á. Löve & D. Löve, Spellenberg, Pippen, Uttal, Fernald, Rydberg, Heller, Kral, Godfrey, Sorrie, Nesom, Strother, Turner, Golden).
     - **Layer 2 (Automatic Protologue Author Extraction):** `extract_protologue_authors()` parses basionym and combination authors directly from Darwin Core `scientificName` strings; sheets annotated by a nomenclatural author of that taxon are classified as `Tier_1_Gold`.
     - **Layer 3 (Corpus-Derived Multi-Herbarium Determiner & Regional Herbarium Scoring):** `discover_corpus_authorities()` promotes any determiner who has identified $\ge 5$ specimens of a *Packera* taxon across $\ge 3$ distinct herbaria (`institutionCode`s)—or annotated $\ge 1$ verified type specimen—to `dynamic_specialists` (`Tier_1_Gold`), and promotes herbaria holding $\ge 1$ type specimen or ranking in the top quartile of holdings ($\ge 10$ sheets with $\ge 70\%$ determiner attribution) to `dynamic_herbaria` (`Tier_2_Silver`). General floristic collections and unverified aggregator records receive `Tier_3_Bronze`.

---

## Phase 2: Stock LeafMachine2 v2.2+ Segmentation, Boundary-Continuity Lobe Rejection & Classical Gatekeeping

- **Scripts:** [`scripts/pipeline/02_segment_and_extract.py`](../scripts/pipeline/02_segment_and_extract.py), [`scripts/vision/configure_leafmachine2.py`](../scripts/vision/configure_leafmachine2.py), [`scripts/vision/lm2_geometry_utils.py`](../scripts/vision/lm2_geometry_utils.py)
- **Empirical Protocol:**
  1. **Stock LeafMachine2 v2.2+ Configuration (`configure_leafmachine2.py`):** Uses stock `PLANT_LeafPriority` (`minimum_confidence_threshold: 0.65`) and compiled Detectron2 PointRend (`leaf_seg_rcnn_v2_1`) weights with `treat_leaflet_as_leaf_whole: False`, `use_efds_for_png_masks: False`, and `do_apply_landmarks: False`.
  2. **Otsu Boundary-Continuity Lobe Rejection (`is_truncated_lobe_or_leaflet`):** Rejects severed lobes or unattached leaf fragments where non-paper plant tissue occupies $> 15\%$ of the bounding-box perimeter (`max_tissue_border_ratio: 0.15`).
  3. **Ruler-Calibrated Physical Area Floor (`min_blade_area_cm2: 0.75`):** Rejects calibrated fragments below $0.75\text{ cm}^2$.
  4. **Classical Dimensionless Descriptors & 100% Label-Blind Botanical Dissection Gates (`classify_and_route_leaf`):** Routes intact basal blades into `Tier_1_Direct` when $\text{Solidity} \ge 0.72$, $\text{Circularity} \ge 0.25$, $\text{Extent} \ge 0.35$, and $\text{UCS} \ge 0.85$. Dissected blades are evaluated via `is_botanical_dissection()` ($0.48 \le \text{Solidity} < 0.72$ pinnatifid; $0.18 \le \text{Solidity} < 0.48$ with $\ge 4$ deep sinuses bipinnatisect) in a 100% label-blind manner (`taxa_with_lyrate_tendency` is retained strictly as a deprecated no-op). Diagnostic figures utilize dynamic species sets (`present_sp`), empirical per-taxon median aspect ratios, and dynamic complex titles (`format_active_complex_title`).
  5. **Resolution-Adaptive Width Profiling & Cauline Contaminant Rejection:** Profiles cross-sectional widths via $w_{\text{span}} = \max(\min(\max(h, w), 350), \lfloor 0.65 \max(h, w) \rfloor)$ and scales petiole thresholds when $\text{axis\_len} > 450\text{ px}$. Rejects upper-stem sessile/clasping cauline bracts (`aspect_ratio > 3.8` and `solidity < 0.78` or displaced $> 2.2 \times$ median rosette radius).
  6. **Crown Point & Petiole Length Estimation:** Estimates rosette crown point $p_{\text{crown}}$ from inward basal midrib ray intersections and Class 9 root clumps. Per-leaf petiole length is $L_{\text{petiole}} = \max(d_{\text{geodesic}}(p_{\text{base}}, p_{\text{crown}}), \|p_{\text{base}} - p_{\text{crown}}\|_2)$. Outer-rosette petiole length is aggregated at the 75th percentile ($P_{75}$).
  7. **Scape Height & Stem Posture:** Computes straight vs. bent scape height: $H_{\text{plant}} = \|p_{\text{bend}} - p_{\text{crown}}\|_2 + \|p_{\text{corymb}} - p_{\text{bend}}\|_2$.
  8. **Post-QC Step 02 Montage Regeneration:** Automatically regenerated after Step 04 outlier pruning via `python scripts/pipeline/02_segment_and_extract.py --regenerate-figures-only`.

---

## Phase 3: Strict Tier-1 Closed-Contour EFA, Inertial Alignment, Symmetrization & OTU Medians

- **Script:** [`scripts/morphometrics/03_fourier_extractor.R`](../scripts/morphometrics/03_fourier_extractor.R)
- **Empirical Protocol:**
  1. **Strict `Tier_1_Direct` Closed-Contour Scope:** Only intact `Tier_1_Direct` basal leaf contours enter `Momocs::efourier` (mirrored hemi-blades or open-curve approximations are strictly excluded).
  2. **Inertial Axis Pre-Alignment (`Momocs::coo_alignxax`) with Suborbicular Guard:** Aligns contours along their major inertial axis when blade circularity $\le 0.75$, protecting suborbicular blades ($>0.75$) from $90^\circ$ rotation artifacts while locking base-to-apex polarity.
  3. **$\ge 99\%$ Cumulative Harmonic Power Verification (`calibrate_harmonicpower_efourier`):** Verifies that $H = 12$ harmonics capture $\ge 99\%$ cumulative harmonic power across the corpus before running `Momocs::efourier(..., nb.h = 12, norm = FALSE, start = FALSE)`.
  4. **Uniform Bilateral Symmetrization (`Momocs::rm_asym`, Claude 2008):** Isolates symmetric harmonic coefficients ($A_n, D_n$) uniformly across all intact contours to remove taphonomic pressing curvature.
  5. **MCD Mahalanobis Outlier Filtering & Unweighted OTU Medians:** Removes contour extraction anomalies via Minimum Covariance Determinant Mahalanobis distance (`MASS::cov.rob`, $p < 0.001$) and `pct_asymmetry > 0.22 | harmonic_roughness > 0.15`, then aggregates multi-leaf vouchers to the Specimen / Operational Taxonomic Unit (`(catalogNumber, plant_individual_id)`) via **unweighted coordinate-wise median** harmonics (`data/tables/leaf_efa_harmonics.csv`) and PCs (`fourier_pcs.csv`, `symmetric_fourier_pcs.csv`) with `morph_tier_source = "Tier_1_Pure"`.

---

## Phase 4: Klingenberg Allometry, Dual CDA & Parametric MANOVA

- **Scripts:** [`scripts/morphometrics/04_gmm_morphotools.R`](../scripts/morphometrics/04_gmm_morphotools.R), [`scripts/morphometrics/authority_weighting_utils.R`](../scripts/morphometrics/authority_weighting_utils.R)
- **Empirical Protocol:**
  1. **Empirical Bayes Weighting & Adaptive Reference Anchor Selection (`authority_weighting_utils.R`):**
     - **`select_adaptive_reference_anchors()`:** Initializes reference anchors with `Tier_1_Gold` and `Tier_2_Silver` `Tier_1_Pure` vouchers pruned at Mahalanobis $\chi^2_{0.95, p}$. When an under-annotated or rare *Packera* taxon has fewer than `min_anchors_per_taxon` (default 8, or $< 3$) clean Gold+Silver anchors, it automatically promotes the most central `Tier_3_Bronze` `Tier_1_Pure` vouchers within $\chi^2_{0.80, p}$ of the taxon medoid (`is_promoted_anchor = TRUE`), ensuring GMM and CDA never drop Gold-poor taxa.
     - **`compute_empirical_bayes_authority_weights()`:** Computes unified Beta-Binomial Empirical Bayes shrinkage weights ($\alpha_0 = 5.0$, prior reliabilities $\pi_{\text{Gold}} = 0.95$, $\pi_{\text{Silver}} = 0.85$, $\pi_{\text{Promoted Bronze}} = 0.65$, $\pi_{\text{Unpromoted Bronze}} = 0.35$):
       $$w_{t, \text{tier}} = \frac{\alpha_0 \pi_{\text{tier}} + \sum_{i \in (t, \text{tier})} (1 - \hat{\epsilon}_i)}{\alpha_0 + n_{t, \text{tier}}}$$
       In `"active_passive"` mode: unpromoted `Tier_3_Bronze` vouchers receive `0.0` (passive projection), while active anchors receive precision weights rescaled to mean `1.0`. Step 04 uses taxon-agnostic discordance strings (`<raw_sp>_Discordant_to_<pred_sp>`) and `"Unassigned_Packera"` fallbacks.
  2. **Klingenberg (2016) Allometric Slope Homogeneity Testing (`test_and_correct_allometry`):** Tests `PC ~ log_centsize * species_standardized` ($\alpha = 0.01$) on ruler-calibrated (`is_calibrated_mm == TRUE`) active anchors (`outputs/reports/allometry_homogeneity_test.csv`, `outputs/figures/shape_allometry_regression.{pdf,png}`). Applies common-slope within-group size correction only when slopes are homogeneous ($p \ge 0.01$).
  3. **Horn's Parallel Analysis & Williams-Titus Cap (`select_informative_pcs`):** Retains informative PCs exceeding the Broken-Stick / Horn threshold, capped for CDA at $P \le \lfloor N_{\min} / 3 \rfloor$ (Williams & Titus 1988).
  4. **Kass & Raftery (1995) GMM (`run_gmm_cluster_analysis`):** Evaluates discrete morphotypes via $2\Delta\text{BIC}$ thresholds ($> 6$ strong, $> 10$ very strong; `outputs/reports/gmm_bayes_factors_summary.csv`).
  5. **Dual Clustering Architecture:**
     - **Track 1: Unsupervised GMM Morphotype CDA:** Fits CDA on Kass-Raftery morphotype clusters (`Morphotype_K1..Kn`) with uniform class priors (`rep(1/g, g)`) and Leave-One-Out Cross-Validation, generating `outputs/figures/fig04_gmm_cda_clustering.{pdf,png}`.
     - **Track 2: Supervised EPA Level III Ecoregion CDA:** Fits CDA on ecoregional ecotypes, exporting `outputs/reports/ecoregion_supervised_cda_summary.csv` and `outputs/figures/fig04c_supervised_ecoregion_cda.{pdf,png}`.
  6. **Parametric Omnibus & Pairwise MANOVA:** Runs Pillai's trace and Wilks' lambda omnibus tests and pairwise Hotelling's $T^2$ contrasts with Benjamini-Hochberg FDR correction (`outputs/reports/manova_group_significance.csv`, `outputs/reports/manova_pairwise_posthoc.csv`).
  7. **Taxon-Adjusted Institutional Batch-Effect Audit (`audit_institutional_batch_effects`):** Fits `lm(PC1 ~ species_standardized + institutionCode)` and `lm(PC2 ~ species_standardized + institutionCode)` when $\ge 2$ taxa are present, confirming institutional batch effect partial $R^2 < 0.05$ (`outputs/reports/institutional_batch_effect_audit.csv`).

---

## Phase 5: Taphonomic Indumentum-Blind Vision XAI, Foreground Patch Pooling & OOD Scoring

- **Script:** [`scripts/analysis/05_cleanlab_vision_xai.py`](../scripts/analysis/05_cleanlab_vision_xai.py)
- **Empirical Protocol:**
  1. **Indumentum-Blind Gaussian Low-Pass Filter ($\sigma = 1.0\text{ px}$ at $224 \times 224$):** `apply_indumentum_blind_lowpass()` attenuates taphonomically abraded trichome/tomentum textures on $224 \times 224$ whole-rosette and basal-leaf crops while preserving macro-architectural venation and lobing.
  2. **5-px Morphologically Dilated Neutral Paper Masking:** Dilates binary `Tier_1_Direct` leaf masks by 5 pixels and composites onto neutral herbarium paper RGB `[242, 239, 230]` to eliminate background mounting cues and label text artifacts.
  3. **Foreground-Masked $16 \times 16$ Spatial Patch Pooling (`dinov2_vitb14_reg`):** Downsamples the binary foreground plant mask to the $16 \times 16$ patch grid (`compute_foreground_patch_mask`) and average-pools `dinov2_vitb14_reg` spatial patch tokens strictly over plant tissue into 768-d L2-normalized feature vectors.
  4. **Active Taxa Derivation & GroupKFold Calibration:** Dynamically extracts `active_taxa` from curated vouchers and maps colorblind-safe Okabe-Ito palettes (`_get_dynamic_okabe_ito_map`). Fits out-of-fold `LogisticRegression(C=0.1, class_weight="balanced")` inside a fold-internal `StandardScaler` + `PCA(0.95)` pipeline on Gold/Silver anchors using `StratifiedGroupKFold(groups=institutionCode)`, calibrates probabilities via Platt/temperature scaling (`calibrate_oof_probabilities`), and projects `Tier_3_Bronze` passively.
  5. **Class-Conditional Cleanlab Noise Pruning, OOD Scoring & Empirical Chromatic Indumentum Audit:** Prunes label issues via `filter_by="prune_by_noise_rate"`, computes Mahalanobis + cosine $k$-NN feature-space distance (`vision_ood_score`, `vision_is_ood`), and enforces Pointing Game `foreground_attention_ratio >= 0.65` (`data/tables/label_noise_audit.csv`). Flags `tomentum_taxon_discordance` on abaxial sheets (`surface_orientation == "abaxial"`) with `given_label != predicted_label` and `c_error > error_threshold` when the predicted taxon has higher empirical median CIE $L^*$ luminance (`med_L[pred] >= med_L[given] + 2.0`) or lower HSV saturation (`med_S[pred] <= med_S[given] - 0.02`), or when `mean_luminance > Q75_L[given]`.

---

## Phase 6: 0–15 cm Rhizosphere/Outcrop ENM Stack, Uncertainty Buffers & Partitioned Niche Equivalency

- **Scripts:** [`scripts/analysis/06_multimodal_spatial_rf.R`](../scripts/analysis/06_multimodal_spatial_rf.R), [`scripts/analysis/edaphic_fine_grain_utils.R`](../scripts/analysis/edaphic_fine_grain_utils.R)
- **Empirical Protocol:**
  1. **$0\text{–}15\text{ cm}$ Trapezoidal Rhizosphere Depth-Weighting:** Integrates $(5 \cdot v_{0\text{–}5} + 10 \cdot v_{5\text{–}15}) / 15$ depth-weighted pedological horizons.
  2. **Strict Non-Blended Edaphic Hierarchy:** CONUS occurrences query 100% SSURGO, 30m POLARIS, and 100m SOLUS100; 250m SoilGrids is restricted strictly to non-CONUS vouchers where POLARIS/SOLUS100 evaluate to `NA`. Raw `soil_*` predictors are excluded from active `env_vars`.
  3. **Fine-Scale Pedological & Lithological Predictors:** Integrates `clr_sand_silt` / `clr_clay_geom` compositional texture, `cfvo_0_15cm`, `awc_0_15cm`, `fine_bedrock_depth_cm`, 50m-buffered SSURGO rock outcrop queries (`ssurgo_rock_outcrop_pct`), 30m 3DEP micro-topography (`topo_tpi_30m`, `topo_twi_30m`), and 6-class USGS SGMC bedrock lithology (`sgmc_lithology_class`).
  4. **Coordinate-Uncertainty Buffered Extraction (`100–1000 m`):** Extracts focal median raster values within circular spatial buffers of $\text{clamp}(\text{coordinateUncertaintyInMeters}, 100, 1000)\text{ m}$.
  5. **Collinearity Pruning, Spatial Block CV & Exhaustive Substrate Contingency (`synthesize_multiscale_edaphics`):** Prunes predictors at $|r| < 0.70$ and $\text{VIF} < 5.0$ (`prune_collinear_predictors`). Evaluates all 6 USGS SGMC lithological classes and SSURGO outcrop/drainage classes via $2 \times 2$ Fisher's Exact Tests with Benjamini-Hochberg FDR correction (`outputs/reports/micro_edaphic_outcrop_contingency.csv`), reporting `Dominant_Lithology_Class` and setting `substrate_specialist_flag` dynamically when outcrop frequency $> 0.30$ or odds ratio $> 2.0$. Trains `ranger::ranger` Random Forests on active anchors across $0.5^\circ \times 0.5^\circ$ spatial blocks (`spatial_block_cv_acc`).
  6. **Organic `PCA-env` Intraspecific Ecotypes & Dynamic Top-Importance Plotting:** In single-taxon intraspecific mode (`--intraspecific-mode`), partitions vouchers in standardized `PCA-env` space into $K = 2..3$ clusters (`Ecotype_1_<TopVar>_Low`, `Ecotype_2_<TopVar>_High`). Panels B and C of `fig06_ecological_niche_phenology.{pdf,png}` dynamically plot the top non-MEM environmental and edaphic predictors ranked by Random Forest permutation importance (`MeanDecreaseAccuracy`) with 85% confidence ellipses (`stat_ellipse(level = 0.85)`).
  7. **Partitioned Climatic vs. Edaphic-Geologic Niche Equivalency:** Computes Schoener's $D$ and Hellinger's $I$ separately in Climatic PCA space ($D_{\text{clim}}, I_{\text{clim}}$) and Edaphic-Geologic PCA space ($D_{\text{edaph}}, I_{\text{edaph}}$) (`outputs/reports/edaphic_niche_equivalency.csv`).

---

## Phase 7 & 8: Reliability-Weighted Consensus, 5-Stage Decision Matrix & Complex Discovery

- **Scripts:** [`scripts/analysis/07_triage_dashboard_synthesis.R`](../scripts/analysis/07_triage_dashboard_synthesis.R), [`scripts/analysis/08_discover_taxon_complexes.R`](../scripts/analysis/08_discover_taxon_complexes.R), [`scripts/morphometrics/pub_figures_utils.R`](../scripts/morphometrics/pub_figures_utils.R)
- **Empirical Protocol:**
  1. **Out-of-Fold Cohen's $\kappa$ Reliability Weighting:** Weights each modality by out-of-fold validation Cohen's $\kappa_m$ ($w_m = \max(0.05, \kappa_m)$), dynamically renormalized over observed modalities $\mathcal{O}_i$ on each specimen (`98.4%` voucher completion across georeferenced and non-georeferenced sheets).
  2. **Hopkins-Adjusted Circular Anthesis Phenology:** Restricted to specimens in active anthesis (`in_anthesis == TRUE`), regresses circular flowering angle $(\sin\theta_i, \cos\theta_i)$ on latitude, elevation, and MAT (`atan2` circular residuals) before computing $|\Delta\text{DOY}_{\text{Hopkins}}|$.
  3. **Taxon-Adjusted Institutional ANCOVA:** Audits institutional batch effects with Benjamini-Hochberg FDR correction.
  4. **Taxon-Agnostic 5-Stage Hierarchical Decision Matrix (`apply_taxonomic_decision_matrix`):**
     - **Stage 1 (`Misidentification_Reassign`):** Evaluated **before** `Hybrid_Introgressant` so sympatric misidentifications are never misclassified as hybrids.
     - **Stage 2 (`Hybrid_Introgressant`):** Morphological intermediate with bi-parental geographic/edaphic sympatry; resolves hybrid parents dynamically via nearest-neighbor Mahalanobis distance (`p_b`).
     - **Stage 3 (`Subspecies_Ecotype`):** Evaluates `is_cohesive_ecotype` (infraspecific rank or secondary GMM cluster with cohesive edaphic/geographic divergence).
     - **Stage 4 (`Ecophenotypic_Plasticity`):** Evaluates `is_taphonomic_plasticity` (edaphic $z$-score $\le 2.0$ with chromatic discordance or late-season $\text{DOY} \ge Q_{75}$ with single-modality divergence).
     - **Stage 5 (`Putative_Cryptic_Species`):** Requires multi-modal **Integration by Congruence** (Padial et al. 2010).
  5. **Dynamic Taxonomic Revision Summaries & Publication Figures:** Generates empirical treatments for every taxon in `active_taxa`, exporting `<Focal_Taxon>_Taxonomic_Revision_Summary.md` and `Figure_Integrative_<Focal_Taxon>_Revision.pdf` when `--focal-taxon` is specified. Publication figures dynamically assign Okabe-Ito colors (`get_dynamic_packera_palette`) and hollow shapes (`get_dynamic_packera_shapes`), formats dynamic complex titles (`format_active_complex_title`), and computes dynamic North American map bounding boxes (`xlim_map`, `ylim_map` padded by 12%).
  6. **Step 08 Unsupervised North American Taxon-Complex Discovery (`scripts/analysis/08_discover_taxon_complexes.R`):** Fuses pairwise morphological ($D_{\text{Mahal}}$), environmental ($1 - D_{\text{Schoener}}$), phenological ($d_{\text{circ}}$), and vision centroid distances across all regional *Packera* species and partitions them via $k$-NN graph modularity (Leiden/Louvain) and Ward.D2 clustering (`outputs/reports/discovered_taxon_complexes.csv`).
  7. **Taxon x GMM Morphotype Subcluster Network:** When $\le 2$ species are analyzed, Step 08 automatically constructs a Taxon x GMM Morphotype subcluster network, resolving cryptic within-taxon structure (`outputs/figures/fig_empirical_taxon_complex_network.{pdf,png}`).

---

## 10. Automated Test Suite Verification

Run the complete unit test suite from the repository root:
```bash
python3 -m unittest discover -s scripts/tests -p "test_*.py" -v
```

---

## 11. Key Literature & Citations

1. Barkley, T. M. 1988. Variation among the Senecioneae (Asteraceae) in North America. *Brittonia* 40(2): 211–221.
2. Broennimann, O., et al. 2012. Measuring ecological niche overlap from occurrence and spatial environmental data. *Global Ecology and Biogeography* 21(4): 481–497.
3. Claude, J. 2008. *Morphometrics with R*. Springer, New York.
4. Collyer, M. L., and D. C. Adams. 2018. RRPP: An R package for fitting linear models to high-dimensional data using residual randomization. *Methods in Ecology and Evolution* 9(7): 1772–1779.
5. Darcet, T., M. Oquab, J. Mairal, and P. Bojanowski. 2024. Vision transformers need registers. *International Conference on Learning Representations (ICLR)*.
6. Dawid, A. P., and A. M. Skene. 1979. Maximum likelihood estimation of observer error-rates using the EM algorithm. *Journal of the Royal Statistical Society: Series C* 28(1): 20–28.
7. de Queiroz, K. 2007. Species concepts and species delimitation. *Systematic Biology* 56(6): 879–886.
8. Kass, R. E., and A. E. Raftery. 1995. Bayes Factors. *Journal of the American Statistical Association* 90(430): 773–795.
9. Klingenberg, C. P. 2016. Size, shape, and form: concepts of allometry in geometric morphometrics. *Development Genes and Evolution* 226(3): 113–137.
10. Kowal, R. R. 1975. Systematics of *Senecio aureus* and allied species on the Gaspé Peninsula, Quebec. *Memoirs of the Torrey Botanical Club* 23(2): 1–113.
11. Kuhl, F. P., and C. R. Giardina. 1982. Elliptic Fourier features of a closed contour. *Computer Graphics and Image Processing* 18(3): 236–258.
12. Northcutt, C. G., L. Jiang, and I. L. Chuang. 2021. Confident Learning: Estimating Uncertainty in Dataset Labels. *Journal of Artificial Intelligence Research* 70: 1373–1411.
13. Padial, J. M., A. Miralles, I. De la Riva, and M. Vences. 2010. The integrative future of taxonomy. *Frontiers in Zoology* 7: 16.
14. Scrucca, L., and A. E. Raftery. 2018. clustvarsel: A package implementing variable selection for Gaussian model-based clustering in R. *Journal of Statistical Software* 84(1): 1–28.
15. Šlenker, M., P. Koutecký, and P. Marhold. 2022. MorphoTools2: an R package for multivariate morphometric analysis. *Bioinformatics* 38(10): 2954–2955.
16. Trock, D. K. 2006. *Packera*. In Flora of North America Editorial Committee (eds.), *Flora of North America North of Mexico*, Vol. 20, 570–602. Oxford University Press, New York.
17. Valavi, R., J. Elith, J. J. Lahoz-Monfort, and G. Guillera-Arroita. 2019. blockCV: An R package for generating spatially or environmentally separated folds for k-fold cross-validation of species distribution models. *Methods in Ecology and Evolution* 10(2): 225–232.
18. Warren, D. L., R. E. Glor, and M. Turelli. 2008. Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. *Evolution* 62(11): 2868–2883.
19. Weaver, W. N., P. S. Ng, and R. LaFrance. 2024. LeafMachine2: Using machine learning to rapidly measure plant traits captured in herbarium specimens. *Applications in Plant Sciences* 12(1): e11545.
20. Williams, B. K., and K. Titus. 1988. Assessment of sampling stability in ecological applications of discriminant analysis. *Ecology* 69(4): 1275–1285.
21. Wright, M. N., and A. Ziegler. 2017. ranger: A fast implementation of random forests for high dimensional data in C++ and R. *Journal of Statistical Software* 77(1): 1–17.
