# Workspace Rules & Non-Negotiable Invariants (`Packera_dubia_morphometrics`)

**Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV4`  
**Affiliation:** University of North Carolina at Chapel Hill Herbarium (NCU)

This file defines the always-on architectural, mathematical, and cross-language data contracts for the `Packera_dubia_morphometrics` repository. Every AI agent and contributor modifying Python (`.py`) or R (`.R`) code in this workspace must enforce the five non-negotiable invariants below without exception.

---

## 1. Strict Archive Isolation

* **Zero `_archive` Coupling:** Never import from, source, read from, or fall back to `scripts/_archive/` or `data/_archive/` (or `outputs/_archive/`) in any active `.py` or `.R` file.
* **No Legacy Fallbacks:** Decommissioned scripts (such as `scripts/_archive/03_fourier_extractor.py`) and archived tables exist solely for historical provenance. Active production modules in `scripts/core/`, `scripts/data_prep/`, `scripts/pipeline/`, `scripts/vision/`, `scripts/morphometrics/`, `scripts/analysis/`, and `scripts/annotation_and_training/` must fail fast with a clear diagnostic error if an expected active input is missing rather than falling back to `_archive` paths.

---

## 2. Zero Silent Synthetic Data

* **No Unflagged Mock Generation:** Never generate random or mock embeddings (`np.random.randn`), static Gaussian attention/Grad-CAM heatmaps, or deterministic trigonometric environmental/pedological (SoilGrids / WorldClim / USDA SSURGO) values during normal pipeline execution.
* **Explicit `--demo` Requirement:** Any synthetic or simulated data generation path must be strictly gated behind an explicit `--demo` CLI flag. If real rasters, SDA network queries, or neural checkpoints are unavailable and `--demo` is not explicitly passed, raise a clear runtime error or mark records as `NA` / uncalibrated rather than fabricating values silently.
* **Never Cache Synthetic Records to Disk:** Even when running under `--demo` or test harnesses, never persist synthetic/fallback records into production cache directories or canonical tables via `saveRDS()`, `write.csv()`, `readr::write_csv()`, or `pd.DataFrame.to_csv()`.

---

## 3. Dual Virtualenv & Submodule Boundaries

* **Immutable `LeafMachine2/` Submodule:** Never modify files inside the git-ignored `LeafMachine2/` directory directly. All custom botanical adapters, configuration builders, PointRend wrappers, and geometric gatekeepers must reside in `scripts/vision/` or `scripts/pipeline/`.
* **Strict Dual-Environment Separation:**
  * **`.venv` (Primary Orchestration & Modern Vision):** Houses Python orchestration (`main.py`), Segment Anything Model 2 (`SAM 2`), `DINOv2` (`ViT-B/14`), `cleanlab` confident learning, and general scientific utilities.
  * **`.venv_LM2` (LeafMachine2 & Detectron2 PointRend Isolation):** Houses LeafMachine2, YOLOv5 (`LeafPriority.pt` and the 38-class Ruler Detector), and compiled Detectron2 PointRend (`leaf_seg_rcnn_v2_1`) with strict ABI pins: `numpy<2.0` and `opencv-contrib-python<=4.10.0.84`. Never upgrade `numpy>=2.0` or `opencv>4.10.0.84` inside `.venv_LM2`.

---

## 4. Canonical Cross-Language CSV Contracts

All Python (`scripts/pipeline/`, `scripts/vision/`, `scripts/analysis/`) and R (`scripts/morphometrics/`, `scripts/analysis/`) stages must adhere to the exact primary keys and shared column schemas below:

| Contract Element | Canonical Column / Pattern | Cross-Language Rule & Semantics |
| :--- | :--- | :--- |
| **Composite Grouping Key** | `(catalogNumber, plant_individual_id)` | Every leaf-level and specimen-level aggregation across Python and R must group by `(catalogNumber, plant_individual_id)`. **Never** collapse multi-plant herbarium sheets to `plant_individual_id = 0` or group by `catalogNumber` alone when aggregating organ traits. |
| **Routing Tier Column** | `assigned_tier` | Canonical closed-silhouette tiers are `"Tier_1_Direct"` and `"Tier_2_Reflected"` (with shorthand `"tier1"` / `"tier2"` in intermediate routing objects). In R, always match case-insensitively via `grepl("^tier_?[12]", assigned_tier, ignore.case = TRUE)`. |
| **Scale Calibration Columns** | `pixels_per_mm`, `is_calibrated_mm` | `pixels_per_mm` stores the empirical ruler conversion factor (`float` or `NA`), and `is_calibrated_mm` (`bool` / `logical`) explicitly flags whether absolute millimeter traits ($\text{CS}_{\text{mm}}$, involucre height/width in mm) are calibrated vs. uncalibrated pixel units. |
| **Morphometric Outlier Flag** | `morphometric_outlier_flag` | Exported by `04_gmm_morphotools.R` (`data/tables/morphometric_flags.csv` / `morphometrics_misidentification_flags.csv`) as `morphometric_outlier_flag` and aliased to `misidentification_flag` when consumed in Phase 6 (`06_multimodal_spatial_rf.R`) and Phase 7 (`07_triage_dashboard_synthesis.R`). |

---

## 5. Atomic Persistence

* **Write-Then-Replace Pattern:** Every tabular export (`.csv`), JSON manifest (`.json`), and model checkpoint write must write to a temporary `.tmp` file in the target directory, flush buffers (`f.flush()` and `os.fsync(f.fileno())` in Python; atomic `file.rename()` from `.tmp` in R), and atomically replace the destination via `Path.replace()` (or `os.replace()`).
* **Zero Partial Artifacts:** Never write directly to final output paths in-place where an interrupted process or OOM signal could leave truncated CSVs or corrupted checkpoints.

---

## Workspace Skills Reference

Activate the repository's specialized skills in `.agents/skills/` whenever working on the corresponding domains:

* **`/morphometric-cv-guardrails`** (`.agents/skills/morphometric-cv-guardrails/SKILL.md`): Mandatory whenever modifying scripts in `scripts/vision/`, `scripts/pipeline/`, `scripts/annotation_and_training/`, or `scripts/morphometrics/`.
* **`/pipeline-contract-audit`** (`.agents/skills/pipeline-contract-audit/SKILL.md`): Mandatory before committing code changes or when debugging cross-stage Python/R pipeline handshakes.
