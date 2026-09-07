# Archived Scripts & Legacy Artifacts — Custom CV/ML Pipeline

**Archived:** 2026-08-30  
**Archived by:** J. Brandon Fuller (University of North Carolina at Chapel Hill / NCU Herbarium)

---

## Reason for Archival

The botanical computer vision workflow transitioned from a fully custom YOLOv8/SAHI/SAM 2 segmentation model to **LeafMachine2 (LM2)** (Weaver et al. 2024). LeafMachine2 provides an actively maintained, pre-trained Plant Component Detector (PCD) purpose-built for herbarium specimen sheets with integrated PointRend boundary refinement and automatic scale ruler isolation.

This architectural shift eliminates the need to maintain custom SAM 2 interactive polygon annotation tools, custom YOLO training loops, SAHI patch inference engines, and custom layout artifact gatekeepers.

All downstream downstream R and Python morphometrics, deep vision XAI, spatial ecological modeling, and triage synthesis modules (`03_fourier_extractor.R`, `04_gmm_morphotools.R`, `05_cleanlab_vision_xai.py`, `06_multimodal_spatial_rf.R`/`.py`, `07_triage_dashboard_synthesis.R`/`.py`) remain active and operate on LM2 structured outputs.

---

## Structure of the Archive

| Subdirectory | Contents | Description |
|:---|:---|:---|
| `core/` | 16 legacy modules | Custom SAHI tiler geometry, legacy artifact filters, augmentation, dataset builders, and leaf spine tracers |
| `vision/` | 4 scripts | Legacy 5-stage leaf extraction (`02_hierarchical_leaf_extractor.py`), SAHI full-sheet inference (`run_sahi_inference.py`), native-DPI tiler (`run_dpi_tiler.py`), and artifact filter gatekeeper |
| `data_prep/` | 4 scripts | Precision SAM 2 annotation GUI (`annotate_with_sam2.py`), 7-class YOLO dataset builder (`build_artifact_robust_dataset.py`), voucher image staging runner (`prepare_lm2_dataset.py`), and symlink staging utility (`staging_utils.py`) |
| `train/` | 6 modules | Custom YOLOv8m-seg fine-tuning engine with mixed-precision AMP, class weighting, and disk caching |
| `tests/` | 5 test suites | Unit and regression tests for legacy YOLO/SAM2 components |
| `analysis/` | 7 modules | Redundant Phase 5 micro-modules (`cleanlab_curator.py`, `dinov2_embeddings.py`, `gradcam_visualizer.py`), duplicate Python Phase 6 (`06_multimodal_spatial_rf.py`, `run_multimodal_spatial_rf.py`), and duplicate Python Phase 7 (`07_triage_dashboard_synthesis.py`, `run_triage_dashboard_synthesis.py`) |
| `root_artifacts/` | Pre-trained weights & models | Legacy YOLO weight checkpoints (`yolov8m-seg.pt`, `yolov8x-seg.pt`), SAM 2 submodule, and training runs |
| `../../data/_archive/` | Legacy datasets | Archived YOLO training sets, SAM2 masks, manual annotations, and tiled patches |
| `../../outputs/_archive/` | Legacy outputs | Historical SAHI detection summaries, tiling logs, synthetic benchmark runs, and GPU profiling logs |

---

## Component Migration Matrix (Legacy vs. Modernized Pipeline)

| Legacy Archived Component | Modernized Replacement Component |
|:---|:---|
| `annotate_with_sam2.py` + `build_artifact_robust_dataset.py` | LeafMachine2 pre-trained `LeafPriority.pt` detector |
| `run_dpi_tiler.py` + `scripts/_archive/core/tiling_geometry.py` | LM2 internal high-resolution multi-crop detector |
| `train/` (`trainer.py`, `train_yolo.py`) | LM2 zero-shot/pre-trained botanical organ detector |
| `run_sahi_inference.py` | LM2 native execution (`LeafMachine2/LeafMachine2.py`) |
| `02_hierarchical_leaf_extractor.py` | `scripts/pipeline/02_segment_and_extract.py` (PointRend + Midrib Reflection + Contour Export) |
| `artifact_filter_gatekeeper.py` | `scripts/vision/lm2_geometry_utils.py` (UCS + Solidity + Midrib Angle) |
| `prepare_lm2_dataset.py` + `staging_utils.py` | Direct ingestion in `scripts/pipeline/02_segment_and_extract.py` reading directly from `data/raw_vouchers/` via `data/tables/curated_vouchers.csv` |
| `cleanlab_curator.py`, `dinov2_embeddings.py`, `gradcam_visualizer.py` | Consolidated [`scripts/analysis/05_cleanlab_vision_xai.py`](../analysis/05_cleanlab_vision_xai.py) |
| `06_multimodal_spatial_rf.py` + `run_multimodal_spatial_rf.py` | Canonical R implementation [`scripts/analysis/06_multimodal_spatial_rf.R`](../analysis/06_multimodal_spatial_rf.R) |
| `07_triage_dashboard_synthesis.py` + `run_triage_dashboard_synthesis.py` | Canonical R implementation [`scripts/analysis/07_triage_dashboard_synthesis.R`](../analysis/07_triage_dashboard_synthesis.R) |

---

## Active Pipeline Entry Points

The modernized 7-stage pipeline entry points are:
1. **Phase 1:** [`scripts/data_prep/01_voucher_harvester.py`](../data_prep/01_voucher_harvester.py)
2. **Phase 2:** [`scripts/pipeline/02_segment_and_extract.py`](../pipeline/02_segment_and_extract.py) (Direct PointRend segmentation & leaf extraction)
3. **Phase 3:** [`scripts/morphometrics/03_fourier_extractor.R`](../morphometrics/03_fourier_extractor.R)
4. **Phase 4:** [`scripts/morphometrics/04_gmm_morphotools.R`](../morphometrics/04_gmm_morphotools.R)
5. **Phase 5:** [`scripts/analysis/05_cleanlab_vision_xai.py`](../analysis/05_cleanlab_vision_xai.py) (Consolidated DINOv2 embeddings, Cleanlab audit, & Grad-CAM XAI)
6. **Phase 6:** [`scripts/analysis/06_multimodal_spatial_rf.R`](../analysis/06_multimodal_spatial_rf.R) (Canonical Spatial Random Forests with MEMs & Warren's Identity)
7. **Phase 7:** [`scripts/analysis/07_triage_dashboard_synthesis.R`](../analysis/07_triage_dashboard_synthesis.R) (Canonical Multi-Evidence Triage Dashboard & Plate Generator)
8. **Pipeline Orchestrator:** [`main.py`](../../main.py) (`harvest`, `segment`, `morphometrics`, `synthesis`, and `run-all`)

Refer to the primary [`README.md`](../../README.md) and [`docs/WORKFLOW_GUIDE.md`](../../docs/WORKFLOW_GUIDE.md) for full execution guides.
