# Archived Scripts — Custom CV/ML Pipeline

**Archived:** 2026-08-30
**Archived by:** J. Brandon Fuller (UNC-CH)

## Reason for Archival

Workflow transitioned from a fully custom computer vision and model-training pipeline to
**LeafMachine2 (LM2)** for automated plant component detection, leaf extraction, and organ
segmentation. LM2 provides a pre-trained, actively maintained detection framework specifically
designed for herbarium specimens, removing the need to maintain custom annotation tooling,
YOLO fine-tuning, SAHI inference, and artifact gatekeeper logic.

The downstream R morphometrics (`03_fourier_extractor.R`, `04_gmm_morphotools.R`) and
multi-modal analysis (`05_cleanlab_vision_xai.py`, `06_multimodal_spatial_rf.R`) stages
are unchanged and remain active in `scripts/`.

---

## What Was Archived

| Subdirectory | Files | Description |
|---|---|---|
| `vision/` | `02_hierarchical_leaf_extractor.py`, `artifact_filter_gatekeeper.py`, `run_dpi_tiler.py`, `run_sahi_inference.py` | Custom 5-stage hierarchical leaf extraction, artifact gatekeeper, native-DPI tiler, and SAHI full-sheet inference runner |
| `data_prep/` | `annotate_with_sam2.py`, `build_artifact_robust_dataset.py` | Interactive SAM 2 ground-truth annotation GUI and 7-class artifact-robust YOLO dataset builder |
| `train/` | `train_yolo.py`, `trainer.py`, `evaluator.py`, `dataset.py`, `config.py`, `__init__.py` | Custom YOLOv8m-seg fine-tuning package with mixed-precision AMP training, botanical class weighting, disk caching, and cross-classification evaluation |
| `tests/` | `test_botanical_topology_classifier.py`, `test_gatekeeper.py`, `test_hierarchical_leaf_extractor.py`, `test_native_dpi_patch_tiler.py`, `test_annotate_with_sam2.py` | Unit tests for all archived pipeline components |
| `root_artifacts/` | `yolov8m-seg.pt`, `yolov8n.pt`, `yolov8x-seg.pt`, `weights/`, `runs/`, `segment-anything-2/` | Pretrained and fine-tuned YOLO weight files, training run logs, and SAM 2 submodule |

---

## Replacement Workflow (LeafMachine2)

| Old Component | Replaced By |
|---|---|
| `annotate_with_sam2.py` + `build_artifact_robust_dataset.py` | LM2 uses pre-trained weights; no custom dataset needed |
| `run_dpi_tiler.py` | LM2 internal tiling |
| `train_yolo.py` / `train/` package | LM2 uses `LeafPriority.pt` pretrained detector |
| `run_sahi_inference.py` | LM2 internal inference engine |
| `02_hierarchical_leaf_extractor.py` | LM2 leaf component extraction pipeline |
| `artifact_filter_gatekeeper.py` | LM2 built-in ruler/label/color-chart detection |

**Active LM2 entry points:**
- `LeafMachine2/LeafMachine2.py` — main LM2 runner
- `LM2_Project/configs/lm2_packera_highperf.yaml` — Packera-optimized LM2 configuration
- `scripts/vision/configure_leafmachine2.py` — programmatic LM2 config generator
- `scripts/data_prep/prepare_lm2_dataset.py` — prepares voucher images for LM2 input
- `.venv_LM2/` — dedicated LeafMachine2 Python virtual environment
- `setup_leafmachine2.sh` — LM2 environment setup script

See the project [`README.md`](../../README.md) for the current end-to-end workflow.
