"""
===============================================================================
Module: lm2_config_template.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Static configuration templates and default hyperparameter constants for
    LeafMachine2 (LM2) inference and Plant Component Detector (PCD) evaluation.
===============================================================================
"""

from __future__ import annotations

from typing import Any, Dict

# Packera-specific defaults for dense rosettes & PointRend
DEFAULT_PCD_WEIGHTS_NAME = "LeafPriority.pt"
DEFAULT_SEG_MODEL_NAME = "Packera_LeafPriority"
DEFAULT_PCD_CONFIDENCE = 0.20
DEFAULT_NMS_THRESH = 0.75
DEFAULT_SUBDIVISION_STEPS = 5
DEFAULT_MIN_LEAF_AREA = 500
DEFAULT_SEG_CONFIDENCE = 0.70
DEFAULT_EFD_ORDER = 40

# Canonical LeafMachine2 template structure matching LeafMachine2.yaml schema
DEFAULT_PACKERA_CONFIG_TEMPLATE: Dict[str, Any] = {
    "leafmachine": {
        "project": {
            "dir_output": "LM2_Project/Data/output",
            "run_name": "Packera_dubia_LM2",
            "dir_images_local": "LM2_Project/Data/images",
            "use_full_component_detector": True,
            "use_single_component_detector": False,
            "use_segmentation": True,
            "use_landmark_detector": False,
            "use_ruler": True,
            "use_archival_component_detector": False,
            "save_overlay_to_pdf": False,
            "save_overlay_to_jpg": True,
            "save_individual_component_images": False,
            "save_individual_ruler_images": False,
            "save_individual_segmentation_images": False,
            "save_individual_landmark_images": False,
            "save_individual_archival_images": False,
            "save_json_data": True,
            "save_csv_data": True,
            "save_segmentation_masks": True,
            "save_ruler_masks": True,
        },
        "plant_component_detector": {
            "detector_type": "Plant_Component_Detector",
            "detector_version": "PCD_v2",
            "detector_weights": "LeafPriority.pt",
            "minimum_confidence_threshold": 0.20,
            "PCD_confidence": 0.20,
            "iou_threshold": 0.75,
            "device": "cuda",
            "batch_size": 50,
            "num_workers": 8,
        },
        "leaf_segmentation": {
            "segmentation_model": "Packera_LeafPriority",
            "segmentation_type": "Detectron2_PointRend",
            "minimum_confidence_threshold": 0.70,
            "NMS_thresh": 0.75,
            "pointrend_subdivision_steps": 5,
            "min_leaf_area": 500,
            "minimum_leaf_area_px": 500,
            "calculate_elliptic_fourier_descriptors": True,
            "elliptic_fourier_descriptor_order": 40,
            "efd_order": 40,
            "device": "cuda",
            "num_workers_seg": 8,
        },
        "overlay": {
            "save_overlay_to_jpg": True,
            "save_overlay_to_pdf": False,
            "line_width": 2,
            "alpha": 0.5,
        },
        "ruler": {
            "ruler_type": "Standard_Scale_Bar",
            "detection_confidence": 0.50,
            "num_workers_ruler": 8,
        },
        "logging": {
            "show_ruler_progress": False,
            "show_segmentation_progress": False,
            "show_landmark_progress": False,
            "show_archival_progress": False,
            "verbose": False,
        },
    }
}
