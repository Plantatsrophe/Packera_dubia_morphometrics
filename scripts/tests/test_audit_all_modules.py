import importlib
import os
import sys
from pathlib import Path
import pytest

CORE_MODULES = [
    "leafmachine2.machine.machine",
    "leafmachine2.machine.machine_fast",
    "leafmachine2.machine.machine_specimencrop",
    "leafmachine2.machine.machine_detect_phenology",
    "leafmachine2.machine.machine_censor_components",
    "leafmachine2.machine.utils_ruler",
    "leafmachine2.machine.binarize_image_ML",
    "leafmachine2.machine.build_custom_overlay",
    "leafmachine2.machine.save_data",
    "leafmachine2.machine.general_utils",
    "leafmachine2.machine.directory_structure",
    "leafmachine2.machine.data_project",
    "leafmachine2.machine.LM2_logger",
    "leafmachine2.machine.fetch_data",
    "leafmachine2.machine.visualize_EFDs",
    "leafmachine2.machine.ruler",
    "leafmachine2.component_detector.component_detector",
    "leafmachine2.component_detector.detect",
    "leafmachine2.component_detector.landmark_processing",
    "leafmachine2.component_detector.armature_processing",
    "leafmachine2.component_detector.models.common",
    "leafmachine2.component_detector.models.yolo",
    "leafmachine2.component_detector.models.experimental",
    "leafmachine2.component_detector.utils.general",
    "leafmachine2.component_detector.utils.metrics",
    "leafmachine2.component_detector.utils.datasets",
    "leafmachine2.component_detector.utils.dataloaders",
    "leafmachine2.segmentation.detectron2.detector",
    "leafmachine2.segmentation.detectron2.predictor_leaf",
    "leafmachine2.segmentation.detectron2.segment_leaves",
    "leafmachine2.segmentation.detectron2.measure_leaf_segmentation",
    "leafmachine2.landmarks.detect",
    "leafmachine2.landmarks.locate",
    "leafmachine2.landmarks.models.utils",
    "leafmachine2.machine.DocEnTR.models.binae",
    "leafmachine2.machine.DocEnTR.predict",
]

@pytest.mark.parametrize("mod_name", CORE_MODULES)
def test_core_module_importable(mod_name):
    # Ensure paths are set
    sys.path.insert(0, os.path.abspath("."))
    sys.path.insert(0, os.path.abspath("LeafMachine2"))
    sys.path.insert(0, os.path.abspath("LeafMachine2/leafmachine2"))
    
    mod = importlib.import_module(mod_name)
    assert mod is not None
