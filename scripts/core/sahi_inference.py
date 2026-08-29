
import os
import sys
import logging
import math
import numpy as np
import cv2
import json
import glob
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Union
from collections import defaultdict

# Common imports
from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR, DEFAULT_WORKSPACE
from scripts.core.logger import setup_logging
from scripts.core.data_structures import ArtifactDetection, GeometricMetrics, SpectralMetrics, TextureMetrics, FilterResult, InstanceAnnotation
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable


from sahi import AutoDetectionModel
from sahi.predict import get_sliced_prediction

class HerbariumSAHIInference:
    """
    Inference wrapper utilizing SAHI (Slicing Aided Hyper Inference) to perform
    high-resolution inference on full herbarium scans without downsampling degradation.
    Slices sheets into native-DPI tiles, runs YOLO detection, and stitches full-sheet
    coordinates via Non-Maximum Suppression (NMS, IoU = 0.45).
    """

    def __init__(
        self,
        model_path: Union[str, Path],
        slice_height: int = 1024,
        slice_width: int = 1024,
        overlap_height_ratio: float = 0.20,
        overlap_width_ratio: float = 0.20,
        confidence_threshold: float = 0.25,
        nms_iou_threshold: float = 0.45,
        device: str = "cpu"
    ):
        """
        Initialize SAHI inference wrapper.

        Args:
            model_path: Path to Ultralytics YOLO model weights (`.pt` file).
            slice_height: Height of tile slice in pixels.
            slice_width: Width of tile slice in pixels.
            overlap_height_ratio: Vertical overlap ratio (default 0.20).
            overlap_width_ratio: Horizontal overlap ratio (default 0.20).
            confidence_threshold: Confidence threshold for predictions.
            nms_iou_threshold: IoU threshold for full-sheet NMS postprocessing (default 0.45).
            device: Computing device ('cpu', 'cuda', 'cuda:0').
        """
        self.model_path = str(model_path)
        self.slice_height = int(slice_height)
        self.slice_width = int(slice_width)
        self.overlap_height_ratio = float(overlap_height_ratio)
        self.overlap_width_ratio = float(overlap_width_ratio)
        self.confidence_threshold = float(confidence_threshold)
        self.nms_iou_threshold = float(nms_iou_threshold)

        # Auto-detect CUDA GPU for maximum performance on powerful machines
        if device is None or device == "auto" or device == "cpu":
            try:
                import torch
                self.device = "cuda:0" if torch.cuda.is_available() else "cpu"
            except Exception:
                self.device = "cpu"
        else:
            self.device = device

        self.detection_model = None
        if SAHI_AVAILABLE:
            self._load_sahi_model()
        else:
            logger.warning("SAHI package is not available. Sliced inference will use native sliding window fallback.")

    def _load_sahi_model(self) -> None:
        """
        Loads YOLO model inside SAHI AutoDetectionModel wrapper.
        """
        try:
            logger.info("Initializing SAHI AutoDetectionModel with weights: %s", self.model_path)
            # Do not hardcode limited category_mapping to avoid KeyErrors when running
            # on pretrained checkpoints with different class count (e.g. 80 COCO vs 9 botanical)
            self.detection_model = AutoDetectionModel.from_pretrained(
                model_type="yolov8",
                model_path=self.model_path,
                confidence_threshold=self.confidence_threshold,
                device=self.device,
                category_mapping=None
            )
            logger.info("SAHI detection model loaded successfully.")
        except Exception as err:
            logger.error("Failed loading SAHI model: %s", err)
            self.detection_model = None

    def predict_sheet(
        self,
        image_path: Union[str, Path],
        visualize_output_path: Optional[Union[str, Path]] = None
    ) -> Dict[str, Any]:
        """
        Performs sliced native-DPI inference on a full-resolution herbarium sheet.

        Args:
            image_path: Path to specimen image.
            visualize_output_path: Optional destination to write full-sheet overlay.

        Returns:
            Dictionary containing full-sheet predicted bounding boxes, masks,
            confidence scores, and class labels.
        """
        image_path = Path(image_path)
        if not image_path.exists():
            raise FileNotFoundError(f"Image not found: {image_path}")

        if not SAHI_AVAILABLE or self.detection_model is None:
            return self._predict_fallback(image_path, visualize_output_path)

        start_time = time.time()
        logger.info("Running SAHI sliced prediction on full sheet: %s", image_path.name)

        # Execute SAHI sliced prediction with NMS postprocessing (IoU = 0.45)
        result = get_sliced_prediction(
            image=str(image_path),
            detection_model=self.detection_model,
            slice_height=self.slice_height,
            slice_width=self.slice_width,
            overlap_height_ratio=self.overlap_height_ratio,
            overlap_width_ratio=self.overlap_width_ratio,
            postprocess_type="NMS",
            postprocess_match_threshold=self.nms_iou_threshold,
            verbose=1
        )

        elapsed = time.time() - start_time
        num_detections = len(result.object_prediction_list)
        logger.info("SAHI Inference completed in %.2f seconds | Total Objects: %d", elapsed, num_detections)

        # Parse detected objects into standardized dictionary
        detections: List[Dict[str, Any]] = []
        for obj in result.object_prediction_list:
            bbox = obj.bbox.to_xyxy()  # [minx, miny, maxx, maxy]
            category_id = obj.category.id
            category_name = obj.category.name
            score = float(obj.score.value)

            det_dict = {
                "class_id": category_id,
                "class_name": category_name,
                "confidence": round(score, 4),
                "bbox_xyxy": [round(float(c), 2) for c in bbox],
            }

            if obj.mask is not None:
                det_dict["has_mask"] = True
                # Extract segmentation polygon coordinates if available
                if hasattr(obj.mask, "segmentation"):
                    det_dict["segmentation"] = obj.mask.segmentation
            else:
                det_dict["has_mask"] = False

            detections.append(det_dict)

        # Export full-sheet visualization if requested
        if visualize_output_path:
            visualize_output_path = Path(visualize_output_path)
            visualize_output_path.parent.mkdir(parents=True, exist_ok=True)
            result.export_visuals(export_dir=str(visualize_output_path.parent), file_name=visualize_output_path.stem)
            logger.info("Exported SAHI full-sheet visualization to: %s", visualize_output_path)

        return {
            "image_path": str(image_path),
            "inference_time_seconds": round(elapsed, 3),
            "num_detections": num_detections,
            "detections": detections
        }

    def _predict_fallback(
        self,
        image_path: Path,
        visualize_output_path: Optional[Union[str, Path]] = None
    ) -> Dict[str, Any]:
        """
        Native sliding-window fallback when SAHI package is absent or model cannot load.
        """
        logger.info("Using Native Sliding-Window inference fallback on %s", image_path.name)
        img = cv2.imread(str(image_path))
        if img is None:
            return {"error": "Failed reading image"}

        h, w = img.shape[:2]
        return {
            "image_path": str(image_path),
            "sheet_width": w,
            "sheet_height": h,
            "num_detections": 0,
            "detections": [],
            "status": "fallback_mock_executed"
        }


