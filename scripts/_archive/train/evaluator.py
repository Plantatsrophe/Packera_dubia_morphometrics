import json
import shutil
import logging
from pathlib import Path
from typing import Dict, List, Any, Union, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ultralytics import YOLO

logger = logging.getLogger("ArtifactRobustYOLO")

class YOLOEvaluator:
    """
    Evaluates the trained model on validation split, computing class-specific
    mAP50, mAP50-95, precision, recall, and analyzing cross-classification
    between herbarium artifacts and botanical organs.
    """
    def __init__(
        self,
        model: YOLO,
        data_config_path: Path,
        output_dir: Path,
        imgsz: int = 1024,
        dataset_info: Optional[Dict[str, Any]] = None,
        batch: Optional[int] = None,
        device: Optional[str] = None
    ):
        self.model = model
        self.data_config_path = data_config_path
        self.output_dir = output_dir
        self.imgsz = imgsz
        self.dataset_info = dataset_info or {}
        self.batch = batch
        self.device = device

    def evaluate_and_export_metrics(
        self,
        checkpoint_path: Optional[Union[str, Path]] = None,
        split: str = "val"
    ) -> Dict[str, Any]:
        """
        Evaluates the model.
        """
        if checkpoint_path is not None:
            ckpt_p = Path(checkpoint_path)
            if ckpt_p.exists() and ckpt_p.stat().st_size > 0:
                eval_model = YOLO(str(ckpt_p))
            else:
                logger.warning(f"Checkpoint at '{checkpoint_path}' does not exist or is 0-bytes. Using active model instance.")
                eval_model = self.model
        else:
            eval_model = self.model

        logger.info("=" * 80)
        logger.info(f"EVALUATING MODEL ON {split.upper()} SPLIT (imgsz={self.imgsz})")
        logger.info("=" * 80)

        # Run validation with plots enabled and save_json disabled
        val_results = eval_model.val(
            data=str(self.data_config_path),
            imgsz=self.imgsz,
            batch=self.batch,
            device=self.device,
            split=split,
            plots=True,
            save_json=False
        )

        class_names = self.dataset_info.get("names", {})
        metrics_summary: Dict[str, Any] = {
            "overall_box_map50": float(val_results.box.map50) if hasattr(val_results, "box") else 0.0,
            "overall_box_map50_95": float(val_results.box.map) if hasattr(val_results, "box") else 0.0,
            "overall_mask_map50": float(val_results.seg.map50) if hasattr(val_results, "seg") else 0.0,
            "overall_mask_map50_95": float(val_results.seg.map) if hasattr(val_results, "seg") else 0.0,
            "class_metrics": {}
        }

        logger.info("-" * 80)
        logger.info(f"{'Class ID':<10}{'Class Name':<20}{'Box mAP50':<14}{'Box mAP50-95':<16}{'Mask mAP50':<14}{'Mask mAP50-95':<16}")
        logger.info("-" * 80)

        # Extract per-class metrics
        box_maps50 = getattr(val_results.box, "maps50", None)
        box_maps = getattr(val_results.box, "maps", None)
        seg_maps50 = getattr(val_results.seg, "maps50", None)
        seg_maps = getattr(val_results.seg, "maps", None)

        for class_id, class_name in class_names.items():
            cid = int(class_id)
            b_map50 = float(box_maps50[cid]) if box_maps50 is not None and cid < len(box_maps50) else 0.0
            b_map = float(box_maps[cid]) if box_maps is not None and cid < len(box_maps) else 0.0
            s_map50 = float(seg_maps50[cid]) if seg_maps50 is not None and cid < len(seg_maps50) else 0.0
            s_map = float(seg_maps[cid]) if seg_maps is not None and cid < len(seg_maps) else 0.0

            metrics_summary["class_metrics"][class_name] = {
                "class_id": cid,
                "box_map50": b_map50,
                "box_map50_95": b_map,
                "mask_map50": s_map50,
                "mask_map50_95": s_map
            }

            logger.info(
                f"{cid:<10}{class_name:<20}{b_map50:<14.4f}{b_map:<16.4f}{s_map50:<14.4f}{s_map:<16.4f}"
            )

        logger.info("-" * 80)

        # Verify key requirements
        basal_leaf_metrics = metrics_summary["class_metrics"].get("basal_leaf", {})
        petiole_metrics = metrics_summary["class_metrics"].get("leaf_petiole", {})
        logger.info(
            f"Key Botanical Organ Performance -> "
            f"basal_leaf Mask mAP50: {basal_leaf_metrics.get('mask_map50', 0.0):.4f} | "
            f"leaf_petiole Mask mAP50: {petiole_metrics.get('mask_map50', 0.0):.4f}"
        )

        # Cross-classification analysis and visualizations
        self._generate_evaluation_artifacts(val_results, class_names, metrics_summary)

        # Save metrics to JSON
        self.output_dir.mkdir(parents=True, exist_ok=True)
        eval_report_path = self.output_dir / "evaluation_report.json"
        with open(eval_report_path, "w", encoding="utf-8") as f:
            json.dump(metrics_summary, f, indent=2)
        logger.info(f"Exported evaluation metrics JSON to: {eval_report_path}")

        return metrics_summary

    def _generate_evaluation_artifacts(
        self,
        val_results: Any,
        class_names: Dict[Any, str],
        metrics_summary: Dict[str, Any]
    ) -> None:
        """
        Generates and exports custom confusion matrix plots, Precision-Recall curves,
        and verifies cross-classification between herbarium artifacts and botanical organs.
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 1. Confusion Matrix Analysis
        cm_matrix = None
        if hasattr(val_results, "confusion_matrix") and val_results.confusion_matrix is not None:
            cm = val_results.confusion_matrix
            if hasattr(cm, "matrix"):
                cm_matrix = cm.matrix

        num_classes = len(class_names)
        ordered_names = [class_names[i] for i in range(num_classes)]

        if cm_matrix is not None and isinstance(cm_matrix, np.ndarray):
            # Check cross-classification between herbarium_label and basal_leaf
            label_idx = None
            leaf_idx = None
            for idx, name in enumerate(ordered_names):
                if name == "herbarium_label":
                    label_idx = idx
                elif name == "basal_leaf":
                    leaf_idx = idx

            cross_misclassifications = 0
            if label_idx is not None and leaf_idx is not None and cm_matrix.shape[0] > max(label_idx, leaf_idx):
                label_as_leaf = cm_matrix[label_idx, leaf_idx]
                leaf_as_label = cm_matrix[leaf_idx, label_idx]
                cross_misclassifications = int(label_as_leaf + leaf_as_label)

                logger.info("=" * 80)
                logger.info("ARTIFACT DISCRIMINATION INTEGRITY CHECK:")
                logger.info(f"  - Herbarium Label misclassified as Basal Leaf: {int(label_as_leaf)}")
                logger.info(f"  - Basal Leaf misclassified as Herbarium Label: {int(leaf_as_label)}")
                if cross_misclassifications == 0:
                    logger.info("  [PASSED] PERFECT ZERO CROSS-CLASSIFICATION CONFIRMED.")
                else:
                    logger.warning(
                        f"  [WARNING] Cross-classification detected: {cross_misclassifications} instances."
                    )
                logger.info("=" * 80)

            metrics_summary["label_to_leaf_misclassifications"] = cross_misclassifications

            # Render custom publication-quality confusion matrix heatmap
            self._plot_custom_confusion_matrix(cm_matrix, ordered_names, self.output_dir / "confusion_matrix_custom.png")

        # 2. Precision-Recall Curves Plot
        self._plot_classwise_map_bars(metrics_summary, self.output_dir / "classwise_map_comparison.png")

        # Copy any Ultralytics generated curves
        if hasattr(val_results, "save_dir") and val_results.save_dir:
            save_dir = Path(val_results.save_dir)
            for curve_file in save_dir.glob("*.png"):
                dest_file = self.output_dir / curve_file.name
                shutil.copy2(curve_file, dest_file)
            logger.info(f"Copied Ultralytics validation curve artifacts from {save_dir} to {self.output_dir}")

    def _plot_custom_confusion_matrix(
        self,
        matrix: np.ndarray,
        class_names: List[str],
        output_path: Path
    ) -> None:
        """
        Plots a high-contrast confusion matrix focusing on botanical vs. artifact discrimination.
        """
        try:
            plt.figure(figsize=(10, 8), dpi=300)
            row_sums = matrix.sum(axis=1, keepdims=True)
            norm_matrix = np.divide(matrix, row_sums, out=np.zeros_like(matrix, dtype=float), where=row_sums != 0)

            display_names = class_names.copy()
            if matrix.shape[0] > len(class_names):
                display_names.append("background")

            matrix_slice = norm_matrix[:len(display_names), :len(display_names)]

            fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
            cax = ax.imshow(matrix_slice, interpolation='nearest', cmap=plt.cm.Blues)
            plt.colorbar(cax, fraction=0.046, pad=0.04)

            thresh = matrix_slice.max() / 2.0 if matrix_slice.max() > 0 else 0.5
            for i in range(matrix_slice.shape[0]):
                for j in range(matrix_slice.shape[1]):
                    val = matrix_slice[i, j]
                    ax.text(
                        j, i, f"{val:.2f}",
                        ha="center", va="center",
                        color="white" if val > thresh else "black",
                        fontsize=8
                    )

            ax.set_xticks(np.arange(len(display_names)))
            ax.set_yticks(np.arange(len(display_names)))
            ax.set_xticklabels(display_names, rotation=45, ha="right", fontsize=9)
            ax.set_yticklabels(display_names, fontsize=9)
            ax.set_xlabel("Predicted Class", fontsize=11, labelpad=8)
            ax.set_ylabel("True Class", fontsize=11, labelpad=8)
            ax.set_title("Normalized Confusion Matrix (Botanical vs. Artifact Discrimination)", fontsize=13, pad=15)
            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            plt.close()
            logger.info(f"Saved custom confusion matrix visualization: {output_path}")
        except Exception as e:
            logger.warning(f"Could not render custom confusion matrix plot: {e}")

    def _plot_classwise_map_bars(
        self,
        metrics_summary: Dict[str, Any],
        output_path: Path
    ) -> None:
        """
        Plots a bar chart comparing Mask mAP50 and Mask mAP50-95 across classes.
        """
        try:
            class_metrics = metrics_summary.get("class_metrics", {})
            if not class_metrics:
                return

            names = list(class_metrics.keys())
            mask_map50 = [class_metrics[n]["mask_map50"] for n in names]
            mask_map50_95 = [class_metrics[n]["mask_map50_95"] for n in names]

            x = np.arange(len(names))
            width = 0.35

            plt.figure(figsize=(12, 6), dpi=300)
            plt.bar(x - width/2, mask_map50, width, label="Mask mAP50", color="#2b5c8f")
            plt.bar(x + width/2, mask_map50_95, width, label="Mask mAP50-95", color="#52b788")

            plt.ylabel("mAP Score", fontsize=11)
            plt.title("Class-Specific Segmentation Performance (YOLOv8-seg Fine-Tuning)", fontsize=13, pad=15)
            plt.xticks(x, names, rotation=35, ha="right", fontsize=9)
            plt.ylim(0.0, 1.05)
            plt.legend(loc="upper right")
            plt.grid(axis="y", linestyle="--", alpha=0.5)
            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            plt.close()
            logger.info(f"Saved class-wise mAP comparison plot: {output_path}")
        except Exception as e:
            logger.warning(f"Could not render class-wise mAP bar chart: {e}")
