import os
import argparse
from pathlib import Path
from typing import Dict, Optional

# Standardized 7-class botanical ontology taxonomy
STANDARD_CLASS_MAPPING: Dict[int, str] = {
    0: "basal_leaf_blade",     # Laminar portion of basal leaves
    1: "leaf_petiole",         # Narrow petiole stalk connecting caudex to blade
    2: "cauline_leaf",         # Sessile/lyrately-pinnatifid leaves on flowering stalk
    3: "cauline_stem",         # Main vertical flowering stalk / scape
    4: "root_rhizome",         # Dark fibrous subterranean roots and caudex
    5: "basal_rosette_clump",  # Dense overlapping basal rosette crown
    6: "capitulum",            # Inflorescence head / involucre
}

def resolve_default_paths(project_root: Optional[Path] = None) -> Dict[str, Path]:
    """
    Resolves standard project paths relative to workspace root.

    Args:
        project_root: Optional custom project root path.

    Returns:
        Dictionary mapping path identifiers to absolute Path objects.
    """
    if project_root is None:
        # Default to parent of scripts/train/ or current working directory
        current_script_dir = Path(__file__).resolve().parent
        if current_script_dir.name == "train" and current_script_dir.parent.name == "scripts":
            project_root = current_script_dir.parent.parent
        else:
            project_root = Path.cwd()

    tiled_dataset_dir = project_root / "data" / "tiled_dataset"
    tiled_config = project_root / "data" / "tiled_dataset_config.yaml"
    legacy_config = project_root / "data" / "dataset_config.yaml"

    # Default to standard dataset config if present, else tiled dataset config
    default_config = legacy_config if legacy_config.exists() else tiled_config

    paths = {
        "root": project_root,
        "data_config": default_config,
        "tiled_config": tiled_config,
        "legacy_config": legacy_config,
        "tiled_dir": tiled_dataset_dir,
        "models_dir": project_root / "models",
        "output_eval_dir": project_root / "outputs" / "training_evaluation",
        "default_weights": project_root / "yolov8x-seg.pt",
        "fallback_weights": project_root / "yolov8m-seg.pt",
        "best_model_export": project_root / "models" / "yolov8_leaf_best.pt"
    }

    # Ensure required output directories exist
    paths["models_dir"].mkdir(parents=True, exist_ok=True)
    paths["output_eval_dir"].mkdir(parents=True, exist_ok=True)

    return paths


def parse_arguments() -> argparse.Namespace:
    """
    Parses command-line arguments for the YOLOv8 fine-tuning script.
    """
    parser = argparse.ArgumentParser(
        description="Fine-tune an artifact-robust YOLOv8-seg model for botanical organ segmentation on sliced native-DPI tiles or full voucher sheets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--weights",
        type=str,
        default="yolov8m-seg.pt",
        help="Initial model backbone weights (e.g. 'yolov8x-seg.pt', 'yolov8m-seg.pt', or path to .pt)."
    )
    parser.add_argument(
        "--data",
        type=str,
        default=None,
        help="Path to YOLO dataset YAML configuration file (default: auto-detects data/tiled_dataset_config.yaml if tiled dataset is present, otherwise data/dataset_config.yaml)."
    )
    parser.add_argument(
        "--split-tiled-dataset",
        action="store_true",
        help="Explicitly (re-)partition data/tiled_dataset into train/val/test splits grouped by parent specimen sheet."
    )
    parser.add_argument(
        "--train-ratio",
        type=float,
        default=0.70,
        help="Fraction of specimen sheets assigned to training split when partitioning tiled dataset."
    )
    parser.add_argument(
        "--val-ratio",
        type=float,
        default=0.15,
        help="Fraction of specimen sheets assigned to validation split when partitioning tiled dataset."
    )
    parser.add_argument(
        "--test-ratio",
        type=float,
        default=0.15,
        help="Fraction of specimen sheets assigned to testing split when partitioning tiled dataset."
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=150,
        help="Total number of training epochs."
    )
    parser.add_argument(
        "--batch",
        type=int,
        default=None,
        help="Training batch size (None for auto-detection based on GPU VRAM)."
    )
    parser.add_argument(
        "--imgsz",
        type=int,
        default=1024,
        help="Input image resolution in pixels (default 1024 for high-resolution botanical tiles)."
    )
    parser.add_argument(
        "--device",
        type=str,
        default=None,
        help="Computing device to use (e.g. '0', '0,1', 'cpu'). Auto-detected if None."
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=16,
        help="Number of DataLoader worker subprocesses (0 recommended on Windows/Google Drive with RAM cache)."
    )
    parser.add_argument(
        "--cache",
        type=str,
        default="disk",
        choices=["ram", "disk", "none"],
        help="Image caching strategy. 'disk' is recommended on Windows to prevent memory exhaustion with high workers."
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume training from existing runs/segment/<name>/weights/last.pt without restarting."
    )
    parser.add_argument(
        "--compile",
        action="store_true",
        help="Enable PyTorch 2.0+ torch.compile(model) with TorchInductor/Triton backend for accelerated training throughput."
    )
    parser.add_argument(
        "--name",
        type=str,
        default="artifact_robust_yolov8_seg",
        help="Experiment run name for checkpoint outputs."
    )
    parser.add_argument(
        "--eval-only",
        action="store_true",
        help="Skip training loop and evaluate existing checkpoint."
    )
    parser.add_argument(
        "--checkpoint",
        type=str,
        default=None,
        help="Path to specific checkpoint to evaluate when --eval-only is set."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run a 1-epoch lightweight trial to verify the training and evaluation pipeline."
    )

    return parser.parse_args()
