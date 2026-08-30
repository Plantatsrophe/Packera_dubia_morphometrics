import os
import json
import shutil
import logging
from pathlib import Path
from typing import Dict, Optional, Any, Union, Tuple

import torch
from ultralytics import YOLO

from .config import resolve_default_paths
from .dataset import verify_dataset_configuration
from .evaluator import YOLOEvaluator

logger = logging.getLogger("ArtifactRobustYOLO")

def detect_optimal_device_and_batch(requested_batch: Optional[int] = None, imgsz: int = 1024) -> Tuple[str, int]:
    """
    Detects hardware accelerator (CUDA GPU / CPU) and computes safe batch size
    based on available VRAM to prevent Out-Of-Memory (OOM) errors during 1024x1024
    instance segmentation training.
    """
    if not torch.cuda.is_available():
        logger.warning("CUDA is not available. Falling back to CPU execution.")
        selected_batch = requested_batch if requested_batch is not None else 2
        return "cpu", selected_batch

    device_name = torch.cuda.get_device_name(0)
    total_memory_gb = torch.cuda.get_device_properties(0).total_memory / (1024 ** 3)
    torch.backends.cudnn.benchmark = True
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True
    if hasattr(torch, "set_float32_matmul_precision"):
        torch.set_float32_matmul_precision("high")
    logger.info(f"CUDA accelerator detected: {device_name} ({total_memory_gb:.2f} GB VRAM) | cuDNN benchmark: Enabled | TF32: Enabled")

    if requested_batch is not None:
        logger.info(f"Using user-specified batch size: {requested_batch}")
        return "0", requested_batch

    # Heuristic batch size selection based on VRAM capacity at 1024x1024
    if total_memory_gb >= 24.0:
        recommended_batch = 16
    elif total_memory_gb >= 12.0:
        recommended_batch = 8
    elif total_memory_gb >= 8.0:
        recommended_batch = 4
    else:
        recommended_batch = 2

    logger.info(
        f"Auto-selected batch size: {recommended_batch} "
        f"(Tailored for {imgsz}x{imgsz} segmentation on {total_memory_gb:.1f} GB VRAM)"
    )
    return "0", recommended_batch


class RobustYOLOTrainer:
    """
    Orchestrates the fine-tuning, custom loss configuration, botanical augmentation,
    and validation evaluation for YOLOv8-seg.
    """
    def __init__(
        self,
        weights: str = "yolov8m-seg.pt",
        data_config: Union[str, Path] = "data/dataset_config.yaml",
        imgsz: int = 1024,
        batch: Optional[int] = None,
        epochs: int = 150,
        device: Optional[str] = None,
        workers: Optional[int] = None,
        cache: Optional[str] = "disk",
        resume: bool = False,
        compile: bool = False,
        project_root: Optional[Path] = None,
        experiment_name: str = "artifact_robust_yolov8_seg"
    ):
        self.paths = resolve_default_paths(project_root)
        self.data_config_path = Path(data_config) if Path(data_config).is_absolute() else self.paths["root"] / data_config
        self.dataset_info = verify_dataset_configuration(self.data_config_path)

        self.imgsz = imgsz
        self.epochs = epochs
        self.experiment_name = experiment_name
        self.resume = resume
        self.compile = compile
        self.cache = cache

        if workers is not None:
            self.workers = workers
        else:
            self.workers = 0 if os.name == "nt" else 4

        if os.name == "nt" and self.cache == "ram" and self.workers > 0:
            logger.warning(
                f"Windows multiprocessing 'spawn' cannot serialize a 20,000+ image RAM cache "
                f"across {self.workers} worker IPC pipes. "
                f"Automatically adjusting workers=0 for ultra-fast zero-copy main-process RAM caching."
            )
            self.workers = 0

        auto_device, auto_batch = detect_optimal_device_and_batch(batch, imgsz=self.imgsz)
        self.device = device if device is not None else auto_device
        self.batch = batch if batch is not None else auto_batch

        last_ckpt_path = self.paths["root"] / "runs" / "segment" / self.experiment_name / "weights" / "last.pt"
        if self.resume:
            if last_ckpt_path.exists() and last_ckpt_path.stat().st_size > 0:
                logger.info(f"Resuming training from existing checkpoint: {last_ckpt_path}")
                self.weights_path = str(last_ckpt_path)
            else:
                logger.warning(f"Resume requested but '{last_ckpt_path}' was not found. Starting from base weights: {weights}")
                self.resume = False
                self.weights_path = self._resolve_model_weights(weights)
        else:
            self.weights_path = self._resolve_model_weights(weights)

        logger.info(f"Initializing YOLO model with backbone: {self.weights_path} | DataLoader workers: {self.workers} | Resume: {self.resume}")
        self.model = YOLO(str(self.weights_path))
        
        self.evaluator = YOLOEvaluator(
            model=self.model,
            data_config_path=self.data_config_path,
            output_dir=self.paths["output_eval_dir"],
            imgsz=self.imgsz,
            dataset_info=self.dataset_info,
            batch=self.batch,
            device=self.device
        )

    def _resolve_model_weights(self, weights_name: str) -> str:
        local_candidate = self.paths["root"] / weights_name
        if local_candidate.exists() and local_candidate.stat().st_size > 0:
            logger.info(f"Found local model weights at: {local_candidate}")
            return str(local_candidate)

        models_dir_candidate = self.paths["models_dir"] / weights_name
        if models_dir_candidate.exists() and models_dir_candidate.stat().st_size > 0:
            logger.info(f"Found weights in models directory: {models_dir_candidate}")
            return str(models_dir_candidate)

        logger.info(f"Using Ultralytics pre-trained weight identifier: {weights_name}")
        return weights_name

    def build_training_hyperparameters(self) -> Dict[str, Any]:
        hyperparameters = {
            "data": str(self.data_config_path),
            "epochs": self.epochs,
            "imgsz": self.imgsz,
            "batch": self.batch,
            "device": self.device,
            "workers": self.workers,
            "cache": self.cache,
            "amp": True,
            "compile": self.compile,

            "optimizer": "AdamW",
            "lr0": 0.001,
            "lrf": 0.01,
            "cos_lr": True,
            "warmup_epochs": 4.0,
            "warmup_momentum": 0.8,
            "warmup_bias_lr": 0.1,
            "weight_decay": 0.0005,
            "momentum": 0.937,

            "deterministic": False,

            "box": 7.5,
            "cls": 1.5,
            "dfl": 1.5,

            "mosaic": 0.0,
            "mixup": 0.0,
            "copy_paste": 0.15,
            "degrees": 25.0,
            "translate": 0.1,
            "scale": 0.2,
            "shear": 2.0,
            "perspective": 0.0,
            "fliplr": 0.5,
            "flipud": 0.0,
            "hsv_h": 0.015,
            "hsv_s": 0.7,
            "hsv_v": 0.4,

            "project": str(self.paths["root"] / "runs" / "segment"),
            "name": self.experiment_name,
            "exist_ok": True,
            "save": True,
            "save_period": -1,
            "patience": 20,
            "val": True,
            "plots": False,
            "verbose": True
        }
        return hyperparameters

    def train(self) -> Any:
        hyperparams = self.build_training_hyperparameters()
        logger.info("=" * 80)
        logger.info(f"STARTING ARTIFACT-ROBUST YOLOV8-SEG MODEL {'RESUMPTION' if self.resume else 'FINE-TUNING'}")
        logger.info("=" * 80)
        logger.info(f"Hyperparameters Summary:\n{json.dumps({k: str(v) for k, v in hyperparams.items()}, indent=2)}")

        if self.resume:
            logger.info("Resuming training from checkpoint without re-running completed epochs.")
            train_results = self.model.train(resume=True, workers=self.workers, cache=self.cache)
        else:
            train_results = self.model.train(**hyperparams)

        logger.info("Training loop completed successfully.")
        return train_results

    def export_best_checkpoint(self) -> Path:
        target_path = self.paths["best_model_export"]
        run_best_weights = self.paths["root"] / "runs" / "segment" / self.experiment_name / "weights" / "best.pt"

        if run_best_weights.exists():
            shutil.copy2(run_best_weights, target_path)
            logger.info(f"Successfully exported best model checkpoint to: {target_path}")
        else:
            logger.warning(
                f"Run best weights not found at: {run_best_weights}. "
                f"Checking for existing model at target path..."
            )
            if not target_path.exists():
                try:
                    self.model.save(str(target_path))
                    logger.info(f"Saved current model weights directly to: {target_path}")
                except Exception as e:
                    logger.error(f"Failed to export model checkpoint: {e}")

        return target_path
