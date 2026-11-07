"""
===============================================================================
Module: artifact_manager.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Artifact and Model Weights Manager. Automates fault-tolerant acquisition,
    streamed downloading with progress feedback, SHA256 integrity verification,
    and atomic filesystem staging for deep learning model checkpoints.
===============================================================================
"""

from __future__ import annotations

import hashlib
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, Optional, Union

from scripts.core.config import PROJECT_ROOT, PipelineConfig
from scripts.core.logger import setup_logging

logger = setup_logging(name="ArtifactManager")


def compute_sha256(file_path: Union[str, Path], chunk_size: int = 1024 * 1024) -> str:
    """Compute SHA256 checksum of a file in chunks."""
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"Cannot compute hash; file does not exist: {path}")
    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest().lower()


def download_file_with_progress(
    url: str,
    dest_path: Union[str, Path],
    expected_sha256: Optional[str] = None,
    chunk_size: int = 1024 * 64,
) -> Path:
    """
    Stream download from a URL to a destination path atomically with a progress indicator.

    Writes to a temporary .tmp file in the destination folder and atomically
    renames upon successful completion to prevent corrupt partial artifacts.
    """
    target = Path(dest_path).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = target.parent / f"{target.name}.tmp.{os.getpid()}"

    req = urllib.request.Request(
        url,
        headers={"User-Agent": "PackeraDubiaPipeline/1.0"},
    )

    logger.info(f"Connecting to download source: {url}")
    try:
        with urllib.request.urlopen(req) as response, open(tmp_path, "wb") as out_f:
            total_size = int(response.headers.get("Content-Length", 0))
            downloaded = 0

            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                out_f.write(chunk)
                downloaded += len(chunk)

                if sys.stderr.isatty():
                    if total_size > 0:
                        pct = (downloaded / total_size) * 100
                        mb_dl = downloaded / (1024 * 1024)
                        mb_tot = total_size / (1024 * 1024)
                        sys.stderr.write(f"\rDownloading {target.name}: {pct:5.1f}% [{mb_dl:.1f}/{mb_tot:.1f} MB]")
                    else:
                        mb_dl = downloaded / (1024 * 1024)
                        sys.stderr.write(f"\rDownloading {target.name}: {mb_dl:.1f} MB")
                    sys.stderr.flush()

            if sys.stderr.isatty():
                sys.stderr.write("\n")
                sys.stderr.flush()

        if expected_sha256 and expected_sha256.strip():
            logger.info("Verifying SHA256 checksum...")
            actual_sha256 = compute_sha256(tmp_path)
            expected_clean = expected_sha256.strip().lower()
            if actual_sha256 != expected_clean:
                raise ValueError(
                    f"Checksum verification failed for {target.name}!\n"
                    f"  Expected: {expected_clean}\n"
                    f"  Actual:   {actual_sha256}"
                )
            logger.info(f"SHA256 verified successfully: {actual_sha256}")

        # Atomic replacement across same filesystem mount
        tmp_path.replace(target)
        logger.info(f"Successfully saved weights to: {target}")
        return target

    except Exception as exc:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass
        logger.error(f"Download failed for {url}: {exc}")
        raise


_CANONICAL_LM3_HF_URLS: Dict[str, tuple[str, str]] = {
    "lm3_archival_detector__yolo26x_det_1280": (
        "https://huggingface.co/phyloforfun/lm3_archival_detector__yolo26x_det_1280/resolve/a3c09585cb952b6388006684f9e2db460224739a/onnx/model.onnx",
        "f87c866941b4cb9e6bf1d3ff801b49d9d5eb222bc48a40da6f49fd371333bf49",
    ),
    "lm3_plant_detector__yolo26x_det_1280": (
        "https://huggingface.co/phyloforfun/lm3_plant_detector__yolo26x_det_1280/resolve/86fcde82e671ca025a9d3d8be09ebe861599851b/onnx/model.onnx",
        "9f7d251cd650debbca7db59a1dbf87377386ab200256baa6f5d78f1a41a80966",
    ),
    "lm3_landmark_detector__yolo26x_pose_640_canonical": (
        "https://huggingface.co/phyloforfun/lm3_landmark_detector__yolo26x_pose_640_canonical/resolve/736d2c8ce354111171a23961dd77239941a30899/onnx/model.onnx",
        "2e840e0869110acdb48d8e0ce733dec72bf5ec4c4f5bf5d07c916075c2a8b227",
    ),
    "lm3_landmark_detector__yolo26x_pose_640": (
        "https://huggingface.co/phyloforfun/lm3_landmark_detector__yolo26x_pose_640_canonical/resolve/736d2c8ce354111171a23961dd77239941a30899/onnx/model.onnx",
        "2e840e0869110acdb48d8e0ce733dec72bf5ec4c4f5bf5d07c916075c2a8b227",
    ),
    "lm3_leaf_segmenter__yolo26x_seg_1024": (
        "https://huggingface.co/phyloforfun/lm3_leaf_segmenter__yolo26x_seg_1024/resolve/b2a612c3e889146e7563666c094e2c0268df4f8f/onnx/model.onnx",
        "7b8398e8a28594c636db3e9c1d591674cd7e5cc5bb8bd0dd432161b85086205c",
    ),
    "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain": (
        "https://huggingface.co/phyloforfun/lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain/resolve/e98803c53bcdd539d3ce02c39124bf986d542b20/onnx/model.onnx",
        "119d4e0e650d7ece9c55460686bbfbc2461c8aa4df4b8996f54692a98eff3b42",
    ),
}


def ensure_model_weights(
    config: Union[PipelineConfig, Dict[str, Any]],
    force: bool = False,
) -> Path:
    """
    Ensure the LeafMachine3 ONNX model weights checkpoint exists and is non-empty.

    Downloads the weights automatically from configured remote URL if missing
    or if force=True.
    """
    if isinstance(config, PipelineConfig):
        weights_path = Path(config.models.pcd_weights_path)
        weights_url = config.models.pcd_weights_url
        weights_sha256 = config.models.pcd_weights_sha256
    else:
        models_dict = config.get("models", {})
        paths_dict = config.get("paths", {})
        raw_path = (
            models_dict.get("pcd_weights_path")
            or paths_dict.get("model_weights")
            or "models/lm3/lm3_plant_detector__yolo26x_det_1280.onnx"
        )
        weights_path = Path(raw_path)
        if not weights_path.is_absolute():
            weights_path = (PROJECT_ROOT / weights_path).resolve()
        weights_url = models_dict.get("pcd_weights_url", "")
        weights_sha256 = models_dict.get("pcd_weights_sha256", "")

    if not force and not weights_path.exists() and weights_path.stem.startswith("lm3_"):
        try:
            from scripts.vision.configure_leafmachine3 import resolve_lm3_weights
            resolved_cand = resolve_lm3_weights(weights_path.name, explicit_path=weights_path)
            if resolved_cand.exists() and resolved_cand.is_file() and resolved_cand.stat().st_size > 0:
                weights_path = resolved_cand
        except Exception:
            pass

    if not force and weights_path.exists() and weights_path.is_file() and weights_path.stat().st_size > 0:
        if weights_sha256 and weights_sha256.strip():
            actual_sha256 = compute_sha256(weights_path)
            if actual_sha256 != weights_sha256.strip().lower():
                logger.warning(
                    f"Existing model weights at {weights_path} have checksum mismatch! Re-downloading..."
                )
                return download_file_with_progress(weights_url, weights_path, expected_sha256=weights_sha256)
        logger.info(
            f"Pre-existing model weights verified: {weights_path} "
            f"({weights_path.stat().st_size / (1024 * 1024):.2f} MB)"
        )
        return weights_path

    if not weights_url and weights_path.stem in _CANONICAL_LM3_HF_URLS:
        weights_url, default_sha = _CANONICAL_LM3_HF_URLS[weights_path.stem]
        if not weights_sha256:
            weights_sha256 = default_sha

    if not weights_url:
        raise ValueError(
            f"Model weights checkpoint not found at '{weights_path}' and no remote "
            "download URL ('models.pcd_weights_url') is specified in configuration."
        )

    logger.info(f"Model weights missing or force refresh requested. Acquiring from: {weights_url}")
    return download_file_with_progress(weights_url, weights_path, expected_sha256=weights_sha256)
