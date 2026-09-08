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


def ensure_model_weights(
    config: Union[PipelineConfig, Dict[str, Any]],
    force: bool = False,
) -> Path:
    """
    Ensure the PointRend model weights checkpoint exists and is non-empty.

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
            or "models/lm2_packera_pcd_finetuned.pth"
        )
        weights_path = Path(raw_path)
        if not weights_path.is_absolute():
            weights_path = (PROJECT_ROOT / weights_path).resolve()
        weights_url = models_dict.get("pcd_weights_url", "")
        weights_sha256 = models_dict.get("pcd_weights_sha256", "")

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

    if not weights_url:
        raise ValueError(
            f"Model weights checkpoint not found at '{weights_path}' and no remote "
            "download URL ('models.pcd_weights_url') is specified in configuration."
        )

    logger.info(f"Model weights missing or force refresh requested. Acquiring from: {weights_url}")
    return download_file_with_progress(weights_url, weights_path, expected_sha256=weights_sha256)
