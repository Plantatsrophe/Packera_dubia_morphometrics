"""
===============================================================================
Module: staging_utils.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Staging utilities for LeafMachine2 project environments: directory scaffold
    setup, image symlinking, and voucher asset audit manifest generation.
===============================================================================
"""

from __future__ import annotations

import csv
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger("LM2Staging")

STANDARD_IMAGE_EXTENSIONS: Set[str] = {
    ".jpg", ".jpeg", ".png", ".tif", ".tiff", ".bmp"
}


def setup_lm2_directories(lm2_project_root: Path) -> Dict[str, Path]:
    """
    Creates standardized directory scaffolding for a LeafMachine2 project workspace.

    Args:
        lm2_project_root: Root directory of the LeafMachine2 project (e.g. 'LM2_Project').

    Returns:
        Dict[str, Path]: Mapping of directory keys to their resolved paths.
    """
    root = Path(lm2_project_root)
    dirs = {
        "root": root,
        "images": root / "Data" / "images",
        "output": root / "Data" / "output",
        "annotations": root / "Data" / "annotations",
        "configs": root / "configs",
    }
    for k, d in dirs.items():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def stage_voucher_symlinks(
    source_dirs: List[Path],
    target_images_dir: Path,
    use_relative_symlinks: bool = True,
    file_extensions: Optional[Set[str]] = None,
) -> Tuple[int, int, List[Tuple[str, str, str]]]:
    """
    Creates filesystem symlinks from raw voucher source directories into the LM2 staging folder.

    Args:
        source_dirs: List of directory paths containing raw specimen images.
        target_images_dir: Destination staging directory (e.g. LM2_Project/Data/images).
        use_relative_symlinks: If True, creates relative symlinks; otherwise absolute.
        file_extensions: Allowed image extensions (default: STANDARD_IMAGE_EXTENSIONS).

    Returns:
        Tuple: (num_created, num_skipped, manifest_records)
    """
    target_dir = Path(target_images_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    allowed_exts = file_extensions or STANDARD_IMAGE_EXTENSIONS
    manifest_records: List[Tuple[str, str, str]] = []

    created_count = 0
    skipped_count = 0

    for s_dir in source_dirs:
        s_path = Path(s_dir)
        if not s_path.exists():
            logger.warning(f"Source image directory not found: {s_path}")
            continue

        for item in sorted(s_path.iterdir()):
            if item.is_file() and item.suffix.lower() in allowed_exts:
                link_dest = target_dir / item.name

                if link_dest.exists() or link_dest.is_symlink():
                    skipped_count += 1
                    manifest_records.append((item.name, str(item.resolve()), str(link_dest.resolve())))
                    continue

                try:
                    if use_relative_symlinks:
                        rel_target = os.path.relpath(item.resolve(), target_dir.resolve())
                        link_dest.symlink_to(rel_target)
                    else:
                        link_dest.symlink_to(item.resolve())

                    created_count += 1
                    manifest_records.append((item.name, str(item.resolve()), str(link_dest.resolve())))
                except OSError as e:
                    logger.error(f"Failed to create symlink {link_dest} -> {item}: {e}")

    logger.info(
        f"Staged symlinks in {target_dir}: {created_count} created, {skipped_count} existing/skipped."
    )
    return created_count, skipped_count, manifest_records


def write_manifest_csv(
    manifest_records: List[Tuple[str, str, str]],
    manifest_output_path: Path
) -> None:
    """
    Writes staged image manifest mapping table to CSV.
    """
    out_path = Path(manifest_output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["filename", "source_path", "staged_symlink_path", "timestamp"])
        now_str = datetime.now().isoformat()
        for r in manifest_records:
            writer.writerow([r[0], r[1], r[2], now_str])

    logger.info(f"Wrote asset manifest to {out_path} ({len(manifest_records)} records)")
