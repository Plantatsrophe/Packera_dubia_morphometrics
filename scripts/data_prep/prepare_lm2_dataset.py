#!/usr/bin/env python3
"""
===============================================================================
Script: prepare_lm2_dataset.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Pre-processor for LeafMachine2 (LM2). Sets up the standardized LM2_Project/Data
    directory structure and creates safe symbolic links (symlinks) to raw herbarium
    sheet images without duplicating large raster files on disk.

    Features:
      - Validates and filters out non-image files (.DS_Store, metadata, CSV, text).
      - Checks for 0-byte corrupt files and validates image file extensions.
      - Handles name collisions and existing symlinks cleanly.
      - Supports --dry-run mode, batch limits (--limit), and optional image verification.
      - Generates an audit manifest CSV of all linked assets.
===============================================================================
"""

import argparse
import csv
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Supported raster image extensions for herbarium processing
VALID_IMAGE_EXTENSIONS: Set[str] = {
    ".jpg",
    ".jpeg",
    ".png",
    ".tif",
    ".tiff",
    ".webp",
    ".bmp",
}

def setup_logging(verbose: bool = False) -> logging.Logger:
    """Configure structured logging output."""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger("LM2_Preprocessor")


def verify_image_file(file_path: Path) -> bool:
    """
    Perform a lightweight header verification using PIL to ensure the file
    is a readable, valid image format.
    """
    try:
        from PIL import Image
        with Image.open(file_path) as img:
            img.verify()
        return True
    except Exception:
        return False


def create_lm2_directories(
    lm2_root: Path,
    images_dir_name: str = "images",
    output_dir_name: str = "output",
    configs_dir_name: str = "configs",
    dry_run: bool = False,
    logger: Optional[logging.Logger] = None,
) -> Tuple[Path, Path, Path]:
    """
    Create the standardized LM2_Project directory structure.
    
    Structure:
      <lm2_root>/
        Data/
          <images_dir_name>/   <- Symlinks to raw herbarium images
          <output_dir_name>/   <- Destination for LM2 analysis output
        <configs_dir_name>/    <- Project-specific configuration files
    """
    data_dir = lm2_root / "Data"
    images_dir = data_dir / images_dir_name
    output_dir = data_dir / output_dir_name
    configs_dir = lm2_root / configs_dir_name

    dirs_to_create = [lm2_root, data_dir, images_dir, output_dir, configs_dir]

    for d in dirs_to_create:
        if not dry_run:
            d.mkdir(parents=True, exist_ok=True)
        if logger:
            logger.debug(f"Ensured directory exists: {d}")

    return images_dir, output_dir, configs_dir


def collect_raw_images(
    input_dirs: List[Path],
    verify_images: bool = False,
    logger: Optional[logging.Logger] = None,
) -> Tuple[List[Path], List[Tuple[Path, str]]]:
    """
    Scan provided input directories and collect valid herbarium sheet images.
    Returns:
      (valid_images, filtered_out_files)
    """
    valid_images: List[Path] = []
    filtered_out: List[Tuple[Path, str]] = []

    for input_dir in input_dirs:
        if not input_dir.exists():
            if logger:
                logger.warning(f"Input directory does not exist: {input_dir}")
            continue

        if logger:
            logger.info(f"Scanning input directory: {input_dir.resolve()}")

        # Iterate through directory items
        for entry in input_dir.iterdir():
            # Skip subdirectories or handle recursively if needed
            if not entry.is_file():
                continue

            # Check for hidden files (e.g. .DS_Store, .gitkeep)
            if entry.name.startswith("."):
                filtered_out.append((entry, "Hidden file"))
                continue

            # Check file size (filter 0-byte corrupt files)
            try:
                if entry.stat().st_size == 0:
                    filtered_out.append((entry, "Zero-byte file"))
                    continue
            except OSError as e:
                filtered_out.append((entry, f"Stat error: {e}"))
                continue

            # Check image extension
            ext = entry.suffix.lower()
            if ext not in VALID_IMAGE_EXTENSIONS:
                filtered_out.append((entry, f"Unsupported extension: '{ext}'"))
                continue

            # Optional deep verification
            if verify_images:
                if not verify_image_file(entry):
                    filtered_out.append((entry, "Corrupt or unreadable image header"))
                    continue

            valid_images.append(entry.resolve())

    return valid_images, filtered_out


def create_symlinks(
    source_images: List[Path],
    target_dir: Path,
    relative: bool = False,
    overwrite: bool = False,
    limit: Optional[int] = None,
    dry_run: bool = False,
    logger: Optional[logging.Logger] = None,
) -> Tuple[List[Dict[str, str]], int, int, int]:
    """
    Safely create symbolic links pointing from target_dir to source_images.

    Returns:
      (manifest_records, created_count, skipped_count, error_count)
    """
    manifest: List[Dict[str, str]] = []
    created_count = 0
    skipped_count = 0
    error_count = 0

    if limit and limit > 0:
        source_images = source_images[:limit]
        if logger:
            logger.info(f"Limiting to first {limit} images as requested.")

    seen_names: Set[str] = set()

    for src_path in source_images:
        dest_filename = src_path.name

        # Handle filename collisions from multiple source directories
        if dest_filename in seen_names:
            if logger:
                logger.warning(f"Filename collision detected for '{dest_filename}'. Skipping duplicate.")
            skipped_count += 1
            continue

        seen_names.add(dest_filename)
        dest_link_path = target_dir / dest_filename

        # Determine link target (relative or absolute)
        if relative:
            try:
                link_target = os.path.relpath(src_path, target_dir)
            except ValueError:
                link_target = str(src_path)
        else:
            link_target = str(src_path)

        link_status = "PENDING"
        
        # Check existing destination
        if dest_link_path.is_symlink() or dest_link_path.exists():
            if overwrite:
                if not dry_run:
                    try:
                        dest_link_path.unlink()
                        os.symlink(link_target, dest_link_path)
                        link_status = "OVERWRITTEN"
                        created_count += 1
                    except OSError as err:
                        if logger:
                            logger.error(f"Failed to overwrite symlink {dest_link_path}: {err}")
                        link_status = f"ERROR: {err}"
                        error_count += 1
                else:
                    link_status = "DRY_RUN_OVERWRITE"
                    created_count += 1
            else:
                link_status = "SKIPPED_EXISTS"
                skipped_count += 1
        else:
            if not dry_run:
                try:
                    os.symlink(link_target, dest_link_path)
                    link_status = "CREATED"
                    created_count += 1
                except OSError as err:
                    if logger:
                        logger.error(f"Failed to create symlink {dest_link_path} -> {link_target}: {err}")
                    link_status = f"ERROR: {err}"
                    error_count += 1
            else:
                link_status = "DRY_RUN_CREATED"
                created_count += 1

        manifest.append({
            "filename": dest_filename,
            "source_path": str(src_path),
            "symlink_path": str(dest_link_path),
            "link_target": str(link_target),
            "file_size_bytes": str(src_path.stat().st_size) if src_path.exists() else "0",
            "status": link_status,
        })

    return manifest, created_count, skipped_count, error_count


def write_manifest_csv(manifest: List[Dict[str, str]], output_csv: Path, logger: Optional[logging.Logger] = None):
    """Write the symlink operation manifest to a CSV file."""
    if not manifest:
        return

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["filename", "source_path", "symlink_path", "link_target", "file_size_bytes", "status"]

    with open(output_csv, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest)

    if logger:
        logger.info(f"Wrote asset manifest to: {output_csv.resolve()}")


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    workspace_root = Path(__file__).resolve().parents[2]
    default_raw_dir = workspace_root / "data" / "raw_vouchers"
    default_lm2_dir = workspace_root / "LM2_Project"

    parser = argparse.ArgumentParser(
        description="Pre-processor for LeafMachine2: Sets up LM2_Project/Data directory and creates safe symlinks.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input-dirs",
        "-i",
        nargs="+",
        type=Path,
        default=[default_raw_dir],
        help="One or more directories containing raw herbarium voucher sheet images.",
    )
    parser.add_argument(
        "--lm2-root",
        "-o",
        type=Path,
        default=default_lm2_dir,
        help="Root directory for the LeafMachine2 project (LM2_Project).",
    )
    parser.add_argument(
        "--images-subdir",
        type=str,
        default="images",
        help="Subdirectory name inside LM2_Project/Data where image symlinks will reside.",
    )
    parser.add_argument(
        "--output-subdir",
        type=str,
        default="output",
        help="Subdirectory name inside LM2_Project/Data for LM2 analysis output.",
    )
    parser.add_argument(
        "--relative",
        action="store_true",
        help="Create relative symlinks instead of absolute symlinks.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing symlinks or files at destination.",
    )
    parser.add_argument(
        "--verify-images",
        action="store_true",
        help="Perform PIL image header verification on each file before symlinking.",
    )
    parser.add_argument(
        "--limit",
        "-n",
        type=int,
        default=None,
        help="Limit number of images to symlink (useful for testing and dry runs).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Simulate the workflow without creating directories or symlinks.",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose debug logging.",
    )

    return parser.parse_args()


def main():
    args = parse_arguments()
    logger = setup_logging(verbose=args.verbose)

    logger.info("=================================================================")
    logger.info(" LeafMachine2 Dataset Pre-Processor & Symlink Manager")
    logger.info("=================================================================")
    logger.info(f"Target LM2 Project Root: {args.lm2_root.resolve()}")
    logger.info(f"Input Directories: {[str(p.resolve()) for p in args.input_dirs]}")
    logger.info(f"Dry Run Mode: {args.dry_run}")
    logger.info(f"Overwrite Mode: {args.overwrite}")
    logger.info(f"Symlink Mode: {'Relative' if args.relative else 'Absolute'}")

    # 1. Create directory structure
    images_dir, output_dir, configs_dir = create_lm2_directories(
        lm2_root=args.lm2_root,
        images_dir_name=args.images_subdir,
        output_dir_name=args.output_subdir,
        dry_run=args.dry_run,
        logger=logger,
    )
    logger.info(f"LM2 Image Input Directory:  {images_dir.resolve()}")
    logger.info(f"LM2 Output Directory:       {output_dir.resolve()}")
    logger.info(f"LM2 Configs Directory:      {configs_dir.resolve()}")

    # 2. Collect and filter images
    valid_images, filtered_out = collect_raw_images(
        input_dirs=args.input_dirs,
        verify_images=args.verify_images,
        logger=logger,
    )

    logger.info(f"Found {len(valid_images)} valid herbarium sheet image(s).")
    if filtered_out:
        logger.info(f"Filtered out {len(filtered_out)} non-image / invalid file(s).")
        for fpath, reason in filtered_out[:10]:
            logger.debug(f"Filtered: {fpath.name} -> Reason: {reason}")
        if len(filtered_out) > 10:
            logger.debug(f"... and {len(filtered_out) - 10} more filtered files.")

    if not valid_images:
        logger.warning("No valid images found to link. Exiting.")
        return

    # 3. Create symlinks
    manifest, created, skipped, errors = create_symlinks(
        source_images=valid_images,
        target_dir=images_dir,
        relative=args.relative,
        overwrite=args.overwrite,
        limit=args.limit,
        dry_run=args.dry_run,
        logger=logger,
    )

    # 4. Write manifest
    manifest_path = args.lm2_root / "Data" / "symlink_manifest.csv"
    if not args.dry_run:
        write_manifest_csv(manifest, manifest_path, logger=logger)

    logger.info("=================================================================")
    logger.info(" Pre-Processing Summary:")
    logger.info(f"   - Total Valid Images Discovered:  {len(valid_images)}")
    logger.info(f"   - Filtered Out Non-Images:        {len(filtered_out)}")
    logger.info(f"   - Symlinks Created/Updated:       {created}")
    logger.info(f"   - Symlinks Skipped (Already Exist):{skipped}")
    logger.info(f"   - Errors Encountered:             {errors}")
    if not args.dry_run:
        logger.info(f"   - Manifest CSV Generated:         {manifest_path.resolve()}")
    logger.info("=================================================================")
    logger.info("Done! LeafMachine2 input dataset is ready.")


if __name__ == "__main__":
    main()
