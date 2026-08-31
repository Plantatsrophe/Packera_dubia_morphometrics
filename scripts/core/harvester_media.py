"""
===============================================================================
Module: harvester_media.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Media extraction, URL optimization, image quality filtering, and
    asynchronous batch image downloading utilities for herbarium vouchers.
===============================================================================
"""

from __future__ import annotations

import asyncio
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import aiohttp
from PIL import Image
from tqdm.asyncio import tqdm as async_tqdm

Image.MAX_IMAGE_PIXELS = None  # Allow decompression of high-resolution botanical sheets (>89 MP)

from scripts.core.config import (
    DEFAULT_MIN_FILE_SIZE_KB,
    DEFAULT_MIN_MEGAPIXELS,
    DEFAULT_MIN_SHARPNESS_LAPLACIAN,
)


def optimize_herbarium_image_url(url: str) -> str:
    """
    Transforms provider-specific URLs to request full-resolution original scans rather than
    dimension-clamped web previews or downscaled thumbnails.

    Args:
        url: Raw image URL string.

    Returns:
        str: Optimized image URL requesting maximum resolution.
    """
    if not url or not isinstance(url, str):
        return ""

    optimized = url.strip()

    # 1. Smithsonian NMNH URLs: strip dimension clamp (e.g. &h=2000) to request native resolution
    if "collections.nmnh.si.edu/media/" in optimized:
        optimized = re.sub(r"[?&][hw]=\d+", "", optimized)
        if "?" not in optimized and "&" in optimized:
            optimized = optimized.replace("&", "?", 1)

    # 2. Symbiota / SERNEC / SEINet / CCH portals: replace web/thumbnail paths with original / large
    if any(k in optimized.lower() for k in ["symbiota", "sernec", "seinet", "cch2", "swbiodiversity"]):
        optimized = re.sub(r"/(?:web|tn|thumbnails?)/", "/orig/", optimized, flags=re.IGNORECASE)
        optimized = re.sub(r"_(?:tn|web|sm)\.(jpe?g|png)", r"_lg.\1", optimized, flags=re.IGNORECASE)

    # 3. IIIF Image API Endpoints: replace constrained dimensions with /full/max/
    if "/full/!" in optimized or "/full/pct:" in optimized or re.search(r"/full/\d+,\d*/", optimized):
        optimized = re.sub(r"/full/(?:!?\d+,\d*|pct:\d+)/", "/full/max/", optimized)

    return optimized


def extract_high_res_image_url(media_list: Optional[List[Dict[str, Any]]]) -> Optional[str]:
    """
    Parses Darwin Core media records, scoring and prioritizing the highest-quality specimen sheet image URL.

    Args:
        media_list: List of media dictionaries from the GBIF occurrence payload.

    Returns:
        Optional[str]: Verified and optimized HTTP/HTTPS image URL, or None if unavailable.
    """
    if not media_list or not isinstance(media_list, list):
        return None

    scored_candidates: List[Tuple[float, str]] = []

    for item in media_list:
        if not isinstance(item, dict):
            continue

        m_type = str(item.get("type", ""))
        m_format = str(item.get("format", "")).lower()
        identifier = str(item.get("identifier", "")).strip()

        if not identifier or not identifier.startswith(("http://", "https://")):
            continue

        score = 0.0

        if m_type == "StillImage" or "image" in m_format:
            score += 50.0
        if identifier.lower().endswith((".jpg", ".jpeg", ".png", ".tif", ".tiff")):
            score += 20.0

        ident_lower = identifier.lower()

        if any(h in ident_lower for h in ["_lg", "_large", "original", "/orig/", "/master/", "/highres/", "/full/", "hires", "high_res"]):
            score += 100.0
        if "max" in ident_lower:
            score += 30.0

        if any(t in ident_lower for t in ["_tn", "_thumb", "thumbnail", "_sm", "_small", "preview", "icon", "mini"]):
            score -= 150.0
        if "detailimages" in ident_lower:
            score -= 25.0
        if re.search(r"[?&]h=(?:[1-9]\d{0,2}|1\d{3}|2000)\b", ident_lower):
            score -= 10.0

        optimized_url = optimize_herbarium_image_url(identifier)
        scored_candidates.append((score, optimized_url))

    if not scored_candidates:
        return None

    scored_candidates.sort(key=lambda x: x[0], reverse=True)
    return scored_candidates[0][1]


def validate_image_quality(
    image_path: Path,
    min_megapixels: float = DEFAULT_MIN_MEGAPIXELS,
    min_file_size_kb: float = DEFAULT_MIN_FILE_SIZE_KB,
    check_sharpness: bool = False,
    min_sharpness: float = DEFAULT_MIN_SHARPNESS_LAPLACIAN,
) -> Tuple[bool, Dict[str, Any]]:
    """
    Evaluates image resolution, byte size, and optical sharpness metrics.

    Args:
        image_path: Path to the local image file.
        min_megapixels: Minimum resolution threshold in Megapixels (W x H / 1e6).
        min_file_size_kb: Minimum compressed file size in Kilobytes.
        check_sharpness: Whether to compute Laplacian variance edge sharpness.
        min_sharpness: Minimum acceptable Laplacian variance score.

    Returns:
        Tuple[bool, Dict[str, Any]]: (is_valid, metrics_dict)
    """
    if not image_path.exists():
        return False, {"valid": False, "reason": "file_not_found"}

    file_size_kb = image_path.stat().st_size / 1024.0
    if file_size_kb < min_file_size_kb:
        return False, {
            "valid": False,
            "reason": "file_size_too_small",
            "file_size_kb": round(file_size_kb, 1),
            "min_file_size_kb": min_file_size_kb,
        }

    try:
        with Image.open(image_path) as img:
            w, h = img.size
            mp = round((w * h) / 1e6, 2)
    except Exception as e:
        return False, {
            "valid": False,
            "reason": f"corrupt_or_unreadable_image: {e}",
            "file_size_kb": round(file_size_kb, 1),
        }

    if mp < min_megapixels:
        return False, {
            "valid": False,
            "reason": "low_resolution",
            "megapixels": mp,
            "min_megapixels": min_megapixels,
            "width": w,
            "height": h,
            "file_size_kb": round(file_size_kb, 1),
        }

    sharpness_score = None
    if check_sharpness:
        try:
            import cv2
            img_cv = cv2.imread(str(image_path), cv2.IMREAD_GRAYSCALE)
            if img_cv is not None:
                sharpness_score = round(float(cv2.Laplacian(img_cv, cv2.CV_64F).var()), 2)
                if sharpness_score < min_sharpness:
                    return False, {
                        "valid": False,
                        "reason": "blurry_or_upscaled",
                        "sharpness": sharpness_score,
                        "min_sharpness": min_sharpness,
                        "megapixels": mp,
                        "width": w,
                        "height": h,
                        "file_size_kb": round(file_size_kb, 1),
                    }
        except Exception:
            pass

    return True, {
        "valid": True,
        "megapixels": mp,
        "width": w,
        "height": h,
        "file_size_kb": round(file_size_kb, 1),
        "sharpness": sharpness_score,
    }


async def download_single_image(
    session: aiohttp.ClientSession,
    image_url: str,
    destination_path: Path,
    semaphore: asyncio.Semaphore,
    max_retries: int = 3
) -> bool:
    """
    Asynchronously downloads a single voucher image file to local storage with retry logic.

    Args:
        session: Active aiohttp ClientSession.
        image_url: Remote image URL string.
        destination_path: Local filesystem destination Path.
        semaphore: Asyncio Semaphore for concurrency throttling.
        max_retries: Number of exponential backoff retry attempts (default: 3).

    Returns:
        bool: True if image downloaded successfully or already existed on disk, False otherwise.
    """
    if destination_path.exists() and destination_path.stat().st_size > 1024:
        return True

    async with semaphore:
        for attempt in range(1, max_retries + 1):
            try:
                timeout = aiohttp.ClientTimeout(total=45, connect=15)
                async with session.get(image_url, timeout=timeout) as response:
                    if response.status == 200:
                        content_type = response.headers.get("Content-Type", "").lower()
                        content = await response.read()
                        if len(content) > 1024 and (not content_type or "image" in content_type or "octet-stream" in content_type or content[:3] == b'\xff\xd8\xff'):
                            destination_path.parent.mkdir(parents=True, exist_ok=True)
                            temp_path = destination_path.with_suffix(".tmp")
                            with open(temp_path, "wb") as f:
                                f.write(content)
                            temp_path.replace(destination_path)
                            return True
                    elif response.status in {404, 410}:
                        return False
            except (aiohttp.ClientError, asyncio.TimeoutError, Exception):
                if attempt == max_retries:
                    return False
                await asyncio.sleep(1.0 * (2 ** (attempt - 1)))
        return False


async def download_all_voucher_images(
    records_to_download: List[Tuple[str, Path]],
    concurrency_limit: int = 15,
    logger: Optional[logging.Logger] = None
) -> Dict[str, int]:
    """
    Coordinates asynchronous batch downloading of voucher images with concurrency control and progress tracking.

    Args:
        records_to_download: List of tuples where each tuple contains (image_url, destination_path).
        concurrency_limit: Maximum concurrent HTTP connections (default: 15).
        logger: Optional Logger instance for diagnostics.

    Returns:
        Dict[str, int]: Download metrics dictionary with 'success', 'skipped', and 'failed' counts.
    """
    semaphore = asyncio.Semaphore(concurrency_limit)
    headers = {
        "User-Agent": "PackeraResearchBot/1.0 (UNC Chapel Hill Herbarium; Evolutionary Morphometrics Lab)"
    }

    stats = {"success": 0, "skipped": 0, "failed": 0}

    pending = []
    for url, dest in records_to_download:
        if dest.exists() and dest.stat().st_size > 1024:
            stats["skipped"] += 1
        else:
            pending.append((url, dest))

    if not pending:
        if logger:
            logger.info(f"All {stats['skipped']} voucher images are already cached locally.")
        return stats

    if logger:
        logger.info(f"Initiating asynchronous download of {len(pending)} pending images (Concurrency: {concurrency_limit})...")

    connector = aiohttp.TCPConnector(limit=concurrency_limit, limit_per_host=5, ssl=False)
    async with aiohttp.ClientSession(headers=headers, connector=connector) as session:
        tasks = [
            download_single_image(session, url, dest, semaphore)
            for url, dest in pending
        ]

        results = await async_tqdm.gather(*tasks, desc="Downloading Voucher Sheets", unit="img")

        for success in results:
            if success:
                stats["success"] += 1
            else:
                stats["failed"] += 1

    return stats
