#!/usr/bin/env python3
"""
===============================================================================
Script: download_environmental_layers.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Idempotent, atomic preflight provisioner for fine-grain environmental layers
    used in Step 05 ecological niche modeling and substrate contingency audits:
      1. USDA-NRCS SOLUS100 (100 m CONUS soil & lithic rasters, Nauman et al. 2024):
         - anylithicdpt_cm_p.tif (depth to lithic/paralithic bedrock, cm; scalar=1)
         - resdept_all_cm_p.tif  (depth to root restriction, cm; scalar=1)
         - fragvol_0_cm_p.tif, fragvol_5_cm_p.tif, fragvol_15_cm_p.tif (% vol; scalar=1)
         - ph1to1h2o_0_cm_p.tif, ph1to1h2o_5_cm_p.tif, ph1to1h2o_15_cm_p.tif (pH*100; scalar=100)
      2. USGS State Geologic Map Compilation (SGMC, Horton et al. 2017, DS-1052):
         - Resolves official bundle from USGS ScienceBase item 5888bf4fe4b05ccb964bab9d
         - Pre-joins SGMC_Geology polygons with SGMC_Units lithology attributes by UNIT_LINK
         - Atomically writes data/environmental/sgmc/sgmc_conus_lithology.gpkg
      3. 1-arc-second (~30 m) Digital Elevation Model tiles (AWS Terrain Skadi / USGS 3DEP):
         - Determines 1x1-degree tiles covering study voucher coordinates
         - Fetches missing 30 m tiles from AWS Open Data Terrain Tiles (skadi .hgt.gz)
           with automatic fallback to USGS 3DEP National Map ImageServer
         - Converts to DEFLATE-compressed GeoTIFFs in data/environmental/dem_3dep_30m/

    Note: 30 m POLARIS (polaris_30m/) is intentionally NOT downloaded by default
    because SSURGO SDA + 100 m SOLUS100 + USGS SGMC + 30 m DEM already provide
    superior lithic, outcrop, and soil chemistry resolution without fetching ~1,800 tiles.
===============================================================================
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import json
import logging
import math
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

logger = logging.getLogger("packera.env_layers")
if not logger.handlers:
    _handler = logging.StreamHandler(sys.stdout)
    _handler.setFormatter(logging.Formatter("[Env Layers] %(message)s"))
    logger.addHandler(_handler)
    logger.setLevel(logging.INFO)

USER_AGENT = "PackeraMorphometricsPipeline/1.0 (UNC-Chapel-Hill-Herbarium-NCU)"

# Official USDA-NRCS SOLUS100 public bucket (Nauman et al. 2024; doi:10.15482/USDA.ADC/25033856)
SOLUS100_BASE_URL = "https://storage.googleapis.com/solus100pub"
SOLUS100_REQUIRED_LAYERS: Tuple[Dict[str, Any], ...] = (
    {
        "filename": "anylithicdpt_cm_p.tif",
        "property": "anylithicdpt",
        "scalar": 1,
        "units": "cm",
        "description": "Depth to lithic or paralithic bedrock contact",
    },
    {
        "filename": "resdept_all_cm_p.tif",
        "property": "resdept",
        "scalar": 1,
        "units": "cm",
        "description": "Depth to any soil root restriction",
    },
    {
        "filename": "fragvol_0_cm_p.tif",
        "property": "fragvol",
        "scalar": 1,
        "units": "percent_volume",
        "description": "Coarse rock fragments (>2 mm) at 0 cm",
    },
    {
        "filename": "fragvol_5_cm_p.tif",
        "property": "fragvol",
        "scalar": 1,
        "units": "percent_volume",
        "description": "Coarse rock fragments (>2 mm) at 5 cm",
    },
    {
        "filename": "fragvol_15_cm_p.tif",
        "property": "fragvol",
        "scalar": 1,
        "units": "percent_volume",
        "description": "Coarse rock fragments (>2 mm) at 15 cm",
    },
    {
        "filename": "ph1to1h2o_0_cm_p.tif",
        "property": "ph1to1h2o",
        "scalar": 100,
        "units": "pH_x100",
        "description": "Soil pH (1:1 H2O) at 0 cm (divide integer raster by 100)",
    },
    {
        "filename": "ph1to1h2o_5_cm_p.tif",
        "property": "ph1to1h2o",
        "scalar": 100,
        "units": "pH_x100",
        "description": "Soil pH (1:1 H2O) at 5 cm (divide integer raster by 100)",
    },
    {
        "filename": "ph1to1h2o_15_cm_p.tif",
        "property": "ph1to1h2o",
        "scalar": 100,
        "units": "pH_x100",
        "description": "Soil pH (1:1 H2O) at 15 cm (divide integer raster by 100)",
    },
    {
        "filename": "sandtotal_0_cm_p.tif",
        "property": "sandtotal",
        "scalar": 1,
        "units": "percent_weight",
        "description": "Total sand (0.05-2.0 mm, %) at 0 cm",
    },
    {
        "filename": "sandtotal_5_cm_p.tif",
        "property": "sandtotal",
        "scalar": 1,
        "units": "percent_weight",
        "description": "Total sand (0.05-2.0 mm, %) at 5 cm",
    },
    {
        "filename": "sandtotal_15_cm_p.tif",
        "property": "sandtotal",
        "scalar": 1,
        "units": "percent_weight",
        "description": "Total sand (0.05-2.0 mm, %) at 15 cm",
    },
    {
        "filename": "claytotal_0_cm_p.tif",
        "property": "claytotal",
        "scalar": 1,
        "units": "percent_weight",
        "description": "Total clay (<0.002 mm, %) at 0 cm",
    },
    {
        "filename": "claytotal_5_cm_p.tif",
        "property": "claytotal",
        "scalar": 1,
        "units": "percent_weight",
        "description": "Total clay (<0.002 mm, %) at 5 cm",
    },
    {
        "filename": "claytotal_15_cm_p.tif",
        "property": "claytotal",
        "scalar": 1,
        "units": "percent_weight",
        "description": "Total clay (<0.002 mm, %) at 15 cm",
    },
    {
        "filename": "silttotal_0_cm_p.tif",
        "property": "silttotal",
        "scalar": 1,
        "units": "percent_weight",
        "description": "Total silt (0.002-0.05 mm, %) at 0 cm",
    },
    {
        "filename": "silttotal_5_cm_p.tif",
        "property": "silttotal",
        "scalar": 1,
        "units": "percent_weight",
        "description": "Total silt (0.002-0.05 mm, %) at 5 cm",
    },
    {
        "filename": "silttotal_15_cm_p.tif",
        "property": "silttotal",
        "scalar": 1,
        "units": "percent_weight",
        "description": "Total silt (0.002-0.05 mm, %) at 15 cm",
    },
    {
        "filename": "cec7_0_cm_p.tif",
        "property": "cec7",
        "scalar": 10,
        "units": "cmol_kg_x10",
        "description": "Cation exchange capacity (pH 7, cmol+/kg x 10) at 0 cm",
    },
    {
        "filename": "cec7_5_cm_p.tif",
        "property": "cec7",
        "scalar": 10,
        "units": "cmol_kg_x10",
        "description": "Cation exchange capacity (pH 7, cmol+/kg x 10) at 5 cm",
    },
    {
        "filename": "cec7_15_cm_p.tif",
        "property": "cec7",
        "scalar": 10,
        "units": "cmol_kg_x10",
        "description": "Cation exchange capacity (pH 7, cmol+/kg x 10) at 15 cm",
    },
)

# Official USGS SGMC ScienceBase catalog item (Horton et al. 2017; doi:10.5066/F7WH2N65)
SGMC_SCIENCEBASE_ITEM_ID = "5888bf4fe4b05ccb964bab9d"
SGMC_SCIENCEBASE_JSON_URL = (
    f"https://www.sciencebase.gov/catalog/item/{SGMC_SCIENCEBASE_ITEM_ID}?format=json"
)
SGMC_SCIENCEBASE_CHILDREN_URL = (
    f"https://www.sciencebase.gov/catalog/items?parentId={SGMC_SCIENCEBASE_ITEM_ID}&format=json"
)
SGMC_TARGET_GPKG = "sgmc_conus_lithology.gpkg"

# AWS Open Data Terrain Tiles (Skadi 1x1 deg 1-arc-second SRTM/3DEP HGT) & USGS 3DEP REST fallback
SKADI_BASE_URL = "https://s3.amazonaws.com/elevation-tiles-prod/skadi"
USGS_3DEP_EXPORT_URL = (
    "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage"
)

MIN_VALID_FILE_BYTES = 1024


def format_dem_tile_id(lat_floor: int, lon_floor: int) -> str:
    """Format integer SW-corner (lat_floor, lon_floor) as SRTM/Skadi tile ID (e.g. N35W079)."""
    ns = "N" if lat_floor >= 0 else "S"
    ew = "E" if lon_floor >= 0 else "W"
    return f"{ns}{abs(int(lat_floor)):02d}{ew}{abs(int(lon_floor)):03d}"


def _is_nonempty_file(path: Path, min_bytes: int = MIN_VALID_FILE_BYTES) -> bool:
    try:
        return path.is_file() and path.stat().st_size >= min_bytes
    except OSError:
        return False


def _is_valid_tiff(path: Path) -> bool:
    if not _is_nonempty_file(path, min_bytes=256):
        return False
    try:
        with path.open("rb") as fh:
            hdr = fh.read(4)
        return hdr in (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")
    except OSError:
        return False


def _atomic_download(
    url: str,
    dest_path: Path,
    timeout: float = 120.0,
    retries: int = 3,
    min_bytes: int = MIN_VALID_FILE_BYTES,
) -> bool:
    """Download URL atomically to dest_path with exponential backoff."""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = dest_path.with_name(f".{dest_path.name}.part.{os.getpid()}")
    last_err: Optional[Exception] = None

    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                status = getattr(resp, "status", 200)
                if status and int(status) >= 400:
                    raise urllib.error.HTTPError(url, int(status), "HTTP error", resp.headers, None)
                with tmp_path.open("wb") as out_f:
                    shutil.copyfileobj(resp, out_f, length=1024 * 256)
                    out_f.flush()
                    os.fsync(out_f.fileno())
            if not _is_nonempty_file(tmp_path, min_bytes=min_bytes):
                raise IOError(f"Downloaded file from {url} was smaller than {min_bytes} bytes")
            tmp_path.replace(dest_path)
            return True
        except urllib.error.HTTPError as exc:
            last_err = exc
            if tmp_path.exists():
                tmp_path.unlink(missing_ok=True)
            if exc.code == 404:
                return False
        except Exception as exc:
            last_err = exc
            if tmp_path.exists():
                tmp_path.unlink(missing_ok=True)

        if attempt < retries:
            time.sleep(min(6.0, 1.2 * (2 ** (attempt - 1))))

    if last_err is not None:
        logger.warning(f"Download failed for {url}: {last_err}")
    return False


def _fetch_json(url: str, timeout: float = 45.0) -> Optional[Dict[str, Any]]:
    try:
        req = urllib.request.Request(
            url,
            headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read().decode("utf-8", errors="replace")
        data = json.loads(raw)
        return data if isinstance(data, dict) else None
    except Exception as exc:
        logger.warning(f"Metadata JSON query failed for {url}: {exc}")
        return None


# =============================================================================
# 1. SOLUS100 (100 m Soil & Bedrock Depth Rasters)
# =============================================================================

def inspect_solus_cache(solus_dir: Path) -> Tuple[bool, List[Dict[str, Any]]]:
    """Return (all_ready, missing_specs) for SOLUS100 rasters in solus_dir."""
    if not solus_dir.is_dir():
        return False, list(SOLUS100_REQUIRED_LAYERS)

    missing: List[Dict[str, Any]] = []
    for spec in SOLUS100_REQUIRED_LAYERS:
        target = solus_dir / str(spec["filename"])
        if not _is_valid_tiff(target):
            missing.append(spec)

    if not missing:
        return True, []

    # Also honor pre-existing custom-named SOLUS rasters only when bedrock, fragvol, pH, texture, and CEC all exist
    tifs = [p for p in solus_dir.glob("*.tif*") if _is_valid_tiff(p)]
    has_bed = any(re.search(r"(anylithicdpt|resdept)", p.name, re.I) for p in tifs)
    has_frg = any(re.search(r"(fragvol|cfvo)", p.name, re.I) for p in tifs)
    has_ph = any(re.search(r"^ph(_|1to1h2o)", p.name, re.I) for p in tifs)
    has_tex = (
        any(re.search(r"sandtotal", p.name, re.I) for p in tifs)
        and any(re.search(r"claytotal", p.name, re.I) for p in tifs)
        and any(re.search(r"silttotal", p.name, re.I) for p in tifs)
    )
    has_cec = any(re.search(r"cec7", p.name, re.I) for p in tifs)
    if has_bed and has_frg and has_ph and has_tex and has_cec:
        return True, []

    return False, missing


def ensure_solus100_layers(
    env_dir: Path,
    force: bool = False,
    max_workers: int = 4,
    timeout: float = 180.0,
) -> Dict[str, Any]:
    """Ensure required USDA-NRCS SOLUS100 100 m GeoTIFFs exist in env_dir/solus_100m."""
    solus_dir = env_dir / "solus_100m"
    solus_dir.mkdir(parents=True, exist_ok=True)

    if not force:
        ready, missing_specs = inspect_solus_cache(solus_dir)
        if ready:
            return {"status": "cached", "downloaded": 0, "dir": str(solus_dir)}
    else:
        missing_specs = list(SOLUS100_REQUIRED_LAYERS)

    logger.info(
        f"Downloading {len(missing_specs)} missing USDA-NRCS SOLUS100 100 m raster(s) into {solus_dir}..."
    )

    def _dl_one(spec: Dict[str, Any]) -> Tuple[str, bool]:
        fname = str(spec["filename"])
        url = f"{SOLUS100_BASE_URL}/{fname}"
        dest = solus_dir / fname
        ok = _atomic_download(url, dest, timeout=timeout, retries=3, min_bytes=MIN_VALID_FILE_BYTES)
        if ok and not _is_valid_tiff(dest):
            dest.unlink(missing_ok=True)
            ok = False
        return fname, ok

    downloaded = 0
    failed: List[str] = []
    workers = max(1, min(int(max_workers), len(missing_specs)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for fname, ok in pool.map(_dl_one, missing_specs):
            if ok:
                downloaded += 1
            else:
                failed.append(fname)

    # Persist scalar metadata manifest alongside rasters for auditability
    manifest_path = solus_dir / "solus100_layers_manifest.json"
    manifest_payload = {
        "source": "USDA-NRCS Soil Landscapes of the United States 100-meter (SOLUS100)",
        "doi": "https://doi.org/10.15482/USDA.ADC/25033856",
        "index_url": f"{SOLUS100_BASE_URL}/index.html",
        "layers": list(SOLUS100_REQUIRED_LAYERS),
    }
    tmp_manifest = manifest_path.with_suffix(".json.tmp")
    tmp_manifest.write_text(json.dumps(manifest_payload, indent=2) + "\n", encoding="utf-8")
    tmp_manifest.replace(manifest_path)

    status = "downloaded" if not failed else ("partial" if downloaded > 0 else "failed")
    return {
        "status": status,
        "downloaded": downloaded,
        "failed": failed,
        "dir": str(solus_dir),
    }


# =============================================================================
# 2. USGS SGMC (State Geologic Map Compilation Pre-Joined Lithology GeoPackage)
# =============================================================================

def inspect_sgmc_cache(sgmc_dir: Path, env_dir: Optional[Path] = None) -> Optional[Path]:
    """Return path to an existing valid SGMC vector file if present, else None."""
    search_dirs = [sgmc_dir]
    if env_dir is not None and env_dir != sgmc_dir:
        search_dirs.append(env_dir)

    for d in search_dirs:
        if not d.is_dir():
            continue
        for p in sorted(d.iterdir()):
            if re.match(r"^sgmc.*\.(gpkg|shp|geojson)$", p.name, flags=re.I) and _is_nonempty_file(p):
                return p
    return None


def resolve_sgmc_download_candidates(item_json: Dict[str, Any]) -> List[Tuple[str, str]]:
    """Extract candidate (filename, download_url) pairs for SGMC zip bundles from ScienceBase JSON."""
    candidates: List[Tuple[str, str]] = []

    def _collect_from_item(obj: Dict[str, Any]) -> None:
        for f_entry in obj.get("files") or []:
            if not isinstance(f_entry, dict):
                continue
            name = str(f_entry.get("name") or "").strip()
            uri = str(f_entry.get("downloadUri") or f_entry.get("url") or "").strip()
            if name.lower().endswith(".zip") and uri.startswith("http"):
                candidates.append((name, uri))
        for facet in obj.get("facets") or []:
            if not isinstance(facet, dict):
                continue
            for f_entry in facet.get("files") or []:
                if not isinstance(f_entry, dict):
                    continue
                name = str(f_entry.get("name") or "").strip()
                uri = str(f_entry.get("downloadUri") or f_entry.get("url") or "").strip()
                if name.lower().endswith(".zip") and uri.startswith("http"):
                    candidates.append((name, uri))
        for wl in obj.get("webLinks") or []:
            if not isinstance(wl, dict):
                continue
            uri = str(wl.get("uri") or "").strip()
            title = str(wl.get("title") or os.path.basename(urllib.parse.urlparse(uri).path)).strip()
            if uri.startswith("http") and (uri.lower().endswith(".zip") or "shapefile" in title.lower() or "gdb" in title.lower()):
                candidates.append((title or "sgmc_bundle.zip", uri))

    _collect_from_item(item_json)

    # Deduplicate while prioritizing Shapefiles zip first (fastest for sf::st_read), then GDB zip
    seen: Set[str] = set()
    uniq: List[Tuple[str, str]] = []
    for name, uri in candidates:
        if uri not in seen:
            seen.add(uri)
            uniq.append((name, uri))

    def _rank_key(item: Tuple[str, str]) -> Tuple[int, str]:
        nm_low = item[0].lower()
        if "shapefile" in nm_low or "shp" in nm_low:
            return (0, nm_low)
        if "geodatabase" in nm_low or "gdb" in nm_low or "sgmc" in nm_low:
            return (1, nm_low)
        return (2, nm_low)

    uniq.sort(key=_rank_key)
    return uniq


def _build_sgmc_join_r_script() -> str:
    """Return R script source that joins SGMC_Geology with SGMC_Units on UNIT_LINK and writes GPKG."""
    return """
args <- commandArgs(trailingOnly = TRUE)
staging_dir <- args[1]
out_gpkg    <- args[2]

suppressPackageStartupMessages({
  library(sf)
})

all_files <- list.files(staging_dir, recursive = TRUE, full.names = TRUE)
gdb_dirs  <- unique(dirname(all_files[grepl("\\\\.gdb(/|\\\\\\\\)", all_files, ignore.case = TRUE)]))
if (length(gdb_dirs) == 0L) {
  cand_dirs <- list.dirs(staging_dir, recursive = TRUE, full.names = TRUE)
  gdb_dirs  <- cand_dirs[grepl("\\\\.gdb$", cand_dirs, ignore.case = TRUE)]
}

shp_geol <- all_files[grepl("SGMC_Geology\\\\.shp$", all_files, ignore.case = TRUE)]
geol_sf  <- NULL
units_df <- NULL

if (length(shp_geol) > 0L) {
  geol_sf <- sf::st_read(shp_geol[1L], quiet = TRUE, stringsAsFactors = FALSE)
  units_cand <- all_files[grepl("SGMC_Units\\\\.(dbf|csv)$|(^|/|\\\\\\\\)Units\\\\.(dbf|csv)$", all_files, ignore.case = TRUE)]
  if (length(units_cand) > 0L) {
    u_path <- units_cand[1L]
    if (grepl("\\\\.csv$", u_path, ignore.case = TRUE)) {
      units_df <- utils::read.csv(u_path, stringsAsFactors = FALSE)
    } else if (requireNamespace("foreign", quietly = TRUE)) {
      units_df <- foreign::read.dbf(u_path, as.is = TRUE)
    } else {
      units_df <- as.data.frame(sf::st_read(u_path, quiet = TRUE, stringsAsFactors = FALSE))
    }
  }
} else if (length(gdb_dirs) > 0L) {
  gdb_path <- gdb_dirs[1L]
  lyrs <- sf::st_layers(gdb_path)$name
  g_lyr <- lyrs[grepl("SGMC_Geology|Geology", lyrs, ignore.case = TRUE)][1L]
  u_lyr <- lyrs[grepl("^SGMC_Units$|^Units$", lyrs, ignore.case = TRUE)][1L]
  if (!is.na(g_lyr)) {
    geol_sf <- sf::st_read(gdb_path, layer = g_lyr, quiet = TRUE, stringsAsFactors = FALSE)
  }
  if (!is.na(u_lyr)) {
    units_df <- as.data.frame(sf::st_read(gdb_path, layer = u_lyr, quiet = TRUE, stringsAsFactors = FALSE))
  }
}

if (is.null(geol_sf) || nrow(geol_sf) == 0L) {
  stop("Failed to locate non-empty SGMC_Geology polygon layer in downloaded bundle.")
}

pick_col <- function(df_obj, candidates) {
  if (is.null(df_obj)) return(rep("", nrow(geol_sf)))
  nms_up <- toupper(names(df_obj))
  for (cand in toupper(candidates)) {
    idx <- match(cand, nms_up)
    if (!is.na(idx)) {
      v <- as.character(df_obj[[idx]])
      v[is.na(v)] <- ""
      return(trimws(v))
    }
  }
  rep("", nrow(df_obj))
}

g_link <- pick_col(geol_sf, c("UNIT_LINK", "UNITLINK"))
g_gen  <- pick_col(geol_sf, c("GENERALIZED_LITH", "GENERALIZE", "GEN_LITH"))
g_orig <- pick_col(geol_sf, c("ORIG_LABEL", "SGMC_LABEL"))

maj1 <- pick_col(geol_sf, c("MAJOR1", "ROCKTYPE1"))
maj2 <- pick_col(geol_sf, c("MAJOR2", "ROCKTYPE2"))
min1 <- pick_col(geol_sf, c("MINOR1"))
unm  <- pick_col(geol_sf, c("UNIT_NAME", "UNITNAME"))

if (!is.null(units_df) && nrow(units_df) > 0L) {
  u_link <- pick_col(units_df, c("UNIT_LINK", "UNITLINK"))
  keep_u <- nzchar(u_link) & !duplicated(u_link)
  if (any(keep_u)) {
    units_u <- units_df[keep_u, , drop = FALSE]
    u_link_u <- u_link[keep_u]
    m_idx <- match(g_link, u_link_u)
    hit <- !is.na(m_idx)
    if (any(hit)) {
      u_maj1 <- pick_col(units_u, c("MAJOR1", "ROCKTYPE1"))[m_idx[hit]]
      u_maj2 <- pick_col(units_u, c("MAJOR2", "ROCKTYPE2"))[m_idx[hit]]
      u_min1 <- pick_col(units_u, c("MINOR1"))[m_idx[hit]]
      u_name <- pick_col(units_u, c("UNIT_NAME", "UNITNAME", "STRAT_UNIT"))[m_idx[hit]]
      maj1[hit] <- ifelse(nzchar(u_maj1), u_maj1, maj1[hit])
      maj2[hit] <- ifelse(nzchar(u_maj2), u_maj2, maj2[hit])
      min1[hit] <- ifelse(nzchar(u_min1), u_min1, min1[hit])
      unm[hit]  <- ifelse(nzchar(u_name), u_name, unm[hit])
    }
  }
}

rk2 <- ifelse(nzchar(maj2), maj2, min1)
geom_col <- sf::st_geometry(geol_sf)
out_sf <- sf::st_sf(
  UNIT_LINK  = g_link,
  ORIG_LABEL = g_orig,
  UNIT_NAME  = unm,
  MAJOR1     = maj1,
  MAJOR2     = maj2,
  ROCKTYPE1  = maj1,
  ROCKTYPE2  = rk2,
  GENERALIZE = g_gen,
  geometry   = geom_col
)

tmp_gpkg <- paste0(out_gpkg, ".tmp.gpkg")
if (file.exists(tmp_gpkg)) unlink(tmp_gpkg)
sf::st_write(out_sf, dsn = tmp_gpkg, layer = "sgmc_lithology", driver = "GPKG", quiet = TRUE, append = FALSE)
if (!file.exists(tmp_gpkg) || file.info(tmp_gpkg)$size < 1024) {
  stop("Failed to write valid sgmc_conus_lithology.gpkg")
}
if (file.exists(out_gpkg)) unlink(out_gpkg)
if (!file.rename(tmp_gpkg, out_gpkg)) {
  file.copy(tmp_gpkg, out_gpkg, overwrite = TRUE)
  unlink(tmp_gpkg)
}
cat(sprintf("SGMC_GPKG_OK rows=%d\\n", nrow(out_sf)))
"""


def ensure_sgmc_layer(
    env_dir: Path,
    force: bool = False,
    rscript_bin: Optional[str] = None,
    timeout: float = 300.0,
) -> Dict[str, Any]:
    """Ensure USGS SGMC lithology polygons (pre-joined with SGMC_Units) exist in env_dir/sgmc."""
    sgmc_dir = env_dir / "sgmc"
    sgmc_dir.mkdir(parents=True, exist_ok=True)

    if not force:
        existing = inspect_sgmc_cache(sgmc_dir, env_dir=env_dir)
        if existing is not None:
            return {"status": "cached", "path": str(existing)}

    r_bin = rscript_bin or shutil.which("Rscript")
    if not r_bin:
        logger.warning("Rscript not found in PATH; skipping USGS SGMC GeoPackage build.")
        return {"status": "skipped_no_rscript", "dir": str(sgmc_dir)}

    item_json = _fetch_json(SGMC_SCIENCEBASE_JSON_URL, timeout=45.0)
    candidates: List[Tuple[str, str]] = []
    if item_json is not None:
        candidates.extend(resolve_sgmc_download_candidates(item_json))
    if not candidates:
        children_json = _fetch_json(SGMC_SCIENCEBASE_CHILDREN_URL, timeout=45.0)
        if children_json is not None:
            for child in children_json.get("items") or []:
                if isinstance(child, dict):
                    candidates.extend(resolve_sgmc_download_candidates(child))

    if not candidates:
        logger.warning(
            f"Could not resolve SGMC bundle download URL from ScienceBase ({SGMC_SCIENCEBASE_JSON_URL})."
        )
        return {"status": "failed_metadata", "dir": str(sgmc_dir)}

    target_gpkg = sgmc_dir / SGMC_TARGET_GPKG
    with tempfile.TemporaryDirectory(prefix="packera_sgmc_") as tmp_dir:
        tmp_root = Path(tmp_dir)
        zip_path = tmp_root / "sgmc_bundle.zip"
        dl_ok = False
        chosen_name = ""
        for cand_name, cand_url in candidates:
            logger.info(f"Downloading USGS SGMC bundle ({cand_name})...")
            if _atomic_download(cand_url, zip_path, timeout=timeout, retries=3, min_bytes=1024 * 64):
                dl_ok = True
                chosen_name = cand_name
                break

        if not dl_ok:
            return {"status": "failed_download", "dir": str(sgmc_dir)}

        extract_dir = tmp_root / "extracted"
        extract_dir.mkdir(parents=True, exist_ok=True)
        try:
            with zipfile.ZipFile(zip_path, "r") as zf:
                zf.extractall(extract_dir)
        except Exception as exc:
            logger.warning(f"Failed to extract SGMC zip bundle ({chosen_name}): {exc}")
            return {"status": "failed_unzip", "dir": str(sgmc_dir)}

        r_script_file = tmp_root / "build_sgmc_gpkg.R"
        r_script_file.write_text(_build_sgmc_join_r_script(), encoding="utf-8")
        proc = subprocess.run(
            [r_bin, str(r_script_file), str(extract_dir), str(target_gpkg)],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0 or not _is_nonempty_file(target_gpkg):
            logger.warning(
                f"USGS SGMC relational join in R failed (exit {proc.returncode}): {proc.stderr.strip()}"
            )
            return {"status": "failed_r_join", "dir": str(sgmc_dir)}

    logger.info(f"Built pre-joined USGS SGMC lithology layer: {target_gpkg}")
    return {"status": "downloaded", "path": str(target_gpkg)}


# =============================================================================
# 3. USGS 3DEP / AWS Terrain 1-Arc-Second (~30 m) Elevation Tiles
# =============================================================================

def collect_voucher_dem_tiles(
    vouchers_csv: Optional[Path] = None,
    workspace_root: Path = PROJECT_ROOT,
    db_path: Optional[Path] = None,
) -> Set[Tuple[int, int]]:
    """Collect unique 1x1-degree SW-corner (floor(lat), floor(lon)) cells across study vouchers."""
    tiles: Set[Tuple[int, int]] = set()

    def _add_coord(lat_raw: Any, lon_raw: Any) -> None:
        try:
            lat = float(lat_raw)
            lon = float(lon_raw)
        except (TypeError, ValueError):
            return
        if not (math.isfinite(lat) and math.isfinite(lon)):
            return
        # Restrict to North American terrestrial bounds covered by 3DEP / CDED / SRTM1
        if 14.0 <= lat <= 60.0 and -168.0 <= lon <= -50.0:
            tiles.add((int(math.floor(lat)), int(math.floor(lon))))

    csv_candidates: List[Path] = []
    if vouchers_csv is not None:
        csv_candidates.append(Path(vouchers_csv))
    canonical_csv = workspace_root / "data" / "tables" / "curated_vouchers.csv"
    if canonical_csv not in csv_candidates:
        csv_candidates.append(canonical_csv)

    for c_path in csv_candidates:
        if not c_path.is_file():
            continue
        try:
            with c_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
                reader = csv.DictReader(fh)
                if not reader.fieldnames:
                    continue
                lat_col = next(
                    (c for c in ("decimalLatitude", "lat") if c in reader.fieldnames),
                    None,
                )
                lon_col = next(
                    (c for c in ("decimalLongitude", "lon") if c in reader.fieldnames),
                    None,
                )
                if not lat_col or not lon_col:
                    continue
                for row in reader:
                    _add_coord(row.get(lat_col), row.get(lon_col))
            if vouchers_csv is not None and c_path == Path(vouchers_csv) and tiles:
                # Active run cohort already resolved all needed 1x1 deg tiles
                return tiles
        except Exception as exc:
            logger.debug(f"Could not read coordinates from {c_path}: {exc}")

    if not tiles:
        sql_path = db_path or (workspace_root / "data" / "packera_pipeline.sqlite")
        if sql_path.is_file():
            try:
                conn = sqlite3.connect(f"file:{sql_path}?mode=ro", uri=True, timeout=5.0)
                try:
                    cur = conn.execute(
                        "SELECT decimalLatitude, decimalLongitude FROM vouchers "
                        "WHERE decimalLatitude IS NOT NULL AND decimalLongitude IS NOT NULL"
                    )
                    for lat_v, lon_v in cur:
                        _add_coord(lat_v, lon_v)
                finally:
                    conn.close()
            except Exception as exc:
                logger.debug(f"SQLite coordinate lookup skipped: {exc}")

    return tiles


def _load_void_tiles_memo(dem_dir: Path) -> Set[str]:
    memo_path = dem_dir / "ocean_void_tiles.json"
    if not memo_path.is_file():
        return set()
    try:
        data = json.loads(memo_path.read_text(encoding="utf-8"))
        if isinstance(data, list):
            return {str(x) for x in data}
    except Exception:
        pass
    return set()


def _save_void_tiles_memo(dem_dir: Path, void_tiles: Set[str]) -> None:
    memo_path = dem_dir / "ocean_void_tiles.json"
    tmp_p = memo_path.with_suffix(".json.tmp")
    tmp_p.write_text(json.dumps(sorted(void_tiles), indent=2) + "\n", encoding="utf-8")
    tmp_p.replace(memo_path)


def _download_single_3dep_tile(
    lat_floor: int,
    lon_floor: int,
    staging_dir: Path,
    dem_dir: Path,
    timeout: float = 90.0,
) -> Tuple[str, str, Optional[Path]]:
    """
    Fetch a single 1x1-degree ~30 m DEM tile.
    Primary: AWS Open Data Skadi SRTM1/3DEP .hgt.gz (1 arc-second, 3601x3601).
    Fallback: USGS 3DEP National Map ImageServer exportImage GeoTIFF.
    Returns (tile_id, outcome, staged_path) where outcome is 'hgt', 'tif', 'void', or 'failed'.
    """
    tile_id = format_dem_tile_id(lat_floor, lon_floor)
    ns_folder = f"{'N' if lat_floor >= 0 else 'S'}{abs(int(lat_floor)):02d}"
    skadi_url = f"{SKADI_BASE_URL}/{ns_folder}/{tile_id}.hgt.gz"
    gz_path = staging_dir / f"{tile_id}.hgt.gz"
    hgt_path = staging_dir / f"{tile_id}.hgt"

    if _atomic_download(skadi_url, gz_path, timeout=timeout, retries=2, min_bytes=512):
        try:
            with gzip.open(gz_path, "rb") as gz_in, hgt_path.open("wb") as hgt_out:
                shutil.copyfileobj(gz_in, hgt_out, length=1024 * 256)
            gz_path.unlink(missing_ok=True)
            # Valid SRTM1 (3601x3601*2 = 25,934,402) or SRTM3 (1201x1201*2 = 2,884,802)
            if _is_nonempty_file(hgt_path, min_bytes=2_000_000):
                return tile_id, "hgt", hgt_path
            hgt_path.unlink(missing_ok=True)
        except Exception:
            gz_path.unlink(missing_ok=True)
            hgt_path.unlink(missing_ok=True)

    # Fallback: USGS 3DEP National Map ImageServer 1x1 deg export (CONUS/US territories)
    if 24 <= lat_floor <= 50 and -126 <= lon_floor <= -65:
        params = urllib.parse.urlencode(
            {
                "bbox": f"{lon_floor},{lat_floor},{lon_floor + 1},{lat_floor + 1}",
                "bboxSR": "4326",
                "imageSR": "4326",
                "size": "3600,3600",
                "format": "tiff",
                "pixelType": "F32",
                "noDataInterpretation": "esriNoDataMatchAny",
                "interpolation": "+RSP_BilinearInterpolation",
                "f": "image",
            }
        )
        usgs_url = f"{USGS_3DEP_EXPORT_URL}?{params}"
        tif_dest = dem_dir / f"{tile_id}.tif"
        if _atomic_download(usgs_url, tif_dest, timeout=timeout, retries=2, min_bytes=4096):
            if _is_valid_tiff(tif_dest):
                return tile_id, "tif", tif_dest
            tif_dest.unlink(missing_ok=True)

    # Check if tile was an HTTP 404 open-ocean cell on Skadi outside land coverage
    try:
        req = urllib.request.Request(skadi_url, method="HEAD", headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=15.0):
            pass
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return tile_id, "void", None
    except Exception:
        pass

    return tile_id, "failed", None


def _convert_hgt_batch_to_geotiff(
    hgt_paths: Sequence[Path],
    dem_dir: Path,
    rscript_bin: Optional[str] = None,
) -> int:
    """Convert staged .hgt tiles into DEFLATE-compressed EPSG:4326 GeoTIFFs via R terra (or keep .hgt)."""
    if not hgt_paths:
        return 0

    r_bin = rscript_bin or shutil.which("Rscript")
    if r_bin:
        r_conv = """
args <- commandArgs(trailingOnly = TRUE)
dem_dir  <- args[1]
hgt_list <- args[-1]

suppressPackageStartupMessages({
  library(terra)
})

n_ok <- 0L
for (h_path in hgt_list) {
  if (!file.exists(h_path)) next
  stem <- sub("\\\\.hgt$", "", basename(h_path), ignore.case = TRUE)
  out_tif <- file.path(dem_dir, paste0(stem, ".tif"))
  tmp_tif <- file.path(dem_dir, paste0(".", stem, ".tmp.tif"))
  ok <- tryCatch({
    r <- terra::rast(h_path)
    if (!nzchar(terra::crs(r))) terra::crs(r) <- "EPSG:4326"
    terra::NAflag(r) <- -32768
    terra::writeRaster(
      r,
      tmp_tif,
      overwrite = TRUE,
      datatype = "INT2S",
      NAflag = -32768,
      gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")
    )
    file.exists(tmp_tif) && file.info(tmp_tif)$size > 1024
  }, error = function(e) FALSE)
  if (isTRUE(ok)) {
    if (file.exists(out_tif)) unlink(out_tif)
    if (!file.rename(tmp_tif, out_tif)) {
      file.copy(tmp_tif, out_tif, overwrite = TRUE)
      unlink(tmp_tif)
    }
    n_ok <- n_ok + 1L
  } else if (file.exists(tmp_tif)) {
    unlink(tmp_tif)
  }
}
cat(sprintf("CONVERTED_TIFS=%d\\n", n_ok))
"""
        with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as rf:
            rf.write(r_conv)
            r_script_path = Path(rf.name)
        try:
            proc = subprocess.run(
                [r_bin, str(r_script_path), str(dem_dir), *[str(p) for p in hgt_paths]],
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode == 0:
                converted = sum(
                    1
                    for p in hgt_paths
                    if _is_valid_tiff(dem_dir / f"{p.stem}.tif")
                )
                if converted == len(hgt_paths):
                    return converted
        finally:
            r_script_path.unlink(missing_ok=True)

    # Fallback: atomically move uncompressed .hgt directly into dem_dir (GDAL/terra reads .hgt natively)
    kept = 0
    for p in hgt_paths:
        tif_p = dem_dir / f"{p.stem}.tif"
        if _is_valid_tiff(tif_p):
            kept += 1
            continue
        dst_hgt = dem_dir / p.name
        shutil.copy2(p, dst_hgt)
        if _is_nonempty_file(dst_hgt, min_bytes=2_000_000):
            kept += 1
    return kept


def inspect_dem_3dep_cache(
    dem_dir: Path,
    required_tiles: Set[Tuple[int, int]],
) -> Tuple[bool, List[Tuple[int, int]]]:
    """Return (all_ready, missing_tiles) for 30 m DEM coverage in dem_dir."""
    if not dem_dir.is_dir():
        return (len(required_tiles) == 0), sorted(required_tiles)

    # Purge legacy 1 km placeholder if present
    stale_placeholder = dem_dir / "dem_3dep_regional.tif"
    if stale_placeholder.exists():
        stale_placeholder.unlink(missing_ok=True)

    # If user placed a custom monolithic regional DEM (not named NxxWyyy.tif), accept it when no tile list exists
    existing_tifs = [
        p
        for p in dem_dir.iterdir()
        if p.suffix.lower() in (".tif", ".tiff", ".hgt")
        and p.name.lower() != "dem_3dep_regional.tif"
        and _is_nonempty_file(p)
    ]
    if not required_tiles:
        return (len(existing_tifs) > 0), []

    tile_names_on_disk = {
        p.stem.upper()
        for p in existing_tifs
        if (p.suffix.lower() == ".hgt" or _is_valid_tiff(p))
    }
    has_monolithic_custom = any(
        not re.match(r"^[NS]\d{2}[EW]\d{3}$", p.stem, flags=re.I) and _is_valid_tiff(p)
        for p in existing_tifs
    )
    if has_monolithic_custom and not tile_names_on_disk:
        return True, []

    void_tiles = _load_void_tiles_memo(dem_dir)
    missing: List[Tuple[int, int]] = []
    for lat_f, lon_f in sorted(required_tiles):
        t_id = format_dem_tile_id(lat_f, lon_f)
        if t_id not in tile_names_on_disk and t_id not in void_tiles:
            missing.append((lat_f, lon_f))

    return (len(missing) == 0), missing


def ensure_dem_3dep_layers(
    env_dir: Path,
    vouchers_csv: Optional[Path] = None,
    workspace_root: Path = PROJECT_ROOT,
    db_path: Optional[Path] = None,
    force: bool = False,
    max_workers: int = 8,
    rscript_bin: Optional[str] = None,
    timeout: float = 90.0,
) -> Dict[str, Any]:
    """Ensure 1-arc-second (~30 m) DEM tiles covering all study voucher coordinates exist."""
    dem_dir = env_dir / "dem_3dep_30m"
    dem_dir.mkdir(parents=True, exist_ok=True)

    required_tiles = collect_voucher_dem_tiles(
        vouchers_csv=vouchers_csv,
        workspace_root=workspace_root,
        db_path=db_path,
    )
    if not required_tiles:
        ready_any, _ = inspect_dem_3dep_cache(dem_dir, set())
        return {
            "status": "cached" if ready_any else "deferred_no_vouchers",
            "downloaded": 0,
            "required_tiles": 0,
            "dir": str(dem_dir),
        }

    if not force:
        ready, missing_tiles = inspect_dem_3dep_cache(dem_dir, required_tiles)
        if ready:
            return {
                "status": "cached",
                "downloaded": 0,
                "required_tiles": len(required_tiles),
                "dir": str(dem_dir),
            }
    else:
        missing_tiles = sorted(required_tiles)

    logger.info(
        f"Downloading {len(missing_tiles)} missing 1x1-deg ~30 m DEM tile(s) "
        f"(out of {len(required_tiles)} cohort cell(s)) into {dem_dir}..."
    )

    void_tiles = _load_void_tiles_memo(dem_dir)
    staged_hgts: List[Path] = []
    direct_tifs = 0
    failed_tiles: List[str] = []

    with tempfile.TemporaryDirectory(prefix="packera_3dep_") as tmp_dir:
        staging_path = Path(tmp_dir)
        workers = max(1, min(int(max_workers), len(missing_tiles)))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(
                    _download_single_3dep_tile,
                    lat_f,
                    lon_f,
                    staging_path,
                    dem_dir,
                    timeout,
                )
                for lat_f, lon_f in missing_tiles
            ]
            for fut in concurrent.futures.as_completed(futures):
                t_id, outcome, p_out = fut.result()
                if outcome == "hgt" and p_out is not None:
                    staged_hgts.append(p_out)
                elif outcome == "tif":
                    direct_tifs += 1
                elif outcome == "void":
                    void_tiles.add(t_id)
                else:
                    failed_tiles.append(t_id)

        converted = _convert_hgt_batch_to_geotiff(
            staged_hgts,
            dem_dir=dem_dir,
            rscript_bin=rscript_bin,
        )

    if void_tiles:
        _save_void_tiles_memo(dem_dir, void_tiles)

    total_dl = converted + direct_tifs
    status = "downloaded" if not failed_tiles else ("partial" if total_dl > 0 else "failed")
    return {
        "status": status,
        "downloaded": total_dl,
        "required_tiles": len(required_tiles),
        "failed": sorted(failed_tiles),
        "dir": str(dem_dir),
    }


# =============================================================================
# Unified Idempotent Entry Point
# =============================================================================

def ensure_environmental_layers(
    env_dir: Path | str = PROJECT_ROOT / "data" / "environmental",
    vouchers_csv: Optional[Path | str] = None,
    workspace_root: Path = PROJECT_ROOT,
    db_path: Optional[Path] = None,
    force: bool = False,
    skip: bool = False,
    demo: bool = False,
    max_workers: int = 6,
    rscript_bin: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Verify that SOLUS100, USGS SGMC, and 30 m DEM layers exist on disk; if all
    present, skip immediately with zero network traffic. Otherwise download
    only the missing layers atomically.
    """
    if (
        demo
        or skip
        or os.environ.get("PACKERA_DEMO_MODE") == "1"
        or os.environ.get("PACKERA_SKIP_ENV_DOWNLOAD") == "1"
    ):
        return {"status": "skipped", "layers": {}}

    env_path = Path(env_dir)
    if not env_path.is_absolute():
        env_path = (workspace_root / env_path).resolve()
    vouchers_path = Path(vouchers_csv) if vouchers_csv is not None else None

    # Fast-path < 50 ms cache check across all 3 layer families before any network or R startup
    if not force:
        solus_ready, _ = inspect_solus_cache(env_path / "solus_100m")
        sgmc_existing = inspect_sgmc_cache(env_path / "sgmc", env_dir=env_path)
        req_tiles = collect_voucher_dem_tiles(
            vouchers_csv=vouchers_path,
            workspace_root=workspace_root,
            db_path=db_path,
        )
        dem_ready, _ = inspect_dem_3dep_cache(env_path / "dem_3dep_30m", req_tiles)

        if solus_ready and sgmc_existing is not None and dem_ready:
            logger.info(
                f"Verified cached SOLUS100, USGS SGMC, and 3DEP 30 m DEM layers in {env_path} "
                f"( Cohort DEM tiles: {len(req_tiles)}; skipping download )."
            )
            return {
                "status": "cached",
                "layers": {
                    "solus_100m": {"status": "cached"},
                    "sgmc": {"status": "cached", "path": str(sgmc_existing)},
                    "dem_3dep_30m": {"status": "cached", "required_tiles": len(req_tiles)},
                },
            }

    env_path.mkdir(parents=True, exist_ok=True)
    solus_res = ensure_solus100_layers(
        env_path,
        force=force,
        max_workers=min(4, max(1, int(max_workers))),
    )
    sgmc_res = ensure_sgmc_layer(
        env_path,
        force=force,
        rscript_bin=rscript_bin,
    )
    dem_res = ensure_dem_3dep_layers(
        env_path,
        vouchers_csv=vouchers_path,
        workspace_root=workspace_root,
        db_path=db_path,
        force=force,
        max_workers=max(1, int(max_workers)),
        rscript_bin=rscript_bin,
    )

    statuses = {solus_res["status"], sgmc_res["status"], dem_res["status"]}
    overall = (
        "cached"
        if statuses == {"cached"}
        else ("ready" if statuses <= {"cached", "downloaded", "deferred_no_vouchers"} else "partial")
    )
    logger.info(
        f"Summary: SOLUS100={solus_res['status']}, "
        f"USGS_SGMC={sgmc_res['status']}, "
        f"DEM_3DEP_30M={dem_res['status']}."
    )
    return {
        "status": overall,
        "layers": {
            "solus_100m": solus_res,
            "sgmc": sgmc_res,
            "dem_3dep_30m": dem_res,
        },
    }


def build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="download_environmental_layers.py",
        description="Idempotent downloader for SOLUS100 (100 m), USGS SGMC, and 3DEP/Skadi 30 m DEM layers.",
    )
    parser.add_argument(
        "--env-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "environmental",
        help="Target environmental root directory (default: data/environmental).",
    )
    parser.add_argument(
        "--vouchers",
        type=Path,
        default=None,
        help="Optional curated_vouchers.csv path used to scope 1x1-degree 30 m DEM tiles.",
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=None,
        help="Optional override path to SQLite database.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force re-download even when cached layers already exist.",
    )
    parser.add_argument(
        "--skip",
        action="store_true",
        default=False,
        help="Skip all downloads immediately (no-op).",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        default=False,
        help="Skip downloads in ephemeral --demo mode.",
    )
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=6,
        help="Concurrent download worker threads (default: 6).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        default=False,
        help="Exit with non-zero status if any missing environmental layer fails to download.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_cli_parser()
    args = parser.parse_args(argv)
    summary = ensure_environmental_layers(
        env_dir=args.env_dir,
        vouchers_csv=args.vouchers,
        workspace_root=PROJECT_ROOT,
        db_path=args.db,
        force=args.force,
        skip=args.skip,
        demo=args.demo,
        max_workers=args.workers,
    )
    if args.strict and summary["status"] not in ("cached", "ready", "skipped"):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
