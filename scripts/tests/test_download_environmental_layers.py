#!/usr/bin/env python3
"""
===============================================================================
Script: test_download_environmental_layers.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration tests for scripts/data_prep/download_environmental_layers.py
    and its startup wiring in run_pipeline.sh and main.py. Verifies:
      1. Fast-skip idempotency (zero HTTP calls when SOLUS100, SGMC, and 3DEP exist).
      2. Automatic download & atomic staging when any layer directory is empty.
      3. Incremental 1x1-degree 3DEP DEM tile provisioning across voucher cohorts.
      4. ScienceBase JSON catalog parsing & R relational join invocation for SGMC.
      5. Demo/skip flag isolation and CLI wiring in run_pipeline.sh and main.py.
===============================================================================
"""

from __future__ import annotations

import gzip
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import MagicMock, patch

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import main as root_main
from scripts.data_prep import download_environmental_layers as env_dl


def _fake_tiff_bytes(size: int = 2048) -> bytes:
    # Little-endian classic TIFF magic header II*\x00 + padding
    return b"II*\x00" + (b"\x00" * max(0, size - 4))


class TestDownloadEnvironmentalLayers(unittest.TestCase):
    """Tests for idempotent SOLUS100, USGS SGMC, and 3DEP 30 m DEM provisioning."""

    def test_format_dem_tile_id_and_voucher_tile_collection(self) -> None:
        self.assertEqual(env_dl.format_dem_tile_id(35, -79), "N35W079")
        self.assertEqual(env_dl.format_dem_tile_id(29, -96), "N29W096")

        with tempfile.TemporaryDirectory() as tmp_d:
            csv_path = Path(tmp_d) / "curated_vouchers.csv"
            csv_path.write_text(
                "catalogNumber,decimalLatitude,decimalLongitude\n"
                "NCU001,35.91,-79.05\n"
                "NCU002,35.12,-79.88\n"
                "GA001,33.42,-82.15\n"
                "BAD_OCEAN,0.0,0.0\n",
                encoding="utf-8",
            )
            tiles = env_dl.collect_voucher_dem_tiles(vouchers_csv=csv_path, workspace_root=Path(tmp_d))
            self.assertEqual(tiles, {(35, -80), (33, -83)})

    def test_fast_skip_zero_network_when_all_layers_exist(self) -> None:
        """When SOLUS100, SGMC, and cohort 3DEP tiles exist on disk, zero HTTP or R calls occur."""
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_root = Path(tmp_d)
            env_dir = tmp_root / "data" / "environmental"
            solus_dir = env_dir / "solus_100m"
            sgmc_dir = env_dir / "sgmc"
            dem_dir = env_dir / "dem_3dep_30m"
            solus_dir.mkdir(parents=True)
            sgmc_dir.mkdir(parents=True)
            dem_dir.mkdir(parents=True)

            for spec in env_dl.SOLUS100_REQUIRED_LAYERS:
                (solus_dir / str(spec["filename"])).write_bytes(_fake_tiff_bytes())

            (sgmc_dir / "sgmc_conus_lithology.gpkg").write_bytes(b"SQLite format 3\x00" + b"\x00" * 2048)

            vouchers_csv = tmp_root / "curated_vouchers.csv"
            vouchers_csv.write_text(
                "catalogNumber,decimalLatitude,decimalLongitude\n"
                "V1,35.8,-79.1\n"
                "V2,34.2,-80.5\n",
                encoding="utf-8",
            )
            (dem_dir / "N35W080.tif").write_bytes(_fake_tiff_bytes())
            (dem_dir / "N34W081.tif").write_bytes(_fake_tiff_bytes())

            with patch.object(env_dl.urllib.request, "urlopen") as mock_http, \
                 patch.object(env_dl.subprocess, "run") as mock_subp:
                res = env_dl.ensure_environmental_layers(
                    env_dir=env_dir,
                    vouchers_csv=vouchers_csv,
                    workspace_root=tmp_root,
                    force=False,
                )
                self.assertEqual(res["status"], "cached")
                mock_http.assert_not_called()
                mock_subp.assert_not_called()

    def test_downloads_missing_layers_and_skips_already_cached_tiles(self) -> None:
        """When folders are empty or partially populated, only missing files are downloaded atomically."""
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_root = Path(tmp_d)
            env_dir = tmp_root / "data" / "environmental"
            dem_dir = env_dir / "dem_3dep_30m"
            dem_dir.mkdir(parents=True)

            # Pre-populate 1 of the 2 needed DEM tiles to verify incremental skip
            (dem_dir / "N35W080.tif").write_bytes(_fake_tiff_bytes())

            vouchers_csv = tmp_root / "curated_vouchers.csv"
            vouchers_csv.write_text(
                "catalogNumber,decimalLatitude,decimalLongitude\n"
                "V1,35.8,-79.1\n"
                "V2,32.6,-83.4\n",
                encoding="utf-8",
            )

            # Build fake SGMC zip bytes
            zip_buf = io.BytesIO()
            with zipfile.ZipFile(zip_buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
                zf.writestr("SGMC_Geology.shp", os.urandom(96 * 1024))
                zf.writestr("SGMC_Units.csv", "UNIT_LINK,MAJOR1,UNIT_NAME\nNC001,Granite,Rolesville\n")
            fake_sgmc_zip = zip_buf.getvalue()

            # Build fake Skadi .hgt.gz bytes (>= 2 MB uncompressed)
            raw_hgt = b"\x00\x64" * 1_100_000
            gz_buf = io.BytesIO()
            with gzip.GzipFile(fileobj=gz_buf, mode="wb") as gz_f:
                gz_f.write(raw_hgt)
            fake_hgt_gz = gz_buf.getvalue()

            requested_urls: list[str] = []

            class _FakeResp:
                def __init__(self, payload: bytes, status: int = 200) -> None:
                    self._bio = io.BytesIO(payload)
                    self.status = status
                    self.headers = {}

                def read(self, n: int = -1) -> bytes:
                    return self._bio.read(n)

                def __enter__(self) -> "_FakeResp":
                    return self

                def __exit__(self, *args: object) -> None:
                    pass

            def _fake_urlopen(req, timeout: float = 60.0):
                url = req.full_url if hasattr(req, "full_url") else str(req)
                requested_urls.append(url)
                if "sciencebase.gov/catalog/item" in url:
                    meta = {
                        "files": [
                            {
                                "name": "USGS_SGMC_Shapefiles.zip",
                                "downloadUri": "https://www.sciencebase.gov/catalog/file/get/5888bf4fe4b05ccb964bab9d?f=shp.zip",
                            }
                        ]
                    }
                    return _FakeResp(json.dumps(meta).encode("utf-8"))
                if "f=shp.zip" in url:
                    return _FakeResp(fake_sgmc_zip)
                if "solus100pub/" in url:
                    return _FakeResp(_fake_tiff_bytes(4096))
                if "elevation-tiles-prod/skadi/" in url:
                    return _FakeResp(fake_hgt_gz)
                raise AssertionError(f"Unexpected URL requested: {url}")

            def _fake_rscript_run(cmd, **kwargs):
                # Emulate Rscript writing sgmc_conus_lithology.gpkg or converting .hgt -> .tif
                r_file = Path(cmd[1]).read_text(encoding="utf-8")
                if "sgmc_lithology" in r_file:
                    out_gpkg = Path(cmd[3])
                    out_gpkg.parent.mkdir(parents=True, exist_ok=True)
                    out_gpkg.write_bytes(b"SQLite format 3\x00" + b"\x00" * 4096)
                elif "CONVERTED_TIFS" in r_file:
                    out_dem_dir = Path(cmd[2])
                    for h_arg in cmd[3:]:
                        hp = Path(h_arg)
                        (out_dem_dir / f"{hp.stem}.tif").write_bytes(_fake_tiff_bytes(4096))
                res = MagicMock()
                res.returncode = 0
                res.stdout = "OK\n"
                res.stderr = ""
                return res

            with patch.object(env_dl.urllib.request, "urlopen", side_effect=_fake_urlopen), \
                 patch.object(env_dl.subprocess, "run", side_effect=_fake_rscript_run):
                summary = env_dl.ensure_environmental_layers(
                    env_dir=env_dir,
                    vouchers_csv=vouchers_csv,
                    workspace_root=tmp_root,
                    rscript_bin="/usr/bin/Rscript",
                )

            self.assertEqual(summary["status"], "ready")
            # Verify cached N35W080 was NOT re-downloaded, while missing N32W084 WAS downloaded
            self.assertFalse(any("N35W080" in u for u in requested_urls))
            self.assertTrue(any("N32W084.hgt.gz" in u for u in requested_urls))
            self.assertTrue((dem_dir / "N32W084.tif").is_file())
            self.assertTrue((env_dir / "sgmc" / "sgmc_conus_lithology.gpkg").is_file())
            self.assertTrue((env_dir / "solus_100m" / "anylithicdpt_cm_p.tif").is_file())
            self.assertTrue((env_dir / "solus_100m" / "ph1to1h2o_5_cm_p.tif").is_file())

            # Second invocation must hit fast-skip with zero HTTP calls
            requested_urls.clear()
            with patch.object(env_dl.urllib.request, "urlopen", side_effect=_fake_urlopen):
                second = env_dl.ensure_environmental_layers(
                    env_dir=env_dir,
                    vouchers_csv=vouchers_csv,
                    workspace_root=tmp_root,
                    rscript_bin="/usr/bin/Rscript",
                )
            self.assertEqual(second["status"], "cached")
            self.assertEqual(requested_urls, [])

    def test_demo_and_skip_flags_bypass_downloads_and_cli_wiring(self) -> None:
        with patch.object(env_dl.urllib.request, "urlopen") as mock_http:
            res_demo = env_dl.ensure_environmental_layers(demo=True)
            res_skip = env_dl.ensure_environmental_layers(skip=True)
            self.assertEqual(res_demo["status"], "skipped")
            self.assertEqual(res_skip["status"], "skipped")
            mock_http.assert_not_called()

        sh_src = (PROJECT_ROOT / "run_pipeline.sh").read_text(encoding="utf-8")
        self.assertIn("run_preflight_env_layers()", sh_src)
        self.assertIn("scripts/data_prep/download_environmental_layers.py", sh_src)
        self.assertIn("--skip-env-download", sh_src)
        self.assertIn("--force-env-download", sh_src)

        parser = root_main.build_parser()
        dl_args = parser.parse_args(["download-env", "--force"])
        self.assertEqual(dl_args.subcommand, "download-env")
        self.assertTrue(dl_args.force)

        synth_args = parser.parse_args(["synthesis", "--skip-env-download"])
        self.assertTrue(synth_args.skip_env_download)


if __name__ == "__main__":
    unittest.main()
