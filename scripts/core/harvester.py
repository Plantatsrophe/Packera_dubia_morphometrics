import os
import sys
import time
import argparse
import asyncio
import aiohttp
import logging
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Union
from collections import defaultdict
import numpy as np
import pandas as pd
import pygbif.occurrences as occ
from tqdm.asyncio import tqdm as async_tqdm

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Configuration path and taxonomy defaults
from scripts.core.config import (
    DEFAULT_WORKSPACE,
    DEFAULT_RAW_DIR,
    DEFAULT_TARGET_TAXA,
    DEFAULT_OUTPUT_CSV,
    DEFAULT_SUMMARY_LOG,
)

# Harvester parsing, scoring, downloading, and logging utilities
from scripts.core.harvester_utils import (
    setup_logger,
    sanitize_filename,
    parse_determiner_tier,
    calculate_circular_phenology,
    infer_regional_group,
    extract_high_res_image_url,
    print_and_log_summary,
    download_all_voucher_images,
)

def harvest_taxa_occurrences(
    taxa_list: List[str],
    max_uncertainty_meters: float = 5000.0,
    max_records_per_taxon: int = 1000,
    logger: Optional[logging.Logger] = None
) -> pd.DataFrame:
    """
    Executes paginated queries to the GBIF Occurrence API and filters records according to DwC standards.
    
    Filters applied:
      - Basis of record: PRESERVED_SPECIMEN
      - Geographic scope: United States (country='US')
      - Spatial validation: Non-null coordinates, coordinateUncertaintyInMeters <= max_uncertainty_meters
      - Temporal validation: Valid collection dates (year, month, day)
      - Media validation: High-resolution sheet image URL present
      
    Args:
        taxa_list: List of botanical species names to query.
        max_uncertainty_meters: Maximum allowed georeferencing uncertainty in meters.
        max_records_per_taxon: Maximum records to harvest per taxon query.
        logger: Logger instance.
        
    Returns:
        pd.DataFrame: Curated dataframe containing parsed metadata and authority scores.
    """
    if logger is None:
        logger = logging.getLogger("VoucherHarvester")

    all_curated_records = []
    seen_catalog_keys = set()

    for taxon in taxa_list:
        logger.info(f"Querying GBIF API for taxon: '{taxon}' (country=US, basisOfRecord=PRESERVED_SPECIMEN)...")
        
        offset = 0
        limit = 300  # GBIF page size limit
        taxon_harvested = 0
        taxon_retained = 0

        while taxon_retained < max_records_per_taxon:
            try:
                # Query GBIF occurrence search endpoint
                response = occ.search(
                    scientificName=taxon,
                    country="US",
                    basisOfRecord="PRESERVED_SPECIMEN",
                    limit=limit,
                    offset=offset
                )
            except Exception as e:
                logger.error(f"GBIF API query error for '{taxon}' at offset {offset}: {e}")
                break

            results = response.get("results", [])
            count = response.get("count", 0)
            if not results:
                break

            for rec in results:
                taxon_harvested += 1

                # 1. Geographic Coordinate & Uncertainty Validation
                lat = rec.get("decimalLatitude")
                lon = rec.get("decimalLongitude")
                
                lat_val = None
                lon_val = None
                if lat is not None and lon is not None:
                    try:
                        lat_val = float(lat)
                        lon_val = float(lon)
                    except (ValueError, TypeError):
                        pass

                # Coordinate uncertainty check (only if coordinates exist)
                uncertainty_val = max_uncertainty_meters
                if lat_val is not None and lon_val is not None:
                    uncertainty_raw = rec.get("coordinateUncertaintyInMeters")
                    if uncertainty_raw is not None:
                        try:
                            parsed_unc = float(uncertainty_raw)
                            if parsed_unc > max_uncertainty_meters:
                                continue
                            uncertainty_val = parsed_unc
                        except (ValueError, TypeError):
                            pass

                # 2. Temporal & Phenology Validation (year, month, day)
                year = rec.get("year")
                month = rec.get("month")
                day = rec.get("day")
                pheno_res = calculate_circular_phenology(year, month, day)
                if pheno_res is None:
                    continue
                doy, pheno_sin, pheno_cos = pheno_res

                # 3. High-Resolution Media Image Validation
                media_list = rec.get("media", [])
                image_url = extract_high_res_image_url(media_list)
                if not image_url:
                    continue

                # 4. Catalog Number & Herbarium Institution Normalization
                raw_catalog = rec.get("catalogNumber")
                inst_code = rec.get("institutionCode") or rec.get("collectionCode") or "UNKNOWN_INST"
                gbif_key = str(rec.get("key", ""))

                if raw_catalog and str(raw_catalog).strip():
                    catalog_number = sanitize_filename(str(raw_catalog).strip())
                else:
                    catalog_number = f"{sanitize_filename(inst_code)}_{gbif_key}"

                # Ensure unique catalog identifier across duplicate uploads
                unique_key = (catalog_number, gbif_key)
                if unique_key in seen_catalog_keys:
                    continue
                seen_catalog_keys.add(unique_key)

                # Local image path destination
                relative_image_path = f"data/raw_vouchers/{catalog_number}.jpg"

                # 5. Taxonomic Authority & Determiner Tier Stratification
                type_status_raw = rec.get("typeStatus")
                identified_by_raw = rec.get("identifiedBy")
                recorded_by_raw = rec.get("recordedBy")
                history_raw = rec.get("verbatimIdentificationHistory")
                locality_raw = rec.get("locality") or rec.get("verbatimLocality")
                habitat_raw = rec.get("habitat")
                species_raw = rec.get("scientificName") or rec.get("species") or taxon

                determiner_tier, type_status, determiner_raw = parse_determiner_tier(
                    type_status_raw=type_status_raw,
                    identified_by_raw=identified_by_raw,
                    recorded_by_raw=recorded_by_raw,
                    history_raw=history_raw,
                    institution_code_raw=inst_code,
                    locality_raw=locality_raw,
                    habitat_raw=habitat_raw
                )

                # 6. Regional Ecological Group Assignment
                state_prov = rec.get("stateProvince")
                regional_group = infer_regional_group(
                    lat=lat_val,
                    lon=lon_val,
                    state_province=state_prov,
                    habitat=habitat_raw,
                    locality=locality_raw
                )

                # Append standardized record
                curated_record = {
                    "catalogNumber": catalog_number,
                    "institutionCode": inst_code,
                    "species_raw": species_raw,
                    "determiner_raw": determiner_raw,
                    "determiner_tier": determiner_tier,
                    "type_status": type_status,
                    "county": rec.get("county") or "",
                    "stateProvince": state_prov or "",
                    "latitude": lat_val,
                    "longitude": lon_val,
                    "coordinateUncertainty": uncertainty_val,
                    "doy": doy,
                    "pheno_sin": pheno_sin,
                    "pheno_cos": pheno_cos,
                    "regional_group": regional_group,
                    "image_path": relative_image_path,
                    "_image_url": image_url,  # Temporary internal column for downloading
                }
                all_curated_records.append(curated_record)
                taxon_retained += 1
                
                if taxon_retained >= max_records_per_taxon:
                    break

            offset += len(results)
            if offset >= count:
                break
            time.sleep(0.1)  # Respectful GBIF rate limiting

        logger.info(f"Taxon '{taxon}': Processed {taxon_harvested} occurrences -> Retained {taxon_retained} curated records meeting quality filters.")

    # Convert to DataFrame
    df = pd.DataFrame(all_curated_records)
    return df


def main() -> None:
    """
    Main entry point for command-line execution of the voucher harvester.
    """
    parser = argparse.ArgumentParser(
        description="Automated GBIF Voucher Harvester & Determiner Authority Scorer for Packera dubia."
    )
    parser.add_argument(
        "--taxa",
        nargs="+",
        default=DEFAULT_TARGET_TAXA,
        help="List of scientific binomials or names to harvest from GBIF."
    )
    parser.add_argument(
        "--max-uncertainty",
        type=float,
        default=5000.0,
        help="Maximum coordinate uncertainty in meters (default: 5000.0)."
    )
    parser.add_argument(
        "--max-records-per-taxon",
        type=int,
        default=5000,
        help="Maximum records to harvest per taxon (default: 5000)."
    )
    parser.add_argument(
        "--download-images",
        action="store_true",
        default=False,
        help="Flag to enable asynchronous high-resolution specimen image downloading."
    )
    parser.add_argument(
        "--concurrency",
        type=int,
        default=15,
        help="Max concurrent asynchronous image downloads (default: 15)."
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default=str(DEFAULT_OUTPUT_CSV),
        help=f"Target path for curated metadata CSV (default: {DEFAULT_OUTPUT_CSV})."
    )
    parser.add_argument(
        "--log-file",
        type=str,
        default=str(DEFAULT_SUMMARY_LOG),
        help=f"Path for summary log file (default: {DEFAULT_SUMMARY_LOG})."
    )

    args = parser.parse_args()

    # Create destination directories
    output_csv_path = Path(args.output_csv)
    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_RAW_DIR.mkdir(parents=True, exist_ok=True)
    log_file_path = Path(args.log_file)
    log_file_path.parent.mkdir(parents=True, exist_ok=True)

    # Initialize logger
    logger = setup_logger(log_file_path)
    logger.info("Starting Packera Voucher Ingestion & Authority Stratification Pipeline...")
    logger.info(f"Target Taxa: {args.taxa}")
    logger.info(f"Max Coordinate Uncertainty Threshold: {args.max_uncertainty} m")

    # Step 1: Harvest and curate metadata records from GBIF
    df_curated = harvest_taxa_occurrences(
        taxa_list=args.taxa,
        max_uncertainty_meters=args.max_uncertainty,
        max_records_per_taxon=args.max_records_per_taxon,
        logger=logger
    )

    download_stats = None

    # Step 2: Asynchronously download specimen images if enabled
    if args.download_images and not df_curated.empty and "_image_url" in df_curated.columns:
        download_queue = []
        for _, row in df_curated.iterrows():
            url = row["_image_url"]
            # Resolve image destination path relative to project workspace
            dest = DEFAULT_WORKSPACE / row["image_path"]
            download_queue.append((url, dest))

        logger.info(f"Starting asynchronous download of {len(download_queue)} voucher sheets...")
        download_stats = asyncio.run(
            download_all_voucher_images(
                records_to_download=download_queue,
                concurrency_limit=args.concurrency,
                logger=logger
            )
        )

    # Step 3: Export standardized CSV (excluding internal temporary columns)
    export_columns = [
        "catalogNumber",
        "institutionCode",
        "species_raw",
        "determiner_raw",
        "determiner_tier",
        "type_status",
        "county",
        "stateProvince",
        "latitude",
        "longitude",
        "coordinateUncertainty",
        "doy",
        "pheno_sin",
        "pheno_cos",
        "regional_group",
        "image_path"
    ]
    
    if not df_curated.empty:
        df_export = df_curated[[col for col in export_columns if col in df_curated.columns]]
        df_export.to_csv(output_csv_path, index=False, encoding="utf-8")
        logger.info(f"Successfully exported {len(df_export)} curated records to: {output_csv_path}")
    else:
        # Create empty table with standardized headers
        pd.DataFrame(columns=export_columns).to_csv(output_csv_path, index=False, encoding="utf-8")
        logger.warning(f"No records met all filtering criteria. Created empty table at: {output_csv_path}")

    # Step 4: Output comprehensive summary log
    print_and_log_summary(df_curated, download_stats, logger)
    logger.info("Pipeline execution completed successfully.")


