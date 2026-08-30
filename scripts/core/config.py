from pathlib import Path
from typing import Dict, List, Tuple

# Standard 7-class botanical ontology mapping for Packera phenotyping
CLASS_NAMES: List[str] = [
    "basal_leaf_blade",     # 0: laminar portion of basal leaves
    "leaf_petiole",         # 1: narrow petiole stalk connecting caudex to blade
    "cauline_leaf",         # 2: sessile/lyrately-pinnatifid leaves on flowering stalk
    "cauline_stem",         # 3: main vertical flowering stalk / scape
    "root_rhizome",         # 4: dark fibrous subterranean roots and caudex
    "basal_rosette_clump",  # 5: dense overlapping basal rosette crown
    "capitulum",            # 6: inflorescence head / involucre
]

CLASS_MAP: Dict[str, int] = {name: idx for idx, name in enumerate(CLASS_NAMES)}

# Distinct RGB/BGR color palette for rendering class overlays during QC
CLASS_COLORS_BGR: Dict[int, Tuple[int, int, int]] = {
    0: (0, 200, 0),      # basal_leaf_blade: Vibrant Green
    1: (50, 255, 150),   # leaf_petiole: Mint Green
    2: (0, 165, 255),    # cauline_leaf: Bright Orange
    3: (0, 215, 255),    # cauline_stem: Gold / Yellow
    4: (50, 100, 200),   # root_rhizome: Rust / Earth Brown
    5: (0, 140, 70),     # basal_rosette_clump: Dark Forest Green
    6: (230, 30, 230),   # capitulum: Vivid Magenta / Purple
}

# Default filesystem paths
DEFAULT_WORKSPACE = Path(__file__).resolve().parent.parent.parent if "__file__" in locals() else Path.cwd()
DEFAULT_RAW_DIR = DEFAULT_WORKSPACE / "data" / "raw_vouchers"
DEFAULT_CURATED_CSV = DEFAULT_WORKSPACE / "data" / "tables" / "curated_vouchers.csv"
DEFAULT_OUTPUT_CSV = DEFAULT_WORKSPACE / "data" / "tables" / "curated_vouchers.csv"
DEFAULT_SUMMARY_LOG = DEFAULT_WORKSPACE / "outputs" / "reports" / "voucher_ingestion_summary.log"
DEFAULT_OUTPUT_DIR = DEFAULT_WORKSPACE / "data" / "yolo_dataset"
DEFAULT_CONFIG_PATH = DEFAULT_WORKSPACE / "data" / "dataset_config.yaml"
DEFAULT_QC_DIR = DEFAULT_WORKSPACE / "outputs" / "dataset_qc"

# Target focal taxa in the Packera dubia complex and allied North American clades
DEFAULT_TARGET_TAXA: List[str] = [
    "Packera dubia",
    "Packera tomentosa",
    "Senecio tomentosus",
    "Packera anonyma",
    "Packera plattensis",
    "Packera paupercula",
]

# Recognized monographic taxonomic specialists for Packera / Senecioneae (Tier 1 Gold)
SPECIALIST_PATTERNS: List[str] = [
    r"\bBarkley\b",
    r"\bT\.?\s*M\.?\s*Barkley\b",
    r"\bTheodore\s+M\.?\s+Barkley\b",
    r"\bTrock\b",
    r"\bD\.?\s*K\.?\s*Trock\b",
    r"\bDebra\s+K\.?\s+Trock\b",
    r"\bKowal\b",
    r"\bR\.?\s*R\.?\s*Kowal\b",
    r"\bRobert\s+R\.?\s+Kowal\b",
    r"\bWeakley\b",
    r"\bA\.?\s*S\.?\s*Weakley\b",
    r"\bAlan\s+S\.?\s+Weakley\b",
    r"\bBain\b",
    r"\bJ\.?\s*F\.?\s*Bain\b",
    r"\bJohn\s+F\.?\s+Bain\b",
    r"\bMahoney\b",
    r"\bA\.?\s*M\.?\s*Mahoney\b",
    r"\bAlison\s+M\.?\s+Mahoney\b",
    r"\bFuller\b",
    r"\bJ\.?\s*B\.?\s*Fuller\b",
    r"\bBrandon\s+Fuller\b",
]

# Major regional research herbaria with high curatorial standards (Tier 2 Silver)
MAJOR_HERBARIA_CODES = {
    "NCU", "GA", "US", "NY", "BRIT", "MO", "WIS", "VDB", "FLAS", "TEX", "LL", "TENN", "F"
}

# Recognized nomenclatural type status designations (Tier 1 Gold)
VALID_TYPE_STATUSES = {
    "HOLOTYPE", "ISOTYPE", "LECTOTYPE", "ISOLECTOTYPE", "SYNTYPE", "ISOSYNTYPE",
    "NEOTYPE", "ISONEOTYPE", "PARATYPE", "ISOPARATYPE", "EPITYPE", "TYPE", "COTYPE"
}
