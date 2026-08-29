from pathlib import Path
from typing import Dict, List, Tuple

# Multi-class schema mapping for botanical vision model
CLASS_NAMES: List[str] = [
    "basal_leaf",       # 0: Leaf blade / intact leaf
    "leaf_petiole",     # 1: Distinct petiole / leaf stalk
    "basal_rosette",    # 2: Clustered basal rosette
    "capitulum",        # 3: Inflorescence / flower head
    "herbarium_label",  # 4: Main specimen metadata label
    "color_chart",      # 5: Calibration color chart / palette
    "ruler_scale",      # 6: Measurement scale / centimeter bar
    "barcode_sticker",  # 7: Digitization barcode / QR sticker
    "mounting_tape",    # 8: Linen, paper, or plastic mounting tape strip
]

CLASS_MAP: Dict[str, int] = {name: idx for idx, name in enumerate(CLASS_NAMES)}

# Distinct RGB color palette for rendering class overlays during QC
CLASS_COLORS_BGR: Dict[int, Tuple[int, int, int]] = {
    0: (0, 200, 0),      # basal_leaf: Vibrant Green
    1: (50, 255, 150),   # leaf_petiole: Mint Green
    2: (0, 140, 70),     # basal_rosette: Dark Forest Green
    3: (0, 215, 255),    # capitulum: Gold / Yellow
    4: (30, 30, 230),    # herbarium_label: Bright Red
    5: (230, 30, 230),   # color_chart: Magenta
    6: (0, 140, 255),    # ruler_scale: Orange
    7: (230, 180, 0),    # barcode_sticker: Cyan
    8: (180, 0, 180),    # mounting_tape: Purple / Violet
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
