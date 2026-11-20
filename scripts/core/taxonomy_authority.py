"""
===============================================================================
Module: taxonomy_authority.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Darwin Core exsiccatae normalization, duplicate stratification, spatial
    balancing/thinning, 3-layer dynamic determiner authority engine, and ICN
    nomenclatural type hierarchy classification.
===============================================================================
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
import re
from typing import Any, Dict, Iterable, List, Mapping, Optional, Set, Tuple, Union
import unicodedata

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]
try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]
try:
    from PIL import Image
except ImportError:
    Image = None  # type: ignore[assignment]

if Image is not None:
    Image.MAX_IMAGE_PIXELS = None

PROJECT_ROOT = Path(__file__).resolve().parents[2]

from scripts.core.config import (
    CIS_MISSISSIPPI_STATES_PROVINCES,
    DEFAULT_CIS_MISSISSIPPI_LON,
    DEFAULT_COUNTRIES,
    DEFAULT_GEOGRAPHIC_SCOPE,
    DEFAULT_HOME_HERBARIUM,
    DEFAULT_MIN_CIS_MISSISSIPPI_FRAC,
    EXCLUDED_WESTERN_STATES,
    GREAT_PLAINS_STRADDLING_STATES,
    MAJOR_HERBARIA_CODES,
    PRIMARY_TYPE_STATUSES,
    SECONDARY_TYPE_STATUSES,
    SPECIALIST_PATTERNS,
    VALID_GEOGRAPHIC_SCOPES,
    WESTERN_LONGITUDE_THRESHOLD,
    normalize_countries_list,
)
from scripts.core.logger import setup_logging
from scripts.core.step01_concept_authority_utils import (
    CORE_PACKERA_GOLD_REGEXES,
    CORE_PACKERA_GOLD_SURNAMES,
    REGIONAL_AND_HISTORICAL_SILVER_REGEXES,
    REGIONAL_AND_HISTORICAL_SILVER_SURNAMES,
    apply_concept_and_segregate_recency_caps,
)

EXPORT_COLUMNS = [
    "catalogNumber",
    "institutionCode",
    "scientificName",
    "species_raw",
    "recordedBy",
    "recordNumber",
    "identifiedBy",
    "determiner_tier",
    "type_status",
    "county",
    "stateProvince",
    "decimalLatitude",
    "decimalLongitude",
    "coordinateUncertainty",
    "year",
    "month",
    "day",
    "eventDate",
    "regional_group",
    "image_path",
    "exsiccatae_key",
    "is_primary_duplicate",
    "duplicate_count",
]


def normalize_collector(recorded_by: Optional[str]) -> str:
    """
    Normalizes botanical collector names for exsiccatae collection event fingerprinting:
    lowercased, non-alphanumeric characters stripped, initials removed, primary surname extracted.
    e.g., 'J. Brandon Fuller' -> 'fuller', 'Fuller, J. B.' -> 'fuller', 'A. Cronquist' -> 'cronquist'.
    """
    if not recorded_by or not isinstance(recorded_by, str):
        return "unknown"
    s = str(recorded_by).strip()
    if not s or s.lower() in {"none", "nan", "null", "unknown", "anonymous", "not evident"}:
        return "unknown"

    # Separate multiple collectors (e.g. 'J. B. Fuller & R. R. Kowal') -> retain primary collector
    primary = re.split(r"\s*(?:&|;|\band\b|\bwith\b|\bet\s+al\.?)\s*", s, flags=re.IGNORECASE)[0].strip()

    # Handle 'Fuller, J. Brandon' format
    if "," in primary:
        surname_part = primary.split(",")[0].strip()
        cleaned = re.sub(r"[^a-zA-Z0-9]", "", surname_part).lower()
        if cleaned:
            return cleaned

    # Tokenize and remove punctuation
    tokens = [re.sub(r"[^a-zA-Z0-9]", "", t).lower() for t in primary.split()]
    tokens = [t for t in tokens if t]

    if not tokens:
        return "unknown"

    # Remove single-letter initials (e.g. 'j', 'b')
    non_initials = [t for t in tokens if len(t) > 1]
    if not non_initials:
        # All tokens were initials (e.g. 'J. B.')
        return "".join(tokens)

    # In 'First Middle Last', the surname is the last token
    return non_initials[-1]


def normalize_record_number(record_number: Optional[Any]) -> str:
    """
    Normalizes collection number: strips boilerplate prefixes ('#', 'No.', 'Coll.'),
    preserves compound/year-prefixed numeric sequences ('2022-014' -> '2022-014') and
    ICN Art. 8.2 split-gathering letter suffixes ('1042b' -> '1042b', '1042-A' -> '1042a').
    """
    if record_number is None:
        return ""
    s = str(record_number).strip()
    if not s or s.lower() in {"none", "nan", "null", "s.n.", "sn", "s. n.", "unnumbered", "sine numero"}:
        return ""
    s_clean = re.sub(r"^(?:coll(?:ection)?\.?\s*)?(?:no\.?|nr\.?|num(?:ber)?\.?|#)\s*", "", s, flags=re.IGNORECASE).strip()
    if re.fullmatch(
        r"[A-Za-z]{2,}\d+(?:[-_./][A-Za-z0-9]+)+|[A-Za-z]{2,}\d+[A-Za-z]{2,}[A-Za-z0-9_./-]*",
        s_clean,
    ):
        return re.sub(r"[./\s_]+", "-", s_clean).lower()
    match = re.search(r"(\d+(?:[-./]\d+)*)(?:[-\s]*([a-zA-Z])(?![a-zA-Z0-9]))?", s_clean)
    if not match:
        return ""
    num_part = re.sub(r"[./]+", "-", match.group(1))
    sfx_part = (match.group(2) or "").lower()
    return f"{num_part}{sfx_part}"


def normalize_event_date(
    event_date_raw: Optional[str] = None,
    year: Optional[Any] = None,
    month: Optional[Any] = None,
    day: Optional[Any] = None,
) -> str:
    """
    Standardizes collection event date to ISO format YYYY-MM-DD.
    Evaluates explicit year, month, day components or parses date string.
    """
    try:
        if year is not None and month is not None and day is not None:
            if not (pd.isna(year) or pd.isna(month) or pd.isna(day)):
                y, m, d = int(year), int(month), int(day)
                if 1700 <= y <= 2100 and 1 <= m <= 12 and 1 <= d <= 31:
                    return f"{y:04d}-{m:02d}-{d:02d}"
    except (ValueError, TypeError):
        pass

    if event_date_raw and isinstance(event_date_raw, str):
        s = event_date_raw.strip()
        match_iso = re.search(r"(\d{4})[-/](\d{1,2})[-/](\d{1,2})", s)
        if match_iso:
            y, m, d = int(match_iso.group(1)), int(match_iso.group(2)), int(match_iso.group(3))
            return f"{y:04d}-{m:02d}-{d:02d}"
        match_y = re.search(r"\b(1[789]\d{2}|20\d{2})\b", s)
        if match_y:
            return f"{match_y.group(1)}-00-00"

    try:
        if year is not None and not pd.isna(year):
            y = int(year)
            if 1700 <= y <= 2100:
                return f"{y:04d}-00-00"
    except (ValueError, TypeError):
        pass

    return "unknown_date"


def generate_exsiccatae_key(
    collector: Optional[str],
    number: Optional[Any],
    date: Optional[str] = None,
    year: Optional[Any] = None,
    month: Optional[Any] = None,
    day: Optional[Any] = None,
) -> str:
    """
    Constructs a standardized collection event signature (exsiccatae fingerprint):
    f"{norm_collector}_{norm_number}_{norm_date}".
    """
    norm_collector = normalize_collector(collector)
    norm_number = normalize_record_number(number)
    norm_date = normalize_event_date(date, year=year, month=month, day=day)
    return f"{norm_collector}_{norm_number}_{norm_date}"


def haversine_distance_meters(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Calculates geodesic distance between two coordinate pairs using the Haversine formula."""
    r_earth = 6371000.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)

    a = math.sin(delta_phi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda / 2.0) ** 2
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
    return r_earth * c


def resolve_voucher_image_path(
    image_path: Union[str, Path],
    workspace_dir: Optional[Path] = None,
) -> Path:
    """
    Resolves a voucher image path on disk by checking the given path, its `.webp` sibling,
    and its `.jpg`/`.jpeg`/`.png` siblings (both relative to `workspace_dir` if provided
    and directly), returning the existing file path on disk (or the `.webp` path by default
    if none exists yet).
    """
    if image_path is None:
        return Path("")
    raw_str = str(image_path).strip()
    if not raw_str:
        return Path("")
    p = Path(raw_str)
    ws_root = Path(workspace_dir) if workspace_dir is not None else None

    base_candidates: List[Path] = []
    if p.is_absolute():
        base_candidates.append(p)
    else:
        if ws_root is not None:
            base_candidates.append(ws_root / p)
        base_candidates.append(p)
        base_candidates.append(PROJECT_ROOT / p)

    for base_p in base_candidates:
        if base_p.is_file():
            return base_p
        for ext in (".webp", ".jpg", ".jpeg", ".png"):
            sib = base_p.with_suffix(ext)
            if sib != base_p and sib.is_file():
                return sib

    default_target = base_candidates[0] if base_candidates else p
    if default_target.suffix.lower() in {".jpg", ".jpeg", ".png", ".webp"}:
        return default_target.with_suffix(".webp")
    return default_target


def get_image_resolution_score(image_path_str: Optional[str], workspace_dir: Optional[Path] = None) -> float:
    """Computes a numeric image resolution score (megapixels or pixel count, fallback to byte size)."""
    if not image_path_str or not isinstance(image_path_str, str) or not image_path_str.strip():
        return 0.0
    resolved_path = resolve_voucher_image_path(image_path_str, workspace_dir=workspace_dir)
    if resolved_path.is_file():
        try:
            with Image.open(resolved_path) as img:
                w, h = img.size
                return float(w * h)
        except Exception:
            try:
                return float(resolved_path.stat().st_size)
            except Exception:
                return 0.0
    return 0.0


def _resolve_row_standardized_packera(row: Mapping[str, Any]) -> str:
    """Resolve standardized canonical 'Packera <epithet>' binomial from a voucher row."""
    from scripts.core.harvester import _extract_canonical_binomial, standardize_packera_taxon
    for col in ("species_standardized", "scientificName", "species_raw"):
        if col not in row:
            continue
        raw_val = row.get(col)
        if raw_val is None or (isinstance(raw_val, float) and math.isnan(raw_val)) or pd.isna(raw_val):
            continue
        s_txt = str(raw_val).strip()
        if not s_txt:
            continue
        std = standardize_packera_taxon(s_txt)
        if std:
            return std
        if s_txt.startswith("Packera "):
            fb = _extract_canonical_binomial(s_txt, {"Packera"})
            if fb:
                return fb
    return ""


def _is_extralimital_status(row_or_status: Any) -> bool:
    """Return True if a voucher row carries an out-of-range or exotic/waif USDA/NatureServe flag."""
    if isinstance(row_or_status, dict) or (pd is not None and isinstance(row_or_status, pd.Series)):
        raw_st = row_or_status.get("usda_county_status")
    else:
        raw_st = row_or_status
    if raw_st is None or (pd is not None and pd.isna(raw_st)):
        return False
    return str(raw_st).strip().lower() in {"out_of_range", "natureserve_exotic_or_waif"}


def _as_bool_flag(val: Any) -> bool:
    if val is None or (pd is not None and pd.isna(val)):
        return False
    if isinstance(val, bool):
        return val
    return str(val).strip().lower() in {"true", "1", "t", "yes"}


def stratify_duplicates(
    df: pd.DataFrame,
    max_distance_meters: float = 50.0,
    workspace_dir: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
    home_herbarium: str = DEFAULT_HOME_HERBARIUM,
) -> pd.DataFrame:
    """
    Identifies duplicate voucher specimens (exsiccatae) representing the same
    field collection event distributed across multiple herbaria.

    Duplicates are resolved by:
      1. Standardized exsiccatae fingerprint (norm_collector + norm_number + norm_date).
      2. If collection number is missing, falls back to geographic coordinate proximity (< 50 m)
         within identical collection dates.

    Designates exactly one primary voucher per duplicate group using authority-first hierarchy:
      Priority 1: Verified Type Status (and ICN Art. 9.6 sub-rank)
      Priority 2: Documented in-range standing
      Priority 3: Determiner tier (`Tier_1_Gold` > `Tier_2_Silver` > `Tier_3_Bronze`)
      Priority 4: Home Herbarium (default: NCU; configurable via home_herbarium)
      Priority 5: Highest specimen image resolution

    Adds/updates columns:
      - exsiccatae_key: Standardized signature string
      - is_primary_duplicate: Boolean flag (True for primary specimen, False for duplicates)
      - duplicate_count: Total sheets in this exsiccatae group
      - exsiccatae_taxonomic_conflict: Boolean flag (True when duplicates disagree in Packera taxon)
      - exsiccatae_gold_taxon: Standardized Packera determination held by Tier_1_Gold duplicate(s)
    """
    if df.empty:
        df["exsiccatae_key"] = pd.Series(dtype=str)
        df["is_primary_duplicate"] = pd.Series(dtype=bool)
        df["duplicate_count"] = pd.Series(dtype=int)
        df["exsiccatae_taxonomic_conflict"] = pd.Series(dtype=bool)
        df["exsiccatae_gold_taxon"] = pd.Series(dtype=str)
        return df

    orig_index = df.index
    df = df.reset_index(drop=True).copy()

    # Ensure metadata columns exist
    if "recordedBy" not in df.columns:
        df["recordedBy"] = df.get("identifiedBy", "")
    if "recordNumber" not in df.columns:
        df["recordNumber"] = ""
    if "image_path" not in df.columns:
        df["image_path"] = ""

    norm_collectors: List[str] = []
    norm_numbers: List[str] = []
    norm_dates: List[str] = []
    has_coords: List[bool] = []
    lats: List[Optional[float]] = []
    lons: List[Optional[float]] = []

    for _, row in df.iterrows():
        rec_by = str(row.get("recordedBy", "") if not pd.isna(row.get("recordedBy", "")) else "")
        rec_num = row.get("recordNumber", "") if not pd.isna(row.get("recordNumber", "")) else ""
        ev_date = str(row.get("eventDate", "") if not pd.isna(row.get("eventDate", "")) else "")
        yr = row.get("year")
        mo = row.get("month")
        dy = row.get("day")

        c_norm = normalize_collector(rec_by)
        n_norm = normalize_record_number(rec_num)
        d_norm = normalize_event_date(ev_date, year=yr, month=mo, day=dy)

        lat_val = None
        lon_val = None
        for col_lat in ["decimalLatitude"]:
            if col_lat in row and not pd.isna(row[col_lat]):
                try:
                    lat_val = float(row[col_lat])
                    break
                except (ValueError, TypeError):
                    pass
        for col_lon in ["decimalLongitude"]:
            if col_lon in row and not pd.isna(row[col_lon]):
                try:
                    lon_val = float(row[col_lon])
                    break
                except (ValueError, TypeError):
                    pass

        norm_collectors.append(c_norm)
        norm_numbers.append(n_norm)
        norm_dates.append(d_norm)
        has_coords.append(lat_val is not None and lon_val is not None)
        lats.append(lat_val)
        lons.append(lon_val)

    n_rows = len(df)
    group_labels: List[Optional[str]] = [None] * n_rows
    row_std_taxa: List[str] = [_resolve_row_standardized_packera(df.iloc[i]) for i in range(n_rows)]

    # Pass 1: Assign keys for records with explicit collection numbers
    for idx in range(n_rows):
        if norm_numbers[idx]:
            group_labels[idx] = f"{norm_collectors[idx]}_{norm_numbers[idx]}_{norm_dates[idx]}"

    # Pass 2: Missing collection number -> fallback to coordinate proximity (< 50m) on identical dates
    missing_num_indices = [idx for idx in range(n_rows) if not norm_numbers[idx]]
    date_to_indices: Dict[str, List[int]] = {}
    for idx in missing_num_indices:
        d = norm_dates[idx]
        date_to_indices.setdefault(d, []).append(idx)

    for d_val, idx_list in date_to_indices.items():
        if d_val == "unknown_date" or d_val.endswith("-00-00"):
            for i in idx_list:
                cat = sanitize_filename(str(df.iloc[i].get("catalogNumber", i)))
                group_labels[i] = f"{norm_collectors[i]}_sn_{cat}_{d_val}"
            continue

        coord_indices: List[int] = []
        no_coord_indices: List[int] = []
        for i in idx_list:
            is_cent = _as_bool_flag(df.iloc[i].get("is_centroid_or_inst_coord"))
            if not has_coords[i] or norm_collectors[i] == "unknown" or is_cent:
                no_coord_indices.append(i)
            else:
                coord_indices.append(i)

        # Records lacking exact coords, named collector, or centroid hygiene stay unique singletons
        for i in no_coord_indices:
            cat = sanitize_filename(str(df.iloc[i].get("catalogNumber", i)))
            group_labels[i] = f"{norm_collectors[i]}_sn_{cat}_{d_val}"

        # Match unnumbered records to any numbered collection on the same date within max_distance_meters
        numbered_date_indices = [
            i for i in range(n_rows)
            if norm_numbers[i]
            and norm_dates[i] == d_val
            and has_coords[i]
            and norm_collectors[i] != "unknown"
            and not _as_bool_flag(df.iloc[i].get("is_centroid_or_inst_coord"))
        ]

        matched_to_numbered = set()
        for i in coord_indices:
            for num_i in numbered_date_indices:
                same_tx = (
                    not row_std_taxa[i]
                    or not row_std_taxa[num_i]
                    or row_std_taxa[i] == row_std_taxa[num_i]
                )
                if same_tx and norm_collectors[i] == norm_collectors[num_i]:
                    dist = haversine_distance_meters(lats[i], lons[i], lats[num_i], lons[num_i])
                    if dist <= max_distance_meters:
                        group_labels[i] = group_labels[num_i]
                        matched_to_numbered.add(i)
                        break

        remaining_coord_indices = [i for i in coord_indices if i not in matched_to_numbered]
        num_rem = len(remaining_coord_indices)
        if num_rem == 0:
            continue

        visited = [False] * num_rem
        for a in range(num_rem):
            if visited[a]:
                continue
            cluster = [a]
            visited[a] = True
            q = [a]
            while q:
                curr = q.pop(0)
                curr_idx = remaining_coord_indices[curr]
                for b in range(num_rem):
                    if not visited[b]:
                        b_idx = remaining_coord_indices[b]
                        same_coll = (
                            norm_collectors[curr_idx] == norm_collectors[b_idx]
                            and norm_collectors[curr_idx] != "unknown"
                        )
                        same_tx = (
                            not row_std_taxa[curr_idx]
                            or not row_std_taxa[b_idx]
                            or row_std_taxa[curr_idx] == row_std_taxa[b_idx]
                        )
                        if same_coll and same_tx:
                            dist = haversine_distance_meters(
                                lats[curr_idx], lons[curr_idx], lats[b_idx], lons[b_idx]
                            )
                            if dist <= max_distance_meters:
                                visited[b] = True
                                cluster.append(b)
                                q.append(b)

            actual_indices = [remaining_coord_indices[c_idx] for c_idx in cluster]
            best_coll = "unknown"
            for act_i in actual_indices:
                if norm_collectors[act_i] != "unknown":
                    best_coll = norm_collectors[act_i]
                    break

            if len(actual_indices) > 1:
                anchor_cat = sanitize_filename(str(df.iloc[actual_indices[0]].get("catalogNumber", actual_indices[0])))
                cluster_key = f"{best_coll}_geo_{anchor_cat}_{d_val}"
            else:
                act_i = actual_indices[0]
                cat = sanitize_filename(str(df.iloc[act_i].get("catalogNumber", act_i)))
                cluster_key = f"{best_coll}_sn_{cat}_{d_val}"

            for act_i in actual_indices:
                group_labels[act_i] = cluster_key

    # Pass 2.5: Merge same-institution portal barcode aliases (e.g., 19398 vs WILLI_19398)
    # when stripped catalog suffixes match and collection metadata is compatible.
    alias_groups: Dict[str, List[int]] = {}
    row_cat_prefixes: List[str] = [""] * n_rows
    row_inst_codes: List[str] = [""] * n_rows
    for idx in range(n_rows):
        row = df.iloc[idx]
        raw_inst = row.get("institutionCode", "")
        inst_clean = "" if (raw_inst is None or pd.isna(raw_inst)) else str(raw_inst).strip().upper()
        row_inst_codes[idx] = inst_clean
        raw_cat = row.get("catalogNumber", "")
        cat_clean = "" if (raw_cat is None or pd.isna(raw_cat)) else str(raw_cat).strip()
        m_pref = re.match(r"^([A-Za-z]{2,12})[-_ ]+(.+)$", cat_clean)
        if m_pref:
            row_cat_prefixes[idx] = m_pref.group(1).strip().upper()
            stripped_cat = m_pref.group(2).strip().upper()
        else:
            stripped_cat = cat_clean.upper()
        if len(stripped_cat) >= 3:
            alias_groups.setdefault(stripped_cat, []).append(idx)

    for _stripped_cat, cand_idxs in alias_groups.items():
        if len(cand_idxs) < 2:
            continue
        for pos_a in range(len(cand_idxs)):
            i = cand_idxs[pos_a]
            for pos_b in range(pos_a + 1, len(cand_idxs)):
                j = cand_idxs[pos_b]
                if group_labels[i] == group_labels[j]:
                    continue
                inst_i = row_inst_codes[i]
                inst_j = row_inst_codes[j]
                pref_i = row_cat_prefixes[i]
                pref_j = row_cat_prefixes[j]
                same_portal = (
                    (inst_i != "" and inst_i == inst_j)
                    or (pref_i != "" and pref_i == inst_j)
                    or (pref_j != "" and pref_j == inst_i)
                    or (pref_i != "" and pref_i == pref_j)
                )
                if not same_portal:
                    continue
                coll_ok = (
                    norm_collectors[i] == norm_collectors[j]
                    or norm_collectors[i] == "unknown"
                    or norm_collectors[j] == "unknown"
                )
                num_ok = (
                    norm_numbers[i] == norm_numbers[j]
                    or not norm_numbers[i]
                    or not norm_numbers[j]
                )
                date_ok = (
                    norm_dates[i] == norm_dates[j]
                    or norm_dates[i] == "unknown_date"
                    or norm_dates[j] == "unknown_date"
                )
                coords_ok = True
                if has_coords[i] and has_coords[j]:
                    dist_ij = haversine_distance_meters(lats[i], lons[i], lats[j], lons[j])
                    coords_ok = dist_ij <= max(float(max_distance_meters), 5000.0)
                if coll_ok and num_ok and date_ok and coords_ok:
                    lbl_i = group_labels[i] or ""
                    lbl_j = group_labels[j] or ""
                    if ("_sn_" in lbl_i and "_sn_" not in lbl_j) or (
                        ("_sn_" in lbl_i) == ("_sn_" in lbl_j) and len(lbl_j) < len(lbl_i)
                    ):
                        canon_lbl, old_lbl = lbl_j, lbl_i
                    else:
                        canon_lbl, old_lbl = lbl_i, lbl_j
                    if canon_lbl and old_lbl and canon_lbl != old_lbl:
                        for k_idx in range(n_rows):
                            if group_labels[k_idx] == old_lbl:
                                group_labels[k_idx] = canon_lbl

    for idx in range(n_rows):
        if group_labels[idx] is None:
            cat = sanitize_filename(str(df.iloc[idx].get("catalogNumber", idx)))
            group_labels[idx] = f"{norm_collectors[idx]}_sn_{cat}_{norm_dates[idx]}"

    df["exsiccatae_key"] = group_labels

    # Pass 3: Cross-duplicate taxonomic conflict auditing and authority-first duplicate ranking
    grouped = df.groupby("exsiccatae_key", sort=False)
    df["duplicate_count"] = grouped["catalogNumber"].transform("count")

    idx_to_pos: Dict[Any, int] = {idx_lbl: pos for pos, idx_lbl in enumerate(df.index)}
    is_primary_list = [False] * n_rows
    pre_conf_col = df["exsiccatae_taxonomic_conflict"] if "exsiccatae_taxonomic_conflict" in df.columns else [False] * n_rows
    pre_gold_col = df["exsiccatae_gold_taxon"] if "exsiccatae_gold_taxon" in df.columns else [""] * n_rows
    conflict_list = [_as_bool_flag(v) for v in pre_conf_col]
    gold_taxon_list = [
        "" if (v is None or (pd is not None and pd.isna(v))) else str(v).strip()
        for v in pre_gold_col
    ]

    for key, group_df in grouped:
        indices = list(group_df.index)
        if len(indices) == 1:
            is_primary_list[idx_to_pos[indices[0]]] = True
            continue

        scored_candidates = []
        idx_to_std_tx: Dict[Any, str] = {}
        idx_to_tier: Dict[Any, str] = {}
        idx_to_score: Dict[Any, Tuple[Any, ...]] = {}
        idx_to_extra: Dict[Any, bool] = {}

        for idx in indices:
            row = group_df.loc[idx]

            # Priority 1: Verified Nomenclatural Type Rank (2 = Primary Type, 1 = Secondary Type, 0 = Non-type)
            # and ICN Art. 9.6 fine-grained float sub-rank (2.2 = Holotype-tier, 2.0 = Isotype-tier,
            # 1.5 = Syntype-tier, 1.0 = Paratype-tier, 0.0 = Non-type)
            raw_ts = row.get("type_status")
            if raw_ts is None or (pd is not None and pd.isna(raw_ts)) or not str(raw_ts).strip():
                raw_ts = row.get("typeStatus")
            type_rank, _ = classify_nomenclatural_type_rank(raw_ts)
            type_subrank = classify_nomenclatural_type_subrank(raw_ts)
            has_type = int(type_rank)

            # Priority 1.5: Documented in-range status ahead of extralimital sheets
            is_extra_i = _is_extralimital_status(row)
            in_range_rank = 0 if is_extra_i else 1

            # Priority 2: Home Herbarium (default: NCU, configurable via home_herbarium)
            target_home = str(home_herbarium).strip().upper() if home_herbarium else ""
            inst = str(row.get("institutionCode", "") if not pd.isna(row.get("institutionCode", "")) else "").strip().upper()
            is_home_herb = 1 if (target_home and inst == target_home) else 0

            # Priority 3: Tier 1 Specialist determinations
            tier = str(row.get("determiner_tier", "") if not pd.isna(row.get("determiner_tier", "")) else "").strip()
            tier_score = 3 if tier == "Tier_1_Gold" else (2 if tier == "Tier_2_Silver" else 1)

            # Priority 4: Highest specimen image resolution
            img_path = str(row.get("image_path", "") if not pd.isna(row.get("image_path", "")) else "")
            resolved_img = resolve_voucher_image_path(img_path, workspace_dir=workspace_dir) if img_path else Path("")
            res_score = get_image_resolution_score(str(resolved_img) if img_path else "", workspace_dir=workspace_dir)

            cat_tiebreak = str(row.get("catalogNumber", idx))

            score_tuple = (has_type, type_subrank, tier_score, in_range_rank, is_home_herb, res_score, cat_tiebreak)
            scored_candidates.append((score_tuple, idx))
            idx_to_std_tx[idx] = _resolve_row_standardized_packera(row)
            idx_to_tier[idx] = tier
            idx_to_score[idx] = score_tuple
            idx_to_extra[idx] = is_extra_i

        distinct_taxa = {tx for tx in idx_to_std_tx.values() if tx}
        inrange_taxa = {idx_to_std_tx[idx] for idx in indices if idx_to_std_tx.get(idx) and not idx_to_extra.get(idx, False)}
        has_conflict = len(distinct_taxa) >= 2

        inrange_gold_counts: Dict[str, int] = {}
        inrange_gold_best_idx: Dict[str, Any] = {}
        all_gold_counts: Dict[str, int] = {}
        all_gold_best_idx: Dict[str, Any] = {}
        for idx in indices:
            tx_i = idx_to_std_tx.get(idx, "")
            if idx_to_tier.get(idx) == "Tier_1_Gold" and tx_i:
                all_gold_counts[tx_i] = all_gold_counts.get(tx_i, 0) + 1
                if tx_i not in all_gold_best_idx or idx_to_score[idx] > idx_to_score[all_gold_best_idx[tx_i]]:
                    all_gold_best_idx[tx_i] = idx
                if not idx_to_extra.get(idx, False) and not conflict_list[idx_to_pos[idx]]:
                    inrange_gold_counts[tx_i] = inrange_gold_counts.get(tx_i, 0) + 1
                    if tx_i not in inrange_gold_best_idx or idx_to_score[idx] > idx_to_score[inrange_gold_best_idx[tx_i]]:
                        inrange_gold_best_idx[tx_i] = idx

        best_gold_taxon = ""
        best_gold_source_idx: Any = None
        best_gold_from_inrange = False
        active_gold_counts = inrange_gold_counts if inrange_gold_counts else all_gold_counts
        active_gold_best = inrange_gold_best_idx if inrange_gold_counts else all_gold_best_idx
        if active_gold_counts:
            best_gold_from_inrange = bool(inrange_gold_counts)
            ranked_gold = sorted(
                active_gold_counts.keys(),
                key=lambda tx: (active_gold_counts[tx], idx_to_score[active_gold_best[tx]]),
                reverse=True,
            )
            best_gold_taxon = ranked_gold[0]
            best_gold_source_idx = active_gold_best[best_gold_taxon]

        multi_gold_disagreement = len(active_gold_counts) >= 2
        unanimous_inrange_tx = next(iter(inrange_taxa)) if len(inrange_taxa) == 1 else ""
        for idx in indices:
            pos_i = idx_to_pos[idx]
            row_tx_i = idx_to_std_tx.get(idx, "")
            if not has_conflict:
                row_conflict = False
            elif best_gold_taxon and not multi_gold_disagreement:
                # Unanimous Tier_1_Gold consensus exists: flag only sheets whose label disagreed with Gold consensus
                row_conflict = bool(row_tx_i and row_tx_i != best_gold_taxon)
            elif unanimous_inrange_tx:
                # Unanimous in-range determination across duplicates: flag only disagreeing extralimital duplicates
                row_conflict = bool(idx_to_extra.get(idx, False) and row_tx_i and row_tx_i != unanimous_inrange_tx)
            else:
                row_conflict = True
            conflict_list[pos_i] = bool(conflict_list[pos_i] or row_conflict)
            if best_gold_taxon and (not multi_gold_disagreement or best_gold_from_inrange or idx_to_extra.get(idx, False)):
                gold_taxon_list[pos_i] = best_gold_taxon

        scored_candidates.sort(key=lambda x: x[0], reverse=True)
        primary_idx = scored_candidates[0][1]
        is_primary_list[idx_to_pos[primary_idx]] = True

    df["is_primary_duplicate"] = is_primary_list
    df["exsiccatae_taxonomic_conflict"] = conflict_list
    df["exsiccatae_gold_taxon"] = gold_taxon_list
    df.index = orig_index

    if logger:
        n_prim = sum(is_primary_list)
        n_dup = n_rows - n_prim
        n_conf = sum(1 for v in conflict_list if v)
        logger.info(
            f"Duplicate Stratification: Processed {n_rows} vouchers -> "
            f"{n_prim} primary duplicate vouchers designated ({n_dup} duplicate sheets flagged, "
            f"{n_conf} sheets in taxonomically conflicted duplicate groups)."
        )

    return df


def reconcile_cross_taxon_exsiccatae_with_db(
    df: pd.DataFrame,
    conn: Any,
    standardize_fn: Optional[Any] = None,
) -> pd.DataFrame:
    """
    Cross-reference `exsiccatae_key` across the full SQLite `vouchers` table so duplicate
    sheets filed under separate online species still detect taxonomic disagreement and
    in-range specialist annotations.
    """
    if pd is None or df is None or df.empty or conn is None or "exsiccatae_key" not in df.columns:
        return df

    try:
        has_tbl = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='vouchers' LIMIT 1"
        ).fetchone()
    except Exception:
        return df
    if not has_tbl:
        return df

    raw_keys = [
        str(v).strip()
        for v in df["exsiccatae_key"].dropna().unique().tolist()
        if str(v).strip() and str(v).strip().lower() not in {"none", "nan"}
    ]
    if not raw_keys:
        return df

    try:
        v_cols = {r[1] for r in conn.execute("PRAGMA table_info(vouchers)").fetchall()}
    except Exception:
        return df

    want_cols = [
        c
        for c in (
            "catalogNumber",
            "exsiccatae_key",
            "species_standardized",
            "scientificName",
            "species_raw",
            "determiner_tier",
            "usda_county_status",
        )
        if c in v_cols
    ]
    if "exsiccatae_key" not in want_cols:
        return df

    db_rows: List[Dict[str, Any]] = []
    sel_sql = ", ".join(want_cols)
    for start in range(0, len(raw_keys), 500):
        chunk = raw_keys[start : start + 500]
        ph = ", ".join(["?"] * len(chunk))
        try:
            cur = conn.execute(
                f"SELECT {sel_sql} FROM vouchers WHERE exsiccatae_key IN ({ph})",
                chunk,
            )
            for row_t in cur.fetchall():
                db_rows.append(dict(zip(want_cols, row_t)))
        except Exception:
            continue

    if not db_rows:
        return df

    def _std_tx(rec: Any) -> str:
        raw_val = ""
        for fld in ("species_standardized", "scientificName", "species_raw"):
            v = rec.get(fld) if isinstance(rec, dict) else rec.get(fld, None)
            if v is not None and not (pd is not None and pd.isna(v)) and str(v).strip():
                raw_val = str(v).strip()
                break
        if not raw_val:
            return ""
        if standardize_fn is not None:
            try:
                res = standardize_fn(raw_val)
                if res:
                    return str(res).strip()
            except Exception:
                pass
        return _resolve_row_standardized_packera(rec if isinstance(rec, dict) else rec)

    # Index DB records grouped by exsiccatae_key
    db_by_key: Dict[str, List[Dict[str, Any]]] = {}
    for drec in db_rows:
        ek = str(drec.get("exsiccatae_key") or "").strip()
        if ek:
            db_by_key.setdefault(ek, []).append(drec)

    out_df = df.copy()
    if "exsiccatae_taxonomic_conflict" not in out_df.columns:
        out_df["exsiccatae_taxonomic_conflict"] = False
    if "exsiccatae_gold_taxon" not in out_df.columns:
        out_df["exsiccatae_gold_taxon"] = ""

    conf_vals = [_as_bool_flag(v) for v in out_df["exsiccatae_taxonomic_conflict"].tolist()]
    gold_vals = [
        "" if (v is None or pd.isna(v)) else str(v).strip()
        for v in out_df["exsiccatae_gold_taxon"].tolist()
    ]

    for key, idx_lbls in out_df.groupby("exsiccatae_key", sort=False).groups.items():
        ek = str(key or "").strip()
        if not ek or ek not in db_by_key:
            continue
        batch_indices = list(idx_lbls)
        batch_cats = {
            str(out_df.at[i_l, "catalogNumber"]).strip()
            for i_l in batch_indices
            if "catalogNumber" in out_df.columns and not pd.isna(out_df.at[i_l, "catalogNumber"])
        }
        siblings: List[Dict[str, Any]] = []
        for i_l in batch_indices:
            r_s = out_df.loc[i_l]
            siblings.append({
                "catalogNumber": str(r_s.get("catalogNumber") or ""),
                "tx": _std_tx(r_s),
                "tier": str(r_s.get("determiner_tier") or "").strip(),
                "extra": _is_extralimital_status(r_s),
                "pre_conf": _as_bool_flag(r_s.get("exsiccatae_taxonomic_conflict")),
            })
        for drec in db_by_key[ek]:
            d_cat = str(drec.get("catalogNumber") or "").strip()
            d_tx = _std_tx(drec)
            # Include DB siblings under other catalogNumbers (or differing taxon snapshots)
            if d_cat and d_cat in batch_cats:
                same_cat_tx = [s["tx"] for s in siblings if s["catalogNumber"].strip() == d_cat]
                if any(t == d_tx for t in same_cat_tx):
                    continue
            siblings.append({
                "catalogNumber": d_cat,
                "tx": d_tx,
                "tier": str(drec.get("determiner_tier") or "").strip(),
                "extra": _is_extralimital_status(drec),
                "pre_conf": False,
            })

        distinct_taxa = {s["tx"] for s in siblings if s["tx"]}
        if len(distinct_taxa) < 2:
            continue

        inrange_taxa = {s["tx"] for s in siblings if s["tx"] and not s["extra"]}
        inrange_gold_counts: Dict[str, int] = {}
        all_gold_counts: Dict[str, int] = {}
        for s in siblings:
            if s["tier"] == "Tier_1_Gold" and s["tx"]:
                all_gold_counts[s["tx"]] = all_gold_counts.get(s["tx"], 0) + 1
                if not s["extra"] and not s["pre_conf"]:
                    inrange_gold_counts[s["tx"]] = inrange_gold_counts.get(s["tx"], 0) + 1

        active_gold = inrange_gold_counts if inrange_gold_counts else all_gold_counts
        best_gold_tx = ""
        best_gold_from_inrange = False
        multi_gold = len(active_gold) >= 2
        if active_gold:
            best_gold_from_inrange = bool(inrange_gold_counts)
            best_gold_tx = sorted(active_gold.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)[0][0]
        unanimous_inrange_tx = next(iter(inrange_taxa)) if len(inrange_taxa) == 1 else ""

        pos_map = {idx_lbl: pos for pos, idx_lbl in enumerate(out_df.index)}
        for i_l in batch_indices:
            pos_i = pos_map[i_l]
            r_s = out_df.loc[i_l]
            r_tx = _std_tx(r_s)
            r_extra = _is_extralimital_status(r_s)
            if best_gold_tx and not multi_gold:
                row_conf = bool(r_tx and r_tx != best_gold_tx)
            elif unanimous_inrange_tx:
                row_conf = bool(r_extra and r_tx and r_tx != unanimous_inrange_tx)
            else:
                row_conf = True
            conf_vals[pos_i] = bool(conf_vals[pos_i] or row_conf)
            if best_gold_tx and not gold_vals[pos_i] and (not multi_gold or best_gold_from_inrange or r_extra):
                gold_vals[pos_i] = best_gold_tx

    out_df["exsiccatae_taxonomic_conflict"] = conf_vals
    out_df["exsiccatae_gold_taxon"] = gold_vals
    return out_df


def spatially_balance_and_thin_vouchers(
    df: pd.DataFrame,
    spatial_thin_km: float = 0.0,
    max_per_state: Optional[int] = None,
    max_per_region: Optional[int] = None,
    primary_duplicates_only: bool = True,
    preserve_all_gold: bool = True,
    manifest_df: Optional[pd.DataFrame] = None,
    workspace_dir: Optional[Path] = None,
    preserve_catalog_numbers: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """
    Stratifies, deduplicates, spatially thins, and computes continuous inverse
    spatial-density weights (`spatial_balance_weight`) for botanical vouchers.
    """
    if df is None:
        return pd.DataFrame(columns=EXPORT_COLUMNS + ["spatial_balance_weight"])
    if df.empty:
        out_empty = df.copy()
        out_empty["spatial_balance_weight"] = pd.Series(dtype=float)
        return out_empty

    from scripts.core.harvester import standardize_packera_taxon

    out_df = df.copy()

    # Ensure standardized taxon column exists
    if "species_standardized" not in out_df.columns:
        src_col = "scientificName" if "scientificName" in out_df.columns else (
            "species_raw" if "species_raw" in out_df.columns else None
        )
        if src_col is not None:
            out_df["species_standardized"] = out_df[src_col].map(
                lambda x: standardize_packera_taxon(str(x)) or str(x).strip() or "Unknown"
                if not pd.isna(x) else "Unknown"
            )
        else:
            out_df["species_standardized"] = "Unknown"
    else:
        src_col = "scientificName" if "scientificName" in out_df.columns else (
            "species_raw" if "species_raw" in out_df.columns else None
        )
        def _fill_sp(row: pd.Series) -> str:
            val = row.get("species_standardized")
            if val is not None and not pd.isna(val) and str(val).strip():
                return str(val).strip()
            if src_col is not None:
                raw_v = row.get(src_col)
                if raw_v is not None and not pd.isna(raw_v) and str(raw_v).strip():
                    return standardize_packera_taxon(str(raw_v)) or str(raw_v).strip()
            return "Unknown"
        out_df["species_standardized"] = out_df.apply(_fill_sp, axis=1)

    preserve_cats_set: Set[str] = {
        str(c).strip() for c in (preserve_catalog_numbers or [])
        if c is not None and str(c).strip()
    }

    # Step A: Exsiccatae Deduplication (unconditionally retaining all nomenclatural type specimens and preserved vouchers)
    if primary_duplicates_only:
        if "is_primary_duplicate" not in out_df.columns:
            if any(c in out_df.columns for c in ("recordedBy", "recordNumber", "exsiccatae_key")):
                out_df = stratify_duplicates(out_df, workspace_dir=workspace_dir)
            else:
                out_df["is_primary_duplicate"] = True
        prim_mask = out_df["is_primary_duplicate"].map(
            lambda v: v is True or str(v).strip().lower() in {"true", "1", "t", "yes"}
        )
        pres_mask_a = (
            out_df["catalogNumber"].astype(str).str.strip().isin(preserve_cats_set)
            if (preserve_cats_set and "catalogNumber" in out_df.columns)
            else pd.Series(False, index=out_df.index)
        )
        ts_col_a = "type_status" if "type_status" in out_df.columns else (
            "typeStatus" if "typeStatus" in out_df.columns else None
        )
        if ts_col_a is not None:
            type_mask_a = out_df[ts_col_a].map(lambda t: classify_nomenclatural_type_rank(t)[0] >= 1)
            keep_mask = prim_mask | type_mask_a | pres_mask_a
        else:
            keep_mask = prim_mask | pres_mask_a
        out_df = out_df.loc[keep_mask].reset_index(drop=True)
        if out_df.empty:
            out_df["spatial_balance_weight"] = pd.Series(dtype=float)
            return out_df

    # Step B: Candidate Priority Ranking
    t1_leaf_counts: Dict[str, int] = {}
    if manifest_df is not None and not manifest_df.empty and "catalogNumber" in manifest_df.columns:
        for cat_k, cnt in manifest_df["catalogNumber"].astype(str).value_counts().items():
            t1_leaf_counts[str(cat_k)] = int(cnt)

    lat_col = "decimalLatitude" if "decimalLatitude" in out_df.columns else None
    lon_col = "decimalLongitude" if "decimalLongitude" in out_df.columns else None

    n_rows = len(out_df)
    lats: List[Optional[float]] = [None] * n_rows
    lons: List[Optional[float]] = [None] * n_rows
    priority_tuples: List[Tuple[int, float, int, int, float, float, str]] = []
    is_gold_flags: List[bool] = []
    is_type_flags: List[bool] = []
    is_preserved_flags: List[bool] = []

    for i in range(n_rows):
        row = out_df.iloc[i]
        cat_str = str(row.get("catalogNumber", i))
        is_preserved_flags.append(cat_str.strip() in preserve_cats_set)

        # 1. Nomenclatural type rank (2 = Primary Type, 1 = Secondary Type, 0 = Non-type) + ICN Art. 9.6 sub-rank
        raw_ts = row.get("type_status")
        if raw_ts is None or (pd is not None and pd.isna(raw_ts)) or not str(raw_ts).strip():
            raw_ts = row.get("typeStatus")
        has_type, _ = classify_nomenclatural_type_rank(raw_ts)
        t_subrank = classify_nomenclatural_type_subrank(raw_ts)
        is_type_flags.append(has_type >= 1)

        # 2. determiner_tier
        tier_str = str(row.get("determiner_tier", "") if not pd.isna(row.get("determiner_tier", "")) else "").strip()
        tier_score = 3 if (tier_str == "Tier_1_Gold" or has_type >= 1) else (2 if tier_str == "Tier_2_Silver" else 1)
        is_gold_flags.append(tier_str == "Tier_1_Gold" or has_type >= 1)

        # 3. Extracted leaf count
        if cat_str in t1_leaf_counts:
            n_t1 = t1_leaf_counts[cat_str]
        elif "leaf_count" in out_df.columns and not pd.isna(row.get("leaf_count")):
            try:
                n_t1 = int(float(row.get("leaf_count", 0)))
            except (ValueError, TypeError):
                n_t1 = 0
        else:
            n_t1 = 0

        # Coordinates & 4. coordinate_precision
        lat_v: Optional[float] = None
        lon_v: Optional[float] = None
        if lat_col and lon_col:
            try:
                lv = float(row.get(lat_col))
                lnv = float(row.get(lon_col))
                if math.isfinite(lv) and math.isfinite(lnv):
                    lat_v, lon_v = lv, lnv
            except (ValueError, TypeError):
                pass
        lats[i] = lat_v
        lons[i] = lon_v

        coord_prec = 0.0
        if lat_v is not None and lon_v is not None:
            unc_raw = row.get("coordinateUncertainty")
            try:
                unc_f = float(unc_raw)
                if math.isfinite(unc_f) and unc_f >= 0.0:
                    coord_prec = 1.0 / (1.0 + unc_f)
                else:
                    coord_prec = 1.0 / (1.0 + 1000.0)
            except (ValueError, TypeError):
                coord_prec = 1.0 / (1.0 + 1000.0)

        # 5. image_resolution_score
        img_p = str(row.get("image_path", "") if not pd.isna(row.get("image_path", "")) else "")
        resolved_img_p = resolve_voucher_image_path(img_p, workspace_dir=workspace_dir) if img_p else Path("")
        res_score = get_image_resolution_score(str(resolved_img_p) if img_p else "", workspace_dir=workspace_dir)

        priority_tuples.append((has_type, t_subrank, tier_score, n_t1, round(coord_prec, 8), round(res_score, 4), cat_str))

    active_indices = list(range(n_rows))

    # Step C: Greedy Priority-Ordered Pairwise Haversine Thinning (spatial_thin_km > 0)
    if spatial_thin_km is not None and float(spatial_thin_km) > 0.0:
        thin_m = float(spatial_thin_km) * 1000.0
        sp_geo_groups: Dict[str, List[int]] = {}
        non_geo_indices: List[int] = []

        for idx in active_indices:
            lat_v = lats[idx]
            lon_v = lons[idx]
            sp_val = str(out_df.iloc[idx]["species_standardized"])
            if lat_v is None or lon_v is None:
                non_geo_indices.append(idx)
            else:
                sp_geo_groups.setdefault(sp_val, []).append(idx)

        thinned_indices: List[int] = []
        for _sp_k, sp_idxs in sp_geo_groups.items():
            sp_sorted = sorted(sp_idxs, key=lambda i: priority_tuples[i], reverse=True)
            kept_sp: List[int] = []
            for cand_i in sp_sorted:
                if is_type_flags[cand_i] or is_preserved_flags[cand_i] or (preserve_all_gold and is_gold_flags[cand_i]):
                    kept_sp.append(cand_i)
                    continue
                c_lat = lats[cand_i]
                c_lon = lons[cand_i]
                too_close = any(
                    haversine_distance_meters(c_lat, c_lon, lats[k_i], lons[k_i]) < thin_m
                    for k_i in kept_sp
                )
                if not too_close:
                    kept_sp.append(cand_i)
            thinned_indices.extend(kept_sp)

        thinned_indices.extend(non_geo_indices)
        active_indices = sorted(set(thinned_indices))

    def _select_spatially_spread_quota(
        cand_pool: List[int],
        quota: int,
        stratum_values: Dict[int, str],
    ) -> List[int]:
        """
        Selects up to `quota` indices from `cand_pool`, unconditionally preserving
        all nomenclatural type specimens (and all Gold anchors when `preserve_all_gold=True`)
        and spreading remaining selections across `stratum_values` (counties or states)
        and farthest Haversine distance.
        """
        if len(cand_pool) <= quota:
            return sorted(cand_pool)

        if preserve_all_gold:
            gold_idxs = [i for i in cand_pool if is_gold_flags[i] or is_preserved_flags[i]]
            rem_pool = [i for i in cand_pool if not (is_gold_flags[i] or is_preserved_flags[i])]
        else:
            gold_idxs = [i for i in cand_pool if is_type_flags[i] or is_preserved_flags[i]]
            rem_pool = [i for i in cand_pool if not (is_type_flags[i] or is_preserved_flags[i])]

        rem_pool.sort(key=lambda i: priority_tuples[i], reverse=True)
        selected: List[int] = list(gold_idxs)
        if len(selected) >= quota:
            return (
                sorted(selected)
                if (preserve_all_gold or any(is_type_flags[i] or is_preserved_flags[i] for i in selected))
                else sorted(selected[:quota])
            )

        stratum_counts: Dict[str, int] = {}
        for i in selected:
            st_key = stratum_values[i]
            stratum_counts[st_key] = stratum_counts.get(st_key, 0) + 1

        while len(selected) < quota and rem_pool:
            min_cnt = min(stratum_counts.get(stratum_values[i], 0) for i in rem_pool)
            eligible = [i for i in rem_pool if stratum_counts.get(stratum_values[i], 0) == min_cnt]

            sel_geo = [i for i in selected if lats[i] is not None and lons[i] is not None]
            if not sel_geo:
                best_i = eligible[0]
            else:
                best_i = eligible[0]
                best_score: Tuple[float, Tuple[int, int, int, float, float, str]] = (-1.0, priority_tuples[best_i])
                for cand_i in eligible:
                    if lats[cand_i] is not None and lons[cand_i] is not None:
                        min_d_km = min(
                            haversine_distance_meters(
                                lats[cand_i], lons[cand_i], lats[s_i], lons[s_i]
                            ) / 1000.0
                            for s_i in sel_geo
                        )
                    else:
                        min_d_km = 0.0
                    cand_score = (round(min_d_km, 2), priority_tuples[cand_i])
                    if cand_score > best_score:
                        best_score = cand_score
                        best_i = cand_i

            selected.append(best_i)
            st_key = stratum_values[best_i]
            stratum_counts[st_key] = stratum_counts.get(st_key, 0) + 1
            rem_pool.remove(best_i)

        return sorted(selected)

    # Step D.1: Per-State Quota Balancing (max_per_state > 0)
    if max_per_state is not None and int(max_per_state) > 0:
        state_quota = int(max_per_state)
        state_groups: Dict[Tuple[str, str], List[int]] = {}
        county_labels: Dict[int, str] = {}
        for idx in active_indices:
            row = out_df.iloc[idx]
            sp_val = str(row["species_standardized"])
            st_raw = row.get("stateProvince", "")
            st_val = str(st_raw if not pd.isna(st_raw) else "").strip().upper() or "UNKNOWN_STATE"
            cnty_raw = row.get("county", "")
            cnty_val = str(cnty_raw if not pd.isna(cnty_raw) else "").strip().upper()
            if not cnty_val:
                if lats[idx] is not None and lons[idx] is not None:
                    cnty_val = f"GRID_{int(round(lats[idx] * 4.0))}_{int(round(lons[idx] * 4.0))}"
                else:
                    cnty_val = f"UNSPECIFIED_{idx}"
            county_labels[idx] = cnty_val
            state_groups.setdefault((sp_val, st_val), []).append(idx)

        state_balanced: List[int] = []
        for _grp_key, grp_idxs in state_groups.items():
            state_balanced.extend(_select_spatially_spread_quota(grp_idxs, state_quota, county_labels))
        active_indices = sorted(set(state_balanced))

    # Step D.2: Per-Biogeographic-Region Quota Balancing (max_per_region > 0)
    if max_per_region is not None and int(max_per_region) > 0:
        region_quota = int(max_per_region)
        region_groups: Dict[Tuple[str, str], List[int]] = {}
        state_labels: Dict[int, str] = {}
        for idx in active_indices:
            row = out_df.iloc[idx]
            sp_val = str(row["species_standardized"])
            reg_raw = row.get("regional_group", "")
            reg_val = str(reg_raw if not pd.isna(reg_raw) else "").strip() or "Unassigned_NA"
            st_raw = row.get("stateProvince", "")
            st_val = str(st_raw if not pd.isna(st_raw) else "").strip().upper() or "UNKNOWN_STATE"
            state_labels[idx] = st_val
            region_groups.setdefault((sp_val, reg_val), []).append(idx)

        region_balanced: List[int] = []
        for _grp_key, grp_idxs in region_groups.items():
            region_balanced.extend(_select_spatially_spread_quota(grp_idxs, region_quota, state_labels))
        active_indices = sorted(set(region_balanced))

    out_df = out_df.iloc[active_indices].reset_index(drop=True)
    if out_df.empty:
        out_df["spatial_balance_weight"] = pd.Series(dtype=float)
        return out_df

    # Step E: Compute Continuous `spatial_balance_weight` strictly over non-quarantined active anchors
    ret_lats = np.array(
        [lats[i] if lats[i] is not None else np.nan for i in active_indices],
        dtype=float,
    )
    ret_lons = np.array(
        [lons[i] if lons[i] is not None else np.nan for i in active_indices],
        dtype=float,
    )
    ret_taxa = out_df["species_standardized"].astype(str).tolist()
    has_auth_wt = "authority_weight" in out_df.columns
    eligible_anchor_mask = np.ones(len(out_df), dtype=bool)
    for k in range(len(out_df)):
        row_k = out_df.iloc[k]
        if _as_bool_flag(row_k.get("is_centroid_or_inst_coord")):
            eligible_anchor_mask[k] = False
            continue
        raw_ts_k = row_k.get("type_status") if not pd.isna(row_k.get("type_status", "")) else row_k.get("typeStatus")
        if classify_nomenclatural_type_rank(raw_ts_k)[0] >= 1:
            continue
        if has_auth_wt and not pd.isna(row_k.get("authority_weight")):
            try:
                if float(row_k.get("authority_weight")) <= 0.0:
                    eligible_anchor_mask[k] = False
                    continue
            except (ValueError, TypeError):
                pass
        tier_k = str(row_k.get("determiner_tier", "") if not pd.isna(row_k.get("determiner_tier", "")) else "").strip()
        if tier_k == "Tier_3_Bronze" and not _as_bool_flag(row_k.get("is_promoted_bronze")):
            eligible_anchor_mask[k] = False

    spatial_weights = np.ones(len(out_df), dtype=float)
    bandwidth_km = 25.0

    taxon_to_rows: Dict[str, List[int]] = {}
    for k, sp_name in enumerate(ret_taxa):
        taxon_to_rows.setdefault(sp_name, []).append(k)

    for _sp_name, t_idxs in taxon_to_rows.items():
        t_arr = np.array(t_idxs, dtype=int)
        valid_geo_mask = np.isfinite(ret_lats[t_arr]) & np.isfinite(ret_lons[t_arr]) & eligible_anchor_mask[t_arr]
        geo_idxs = t_arr[valid_geo_mask]
        if len(geo_idxs) == 0:
            continue

        # Continuous geodesic 25 km + 75 km multi-scale Gaussian kernel specimen density d_local
        # (free of arbitrary US state political-area penalties)
        rad_lat = np.radians(ret_lats[geo_idxs])
        rad_lon = np.radians(ret_lons[geo_idxs])
        dlat_mat = rad_lat[:, None] - rad_lat[None, :]
        dlon_mat = rad_lon[:, None] - rad_lon[None, :]
        a_mat = (
            np.sin(dlat_mat / 2.0) ** 2
            + np.cos(rad_lat[:, None]) * np.cos(rad_lat[None, :]) * (np.sin(dlon_mat / 2.0) ** 2)
        )
        dist_km_mat = 2.0 * 6371.0 * np.arcsin(np.clip(np.sqrt(a_mat), 0.0, 1.0))
        kern_fine = np.exp(-(dist_km_mat ** 2) / (2.0 * (bandwidth_km ** 2)))
        kern_regional = 0.35 * np.exp(-(dist_km_mat ** 2) / (2.0 * ((3.0 * bandwidth_km) ** 2)))
        kern_mat = kern_fine + kern_regional
        np.fill_diagonal(kern_mat, 0.0)
        d_local = 1.0 + np.sum(kern_mat, axis=1)

        w_raw = 1.0 / np.sqrt(np.maximum(d_local, 1e-8))
        mean_raw = float(np.mean(w_raw))
        if mean_raw > 0.0:
            w_norm = w_raw / mean_raw
        else:
            w_norm = np.ones_like(w_raw)

        w_spatial = np.clip(w_norm, 0.20, 3.00)
        for _ in range(3):
            m_cur = float(np.mean(w_spatial))
            if m_cur > 0.0:
                w_spatial = np.clip(w_spatial / m_cur, 0.20, 3.00)

        spatial_weights[geo_idxs] = np.round(w_spatial, 6)

    out_df["spatial_balance_weight"] = spatial_weights
    return out_df



def setup_logger(log_file_path: Optional[Path] = None, verbose: bool = False) -> logging.Logger:
    """Configures or retrieves a formatted logger for voucher harvesting."""
    return setup_logging(log_file=log_file_path, verbose=verbose, name="VoucherHarvester")


def sanitize_filename(name: str) -> str:
    """Sanitizes arbitrary strings into safe, valid filesystem filenames across platforms."""
    clean = re.sub(r'[\\/*?:"<>|\s]+', "_", str(name).strip())
    clean = clean.strip("._")
    return clean if clean else "voucher_unknown"


GENERIC_DETERMINER_TOKENS: Set[str] = {
    "",
    "unknown",
    "anonymous",
    "none",
    "nan",
    "null",
    "na",
    "n/a",
    "staff",
    "student",
    "students",
    "herbarium",
    "curator",
    "volunteer",
    "unassigned",
    "unspecified",
    "not",
    "evident",
    "indet",
    "indetermined",
    "undetermined",
    "various",
    "class",
    "admin",
    "database",
    "digitizer",
    "illegible",
    "unrecorded",
}

_INVALID_INSTITUTION_CODES: Set[str] = {
    "",
    "UNKNOWN",
    "UNKNOWN_INST",
    "NONE",
    "NAN",
    "NULL",
    "NA",
}


def _fold_ascii(text: str) -> str:
    """Fold Unicode diacritics to base ASCII lowercase characters (e.g., 'Löve' -> 'love')."""
    norm = unicodedata.normalize("NFKD", text)
    return "".join(ch for ch in norm if not unicodedata.combining(ch)).lower()


_MISSING_IDENTIFIED_BY_TOKENS: Set[str] = {
    "", "none", "nan", "null", "na", "n/a", "<na>",
    "unknown", "unassigned", "unspecified", "unrecorded", "not recorded", "not evident",
}
_NON_PERSON_DETERMINATION_TOKENS: Tuple[str, ...] = (
    "nomenclatural adjustment",
    "nomenclature update",
    "nomenclatural update",
    "taxonomic update",
    "taxonomic adjustment",
    "synonymy update",
    "synonymy",
    "folder update",
    "herbarium folder",
    "database import",
    "batch update",
    "flora of ",
    "fna vol",
    "fna ",
    "usda plants",
    "bonap",
    "checklist",
    "kartesz",
    "[exsiccatae ",
    "sernec",
    "gbif",
    "idigbio",
    "unverified",
    "unknown",
    "anonymous",
    "leg. ign",
    "s.n.",
)
_FLORA_CITATION_REGEX = re.compile(
    r"\b(?:fide|sensu|sec\.|according\s+to|per|after|following)\s+[A-Z]"
    r"|\b(?:FNA|BONAP|USDA)\b"
    r"|^[A-Z][a-z]+(?:\s+(?:et\s+al\.?|&\s*[A-Z][a-z]+))?\s*\(?\b(?:19|20)\d{2}[a-z]?\)?$",
    re.IGNORECASE,
)
_AMBIGUOUS_COMMON_SURNAMES: Set[str] = {
    "anderson", "atwood", "baker", "bell", "britton", "brooks", "brown", "bush", "campbell", "carter",
    "chapman", "clark", "collins", "cook", "cooper", "cox", "davis", "diamond", "dorn",
    "duncan", "estes", "evans", "fisher", "freeman", "fuller", "golden", "gray", "green", "hall", "heller", "hill",
    "horn", "howard", "jackson", "james", "jeffrey", "johnson", "jones", "kelly", "king",
    "lee", "lewis", "long", "love", "martin", "miller", "mitchell", "moore", "morgan",
    "morris", "murphy", "nelson", "packer", "parker", "perry", "phillips", "powell", "price",
    "reed", "richardson", "roberts", "robinson", "rogers", "ross", "russell",
    "scott", "smith", "stewart", "taylor", "thomas", "thompson", "turner", "walker",
    "ward", "watson", "weber", "welsh", "white", "williams", "wilson", "wood", "wright", "young",
}
CORE_PACKERA_MONOGRAPHERS: Set[str] = set(CORE_PACKERA_GOLD_SURNAMES)
REGIONAL_FLORA_SPECIALISTS: Set[str] = set(REGIONAL_AND_HISTORICAL_SILVER_SURNAMES)
PACKERA_SPECIALISTS: Set[str] = CORE_PACKERA_MONOGRAPHERS | REGIONAL_FLORA_SPECIALISTS
_CORE_QUALIFIED_SPECIALIST_REGEXES: Tuple[re.Pattern, ...] = CORE_PACKERA_GOLD_REGEXES
_REGIONAL_QUALIFIED_SPECIALIST_REGEXES: Tuple[re.Pattern, ...] = REGIONAL_AND_HISTORICAL_SILVER_REGEXES
_QUALIFIED_SPECIALIST_REGEXES: Tuple[re.Pattern, ...] = (
    _CORE_QUALIFIED_SPECIALIST_REGEXES + _REGIONAL_QUALIFIED_SPECIALIST_REGEXES
)
_COLLECTOR_REFERENCE_REGEX = re.compile(
    r"^(?:(?:id(?:entified)?|det(?:ermined)?|annot(?:ated)?|ver(?:ified)?|conf(?:irmed)?)\.?\s*(?:by\s+)?)?"
    r"(?:the\s+)?(?:collector(?:s)?|same(?:\s+as\s+collector(?:s)?)?)\.?$",
    re.IGNORECASE,
)


def _has_core_monographer_match(text: str) -> bool:
    """Return True if `text` matches a core Packera / Senecioneae monographer."""
    if not text or not str(text).strip():
        return False
    s_str = str(text).strip()
    if _has_non_person_determination_token(s_str):
        return False
    surnames = set(normalize_determiner_surnames(s_str))
    if surnames & CORE_PACKERA_MONOGRAPHERS:
        return True
    return any(q_rx.search(s_str) for q_rx in _CORE_QUALIFIED_SPECIALIST_REGEXES)


def _has_regional_or_historical_specialist_match(text: str) -> bool:
    """Return True if `text` matches a regional flora or historical Senecioneae specialist."""
    if not text or not str(text).strip():
        return False
    s_str = str(text).strip()
    if _has_non_person_determination_token(s_str):
        return False
    surnames = set(normalize_determiner_surnames(s_str))
    if surnames & REGIONAL_FLORA_SPECIALISTS:
        return True
    return any(q_rx.search(s_str) for q_rx in _REGIONAL_QUALIFIED_SPECIALIST_REGEXES)


def _has_specialist_match(text: str, extra_patterns: Optional[Iterable[str]] = None) -> bool:
    """
    Returns True if `text` matches a recognized North American Packera / Senecioneae
    taxonomic specialist via normalized surnames (`PACKERA_SPECIALISTS`), initial-qualified
    specialist regexes (`_QUALIFIED_SPECIALIST_REGEXES`), `SPECIALIST_PATTERNS`, or optional
    `extra_patterns`.
    """
    if not text or not str(text).strip():
        return False
    s_str = str(text).strip()
    if _has_non_person_determination_token(s_str):
        return False
    surnames = set(normalize_determiner_surnames(s_str))
    if surnames & PACKERA_SPECIALISTS:
        return True
    for q_rx in _QUALIFIED_SPECIALIST_REGEXES:
        if q_rx.search(s_str):
            return True
    for pat in SPECIALIST_PATTERNS:
        if re.search(pat, s_str, re.IGNORECASE):
            return True
    if extra_patterns:
        for pat in extra_patterns:
            if pat and re.search(pat, s_str, re.IGNORECASE):
                return True
    return False


def _has_non_person_determination_token(val: Optional[Any]) -> bool:
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return False
    if pd is not None and pd.isna(val):
        return False
    s_raw = str(val).strip()
    s_low = s_raw.lower()
    if not s_low:
        return False
    if any(tok in s_low for tok in _NON_PERSON_DETERMINATION_TOKENS):
        return True
    return bool(_FLORA_CITATION_REGEX.search(s_raw))


def _is_missing_identified_by(val: Optional[Any]) -> bool:
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return True
    return str(val).strip().lower() in _MISSING_IDENTIFIED_BY_TOKENS


def _is_explicit_collector_reference(val: Optional[Any]) -> bool:
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return False
    s = re.sub(r"\s+", " ", str(val).strip())
    return bool(s and _COLLECTOR_REFERENCE_REGEX.match(s))


def _should_fallback_to_collector(identified_by_val: Optional[Any]) -> bool:
    if _has_non_person_determination_token(identified_by_val):
        return False
    return _is_missing_identified_by(identified_by_val) or _is_explicit_collector_reference(identified_by_val)


def normalize_determiner_surnames(identified_by_raw: Optional[Any]) -> List[str]:
    """
    Normalizes Darwin Core `identifiedBy` strings into clean ASCII-folded determiner
    surnames, filtering out dates, prefixes, initials, and generic non-person strings.
    """
    if (
        identified_by_raw is None
        or (isinstance(identified_by_raw, float) and math.isnan(identified_by_raw))
        or _is_missing_identified_by(identified_by_raw)
        or _is_explicit_collector_reference(identified_by_raw)
        or _has_non_person_determination_token(identified_by_raw)
    ):
        return []
    raw = str(identified_by_raw).strip()
    if not raw:
        return []
    raw_low = raw.lower()
    if raw_low in GENERIC_DETERMINER_TOKENS or any(
        bad in raw_low for bad in ("herbarium staff", "not evident", "unrecorded", "anonymous", "digitiz")
    ):
        return []

    # Strip annotation prefixes, affirmation marks, parenthetical notes, dates, and years
    cleaned = re.sub(r"\([^)]*\)", " ", raw)
    cleaned = re.sub(
        r"\b(?:det\.?\s*(?:by)?|determined\s+by|annot\.?\s*(?:by)?|annotated\s+by|ver\.?\s*(?:by)?|verified\s+by|conf\.?\s*(?:by)?|confirmed\s+by|rev\.?\s*(?:by)?)\b",
        " ",
        cleaned,
        flags=re.IGNORECASE,
    )
    cleaned = re.sub(
        r"[,\s]*(?:\d{1,2}\s+)?\b(?:jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|jun(?:e)?|jul(?:y)?|aug(?:ust)?|sep(?:t(?:ember)?)?|oct(?:ober)?|nov(?:ember)?|dec(?:ember)?)\.?\s*\d{0,4}\b",
        " ",
        cleaned,
        flags=re.IGNORECASE,
    )
    cleaned = re.sub(r"[,\s]*\b\d{4}(?:[-/]\d{1,2}(?:[-/]\d{1,2})?)?\b", " ", cleaned)
    cleaned = re.sub(r"\b\d+\b", " ", cleaned)
    cleaned = cleaned.replace("!", " ").strip(" ,;|-")
    if not cleaned:
        return []

    parts = re.split(r"\s*(?:&|;|\band\b|\bwith\b|\bet\s+al\.?)\s*", cleaned, flags=re.IGNORECASE)
    surnames: List[str] = []
    for part in parts:
        p = part.strip(" ,;|-.")
        if not p:
            continue
        p = re.sub(r"\b(?:dr|prof|mr|mrs|ms|jr|sr|ii|iii|iv)\b\.?", " ", p, flags=re.IGNORECASE).strip(" ,.")
        if not p:
            continue
        if "," in p:
            cand_raw = p.split(",", 1)[0].strip()
            cand_tokens = [re.sub(r"[^A-Za-zÀ-ÿ-]", "", t) for t in cand_raw.split()]
            cand_tokens = [t for t in cand_tokens if len(t) > 1]
            surname = _fold_ascii(cand_tokens[-1]) if cand_tokens else ""
        else:
            # Handle compact initials like 'Á.Löve' or 'D.K.Trock'
            p_spaced = re.sub(r"\.", ". ", p)
            tokens = [re.sub(r"[^A-Za-zÀ-ÿ-]", "", t) for t in p_spaced.split()]
            non_init = [t for t in tokens if len(t) > 1]
            surname = _fold_ascii(non_init[-1]) if non_init else ""

        surname = surname.strip("-")
        if (
            len(surname) >= 3
            and surname not in GENERIC_DETERMINER_TOKENS
            and surname not in {"collector", "collectors", "same"}
            and not any(tok in GENERIC_DETERMINER_TOKENS for tok in _fold_ascii(p).split() if tok in {"staff", "student", "students", "herbarium", "unknown", "anonymous"})
        ):
            if surname not in surnames:
                surnames.append(surname)

    return surnames


_HOLOTYPE_TIER_STATUSES: Set[str] = {"HOLOTYPE", "LECTOTYPE", "NEOTYPE", "EPITYPE"}
_ISOTYPE_TIER_STATUSES: Set[str] = {"ISOTYPE", "ISOLECTOTYPE", "ISONEOTYPE"}
_SYNTYPE_TIER_STATUSES: Set[str] = {"SYNTYPE", "ISOSYNTYPE"}
_PARATYPE_TIER_STATUSES: Set[str] = {"PARATYPE", "ISOPARATYPE", "TYPE", "COTYPE"}
_NEGATED_TYPE_RE = re.compile(r"\b(?:not\s+a?\s*type|non[-\s]?type)\b", re.IGNORECASE)


def classify_nomenclatural_type_rank(type_status_raw: Optional[Any]) -> Tuple[int, str]:
    """
    Classify a Darwin Core `typeStatus` value into a 3-tier nomenclatural type rank:
      - (2, cleaned_label) for PRIMARY_TYPE_STATUSES (Holotype, Lectotype, Neotype,
        Epitype, Isotype, Isolectotype, Isoneotype; `Tier_0_Primary_Type`)
      - (1, cleaned_label) for SECONDARY_TYPE_STATUSES (Syntype, Isosyntype,
        Paratype, Isoparatype, Type, Cotype; `Tier_0_Secondary_Type`)
      - (0, "None") for non-types or negated strings.
    """
    if type_status_raw is None or (isinstance(type_status_raw, float) and math.isnan(type_status_raw)):
        return 0, "None"
    if pd is not None and pd.isna(type_status_raw):
        return 0, "None"
    type_status = str(type_status_raw).strip()
    if not type_status or type_status.upper() in {
        "NONE", "NOT A TYPE", "NOTATYPE", "NON-TYPE", "NONTYPE",
        "UNSPECIFIED", "NULL", "NAN", "NA", "N/A", "FALSE", "0", "UNKNOWN",
    }:
        return 0, "None"
    if _NEGATED_TYPE_RE.search(type_status):
        return 0, "None"
    for prim_type in sorted(PRIMARY_TYPE_STATUSES, key=len, reverse=True):
        if re.search(rf"\b{prim_type}\b", type_status, re.IGNORECASE):
            return 2, prim_type.title()
    for sec_type in sorted(SECONDARY_TYPE_STATUSES, key=len, reverse=True):
        if re.search(rf"\b{sec_type}\b", type_status, re.IGNORECASE):
            return 1, sec_type.title()
    return 0, "None"


def classify_nomenclatural_type_subrank(type_status_raw: Optional[Any]) -> float:
    """
    Classify a Darwin Core `typeStatus` value into an ICN Art. 9.6 fine-grained float sub-rank:
      - 2.2: Holotype-tier (HOLOTYPE, LECTOTYPE, NEOTYPE, EPITYPE)
      - 2.0: Isotype-tier (ISOTYPE, ISOLECTOTYPE, ISONEOTYPE)
      - 1.5: Syntype-tier (SYNTYPE, ISOSYNTYPE; jointly name-bearing under ICN Art. 9.6)
      - 1.0: Paratype-tier (PARATYPE, ISOPARATYPE, TYPE, COTYPE; non-name-bearing)
      - 0.0: Non-types or negated strings.
    """
    if type_status_raw is None or (isinstance(type_status_raw, float) and math.isnan(type_status_raw)):
        return 0.0
    if pd is not None and pd.isna(type_status_raw):
        return 0.0
    type_status = str(type_status_raw).strip()
    if not type_status or type_status.upper() in {
        "NONE", "NOT A TYPE", "NOTATYPE", "NON-TYPE", "NONTYPE",
        "UNSPECIFIED", "NULL", "NAN", "NA", "N/A", "FALSE", "0", "UNKNOWN",
    }:
        return 0.0
    if _NEGATED_TYPE_RE.search(type_status):
        return 0.0
    for tier_val, tier_set in (
        (2.2, _HOLOTYPE_TIER_STATUSES),
        (2.0, _ISOTYPE_TIER_STATUSES),
        (1.5, _SYNTYPE_TIER_STATUSES),
        (1.0, _PARATYPE_TIER_STATUSES),
    ):
        for tok in sorted(tier_set, key=len, reverse=True):
            if re.search(rf"\b{tok}\b", type_status, re.IGNORECASE):
                return tier_val
    return 0.0


def _is_verified_type_status(type_status_raw: Optional[Any]) -> Tuple[bool, str]:
    """Return (is_type, cleaned_type_label) for a Darwin Core typeStatus value."""
    t_rank, t_lbl = classify_nomenclatural_type_rank(type_status_raw)
    return (t_rank >= 1), t_lbl


def discover_corpus_authorities(
    df: pd.DataFrame,
    min_gold_determinations: Optional[int] = None,
    min_gold_herbaria: int = 3,
    min_herbarium_vouchers: int = 10,
) -> Tuple[Set[str], Set[str]]:
    """
    Layer 3 Authority Engine: Discovers empirical taxonomic specialists and regional
    research herbaria directly from the assembled Darwin Core voucher corpus.

    Rules:
      - Dynamic Specialist Discovery (`dynamic_specialists`):
        1. Determiner (`identifiedBy` only, never falling back to `recordedBy` or
           matching `_NON_PERSON_DETERMINATION_TOKENS`) has identified >= `min_gold_determinations`
           specimens across >= `min_gold_herbaria` distinct `institutionCode`s (and in full
           corpora with >= 100 rows, either >= 4 distinct Packera species or >= 1 type specimen), OR
        2. Determiner has explicitly annotated >= 1 verified type specimen (giving +2
           weight to type determinations, excluding collector-reference fallbacks) AND
           reaches weighted count >= 3 across >= 2 distinct `institutionCode`s, AND
        3. Bare surnames in `_AMBIGUOUS_COMMON_SURNAMES` are excluded unless already in
           `PACKERA_SPECIALISTS` or backed by >= 1 verified type determination, AND
        4. Determiner has cross-duplicate disagreement rate <= 15% across shared
           `exsiccatae_key` duplicate groups when >= 2 duplicate groups exist.
      - Dynamic Regional Research Herbarium Discovery (`dynamic_herbaria`):
        1. `institutionCode` holds >= 1 verified type specimen for target Packera taxa
           AND >= 5 total sheets in the corpus, OR
        2. `institutionCode` ranks in the top quartile of institutional holdings
           (>= `min_herbarium_vouchers` curated sheets) with >= 70% determiner-attributed
           records for any harvested taxon.

    Returns:
        Tuple[Set[str], Set[str]]: (dynamic_specialists, dynamic_herbaria)
    """
    dynamic_specialists: Set[str] = set()
    dynamic_herbaria: Set[str] = set()
    if df is None or df.empty:
        return dynamic_specialists, dynamic_herbaria

    is_large_corpus = len(df) >= 100
    eff_min_det = (15 if is_large_corpus else 5) if min_gold_determinations is None else int(min_gold_determinations)
    eff_min_herbs = max(1, int(min_gold_herbaria))
    eff_min_herb_vouchers = max(1, int(min_herbarium_vouchers))

    id_col = (
        "_raw_identified_by"
        if "_raw_identified_by" in df.columns
        else ("identifiedBy" if "identifiedBy" in df.columns else None)
    )
    sp_col = (
        "species_standardized"
        if "species_standardized" in df.columns
        else ("scientificName" if "scientificName" in df.columns else ("species_raw" if "species_raw" in df.columns else None))
    )
    type_col = (
        "_raw_type_status"
        if "_raw_type_status" in df.columns
        else ("type_status" if "type_status" in df.columns else ("typeStatus" if "typeStatus" in df.columns else None))
    )
    inst_col = "institutionCode" if "institutionCode" in df.columns else None
    has_exsiccatae_col = "exsiccatae_key" in df.columns

    # Track (determiner_surname, standardized_taxon) -> [count, set_of_institutions]
    det_taxon_counts: Dict[Tuple[str, str], int] = {}
    det_taxon_insts: Dict[Tuple[str, str], Set[str]] = {}
    det_total_counts: Dict[str, int] = {}
    det_total_insts: Dict[str, Set[str]] = {}
    det_species: Dict[str, Set[str]] = {}
    det_type_hits: Dict[str, int] = {}
    det_type_annotated: Set[str] = set()

    # Track institutional holdings, type holdings, and determiner-attribution rates
    inst_total_sheets: Dict[str, int] = {}
    inst_type_sheets: Dict[str, int] = {}
    inst_taxon_totals: Dict[Tuple[str, str], int] = {}
    inst_taxon_attributed: Dict[Tuple[str, str], int] = {}

    from scripts.core.harvester import _extract_canonical_binomial, standardize_packera_taxon

    # Track duplicate groups for cross-duplicate disagreement auditing (<= 15% threshold)
    group_entries: Dict[str, List[Tuple[str, List[str], bool]]] = {}

    for _, row in df.iterrows():
        raw_sp = row.get(sp_col) if sp_col else None
        std_sp = standardize_packera_taxon(raw_sp)
        if not std_sp and raw_sp is not None and not pd.isna(raw_sp):
            s_tmp = str(raw_sp).strip()
            if s_tmp.startswith("Packera "):
                std_sp = _extract_canonical_binomial(s_tmp, {"Packera"})
        if not std_sp:
            std_sp = _resolve_row_standardized_packera(row)

        raw_inst = row.get(inst_col) if inst_col else ""
        inst = str(raw_inst if not pd.isna(raw_inst) else "").strip().upper()
        valid_inst = inst not in _INVALID_INSTITUTION_CODES

        raw_id = row.get(id_col) if id_col else ""
        is_explicit_col_ref = _is_explicit_collector_reference(raw_id)
        is_non_person_id = _has_non_person_determination_token(raw_id)
        has_explicit_det_str = (
            raw_id is not None
            and not (pd is not None and pd.isna(raw_id))
            and bool(str(raw_id).strip())
            and not is_explicit_col_ref
            and not is_non_person_id
        )
        # Never fall back to recordedBy when counting determinations in discover_corpus_authorities
        surnames = normalize_determiner_surnames(raw_id) if has_explicit_det_str else []
        has_attributed_det = len(surnames) > 0
        raw_rec_by = row.get("recordedBy", "") if "recordedBy" in row else ""
        rec_col_norm = normalize_collector("" if (raw_rec_by is None or pd.isna(raw_rec_by)) else str(raw_rec_by))
        cross_collector_surnames = [
            s_n for s_n in surnames
            if not (rec_col_norm and rec_col_norm != "unknown" and s_n == rec_col_norm)
        ]

        raw_type = row.get(type_col) if type_col else None
        is_type, _ = _is_verified_type_status(raw_type)

        if is_type:
            if valid_inst and inst != "UNKNOWN" and std_sp:
                inst_type_sheets[inst] = inst_type_sheets.get(inst, 0) + 1
            if has_explicit_det_str:
                for s_name in cross_collector_surnames:
                    det_type_annotated.add(s_name)
                    det_type_hits[s_name] = det_type_hits.get(s_name, 0) + 1

        if valid_inst:
            inst_total_sheets[inst] = inst_total_sheets.get(inst, 0) + 1
            if std_sp:
                it_key = (inst, std_sp)
                inst_taxon_totals[it_key] = inst_taxon_totals.get(it_key, 0) + 1
                if has_attributed_det:
                    inst_taxon_attributed[it_key] = inst_taxon_attributed.get(it_key, 0) + 1

        if std_sp and cross_collector_surnames:
            det_weight = 2 if (is_type and has_explicit_det_str) else 1
            for s_name in cross_collector_surnames:
                dt_key = (s_name, std_sp)
                det_taxon_counts[dt_key] = det_taxon_counts.get(dt_key, 0) + det_weight
                det_total_counts[s_name] = det_total_counts.get(s_name, 0) + det_weight
                det_species.setdefault(s_name, set()).add(std_sp)
                if valid_inst:
                    det_taxon_insts.setdefault(dt_key, set()).add(inst)
                    det_total_insts.setdefault(s_name, set()).add(inst)

        dup_key = ""
        if has_exsiccatae_col:
            raw_ek = row.get("exsiccatae_key")
            if raw_ek is not None and not pd.isna(raw_ek):
                dup_key = str(raw_ek).strip()
        if not dup_key and "recordNumber" in row:
            r_num = row.get("recordNumber", "")
            norm_num = normalize_record_number("" if pd.isna(r_num) else r_num)
            if norm_num:
                r_by = row.get("recordedBy", "")
                norm_col = normalize_collector("" if pd.isna(r_by) else str(r_by))
                ev_dt = row.get("eventDate", "")
                norm_dt = normalize_event_date(
                    "" if pd.isna(ev_dt) else str(ev_dt),
                    year=row.get("year"),
                    month=row.get("month"),
                    day=row.get("day"),
                )
                dup_key = f"{norm_col}_{norm_num}_{norm_dt}"
        if dup_key and std_sp:
            pre_conf = bool(
                "exsiccatae_taxonomic_conflict" in row
                and not pd.isna(row.get("exsiccatae_taxonomic_conflict"))
                and bool(row.get("exsiccatae_taxonomic_conflict"))
            )
            group_entries.setdefault(dup_key, []).append((std_sp, surnames, pre_conf))

    # Build cross-duplicate conflict counts per determiner before promotion
    det_dup_groups: Dict[str, Set[str]] = {}
    det_conflicted_groups: Dict[str, Set[str]] = {}
    for dup_key, entries in group_entries.items():
        if len(entries) < 2:
            continue
        distinct_group_taxa = {e[0] for e in entries if e[0]}
        is_conflicted_group = (len(distinct_group_taxa) >= 2) or any(e[2] for e in entries)
        for _, entry_surnames, _ in entries:
            for s_name in entry_surnames:
                det_dup_groups.setdefault(s_name, set()).add(dup_key)
                if is_conflicted_group:
                    det_conflicted_groups.setdefault(s_name, set()).add(dup_key)

    def _passes_specialist_gates(s_name: str) -> bool:
        type_hits = det_type_hits.get(s_name, 0)
        if s_name in _AMBIGUOUS_COMMON_SURNAMES and s_name not in PACKERA_SPECIALISTS and type_hits < 1:
            return False
        dup_groups = len(det_dup_groups.get(s_name, set()))
        conflicted_groups = len(det_conflicted_groups.get(s_name, set()))
        if dup_groups >= 2 and (conflicted_groups / float(dup_groups)) > 0.15:
            return False
        if is_large_corpus and type_hits < 1 and len(det_species.get(s_name, set())) < 4:
            return False
        return True

    for dt_key, count in det_taxon_counts.items():
        s_name = dt_key[0]
        if not _passes_specialist_gates(s_name):
            continue
        inst_set = det_taxon_insts.get(dt_key, set())
        if (count >= eff_min_det and len(inst_set) >= eff_min_herbs) or (
            s_name in det_type_annotated and count >= 3 and len(inst_set) >= 2
        ):
            dynamic_specialists.add(s_name)

    for s_name, tot_count in det_total_counts.items():
        if s_name in dynamic_specialists or not _passes_specialist_gates(s_name):
            continue
        tot_insts = det_total_insts.get(s_name, set())
        if (tot_count >= eff_min_det and len(tot_insts) >= eff_min_herbs) or (
            s_name in det_type_annotated and tot_count >= 3 and len(tot_insts) >= 2
        ):
            dynamic_specialists.add(s_name)

    if inst_total_sheets:
        for inst, n_types in inst_type_sheets.items():
            if n_types >= 1 and inst_total_sheets.get(inst, 0) >= 5:
                dynamic_herbaria.add(inst)
        q75_holdings = float(np.percentile(list(inst_total_sheets.values()), 75))
        for inst, total_n in inst_total_sheets.items():
            if total_n >= eff_min_herb_vouchers and total_n >= q75_holdings:
                min_sp_sheets = min(3, total_n)
                for (i_code, sp_name), sp_total in inst_taxon_totals.items():
                    if i_code == inst and sp_total >= min_sp_sheets:
                        n_attr = inst_taxon_attributed.get((i_code, sp_name), 0)
                        if (n_attr / float(sp_total)) >= 0.70:
                            dynamic_herbaria.add(inst)
                            break

    return dynamic_specialists, dynamic_herbaria


def parse_determiner_tier(
    type_status_raw: Optional[str],
    identified_by_raw: Optional[str],
    recorded_by_raw: Optional[str],
    history_raw: Optional[str],
    institution_code_raw: Optional[str],
    locality_raw: Optional[str],
    habitat_raw: Optional[str],
    dynamic_specialists: Optional[Iterable[str]] = None,
    dynamic_herbaria: Optional[Iterable[str]] = None,
    is_open_or_hybrid: bool = False,
) -> Tuple[str, str, str]:
    """
    Evaluates Darwin Core fields to assign a taxonomic determination credibility tier
    using the 3-Layer Dynamic Expert & Herbarium Authority Engine:
      - Layer 1: Genus-wide Packera / Senecioneae monographer patterns (`SPECIALIST_PATTERNS`).
      - Layer 2 & 3: Dynamic protologue authors and corpus-derived multi-institution
        specialists (`dynamic_specialists`) and regional research herbaria (`dynamic_herbaria`).

    Tiers:
      - Tier_1_Gold: Primary/secondary nomenclatural types, monographers, protologue authors,
        or corpus-verified multi-herbarium specialists.
      - Tier_2_Silver: Determinations from major or corpus-qualified research herbaria with
        complete ecological locality and attributed determiner.
      - Tier_3_Bronze: Unverified candidate collections, general collectors, or missing determiner.

    Returns:
        Tuple[str, str, str]: (determiner_tier, type_status_clean, identified_by_combined)
    """
    type_status = "" if (type_status_raw is None or (isinstance(type_status_raw, float) and math.isnan(type_status_raw))) else str(type_status_raw).strip()
    identified_by = "" if (identified_by_raw is None or (isinstance(identified_by_raw, float) and math.isnan(identified_by_raw))) else str(identified_by_raw).strip()
    recorded_by = "" if (recorded_by_raw is None or (isinstance(recorded_by_raw, float) and math.isnan(recorded_by_raw))) else str(recorded_by_raw).strip()
    history = "" if (history_raw is None or (isinstance(history_raw, float) and math.isnan(history_raw))) else str(history_raw).strip()
    institution = "" if (institution_code_raw is None or (isinstance(institution_code_raw, float) and math.isnan(institution_code_raw))) else str(institution_code_raw).strip().upper()
    locality = "" if (locality_raw is None or (isinstance(locality_raw, float) and math.isnan(locality_raw))) else str(locality_raw).strip()
    habitat = "" if (habitat_raw is None or (isinstance(habitat_raw, float) and math.isnan(habitat_raw))) else str(habitat_raw).strip()

    type_rank, type_status_clean = classify_nomenclatural_type_rank(type_status)
    is_type = type_rank > 0

    is_non_person_det = _has_non_person_determination_token(identified_by)
    has_collector_ref = False if is_non_person_det else _is_explicit_collector_reference(identified_by)
    has_explicit_det = bool(
        not is_non_person_det
        and identified_by
        and not _is_missing_identified_by(identified_by)
        and not has_collector_ref
    )
    fallback_to_collector = (
        False
        if (is_non_person_det and type_rank <= 0)
        else _should_fallback_to_collector(identified_by)
    )

    if has_collector_ref and recorded_by:
        identified_by_combined = recorded_by
    elif has_explicit_det:
        identified_by_combined = identified_by
    else:
        identified_by_combined = history if history else (recorded_by if recorded_by else identified_by)

    # Quarantine non-type open-nomenclature and hybrid-formula sheets strictly to Tier_3_Bronze
    if is_open_or_hybrid and not is_type:
        return "Tier_3_Bronze", type_status_clean, identified_by_combined

    # Match specialists strictly on current determination (`identifiedBy`), only consulting
    # `verbatimIdentificationHistory` when `identifiedBy` is absent/collector-ref (so an older
    # superseded specialist annotation in history never promotes a newer non-specialist override).
    eff_id_for_specialist = "" if is_non_person_det else identified_by
    eff_hist_for_specialist = (
        ""
        if (has_explicit_det or _has_non_person_determination_token(history))
        else history
    )

    is_core_monographer = False
    is_regional_specialist_silver = False
    if not (is_non_person_det and not eff_hist_for_specialist and type_rank <= 0):
        explicit_det_text = f"{eff_id_for_specialist} | {eff_hist_for_specialist}".strip(" |")
        if explicit_det_text and _has_core_monographer_match(explicit_det_text):
            is_core_monographer = True
        elif explicit_det_text and (
            _has_regional_or_historical_specialist_match(explicit_det_text)
            or _has_specialist_match(explicit_det_text)
        ):
            is_regional_specialist_silver = True
        elif fallback_to_collector and recorded_by:
            if _has_core_monographer_match(recorded_by):
                is_core_monographer = True
            elif _has_regional_or_historical_specialist_match(recorded_by) or _has_specialist_match(recorded_by):
                # Regional flora collectors on blank-identifiedBy field labels route to Silver
                is_regional_specialist_silver = True

    # Corpus-discovered multi-herbarium specialists cap at Tier_2_Silver (never Tier_1_Gold),
    # and apply ONLY to explicit determinations/unoverridden history (never collector fallback).
    is_dynamic_corpus_specialist_silver = False
    if not is_core_monographer and not is_regional_specialist_silver and dynamic_specialists and (has_explicit_det or eff_hist_for_specialist):
        dyn_spec_set: Set[str] = set()
        for item in dynamic_specialists:
            if item is None:
                continue
            s_tok = str(item).strip().lower()
            if len(s_tok) >= 3 and s_tok not in _AMBIGUOUS_COMMON_SURNAMES:
                dyn_spec_set.add(s_tok)
                dyn_spec_set.add(_fold_ascii(s_tok))
        if dyn_spec_set:
            det_sources = [eff_id_for_specialist if has_explicit_det else "", eff_hist_for_specialist]
            det_surnames: Set[str] = set()
            for src in det_sources:
                det_surnames.update(normalize_determiner_surnames(src))
            if det_surnames & dyn_spec_set:
                is_dynamic_corpus_specialist_silver = True
            else:
                det_audit_raw = " | ".join(det_sources).lower()
                det_audit_ascii = _fold_ascii(det_audit_raw)
                for spec_tok in dyn_spec_set:
                    pat = rf"(?<![a-zA-ZÀ-ÿ]){re.escape(spec_tok)}(?![a-zA-ZÀ-ÿ])"
                    if re.search(pat, det_audit_raw) or re.search(pat, det_audit_ascii):
                        is_dynamic_corpus_specialist_silver = True
                        break

    if is_type or is_core_monographer:
        return "Tier_1_Gold", type_status_clean, identified_by_combined

    if is_regional_specialist_silver or is_dynamic_corpus_specialist_silver:
        return "Tier_2_Silver", type_status_clean, identified_by_combined

    if "[exsiccatae " in identified_by.lower():
        return "Tier_2_Silver", type_status_clean, identified_by_combined

    dyn_herb_set = {str(h).strip().upper() for h in dynamic_herbaria if h} if dynamic_herbaria else set()
    is_major_herbarium = (institution in MAJOR_HERBARIA_CODES) or (institution in dyn_herb_set)
    has_rich_locality = len(locality) > 5 or len(habitat) > 3
    eff_silver_det = "" if is_non_person_det else (recorded_by if (has_collector_ref and recorded_by) else identified_by)
    silver_low = eff_silver_det.lower()
    is_student_or_intern = bool(re.search(r"\b(?:student|undergrad(?:uate)?|intern|class|volunteer)\b", silver_low))
    has_determiner = bool(
        eff_silver_det
        and not is_student_or_intern
        and silver_low not in GENERIC_DETERMINER_TOKENS
        and not _has_non_person_determination_token(eff_silver_det)
        and len(normalize_determiner_surnames(eff_silver_det)) > 0
    )

    if is_major_herbarium and has_rich_locality and has_determiner:
        return "Tier_2_Silver", type_status_clean, identified_by_combined

    return "Tier_3_Bronze", type_status_clean, identified_by_combined


def apply_corpus_authority_tiers(
    df: pd.DataFrame,
    dynamic_specialists: Optional[Set[str]] = None,
    dynamic_herbaria: Optional[Set[str]] = None,
    db_path: Optional[Union[str, Path]] = None,
    workspace_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Applies the 3-Layer Dynamic Authority Engine across an assembled DataFrame:
      - Discovers corpus-level multi-herbarium specialists and dynamic regional herbaria.
      - Collects per-taxon protologue authors from Darwin Core `scientificName` / `species_raw`.
      - Detects undated institutional batch transfers (>= 25 sheets at 1 institution across >= 10
        collectors) and demotes non-type historical sheets collected by others from `Tier_1_Gold`
        to `Tier_2_Silver` for out-of-sample screening, while exempting recognized *Packera*
        specialists who annotate across multiple herbaria or multiple species in a batch.
      - Re-evaluates `determiner_tier` across all rows cleanly without downgrading pre-scored
        Silver rows whose `locality` / `habitat` columns were already omitted upon CSV persistence.
    """
    if df is None or df.empty:
        return df

    out = df.copy()
    if dynamic_specialists is None or dynamic_herbaria is None:
        disc_specs, disc_herbs = discover_corpus_authorities(out)
        eff_specs = disc_specs if dynamic_specialists is None else (set(dynamic_specialists) | disc_specs)
        eff_herbs = disc_herbs if dynamic_herbaria is None else (set(dynamic_herbaria) | disc_herbs)
    else:
        eff_specs = set(dynamic_specialists)
        eff_herbs = set(dynamic_herbaria)

    from scripts.core.harvester import standardize_packera_taxon

    # Tally (institutionCode, identifiedBy) batch statistics to detect undated institutional batch transfers
    batch_voucher_counts: Dict[Tuple[str, str], int] = {}
    batch_collectors: Dict[Tuple[str, str], Set[str]] = {}
    batch_dates: Dict[Tuple[str, str], Set[str]] = {}
    det_institutions: Dict[str, Set[str]] = {}
    surname_institutions: Dict[str, Set[str]] = {}
    key_species_counts: Dict[Tuple[str, str], Dict[str, int]] = {}

    for _, row in out.iterrows():
        raw_names_row = [
            row.get("scientificName"),
            row.get("species_raw"),
        ]
        std_sp = None
        for rn in raw_names_row:
            std_sp = standardize_packera_taxon(rn)
            if std_sp:
                break
        if not std_sp:
            std_sp = _resolve_row_standardized_packera(row)

        inst_r = row.get("institutionCode")
        inst_key = "" if (inst_r is None or pd.isna(inst_r)) else str(inst_r).strip().upper()
        id_r = row.get("_raw_identified_by") if "_raw_identified_by" in out.columns else row.get("identifiedBy")
        id_str_clean = "" if (id_r is None or pd.isna(id_r)) else str(id_r).strip()
        id_key = re.sub(r"\s+", " ", id_str_clean.lower()) if id_str_clean else ""
        if (
            inst_key
            and inst_key not in _INVALID_INSTITUTION_CODES
            and id_key
            and not _is_missing_identified_by(id_key)
            and not _is_explicit_collector_reference(id_key)
            and not _has_non_person_determination_token(id_key)
        ):
            b_key = (inst_key, id_key)
            batch_voucher_counts[b_key] = batch_voucher_counts.get(b_key, 0) + 1
            det_institutions.setdefault(id_key, set()).add(inst_key)
            for s_n in normalize_determiner_surnames(id_str_clean):
                surname_institutions.setdefault(s_n, set()).add(inst_key)
            sp_val = str(std_sp or "").strip().lower()
            if sp_val:
                sp_map = key_species_counts.setdefault(b_key, {})
                sp_map[sp_val] = sp_map.get(sp_val, 0) + 1
            rec_r = row.get("recordedBy")
            col_norm = normalize_collector("" if (rec_r is None or pd.isna(rec_r)) else str(rec_r))
            if col_norm and col_norm != "unknown":
                batch_collectors.setdefault(b_key, set()).add(col_norm)
            dt_id_r = row.get("_raw_date_identified") if "_raw_date_identified" in out.columns else row.get("dateIdentified")
            dt_id_str = "" if (dt_id_r is None or pd.isna(dt_id_r)) else str(dt_id_r).strip()
            if dt_id_str and dt_id_str.lower() not in {"none", "nan", "null", "na", "n/a", "unknown"}:
                batch_dates.setdefault(b_key, set()).add(dt_id_str)

    flagged_batch_keys: Set[Tuple[str, str]] = set()
    for b_key, v_cnt in batch_voucher_counts.items():
        if v_cnt < 25 or len(batch_collectors.get(b_key, set())) < 10:
            continue
        _inst_k, _id_k = b_key
        if _has_specialist_match(_id_k):
            inst_set_for_det = set(det_institutions.get(_id_k, set()))
            for _sn in normalize_determiner_surnames(_id_k):
                if _sn in PACKERA_SPECIALISTS:
                    inst_set_for_det |= surname_institutions.get(_sn, set())
            n_inst_for_det = len(inst_set_for_det)
            sp_counts = key_species_counts.get(b_key, {})
            n_species_in_batch = len(sp_counts)
            max_sp_share = (
                (max(sp_counts.values()) / max(1, sum(sp_counts.values())))
                if sp_counts
                else 1.0
            )
            if n_inst_for_det >= 2 or (n_species_in_batch >= 2 and max_sp_share < 0.85):
                continue
        flagged_batch_keys.add(b_key)

    has_determiner_cols = any(
        c in out.columns
        for c in (
            "_raw_identified_by",
            "identifiedBy",
            "recordedBy",
            "_raw_type_status",
            "type_status",
            "typeStatus",
            "_raw_history",
            "verbatimIdentificationHistory",
        )
    )
    if not has_determiner_cols and "determiner_tier" in out.columns:
        return out

    has_locality_cols = any(
        c in out.columns for c in ("_raw_locality", "locality", "verbatimLocality", "_raw_habitat", "habitat")
    )

    new_tiers: List[str] = []
    new_types: List[str] = []
    update_type_col = "type_status" in out.columns

    for _, row in out.iterrows():
        t_raw = row.get("_raw_type_status") if "_raw_type_status" in out.columns else row.get("type_status", row.get("typeStatus"))
        id_raw = row.get("_raw_identified_by") if "_raw_identified_by" in out.columns else row.get("identifiedBy")
        rec_raw = row.get("recordedBy")
        hist_raw = row.get("_raw_history") if "_raw_history" in out.columns else row.get("verbatimIdentificationHistory")
        inst_raw = row.get("institutionCode")
        loc_raw = (
            row.get("_raw_locality")
            if "_raw_locality" in out.columns
            else (row.get("locality") if "locality" in out.columns else row.get("verbatimLocality"))
        )
        hab_raw = row.get("_raw_habitat") if "_raw_habitat" in out.columns else row.get("habitat")

        raw_sci = row.get("scientificName") or row.get("species_raw")
        std_sp = standardize_packera_taxon(raw_sci) or _resolve_row_standardized_packera(row)
        c_stat_r = str(row.get("concept_synonymy_status", "") if not pd.isna(row.get("concept_synonymy_status", "")) else "").strip()
        row_open_hyb = _as_bool_flag(row.get("is_open_or_hybrid_concept")) or (c_stat_r in {"Open_Nomenclature", "Hybrid_Formula"})

        tier_eval, type_clean, _ = parse_determiner_tier(
            type_status_raw=t_raw,
            identified_by_raw=id_raw,
            recorded_by_raw=rec_raw,
            history_raw=hist_raw,
            institution_code_raw=inst_raw,
            locality_raw=loc_raw,
            habitat_raw=hab_raw,
            dynamic_specialists=eff_specs,
            dynamic_herbaria=eff_herbs,
            is_open_or_hybrid=row_open_hyb,
        )

        prev_tier = str(row.get("determiner_tier", "") if not pd.isna(row.get("determiner_tier", "")) else "").strip()
        ex_gold = str(row.get("exsiccatae_gold_taxon", "") if not pd.isna(row.get("exsiccatae_gold_taxon", "")) else "").strip()
        is_exsiccatae_transfer = "[exsiccatae " in str(id_raw or "").lower()
        t_rank_row, _ = classify_nomenclatural_type_rank(t_raw)
        if row_open_hyb and t_rank_row == 0:
            tier_eval = "Tier_3_Bronze"
        elif is_exsiccatae_transfer:
            tier_eval = "Tier_2_Silver"
        elif prev_tier in {"Tier_1_Gold", "Tier_2_Silver"} and ex_gold and std_sp == ex_gold and tier_eval == "Tier_3_Bronze":
            tier_eval = prev_tier
        elif not has_locality_cols and prev_tier == "Tier_2_Silver" and tier_eval == "Tier_3_Bronze":
            tier_eval = "Tier_2_Silver"

        # Demote undated (or batch-identical-dated) institutional batch transfers on historical sheets
        # collected by others (type_rank == 0) from Tier_1_Gold to Tier_2_Silver
        if tier_eval == "Tier_1_Gold" and not (ex_gold and std_sp == ex_gold):
            t_rank, _ = classify_nomenclatural_type_rank(t_raw)
            if t_rank == 0:
                inst_key = "" if (inst_raw is None or pd.isna(inst_raw)) else str(inst_raw).strip().upper()
                id_str_row = "" if (id_raw is None or pd.isna(id_raw)) else str(id_raw).strip()
                id_key = re.sub(r"\s+", " ", id_str_row.lower()) if id_str_row else ""
                b_key = (inst_key, id_key)
                if b_key in flagged_batch_keys:
                    dt_id_r = row.get("_raw_date_identified") if "_raw_date_identified" in out.columns else row.get("dateIdentified")
                    dt_id_str = "" if (dt_id_r is None or pd.isna(dt_id_r)) else str(dt_id_r).strip()
                    has_embedded_year = bool(id_str_row and re.search(r"\b(19|20)\d{2}\b", id_str_row))
                    has_date = bool(dt_id_str and re.search(r"\b(19|20)\d{2}\b", dt_id_str)) or has_embedded_year
                    has_blank_or_batch_date = (
                        not has_date
                        or (not has_embedded_year and len(batch_dates.get(b_key, set())) <= 1)
                    )
                    if has_blank_or_batch_date:
                        det_surnames = set(normalize_determiner_surnames(id_raw))
                        rec_str = "" if (rec_raw is None or pd.isna(rec_raw)) else str(rec_raw)
                        rec_surnames = set(normalize_determiner_surnames(rec_str))
                        rec_primary = normalize_collector(rec_str)
                        if rec_primary and rec_primary != "unknown":
                            rec_surnames.add(rec_primary)
                        rec_ascii = _fold_ascii(rec_str)
                        collector_matches_det = bool(det_surnames & rec_surnames) or any(
                            bool(re.search(rf"(?<![a-zA-ZÀ-ÿ]){re.escape(ds)}(?![a-zA-ZÀ-ÿ])", rec_ascii))
                            for ds in det_surnames
                        )
                        if not collector_matches_det:
                            tier_eval = "Tier_2_Silver"

        new_tiers.append(tier_eval)
        prev_type = str(row.get("type_status", "") if not pd.isna(row.get("type_status", "")) else "").strip()
        new_types.append(type_clean if (type_clean != "None" or not prev_type) else prev_type)

    out["determiner_tier"] = new_tiers
    if update_type_col:
        out["type_status"] = new_types
    return apply_concept_and_segregate_recency_caps(out, db_path=db_path, workspace_dir=workspace_dir)


def resolve_geographic_scope_settings(
    geographic_scope: Optional[str] = DEFAULT_GEOGRAPHIC_SCOPE,
    exclude_western: bool = True,
    include_western: bool = False,
    countries: Optional[Union[str, List[str]]] = None,
) -> Tuple[str, bool, List[str]]:
    """
    Resolves effective `(geographic_scope, exclude_western, countries)` settings:
      - `eastern_na` (default): `exclude_western=True` (unless `include_western=True` or
        `exclude_western=False`), `countries=["US", "CA"]` from `DEFAULT_COUNTRIES`
        (or `["US", "CA", "MX"]` when `include_western=True`).
      - `north_america`: automatically sets `exclude_western=False` and expands default
        `countries` to `["US", "CA", "MX"]` unless explicit `countries` were supplied.
      - `global`: automatically sets `exclude_western=False` and sets `countries=[]` (unrestricted
        global GBIF search) unless explicit `countries` were supplied.
    """
    scope = str(geographic_scope or DEFAULT_GEOGRAPHIC_SCOPE).strip().lower()
    if scope not in VALID_GEOGRAPHIC_SCOPES:
        raise ValueError(f"Unsupported geographic_scope={geographic_scope!r}; expected one of {sorted(VALID_GEOGRAPHIC_SCOPES)}.")

    eff_exclude_western = bool(exclude_western) and (not bool(include_western)) and (scope == "eastern_na")

    if countries is not None:
        eff_countries = normalize_countries_list(countries, default=[])
    elif scope == "global":
        eff_countries = []
    elif scope == "north_america" or bool(include_western):
        eff_countries = ["US", "CA", "MX"]
    else:
        eff_countries = normalize_countries_list(DEFAULT_COUNTRIES, default=["US", "CA"])

    return scope, eff_exclude_western, eff_countries


def is_excluded_western_region(
    state_province: Optional[str],
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    western_longitude_threshold: float = WESTERN_LONGITUDE_THRESHOLD,
    exclude_western: bool = True,
    is_type_specimen: bool = False,
) -> bool:
    """
    Determines whether a record originates from a western state/province (`EXCLUDED_WESTERN_STATES`)
    or west of `western_longitude_threshold` (default `-104.05`), while exempting Eastern Canadian
    provinces, Washington D.C., and verified nomenclatural type specimens (`is_type_specimen=True`).
    """
    if not exclude_western:
        return False

    has_valid_state = False
    if state_province is not None and not (isinstance(state_province, float) and math.isnan(state_province)):
        if pd is None or not pd.isna(state_province):
            raw_state = str(state_province).strip()
            cleaned = re.sub(r"\((?:state|province|territory)\)", "", raw_state, flags=re.IGNORECASE).strip(" ._,-")
            upper_state = cleaned.upper()
            if upper_state and upper_state not in {"NONE", "NAN", "NA", "N/A", "NULL", "UNKNOWN", "UNSPECIFIED", "UNRECORDED"}:
                has_valid_state = True
                if upper_state in {
                    "WASHINGTON, D.C.", "WASHINGTON D.C.", "WASHINGTON DC", "DISTRICT OF COLUMBIA", "DC",
                    "ON", "ONTARIO", "QC", "QUEBEC", "QUÉBEC", "NB", "NEW BRUNSWICK",
                    "NS", "NOVA SCOTIA", "PE", "PEI", "PRINCE EDWARD ISLAND",
                    "NL", "NEWFOUNDLAND", "NEWFOUNDLAND AND LABRADOR",
                }:
                    return False

                if upper_state in EXCLUDED_WESTERN_STATES:
                    return True

    has_valid_lon = False
    if lon is not None:
        try:
            lon_val = float(lon)
            if math.isfinite(lon_val):
                has_valid_lon = True
                if lon_val < float(western_longitude_threshold):
                    return True
        except (ValueError, TypeError):
            pass

    if not has_valid_lon and not has_valid_state and not is_type_specimen:
        return True

    return False


def _clean_upper_state_province(state_province: Optional[Any]) -> str:
    if state_province is None or (isinstance(state_province, float) and math.isnan(state_province)):
        return ""
    if pd is not None and pd.isna(state_province):
        return ""
    raw_state = str(state_province).strip()
    cleaned = re.sub(r"\((?:state|province|territory)\)", "", raw_state, flags=re.IGNORECASE).strip(" ._,-")
    up = cleaned.upper()
    if up in {"NONE", "NAN", "NA", "N/A", "NULL", "UNKNOWN", "UNSPECIFIED", "UNRECORDED"}:
        return ""
    return up


def is_cis_mississippi_record(
    state_province: Optional[str],
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    cis_longitude_threshold: float = DEFAULT_CIS_MISSISSIPPI_LON,
) -> bool:
    """
    Returns True strictly when an occurrence lies east of the Mississippi River
    (`CIS_MISSISSIPPI_STATES_PROVINCES`, or `lon >= cis_longitude_threshold` when state is missing).
    Trans-Mississippi states (MN, IA, MO, AR, LA, ND, SD, NE, KS, OK, TX, MB) return False so
    Cordilleran and Great Plains species cannot pass Eastern NA species-residency discovery.
    """
    up_state = _clean_upper_state_province(state_province)
    lon_val: Optional[float] = None
    if lon is not None:
        try:
            cand_lon = float(lon)
            if math.isfinite(cand_lon):
                lon_val = cand_lon
        except (ValueError, TypeError):
            lon_val = None

    if up_state:
        if up_state in CIS_MISSISSIPPI_STATES_PROVINCES:
            # Guard against grossly transposed western coordinates inside an eastern state label
            if lon_val is not None and lon_val < -92.5 and up_state not in {"ONTARIO", "ON"}:
                return False
            return True
        return False

    if lon_val is not None and float(cis_longitude_threshold) <= lon_val <= -50.0:
        if lat is not None:
            try:
                lat_f = float(lat)
                if math.isfinite(lat_f) and lat_f < 24.0:
                    return False
            except (ValueError, TypeError):
                pass
        return True
    return False


def evaluate_eastern_na_cohort_filter(
    df: "pd.DataFrame",
    standardize_taxon_fn: Any,
    *,
    western_longitude_threshold: float = WESTERN_LONGITUDE_THRESHOLD,
    transcontinental_longitude_threshold: Optional[float] = None,
    cis_longitude_threshold: float = DEFAULT_CIS_MISSISSIPPI_LON,
    min_cis_records: int = 10,
    min_cis_fraction: float = DEFAULT_MIN_CIS_MISSISSIPPI_FRAC,
    strict_species_whitelist: bool = False,
) -> Tuple["pd.Series", Dict[str, Any]]:
    """
    Two-stage Cis-Mississippi Eastern North America + Ecoregion-Bounded cohort filter:
      1. Species Gate: Retains genuinely cis-Mississippi Eastern Packera species from
         `resolve_eastern_resident_whitelist()` (YAML + SQLite TDWG/NatureServe) or
         meeting `cis_count >= eff_min_cis` and `cis_fraction >= min_cis_fraction`.
      2. Voucher Gate: Retains cis-Mississippi vouchers plus trans-Mississippi disjuncts
         up to `western_longitude_threshold` (default `-100.0W`), screens transcontinental
         boreal taxa by EPA/CEC Level III ecoregion via `evaluate_transcontinental_ecoregion_gate()`,
         and rescues verified nomenclatural types of retained Eastern species.
    """
    if df is None or df.empty:
        empty_mask = pd.Series(dtype=bool, index=df.index if df is not None else None)
        return empty_mask, {
            "retained_species": [],
            "dropped_species": [],
            "species_breakdown": {},
            "biogeographic_boundary_status": {},
        }

    from scripts.core.step01_biogeography_utils import (
        _canonical_taxon_equiv_key,
        evaluate_transcontinental_ecoregion_gate,
        resolve_eastern_resident_whitelist,
        resolve_transcontinental_taxa,
    )
    _ = transcontinental_longitude_threshold
    whitelist_keys = {_canonical_taxon_equiv_key(s): s for s in resolve_eastern_resident_whitelist()}
    transcontinental_keys = {_canonical_taxon_equiv_key(s) for s in resolve_transcontinental_taxa()}

    def _resolve_sp(row: Mapping[str, Any]) -> str:
        for col in ("species_standardized", "scientificName", "species_raw"):
            val = row.get(col)
            if val is not None and str(val).strip() and str(val).strip().lower() not in {"nan", "none"}:
                std = standardize_taxon_fn(val) if callable(standardize_taxon_fn) else None
                if std:
                    return str(std).strip()
        raw_fb = str(row.get("scientificName") or row.get("species_raw") or "").strip()
        return raw_fb

    lat_col = "decimalLatitude" if "decimalLatitude" in df.columns else ("latitude" if "latitude" in df.columns else None)
    lon_col = "decimalLongitude" if "decimalLongitude" in df.columns else ("longitude" if "longitude" in df.columns else None)

    n_rows = len(df)
    sp_list: List[str] = []
    is_type_list: List[bool] = []
    is_cis_list: List[bool] = []
    is_west_excl_list: List[bool] = []
    is_transcont_west_list: List[bool] = []
    bio_status_list: List[str] = []

    records = df.to_dict(orient="records")
    for rec in records:
        sp = _resolve_sp(rec)
        sp_list.append(sp)
        if "_is_type_rec" in rec and rec["_is_type_rec"] is not None and not pd.isna(rec["_is_type_rec"]):
            is_t = bool(rec["_is_type_rec"])
        else:
            is_t, _ = _is_verified_type_status(rec.get("type_status") or rec.get("_raw_type_status"))
        is_type_list.append(is_t)

        st_val = rec.get("stateProvince")
        lat_val = rec.get(lat_col) if lat_col else None
        lon_val = rec.get(lon_col) if lon_col else None

        cis_ok = is_cis_mississippi_record(
            st_val, lat=lat_val, lon=lon_val, cis_longitude_threshold=cis_longitude_threshold
        )
        is_cis_list.append(cis_ok)

        west_excl = is_excluded_western_region(
            st_val,
            lat=lat_val,
            lon=lon_val,
            western_longitude_threshold=western_longitude_threshold,
            exclude_western=True,
            is_type_specimen=is_t,
        )
        is_west_excl_list.append(west_excl)

        tc_west = False
        sp_key = _canonical_taxon_equiv_key(sp)
        if not is_t and sp_key in transcontinental_keys:
            tc_west, bio_status = evaluate_transcontinental_ecoregion_gate(
                st_val,
                lat=lat_val,
                lon=lon_val,
                regional_group=rec.get("regional_group"),
                is_type_specimen=is_t,
                western_longitude_threshold=western_longitude_threshold,
            )
        elif is_t:
            bio_status = "type_exempt"
        elif west_excl:
            bio_status = "excluded_plains_cordilleran"
        elif cis_ok:
            bio_status = "core_eastern"
        else:
            bio_status = "disjunct_western_extension"
        is_transcont_west_list.append(tc_west)
        bio_status_list.append(bio_status)

    eff_min_cis = int(min_cis_records) if n_rows >= 50 else min(3, max(1, int(min_cis_records)))

    sp_totals: Dict[str, int] = {}
    sp_non_type_totals: Dict[str, int] = {}
    sp_cis_counts: Dict[str, int] = {}
    for sp, is_t, is_cis in zip(sp_list, is_type_list, is_cis_list):
        if not sp:
            continue
        sp_totals[sp] = sp_totals.get(sp, 0) + 1
        if not is_t:
            sp_non_type_totals[sp] = sp_non_type_totals.get(sp, 0) + 1
            if is_cis:
                sp_cis_counts[sp] = sp_cis_counts.get(sp, 0) + 1

    retained_species: Set[str] = set()
    for sp, tot_n in sp_totals.items():
        cis_n = sp_cis_counts.get(sp, 0)
        nt_n = max(1, sp_non_type_totals.get(sp, 0))
        cis_frac = float(cis_n) / float(nt_n)
        sp_key = _canonical_taxon_equiv_key(sp)
        in_whitelist = sp_key in whitelist_keys
        if strict_species_whitelist:
            if in_whitelist:
                retained_species.add(sp)
        else:
            if in_whitelist and (cis_n >= 1 or tot_n == sp_totals[sp] - sp_non_type_totals.get(sp, 0) or n_rows < 50):
                if cis_n >= 1 or sp_non_type_totals.get(sp, 0) == 0:
                    retained_species.add(sp)
                elif cis_frac >= float(min_cis_fraction):
                    retained_species.add(sp)
                else:
                    if sp_key not in transcontinental_keys:
                        retained_species.add(sp)
            elif cis_n >= eff_min_cis and cis_frac >= float(min_cis_fraction):
                retained_species.add(sp)

    keep_flags: List[bool] = []
    species_breakdown: Dict[str, Dict[str, Any]] = {}
    for sp in sorted(sp_totals.keys()):
        species_breakdown[sp] = {
            "total": sp_totals[sp],
            "cis_mississippi": sp_cis_counts.get(sp, 0),
            "retained": 0,
            "disjunct_retained": 0,
            "dropped": 0,
            "resident_accepted": sp in retained_species,
        }

    bio_status_counts: Dict[str, int] = {}
    for sp, is_t, is_cis, west_excl, tc_west, b_st in zip(
        sp_list, is_type_list, is_cis_list, is_west_excl_list, is_transcont_west_list, bio_status_list
    ):
        if not sp or sp not in retained_species:
            keep_flags.append(False)
            bio_status_counts["excluded_non_resident_taxon"] = bio_status_counts.get("excluded_non_resident_taxon", 0) + 1
            if sp in species_breakdown:
                species_breakdown[sp]["dropped"] += 1
            continue
        if is_t:
            keep_flags.append(True)
            bio_status_counts["type_exempt"] = bio_status_counts.get("type_exempt", 0) + 1
            species_breakdown[sp]["retained"] += 1
            continue
        if west_excl or tc_west:
            keep_flags.append(False)
            bio_status_counts["excluded_plains_cordilleran"] = bio_status_counts.get("excluded_plains_cordilleran", 0) + 1
            species_breakdown[sp]["dropped"] += 1
            continue
        keep_flags.append(True)
        bio_status_counts[b_st] = bio_status_counts.get(b_st, 0) + 1
        species_breakdown[sp]["retained"] += 1
        if not is_cis:
            species_breakdown[sp]["disjunct_retained"] += 1

    keep_series = pd.Series(keep_flags, index=df.index, dtype=bool)
    dropped_species = sorted(set(sp_totals.keys()) - retained_species)
    stats = {
        "initial_rows": n_rows,
        "retained_rows": int(keep_series.sum()),
        "dropped_rows": int(n_rows - int(keep_series.sum())),
        "retained_species": sorted(retained_species),
        "dropped_species": dropped_species,
        "species_breakdown": species_breakdown,
        "biogeographic_boundary_status": bio_status_counts,
    }
    return keep_series, stats


def select_cascading_outgroup_anchors(
    df: pd.DataFrame,
    outgroup_taxa: List[str],
    max_anchors: int = 30,
    workspace_dir: Optional[Path] = None,
    focal_df: Optional[pd.DataFrame] = None,
) -> Tuple[List[pd.DataFrame], Dict[str, int], Dict[str, Dict[str, int]]]:
    """
    Selects up to `max_anchors` reference sheets per outgroup taxon by cascading across
    nomenclatural types, range-clean standing, and determiner authority tiers
    (`Tier_1_Gold` -> `Tier_2_Silver` -> `Tier_3_Bronze`), with geographic round-robin
    stratification across `(stateProvince, regional_group)` within each stratum
    prioritizing regions overlapping `focal_df`.
    """
    outgroup_dfs: List[pd.DataFrame] = []
    outgroup_counts: Dict[str, int] = {}
    outgroup_tier_breakdown: Dict[str, Dict[str, int]] = {}

    if df is None or df.empty or not outgroup_taxa or max_anchors <= 0:
        for tx in outgroup_taxa or []:
            outgroup_counts[tx] = 0
            outgroup_tier_breakdown[tx] = {}
        return outgroup_dfs, outgroup_counts, outgroup_tier_breakdown

    tier_rank = {"Tier_1_Gold": 3, "Tier_2_Silver": 2, "Tier_3_Bronze": 1}
    ws_root = Path(workspace_dir).resolve() if workspace_dir is not None else PROJECT_ROOT

    def _norm_tok(v: Any) -> str:
        if v is None or (pd is not None and pd.isna(v)):
            return ""
        return str(v).strip().lower()

    focal_states: Set[str] = set()
    focal_regions: Set[str] = set()
    if focal_df is not None and not focal_df.empty:
        if "stateProvince" in focal_df.columns:
            focal_states = {_norm_tok(v) for v in focal_df["stateProvince"].tolist()} - {""}
        if "regional_group" in focal_df.columns:
            focal_regions = {_norm_tok(v) for v in focal_df["regional_group"].tolist()} - {""}

    for tx in outgroup_taxa:
        tx_df = df[df["species_standardized"] == tx].copy()
        if tx_df.empty:
            outgroup_counts[tx] = 0
            outgroup_tier_breakdown[tx] = {}
            continue

        def _has_local(row: pd.Series) -> int:
            img_p = str(row.get("image_path", "") if not pd.isna(row.get("image_path", "")) else "").strip()
            cat_str = str(row.get("catalogNumber", "") if not pd.isna(row.get("catalogNumber", "")) else "").strip()
            cat_rel = f"data/raw_vouchers/{cat_str}.webp" if cat_str else ""
            if img_p and resolve_voucher_image_path(img_p, workspace_dir=ws_root).is_file():
                return 1
            if cat_rel and resolve_voucher_image_path(cat_rel, workspace_dir=ws_root).is_file():
                return 1
            return 0

        def _has_url_or_path(row: pd.Series) -> bool:
            for c in ("image_path", "_image_url", "image_url", "identifier", "good_quality_AccessURI"):
                v = row.get(c, None)
                if v is not None and not pd.isna(v) and str(v).strip():
                    return True
            return False

        def _is_prim(row: pd.Series) -> int:
            if "is_primary_duplicate" not in row.index or pd.isna(row.get("is_primary_duplicate")):
                return 1
            v = row.get("is_primary_duplicate")
            return 1 if (v is True or str(v).strip().lower() in {"true", "1", "t", "yes"}) else 0

        tx_df["_has_local"] = tx_df.apply(_has_local, axis=1)
        tx_df["_has_valid"] = (tx_df["_has_local"] > 0) | tx_df.apply(_has_url_or_path, axis=1)
        tx_df = tx_df[tx_df["_has_valid"]].copy()
        if tx_df.empty:
            outgroup_counts[tx] = 0
            outgroup_tier_breakdown[tx] = {}
            continue

        ts_col = "type_status" if "type_status" in tx_df.columns else ("typeStatus" if "typeStatus" in tx_df.columns else None)
        tx_df["_type_rank"] = (
            tx_df[ts_col].map(lambda t: classify_nomenclatural_type_rank(t)[0])
            if ts_col
            else 0
        )
        tx_df["_type_subrank"] = (
            tx_df[ts_col].map(classify_nomenclatural_type_subrank)
            if ts_col
            else 0.0
        )
        tx_df["_tier_rank"] = (
            tx_df["determiner_tier"].astype(str).str.strip().map(lambda t: tier_rank.get(t, 0))
            if "determiner_tier" in tx_df.columns
            else 0
        )
        tx_df["_is_prim"] = tx_df.apply(_is_prim, axis=1)
        tx_df["_range_clean_rank"] = [
            0 if (_is_extralimital_status(r) or _as_bool_flag(r.get("exsiccatae_taxonomic_conflict"))) else 1
            for _, r in tx_df.iterrows()
        ]
        tx_df["_img_res"] = tx_df["image_path"].map(
            lambda p: get_image_resolution_score(str(p) if not pd.isna(p) else "", workspace_dir=ws_root)
        ) if "image_path" in tx_df.columns else 0.0

        tx_df = tx_df.sort_values(
            by=["_type_rank", "_type_subrank", "_range_clean_rank", "_tier_rank", "_is_prim", "_has_local", "_img_res"],
            ascending=[False, False, False, False, False, False, False],
            kind="mergesort",
        )

        # Within each (type_rank, range_clean_rank, tier_rank) cascade stratum, round-robin
        # across (stateProvince, regional_group) buckets with focal-overlapping buckets first
        selected_positions: List[int] = []
        cap_n = int(max_anchors)
        strata_groups: Dict[Tuple[int, int, int], List[int]] = {}
        strata_order: List[Tuple[int, int, int]] = []
        for pos_i, (_, row_i) in enumerate(tx_df.iterrows()):
            s_key = (int(row_i["_type_rank"]), int(row_i["_range_clean_rank"]), int(row_i["_tier_rank"]))
            if s_key not in strata_groups:
                strata_groups[s_key] = []
                strata_order.append(s_key)
            strata_groups[s_key].append(pos_i)

        for s_key in strata_order:
            if len(selected_positions) >= cap_n:
                break
            stratum_positions = strata_groups[s_key]
            buckets: Dict[Tuple[str, str], List[int]] = {}
            bucket_order: List[Tuple[str, str]] = []
            for p_i in stratum_positions:
                r_p = tx_df.iloc[p_i]
                b_key = (_norm_tok(r_p.get("stateProvince")), _norm_tok(r_p.get("regional_group")))
                if b_key not in buckets:
                    buckets[b_key] = []
                    bucket_order.append(b_key)
                buckets[b_key].append(p_i)
            if focal_states or focal_regions:
                focal_buckets = [
                    bk for bk in bucket_order
                    if (bk[0] and bk[0] in focal_states) or (bk[1] and bk[1] in focal_regions)
                ]
                other_buckets = [bk for bk in bucket_order if bk not in focal_buckets]
                bucket_order = focal_buckets + other_buckets
            wave_idx = 0
            while len(selected_positions) < cap_n:
                added_wave = False
                for bk in bucket_order:
                    b_list = buckets[bk]
                    if wave_idx < len(b_list):
                        selected_positions.append(b_list[wave_idx])
                        added_wave = True
                        if len(selected_positions) >= cap_n:
                            break
                if not added_wave:
                    break
                wave_idx += 1

        tx_head = tx_df.iloc[selected_positions]
        n_type_sel = int((tx_head["_type_rank"] >= 1).sum())
        n_prim_type_sel = int((tx_head["_type_rank"] == 2).sum())
        n_sec_type_sel = int((tx_head["_type_rank"] == 1).sum())
        tx_sel = tx_head.drop(
            columns=["_has_local", "_has_valid", "_type_rank", "_type_subrank", "_range_clean_rank", "_tier_rank", "_is_prim", "_img_res"]
        )
        outgroup_dfs.append(tx_sel)
        outgroup_counts[tx] = len(tx_sel)
        if "determiner_tier" in tx_sel.columns:
            breakdown = {
                str(k): int(v) for k, v in tx_sel["determiner_tier"].astype(str).value_counts().items()
            }
        else:
            breakdown = {"Unspecified": len(tx_sel)}
        breakdown["Type_Specimens"] = n_type_sel
        if n_prim_type_sel > 0:
            breakdown["Tier_0_Primary_Type"] = n_prim_type_sel
        if n_sec_type_sel > 0:
            breakdown["Tier_0_Secondary_Type"] = n_sec_type_sel
        outgroup_tier_breakdown[tx] = breakdown

    return outgroup_dfs, outgroup_counts, outgroup_tier_breakdown
