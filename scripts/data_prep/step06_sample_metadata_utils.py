#!/usr/bin/env python3
"""
Module: scripts/data_prep/step06_sample_metadata_utils.py
Project: Plantmatics - Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  Sample-ID parsing, staged genomic asset discovery, field/user CSV metadata
  parsing, SNPsaurus stats/key readers, and G3 sample coordinate hygiene.
"""

from __future__ import annotations

import csv, math, os, re, unicodedata
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from scripts.core.coord_hygiene_utils import (
    _CA_PROV_ABBR_TO_FULL, _US_STATE_ABBR_TO_FULL, _normalize_state,
    evaluate_coordinate_centroid_quarantine,
)
from scripts.data_prep.step06_synonym_utils import canonical_taxon_equiv_key, normalize_step06_taxon_name


def _resolve_target_filter(attr: str, env_key: str, default: str) -> str:
    try:
        from scripts.core import config as _cfg
        fn = getattr(_cfg, attr, None)
        if callable(fn) and (val := str(fn()).strip()):
            return val.lower()
    except Exception:
        pass
    return (os.environ.get(f"PLANTMATICS_{env_key}") or os.environ.get(f"PACKERA_{env_key}") or default).strip().lower()


PAC_CODE_TO_TAXON: Dict[str, Tuple[str, str]] = {
    "TOM": ("Packera tomentosa", "Packera dubia"), "DUB": ("Packera dubia", "Packera dubia"),
    "SER": ("Packera serpenticola", "Packera serpenticola"), "SCH": ("Packera schweinitziana", "Packera schweinitziana"),
    "AUR": ("Packera aurea", "Packera aurea"), "NAUR": ("Packera aurea", "Packera aurea"),
    "PAU": ("Packera paupercula", "Packera paupercula"), "MIL": ("Packera millefolium", "Packera millefolium"),
    "ANT": ("Packera antennariifolia", "Packera antennariifolia"), "OBO": ("Packera obovata", "Packera obovata"),
    "ANO": ("Packera anonyma", "Packera anonyma"), "GLA": ("Packera glabella", "Packera glabella"),
    "PSE": ("Packera pseudaurea", "Packera pseudaurea"), "PLATT": ("Packera plattensis", "Packera plattensis"),
    "NOV": ("Packera sp. nov.", "Packera sp. nov."),
}
OUTGROUP_CODE_TO_TAXON: Dict[str, Tuple[str, str]] = {"EREHIE": ("Erechtites hieraciifolius", "Erechtites hieraciifolius")}
_MONTH_MAP: Dict[str, int] = {
    "jan": 1, "january": 1, "feb": 2, "february": 2, "mar": 3, "march": 3, "apr": 4, "april": 4,
    "may": 5, "jun": 6, "june": 6, "jul": 7, "july": 7, "aug": 8, "august": 8,
    "sep": 9, "sept": 9, "september": 9, "oct": 10, "october": 10, "nov": 11, "november": 11, "dec": 12, "december": 12,
}
_US_CA_STATE_SET = set(_US_STATE_ABBR_TO_FULL.values()) | set(_CA_PROV_ABBR_TO_FULL.values()) | {"usa", "united states", "canada"}

VARIANT_SUFFIXES = (".vcf.gz", ".vcf")
DOSAGE_HINTS = ("dosage", "polyrad", "updog")
TREE_PREFIXES = ("RAxML_bestTree", "RAxML_bipartitions")
TREE_SUFFIXES = (".treefile", ".tre", ".nwk", ".tree")
ALIGNMENT_SUFFIXES = (".phy", ".phylip", ".nex", ".nexus", ".fasta", ".fa", ".fas", ".fna")
_TRANSIENT_DOWNLOAD_SUFFIXES = (".crdownload", ".part", ".tmp")
_RAW_SEQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fq.g", ".fastq.g", ".fastqsanger", ".bam", ".bai", ".sam", ".cram")
_UNNUMBERED_RECORD_TOKENS = frozenset({"sn", "unnumbered", "none", "na", "n/a", "unknown", "null"})
_NON_COLLECTOR_TOKENS = frozenset({"unknown", "anonymous", "anon", "na", "none", "null", "collector", "unattributed"})
_NON_TYPE_STATUS_TOKENS = frozenset({"NONE", "NOT A TYPE", "NOTATYPE", "NON-TYPE", "NONTYPE", "UNSPECIFIED", "NULL", "NAN", "NA", "N/A", "FALSE", "0", "UNKNOWN"})


@dataclass
class StagedGenomicAssets:
    variant_files: List[Path]
    dosage_files: List[Path]
    stats_files: List[Path]
    key_files: List[Path]
    metadata_csvs: List[Path]
    tree_files: List[Path]
    alignment_files: List[Path]
    admixture_q_files: List[Path]
    organellar_files: List[Path] = field(default_factory=list)
    public_reference_files: List[Path] = field(default_factory=list)

    @property
    def has_genetic_payload(self) -> bool:
        return bool(self.variant_files or self.dosage_files or self.tree_files or self.admixture_q_files or self.alignment_files)

    @property
    def has_public_or_organellar_payload(self) -> bool:
        return bool(self.public_reference_files or self.organellar_files)


def _canonical_taxon_equiv_key(taxon: Any) -> str:
    return canonical_taxon_equiv_key(taxon)


def normalize_taxon_name(raw_taxon: Optional[str]) -> Tuple[Optional[str], Optional[str], bool]:
    """Return `(taxon_raw, species_standardized, is_novel_or_hybrid)` preserving legacy synonymy."""
    return normalize_step06_taxon_name(raw_taxon)


def strip_smart_quotes(text: Any) -> str:
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return ""
    s = str(text)
    for ch in ('"', "\u201c", "\u201d", "\u2018", "\u2019", "'"):
        s = s.replace(ch, "")
    return s.strip()


def clean_sample_id(raw_id: Any) -> str:
    s = strip_smart_quotes(raw_id)
    if s.endswith("_sorted"):
        s = s[: -len("_sorted")]
    return s.strip()


def _safe_text(val: Any) -> str:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return ""
    s = strip_smart_quotes(val)
    return "" if s.lower() in {"nan", "none", "null"} else s


def _extract_collector_surname(collector: Any) -> str:
    if not (raw := _safe_text(collector)):
        return ""
    if not (primary := re.split(r"\s*(?:;|&|\||\band\b|\bet\s+al\.?)\s*", raw, maxsplit=1, flags=re.IGNORECASE)[0].strip()):
        return ""
    ascii_name = unicodedata.normalize("NFKD", primary).encode("ascii", "ignore").decode("ascii")
    ascii_name = re.sub(r"(?:,\s*|\s+)\b(?:jr|sr|ii|iii|iv)\.?\s*$", "", ascii_name, flags=re.IGNORECASE).strip()
    if "," in ascii_name:
        before, _ = ascii_name.split(",", 1)
        tokens = before.strip().split()
        raw_surname = tokens[-1] if (len(tokens) >= 2 and (tokens[0].endswith(".") or len(tokens[0].rstrip(".")) == 1)) else before
    else:
        tokens = ascii_name.strip().split()
        non_init = [t for t in tokens if len(t.rstrip(".")) > 1]
        raw_surname = non_init[-1] if non_init else (tokens[-1] if tokens else "")
    surname = re.sub(r"[^a-z0-9]+", "", raw_surname.lower())
    return "" if (not surname or surname in _NON_COLLECTOR_TOKENS) else surname


def _normalize_record_number(record_number: Any) -> str:
    if not (raw := _safe_text(record_number)):
        return ""
    low = re.sub(r"^(?:no\.?|#)\s*", "", raw, flags=re.IGNORECASE).strip().lower()
    compact = re.sub(r"[\s._-]+", "", low)
    if not compact or low in {"s.n.", "s.n", "s. n.", "s. n"} or compact in _UNNUMBERED_RECORD_TOKENS or bool(re.fullmatch(r"0+", compact)):
        return ""
    return re.sub(r"^0+(?=\d)", "", re.sub(r"[^a-z0-9]+", "", low))


def _normalize_event_date_part(event_date: Any) -> str:
    if not (raw := _safe_text(event_date)):
        return ""
    if m := re.search(r"(\d{4})[-/](\d{1,2})[-/](\d{1,2})", raw):
        return f"{m.group(1)}-{int(m.group(2)):02d}-{int(m.group(3)):02d}"
    if (m := re.search(r"\b(\d{1,2})(?:st|nd|rd|th)?\s+([A-Za-z]+)\.?\s*,?\s*(1[7-9]\d{2}|20\d{2})\b", raw, re.I)) and (mon := _MONTH_MAP.get(m.group(2).lower())):
        return f"{m.group(3)}-{mon:02d}-{int(m.group(1)):02d}"
    if (m := re.search(r"\b([A-Za-z]+)\.?\s+(\d{1,2})(?:st|nd|rd|th)?\s*,?\s*(1[7-9]\d{2}|20\d{2})\b", raw, re.I)) and (mon := _MONTH_MAP.get(m.group(1).lower())):
        return f"{m.group(3)}-{mon:02d}-{int(m.group(2)):02d}"
    if m := re.search(r"(\d{4})[-/](\d{1,2})\b", raw):
        return f"{m.group(1)}-{int(m.group(2)):02d}"
    return m.group(1) if (m := re.search(r"\b(1[7-9]\d{2}|20\d{2})\b", raw)) else ""


def _normalize_exsiccatae_key(collector: str, record_number: str, event_date: str = "") -> str:
    """Build deterministic `<surname>|<number>|<date>` key for numbered exsiccatae; return `""` when unnumbered."""
    surname, norm_num = _extract_collector_surname(collector), _normalize_record_number(record_number)
    if not surname or not norm_num:
        return ""
    date_part = _normalize_event_date_part(event_date)
    return f"{surname}|{norm_num}|{date_part}" if date_part else f"{surname}|{norm_num}"


def parse_sample_code(sample_id: str) -> Dict[str, Optional[str]]:
    """Parse SNPsaurus run_5563 sample IDs and hyphenated Project_ID codes cleanly."""
    cid = clean_sample_id(sample_id)
    empty = {"taxon_code": None, "taxon_raw": None, "species_standardized": None, "population_id": None, "replicate_pair_id": None, "replicate_suffix": None}
    if not cid:
        return empty
    if m := re.match(r"^PAC(NAUR|[A-Z]{3,5})\d*-([A-Z0-9]+)(?:-([A-Z0-9-]+))?$", cid, re.I):
        code, pop = m.group(1).upper(), m.group(2).upper()
        t_raw, sp_std = PAC_CODE_TO_TAXON.get(code, (None, None))
        return {"taxon_code": code, "taxon_raw": t_raw, "species_standardized": sp_std, "population_id": pop, "replicate_pair_id": f"PAC{code}{pop}", "replicate_suffix": None}
    if m := re.match(r"^PAC(NAUR|PLATT|[A-Z]{3})([A-Z0-9]+?)([MP])$", cid, re.I):
        code, pop, rep = m.group(1).upper(), m.group(2).upper(), m.group(3).upper()
        t_raw, sp_std = PAC_CODE_TO_TAXON.get(code, (None, None))
        return {"taxon_code": code, "taxon_raw": t_raw, "species_standardized": sp_std, "population_id": pop, "replicate_pair_id": f"PAC{code}{pop}", "replicate_suffix": rep}
    if m := re.match(r"^PAC(NAUR|PLATT|[A-Z]{3})([A-Z0-9]+)$", cid, re.I):
        code, pop = m.group(1).upper(), m.group(2).upper()
        t_raw, sp_std = PAC_CODE_TO_TAXON.get(code, (None, None))
        return {"taxon_code": code, "taxon_raw": t_raw, "species_standardized": sp_std, "population_id": pop, "replicate_pair_id": cid, "replicate_suffix": None}
    m_og_hy = re.match(r"^EREHIE\d*[-_]([A-Z0-9]+?)(?:[-_]?([MP])|[-_]([A-Z0-9]+))?$", cid, re.I)
    m_og_flat = None if m_og_hy else re.match(r"^EREHIE([A-Z0-9]+?)([MP])?$", cid, re.I)
    if m_og_hy or m_og_flat:
        tax_raw, sp_std = OUTGROUP_CODE_TO_TAXON["EREHIE"]
        pop = (m_og_hy.group(1) if m_og_hy else m_og_flat.group(1)).upper()
        sfx = ((m_og_hy.group(2) or m_og_hy.group(3) or "") if m_og_hy else (m_og_flat.group(2) or "")).upper() or None
        return {"taxon_code": "EREHIE", "taxon_raw": tax_raw, "species_standardized": sp_std, "population_id": pop, "replicate_pair_id": f"EREHIE{pop}", "replicate_suffix": sfx}
    return {**empty, "replicate_pair_id": cid}


def discover_genomic_inputs(genomics_dir: Path) -> StagedGenomicAssets:
    """Discover staged variant, dosage, metadata, alignment, organellar, public reference, admixture, and tree inputs."""
    assets = StagedGenomicAssets([], [], [], [], [], [], [], [], [], [])
    if not genomics_dir.exists() or not genomics_dir.is_dir():
        return assets
    skip_names = {"parsed_genomic_samples.csv", "allele_dosage_matrix.tsv", "genomic_manifest.json", "hdplot_locus_diagnostics.tsv", "locus_weights.tsv"}
    for p in sorted(genomics_dir.rglob("*")):
        if not p.is_file() or any(part.startswith(".") for part in p.relative_to(genomics_dir).parts):
            continue
        name, low = p.name, p.name.lower()
        if low.endswith(_TRANSIENT_DOWNLOAD_SUFFIXES) or low.endswith(_RAW_SEQ_SUFFIXES) or low in skip_names:
            continue
        if low.startswith("public_packera_"):
            assets.public_reference_files.append(p)
        elif low.startswith("skimmed_") or low in ("bain_jansen_sequences.fasta", "cytonuclear_discordance.csv", "cytonuclear_discordance.tsv"):
            assets.organellar_files.append(p)
        elif low.endswith(VARIANT_SUFFIXES):
            assets.variant_files.append(p)
        elif low.endswith((".tsv", ".csv")) and any(h in low for h in DOSAGE_HINTS):
            assets.dosage_files.append(p)
        elif low.endswith(".csv") and any(h in low for h in ("field", "specimen", "user_genetic_metadata", "metadata")):
            assets.metadata_csvs.append(p)
        elif low.endswith(".txt") and "stats" in low:
            assets.stats_files.append(p)
        elif low.endswith(".txt") and "key" in low:
            assets.key_files.append(p)
        elif name.startswith(TREE_PREFIXES) or low.endswith(TREE_SUFFIXES):
            assets.tree_files.append(p)
        elif low.endswith(".q"):
            assets.admixture_q_files.append(p)
        elif low.endswith(ALIGNMENT_SUFFIXES) and not low.startswith(("vcfref_", "ref_packera")):
            assets.alignment_files.append(p)
    return assets


def _read_text_rows(path: Path) -> pd.DataFrame:
    for enc in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            return pd.read_csv(path, dtype=str, encoding=enc, keep_default_na=False)
        except UnicodeDecodeError:
            continue
    return pd.read_csv(path, dtype=str, encoding="latin-1", errors="replace", keep_default_na=False)


def _parse_lat_lon(val: Any) -> Tuple[Optional[float], Optional[float]]:
    if val is None or not (s := str(val).strip()):
        return None, None
    if len(m := re.findall(r"[-+]?\d+(?:\.\d+)?", s)) >= 2:
        try:
            lat, lon = float(m[0]), float(m[1])
            if -90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0:
                return lat, lon
        except ValueError:
            pass
    return None, None


def sanitize_step06_sample_coordinates(
    lat: Any, lon: Any, state_prov: Any = "", county: Any = "", locality: Any = "", *, is_vouchered: bool = False,
) -> Tuple[Optional[float], Optional[float], bool]:
    """Nullify out-of-bounds, (0,0), US/CA sign-flipped, and unvouchered rooftop coordinates."""
    if lat is None or lon is None or isinstance(lat, bool) or isinstance(lon, bool):
        return None, None, False
    try:
        lat_f, lon_f = float(lat), float(lon)
    except (TypeError, ValueError):
        return None, None, True
    if not (math.isfinite(lat_f) and math.isfinite(lon_f)) or not (-90.0 <= lat_f <= 90.0 and -180.0 <= lon_f <= 180.0):
        return None, None, True
    if abs(lat_f) < 1e-4 and abs(lon_f) < 1e-4:
        return None, None, True
    if _normalize_state(state_prov) in _US_CA_STATE_SET and lat_f > 15.0 and lon_f > 40.0:
        return None, None, True
    if not is_vouchered:
        flagged, reason = evaluate_coordinate_centroid_quarantine(lat_f, lon_f, locality=str(locality or ""), county=str(county or ""), state_prov=str(state_prov or ""))
        if flagged and reason == "herbarium_institution_rooftop":
            return None, None, True
    return lat_f, lon_f, False


def _parse_field_metadata(metadata_csvs: Sequence[Path]) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    """Parse Field Specimen Database CSVs and optional `user_genetic_metadata.csv`."""
    field_rows: List[Dict[str, Any]] = []
    direct_user_meta: Dict[str, Dict[str, Any]] = {}
    target_family = _resolve_target_filter("resolve_target_family", "TARGET_FAMILY", "asteraceae")
    target_genus = _resolve_target_filter("resolve_target_genus", "TARGET_GENUS", "packera")

    for csv_path in metadata_csvs:
        df = _read_text_rows(csv_path)
        if df.empty:
            continue
        cols_low = {c.lower().strip(): c for c in df.columns}
        if (sid_col := cols_low.get("genetic_sample_id") or cols_low.get("sample_id")) is not None:
            ex_col = cols_low.get("exsiccatae_key")
            for _, row in df.iterrows():
                if not (sid := clean_sample_id(row.get(sid_col, ""))):
                    continue
                meta_row = {str(k).strip(): (str(v).strip() if str(v).strip() != "" else None) for k, v in row.items()}
                if ex_col and "exsiccatae_key" not in meta_row:
                    meta_row["exsiccatae_key"] = str(row.get(ex_col, "")).strip() or None
                direct_user_meta[sid] = meta_row
            continue

        fam_col, tax_col = cols_low.get("family"), cols_low.get("taxon") or cols_low.get("scientificname")
        pid_col, latlon_col = cols_low.get("project_id"), cols_low.get("lat_lon")
        lat_col, lon_col = cols_low.get("decimallatitude") or cols_low.get("latitude"), cols_low.get("decimallongitude") or cols_low.get("longitude")
        state_col, county_col = cols_low.get("state") or cols_low.get("stateprovince"), cols_low.get("county")
        date_col, notes_col, ex_col = cols_low.get("date_collected") or cols_low.get("eventdate"), cols_low.get("notes"), cols_low.get("exsiccatae_key")
        rec_by_col = cols_low.get("recordedby") or cols_low.get("collector")
        rec_num_col = cols_low.get("recordnumber") or cols_low.get("record_number") or cols_low.get("collector_number") or cols_low.get("collection_number") or cols_low.get("collector_id")
        alias_col, q_col = cols_low.get("population_alias"), cols_low.get("is_putative_novel_or_hybrid_query")

        for _, row in df.iterrows():
            taxon_str = str(row.get(tax_col, "")).strip() if tax_col else ""
            fam_str = str(row.get(fam_col, "")).strip() if fam_col else ""
            if (fam_str and target_family and fam_str.lower() != target_family) or (target_genus and not taxon_str.lower().startswith(target_genus)):
                continue
            proj_clean = strip_smart_quotes(row.get(pid_col, "")) if pid_col else ""
            parsed_proj = parse_sample_code(proj_clean)
            lat, lon = (None, None)
            if lat_col and lon_col:
                try:
                    lat, lon = float(row.get(lat_col, "")), float(row.get(lon_col, ""))
                except (TypeError, ValueError):
                    lat, lon = (None, None)
            if (lat is None or lon is None) and latlon_col:
                lat, lon = _parse_lat_lon(row.get(latlon_col, ""))
            t_raw, sp_std, is_query = normalize_taxon_name(taxon_str)
            if q_col and str(row.get(q_col, "")).strip().lower() in ("1", "true", "yes"):
                is_query = True
            ev_date_str = str(row.get(date_col, "")).strip() or None if date_col else None
            rec_by = str(row.get(rec_by_col, "")).strip() or None if rec_by_col else None
            rec_num = str(row.get(rec_num_col, "")).strip() or None if rec_num_col else None
            ex_key = (str(row.get(ex_col, "")).strip() or None if ex_col else None) or (_normalize_exsiccatae_key(rec_by, rec_num, ev_date_str or "") or None if (rec_by and rec_num) else None)
            field_rows.append({
                "project_id": proj_clean or None, "taxon_code": parsed_proj["taxon_code"], "population_id": parsed_proj["population_id"],
                "population_alias": str(row.get(alias_col, "")).strip() or None if alias_col else None,
                "taxon_raw": t_raw, "species_standardized": sp_std, "is_putative_novel_or_hybrid_query": is_query,
                "decimalLatitude": lat, "decimalLongitude": lon,
                "stateProvince": str(row.get(state_col, "")).strip() or None if state_col else None,
                "county": str(row.get(county_col, "")).strip() or None if county_col else None,
                "eventDate": ev_date_str, "recordedBy": rec_by, "recordNumber": rec_num, "exsiccatae_key": ex_key,
                "notes": (str(row.get(notes_col, "")).strip() or None) if notes_col else None,
            })
    return field_rows, direct_user_meta


def _parse_stats_file(stats_path: Path) -> Dict[str, Dict[str, Any]]:
    """Parse SNPsaurus `packera_stats.txt` cleanly, returning per-sample QC metrics."""
    out: Dict[str, Dict[str, Any]] = {}
    with stats_path.open("r", encoding="utf-8-sig", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not (sid := clean_sample_id(row.get("sample") or row.get("Sample") or "")):
                continue
            def _num(keys: Sequence[str], as_float: bool = False) -> Any:
                for k in keys:
                    if k in row and (s := str(row[k]).replace("%", "").strip()):
                        try: return float(s) if as_float else int(float(s))
                        except ValueError: return None
                return None
            n_miss, n_00, n_01, n_11 = (_num([k]) or 0 for k in ("# ./. genotypes", "# 0/0 genotypes", "# 0/1 genotypes", "# 1/1 genotypes"))
            n_called, n_total = n_00 + n_01 + n_11, n_00 + n_01 + n_11 + n_miss
            out[sid] = {
                "hq_reads": _num(["HQ reads", "hq_reads"]), "mapped_reads": _num(["HQ reads mapped to reference", "mapped_reads"]),
                "pct_repetitive": _num(["percent repetitive reads", "pct_repetitive"], True),
                "n_called_loci": n_called if n_total > 0 else None,
                "missing_genotype_rate": float(n_miss / n_total) if n_total > 0 else None,
                "observed_het": float(n_01 / n_called) if n_called > 0 else None,
            }
    return out


def _parse_key_samples(key_files: Sequence[Path]) -> List[str]:
    samples: List[str] = []
    for kp in key_files:
        with kp.open("r", encoding="utf-8-sig", errors="replace") as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) >= 2 and "-" in parts[0] and (sid := clean_sample_id(parts[1])) and sid not in samples:
                    samples.append(sid)
    return samples

