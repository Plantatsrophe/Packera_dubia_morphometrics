#!/usr/bin/env python3
"""
Module: scripts/data_prep/06_ingest_genetic_data.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  Independent Step 06 Tier A Genomic Evidence Stream ingestion, HDplot paralog
  filter, Empirical Bayes polyploid dosage estimation, cytotype mixture scoring,
  organellar/cytonuclear ingestion, and 4-tier voucher linker.
  Auto-discovers flat staged files inside `data/genomics/`.

  Invariants:
    - Executes cleanly as a zero-overhead no-op when `data/genomics/` contains no
      supported genomic payloads (empty or `.gitkeep` only).
    - Strictly independent of Step 03/04 morphometrics and Step 05 ecological niche
      modeling (reads only taxonomy/geography metadata from `vouchers` and staged CSVs).
    - Gates high-dropout libraries (`hq_reads < 50000` or `missing_genotype_rate > 0.80`)
      with `geno_quality_pass = 0` for audit reporting while excluding them from
      active ordination and dosage export.
    - Masks excess-heterozygosity paralog loci (`Ho > 0.65` across reference samples)
      without applying global minor-allele-frequency pruning (`MAF >= 0.05` disabled).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sqlite3
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import get_db_connection, upsert_genomic_samples
from scripts.core.logger import setup_logging

logger = setup_logging(name="packera.genomics.ingest")

# Taxon code table used in Fuller / SNPsaurus field and sequencing labels
PAC_CODE_TO_TAXON: Dict[str, Tuple[str, str]] = {
    "TOM": ("Packera tomentosa", "Packera dubia"),
    "DUB": ("Packera dubia", "Packera dubia"),
    "SER": ("Packera serpenticola", "Packera serpenticola"),
    "SCH": ("Packera schweinitziana", "Packera schweinitziana"),
    "AUR": ("Packera aurea", "Packera aurea"),
    "NAUR": ("Packera aurea", "Packera aurea"),
    "PAU": ("Packera paupercula", "Packera paupercula"),
    "MIL": ("Packera millefolium", "Packera millefolium"),
    "ANT": ("Packera antennariifolia", "Packera antennariifolia"),
    "OBO": ("Packera obovata", "Packera obovata"),
    "ANO": ("Packera anonyma", "Packera anonyma"),
    "GLA": ("Packera glabella", "Packera glabella"),
    "PSE": ("Packera pseudaurea", "Packera pseudaurea"),
    "PLATT": ("Packera plattensis", "Packera plattensis"),
    "NOV": ("Packera sp. nov.", "Packera sp. nov."),
}

KNOWN_NOVEL_OR_HYBRID_POPS = {"VAL1", "VAL2", "PAR1"}
KNOWN_NOVEL_OR_HYBRID_SAMPLE_PAIRS = {
    ("MIL", "VAL1"),
    ("MIL", "VAL2"),
    ("NOV", "VAL1"),
    ("NOV", "VAL2"),
    ("PAU", "PAR1"),
    ("NOV", "PAR1"),
}


def _is_novel_or_hybrid_sample_pair(code: Optional[str], pop: Optional[str]) -> bool:
    c = (code or "").upper()
    p = (pop or "").upper()
    if c == "NOV":
        return True
    if not c and p in KNOWN_NOVEL_OR_HYBRID_POPS:
        return True
    return (c, p) in KNOWN_NOVEL_OR_HYBRID_SAMPLE_PAIRS


POPULATION_ALIAS_MAP: Dict[str, Tuple[str, ...]] = {
    "NCS1": ("NCS1", "NCSP"),
    "NCSP": ("NCSP", "NCS1"),
    "TX2": ("TX2", "TX1"),
}

VARIANT_SUFFIXES = (".vcf.gz", ".vcf")
DOSAGE_HINTS = ("dosage", "polyrad", "updog")
TREE_PREFIXES = ("RAxML_bestTree", "RAxML_bipartitions")
TREE_SUFFIXES = (".treefile", ".tre", ".nwk", ".tree")
ALIGNMENT_SUFFIXES = (".phy", ".phylip", ".nex", ".nexus", ".fasta", ".fa", ".fas", ".fna")


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

    @property
    def has_genetic_payload(self) -> bool:
        return bool(
            self.variant_files
            or self.dosage_files
            or self.tree_files
            or self.admixture_q_files
            or self.alignment_files
        )


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


def normalize_taxon_name(raw_taxon: Optional[str]) -> Tuple[Optional[str], Optional[str], bool]:
    """Return `(taxon_raw, species_standardized, is_novel_or_hybrid)` preserving legacy synonymy."""
    if raw_taxon is None:
        return None, None, False
    t_raw = " ".join(str(raw_taxon).strip().split())
    if not t_raw:
        return None, None, False

    low = t_raw.lower()
    is_query = bool(
        re.search(r"\bsp\.\s*nov\b", low)
        or re.search(r"\bcf\.", low)
        or re.search(r"\baff\.", low)
        or re.search(r"\bhybrid\b", low)
        or " x " in low
        or "\u00d7" in low
    )
    if is_query:
        return t_raw, t_raw, True

    if low in ("packera tomentosa", "p. tomentosa", "senecio tomentosus"):
        return t_raw, "Packera dubia", False

    # Lazy import keeps ingest clean even if pygbif is unavailable in light test envs
    try:
        from scripts.core.harvester import standardize_packera_taxon

        std = standardize_packera_taxon(t_raw)
        if std:
            return t_raw, std, False
    except Exception:
        pass
    return t_raw, t_raw, False


def parse_sample_code(sample_id: str) -> Dict[str, Optional[str]]:
    """
    Parse SNPsaurus run_5563 sample IDs (e.g. `PACTOMSCG1M`, `PACPAUVAL1M`, `PACSERNCS1M`)
    and hyphenated Project_ID codes (e.g. `PAC"MIL"01-VAL1-20200704`, `PACNOV01-PAR1-20200702`).
    """
    cid = clean_sample_id(sample_id)
    if not cid:
        return {
            "taxon_code": None,
            "taxon_raw": None,
            "species_standardized": None,
            "population_id": None,
            "replicate_pair_id": None,
            "replicate_suffix": None,
        }

    # Hyphenated Project_ID pattern: PAC[N]CODE[num]-POP[-EXTRA]
    m_hyph = re.match(r"^PAC(NAUR|[A-Z]{3,5})\d*-([A-Z0-9]+)(?:-([A-Z0-9-]+))?$", cid, re.IGNORECASE)
    if m_hyph:
        code = m_hyph.group(1).upper()
        pop = m_hyph.group(2).upper()
        t_raw, sp_std = PAC_CODE_TO_TAXON.get(code, (None, None))
        return {
            "taxon_code": code,
            "taxon_raw": t_raw,
            "species_standardized": sp_std,
            "population_id": pop,
            "replicate_pair_id": f"PAC{code}{pop}",
            "replicate_suffix": None,
        }

    # Unhyphenated SNPsaurus pattern: PAC<3-letter><POP><M|P>
    m_snp = re.match(r"^PAC(NAUR|PLATT|[A-Z]{3})([A-Z0-9]+?)([MP])$", cid, re.IGNORECASE)
    if m_snp:
        code = m_snp.group(1).upper()
        pop = m_snp.group(2).upper()
        rep = m_snp.group(3).upper()
        t_raw, sp_std = PAC_CODE_TO_TAXON.get(code, (None, None))
        return {
            "taxon_code": code,
            "taxon_raw": t_raw,
            "species_standardized": sp_std,
            "population_id": pop,
            "replicate_pair_id": f"PAC{code}{pop}",
            "replicate_suffix": rep,
        }

    # Fallback: PAC<3-letter><POP> without M/P suffix
    m_nomp = re.match(r"^PAC(NAUR|PLATT|[A-Z]{3})([A-Z0-9]+)$", cid, re.IGNORECASE)
    if m_nomp:
        code = m_nomp.group(1).upper()
        pop = m_nomp.group(2).upper()
        t_raw, sp_std = PAC_CODE_TO_TAXON.get(code, (None, None))
        return {
            "taxon_code": code,
            "taxon_raw": t_raw,
            "species_standardized": sp_std,
            "population_id": pop,
            "replicate_pair_id": cid,
            "replicate_suffix": None,
        }

    return {
        "taxon_code": None,
        "taxon_raw": None,
        "species_standardized": None,
        "population_id": None,
        "replicate_pair_id": cid,
        "replicate_suffix": None,
    }


_TRANSIENT_DOWNLOAD_SUFFIXES = (".crdownload", ".part", ".tmp")
_RAW_SEQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fastqsanger", ".bam", ".bai", ".sam", ".cram")


def discover_genomic_inputs(genomics_dir: Path) -> StagedGenomicAssets:
    """Discover flat staged variant, dosage, metadata, alignment, organellar, admixture, and phylogeny inputs."""
    assets = StagedGenomicAssets([], [], [], [], [], [], [], [], [])
    if not genomics_dir.exists() or not genomics_dir.is_dir():
        return assets

    for p in sorted(genomics_dir.iterdir()):
        if not p.is_file():
            continue
        name = p.name
        low = name.lower()
        if name.startswith(".") or low.endswith(_TRANSIENT_DOWNLOAD_SUFFIXES) or low.endswith(_RAW_SEQ_SUFFIXES):
            continue
        if low in ("parsed_genomic_samples.csv", "allele_dosage_matrix.tsv", "genomic_manifest.json", "hdplot_locus_diagnostics.tsv"):
            continue
        if low.startswith("skimmed_") or low in ("bain_jansen_sequences.fasta", "cytonuclear_discordance.csv", "cytonuclear_discordance.tsv"):
            assets.organellar_files.append(p)
            continue
        if low.endswith(".vcf.gz") or low.endswith(".vcf"):
            assets.variant_files.append(p)
        elif low.endswith((".tsv", ".csv")) and any(h in low for h in DOSAGE_HINTS):
            assets.dosage_files.append(p)
        elif low.endswith(".csv") and ("field" in low or "specimen" in low or "user_genetic_metadata" in low or "metadata" in low):
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
    if val is None:
        return None, None
    s = str(val).strip()
    if not s:
        return None, None
    m = re.findall(r"[-+]?\d+(?:\.\d+)?", s)
    if len(m) >= 2:
        try:
            lat = float(m[0])
            lon = float(m[1])
            if -90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0:
                return lat, lon
        except ValueError:
            pass
    return None, None


def _parse_field_metadata(metadata_csvs: Sequence[Path]) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    """
    Parse Field Specimen Database CSVs and optional `user_genetic_metadata.csv`.
    Returns `(field_packera_rows, direct_user_metadata_by_sample_id)`.
    """
    field_rows: List[Dict[str, Any]] = []
    direct_user_meta: Dict[str, Dict[str, Any]] = {}

    for csv_path in metadata_csvs:
        df = _read_text_rows(csv_path)
        if df.empty:
            continue
        cols_low = {c.lower().strip(): c for c in df.columns}
        # Explicit user_genetic_metadata schema if sample_id / genetic_sample_id column present
        sid_col = cols_low.get("genetic_sample_id") or cols_low.get("sample_id")
        if sid_col is not None:
            for _, row in df.iterrows():
                sid = clean_sample_id(row.get(sid_col, ""))
                if not sid:
                    continue
                direct_user_meta[sid] = {str(k).strip(): (str(v).strip() if str(v).strip() != "" else None) for k, v in row.items()}
            continue

        fam_col = cols_low.get("family")
        tax_col = cols_low.get("taxon") or cols_low.get("scientificname")
        pid_col = cols_low.get("project_id")
        latlon_col = cols_low.get("lat_lon")
        lat_col = cols_low.get("decimallatitude") or cols_low.get("latitude")
        lon_col = cols_low.get("decimallongitude") or cols_low.get("longitude")
        state_col = cols_low.get("state") or cols_low.get("stateprovince")
        county_col = cols_low.get("county")
        date_col = cols_low.get("date_collected") or cols_low.get("eventdate")
        notes_col = cols_low.get("notes")

        for _, row in df.iterrows():
            taxon_str = str(row.get(tax_col, "")).strip() if tax_col else ""
            fam_str = str(row.get(fam_col, "")).strip() if fam_col else ""
            if fam_str and fam_str.lower() != "asteraceae":
                continue
            if not taxon_str.lower().startswith("packera"):
                continue

            proj_clean = strip_smart_quotes(row.get(pid_col, "")) if pid_col else ""
            parsed_proj = parse_sample_code(proj_clean)
            lat, lon = (None, None)
            if lat_col and lon_col:
                try:
                    lat = float(row.get(lat_col, ""))
                    lon = float(row.get(lon_col, ""))
                except (TypeError, ValueError):
                    lat, lon = (None, None)
            if (lat is None or lon is None) and latlon_col:
                lat, lon = _parse_lat_lon(row.get(latlon_col, ""))

            t_raw, sp_std, is_query = normalize_taxon_name(taxon_str)
            notes_str = str(row.get(notes_col, "")).strip() if notes_col else ""
            field_rows.append(
                {
                    "project_id": proj_clean or None,
                    "taxon_code": parsed_proj["taxon_code"],
                    "population_id": parsed_proj["population_id"],
                    "taxon_raw": t_raw,
                    "species_standardized": sp_std,
                    "is_putative_novel_or_hybrid_query": is_query,
                    "decimalLatitude": lat,
                    "decimalLongitude": lon,
                    "stateProvince": str(row.get(state_col, "")).strip() or None if state_col else None,
                    "county": str(row.get(county_col, "")).strip() or None if county_col else None,
                    "eventDate": str(row.get(date_col, "")).strip() or None if date_col else None,
                    "notes": notes_str or None,
                }
            )

    return field_rows, direct_user_meta


def _parse_stats_file(stats_path: Path) -> Dict[str, Dict[str, Any]]:
    """Parse SNPsaurus `packera_stats.txt` cleanly, returning per-sample QC metrics."""
    out: Dict[str, Dict[str, Any]] = {}
    with stats_path.open("r", encoding="utf-8-sig", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            raw_sample = row.get("sample") or row.get("Sample") or ""
            sid = clean_sample_id(raw_sample)
            if not sid:
                continue

            def _int_val(keys: Sequence[str]) -> Optional[int]:
                for k in keys:
                    if k in row and str(row[k]).strip() != "":
                        try:
                            return int(float(str(row[k]).strip()))
                        except ValueError:
                            return None
                return None

            def _pct_val(keys: Sequence[str]) -> Optional[float]:
                for k in keys:
                    if k in row and str(row[k]).strip() != "":
                        s = str(row[k]).replace("%", "").strip()
                        try:
                            return float(s)
                        except ValueError:
                            return None
                return None

            hq = _int_val(["HQ reads", "hq_reads"])
            mapped = _int_val(["HQ reads mapped to reference", "mapped_reads"])
            pct_rep = _pct_val(["percent repetitive reads", "pct_repetitive"])
            n_miss = _int_val(["# ./. genotypes"]) or 0
            n_00 = _int_val(["# 0/0 genotypes"]) or 0
            n_01 = _int_val(["# 0/1 genotypes"]) or 0
            n_11 = _int_val(["# 1/1 genotypes"]) or 0
            n_called = n_00 + n_01 + n_11
            n_total = n_called + n_miss
            miss_rate = float(n_miss / n_total) if n_total > 0 else None
            obs_het = float(n_01 / n_called) if n_called > 0 else None
            out[sid] = {
                "hq_reads": hq,
                "mapped_reads": mapped,
                "pct_repetitive": pct_rep,
                "n_called_loci": n_called if n_total > 0 else None,
                "missing_genotype_rate": miss_rate,
                "observed_het": obs_het,
            }
    return out


def _parse_key_samples(key_files: Sequence[Path]) -> List[str]:
    samples: List[str] = []
    for kp in key_files:
        with kp.open("r", encoding="utf-8-sig", errors="replace") as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) >= 2 and "-" in parts[0]:
                    sid = clean_sample_id(parts[1])
                    if sid and sid not in samples:
                        samples.append(sid)
    return samples


def _gt_to_dosage(gt_token: str) -> float:
    """Convert diploid or polyploid GT (`0/0`, `0/1`, `1/1`, `0/0/1/1`, `1|0`) to 0..1 dosage."""
    if not gt_token or "." in gt_token:
        return np.nan
    alleles = re.split(r"[/|]", gt_token)
    if not alleles:
        return np.nan
    try:
        vals = [int(a) for a in alleles]
    except ValueError:
        return np.nan
    if not vals:
        return np.nan
    return float(sum(1 for a in vals if a > 0)) / float(len(vals))


def _safe_nonneg_int(tok: str) -> Optional[int]:
    if not tok or tok in (".", "./.", ".|."):
        return None
    try:
        v = int(float(tok))
        return v if v >= 0 else None
    except ValueError:
        return None


def _extract_cell_depths(
    fields: Sequence[str],
    fmt_idx: Mapping[str, int],
    gt_dosage: float,
) -> Tuple[float, float]:
    """Extract `(n_ref, n_alt)` supporting FreeBayes RO/AO, GATK `ref,alt` AD, and BBMap single-int AD."""
    ro_i = fmt_idx.get("RO", -1)
    ao_i = fmt_idx.get("AO", -1)
    if ro_i >= 0 and ao_i >= 0 and ro_i < len(fields) and ao_i < len(fields):
        ro_v = _safe_nonneg_int(fields[ro_i])
        ao_vals = [_safe_nonneg_int(t) for t in fields[ao_i].split(",")]
        ao_clean = [x for x in ao_vals if x is not None]
        if ro_v is not None and ao_clean:
            return float(ro_v), float(sum(ao_clean))

    ad_i = fmt_idx.get("AD", -1)
    dp_i = fmt_idx.get("DP", -1)
    dp_v = _safe_nonneg_int(fields[dp_i]) if (0 <= dp_i < len(fields)) else None
    if 0 <= ad_i < len(fields):
        ad_tok = fields[ad_i].strip()
        if ad_tok and ad_tok != ".":
            if "," in ad_tok:
                ad_parts = [_safe_nonneg_int(t) for t in ad_tok.split(",")]
                if ad_parts[0] is not None and any(x is not None for x in ad_parts[1:]):
                    n_ref = float(ad_parts[0])
                    n_alt = float(sum(x for x in ad_parts[1:] if x is not None))
                    return n_ref, n_alt
            else:
                n_alt_i = _safe_nonneg_int(ad_tok)
                if n_alt_i is not None:
                    # BBMap v38.69 records a single integer in AD equal to alt allele read count
                    n_ref_f = float(max(0, dp_v - n_alt_i)) if dp_v is not None else 0.0
                    return n_ref_f, float(n_alt_i)

    if dp_v is not None and dp_v > 0 and np.isfinite(gt_dosage):
        n_alt_f = float(round(dp_v * float(gt_dosage)))
        return float(max(0.0, float(dp_v) - n_alt_f)), n_alt_f
    return np.nan, np.nan


def _parse_vcf_with_depths(
    vcf_path: Path,
) -> Tuple[List[str], List[str], np.ndarray, np.ndarray, np.ndarray]:
    """Stream VCF / `.vcf.gz` into `(sample_ids, locus_ids, dosage_mat, ref_depth_mat, alt_depth_mat)`."""
    opener = gzip.open if vcf_path.name.lower().endswith(".gz") else open
    sample_ids: List[str] = []
    locus_ids: List[str] = []
    rows_loci: List[np.ndarray] = []
    rows_ref: List[np.ndarray] = []
    rows_alt: List[np.ndarray] = []

    with opener(vcf_path, "rt", encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.rstrip("\r\n").split("\t")
                sample_ids = [clean_sample_id(s) for s in parts[9:]]
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 10 or not sample_ids:
                continue
            chrom, pos, vid = parts[0], parts[1], parts[2]
            locus_id = vid if (vid and vid != ".") else f"{chrom}_{pos}"
            if locus_ids and locus_id == locus_ids[-1]:
                locus_id = f"{chrom}_{pos}_{len(locus_ids)}"
            fmt = parts[8].split(":") if len(parts) > 8 else ["GT"]
            fmt_idx = {tag: i for i, tag in enumerate(fmt)}
            gt_idx = fmt_idx.get("GT", 0)
            dos_idx = fmt_idx.get("DS", -1)

            vec = np.full(len(sample_ids), np.nan, dtype=np.float32)
            rvec = np.full(len(sample_ids), np.nan, dtype=np.float32)
            avec = np.full(len(sample_ids), np.nan, dtype=np.float32)
            for idx, cell in enumerate(parts[9 : 9 + len(sample_ids)]):
                fields = cell.split(":")
                gt_tok = fields[gt_idx] if gt_idx < len(fields) else "."
                gt_d = _gt_to_dosage(gt_tok)
                if 0 <= dos_idx < len(fields) and fields[dos_idx] not in ("", "."):
                    try:
                        dv = float(fields[dos_idx])
                        vec[idx] = np.float32(dv if dv <= 1.0 else dv / 2.0)
                    except ValueError:
                        vec[idx] = np.float32(gt_d)
                else:
                    vec[idx] = np.float32(gt_d)
                if np.isfinite(vec[idx]):
                    nr, na = _extract_cell_depths(fields, fmt_idx, float(vec[idx]))
                    rvec[idx] = np.float32(nr)
                    avec[idx] = np.float32(na)

            locus_ids.append(locus_id)
            rows_loci.append(vec)
            rows_ref.append(rvec)
            rows_alt.append(avec)

    if not sample_ids or not rows_loci:
        empty = np.empty((len(sample_ids), 0), dtype=np.float32)
        return sample_ids, [], empty, empty.copy(), empty.copy()
    return (
        sample_ids,
        locus_ids,
        np.vstack(rows_loci).T.astype(np.float32),
        np.vstack(rows_ref).T.astype(np.float32),
        np.vstack(rows_alt).T.astype(np.float32),
    )


def _parse_vcf_to_dosage(
    vcf_path: Path,
) -> Tuple[List[str], List[str], np.ndarray]:
    """Stream VCF / `.vcf.gz` into `(sample_ids, locus_ids, dosage_matrix[n_samples, n_loci])`."""
    sids, lids, dmat, _, _ = _parse_vcf_with_depths(vcf_path)
    return sids, lids, dmat


def _parse_dosage_table(dosage_path: Path) -> Tuple[List[str], List[str], np.ndarray]:
    sep = "\t" if dosage_path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(dosage_path, sep=sep)
    if df.empty or df.shape[1] < 2:
        return [], [], np.empty((0, 0), dtype=np.float32)
    id_col = df.columns[0]
    sample_ids = [clean_sample_id(s) for s in df[id_col].tolist()]
    locus_cols = [str(c) for c in df.columns[1:]]
    mat = df[locus_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)
    # Normalize 0..2 diploid allele counts cleanly into 0..1 dosage when max > 1.05
    finite = mat[np.isfinite(mat)]
    if finite.size > 0 and float(np.max(finite)) > 1.05:
        max_p = 2.0 if float(np.max(finite)) <= 2.05 else float(np.max(finite))
        mat = mat / np.float32(max_p)
    return sample_ids, locus_cols, mat


def _load_db_vouchers(conn: Optional[sqlite3.Connection]) -> List[Dict[str, Any]]:
    if conn is None:
        return []
    try:
        cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(vouchers);").fetchall()}
    except sqlite3.Error:
        return []
    if not cols:
        return []
    select_cols = [
        c
        for c in (
            "catalogNumber", "recordNumber", "recordedBy", "scientificName",
            "species_standardized", "decimalLatitude", "decimalLongitude",
            "stateProvince", "county", "eventDate", "determiner_tier",
        )
        if c in cols
    ]
    rows = conn.execute(f"SELECT {', '.join(select_cols)} FROM vouchers;").fetchall()
    out: List[Dict[str, Any]] = []
    for r in rows:
        out.append({c: r[idx] for idx, c in enumerate(select_cols)})
    return out


def _match_field_record(
    sample_id: str,
    parsed: Dict[str, Optional[str]],
    field_rows: Sequence[Mapping[str, Any]],
) -> Optional[Mapping[str, Any]]:
    pop = (parsed.get("population_id") or "").upper()
    code = (parsed.get("taxon_code") or "").upper()
    sp_std = parsed.get("species_standardized")
    if not pop:
        return None

    exact_pid = [
        r
        for r in field_rows
        if r.get("project_id") and strip_smart_quotes(r["project_id"]).upper() == sample_id.upper()
    ]
    if exact_pid:
        return exact_pid[0]

    candidates = [r for r in field_rows if (r.get("population_id") or "").upper() == pop]
    if not candidates:
        for alt_pop in POPULATION_ALIAS_MAP.get(pop, (pop,))[1:]:
            alt_cands = [r for r in field_rows if (r.get("population_id") or "").upper() == alt_pop]
            if alt_cands:
                candidates = alt_cands
                break
    if not candidates:
        return None

    # Prefer novel-query records only for novel/hybrid pairs (avoid collapsing sympatric PACPAUVAL1 to sp. nov.)
    if _is_novel_or_hybrid_sample_pair(code, pop):
        novel_cands = [r for r in candidates if r.get("is_putative_novel_or_hybrid_query")]
        if novel_cands:
            return novel_cands[0]

    same_code = [r for r in candidates if (r.get("taxon_code") or "").upper() == code and code]
    if same_code:
        return same_code[0]

    same_sp = [r for r in candidates if sp_std and r.get("species_standardized") == sp_std]
    if same_sp:
        return same_sp[0]

    return candidates[0]


def _score_voucher_match(
    voucher: Mapping[str, Any],
    code: str,
    sp_std: Optional[str],
    pop_tok: str,
    is_alias_fallback: bool,
) -> Optional[Tuple[int, int, str]]:
    """Score compatible voucher match using anchored token regex + taxon compatibility."""
    rec = strip_smart_quotes(voucher.get("recordNumber") or "")
    cat = strip_smart_quotes(voucher.get("catalogNumber") or "")
    v_sp = str(voucher.get("species_standardized") or "").strip()

    pat = re.compile(rf"(?:^|[-_]){re.escape(pop_tok)}(?:[-_]|$)", re.IGNORECASE)
    rec_hit = bool(rec and pat.search(rec))
    cat_hit = bool(cat and pat.search(cat))
    if not (rec_hit or cat_hit):
        return None

    rec_up = rec.upper()
    code_up = code.upper()
    # Collector tag taxon compatibility check prevents e.g. PACSER matching PACAUR at NCSP
    prefix_compatible = False
    prefix_incompatible = False
    if rec_up.startswith("PAC") and len(rec_up) >= 6:
        rec_code = "NAUR" if rec_up.startswith("PACNAUR") else rec_up[3:6]
        if rec_code == code_up:
            prefix_compatible = True
        elif {rec_code, code_up}.issubset({"TOM", "DUB"}):
            prefix_compatible = True
        elif rec_code == "TX" and code_up in ("TOM", "DUB"):
            prefix_compatible = True
        else:
            prefix_incompatible = True

    if prefix_incompatible:
        return None

    sp_match = bool(sp_std and v_sp == sp_std)
    if not prefix_compatible and sp_std and v_sp and v_sp != sp_std:
        return None

    tier_rank = 0 if (voucher.get("determiner_tier") == "Tier_1_Gold") else (1 if voucher.get("determiner_tier") == "Tier_2_Silver" else 2)
    has_coords = 0 if (voucher.get("decimalLatitude") is not None and voucher.get("decimalLongitude") is not None) else 1
    alias_pen = 1 if is_alias_fallback else 0
    return (alias_pen, 0 if prefix_compatible else 1, 0 if sp_match else 1, tier_rank, has_coords, str(voucher.get("catalogNumber") or ""))


def _link_sample_to_voucher(
    sample_id: str,
    parsed: Dict[str, Optional[str]],
    user_meta: Optional[Mapping[str, Any]],
    field_row: Optional[Mapping[str, Any]],
    vouchers: Sequence[Mapping[str, Any]],
) -> Dict[str, Any]:
    """Resolve 4-tier voucher linkage (`G1_Direct_Voucher` .. `G4_Ungeoreferenced_Tip`) without fabricating IDs."""
    voucher_by_cat = {str(v.get("catalogNumber") or ""): v for v in vouchers if v.get("catalogNumber")}

    # 1. Explicit user_genetic_metadata catalogNumber -> G1_Direct_Voucher
    if user_meta and user_meta.get("catalogNumber"):
        u_cat = str(user_meta["catalogNumber"]).strip()
        matched = voucher_by_cat.get(u_cat)
        if matched:
            return {"catalogNumber": u_cat, "voucher": matched, "linkage_tier": "G1_Direct_Voucher"}

    # 2. Exact Project_ID / sample_id == catalogNumber or recordNumber -> G1_Direct_Voucher
    probe_ids = [sample_id]
    if user_meta and user_meta.get("project_id"):
        probe_ids.append(str(user_meta["project_id"]))
    for pid in probe_ids:
        pid_clean = strip_smart_quotes(pid).upper()
        if not pid_clean:
            continue
        for v in vouchers:
            cat_c = strip_smart_quotes(v.get("catalogNumber") or "").upper()
            rec_c = strip_smart_quotes(v.get("recordNumber") or "").upper()
            if pid_clean == cat_c or (rec_c and (pid_clean == rec_c or pid_clean.startswith(rec_c + "-"))):
                return {
                    "catalogNumber": str(v.get("catalogNumber")),
                    "voucher": v,
                    "linkage_tier": "G1_Direct_Voucher",
                }

    # 3. Population proxy match -> G2_Population_Voucher_Proxy
    pop = (parsed.get("population_id") or "").upper()
    code = (parsed.get("taxon_code") or "").upper()
    sp_std = parsed.get("species_standardized")
    if pop and vouchers:
        aliases = POPULATION_ALIAS_MAP.get(pop, (pop,))
        scored: List[Tuple[Tuple[Any, ...], Mapping[str, Any]]] = []
        for idx, tok in enumerate(aliases):
            is_fb = idx > 0
            tok_hits: List[Tuple[Tuple[Any, ...], Mapping[str, Any]]] = []
            for v in vouchers:
                sc = _score_voucher_match(v, code=code, sp_std=sp_std, pop_tok=tok, is_alias_fallback=is_fb)
                if sc is not None:
                    tok_hits.append((sc, v))
            if tok_hits:
                scored.extend(tok_hits)
                break
        if scored:
            scored.sort(key=lambda item: item[0])
            best_v = scored[0][1]
            return {
                "catalogNumber": str(best_v.get("catalogNumber")),
                "voucher": best_v,
                "linkage_tier": "G2_Population_Voucher_Proxy",
            }

    # 4. Georeferenced field sample vs ungeoreferenced tip -> G3 or G4
    has_user_coords = bool(
        user_meta
        and user_meta.get("decimalLatitude") is not None
        and user_meta.get("decimalLongitude") is not None
    )
    has_field_coords = bool(
        field_row
        and field_row.get("decimalLatitude") is not None
        and field_row.get("decimalLongitude") is not None
    )
    if has_user_coords or has_field_coords:
        return {"catalogNumber": None, "voucher": None, "linkage_tier": "G3_Georeferenced_Unvouchered"}
    return {"catalogNumber": None, "voucher": None, "linkage_tier": "G4_Ungeoreferenced_Tip"}


def _infer_sample_cytotype(
    ref_vec: np.ndarray,
    alt_vec: np.ndarray,
    dosage_vec: np.ndarray,
    user_ploidy: Optional[int],
) -> Tuple[str, Optional[float]]:
    """Infer cytotype (`2x`/`4x`/`6x`/`mixed/uncertain`) and delta-logLik from het allele balance."""
    valid = np.isfinite(ref_vec) & np.isfinite(alt_vec) & np.isfinite(dosage_vec)
    dp = ref_vec + alt_vec
    het_mask = valid & (dp >= 10.0) & (dosage_vec > 0.05) & (dosage_vec < 0.95)
    ratios = (alt_vec[het_mask] / np.maximum(1.0, dp[het_mask])).astype(np.float64)
    ratios = ratios[(ratios >= 0.10) & (ratios <= 0.90)]
    if ratios.size < 15:
        if user_ploidy in (2, 4, 6):
            return f"{user_ploidy}x", 0.0
        return "mixed/uncertain", None

    sigma = 0.08
    inv_var = 1.0 / (2.0 * sigma * sigma)

    def _mix_ll(centers: Sequence[float]) -> float:
        c = np.asarray(centers, dtype=np.float64)[None, :]
        diffs = ratios[:, None] - c
        comp = -0.5 * np.log(2.0 * np.pi * sigma * sigma) - (diffs * diffs) * inv_var - np.log(len(centers))
        m = np.max(comp, axis=1)
        return float(np.sum(m + np.log(np.sum(np.exp(comp - m[:, None]), axis=1))))

    ll_2x = _mix_ll([0.50])
    ll_4x = _mix_ll([0.25, 0.50, 0.75]) - 1.5 * float(np.log(max(2, ratios.size)))
    ll_6x = _mix_ll([1.0 / 6.0, 2.0 / 6.0, 0.50, 4.0 / 6.0, 5.0 / 6.0]) - 3.0 * float(np.log(max(2, ratios.size)))
    delta = float(max(ll_4x, ll_6x) - ll_2x)
    if delta <= 2.0:
        return "2x", round(delta, 4)
    if (ll_6x - ll_4x) > 3.0:
        return "6x", round(delta, 4)
    return "4x", round(delta, 4)


def _evaluate_hdplot_loci(
    sub_pass: np.ndarray,
    sub_ref: np.ndarray,
    ref_d_ref: np.ndarray,
    alt_d_ref: np.ndarray,
    locus_ids: Sequence[str],
    h_max: float = 0.60,
    d_max: float = 4.0,
) -> Tuple[np.ndarray, pd.DataFrame]:
    """Apply McKinney et al. (2017) HDplot paralog filter (`H <= 0.60`, `|D| <= 4.0`) without MAF pruning."""
    n_loci = len(locus_ids)
    keep_mask = np.ones(n_loci, dtype=bool)
    diag_rows: List[Dict[str, Any]] = []

    has_ref_depths = ref_d_ref.size > 0 and alt_d_ref.size > 0
    for j in range(n_loci):
        col_p = sub_pass[:, j]
        fin_p = col_p[np.isfinite(col_p)]
        col_r = sub_ref[:, j]
        fin_mask_r = np.isfinite(col_r)
        fin_r = col_r[fin_mask_r]
        n_called_r = int(fin_r.size)
        het_mask_r = fin_mask_r & (col_r > 0.05) & (col_r < 0.95)
        n_hets_r = int(np.sum(het_mask_r))
        h_obs = float(n_hets_r / float(n_called_r)) if n_called_r > 0 else 0.0

        na_sum = 0.0
        nb_sum = 0.0
        d_stat = 0.0
        if has_ref_depths and n_hets_r > 0:
            rr = ref_d_ref[het_mask_r, j]
            aa = alt_d_ref[het_mask_r, j]
            ok_d = np.isfinite(rr) & np.isfinite(aa)
            if np.any(ok_d):
                na_sum = float(np.sum(rr[ok_d]))
                nb_sum = float(np.sum(aa[ok_d]))
                n_tot = na_sum + nb_sum
                if n_tot >= 10.0:
                    d_stat = float((na_sum - 0.5 * n_tot) / np.sqrt(0.25 * n_tot))

        reason = "PASS"
        if fin_p.size < 2:
            keep_mask[j] = False
            reason = "monomorphic_or_sparse"
        elif n_called_r >= 4 and h_obs > float(h_max):
            keep_mask[j] = False
            reason = f"HDplot_H>{h_max:.2f}"
        elif (na_sum + nb_sum) >= 10.0 and abs(d_stat) > float(d_max):
            keep_mask[j] = False
            reason = f"HDplot_|D|>{d_max:.1f}"
        elif float(np.std(fin_p)) <= 1e-8:
            keep_mask[j] = False
            reason = "monomorphic_or_sparse"

        diag_rows.append({
            "locus_id": locus_ids[j],
            "n_called_ref": n_called_r,
            "n_hets_ref": n_hets_r,
            "H_obs": round(h_obs, 4) if n_called_r > 0 else np.nan,
            "het_ref_reads": round(na_sum, 1),
            "het_alt_reads": round(nb_sum, 1),
            "het_total_reads": round(na_sum + nb_sum, 1),
            "D_read_ratio": round(d_stat, 4),
            "hdplot_keep": 1 if keep_mask[j] else 0,
            "filter_reason": reason,
        })

    return keep_mask, pd.DataFrame(diag_rows)


def _compute_eb_polyrad_dosage(
    dosage_sub: np.ndarray,
    ref_sub: np.ndarray,
    alt_sub: np.ndarray,
    cytotypes: Sequence[str],
    eps: float = 0.005,
) -> np.ndarray:
    """Estimate Empirical Bayes discrete polyploid posterior dosage in `[0, 1]` from allele read counts."""
    out = dosage_sub.copy().astype(np.float32)
    if ref_sub.size == 0 or alt_sub.size == 0 or out.size == 0:
        return out

    from math import comb

    p_hat = np.nanmean(dosage_sub, axis=0)
    p_hat = np.where(np.isfinite(p_hat), np.clip(p_hat, 0.01, 0.99), 0.5).astype(np.float64)

    for i, cy in enumerate(cytotypes):
        ploidy_p = 2 if cy == "2x" else (6 if cy == "6x" else 4)
        ks = np.arange(ploidy_p + 1, dtype=np.float64)
        ratios_k = ks / float(ploidy_p)
        ek = (1.0 - eps) * ratios_k + eps * (1.0 - ratios_k)
        ek = np.clip(ek, 1e-6, 1.0 - 1e-6)
        log_comb = np.log(np.asarray([comb(ploidy_p, int(k)) for k in ks], dtype=np.float64))

        nr = ref_sub[i, :].astype(np.float64)
        na = alt_sub[i, :].astype(np.float64)
        called = np.isfinite(dosage_sub[i, :]) & np.isfinite(nr) & np.isfinite(na) & ((nr + na) > 0.0)
        if not np.any(called):
            continue

        p_c = p_hat[called][:, None]
        prior = log_comb[None, :] + ks[None, :] * np.log(p_c) + (ploidy_p - ks[None, :]) * np.log(1.0 - p_c)
        ll = na[called][:, None] * np.log(ek[None, :]) + nr[called][:, None] * np.log(1.0 - ek[None, :])
        log_post = prior + ll
        log_post -= np.max(log_post, axis=1, keepdims=True)
        w = np.exp(log_post)
        w /= np.sum(w, axis=1, keepdims=True)
        out[i, called] = np.sum(w * ratios_k[None, :], axis=1).astype(np.float32)

    return np.clip(out, 0.0, 1.0)


def _read_fasta_headers_seqs(path: Path) -> List[Tuple[str, str]]:
    entries: List[Tuple[str, str]] = []
    cur_h: Optional[str] = None
    buf: List[str] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if cur_h is not None:
                    entries.append((cur_h, "".join(buf).upper()))
                cur_h = s[1:].split()[0]
                buf = []
            else:
                buf.append(s)
    if cur_h is not None:
        entries.append((cur_h, "".join(buf).upper()))
    return entries


def _load_organellar_metadata(organellar_files: Sequence[Path]) -> Dict[str, Dict[str, Any]]:
    """Read staged cytonuclear table and skimmed plastid/ITS consensus FASTAs."""
    meta: Dict[str, Dict[str, Any]] = {}
    for p in organellar_files:
        low = p.name.lower()
        if low.startswith("cytonuclear_discordance") and low.endswith((".csv", ".tsv")):
            sep = "\t" if low.endswith(".tsv") else ","
            try:
                df = pd.read_csv(p, sep=sep, dtype=str, keep_default_na=False)
            except Exception:
                continue
            for _, row in df.iterrows():
                sid = clean_sample_id(row.get("genetic_sample_id", ""))
                if not sid:
                    continue
                rec = meta.setdefault(sid, {})
                if row.get("plastid_clade"):
                    rec["plastid_clade"] = str(row["plastid_clade"]).strip()
                if row.get("its_ribotype"):
                    rec["its_ribotype"] = str(row["its_ribotype"]).strip()
                disc = str(row.get("cytonuclear_discordant", "")).strip().lower()
                if disc in ("1", "true", "yes"):
                    rec["cytonuclear_discordant"] = 1
                elif disc in ("0", "false", "no"):
                    rec["cytonuclear_discordant"] = 0
        elif low.endswith((".fasta", ".fa", ".fas")):
            target_col = "plastid_clade" if "plastid" in low else ("its_ribotype" if "its" in low else None)
            if not target_col:
                continue
            seq_map: Dict[str, str] = {}
            for hdr, seq in _read_fasta_headers_seqs(p):
                sid = clean_sample_id(hdr.split("|")[0])
                clade_tag = None
                for tok in hdr.split("|")[1:]:
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        if k.lower() in ("clade", "ribotype", target_col):
                            clade_tag = v.strip()
                if sid:
                    rec = meta.setdefault(sid, {})
                    if not rec.get(target_col):
                        if clade_tag:
                            rec[target_col] = clade_tag
                        else:
                            # Deterministic 4-mer signature hash if clade header metadata absent
                            if seq not in seq_map:
                                seq_map[seq] = f"{'cp' if target_col == 'plastid_clade' else 'ITS'}_H{len(seq_map) + 1}"
                            rec[target_col] = seq_map[seq]
    return meta


def run_genomic_ingestion(
    genomics_dir: Path,
    out_dir: Path,
    db_path: Optional[Path] = None,
    min_hq_reads: int = 50000,
    max_missing_rate: float = 0.80,
    paralog_ho_threshold: float = 0.60,
    hdplot_d_threshold: float = 4.0,
    demo: bool = False,
) -> Dict[str, Any]:
    """Ingest staged genomic assets, apply HDplot & EB dosage estimation, link vouchers, and upsert Tier A SQLite."""
    gen_res = genomics_dir.resolve()
    out_res = out_dir.resolve()
    if out_res != gen_res and gen_res in out_res.parents:
        raise ValueError(f"Refusing to create subdirectory {out_dir} inside flat genomics directory {genomics_dir}")

    assets = discover_genomic_inputs(genomics_dir)
    if not assets.has_genetic_payload:
        logger.info("No active genomic payloads found in %s; skipping genomic ingestion (no-op).", genomics_dir)
        return {
            "status": "skipped_empty",
            "n_samples_total": 0,
            "n_samples_pass": 0,
            "n_loci_raw": 0,
            "n_loci_retained": 0,
        }

    field_rows, direct_user_meta = _parse_field_metadata(assets.metadata_csvs)
    organellar_map = _load_organellar_metadata(assets.organellar_files)
    stats_map: Dict[str, Dict[str, Any]] = {}
    for sp in assets.stats_files:
        stats_map.update(_parse_stats_file(sp))

    dosage_samples: List[str] = []
    locus_ids: List[str] = []
    dosage_mat = np.empty((0, 0), dtype=np.float32)
    ref_depth_mat = np.empty((0, 0), dtype=np.float32)
    alt_depth_mat = np.empty((0, 0), dtype=np.float32)
    from_precomputed_dosage = False
    marker_platform = "nextRAD"
    if assets.variant_files:
        primary_vcf = assets.variant_files[0]
        if "hifi" in primary_vcf.name.lower() or "pacbio" in primary_vcf.name.lower():
            marker_platform = "PacBio_HiFi"
        elif "ddrad" in primary_vcf.name.lower():
            marker_platform = "ddRAD"
        dosage_samples, locus_ids, dosage_mat, ref_depth_mat, alt_depth_mat = _parse_vcf_with_depths(primary_vcf)
    if assets.dosage_files:
        d_sids, d_loci, d_mat = _parse_dosage_table(assets.dosage_files[0])
        if d_sids and d_loci:
            marker_platform = "PacBio_HiFi" if not assets.variant_files else marker_platform
            dosage_samples, locus_ids, dosage_mat = d_sids, d_loci, d_mat
            ref_depth_mat = np.empty((0, 0), dtype=np.float32)
            alt_depth_mat = np.empty((0, 0), dtype=np.float32)
            from_precomputed_dosage = True

    all_sample_ids: List[str] = []
    for src_ids in (dosage_samples, list(stats_map.keys()), _parse_key_samples(assets.key_files), list(direct_user_meta.keys())):
        for sid in src_ids:
            if sid and sid not in all_sample_ids:
                all_sample_ids.append(sid)

    conn: Optional[sqlite3.Connection] = None
    vouchers: List[Dict[str, Any]] = []
    try:
        conn = get_db_connection(
            db_path=str(db_path) if db_path else None,
            read_only=demo,
            demo=False,
            init_schema=(not demo),
            workspace_root=PROJECT_ROOT,
        )
        vouchers = _load_db_vouchers(conn)
    except Exception as exc:
        logger.warning("SQLite voucher lookup unavailable (%s); continuing with standalone linking.", exc)
        vouchers = []

    sample_idx_map = {sid: i for i, sid in enumerate(dosage_samples)}
    n_loci_raw = int(dosage_mat.shape[1]) if dosage_mat.ndim == 2 else 0

    records: List[Dict[str, Any]] = []
    for sid in all_sample_ids:
        parsed = parse_sample_code(sid)
        u_meta = direct_user_meta.get(sid)
        f_row = _match_field_record(sid, parsed, field_rows)
        st = stats_map.get(sid, {})
        org = organellar_map.get(sid, {})

        raw_t = parsed.get("taxon_raw")
        sp_std = parsed.get("species_standardized")
        pop_up = (parsed.get("population_id") or "").upper()
        code_up = (parsed.get("taxon_code") or "").upper()
        force_novel_pair = _is_novel_or_hybrid_sample_pair(code_up, pop_up)
        f_same_taxon = bool(f_row is not None and (not code_up or (f_row.get("taxon_code") or "").upper() == code_up))
        is_query = False
        if f_row is not None:
            if force_novel_pair or (f_same_taxon and f_row.get("is_putative_novel_or_hybrid_query")):
                raw_t = f_row.get("taxon_raw") or raw_t
                sp_std = f_row.get("species_standardized") or sp_std
                is_query = True
            elif not sp_std and f_same_taxon:
                raw_t = f_row.get("taxon_raw")
                sp_std = f_row.get("species_standardized")
                is_query = bool(f_row.get("is_putative_novel_or_hybrid_query"))
        if force_novel_pair:
            is_query = True
            if not raw_t:
                raw_t = (f_row or {}).get("taxon_raw") or "Packera sp. nov."
            if not sp_std or (f_row is not None and f_row.get("is_putative_novel_or_hybrid_query")):
                sp_std = (f_row or {}).get("species_standardized") or raw_t
        if u_meta is not None:
            if u_meta.get("taxon_raw") or u_meta.get("species_standardized"):
                ut_raw, ut_std, ut_q = normalize_taxon_name(u_meta.get("taxon_raw") or u_meta.get("species_standardized"))
                raw_t = ut_raw or raw_t
                sp_std = ut_std or sp_std
                is_query = is_query or ut_q
            if str(u_meta.get("is_putative_novel_or_hybrid_query", "")).strip().lower() in ("1", "true", "yes"):
                is_query = True

        hq_reads = st.get("hq_reads")
        mapped_reads = st.get("mapped_reads")
        pct_rep = st.get("pct_repetitive")
        n_called = st.get("n_called_loci")
        miss_rate = st.get("missing_genotype_rate")
        obs_het = st.get("observed_het")

        ploidy = None
        if u_meta and u_meta.get("ploidy_level"):
            try:
                ploidy = int(float(u_meta["ploidy_level"]))
            except ValueError:
                ploidy = None

        inf_cyto = "mixed/uncertain" if ploidy not in (2, 4, 6) else f"{ploidy}x"
        nq_delta: Optional[float] = None
        if sid in sample_idx_map and n_loci_raw > 0:
            sidx = sample_idx_map[sid]
            vec = dosage_mat[sidx, :]
            called_mask = np.isfinite(vec)
            n_called = int(np.sum(called_mask))
            miss_rate = float(1.0 - (n_called / float(n_loci_raw)))
            obs_het = float(np.mean((vec[called_mask] > 0.05) & (vec[called_mask] < 0.95))) if n_called > 0 else None
            if ref_depth_mat.size > 0 and alt_depth_mat.size > 0:
                inf_cyto, nq_delta = _infer_sample_cytotype(ref_depth_mat[sidx, :], alt_depth_mat[sidx, :], vec, ploidy)

        if hq_reads is not None:
            qc_pass = (hq_reads >= int(min_hq_reads)) and (miss_rate is None or miss_rate <= float(max_missing_rate))
        elif n_loci_raw > 0:
            qc_pass = (miss_rate is not None and miss_rate <= float(max_missing_rate)) and ((n_called or 0) > 0)
        else:
            qc_pass = True

        link = _link_sample_to_voucher(sid, parsed, u_meta, f_row, vouchers)
        v_rec = link.get("voucher") or {}

        lat = lon = state = county = ev_date = None
        if u_meta and u_meta.get("decimalLatitude") and u_meta.get("decimalLongitude"):
            try:
                lat, lon = float(u_meta["decimalLatitude"]), float(u_meta["decimalLongitude"])
            except ValueError:
                lat = lon = None
        if lat is None and f_row is not None:
            lat, lon = f_row.get("decimalLatitude"), f_row.get("decimalLongitude")
            state, county, ev_date = f_row.get("stateProvince"), f_row.get("county"), f_row.get("eventDate")
        if lat is None and v_rec:
            lat, lon = v_rec.get("decimalLatitude"), v_rec.get("decimalLongitude")
        if state is None and v_rec:
            state = v_rec.get("stateProvince")
        if county is None and v_rec:
            county = v_rec.get("county")
        if ev_date is None and v_rec:
            ev_date = v_rec.get("eventDate")
        if u_meta:
            state = u_meta.get("stateProvince") or state
            county = u_meta.get("county") or county
            ev_date = u_meta.get("eventDate") or ev_date

        records.append({
            "genetic_sample_id": sid,
            "catalogNumber": link.get("catalogNumber"),
            "project_id": (
                (u_meta or {}).get("project_id")
                or ((f_row or {}).get("project_id") if (f_same_taxon or force_novel_pair) else None)
                or parsed.get("replicate_pair_id")
            ),
            "population_id": (u_meta or {}).get("population_id") or parsed.get("population_id"),
            "replicate_pair_id": (u_meta or {}).get("replicate_pair_id") or parsed.get("replicate_pair_id"),
            "taxon_raw": raw_t,
            "species_standardized": sp_std,
            "is_putative_novel_or_hybrid_query": 1 if is_query else 0,
            "ploidy_level": ploidy,
            "marker_platform": (u_meta or {}).get("marker_platform") or marker_platform,
            "decimalLatitude": lat,
            "decimalLongitude": lon,
            "stateProvince": state,
            "county": county,
            "eventDate": ev_date,
            "hq_reads": hq_reads,
            "mapped_reads": mapped_reads,
            "pct_repetitive": pct_rep,
            "n_called_loci": n_called,
            "missing_genotype_rate": miss_rate,
            "observed_het": obs_het,
            "geno_quality_pass": 1 if qc_pass else 0,
            "linkage_tier": link.get("linkage_tier"),
            "hdplot_retained_loci": 0,
            "inferred_cytotype": (u_meta or {}).get("inferred_cytotype") or inf_cyto,
            "nquire_delta_logl": nq_delta,
            "plastid_clade": (u_meta or {}).get("plastid_clade") or org.get("plastid_clade"),
            "its_ribotype": (u_meta or {}).get("its_ribotype") or org.get("its_ribotype"),
            "cytonuclear_discordant": org.get("cytonuclear_discordant"),
        })

    out_dir.mkdir(parents=True, exist_ok=True)
    n_loci_retained = 0
    dosage_out_path = out_dir / "allele_dosage_matrix.tsv"
    hdplot_out_path = out_dir / "hdplot_locus_diagnostics.tsv"

    if dosage_samples and n_loci_raw > 0:
        pass_recs = [r for r in records if r["geno_quality_pass"] == 1 and r["genetic_sample_id"] in sample_idx_map]
        pass_ids = [r["genetic_sample_id"] for r in pass_recs]
        ref_pass_ids = [r["genetic_sample_id"] for r in pass_recs if r["is_putative_novel_or_hybrid_query"] == 0]
        pass_rows = [sample_idx_map[s] for s in pass_ids]
        ref_rows = [sample_idx_map[s] for s in ref_pass_ids] if len(ref_pass_ids) >= 4 else pass_rows

        if pass_rows:
            sub_pass = dosage_mat[pass_rows, :]
            sub_ref = dosage_mat[ref_rows, :]
            rd_ref = ref_depth_mat[ref_rows, :] if ref_depth_mat.size > 0 else np.empty((0, 0), dtype=np.float32)
            ad_ref = alt_depth_mat[ref_rows, :] if alt_depth_mat.size > 0 else np.empty((0, 0), dtype=np.float32)

            keep_mask, hd_df = _evaluate_hdplot_loci(
                sub_pass, sub_ref, rd_ref, ad_ref, locus_ids,
                h_max=paralog_ho_threshold, d_max=hdplot_d_threshold,
            )
            hd_df.to_csv(hdplot_out_path, sep="\t", index=False, na_rep="NA")

            kept_indices = np.where(keep_mask)[0]
            n_loci_retained = int(kept_indices.size)
            kept_loci = [locus_ids[int(j)] for j in kept_indices]
            kept_dosage = sub_pass[:, kept_indices]
            if not from_precomputed_dosage and ref_depth_mat.size > 0 and n_loci_retained > 0:
                rd_pass = ref_depth_mat[pass_rows, :][:, kept_indices]
                ad_pass = alt_depth_mat[pass_rows, :][:, kept_indices]
                cytos = [str(r["inferred_cytotype"]) for r in pass_recs]
                kept_dosage = _compute_eb_polyrad_dosage(kept_dosage, rd_pass, ad_pass, cytos)

            for pos_i, r in enumerate(pass_recs):
                r["hdplot_retained_loci"] = int(np.sum(np.isfinite(kept_dosage[pos_i, :]))) if n_loci_retained > 0 else 0

            dos_df = pd.DataFrame(np.round(kept_dosage, 4), columns=kept_loci)
            dos_df.insert(0, "genetic_sample_id", pass_ids)
            dos_df.to_csv(dosage_out_path, sep="\t", index=False, na_rep="NA")

    samples_df = pd.DataFrame(records)
    samples_csv = out_dir / "parsed_genomic_samples.csv"
    samples_df.to_csv(samples_csv, index=False)

    manifest = {
        "genomics_dir": str(genomics_dir),
        "out_dir": str(out_dir),
        "parsed_samples_csv": str(samples_csv),
        "allele_dosage_tsv": str(dosage_out_path) if dosage_out_path.exists() else None,
        "hdplot_diagnostics_tsv": str(hdplot_out_path) if hdplot_out_path.exists() else None,
        "tree_files": [str(p) for p in assets.tree_files],
        "alignment_files": [str(p) for p in assets.alignment_files],
        "organellar_files": [str(p) for p in assets.organellar_files],
        "admixture_q_files": [str(p) for p in assets.admixture_q_files],
        "n_samples_total": len(records),
        "n_samples_pass": int(sum(1 for r in records if r["geno_quality_pass"] == 1)),
        "n_loci_raw": n_loci_raw,
        "n_loci_retained": n_loci_retained,
    }
    manifest_path = out_dir / "genomic_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    if conn is not None:
        try:
            if not demo:
                upsert_genomic_samples(conn, records, preserve_existing_on_null=False)
        finally:
            conn.close()

    logger.info(
        "Ingested %d genomic samples (%d QC-passing, %d/%d singleton loci retained after HDplot H<=%.2f, |D|<=%.1f).",
        manifest["n_samples_total"],
        manifest["n_samples_pass"],
        n_loci_retained,
        n_loci_raw,
        paralog_ho_threshold,
        hdplot_d_threshold,
    )
    return {"status": "ingested", **manifest}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Ingest staged genomic data into Packera Tier A SQLite.")
    parser.add_argument("--genomics-dir", type=Path, default=PROJECT_ROOT / "data" / "genomics")
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "outputs" / "reports")
    parser.add_argument("--db-path", "--db", dest="db_path", type=Path, default=None)
    parser.add_argument("--min-hq-reads", type=int, default=50000)
    parser.add_argument("--max-missing", type=float, default=0.80)
    parser.add_argument("--paralog-ho", type=float, default=0.60)
    parser.add_argument("--hdplot-d", type=float, default=4.0)
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--check-only", action="store_true", help="Print HAS_GENOMIC_INPUTS=0/1 and exit 0 cleanly.")
    args = parser.parse_args(list(argv) if argv is not None else None)

    assets = discover_genomic_inputs(args.genomics_dir)
    if args.check_only:
        print(f"HAS_GENOMIC_INPUTS={1 if assets.has_genetic_payload else 0}")
        return 0

    run_genomic_ingestion(
        genomics_dir=args.genomics_dir,
        out_dir=args.out_dir,
        db_path=args.db_path,
        min_hq_reads=args.min_hq_reads,
        max_missing_rate=args.max_missing,
        paralog_ho_threshold=args.paralog_ho,
        hdplot_d_threshold=args.hdplot_d,
        demo=bool(args.demo),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
