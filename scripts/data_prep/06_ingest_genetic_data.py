#!/usr/bin/env python3
"""
Module: scripts/data_prep/06_ingest_genetic_data.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  Independent Step 06 Tier A Genomic Evidence Stream ingestion, HDplot paralog
  filter, taxon-stratified Empirical Bayes polyploid dosage estimation, nQuire
  cytotype mixture scoring, organellar/cytonuclear ingestion, and 4-tier voucher linker.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import sqlite3
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import resolve_genomics_config, resolve_max_workers
from scripts.core.db import get_db_connection, upsert_genomic_samples
from scripts.core.logger import setup_logging
from scripts.data_prep.step06_public_cohort_ingestor import has_valid_public_genomic_records, resolve_and_stage_dynamic_genomic_cohort
from scripts.data_prep.step06_ploidy_hdplot_utils import (
    _CYTO_TO_PLOIDY, _compute_eb_polyrad_dosage, _em_nquire_component, _evaluate_hdplot_loci,
    _fit_nquire_em_models, _infer_sample_cytotype, _load_paraloss1272_catalog, _lookup_paraloss1272_locus,
)
from scripts.data_prep.step06_sample_metadata_utils import (
    ALIGNMENT_SUFFIXES, DOSAGE_HINTS, OUTGROUP_CODE_TO_TAXON, PAC_CODE_TO_TAXON, TREE_PREFIXES, TREE_SUFFIXES, VARIANT_SUFFIXES,
    StagedGenomicAssets, _canonical_taxon_equiv_key, _normalize_event_date_part, _normalize_exsiccatae_key, _parse_field_metadata,
    _parse_key_samples, _parse_stats_file, clean_sample_id, discover_genomic_inputs, normalize_taxon_name, parse_sample_code,
    sanitize_step06_sample_coordinates, strip_smart_quotes,
)
from scripts.data_prep.step06_vcf_dosage_utils import (
    _extract_cell_depths, _load_organellar_metadata, _parse_dosage_table, _parse_vcf_with_depths, align_vcf_depths_to_dosage_grid,
)
from scripts.data_prep.step06_voucher_linker_utils import (
    _build_exsiccatae_index, _extract_meta_exsiccatae_keys, _link_sample_to_voucher, _load_db_cytology_priors,
    _load_db_vouchers, _match_field_record, _normalize_exsiccatae_lookup_key, _parse_int_ploidy,
)

logger = setup_logging(name="plantmatics.genomics.ingest")


def _resolve_sample_record(
    sid: str,
    vouchers: Sequence[Dict[str, Any]],
    exsiccatae_index: Dict[str, List[Dict[str, Any]]],
    field_rows: Sequence[Dict[str, Any]],
    direct_user_meta: Dict[str, Dict[str, Any]],
    stats_map: Dict[str, Dict[str, Any]],
    organellar_map: Dict[str, Dict[str, Any]],
    db_priors: Dict[str, Dict[str, Any]],
    ex_prior_lc: Dict[str, int],
    dosage_mat: np.ndarray,
    ref_depth_mat: np.ndarray,
    alt_depth_mat: np.ndarray,
    sample_idx_map: Dict[str, int],
    n_loci_raw: int,
    marker_platform: str,
    min_hq_reads: int,
    max_missing_rate: float,
    cyto_ctx_out: Optional[Dict[str, Tuple[Optional[int], Optional[int], bool, bool]]] = None,
) -> Tuple[Dict[str, Any], Optional[int]]:
    parsed = parse_sample_code(sid)
    u_meta = direct_user_meta.get(sid)
    f_row = _match_field_record(sid, parsed, field_rows)
    st = stats_map.get(sid, {})
    org = organellar_map.get(sid, {})

    raw_t = parsed.get("taxon_raw")
    sp_std = parsed.get("species_standardized")
    code_up = (parsed.get("taxon_code") or "").upper()
    is_nov_code = code_up == "NOV"
    f_same_taxon = bool(f_row is not None and (not code_up or (f_row.get("taxon_code") or "").upper() == code_up))
    is_query = is_nov_code
    if f_row is not None:
        if is_nov_code or (f_same_taxon and f_row.get("is_putative_novel_or_hybrid_query")):
            raw_t = f_row.get("taxon_raw") or raw_t
            sp_std = f_row.get("species_standardized") or sp_std
            is_query = is_query or bool(f_row.get("is_putative_novel_or_hybrid_query"))
        elif not sp_std and f_same_taxon:
            raw_t = f_row.get("taxon_raw")
            sp_std = f_row.get("species_standardized")
            is_query = bool(f_row.get("is_putative_novel_or_hybrid_query"))
    if is_nov_code:
        raw_t = raw_t or (f_row or {}).get("taxon_raw") or "Packera sp. nov."
        sp_std = sp_std or (f_row or {}).get("species_standardized") or raw_t
    if u_meta is not None:
        if u_meta.get("taxon_raw") or u_meta.get("species_standardized"):
            ut_raw, ut_std, ut_q = normalize_taxon_name(u_meta.get("taxon_raw") or u_meta.get("species_standardized"))
            raw_t = ut_raw or raw_t
            sp_std = ut_std or sp_std
            is_query = is_query or ut_q
        if str(u_meta.get("is_putative_novel_or_hybrid_query", "")).strip().lower() in ("1", "true", "yes"):
            is_query = True

    hq_reads, mapped_reads, pct_rep = st.get("hq_reads"), st.get("mapped_reads"), st.get("pct_repetitive")
    n_called, miss_rate, obs_het = st.get("n_called_loci"), st.get("missing_genotype_rate"), st.get("observed_het")
    link = _link_sample_to_voucher(sid, parsed, u_meta, f_row, vouchers, exsiccatae_index=exsiccatae_index)
    v_rec = link.get("voucher") or {}

    ploidy = _parse_int_ploidy(u_meta.get("ploidy_level")) if u_meta else None
    if ploidy is None and (cat_link := str(link.get("catalogNumber") or "").strip()):
        ploidy = db_priors["by_catalog"].get(cat_link)
    if ploidy is None and (ex_link := _normalize_exsiccatae_lookup_key(link.get("exsiccatae_key")).lower()):
        ploidy = ex_prior_lc.get(ex_link)
    if ploidy is None and (gb_acc := str((u_meta or {}).get("genbank_accession") or parsed.get("genbank_accession") or (f_row or {}).get("genbank_accession") or "").strip()):
        ploidy = db_priors["by_genbank"].get(gb_acc)

    tax_sp = str(sp_std or v_rec.get("species_standardized") or "").strip()
    tax_info = None
    if not is_query and tax_sp:
        tax_info = db_priors["by_taxon"].get(tax_sp)
        if tax_info is None:
            eq = _canonical_taxon_equiv_key(tax_sp)
            tax_info = next((v for k, v in db_priors["by_taxon"].items() if _canonical_taxon_equiv_key(k) == eq), None)
    sp_modal = tax_info.get("modal_ploidy_level") if tax_info else None
    sp_multi = bool(tax_info.get("is_multi_cytotype_taxon", 0)) if tax_info else False

    inf_cyto = f"{ploidy}x" if ploidy in (2, 3, 4, 6, 8) else (f"{sp_modal}x" if (sp_modal in (2, 3, 4, 6, 8) and not sp_multi) else "mixed/uncertain")
    nq_delta: Optional[float] = None
    has_user_cyto = bool(str((u_meta or {}).get("inferred_cytotype") or "").strip())
    if cyto_ctx_out is not None:
        cyto_ctx_out[sid] = (ploidy, sp_modal, sp_multi, has_user_cyto)
    if sid in sample_idx_map and n_loci_raw > 0:
        sidx = sample_idx_map[sid]
        vec = dosage_mat[sidx, :]
        called_mask = np.isfinite(vec)
        n_called = int(np.sum(called_mask))
        miss_rate = float(1.0 - (n_called / float(n_loci_raw)))
        obs_het = float(np.mean((vec[called_mask] > 0.05) & (vec[called_mask] < 0.95))) if n_called > 0 else None

    fb_ploidy = ploidy if ploidy in (2, 3, 4, 6, 8) else (sp_modal if sp_modal in (2, 3, 4, 6, 8) else None)
    if hq_reads is not None:
        qc_pass = (hq_reads >= int(min_hq_reads)) and (miss_rate is None or miss_rate <= float(max_missing_rate))
    elif n_loci_raw > 0:
        qc_pass = (miss_rate is not None and miss_rate <= float(max_missing_rate)) and ((n_called or 0) > 0)
    else:
        qc_pass = True

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
    if state is None and v_rec: state = v_rec.get("stateProvince")
    if county is None and v_rec: county = v_rec.get("county")
    if ev_date is None and v_rec: ev_date = v_rec.get("eventDate")
    if u_meta:
        state = u_meta.get("stateProvince") or state
        county = u_meta.get("county") or county
        ev_date = u_meta.get("eventDate") or ev_date

    is_vouchered = bool(link.get("catalogNumber"))
    loc_txt = (u_meta or {}).get("locality") or (f_row or {}).get("notes") or ""
    lat, lon, was_nullified = sanitize_step06_sample_coordinates(
        lat, lon, state_prov=state or "", county=county or "", locality=loc_txt, is_vouchered=is_vouchered,
    )
    if not is_vouchered and was_nullified and link.get("linkage_tier") == "G3_Georeferenced_Unvouchered":
        link["linkage_tier"] = "G4_Ungeoreferenced_Tip"

    rec = {
        "genetic_sample_id": sid,
        "catalogNumber": link.get("catalogNumber"),
        "exsiccatae_key": link.get("exsiccatae_key"),
        "project_id": (u_meta or {}).get("project_id") or ((f_row or {}).get("project_id") if (f_same_taxon or is_nov_code) else None) or parsed.get("replicate_pair_id"),
        "population_id": (u_meta or {}).get("population_id") or parsed.get("population_id"),
        "replicate_pair_id": (u_meta or {}).get("replicate_pair_id") or parsed.get("replicate_pair_id"),
        "taxon_raw": raw_t, "species_standardized": sp_std,
        "is_putative_novel_or_hybrid_query": 1 if is_query else 0,
        "ploidy_level": ploidy, "marker_platform": (u_meta or {}).get("marker_platform") or marker_platform,
        "decimalLatitude": lat, "decimalLongitude": lon,
        "stateProvince": state, "county": county, "eventDate": ev_date,
        "hq_reads": hq_reads, "mapped_reads": mapped_reads, "pct_repetitive": pct_rep,
        "n_called_loci": n_called, "missing_genotype_rate": miss_rate, "observed_het": obs_het,
        "geno_quality_pass": 1 if qc_pass else 0, "linkage_tier": link.get("linkage_tier"),
        "hdplot_retained_loci": 0, "inferred_cytotype": (u_meta or {}).get("inferred_cytotype") or inf_cyto,
        "nquire_delta_logl": nq_delta,
        "plastid_clade": (u_meta or {}).get("plastid_clade") or org.get("plastid_clade"),
        "its_ribotype": (u_meta or {}).get("its_ribotype") or org.get("its_ribotype"),
        "cytonuclear_discordant": org.get("cytonuclear_discordant"),
        "public_plastome_match": org.get("public_plastome_match"),
        "public_its_match": org.get("public_its_match"),
        "waxy_lineage_match": org.get("waxy_lineage_match"),
    }
    return rec, fb_ploidy


_EMPTY_SAMPLE_COLS = (
    "genetic_sample_id", "catalogNumber", "exsiccatae_key", "project_id", "population_id", "replicate_pair_id",
    "taxon_raw", "species_standardized", "is_putative_novel_or_hybrid_query", "ploidy_level", "marker_platform",
    "decimalLatitude", "decimalLongitude", "stateProvince", "county", "eventDate", "hq_reads", "mapped_reads",
    "pct_repetitive", "n_called_loci", "missing_genotype_rate", "observed_het", "geno_quality_pass", "linkage_tier",
    "hdplot_retained_loci", "inferred_cytotype", "nquire_delta_logl", "plastid_clade", "its_ribotype",
    "cytonuclear_discordant", "public_plastome_match", "public_its_match", "waxy_lineage_match", "role_in_run",
)


def _global_vcf_cache_path(
    assets: StagedGenomicAssets, vouchers: Sequence[Dict[str, Any]], db_priors: Dict[str, Dict[str, Any]],
    min_hq_reads: int, max_missing_rate: float, paralog_ho_threshold: float, hdplot_d_threshold: float,
) -> Optional[Path]:
    sig_files = [*assets.variant_files[:1], *assets.dosage_files[:1], *assets.stats_files, *assets.key_files, *assets.metadata_csvs, *assets.organellar_files]
    if not (assets.variant_files or assets.dosage_files):
        return None
    h = hashlib.sha256()
    for p in sig_files:
        try:
            st = p.stat()
            h.update(f"{p.resolve()}|{st.st_mtime_ns}|{st.st_size}|".encode("utf-8"))
        except OSError:
            return None
    v_sig = "|".join(f"{v.get('catalogNumber','')}:{v.get('species_standardized','')}" for v in vouchers[:200])
    h.update(f"{min_hq_reads}|{max_missing_rate}|{paralog_ho_threshold}|{hdplot_d_threshold}|{len(vouchers)}|{len(db_priors.get('by_taxon', {}))}|{v_sig}".encode("utf-8"))
    return Path(tempfile.gettempdir()) / "plantmatics_step06_vcf_cache" / f"{h.hexdigest()[:24]}.pkl"


def run_genomic_ingestion(
    genomics_dir: Path,
    out_dir: Path,
    db_path: Optional[Path] = None,
    min_hq_reads: int = 50000,
    max_missing_rate: float = 0.80,
    paralog_ho_threshold: float = 0.60,
    hdplot_d_threshold: float = 4.0,
    demo: bool = False,
    vouchers_csv: Optional[Path] = None,
    focal_taxon: Optional[str] = None,
    intraspecific_mode: bool = False,
    genomic_tier: str = "legacy",
) -> Dict[str, Any]:
    """Ingest staged genomic assets, apply HDplot & taxon-stratified EB dosage estimation, link vouchers, and upsert Tier A SQLite."""
    gen_res, out_res = genomics_dir.resolve(), out_dir.resolve()
    if out_res != gen_res and gen_res in out_res.parents:
        raise ValueError(f"Refusing to write Tier B parsed outputs {out_dir} inside Tier A genomics directory {genomics_dir}")

    assets = discover_genomic_inputs(genomics_dir)
    eff_tier_req = str(genomic_tier or "legacy").strip().lower()
    if eff_tier_req == "legacy" and (vouchers_csv is not None or (focal_taxon and str(focal_taxon).strip())):
        eff_tier_req = "auto"
    has_dyn_inputs = assets.has_genetic_payload or (
        eff_tier_req != "legacy" and (has_valid_public_genomic_records(assets) or bool(focal_taxon and str(focal_taxon).strip()) or (vouchers_csv is not None and Path(vouchers_csv).is_file()))
    )
    if not has_dyn_inputs:
        logger.info("No active genomic payloads found in %s; skipping genomic ingestion (no-op).", genomics_dir)
        return {"status": "skipped_empty", "n_samples_total": 0, "n_samples_pass": 0, "n_loci_raw": 0, "n_loci_retained": 0}

    conn: Optional[sqlite3.Connection] = None
    vouchers: List[Dict[str, Any]] = []
    db_priors = _load_db_cytology_priors(None)
    try:
        conn = get_db_connection(db_path=str(db_path) if db_path else None, read_only=demo, demo=False, init_schema=(not demo), workspace_root=PROJECT_ROOT)
        vouchers = _load_db_vouchers(conn)
        db_priors = _load_db_cytology_priors(conn)
    except Exception as exc:
        logger.warning("SQLite voucher lookup unavailable (%s); continuing with standalone linking.", exc)
        vouchers, db_priors = [], _load_db_cytology_priors(None)

    out_dir.mkdir(parents=True, exist_ok=True)
    dosage_out_path, hdplot_out_path, weights_out_path = out_dir / "allele_dosage_matrix.tsv", out_dir / "hdplot_locus_diagnostics.tsv", out_dir / "locus_weights.tsv"
    exsiccatae_index = _build_exsiccatae_index(vouchers)
    records: List[Dict[str, Any]] = []
    dos_df: Optional[pd.DataFrame] = None
    n_loci_raw, n_loci_retained = 0, 0
    cache_p = _global_vcf_cache_path(assets, vouchers, db_priors, min_hq_reads, max_missing_rate, paralog_ho_threshold, hdplot_d_threshold) if assets.has_genetic_payload else None
    loaded_cache = False
    if cache_p is not None and cache_p.is_file():
        try:
            with cache_p.open("rb") as fh:
                c_obj = pickle.load(fh)
            records = [dict(r) for r in c_obj["records"]]
            dos_df, hd_df, weights_df = c_obj["dos_df"], c_obj["hd_df"], c_obj["weights_df"]
            n_loci_raw, n_loci_retained = int(c_obj["n_loci_raw"]), int(c_obj["n_loci_retained"])
            if hd_df is not None: hd_df.to_csv(hdplot_out_path, sep="\t", index=False, na_rep="NA")
            if weights_df is not None: weights_df.to_csv(weights_out_path, sep="\t", index=False, na_rep="NA")
            loaded_cache = True
        except Exception:
            loaded_cache = False

    if not loaded_cache and assets.has_genetic_payload:
        field_rows, direct_user_meta = _parse_field_metadata(assets.metadata_csvs)
        organellar_map = _load_organellar_metadata(assets.organellar_files)
        stats_map: Dict[str, Dict[str, Any]] = {}
        for sp in assets.stats_files:
            stats_map.update(_parse_stats_file(sp))
        dosage_samples, locus_ids, dosage_mat, ref_depth_mat, alt_depth_mat = [], [], np.empty((0, 0), dtype=np.float32), np.empty((0, 0), dtype=np.float32), np.empty((0, 0), dtype=np.float32)
        from_precomputed_dosage, marker_platform = False, "nextRAD"
        if assets.variant_files:
            primary_vcf = assets.variant_files[0]
            if "hifi" in primary_vcf.name.lower() or "pacbio" in primary_vcf.name.lower(): marker_platform = "PacBio_HiFi"
            elif "ddrad" in primary_vcf.name.lower(): marker_platform = "ddRAD"
            dosage_samples, locus_ids, dosage_mat, ref_depth_mat, alt_depth_mat = _parse_vcf_with_depths(primary_vcf)
        if assets.dosage_files:
            d_sids, d_loci, d_mat = _parse_dosage_table(assets.dosage_files[0])
            if d_sids and d_loci:
                marker_platform = "PacBio_HiFi" if not assets.variant_files else marker_platform
                if ref_depth_mat.size > 0 and alt_depth_mat.size > 0:
                    ref_depth_mat, alt_depth_mat = align_vcf_depths_to_dosage_grid(d_sids, d_loci, dosage_samples, locus_ids, ref_depth_mat, alt_depth_mat)
                else:
                    ref_depth_mat, alt_depth_mat = np.empty((0, 0), dtype=np.float32), np.empty((0, 0), dtype=np.float32)
                dosage_samples, locus_ids, dosage_mat, from_precomputed_dosage = d_sids, d_loci, d_mat, True
        all_sample_ids = list(dict.fromkeys(sid for src in (dosage_samples, list(stats_map.keys()), _parse_key_samples(assets.key_files), list(direct_user_meta.keys())) for sid in src if sid))
        gcfg = resolve_genomics_config() or {}
        min_het_cfg, strong_bic_cfg, min_logl_per_locus_cfg = int(gcfg.get("nquire_min_het_loci", 60)), float(gcfg.get("nquire_strong_delta_bic", 25.0)), float(gcfg.get("nquire_min_delta_logl_per_locus", 0.15))
        sample_idx_map, n_loci_raw = {sid: i for i, sid in enumerate(dosage_samples)}, (int(dosage_mat.shape[1]) if dosage_mat.ndim == 2 else 0)
        ex_prior_lc = {_normalize_exsiccatae_lookup_key(k).lower(): v for k, v in db_priors["by_exsiccatae"].items()}
        fallback_ploidy_by_sid, cyto_ctx_by_sid = {}, {}
        for sid in all_sample_ids:
            rec, fb_pl = _resolve_sample_record(sid, vouchers, exsiccatae_index, field_rows, direct_user_meta, stats_map, organellar_map, db_priors, ex_prior_lc, dosage_mat, ref_depth_mat, alt_depth_mat, sample_idx_map, n_loci_raw, marker_platform, min_hq_reads, max_missing_rate, cyto_ctx_out=cyto_ctx_by_sid)
            records.append(rec); fallback_ploidy_by_sid[sid] = fb_pl
        paraloss_catalog = _load_paraloss1272_catalog(genomics_dir)
        hd_df, weights_df = None, None
        if dosage_samples and n_loci_raw > 0:
            pass_recs = [r for r in records if r["geno_quality_pass"] == 1 and r["genetic_sample_id"] in sample_idx_map]
            pass_ids = [r["genetic_sample_id"] for r in pass_recs]
            ref_pass_ids = [r["genetic_sample_id"] for r in pass_recs if r["is_putative_novel_or_hybrid_query"] == 0]
            pass_rows = [sample_idx_map[s] for s in pass_ids]
            ref_rows = [sample_idx_map[s] for s in ref_pass_ids] if len(ref_pass_ids) >= 4 else pass_rows
            if pass_rows:
                rd_ref = ref_depth_mat[ref_rows, :] if ref_depth_mat.size > 0 else np.empty((0, 0), dtype=np.float32)
                ad_ref = alt_depth_mat[ref_rows, :] if alt_depth_mat.size > 0 else np.empty((0, 0), dtype=np.float32)
                keep_mask, hd_df = _evaluate_hdplot_loci(dosage_mat[pass_rows, :], dosage_mat[ref_rows, :], rd_ref, ad_ref, locus_ids, h_max=paralog_ho_threshold, d_max=hdplot_d_threshold, paraloss_catalog=paraloss_catalog)
                hd_df.to_csv(hdplot_out_path, sep="\t", index=False, na_rep="NA")
                kept_indices = np.where(keep_mask)[0]
                n_loci_retained = int(kept_indices.size)
                kept_loci, kept_dosage = [locus_ids[int(j)] for j in kept_indices], dosage_mat[pass_rows, :][:, kept_indices]
                if ref_depth_mat.size > 0 and alt_depth_mat.size > 0 and n_loci_retained > 0:
                    sc_indices = np.where(keep_mask & (hd_df["is_paraloss1272_single_copy"].to_numpy(dtype=int) == 1))[0]
                    def _fit_cyto_one(r_item: Dict[str, Any]) -> Tuple[Dict[str, Any], str, Any, str, bool]:
                        sid = r_item["genetic_sample_id"]; sidx = sample_idx_map[sid]; use_idx = kept_indices
                        if sc_indices.size >= min_het_cfg:
                            rv_sc, av_sc, dv_sc = ref_depth_mat[sidx, sc_indices], alt_depth_mat[sidx, sc_indices], dosage_mat[sidx, sc_indices]; dp_sc = rv_sc + av_sc
                            hm_sc = np.isfinite(rv_sc) & np.isfinite(av_sc) & np.isfinite(dv_sc) & (dp_sc >= 10.0) & (dv_sc > 0.05) & (dv_sc < 0.95)
                            if np.any(hm_sc) and int(np.sum(((av_sc[hm_sc] / np.maximum(1.0, dp_sc[hm_sc])) >= 0.10) & ((av_sc[hm_sc] / np.maximum(1.0, dp_sc[hm_sc])) <= 0.90))) >= min_het_cfg:
                                use_idx = sc_indices
                        u_pl, sp_m, sp_mu, has_u_cyto = cyto_ctx_by_sid.get(sid, (r_item.get("ploidy_level"), None, False, False))
                        inf_c, nq_d = _infer_sample_cytotype(ref_depth_mat[sidx, use_idx], alt_depth_mat[sidx, use_idx], dosage_mat[sidx, use_idx], u_pl, species_modal_ploidy=sp_m, is_multi_cytotype=sp_mu, min_het_loci=min_het_cfg, strong_delta_bic=strong_bic_cfg, min_delta_logl_per_locus=min_logl_per_locus_cfg)
                        return r_item, sid, inf_c, nq_d, has_u_cyto
                    n_cw = min(resolve_max_workers(), len(pass_recs)) if len(pass_recs) >= 8 else 1
                    cyto_results = list(ThreadPoolExecutor(max_workers=n_cw).map(_fit_cyto_one, pass_recs)) if n_cw > 1 else [_fit_cyto_one(r_item) for r_item in pass_recs]
                    for r_item, sid, inf_c, nq_d, has_u_cyto in cyto_results:
                        if not has_u_cyto: r_item["inferred_cytotype"] = inf_c
                        r_item["nquire_delta_logl"] = nq_d
                        if str(r_item.get("inferred_cytotype") or "") in _CYTO_TO_PLOIDY: fallback_ploidy_by_sid[sid] = _CYTO_TO_PLOIDY[str(r_item["inferred_cytotype"])]
                if not from_precomputed_dosage and ref_depth_mat.size > 0 and n_loci_retained > 0:
                    kept_dosage = _compute_eb_polyrad_dosage(kept_dosage, ref_depth_mat[pass_rows, :][:, kept_indices], alt_depth_mat[pass_rows, :][:, kept_indices], [str(r["inferred_cytotype"]) for r in pass_recs], fallback_ploidies=[fallback_ploidy_by_sid.get(r["genetic_sample_id"]) for r in pass_recs], taxa_labels=[None if r["is_putative_novel_or_hybrid_query"] else r.get("species_standardized") for r in pass_recs])
                for pos_i, r in enumerate(pass_recs):
                    r["hdplot_retained_loci"] = int(np.sum(np.isfinite(kept_dosage[pos_i, :]))) if n_loci_retained > 0 else 0
                dos_df = pd.DataFrame(np.round(kept_dosage, 4), columns=kept_loci); dos_df.insert(0, "genetic_sample_id", pass_ids)
                if paraloss_catalog and any(_lookup_paraloss1272_locus(lid, paraloss_catalog)[0] for lid in locus_ids):
                    weights_df = hd_df.loc[hd_df["hdplot_keep"] == 1, ["locus_id", "locus_paralogy_weight", "is_paraloss1272_single_copy"]]
                    weights_df.to_csv(weights_out_path, sep="\t", index=False, na_rep="NA")
        if cache_p is not None and not demo:
            try:
                cache_p.parent.mkdir(parents=True, exist_ok=True)
                with cache_p.open("wb") as fh:
                    pickle.dump({"records": records, "dos_df": dos_df, "hd_df": hd_df, "weights_df": weights_df, "n_loci_raw": n_loci_raw, "n_loci_retained": n_loci_retained}, fh)
            except Exception:
                pass

    active_tier = "Tier_1_NextRAD_PopGen"
    tier_manifest_path: Optional[Path] = None
    if eff_tier_req != "legacy":
        cohort = resolve_and_stage_dynamic_genomic_cohort(
            genomics_dir=genomics_dir, out_dir=out_dir, assets=assets, global_records=records, global_dosage_df=dos_df,
            vouchers=vouchers, exsiccatae_index=exsiccatae_index, db_priors=db_priors, vouchers_csv=vouchers_csv,
            focal_taxon=focal_taxon, intraspecific_mode=bool(intraspecific_mode), genomic_tier=eff_tier_req,
        )
        records, dos_df = cohort["records"], cohort["dosage_df"]
        active_tier = str(cohort["manifest"]["active_genomic_tier"])
        tier_manifest_path = out_dir / "genomic_tier_manifest.json"
    if dos_df is not None and not dos_df.empty:
        dos_df.to_csv(dosage_out_path, sep="\t", index=False, na_rep="NA")

    samples_csv = out_dir / "parsed_genomic_samples.csv"
    (pd.DataFrame(records) if records else pd.DataFrame(columns=list(_EMPTY_SAMPLE_COLS))).to_csv(samples_csv, index=False)
    manifest = {
        "genomics_dir": str(genomics_dir), "out_dir": str(out_dir), "parsed_samples_csv": str(samples_csv),
        "allele_dosage_tsv": str(dosage_out_path) if dosage_out_path.exists() else None,
        "hdplot_diagnostics_tsv": str(hdplot_out_path) if hdplot_out_path.exists() else None,
        "locus_weights_tsv": str(weights_out_path) if weights_out_path.exists() else None,
        "tier_manifest_json": str(tier_manifest_path) if tier_manifest_path and tier_manifest_path.exists() else None,
        "active_genomic_tier": active_tier,
        "tree_files": [str(p) for p in assets.tree_files], "alignment_files": [str(p) for p in assets.alignment_files],
        "organellar_files": [str(p) for p in assets.organellar_files], "public_reference_files": [str(p) for p in assets.public_reference_files],
        "admixture_q_files": [str(p) for p in assets.admixture_q_files],
        "n_samples_total": len(records), "n_samples_pass": int(sum(1 for r in records if r.get("geno_quality_pass") == 1)),
        "n_loci_raw": n_loci_raw, "n_loci_retained": n_loci_retained,
    }
    (out_dir / "genomic_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    if conn is not None:
        try:
            if not demo and records:
                upsert_genomic_samples(conn, records, preserve_existing_on_null=True)
        finally:
            conn.close()

    logger.info("Ingested %d genomic samples (%d QC-passing, %d/%d singleton loci retained after HDplot H<=%.2f, |D|<=%.1f).", manifest["n_samples_total"], manifest["n_samples_pass"], n_loci_retained, n_loci_raw, paralog_ho_threshold, hdplot_d_threshold)
    return {"status": "ingested", **manifest}


def main(argv: Optional[Sequence[str]] = None) -> int:
    gcfg = resolve_genomics_config()
    parser = argparse.ArgumentParser(description="Ingest staged genomic data into Packera Tier A SQLite.")
    parser.add_argument("--genomics-dir", type=Path, default=PROJECT_ROOT / "data" / "genomics")
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "outputs" / "reports")
    parser.add_argument("--db-path", "--db", dest="db_path", type=Path, default=None)
    parser.add_argument("--vouchers", type=Path, default=None)
    parser.add_argument("--focal-taxon", type=str, default=None)
    parser.add_argument("--intraspecific-mode", action="store_true", default=False)
    parser.add_argument("--genomic-tier", default="legacy", choices=("legacy", "auto", "tier1", "tier2", "tier3", "tier4"))
    parser.add_argument("--min-hq-reads", type=int, default=int(gcfg["min_hq_reads"]))
    parser.add_argument("--max-missing", type=float, default=float(gcfg["max_missing_rate"]))
    parser.add_argument("--paralog-ho", type=float, default=float(gcfg["paralog_ho_threshold"]))
    parser.add_argument("--hdplot-d", type=float, default=float(gcfg["hdplot_d_threshold"]))
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--check-only", action="store_true", help="Print HAS_GENOMIC_INPUTS=0/1 and exit 0 cleanly.")
    parser.add_argument("--check-dynamic", action="store_true", help="Print HAS_DYNAMIC_GENOMIC_INPUTS=0/1 and exit 0 cleanly.")
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.check_only:
        print(f"HAS_GENOMIC_INPUTS={1 if discover_genomic_inputs(args.genomics_dir).has_genetic_payload else 0}")
        return 0
    if args.check_dynamic:
        assets = discover_genomic_inputs(args.genomics_dir)
        ok = assets.has_genetic_payload or has_valid_public_genomic_records(assets) or bool(args.focal_taxon) or bool(args.vouchers and Path(args.vouchers).is_file())
        print(f"HAS_DYNAMIC_GENOMIC_INPUTS={1 if ok else 0}")
        return 0
    run_genomic_ingestion(
        args.genomics_dir, args.out_dir, db_path=args.db_path, min_hq_reads=args.min_hq_reads, max_missing_rate=args.max_missing,
        paralog_ho_threshold=args.paralog_ho, hdplot_d_threshold=args.hdplot_d, demo=bool(args.demo),
        vouchers_csv=args.vouchers, focal_taxon=args.focal_taxon, intraspecific_mode=bool(args.intraspecific_mode), genomic_tier=args.genomic_tier,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
