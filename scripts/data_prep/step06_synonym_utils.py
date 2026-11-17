#!/usr/bin/env python3
"""
Module: scripts/data_prep/step06_synonym_utils.py
Project: Plantmatics - Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  Canonical Senecio-to-Packera synonymy, orthographic equivalence keys, and
  FASTA header binomial extractors shared across Step 06 ingestion and
  organellar/public-sequence modules.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Optional, Tuple

CANONICAL_SENECIO_PACKERA_SYNONYMS: Dict[str, str] = {
    "senecio tomentosus": "Packera dubia",
    "packera tomentosa": "Packera dubia",
    "p. tomentosa": "Packera dubia",
    "senecio aureus": "Packera aurea",
    "senecio anonymus": "Packera anonyma",
    "senecio smallii": "Packera anonyma",
    "senecio pauperculus": "Packera paupercula",
    "senecio schweinitzianus": "Packera schweinitziana",
    "senecio robbinsii": "Packera schweinitziana",
    "senecio millefolium": "Packera millefolium",
    "packera millefolia": "Packera millefolium",
    "packera memmingeri": "Packera x memmingeri",
    "packera xmemmingeri": "Packera x memmingeri",
    "packera x memmingeri": "Packera x memmingeri",
    "packera \u00d7memmingeri": "Packera x memmingeri",
    "packera \u00d7 memmingeri": "Packera x memmingeri",
    "senecio obovatus": "Packera obovata",
    "senecio plattensis": "Packera plattensis",
    "senecio antennariifolius": "Packera antennariifolia",
    "senecio glabellus": "Packera glabella",
    "senecio crawfordii": "Packera crawfordii",
    "senecio loratifolius": "Packera loratifolia",
    "senecio pseudaureus": "Packera pseudaurea",
    "senecio indecorus": "Packera indecora",
    "senecio canus": "Packera cana",
    "senecio strepthanthifolius": "Packera streptanthifolia",
    "senecio streptanthifolius": "Packera streptanthifolia",
    "senecio fendleri": "Packera fendleri",
    "senecio multilobatus": "Packera multilobata",
    "senecio neomexicanus": "Packera neomexicana",
    "senecio debilis": "Packera debilis",
    "senecio werneriifolius": "Packera werneriifolia",
    "senecio hartianus": "Packera hartiana",
}

_NON_EPITHET_TOKENS = frozenset({
    "rbcl", "trnl", "its", "its1", "its2", "58s", "waxy", "gbssi", "core",
    "plastid", "chloroplast", "plastome", "ribosomal", "nuclear", "sp",
    "unk", "unknown", "var", "subsp", "cf", "aff",
})

_CANON_EQUIV_OVERRIDES: Dict[str, str] = {
    "packera millefolia": "packera millefolium",
    "packera millefolium": "packera millefolium",
    "packera memmingeri": "packera xmemmingeri",
    "packera x memmingeri": "packera xmemmingeri",
    "packera xmemmingeri": "packera xmemmingeri",
    "packera tomentosa": "packera dubia",
    "p. tomentosa": "packera dubia",
    "senecio tomentosus": "packera dubia",
    "packera dubia": "packera dubia",
}


def _resolve_active_genus(target_genus: Optional[str] = None) -> str:
    if target_genus and str(target_genus).strip():
        return str(target_genus).strip().capitalize()
    try:
        from scripts.core.config import resolve_target_genus
        if g := str(resolve_target_genus() or "").strip():
            return g.capitalize()
    except Exception:
        pass
    return "Packera"


def canonical_taxon_equiv_key(taxon: Any, target_genus: Optional[str] = None) -> str:
    """Normalize whitespace, underscores, hybrid signs, and synonyms to a canonical lowercase comparison key."""
    if taxon is None:
        return ""
    if isinstance(taxon, float) and taxon != taxon:
        return ""
    raw = re.sub(r"[\s_]+", " ", str(taxon)).strip().lower().replace("\u00d7", "x")
    if not raw or raw in {"nan", "none", "null"}:
        return ""
    if _resolve_active_genus(target_genus).lower() == "packera":
        raw = re.sub(r"\bpackera\s+x\s*memmingeri\b", "packera xmemmingeri", raw)
        if raw in _CANON_EQUIV_OVERRIDES:
            return _CANON_EQUIV_OVERRIDES[raw]
        if raw in CANONICAL_SENECIO_PACKERA_SYNONYMS:
            canon = CANONICAL_SENECIO_PACKERA_SYNONYMS[raw].lower()
            return "packera xmemmingeri" if canon == "packera x memmingeri" else canon
    return raw


def normalize_step06_taxon_name(
    raw_taxon: Optional[str], target_genus: Optional[str] = None,
) -> Tuple[Optional[str], Optional[str], bool]:
    """Return `(taxon_raw, species_standardized, is_novel_or_hybrid)` using active genus synonymy."""
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

    active_genus = _resolve_active_genus(target_genus)
    is_packera = active_genus.lower() == "packera"
    if is_packera and low in CANONICAL_SENECIO_PACKERA_SYNONYMS:
        return t_raw, CANONICAL_SENECIO_PACKERA_SYNONYMS[low], False

    try:
        import scripts.core.harvester as _harv
        std_fn = getattr(_harv, "standardize_taxon_name", None)
        if callable(std_fn):
            std = std_fn(t_raw, target_genus=active_genus)
        elif is_packera:
            std = _harv.standardize_packera_taxon(t_raw)
        elif (m := re.match(rf"^{re.escape(active_genus.lower())}\s+([a-z][a-z-]*)\b", low)) and m.group(1) not in _NON_EPITHET_TOKENS:
            std = f"{active_genus} {m.group(1)}"
        else:
            std = None
        if std:
            return t_raw, std, False
    except Exception:
        pass
    return t_raw, t_raw, False


def normalize_public_header_binomial(text: str, target_genus: Optional[str] = None) -> str:
    """Normalize a FASTA header or taxon string into `<Genus>_<epithet>` (or `<Genus>_sp`)."""
    active_genus = _resolve_active_genus(target_genus)
    g_low = active_genus.lower()
    hay = str(text or "").replace("_", " ").lower()
    if g_low == "packera":
        for syn, canon in CANONICAL_SENECIO_PACKERA_SYNONYMS.items():
            if syn in hay:
                return "Packera_xmemmingeri" if canon == "Packera x memmingeri" else canon.replace(" ", "_")
    if m := re.search(rf"\b{re.escape(g_low)}\s+[0-9]*([a-z-]+)\b", hay):
        ep = m.group(1).strip("-")
        if g_low == "packera":
            if ep == "tomentosa":
                return "Packera_dubia"
            if ep == "millefolia":
                return "Packera_millefolium"
            if ep in ("memmingeri", "xmemmingeri"):
                return "Packera_xmemmingeri"
        return f"{active_genus}_{ep}"
    return f"{active_genus}_sp"


def extract_species_from_seed_key(seed_key: str, target_genus: Optional[str] = None) -> str:
    """Extract canonical `<Genus> <epithet>` from a reference FASTA header without accession hardcoding."""
    active_genus = _resolve_active_genus(target_genus)
    g_low = active_genus.lower()
    norm = re.sub(r"[|_/\t]+", " ", str(seed_key or "")).lower()
    if g_low == "packera":
        for syn, canon in CANONICAL_SENECIO_PACKERA_SYNONYMS.items():
            if syn in norm:
                return canon
    if m := re.search(rf"\b{re.escape(g_low)}\s+[0-9]*([a-z]+(?:-[a-z]+)?)(?=\s|$|[0-9])", norm):
        ep = m.group(1).strip("-")
        if ep and ep not in _NON_EPITHET_TOKENS:
            if g_low == "packera":
                if ep == "tomentosa":
                    return "Packera dubia"
                if ep == "millefolia":
                    return "Packera millefolium"
                if ep in ("memmingeri", "xmemmingeri"):
                    return "Packera x memmingeri"
            return f"{active_genus} {ep}"
    return ""

