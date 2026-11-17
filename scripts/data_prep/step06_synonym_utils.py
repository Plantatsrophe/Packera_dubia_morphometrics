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


def canonical_taxon_equiv_key(taxon: Any) -> str:
    """Normalize whitespace, underscores, hybrid signs, and synonyms to a canonical lowercase comparison key."""
    if taxon is None:
        return ""
    if isinstance(taxon, float) and taxon != taxon:
        return ""
    raw = re.sub(r"[\s_]+", " ", str(taxon)).strip().lower().replace("\u00d7", "x")
    if not raw or raw in {"nan", "none", "null"}:
        return ""
    raw = re.sub(r"\bpackera\s+x\s*memmingeri\b", "packera xmemmingeri", raw)
    if raw in _CANON_EQUIV_OVERRIDES:
        return _CANON_EQUIV_OVERRIDES[raw]
    if raw in CANONICAL_SENECIO_PACKERA_SYNONYMS:
        canon = CANONICAL_SENECIO_PACKERA_SYNONYMS[raw].lower()
        return "packera xmemmingeri" if canon == "packera x memmingeri" else canon
    return raw


def normalize_step06_taxon_name(raw_taxon: Optional[str]) -> Tuple[Optional[str], Optional[str], bool]:
    """Return `(taxon_raw, species_standardized, is_novel_or_hybrid)` using canonical Packera synonymy."""
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

    if low in CANONICAL_SENECIO_PACKERA_SYNONYMS:
        return t_raw, CANONICAL_SENECIO_PACKERA_SYNONYMS[low], False

    try:
        from scripts.core.harvester import standardize_packera_taxon

        std = standardize_packera_taxon(t_raw)
        if std:
            return t_raw, std, False
    except Exception:
        pass
    return t_raw, t_raw, False


def normalize_public_header_binomial(text: str) -> str:
    """Normalize a FASTA header or taxon string into `Packera_<epithet>` (or `Packera_sp`)."""
    hay = str(text or "").replace("_", " ").lower()
    for syn, canon in CANONICAL_SENECIO_PACKERA_SYNONYMS.items():
        if syn in hay:
            return "Packera_xmemmingeri" if canon == "Packera x memmingeri" else canon.replace(" ", "_")
    if m := re.search(r"\bpackera\s+[0-9]*([a-z-]+)\b", hay):
        ep = m.group(1).strip("-")
        if ep == "tomentosa":
            return "Packera_dubia"
        if ep == "millefolia":
            return "Packera_millefolium"
        if ep in ("memmingeri", "xmemmingeri"):
            return "Packera_xmemmingeri"
        return f"Packera_{ep}"
    return "Packera_sp"


def extract_species_from_seed_key(seed_key: str) -> str:
    """Extract canonical `Packera <epithet>` from a reference FASTA header without accession hardcoding."""
    norm = re.sub(r"[|_/\t]+", " ", str(seed_key or "")).lower()
    for syn, canon in CANONICAL_SENECIO_PACKERA_SYNONYMS.items():
        if syn in norm:
            return canon
    if m := re.search(r"\bpackera\s+[0-9]*([a-z]+(?:-[a-z]+)?)(?=\s|$|[0-9])", norm):
        ep = m.group(1).strip("-")
        if ep and ep not in _NON_EPITHET_TOKENS:
            if ep == "tomentosa":
                return "Packera dubia"
            if ep == "millefolia":
                return "Packera millefolium"
            if ep in ("memmingeri", "xmemmingeri"):
                return "Packera x memmingeri"
            return f"Packera {ep}"
    return ""

