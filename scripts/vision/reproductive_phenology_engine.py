#!/usr/bin/env python3
"""
===============================================================================
Module: reproductive_phenology_engine.py
Project: Plantmatics Genus-Agnostic Morphometrics & Specimen Triage Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Universal optical reproductive phenology engine routing across:
      1. `asteraceae_radiate` -> `capitulum_phenology_classifier.py` (preserving
         peduncle-oriented yellow/aged-ochre corolla & pappus plume biomarkers).
      2. `generic_lm3` -> family-agnostic reproductive classification directly
         from LeafMachine3 / detector organ classes and labels without penalizing
         non-yellow corollas or absence of pappus bristles.
===============================================================================
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.capitulum_phenology_classifier import (
    VOUCHER_ANTHESIS,
    VOUCHER_BUD,
    VOUCHER_FRUIT,
    VOUCHER_STERILE,
    VoucherPhenologyResult,
    classify_asteraceae_voucher_crops,
)

_NON_REPRO_TOKENS = frozenset({
    "leaf", "leaf_whole", "leaf_partial", "leaflet", "leaflet_whole", "leaflet_one",
    "roots", "root", "wood", "stem", "petiole", "hole", "specimen", "ruler", "label",
})
_SENESCENT_TOKENS = ("senescent", "senescence", "withered")
_BUD_TOKENS = ("bud_whole", "bud_one", "bud", "buds", "flower_bud", "floral_bud")
_ANTHESIS_TOKENS = (
    "flower_whole", "flower_one", "flower_many", "flower", "flowers",
    "anthesis", "catkin", "catkins",
)
_FRUIT_TOKENS = (
    "fruit_whole", "fruit_one", "fruit_many", "seed_fruit_one", "seed_fruit_many",
    "seed_fruit", "fruit", "fruits", "acorn", "capsule", "achene", "pod",
    "samara", "silique", "follicle", "nutlet", "cone", "drupe", "berry", "seed",
)
_LM3_CLASS_STAGE_MAP: Dict[int, Optional[str]] = {
    0: None, 1: None, 2: "fruit", 3: "fruit",
    4: "anthesis", 5: "anthesis", 6: "bud",
    7: None, 8: None, 9: None, 10: None,
}
_LABEL_KEYS = ("label", "class_name", "name", "class", "category", "organ_class", "cls_name")
_CID_KEYS = ("class_id", "cls_id", "cls", "category_id")
_CONF_KEYS = ("confidence", "score", "conf", "prob")


def resolve_phenology_plugin(
    phenology_plugin: str = "generic_lm3",
    target_family: Optional[str] = None,
) -> str:
    """Resolve effective phenology plugin (`asteraceae_radiate` vs `generic_lm3`)."""
    fam_norm = str(target_family or "").strip().lower()
    plug_norm = str(phenology_plugin or "generic_lm3").strip().lower()
    if plug_norm == "auto":
        return "asteraceae_radiate" if fam_norm == "asteraceae" else "generic_lm3"
    if fam_norm == "asteraceae" or plug_norm in ("asteraceae_radiate", "asteraceae", "capitulum"):
        return "asteraceae_radiate"
    return "generic_lm3"


def _first_field(obj: Any, keys: Tuple[str, ...]) -> Any:
    """Return the first non-None mapping entry or attribute matching `keys`."""
    is_map = isinstance(obj, Mapping)
    for k in keys:
        val = obj.get(k) if is_map else getattr(obj, k, None)
        if val is not None:
            return val
    return None


def _coerce_items_and_catalog(
    voucher_crops_or_boxes: Any,
    catalog_number: str = "",
) -> Tuple[List[Any], str]:
    """Normalize input crops/boxes/dicts/tuples into a list of detection items."""
    cat = str(catalog_number or "").strip()
    if voucher_crops_or_boxes is None:
        return [], cat

    if isinstance(voucher_crops_or_boxes, Mapping):
        if not cat:
            cat = str(_first_field(voucher_crops_or_boxes, ("catalogNumber", "catalog_number")) or "").strip()
        for key in ("detections", "boxes", "crops", "items", "organs", "head_classifications"):
            val = voucher_crops_or_boxes.get(key)
            if val is not None and not isinstance(val, (str, bytes)):
                return list(val), cat
        return [voucher_crops_or_boxes], cat

    if isinstance(voucher_crops_or_boxes, np.ndarray):
        return ([voucher_crops_or_boxes] if voucher_crops_or_boxes.ndim == 3 else []), cat

    if isinstance(voucher_crops_or_boxes, tuple):
        is_single_det = len(voucher_crops_or_boxes) in (2, 3, 4) and (
            isinstance(voucher_crops_or_boxes[0], np.ndarray)
            or (
                isinstance(voucher_crops_or_boxes[0], (tuple, list))
                and len(voucher_crops_or_boxes[0]) == 4
                and all(isinstance(v, (int, float, np.number)) for v in voucher_crops_or_boxes[0])
            )
        )
        return ([voucher_crops_or_boxes] if is_single_det else list(voucher_crops_or_boxes)), cat

    if isinstance(voucher_crops_or_boxes, str):
        return ([voucher_crops_or_boxes] if voucher_crops_or_boxes.strip() else []), cat

    return list(voucher_crops_or_boxes), cat


def _classify_label_string(raw_label: Any) -> Optional[str]:
    """Map an organ label string to `anthesis`, `fruit`, `bud`, `senescent`, or None."""
    if raw_label is None:
        return None
    norm = re.sub(r"[^a-z0-9]+", "_", str(raw_label).strip().lower()).strip("_")
    if not norm or norm in _NON_REPRO_TOKENS:
        return None
    parts = set(norm.split("_"))

    if any(tok in norm for tok in _SENESCENT_TOKENS):
        return "senescent"
    if norm in _BUD_TOKENS or bool(parts & {"bud", "buds"}):
        return "bud"
    if any(tok in norm for tok in _ANTHESIS_TOKENS) or bool(parts & {"flower", "flowers", "anthesis", "catkin", "catkins"}):
        return "anthesis"
    if any(tok in norm for tok in _FRUIT_TOKENS) or bool(parts & {"fruit", "fruits", "acorn", "capsule", "achene", "pod", "samara", "seed"}):
        return "fruit"
    return None


def _extract_generic_stage_and_conf(item: Any) -> Tuple[Optional[str], Optional[float]]:
    """Extract `(stage, confidence)` from a dict, tuple, object, or label string."""
    raw_label = raw_cid = raw_conf = head_state = None

    if isinstance(item, str):
        raw_label = item
    elif isinstance(item, (tuple, list)):
        seq = list(item)
        if len(seq) == 2:
            if isinstance(seq[1], (int, np.integer)) and not isinstance(seq[1], bool):
                raw_cid = int(seq[1])
            else:
                raw_label = seq[1]
        elif len(seq) >= 3:
            c_idx, l_idx = (1, 2) if len(seq) == 3 else (2, 3)
            if isinstance(seq[c_idx], (int, float, np.number)) and not isinstance(seq[c_idx], bool):
                raw_conf = float(seq[c_idx])
            if isinstance(seq[l_idx], (int, np.integer)) and not isinstance(seq[l_idx], bool):
                raw_cid = int(seq[l_idx])
            else:
                raw_label = seq[l_idx]
    elif not isinstance(item, np.ndarray):
        raw_label = _first_field(item, _LABEL_KEYS)
        raw_cid = _first_field(item, _CID_KEYS)
        raw_conf = _first_field(item, _CONF_KEYS)
        head_state = item.get("head_state") if isinstance(item, Mapping) else getattr(item, "head_state", None)
        if head_state is None and isinstance(item, Mapping) and isinstance(item.get("head_pheno"), Mapping):
            head_state = item["head_pheno"].get("head_state")

    conf_val: Optional[float] = None
    if raw_conf is not None:
        try:
            cf = float(raw_conf)
            if math.isfinite(cf):
                if cf < 0.0:
                    return None, None
                conf_val = max(0.0, min(1.0, cf / 100.0 if cf > 1.0 else cf))
        except (TypeError, ValueError):
            conf_val = None

    stage = _classify_label_string(raw_label)
    if stage is None and raw_label is not None:
        norm_l = re.sub(r"[^a-z0-9]+", "_", str(raw_label).strip().lower()).strip("_")
        if norm_l in _NON_REPRO_TOKENS:
            return None, None
    if stage is None and raw_cid is not None:
        try:
            stage = _LM3_CLASS_STAGE_MAP.get(int(raw_cid))
        except (TypeError, ValueError):
            stage = None
    if stage is None and head_state is not None:
        hs = str(head_state).strip().lower()
        if hs in ("anthesis", "fruit", "bud", "senescent"):
            stage = hs
    return stage, conf_val


def classify_voucher_reproductive_phenology(
    voucher_crops_or_boxes: Any,
    sheet_rgb: Optional[Any] = None,
    phenology_plugin: str = "generic_lm3",
    target_family: Optional[str] = None,
    *,
    catalog_number: str = "",
    min_anthesis_heads_for_flowering: int = 1,
    capitulum_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Classify voucher-level reproductive phenology across botanical families.

    Routes Asteraceae (`target_family == "Asteraceae"` or `phenology_plugin == "asteraceae_radiate"`,
    or `"auto"` when `target_family` is `"Asteraceae"`) to `capitulum_phenology_classifier.py`,
    and classifies all other families using family-agnostic `generic_lm3` organ detections/labels.
    """
    items, cat = _coerce_items_and_catalog(voucher_crops_or_boxes, catalog_number=catalog_number)
    resolved_plugin = resolve_phenology_plugin(phenology_plugin=phenology_plugin, target_family=target_family)
    if resolved_plugin == "asteraceae_radiate":
        return classify_asteraceae_voucher_crops(
            items=items,
            sheet_rgb=sheet_rgb,
            catalog_number=cat,
            min_anthesis_heads_for_flowering=min_anthesis_heads_for_flowering,
            capitulum_kwargs=capitulum_kwargs,
        )

    n_anthesis = n_fruiting = n_bud = n_senescent = 0
    confs: List[float] = []
    for item in items:
        stage, conf_val = _extract_generic_stage_and_conf(item)
        if stage == "anthesis":
            n_anthesis += 1
        elif stage == "fruit":
            n_fruiting += 1
        elif stage == "bud":
            n_bud += 1
        elif stage == "senescent":
            n_senescent += 1
        else:
            continue
        if conf_val is not None:
            confs.append(conf_val)

    min_anth = max(1, int(min_anthesis_heads_for_flowering))
    if n_anthesis >= min_anth:
        voucher_state = VOUCHER_ANTHESIS
    elif n_fruiting >= 1:
        voucher_state = VOUCHER_FRUIT
    elif n_bud >= 1:
        voucher_state = VOUCHER_BUD
    else:
        voucher_state = VOUCHER_STERILE

    total_repro = int(n_anthesis + n_fruiting + n_bud + n_senescent)
    repro_conf = round(float(np.mean(confs)), 4) if confs else (1.0 if total_repro > 0 else 0.0)
    return VoucherPhenologyResult({
        "catalogNumber": cat,
        "total_capitula": total_repro,
        "total_reproductive_organs": total_repro,
        "n_anthesis": int(n_anthesis),
        "n_fruiting": int(n_fruiting),
        "n_bud": int(n_bud),
        "n_senescent": int(n_senescent),
        "voucher_phenological_state": voucher_state,
        "is_verified_anthesis": bool(n_anthesis >= 1),
        "reproductive_stage_confidence": repro_conf,
        "dominant_pappus_ratio": 0.0,
        "phenology_plugin": "generic_lm3",
    })


__all__ = ["classify_voucher_reproductive_phenology", "resolve_phenology_plugin"]
