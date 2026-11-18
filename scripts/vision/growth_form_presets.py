#!/usr/bin/env python3
"""
===============================================================================
Module: growth_form_presets.py
Project: Plantmatics Genus-Agnostic Morphometrics & Specimen Triage Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Declarative growth-form preset registry for Step 02 LeafMachine3 organ
    clustering, leaflet promotion, non-laminate bypass gating, and optical
    reproductive phenology plugin routing.
===============================================================================
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


@dataclass(frozen=True)
class GrowthFormSpec:
    """Declarative behavior specification for a botanical growth-form preset."""

    name: str
    consolidate_rosette_mm: Optional[float]
    allow_cauline_primary: bool
    treat_leaflet_as_leaf_whole: bool
    downweight_efa_outline: bool
    bypass_planar_leaf_gate: bool
    phenology_plugin: str


GROWTH_FORM_PRESETS: Dict[str, GrowthFormSpec] = {
    "rosette_herb": GrowthFormSpec(
        name="rosette_herb",
        consolidate_rosette_mm=35.0,
        allow_cauline_primary=False,
        treat_leaflet_as_leaf_whole=False,
        downweight_efa_outline=False,
        bypass_planar_leaf_gate=False,
        phenology_plugin="auto",
    ),
    "cauline_herb": GrowthFormSpec(
        name="cauline_herb",
        consolidate_rosette_mm=20.0,
        allow_cauline_primary=True,
        treat_leaflet_as_leaf_whole=False,
        downweight_efa_outline=False,
        bypass_planar_leaf_gate=False,
        phenology_plugin="generic_lm3",
    ),
    "woody_branch": GrowthFormSpec(
        name="woody_branch",
        consolidate_rosette_mm=None,
        allow_cauline_primary=True,
        treat_leaflet_as_leaf_whole=False,
        downweight_efa_outline=False,
        bypass_planar_leaf_gate=False,
        phenology_plugin="generic_lm3",
    ),
    "compound_leaf": GrowthFormSpec(
        name="compound_leaf",
        consolidate_rosette_mm=None,
        allow_cauline_primary=True,
        treat_leaflet_as_leaf_whole=True,
        downweight_efa_outline=False,
        bypass_planar_leaf_gate=False,
        phenology_plugin="generic_lm3",
    ),
    "graminoid_linear": GrowthFormSpec(
        name="graminoid_linear",
        consolidate_rosette_mm=25.0,
        allow_cauline_primary=True,
        treat_leaflet_as_leaf_whole=False,
        downweight_efa_outline=True,
        bypass_planar_leaf_gate=False,
        phenology_plugin="generic_lm3",
    ),
    "non_laminate_bypass": GrowthFormSpec(
        name="non_laminate_bypass",
        consolidate_rosette_mm=None,
        allow_cauline_primary=True,
        treat_leaflet_as_leaf_whole=False,
        downweight_efa_outline=True,
        bypass_planar_leaf_gate=True,
        phenology_plugin="generic_lm3",
    ),
}

DEFAULT_GROWTH_FORM_PRESET: str = "rosette_herb"


def _resolve_preset_from_config(cfg: Optional[Any] = None) -> Optional[str]:
    """Resolve growth-form preset name from core config when available."""
    try:
        from scripts.core.config import resolve_growth_form_preset

        try:
            raw = resolve_growth_form_preset(cfg=cfg)
        except TypeError:
            raw = resolve_growth_form_preset()
        val = str(raw or "").strip()
        if val:
            return val
    except ImportError:
        pass
    except Exception:
        pass

    if cfg is not None:
        for sec_name in ("taxa", "segmentation"):
            sec = getattr(cfg, sec_name, None)
            if sec is not None:
                val = str(getattr(sec, "growth_form_preset", "") or "").strip()
                if val:
                    return val
    return None


def get_growth_form_spec(
    preset_name: Optional[str] = None,
    cfg: Optional[Any] = None,
) -> GrowthFormSpec:
    """
    Resolve a `GrowthFormSpec` from an explicit name, `PLANTMATICS_GROWTH_FORM_PRESET`,
    core `resolve_growth_form_preset()`, or the `"rosette_herb"` default.

    Raises `ValueError` for unknown preset names.
    """
    if preset_name is not None:
        key = str(preset_name).strip().lower()
        if not key or key not in GROWTH_FORM_PRESETS:
            raise ValueError(
                f"Unknown growth-form preset {preset_name!r}; "
                f"expected one of {sorted(GROWTH_FORM_PRESETS)}."
            )
        return GROWTH_FORM_PRESETS[key]

    env_val = os.environ.get("PLANTMATICS_GROWTH_FORM_PRESET", "").strip()
    if env_val:
        key = env_val.lower()
        if key not in GROWTH_FORM_PRESETS:
            raise ValueError(
                f"Unknown PLANTMATICS_GROWTH_FORM_PRESET {env_val!r}; "
                f"expected one of {sorted(GROWTH_FORM_PRESETS)}."
            )
        return GROWTH_FORM_PRESETS[key]

    cfg_val = _resolve_preset_from_config(cfg=cfg)
    if cfg_val:
        key = cfg_val.lower()
        if key not in GROWTH_FORM_PRESETS:
            raise ValueError(
                f"Unknown configured growth-form preset {cfg_val!r}; "
                f"expected one of {sorted(GROWTH_FORM_PRESETS)}."
            )
        return GROWTH_FORM_PRESETS[key]

    return GROWTH_FORM_PRESETS[DEFAULT_GROWTH_FORM_PRESET]


def resolve_effective_target_family(
    target_family: Optional[str] = None,
    cfg: Optional[Any] = None,
) -> str:
    """Resolve target botanical family cleanly even when vision modules run standalone."""
    if target_family is not None and str(target_family).strip():
        return str(target_family).strip()
    env_fam = os.environ.get("PLANTMATICS_TARGET_FAMILY", "").strip()
    if env_fam:
        return env_fam
    try:
        from scripts.core.config import resolve_target_family

        val = str(resolve_target_family(cfg=cfg) or "").strip()
        if val:
            return val
    except ImportError:
        pass
    except Exception:
        pass
    if cfg is not None:
        taxa_sec = getattr(cfg, "taxa", None)
        if taxa_sec is not None:
            val = str(getattr(taxa_sec, "target_family", "") or "").strip()
            if val:
                return val
    return "Asteraceae"


__all__ = [
    "DEFAULT_GROWTH_FORM_PRESET",
    "GROWTH_FORM_PRESETS",
    "GrowthFormSpec",
    "get_growth_form_spec",
    "resolve_effective_target_family",
]
