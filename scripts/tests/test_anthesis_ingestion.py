#!/usr/bin/env python3
"""
===============================================================================
Script: test_anthesis_ingestion.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration test suite validating:
    1. Configuration synchronization for thresholds.capitulum_phenology.
    2. The 3-Tier Anthesis Verification Filter:
       - Tier 1: Calendar gate (50 <= DOY <= 220).
       - Tier 2: Darwin Core text check (fruit/sterile exclusion without 'fl').
       - Tier 3: Computer Vision check (voucher_phenological_state == 'anthesis').
       - Graceful fallback for unreadable crops or missing CV records.
    3. Verification telemetry formatting:
       [INFO] Filtered phenology vouchers: X verified anthesis, Y fruiting-only (excluded), Z sterile (excluded).
    4. Verified Latitudinal Baseline Regression (lm(DOY ~ decimalLatitude)).
    5. Structural syntax and bracket balance of scripts/analysis/06_multimodal_spatial_rf.R.

Execution:
    python -m unittest scripts/tests/test_anthesis_ingestion.py
===============================================================================
"""

from __future__ import annotations

import re
import sys
import unittest
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import PipelineConfig, CapitulumPhenologyConfig


def simulate_3tier_anthesis_filter(df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, int], str]:
    """
    Python reference simulation of the 3-Tier Anthesis Verification Filter
    matching scripts/analysis/06_multimodal_spatial_rf.R.
    """
    df = df.copy()

    # Coordinates
    lat_val = df["decimalLatitude"] if "decimalLatitude" in df.columns else df["latitude"]
    df["decimalLatitude"] = pd.to_numeric(lat_val, errors="coerce")
    df["doy"] = pd.to_numeric(df["doy"], errors="coerce")

    # Tier 1: Calendar gate: 50 <= DOY <= 220 & 24 <= Latitude <= 55
    pass_calendar = (
        df["doy"].between(50, 220)
        & df["decimalLatitude"].between(24.0, 55.0)
    )

    # Tier 2: Darwin Core check
    if "reproductiveCondition" in df.columns:
        rc = df["reproductiveCondition"].fillna("").astype(str).str.lower().str.strip()
        has_fruit = rc.str.contains("fruit", regex=False)
        has_sterile = rc.str.contains("sterile", regex=False)
        has_fl = rc.str.contains("fl", regex=False)
        dwc_fruit_only = has_fruit & (~has_fl)
        dwc_sterile_only = has_sterile & (~has_fl)
    else:
        dwc_fruit_only = pd.Series(False, index=df.index)
        dwc_sterile_only = pd.Series(False, index=df.index)

    pass_dwc = ~(dwc_fruit_only | dwc_sterile_only)

    # Tier 3: Computer Vision check
    if "voucher_phenological_state" in df.columns:
        v_state = df["voucher_phenological_state"].fillna("").astype(str).str.lower().str.strip()
    else:
        v_state = pd.Series("", index=df.index)

    if "total_capitula" in df.columns:
        tot_caps = pd.to_numeric(df["total_capitula"], errors="coerce").fillna(0)
    else:
        tot_caps = pd.Series(0, index=df.index)

    has_valid_cv = (v_state != "") & (~v_state.isin(["unknown", "na"])) & (tot_caps > 0)

    cv_anthesis = pd.Series(False, index=df.index)
    cv_fruit = pd.Series(False, index=df.index)
    cv_sterile = pd.Series(False, index=df.index)

    cv_anthesis[has_valid_cv & (v_state == "anthesis")] = True
    cv_fruit[has_valid_cv & (v_state == "fruit")] = True
    cv_sterile[has_valid_cv & v_state.isin(["sterile", "bud"])] = True

    # Graceful fallback: unreadable crops / missing CV records
    missing_cv = ~has_valid_cv
    cv_anthesis[missing_cv & pass_calendar & pass_dwc] = True
    cv_fruit[missing_cv & dwc_fruit_only] = True
    cv_sterile[missing_cv & dwc_sterile_only] = True

    is_verified_anthesis = pass_calendar & pass_dwc & cv_anthesis
    is_fruit_excluded = pass_calendar & (cv_fruit | dwc_fruit_only) & (~is_verified_anthesis)
    is_sterile_excluded = pass_calendar & (cv_sterile | dwc_sterile_only) & (~is_verified_anthesis) & (~is_fruit_excluded)

    df["is_verified_anthesis"] = is_verified_anthesis
    df["is_fruit_excluded"] = is_fruit_excluded
    df["is_sterile_excluded"] = is_sterile_excluded

    n_verified = int(is_verified_anthesis.sum())
    n_fruit = int(is_fruit_excluded.sum())
    n_sterile = int(is_sterile_excluded.sum())

    counts = {
        "verified_anthesis": n_verified,
        "fruiting_only": n_fruit,
        "sterile": n_sterile,
    }

    telemetry = (
        f"[INFO] Filtered phenology vouchers: {n_verified} verified anthesis, "
        f"{n_fruit} fruiting-only (excluded), {n_sterile} sterile (excluded)."
    )

    return df, counts, telemetry


class TestAnthesisIngestion(unittest.TestCase):
    """Test suite for anthesis ingestion, 3-tier filtering, and latitudinal cline modeling."""

    def test_config_capitulum_phenology_synchronization(self):
        """Verify that PipelineConfig and config.yaml properly synchronize capitulum_phenology."""
        cfg = PipelineConfig.from_yaml()
        self.assertTrue(hasattr(cfg.thresholds, "capitulum_phenology"))
        pheno_cfg = cfg.thresholds.capitulum_phenology

        self.assertIsInstance(pheno_cfg, CapitulumPhenologyConfig)
        self.assertEqual(pheno_cfg.min_pappus_luminance, 75.0)
        self.assertEqual(pheno_cfg.max_pappus_saturation, 0.18)
        self.assertEqual(pheno_cfg.min_pappus_area_ratio, 0.18)
        self.assertEqual(pheno_cfg.min_floret_yellow_saturation, 0.32)
        self.assertEqual(pheno_cfg.min_floret_yellow_ratio, 0.12)
        self.assertEqual(pheno_cfg.min_anthesis_heads_for_flowering, 1)

    def test_calendar_gate_tier1(self):
        """Verify Tier 1 calendar gate bounds (50 <= DOY <= 220)."""
        mock_data = pd.DataFrame({
            "catalogNumber": ["V01", "V02", "V03", "V04"],
            "decimalLatitude": [36.0, 36.0, 36.0, 36.0],
            "doy": [45.0, 50.0, 220.0, 225.0],
            "voucher_phenological_state": ["anthesis", "anthesis", "anthesis", "anthesis"],
            "total_capitula": [5, 5, 5, 5],
        })
        filtered_df, counts, _ = simulate_3tier_anthesis_filter(mock_data)
        self.assertFalse(filtered_df.loc[0, "is_verified_anthesis"])  # DOY 45 < 50
        self.assertTrue(filtered_df.loc[1, "is_verified_anthesis"])   # DOY 50 == 50
        self.assertTrue(filtered_df.loc[2, "is_verified_anthesis"])   # DOY 220 == 220
        self.assertFalse(filtered_df.loc[3, "is_verified_anthesis"])  # DOY 225 > 220
        self.assertEqual(counts["verified_anthesis"], 2)

    def test_darwin_core_tier2(self):
        """Verify Tier 2 Darwin Core check: fruit/sterile exclusion without 'fl'."""
        mock_data = pd.DataFrame({
            "catalogNumber": ["V01", "V02", "V03", "V04", "V05"],
            "decimalLatitude": [36.0, 36.0, 36.0, 36.0, 36.0],
            "doy": [120.0, 120.0, 120.0, 120.0, 120.0],
            "reproductiveCondition": [
                "fl",
                "fl & fruit",
                "fruit",
                "sterile rosette",
                np.nan
            ],
            "voucher_phenological_state": ["anthesis", "anthesis", "anthesis", "anthesis", "anthesis"],
            "total_capitula": [5, 5, 5, 5, 5],
        })
        filtered_df, counts, _ = simulate_3tier_anthesis_filter(mock_data)
        self.assertTrue(filtered_df.loc[0, "is_verified_anthesis"])   # "fl" passes
        self.assertTrue(filtered_df.loc[1, "is_verified_anthesis"])   # "fl & fruit" contains "fl", passes
        self.assertFalse(filtered_df.loc[2, "is_verified_anthesis"])  # "fruit" without "fl" excluded
        self.assertTrue(filtered_df.loc[2, "is_fruit_excluded"])
        self.assertFalse(filtered_df.loc[3, "is_verified_anthesis"])  # "sterile rosette" excluded
        self.assertTrue(filtered_df.loc[3, "is_sterile_excluded"])
        self.assertTrue(filtered_df.loc[4, "is_verified_anthesis"])   # NaN passes DwC

    def test_computer_vision_tier3_and_fallback(self):
        """Verify Tier 3 CV check and graceful fallback for missing CV records."""
        mock_data = pd.DataFrame({
            "catalogNumber": ["V01", "V02", "V03", "V04", "V05"],
            "decimalLatitude": [36.0, 36.0, 36.0, 36.0, 36.0],
            "doy": [120.0, 120.0, 120.0, 120.0, 120.0],
            "voucher_phenological_state": [
                "anthesis",
                "fruit",
                "sterile",
                np.nan,       # Missing CV record
                "unknown",    # Unclassified CV record
            ],
            "total_capitula": [5, 5, 5, 0, 0],
        })
        filtered_df, counts, telemetry = simulate_3tier_anthesis_filter(mock_data)
        self.assertTrue(filtered_df.loc[0, "is_verified_anthesis"])   # CV anthesis
        self.assertFalse(filtered_df.loc[1, "is_verified_anthesis"])  # CV fruit excluded
        self.assertTrue(filtered_df.loc[1, "is_fruit_excluded"])
        self.assertFalse(filtered_df.loc[2, "is_verified_anthesis"])  # CV sterile excluded
        self.assertTrue(filtered_df.loc[2, "is_sterile_excluded"])
        self.assertTrue(filtered_df.loc[3, "is_verified_anthesis"])   # Fallback: calendar + DwC pass
        self.assertTrue(filtered_df.loc[4, "is_verified_anthesis"])   # Fallback: calendar + DwC pass

        self.assertEqual(counts["verified_anthesis"], 3)
        self.assertEqual(counts["fruiting_only"], 1)
        self.assertEqual(counts["sterile"], 1)
        self.assertIn("Filtered phenology vouchers: 3 verified anthesis, 1 fruiting-only (excluded), 1 sterile (excluded).", telemetry)

    def test_latitudinal_cline_fit_on_verified_cohort(self):
        """Verify that fitting lm(DOY ~ decimalLatitude) strictly on verified anthesis yields expected slope."""
        from scipy import stats

        # Synthesize anthesis cohort with ~4.0 days/degree N slope + fruit outliers
        np.random.seed(42)
        lats = np.random.uniform(30.0, 42.0, 50)
        # Anthesis vouchers: DOY = 10 + 3.5 * lat
        anthesis_doy = 10.0 + 3.5 * lats + np.random.normal(0, 2.0, 50)
        # Late fruiting vouchers: DOY = 180 + 1.0 * lat (distorts regression if included)
        fruit_lats = np.random.uniform(30.0, 42.0, 20)
        fruit_doy = 180.0 + 1.0 * fruit_lats

        cohort_df = pd.DataFrame({
            "catalogNumber": [f"ANTH_{i:02d}" for i in range(50)] + [f"FRUIT_{i:02d}" for i in range(20)],
            "decimalLatitude": np.concatenate([lats, fruit_lats]),
            "doy": np.concatenate([anthesis_doy, fruit_doy]),
            "voucher_phenological_state": ["anthesis"] * 50 + ["fruit"] * 20,
            "total_capitula": [4] * 70,
        })

        filtered_df, _, _ = simulate_3tier_anthesis_filter(cohort_df)
        verified = filtered_df[filtered_df["is_verified_anthesis"]]

        # Unfiltered regression slope (distorted by fruiting vouchers)
        unfiltered_reg = stats.linregress(cohort_df["decimalLatitude"], cohort_df["doy"])
        # Verified anthesis regression slope (clean Hopkins' Law progression)
        verified_reg = stats.linregress(verified["decimalLatitude"], verified["doy"])

        self.assertEqual(len(verified), 50)
        self.assertAlmostEqual(verified_reg.slope, 3.5, delta=0.5)
        self.assertGreater(verified_reg.rvalue ** 2, unfiltered_reg.rvalue ** 2)

    def test_r_script_syntax_and_balance(self):
        """Verify that 06_multimodal_spatial_rf.R has balanced braces, brackets, and contains anthesis filter."""
        r_file = PROJECT_ROOT / "scripts" / "analysis" / "06_multimodal_spatial_rf.R"
        self.assertTrue(r_file.exists(), f"Missing R script at {r_file}")

        content = r_file.read_text(encoding="utf-8")

        # Check keyword integration
        self.assertIn("--pheno-states", content)
        self.assertIn("voucher_phenological_states.csv", content)
        self.assertIn("Filtered phenology vouchers:", content)
        self.assertIn("verified anthesis", content)
        self.assertIn("fruiting-only (excluded)", content)
        self.assertIn("is_verified_anthesis", content)
        self.assertIn("model_latitudinal_spring_baseline", content)

        # Robust parser that strips R comments and string literals
        out = []
        in_str = False
        str_char = None
        escape = False
        in_comment = False

        for ch in content:
            if in_comment:
                if ch == "\n":
                    in_comment = False
                    out.append(ch)
                continue

            if in_str:
                if escape:
                    escape = False
                elif ch == "\\":
                    escape = True
                elif ch == str_char:
                    in_str = False
                continue

            if ch == "#":
                in_comment = True
                continue
            elif ch in ('"', "'"):
                in_str = True
                str_char = ch
                continue
            out.append(ch)

        stripped = "".join(out)

        open_curly = stripped.count("{")
        close_curly = stripped.count("}")
        open_paren = stripped.count("(")
        close_paren = stripped.count(")")
        open_bracket = stripped.count("[")
        close_bracket = stripped.count("]")

        self.assertEqual(open_curly, close_curly, f"Curly brace mismatch: {open_curly} vs {close_curly}")
        self.assertEqual(open_paren, close_paren, f"Parenthesis mismatch: {open_paren} vs {close_paren}")
        self.assertEqual(open_bracket, close_bracket, f"Square bracket mismatch: {open_bracket} vs {close_bracket}")


if __name__ == "__main__":
    unittest.main()
