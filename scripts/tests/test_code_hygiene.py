#!/usr/bin/env python3
"""
===============================================================================
Script: test_code_hygiene.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Static AST and text/grep hygiene test suite verifying that all active .py,
    .R, and .sh files in the repository (excluding archived script folders,
    .git, .venv, .venv_LM3, and LeafMachine3) are 100% isolated from legacy
    archived code/data directories and contain zero legacy LeafMachine2 or
    Detectron2 PointRend remnants.
===============================================================================
"""

from __future__ import annotations

import ast
import sys
import unittest
from pathlib import Path
from typing import List

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

ARCHIVE_MARKER = "_" + "archive"
SELF_PATH = Path(__file__).resolve()

EXCLUDED_DIRS = {
    ARCHIVE_MARKER,
    ".git",
    ".venv",
    ".venv_LM3",
    "LeafMachine3",
    "__pycache__",
}

FORBIDDEN_LEGACY_LM2_TOKENS = (
    "LeafMachine2",
    "leafmachine2",
    ".venv_LM2",
    "setup_leafmachine2",
    "configure_leafmachine2",
    "lm2_data_loader",
    "lm2_inference_engines",
    "lm2_geometry_utils",
    "finetune_lm2_pcd",
    "PointRendInferenceEngine",
    "LM2Stage1Detector",
    "detectron2",
)

FORBIDDEN_OBSOLETE_SAM2_AND_PRUNED_HELPER_TOKENS = (
    # SAM-2 & tuning tokens
    "sam2",
    "SAM2",
    "annotate_with_sam2",
    "sam2_annotator_utils",
    "SAM2BoundaryRefiner",
    "SAM2PromptBundle",
    "refine_mask_with_sam2",
    "tune_geometry_parameters",
    "run_calibrate_geometry",
    # Pruned helper tokens
    "detect_ruler_scale_hough",
    "trim_petiole_tail",
    "_sample_cross_sectional_widths",
    "_disambiguate_base_apex_endpoints",
    "homologize_contour_starting_point",
    "_align_and_normalize_contour",
    "align_contour_basal_apical",
    "is_non_lamina_petiole_ribbon",
    # Pre-LeafMachine3 legacy Python QC, mask-bridging & Tier-2 reflection tokens
    "leaf_gatekeeper",
    "LeafRoutingResult",
    "TIER_1_DIRECT",
    "TIER_2_REFLECTED",
    "Tier_1_Direct",
    "Tier_2_Reflected",
    "Tier_1_Pure",
    "Tier_2_Fallback",
    "tier1_intact",
    "tier2_reflected",
    "assigned_tier",
    "tier_assignment",
    "ucs_score",
    "min_ucs",
    "min_solidity",
    "enable_geometric_qc",
    "is_botanical_dissection",
    "compute_geometric_metrics",
    "compute_classical_shape_descriptors",
    "_prepare_lamina_mask",
    "extract_tier1_pristine",
    "compute_sheet_adaptive_min_blade_area_cm2",
    "FoldDetectionConfig",
    "DissectionConfig",
    "n_tier1_leaves",
    "n_tier2_leaves",
    "is_folded",
    "reflection_applied",
    "render_tier1_vs_tier2_2d_pdf",
    "tier_symmetry_validation",
    "tier1_vs_tier2_density_overlay",
    "find_lamina_midrib_endpoints",
    # LeafMachine2 tokens
    *FORBIDDEN_LEGACY_LM2_TOKENS,
)


def iter_active_source_files(include_shell: bool = False) -> List[Path]:
    """Return all active .py, .R (and optionally .sh) files under scripts/ and root, excluding archive and submodule/venv directories."""
    exts = ("*.py", "*.R", "*.sh") if include_shell else ("*.py", "*.R")
    files: List[Path] = []
    for ext in exts:
        for root_file in sorted(PROJECT_ROOT.glob(ext)):
            resolved = root_file.resolve()
            if resolved.is_file():
                files.append(resolved)

    scripts_dir = PROJECT_ROOT / "scripts"
    for ext in exts:
        for path in sorted(scripts_dir.rglob(ext)):
            resolved = path.resolve()
            if not resolved.is_file():
                continue
            rel_parts = set(resolved.relative_to(PROJECT_ROOT).parts)
            if rel_parts & EXCLUDED_DIRS:
                continue
            files.append(resolved)
    return files


class TestCodeHygiene(unittest.TestCase):
    """Static AST and grep hygiene checks for archive isolation and zero legacy LM2 remnants across active files."""

    def test_no_active_files_reference_legacy_folders_via_grep(self):
        """Scan all active .py and .R files and fail if any references the archive folder."""
        violations: List[str] = []
        for src_path in iter_active_source_files():
            content = src_path.read_text(encoding="utf-8", errors="ignore")
            for lineno, line in enumerate(content.splitlines(), start=1):
                if ARCHIVE_MARKER in line:
                    rel = src_path.relative_to(PROJECT_ROOT)
                    violations.append(f"{rel}:{lineno}: {line.strip()}")

        self.assertEqual(
            violations,
            [],
            f"Active source files must not reference '{ARCHIVE_MARKER}':\n" + "\n".join(violations),
        )

    def test_no_active_python_ast_imports_or_constants_reference_legacy_folders(self):
        """Parse AST of all active .py files and verify zero imports or string literals reference archive."""
        ast_violations: List[str] = []
        for src_path in iter_active_source_files():
            if src_path.suffix != ".py":
                continue
            source = src_path.read_text(encoding="utf-8")
            tree = ast.parse(source, filename=str(src_path))
            rel = src_path.relative_to(PROJECT_ROOT)

            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        if ARCHIVE_MARKER in alias.name:
                            ast_violations.append(f"{rel}:{node.lineno} import {alias.name}")
                elif isinstance(node, ast.ImportFrom):
                    mod_name = node.module or ""
                    if ARCHIVE_MARKER in mod_name:
                        ast_violations.append(f"{rel}:{node.lineno} from {mod_name} import ...")
                elif isinstance(node, ast.Constant) and isinstance(node.value, str):
                    if ARCHIVE_MARKER in node.value:
                        ast_violations.append(f"{rel}:{node.lineno} string constant '{node.value}'")

        self.assertEqual(
            ast_violations,
            [],
            f"AST inspection detected '{ARCHIVE_MARKER}' references in active Python files:\n"
            + "\n".join(ast_violations),
        )

    def test_no_active_files_reference_leafmachine2_or_pointrend(self):
        """Scan all active .py, .R, and .sh files (excluding test_code_hygiene.py, archived folders, .git, .venv, .venv_LM3, LeafMachine3) and assert zero forbidden LM2/PointRend tokens."""
        violations: List[str] = []
        for src_path in iter_active_source_files(include_shell=True):
            if src_path == SELF_PATH:
                continue
            content = src_path.read_text(encoding="utf-8", errors="ignore")
            rel = src_path.relative_to(PROJECT_ROOT)
            for lineno, line in enumerate(content.splitlines(), start=1):
                for token in FORBIDDEN_LEGACY_LM2_TOKENS:
                    if token in line:
                        violations.append(f"{rel}:{lineno} [{token}]: {line.strip()}")

        self.assertEqual(
            violations,
            [],
            "Active source files must not contain any legacy LeafMachine2 / PointRend tokens:\n"
            + "\n".join(violations),
        )

    def test_no_files_or_folders_named_lm2_or_leafmachine2(self):
        """Verify that zero files or folders in the workspace use 'lm2' or 'leafmachine2' in their naming."""
        skip_roots = {".git", ".venv", ".venv_LM3", "LeafMachine3", "data", "runs", "outputs"}
        bad_paths: List[str] = []
        for entry in sorted(PROJECT_ROOT.iterdir()):
            name_low = entry.name.lower()
            if "lm2" in name_low or "leafmachine2" in name_low:
                bad_paths.append(str(entry.relative_to(PROJECT_ROOT)))
            if entry.name in skip_roots or not entry.is_dir():
                continue
            for sub in entry.rglob("*"):
                sub_low = sub.name.lower()
                if "lm2" in sub_low or "leafmachine2" in sub_low:
                    bad_paths.append(str(sub.relative_to(PROJECT_ROOT)))

        self.assertEqual(
            bad_paths,
            [],
            "No files or folders in the workspace may use 'lm2' or 'leafmachine2' in their name:\n"
            + "\n".join(bad_paths),
        )

    def test_no_active_files_reference_obsolete_sam2_lm2_or_pruned_helpers(self):
        """Scan all active .py, .R, and .sh files (excluding test_code_hygiene.py) and assert zero forbidden SAM-2, LM2, or pruned helper tokens."""
        violations: List[str] = []
        for src_path in iter_active_source_files(include_shell=True):
            if src_path == SELF_PATH:
                continue
            content = src_path.read_text(encoding="utf-8", errors="ignore")
            rel = src_path.relative_to(PROJECT_ROOT)
            for lineno, line in enumerate(content.splitlines(), start=1):
                for token in FORBIDDEN_OBSOLETE_SAM2_AND_PRUNED_HELPER_TOKENS:
                    if token in line:
                        violations.append(f"{rel}:{lineno} [{token}]: {line.strip()}")

        self.assertEqual(
            violations,
            [],
            "Active source files must not contain any obsolete SAM-2, LeafMachine2, or pruned helper tokens:\n"
            + "\n".join(violations),
        )


if __name__ == "__main__":
    unittest.main()
