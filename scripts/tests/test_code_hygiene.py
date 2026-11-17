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

import re
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
    "sam2", "SAM2", "annotate_with_sam2", "sam2_annotator_utils", "SAM2BoundaryRefiner",
    "SAM2PromptBundle", "refine_mask_with_sam2", "tune_geometry_parameters", "run_calibrate_geometry",
    # Pruned helper tokens
    "detect_ruler_scale_hough", "trim_petiole_tail", "_sample_cross_sectional_widths",
    "_disambiguate_base_apex_endpoints", "homologize_contour_starting_point",
    "_align_and_normalize_contour", "align_contour_basal_apical", "is_non_lamina_petiole_ribbon",
    "extract_protologue_authors",
    # Pre-LeafMachine3 legacy Python QC, mask-bridging & Tier-2 reflection tokens
    "leaf_gatekeeper", "LeafRoutingResult", "TIER_1_DIRECT", "TIER_2_REFLECTED",
    "Tier_1_Direct", "Tier_2_Reflected", "Tier_1_Pure", "Tier_2_Fallback",
    "tier1_intact", "tier2_reflected", "assigned_tier", "tier_assignment",
    "ucs_score", "min_ucs", "min_solidity", "enable_geometric_qc",
    "is_botanical_dissection", "compute_geometric_metrics", "compute_classical_shape_descriptors",
    "_prepare_lamina_mask", "extract_tier1_pristine", "compute_sheet_adaptive_min_blade_area_cm2",
    "FoldDetectionConfig", "DissectionConfig", "n_tier1_leaves", "n_tier2_leaves",
    "is_folded", "reflection_applied", "render_tier1_vs_tier2_2d_pdf",
    "tier_symmetry_validation", "tier1_vs_tier2_density_overlay", "find_lamina_midrib_endpoints",
    # Purged Steps 1-4 twin columns & circular helpers (Zero-Trust Variable Overhaul)
    "determiner_raw", "pheno_sin", "pheno_cos", "calculate_circular_phenology",
    "capitulum_aspect_ratio", "blade_length_to_width_ratio", "blade_area_cm2",
    "x_norm", "y_norm", "Sym_PC", "is_promoted_anchor", "heteroblasty_index",
    "whiten_intraspecific_deme_plasticity", "compute_cluster_spatial_segregation_p",
    "residualize_within_deme_phenophase_plasticity", ".augment_foliar_architecture_features",
    ".compute_soft_geo_range_log_prior",
    # Purged Step 6 twin columns, obsolete vision hooks, and misnumbered step07/08 helpers
    "w_vision_kappa", "dinov2_predicted_taxon", "vision_conflict_flag", "final_curation_triage",
    "morphometrics_misidentification_flags", "probabilistic_anomaly_score", "composite_conflict_score",
    "integrative_consensus_prob", "n_active_modalities", "render_step07_synthesis_plate_pdf",
    "render_step08_complex_discovery_pdf", "load_step08_merged_inputs",
    "fig07_multimodal_discordance_summary", "fig07b_taxonomic_treatment_summary",
    # Retired Step 06 synthesis / Step 06b / Step 08 helper names and complex-table twin columns
    "06_triage_dashboard_synthesis", "step06b", "Step 06b", "fig06a_integrative_taxonomic_synthesis",
    "fig06b_taxon_complex_delimitation_network", "fig06c_", "ingest_genetic_data",
    "07_discover_taxon_complexes", "render_step07_complex_discovery_pdf",
    ".discover_step08_eco_cols", ".discover_step08_dynamic_morph_cols",
    "resolve_step08_morph_features", "resolve_step08_eco_predictors",
    ".parse_step08_semicolon", ".sort_step08_pcs", "n_Pure_Gold_Silver",
    "leiden_community", "ward_silhouette_cluster", "top_neighbor_1", "top_neighbor_2",
    "kappa_morph_weight", "kappa_eco_weight", "fig_empirical_taxon_complex_network",
    "vision_xai_audit", "dinov2_feature_cache", "cda_biplot_pdf",
    "generate_cda_biplot_pdf", "render_allometry_diagnostic_pdf", "render_cda_biplot_4panel_pdf",
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
                    if token in ("x_norm", "y_norm", "ingest_genetic_data"):
                        hit = bool(re.search(rf"(?<![A-Za-z0-9_]){re.escape(token)}(?![A-Za-z0-9_])", line))
                    elif token == "fig06c_":
                        hit = bool(re.search(r"fig06c_(?!cytonuclear_paralog_cytotype)", line))
                    else:
                        hit = token in line
                    if hit:
                        violations.append(f"{rel}:{lineno} [{token}]: {line.strip()}")

        self.assertEqual(
            violations,
            [],
            "Active source files must not contain any obsolete SAM-2, LeafMachine2, or pruned helper tokens:\n"
            + "\n".join(violations),
        )

    def test_new_public_genomics_modules_respect_500_line_ceiling(self):
        """Verify public Packera genomics, cytology, taxonomy, geochemistry, monograph, and hygiene modules stay within line ceilings."""
        strict_ceilings = {
            PROJECT_ROOT / "scripts" / "data_prep" / "download_public_packera_sequences.py": 500,
            PROJECT_ROOT / "scripts" / "data_prep" / "ena_paraloss_discovery.py": 300,
            PROJECT_ROOT / "scripts" / "data_prep" / "download_public_cytology.py": 450,
            PROJECT_ROOT / "scripts" / "analysis" / "step07_cytology_integration.R": 380,
            PROJECT_ROOT / "scripts" / "tests" / "test_public_cytology_integration.py": 380,
            PROJECT_ROOT / "scripts" / "data_prep" / "download_taxonomic_backbone.py": 450,
            PROJECT_ROOT / "scripts" / "data_prep" / "download_geochem_wetland_layers.py": 450,
            PROJECT_ROOT / "scripts" / "analysis" / "edaphic_geochem_wetland_utils.R": 320,
            PROJECT_ROOT / "scripts" / "tests" / "test_external_data_sources.py": 420,
            PROJECT_ROOT / "scripts" / "figures" / "step06c_figures.R": 350,
            PROJECT_ROOT / "scripts" / "analysis" / "step06_coalescent_models.py": 380,
            PROJECT_ROOT / "scripts" / "analysis" / "step06_model_selection.R": 380,
            PROJECT_ROOT / "scripts" / "analysis" / "step06_spatial_landscape_genomics.R": 380,
            PROJECT_ROOT / "scripts" / "tests" / "test_step06_model_selection.py": 420,
            PROJECT_ROOT / "scripts" / "morphometrics" / "gmm_icl_bootstrap_utils.R": 280,
            PROJECT_ROOT / "scripts" / "analysis" / "step05_circular_hopkins_models.R": 320,
            PROJECT_ROOT / "scripts" / "tests" / "test_step04_05_model_selection.py": 400,
            PROJECT_ROOT / "scripts" / "core" / "coord_hygiene_utils.py": 260,
            PROJECT_ROOT / "scripts" / "analysis" / "step05_env_outlier_utils.R": 260,
            PROJECT_ROOT / "scripts" / "analysis" / "step06_rogue_sequence_qc.R": 280,
            PROJECT_ROOT / "scripts" / "tests" / "test_outlier_quarantine_gates.py": 400,
            PROJECT_ROOT / "scripts" / "analysis" / "step06_public_backbone.R": 475,
            PROJECT_ROOT / "scripts" / "analysis" / "step07_deme_propagation_utils.R": 260,
            PROJECT_ROOT / "scripts" / "analysis" / "step07_genomic_propagation.R": 400,
            PROJECT_ROOT / "scripts" / "tests" / "test_step07_deme_propagation.py": 300,
            PROJECT_ROOT / "scripts" / "data_prep" / "organellar_cytonuclear_utils.py": 220,
            PROJECT_ROOT / "scripts" / "data_prep" / "skim_offtarget_organellar.py": 460,
            PROJECT_ROOT / "scripts" / "analysis" / "step07_organellar_triage_utils.R": 260,
            PROJECT_ROOT / "scripts" / "tests" / "test_step07_organellar_triage.py": 280,
            PROJECT_ROOT / "scripts" / "analysis" / "step05_layer_hygiene_utils.R": 280,
            PROJECT_ROOT / "scripts" / "tests" / "test_step05_layer_hygiene.py": 380,
            PROJECT_ROOT / "scripts" / "core" / "step01_biogeography_utils.py": 300,
            PROJECT_ROOT / "scripts" / "data_prep" / "prune_eastern_cohort.py": 500,
            PROJECT_ROOT / "scripts" / "tests" / "test_step01_robustness_audit.py": 380,
            PROJECT_ROOT / "scripts" / "vision" / "step02_photometry_corymb_utils.py": 360,
            PROJECT_ROOT / "scripts" / "tests" / "test_step02_vision_hygiene.py": 380,
            PROJECT_ROOT / "scripts" / "morphometrics" / "step03_config_scale_utils.R": 240,
            PROJECT_ROOT / "scripts" / "morphometrics" / "03_fourier_extractor.R": 500,
            PROJECT_ROOT / "scripts" / "morphometrics" / "fourier_parallel_utils.R": 500,
            PROJECT_ROOT / "scripts" / "morphometrics" / "plasticity_mitigation_utils.R": 470,
            PROJECT_ROOT / "scripts" / "figures" / "step03_figures.R": 385,
            PROJECT_ROOT / "scripts" / "tests" / "test_step03_audit_contracts.py": 360,
            PROJECT_ROOT / "scripts" / "morphometrics" / "04_gmm_morphotools.R": 500,
            PROJECT_ROOT / "scripts" / "morphometrics" / "step04_allometry_merger_utils.R": 380,
            PROJECT_ROOT / "scripts" / "morphometrics" / "pub_synthesis_plates.R": 120,
            PROJECT_ROOT / "scripts" / "figures" / "step04_figures.R": 500,
            PROJECT_ROOT / "scripts" / "tests" / "test_step04_audit_hygiene.py": 380,
            PROJECT_ROOT / "scripts" / "analysis" / "step05_config_geo_utils.R": 320,
            PROJECT_ROOT / "scripts" / "figures" / "step05_figures.R": 500,
            PROJECT_ROOT / "scripts" / "figures" / "step05_sdm_figures.R": 320,
            PROJECT_ROOT / "scripts" / "tests" / "test_step05_audit_hygiene.py": 380,
        }
        bounded_paths = (
            *strict_ceilings.keys(),
            PROJECT_ROOT / "scripts" / "analysis" / "genomic_triage_utils.R",
            PROJECT_ROOT / "scripts" / "pipeline" / "resolve_eastern_na_packera.py",
            PROJECT_ROOT / "scripts" / "tests" / "test_eastern_na_monograph_cli.py",
        )
        for mod_path in bounded_paths:
            self.assertTrue(mod_path.exists(), f"Missing expected bounded module: {mod_path}")
            line_count = len(mod_path.read_text(encoding="utf-8").splitlines())
            self.assertLess(line_count, 500, f"{mod_path.relative_to(PROJECT_ROOT)} exceeds 500-line ceiling ({line_count} lines)")
            if mod_path in strict_ceilings:
                ceiling = strict_ceilings[mod_path]
                self.assertLess(line_count, ceiling, f"{mod_path.relative_to(PROJECT_ROOT)} exceeds {ceiling}-line ceiling ({line_count} lines)")

    def test_data_cytology_directory_stays_flat(self):
        """Verify `data/cytology/` exists and contains zero nested subdirectories."""
        cytology_dir = PROJECT_ROOT / "data" / "cytology"
        self.assertTrue(cytology_dir.is_dir(), f"Missing expected cytology directory: {cytology_dir}")
        subdirs = [p.name for p in cytology_dir.iterdir() if p.is_dir()]
        self.assertEqual(subdirs, [], f"data/cytology/ must remain flat with zero subdirectories: {subdirs}")

    def test_data_taxonomy_directory_stays_flat(self):
        """Verify `data/taxonomy/` contains zero nested subdirectories when present."""
        taxonomy_dir = PROJECT_ROOT / "data" / "taxonomy"
        if taxonomy_dir.exists():
            self.assertTrue(taxonomy_dir.is_dir(), f"Expected directory at {taxonomy_dir}")
            subdirs = [p.name for p in taxonomy_dir.iterdir() if p.is_dir()]
            self.assertEqual(subdirs, [], f"data/taxonomy/ must remain flat with zero subdirectories: {subdirs}")

    def test_no_hardcoded_workstation_legacy_repo_path(self):
        """Scan all active .py, .R, and .sh files and assert zero hardcoded workstation legacy repo paths."""
        legacy_abs_path = "/" + "home/brandon/" + "Packera_dubia_morphometrics"
        violations: List[str] = []
        for src_path in iter_active_source_files(include_shell=True):
            if src_path == SELF_PATH:
                continue
            content = src_path.read_text(encoding="utf-8", errors="ignore")
            rel = src_path.relative_to(PROJECT_ROOT)
            for lineno, line in enumerate(content.splitlines(), start=1):
                if legacy_abs_path in line:
                    violations.append(f"{rel}:{lineno}: {line.strip()}")

        self.assertEqual(
            violations,
            [],
            "Active source files must not hardcode the legacy workstation repository path:\n"
            + "\n".join(violations),
        )

    def test_no_duplicate_python_defs_in_same_scope(self):
        """Assert zero shadowed or duplicate def/async def/class definitions within any module or class body."""
        def _is_overload(node: ast.AST) -> bool:
            for dec in getattr(node, "decorator_list", []):
                if isinstance(dec, ast.Name) and dec.id == "overload":
                    return True
                if isinstance(dec, ast.Attribute) and dec.attr == "overload":
                    return True
            return False

        def _check_body(body: List[ast.stmt], scope_label: str, dupes: List[str]) -> None:
            seen: dict[str, int] = {}
            for stmt in body:
                if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    if _is_overload(stmt):
                        continue
                    if stmt.name in seen:
                        dupes.append(f"{scope_label}:{stmt.lineno} duplicates {stmt.name} (first at L{seen[stmt.name]})")
                    else:
                        seen[stmt.name] = stmt.lineno
                    if isinstance(stmt, ast.ClassDef):
                        _check_body(stmt.body, f"{scope_label}.{stmt.name}", dupes)

        dupes: List[str] = []
        py_files = [PROJECT_ROOT / "main.py"] + sorted(
            p for p in (PROJECT_ROOT / "scripts").rglob("*.py")
            if ARCHIVE_MARKER not in p.parts
        )
        for py_path in py_files:
            if not py_path.exists():
                continue
            tree = ast.parse(py_path.read_text(encoding="utf-8"), filename=str(py_path))
            _check_body(tree.body, str(py_path.relative_to(PROJECT_ROOT)), dupes)

        self.assertEqual(dupes, [], "Duplicate Python definitions found in same scope:\n" + "\n".join(dupes))

    def test_no_broken_symlinks(self):
        """Assert zero broken symlinks across models/, scripts/, and runs/latest*."""
        broken: List[str] = []
        for root_dir in (PROJECT_ROOT / "models", PROJECT_ROOT / "scripts"):
            if root_dir.exists():
                for entry in root_dir.rglob("*"):
                    if "__pycache__" in entry.parts:
                        continue
                    if entry.is_symlink() and not entry.exists():
                        broken.append(str(entry.relative_to(PROJECT_ROOT)))
        runs_dir = PROJECT_ROOT / "runs"
        if runs_dir.exists():
            for entry in runs_dir.glob("latest*"):
                if entry.is_symlink() and not entry.exists():
                    broken.append(str(entry.relative_to(PROJECT_ROOT)))
        self.assertEqual(broken, [], "Broken symlinks detected:\n" + "\n".join(broken))

    def test_no_pdf_figure_artifacts(self):
        """Assert zero Rplots.pdf or figure .pdf artifacts at repo root or figure directories."""
        found_pdfs: List[str] = []
        rplots = PROJECT_ROOT / "Rplots.pdf"
        if rplots.exists() or rplots.is_symlink():
            found_pdfs.append("Rplots.pdf")
        out_figs = PROJECT_ROOT / "outputs" / "figures"
        if out_figs.exists():
            for pdf_p in out_figs.rglob("*.pdf"):
                found_pdfs.append(str(pdf_p.relative_to(PROJECT_ROOT)))
        runs_dir = PROJECT_ROOT / "runs"
        if runs_dir.exists():
            for fig_dir in runs_dir.glob("*/outputs/figures"):
                if fig_dir.is_symlink():
                    continue
                for pdf_p in fig_dir.rglob("*.pdf"):
                    found_pdfs.append(str(pdf_p.relative_to(PROJECT_ROOT)))
        self.assertEqual(found_pdfs, [], "Unexpected PDF figure artifacts found:\n" + "\n".join(sorted(found_pdfs)))


if __name__ == "__main__":
    unittest.main()
