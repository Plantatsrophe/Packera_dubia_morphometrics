---
name: pipeline-contract-audit
description: Enforces LeafMachine3 computer vision and morphometric guardrails, audits cross-stage Python and R SQLite/CSV schema contracts, and runs the test suite.
---

# Pipeline Contract Audit & CV Guardrails (`UNC-BOT-SOP-2026-04-REV5`)

Use this skill before committing code changes, after modifying pipeline stages (`01` through `08`), or when diagnosing missing-column, coordinate unscaling, or join-key discrepancies across Python and R scripts.

---

## Part 1: Computer Vision & Morphometric Guardrails (LeafMachine3)

Every CV, extraction, and morphometrics module in `scripts/vision/`, `scripts/pipeline/`, and `scripts/morphometrics/` must enforce these invariants:

### 1. Working-Frame Coordinate Unscaling (`work_scale`)
* When LeafMachine3 processes downscaled/letterboxed representations (`work_scale`), all bounding boxes, polygon vertices, landmark coordinates, and binary masks saved to database tables or exported manifests must be divided by `work_scale` to return to native full-sheet coordinates:
  $$\text{coord}_{\text{native}} = \frac{\text{coord}_{\text{work}}}{\text{work\_scale}}$$
* Columns with suffix `_px` (e.g., `petiole_length_px`, `plant_height_px`, contour vertices) must strictly represent coordinates on the original unscaled image.

### 2. 3-Class Organ Mask Handling & UNet++ v9 Rachis Bridging
* Module 8 instance segmentation (`lm3_leaf_segmenter__yolo26x_seg_1024`) outputs 3 semantic classes: `Leaf` (lamina, class 0), `Petiole` (class 1), and `Hole` (sinuses/perforations, class 2).
* Hole pixels falling within lamina polygons must be subtracted from the blade mask.
* In deeply dissected or pinnatifid leaves where narrow rachis segments experience minor breaks, Module 4 UNet++ v9 `specimen_paper_cleaned_mask` (`LLR_Hyst_KeepBlobs: true` color likelihood) acts as an anatomical bridge constraint (`bridge_dissected_rachis_components`), reconnecting genuine pinnae without bridging across mounting paper.

### 3. Tip-Up Clockwise `lamina_base`-Anchored Contour Geometry
* Normalized 120-vertex contours in `extracted_leaves.contour_xy_json` must be oriented vertically with apex pointing upward ($y$-axis toward apex) and starting homologously at the anatomical lamina-petiole junction (`lamina_base` anchor $p_{\text{base}}$), traversing clockwise around the perimeter.
* Because coordinates are standardized and homologized in Step 02, `Momocs::efourier` in Step 03 must be executed with `norm = FALSE` and `start = FALSE`:
  ```r
  Momocs::efourier(out_obj, nb.h = H_opt, norm = FALSE, start = FALSE)
  ```
  preserving biologically meaningful absolute orientation, rotation, and starting-point homology.

### 4. Strict Active vs. Passive Determiner Authority Split
* **Active Reference Anchors (`Tier_1_Gold`, `Tier_2_Silver`, promoted `Tier_3_Bronze`, `Tier_0_Type`):** Used exclusively for model fitting, Mahalanobis scatter estimation, CDA discriminant axes, and Spatial RF training.
* **Passively Projected Vouchers (`Tier_4_Unverified` / unpromoted `Tier_3_Bronze`):** Assigned weight $0.0$ during training and projected into fitted spaces without altering reference geometry.

---

## Part 2: Automated Cross-Stage SQLite/CSV Contract & Hygiene Audit

### 1. Zero `_archive` Coupling & Legacy Token Scan
Active files must never reference `_archive` paths and must contain zero legacy framework remnants:

```bash
python3 - <<'EOF'
import sys
from pathlib import Path

root = Path(".").resolve()
active_files = [
    p for p in list(root.rglob("*.py")) + list(root.rglob("*.R"))
    if "_archive" not in p.parts and ".venv" not in p.parts
    and ".venv_LM3" not in p.parts and "LeafMachine3" not in p.parts
]

violations = []
for f in sorted(active_files):
    rel = f.relative_to(root)
    if rel.parts[:2] == ("scripts", "tests"):
        continue
    for lineno, line in enumerate(f.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if "_archive" in line:
            violations.append(f"{rel}:{lineno}: {stripped}")

if violations:
    print("[FAIL] Forbidden _archive references found:")
    for v in violations:
        print("  -", v)
    sys.exit(1)
print(f"[PASS] Scanned {len(active_files)} active files: 0 _archive leaks.")
EOF
```

### 2. Cross-Stage Handshake Contract Verification Script
Verifies column agreements, composite keys `(catalogNumber, plant_individual_id)`, and handshakes across stages:

```bash
python3 - <<'EOF'
import sys
from pathlib import Path

scripts = [
    "scripts/core/db.py", "scripts/core/db_utils.R", "scripts/core/harvester.py",
    "scripts/pipeline/02_segment_and_extract.py", "scripts/morphometrics/03_fourier_extractor.R",
    "scripts/morphometrics/04_gmm_morphotools.R", "scripts/morphometrics/authority_weighting_utils.R",
    "scripts/analysis/05_cleanlab_vision_xai.py", "scripts/analysis/vision_dissection_utils.py",
    "scripts/analysis/06_multimodal_spatial_rf.R", "scripts/analysis/07_triage_dashboard_synthesis.R",
    "scripts/analysis/08_discover_taxon_complexes.R",
]

missing = [s for s in scripts if not Path(s).exists()]
if missing:
    print(f"[FAIL] Missing pipeline scripts: {missing}")
    sys.exit(1)

s_db = Path("scripts/core/db.py").read_text()
s_db_r = Path("scripts/core/db_utils.R").read_text()
s_harv = Path("scripts/core/harvester.py").read_text()
s02 = Path("scripts/pipeline/02_segment_and_extract.py").read_text()
s03 = Path("scripts/morphometrics/03_fourier_extractor.R").read_text()
s04 = Path("scripts/morphometrics/04_gmm_morphotools.R").read_text()
s_auth = Path("scripts/morphometrics/authority_weighting_utils.R").read_text()
s05 = Path("scripts/analysis/05_cleanlab_vision_xai.py").read_text()
s_vis = Path("scripts/analysis/vision_dissection_utils.py").read_text()
s06 = Path("scripts/analysis/06_multimodal_spatial_rf.R").read_text()
s07 = Path("scripts/analysis/07_triage_dashboard_synthesis.R").read_text()
s08 = Path("scripts/analysis/08_discover_taxon_complexes.R").read_text()

checks = [
    ("db.py defines schema, WAL mode, tables, and master view", "v_integrative_triage_master" in s_db and "PRAGMA journal_mode = WAL" in s_db),
    ("db_utils.R exposes R connection and master view helpers", "v_integrative_triage_master" in s_db_r and "RSQLite::SQLite()" in s_db_r),
    (".discover_step08_eco_cols extracts type ecological baseline", ".discover_step08_eco_cols" in s_db_r and "type_edaphic_z_score" in s_db_r),
    ("01 parses determiner tier with collector isolation & duplicate guard", "_should_fallback_to_collector" in s_harv and "exsiccatae_taxonomic_conflict" in s_harv),
    ("01 defines nomenclatural type hierarchy, subranks & cascading priority", "PRIMARY_TYPE_STATUSES" in s_harv and "classify_nomenclatural_type_rank" in s_harv and "_type_rank" in s_harv),
    ("01 historical types exempt from uncertainty filters and retained in thinning", "is_rec_type" in s_harv and "spatially_balance_and_thin_vouchers" in s_harv),
    ("02 exports assigned_tier, pixels_per_mm & rosette centroids", "assigned_tier" in s02 and "pixels_per_mm" in s02 and "rosette_centroids_by_pid" in s02),
    ("02 supports Class 0/1 candidates & is_cauline_or_stem_contaminant", "s1_cls not in (0, 1)" in s02 and "is_cauline_or_stem_contaminant" in s02),
    ("03 groups by (catalogNumber, plant_individual_id), uses cov.wt & exports morph_tier_source", "plant_individual_id" in s03 and "catalogNumber" in s03 and "cov.wt" in s03),
    ("03 matches tier via case-insensitive regex", "tier_?1" in s03.lower() or "tier_?[12]" in s03.lower()),
    ("03 exports harmonic_roughness, is_contour_outlier & leaf_quality_weight", "harmonic_roughness" in s03 and "is_contour_outlier" in s03),
    ("03 enforces nomenclatural type MCD outlier immunity & precision weighting", "mcd_out" in s03 and "is_type_specimen" in s03 and "compute_authority_precision_weights" in s03),
    ("04 enforces 2-step Dawid-Skene / EB authority update & full-rank can3", "compute_empirical_bayes_authority_weights" in s_auth and "can3" in s04),
    ("04 exports morphometric flags, loocv_accuracy & cda_posterior_margin", "morphometric_outlier_flag" in s04 and "loocv_accuracy" in s04 and "cda_posterior_margin" in s04),
    ("04 enforces type-medoid shrinkage & type morphological control", "compute_type_specimen_morph_control" in s04 and ".identify_type_specimen_vec" in s_auth),
    ("05 uses promoted bronze anchor parity & abaxial artifact gate", "is_promoted_bronze" in s05 and "filter_abaxial_artifacts" in s05),
    ("05 uses foreground patch pooling, OOD scoring & attention ratio", ("foreground_masked_patch_pooling" in s05 or "compute_foreground_tissue_mask" in s05) and "vision_ood_score" in s05),
    ("05 enforces type visual prototype control & Cleanlab corruption immunity", "compute_type_visual_control" in s_vis and "fit_calibrated_oof_fold_classifier" in s05),
    ("06 uses 0-15cm depth weighting, VIF pruning & spatial_block_cv_acc", "left_join" in s06 and "prune_collinear_predictors" in s06 and "spatial_block_cv_acc" in s06),
    ("06 pools Spatial Block-CV RF posteriors and computes edaphic_z_score", "spatial_block_cv_acc" in s06 and "edaphic_z_score" in s06),
    ("06 derives locus classicus reference set & type_edaphic_z_score", "compute_type_locality_ecological_control" in s06 and "type_edaphic_z_score" in s06),
    ("06 enforces resolve_step06_authority_weights fallback pinning", "resolve_step06_authority_weights" in s06 and "3.0" in s06),
    ("07 validates contracts & coalesces non-georeferenced vouchers", "validate_upstream_contracts" in s07 and "taxonomic_treatment_broad" in s07 and "is_completed_voucher" in s07),
    ("07 enforces cross-modal non-contradiction and precision weighting", "morph_or_vis_corrupt_reassign" in s07 and ("prec_morph" in s07 or "lookup_prec" in s07)),
    ("07 enforces type control Rules A, B, and C", "Confirmed_Type_Control" in s07 and "Type_Synonymy_Candidate" in s07 and "is_sensu_stricto_locked_i" in s07),
    ("08 enforces type-anchored complex discovery and delimitation", "select_taxon_anchors" in s08 and "is_type_spec" in s08 and "gamma_k" in s08),
]

failed = False
for label, ok in checks:
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {label}")
    if not ok:
        failed = True
if failed:
    sys.exit(1)
EOF
```

### 3. Unit & Integration Test Suite Verification Gate
Always execute the comprehensive test suite to confirm zero regressions:

```bash
pytest scripts/tests/ -q
pytest scripts/tests/test_code_hygiene.py -q
```
