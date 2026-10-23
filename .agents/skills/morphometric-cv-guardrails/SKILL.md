---
name: morphometric-cv-guardrails
description: Enforces scientific, mathematical, and computer vision protocols for LeafMachine2, Detectron2 PointRend, OpenCV geometry, and R Momocs/MorphoTools2 morphometrics. Use whenever modifying scripts in scripts/vision/, scripts/pipeline/, scripts/annotation_and_training/, or scripts/morphometrics/.
---

# Morphometric & Computer Vision Guardrails (`UNC-BOT-SOP-2026-04-REV4`)

Activate and follow these scientific, mathematical, and computer vision guardrails whenever creating, refactoring, or debugging code in:
* `scripts/vision/` (`lm2_geometry_utils.py`, `configure_leafmachine2.py`, `capitulum_phenology_classifier.py`, `tune_geometry_parameters.py`)
* `scripts/pipeline/` (`02_segment_and_extract.py`)
* `scripts/annotation_and_training/` (`annotate_with_sam2.py`, `sam2_annotator_utils.py`, `finetune_lm2_pcd.py`)
* `scripts/morphometrics/` (`03_fourier_extractor.R`, `04_gmm_morphotools.R`)

---

## 1. Two-Stage LeafMachine2 Handshake

Never run Detectron2 PointRend directly on full-resolution herbarium sheets, as PointRend is trained on cropped single-organ bounding boxes and fails when presented with whole-sheet scale variance and background clutter.

1. **Stage 1 (Full-Sheet Component & Ruler Detection):**
   * Run LeafMachine2's stock YOLOv5 Plant Component Detector (`LeafPriority.pt`) and the **LM2 38-class Ruler Detector** on full herbarium sheets to localize organ bounding boxes (`bbox`) and extract empirical scale bars (`pixels_per_mm`, `is_calibrated_mm`).
   * Keep `detect_ruler_scale_hough()` only as an asynchronous secondary fallback if the 38-class Ruler Detector finds no ruler block.
2. **Stage 2 (Single-Leaf Native-DPI PointRend Segmentation):**
   * For each candidate leaf bounding box from Stage 1, extract a **12%-padded native-DPI single-leaf crop** from the original high-resolution voucher image.
   * Run the fine-tuned Detectron2 PointRend checkpoint (`leaf_seg_rcnn_v2_1` / `models/lm2_packera_pcd_finetuned.pth`) **strictly on these 12%-padded native-DPI single-leaf crops** (never on full sheets).
   * Project the segmented binary crop mask back into full-sheet coordinates only for spatial clustering (`cluster_plant_individuals`), while preserving native-DPI crop resolution for contour extraction.

---

## 2. Label-Blind Gatekeeping

To prevent taxonomic circularity under the Unified Species Concept (de Queiroz 2007), all geometric quality gatekeeping and tier routing in `classify_and_route_leaf()` must be **100% blind to the specimen's species determination (`taxon` / `scientificName` / `species_raw`)**.

* **Pristine Tier 1 (`Tier_1_Direct`):** Requires $\text{UCS} \ge 0.85$ and $\text{Solidity} \ge 0.72$.
* **Dissected / Lyrate-Pinnatifid Tier 1 (`Tier 1 (Dissected)`):**
  * Low-solidity leaves in the range $0.50 \le \text{Solidity} < 0.72$ must be routed into `Tier 1 (Dissected)` **solely** on the basis of intrinsic geometric bilateral sinus topology via `is_botanical_dissection()` (evaluating paired convexity defects, sinus depth ratios $\ge 0.08$, and bilateral lobe symmetry across the midrib).
  * **Forbidden Pattern:** Never use `taxon_has_lyrate_tendency` or any species-name lookup (`Packera paupercula`, `Packera plattensis`, etc.) to bypass or override `is_botanical_dissection()`. A damaged *P. paupercula* leaf with unilateral herbivory must not be promoted to Tier 1 because of its taxonomic label, and a genuinely lyrate-pinnatifid *P. dubia* or *P. anonyma* leaf must not be rejected because of its label.

---

## 3. OpenCV Affine & Contour Math

Every geometric operation in `scripts/vision/lm2_geometry_utils.py` and `scripts/pipeline/02_segment_and_extract.py` must obey these four exact coordinate conventions:

1. **Crop with Padding Before Rotating:**
   * Before calling `cv2.warpAffine()`, always crop the binary mask to its tight foreground bounding box plus sufficient zero-border padding (e.g., `pad = int(math.ceil(0.5 * math.hypot(h, w)))` via `cv2.copyMakeBorder`) so that affine rotation never clips apices, basal lobes, or reflected lamina margins against image canvas boundaries.
2. **Sign of Rotation Angle (`+angle` in OpenCV):**
   * Given the directed midrib vector $\overrightarrow{p_{\text{base}} p_{\text{apex}}} = (\Delta x, \Delta y) = (x_{\text{apex}} - x_{\text{base}},\, y_{\text{apex}} - y_{\text{base}})$ with $\text{angle} = \text{degrees}(\text{atan2}(\Delta y, \Delta x))$ in image coordinates ($y$ pointing downward), OpenCV's `cv2.getRotationMatrix2D(center, theta, 1.0)` rotates **counter-clockwise** in screen space (which equals clockwise in standard Cartesian space) by `theta`.
   * Therefore, always pass **`+angle`** (`cv2.getRotationMatrix2D((cx, cy), +angle, 1.0)`) to rotate the directed vector $\overrightarrow{p_{\text{base}} p_{\text{apex}}}$ onto the positive horizontal $+x$ axis ($0^\circ$). Passing `-angle` doubles the inclination to $2 \cdot \text{angle}$ and breaks horizontal midrib alignment.
3. **Disambiguate `p_base` vs. `p_apex` via Terminal Stalk Width Tapering:**
   * Principal axis projection (`cv2.fitLine`) identifies the two longitudinal extrema along the midrib but is sign-ambiguous ($180^\circ$ arbitrary).
   * Always disambiguate `p_base` (petiole insertion) from `p_apex` (distal blade tip) by measuring orthogonal cross-sectional lamina width across the terminal $5\%\text{–}15\%$ zones of both extrema: the petiole base (`p_base`) exhibits a narrower, high-aspect-ratio stalk taper prior to blade inflection (`trim_petiole_tail`), whereas `p_apex` exhibits a broader, convex distal lamina termination.
4. **Isotropic Contour Normalization (`max(span_x, span_y)`):**
   * When exporting normalized contour coordinates (`x_norm`, `y_norm`) in `export_standardized_contour()` and `resample_normalized_contour()`, always scale **both** axes isotropically by the single scalar `max_span = max(span_x, span_y, 1e-6)`:
     ```python
     max_span = max(span_x, span_y, 1e-6)
     x_norm = (x_resamp - min_x) / max_span
     y_norm = (y_resamp - min_y) / max_span
     ```
   * **Forbidden Pattern:** Never divide `x` by `span_x` and `y` by `span_y` independently, which stretches every leaf into a $1 \times 1$ unit square and destroys the diagnostic length-to-width aspect ratio (e.g., collapsing narrowly oblanceolate *P. anonyma* with $L/W > 4:1$ into the same aspect ratio as suborbicular *P. paupercula*).

---

## 4. R Elliptic Fourier Analysis (EFA) & Allometry

When modifying `scripts/morphometrics/03_fourier_extractor.R` and `scripts/morphometrics/04_gmm_morphotools.R`, enforce the following mathematical protocols:

1. **Oriented Aspect Ratio via `Momocs::coo_lw`:**
   * Compute leaf length-to-width aspect ratio from the principal-axis-aligned contour via `Momocs::coo_lw(coo)` (or oriented minimum bounding box / SVD principal axes fallback), never from unaligned canvas bounding boxes (`diff(range(x)) / diff(range(y))`), which vary arbitrarily with mounting angle on the herbarium sheet.
2. **Pre-Aligned EFA without First-Harmonic $180^\circ$ Flips (`norm = FALSE, start = FALSE`):**
   * Because contours exported by Phase 2 are already homologized to start at `p_base`, aligned along the directed $\overrightarrow{p_{\text{base}} p_{\text{apex}}}$ axis, centered, and scaled isotropically by centroid size / `max_span`, run:
     ```r
     Momocs::efourier(out_obj, nb.h = opts$harmonics, norm = FALSE, start = FALSE)
     ```
     (after explicit baseline alignment `coo_center()`, `coo_scale()`, and base-to-apex horizontal alignment).
   * **Why:** Setting `norm = TRUE` in `Momocs::efourier()` aligns outlines using the major axis of the first harmonic ellipse ($E_1$), which is invariant to $180^\circ$ reversals on nearly elliptical blades and randomly flips obovate/spatulate leaves (*P. anonyma*, *P. plattensis*) into ovate orientations, corrupting odd-numbered symmetric harmonics.
3. **Millimeter-Calibrated Centroid Size ($\text{CS}_{\text{mm}}$):**
   * Always propagate `pixels_per_mm` and `is_calibrated_mm` from the extracted leaf manifest into `03_fourier_extractor.R` and `04_gmm_morphotools.R`.
   * For calibrated specimens (`is_calibrated_mm == TRUE` and `pixels_per_mm > 0`), compute true physical Centroid Size in millimeters:
     $$\text{CS}_{\text{mm}} = \frac{\text{CS}_{\text{px}}}{\text{pixels\_per\_mm}}, \qquad \text{log\_centsize} = \ln(\text{CS}_{\text{mm}})$$
   * Never mix raw uncalibrated pixel centroid sizes ($\text{CS}_{\text{px}}$, which conflates camera sensor DPI with biological leaf size) with millimeter-calibrated centroid sizes in allometric regressions.
4. **Pooled Within-Species Allometric Regression on Calibrated Tier 1 Gold Anchors:**
   * In `test_and_correct_allometry()` (`04_gmm_morphotools.R`), fit the allometric model on **calibrated (`is_calibrated_mm == TRUE`) Tier 1 Gold (`determiner_tier == "Tier_1_Gold"`)** specimens using the additive pooled within-species model:
     ```r
     fit_allom <- stats::lm(Y_gold ~ log_centsize + scientificName, data = gold_calibrated_df)
     ```
   * Extract the pooled within-species allometric slope vector $\hat{\beta}_{\text{size}}$ (the coefficient for `log_centsize`) and subtract the size-dependent component $\hat{\beta}_{\text{size}} \cdot (\text{log\_centsize}_i - \overline{\text{log\_centsize}})$ from all calibrated specimens, preserving the between-species intercept differences ($\alpha_{\text{species}}$).
   * **Forbidden Pattern:** Never take raw residuals from a global `lm(Y ~ log_centsize)` across all pooled species without controlling for `scientificName`; when species naturally differ in mean leaf size, a naive global regression strips out true interspecific taxonomic divergence along with ontogenetic allometry.

---

## Verification Checklist Before Completing Changes

1. Confirm no changes were made inside `LeafMachine2/`.
2. Run unit tests covering geometry, orientation, and Fourier decomposition:
   ```bash
   python scripts/tests/run_audit.py
   ```
3. Verify that `test_geometry_edge_cases.py`, `test_symmetric_fourier_validation.py`, and `test_allometry_and_phenology.py` pass with 0 failures.
