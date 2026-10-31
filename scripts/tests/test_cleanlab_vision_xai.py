#!/usr/bin/env python3
"""
scripts/tests/test_cleanlab_vision_xai.py
=========================================
Unit and regression tests for Indumentum-Blind 224x224 Macro-Architecture Pooling,
Foreground-Masked DINOv2 Token Averaging, Authority-Anchored Institution-Blocked
Out-of-Fold Calibration, Confident Learning + Feature-Space OOD, and Pointing Game
Foreground Attribution Auditing.
"""

import importlib
import inspect
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch as mock_patch
import numpy as np
import pandas as pd
from PIL import Image

mod_xai = importlib.import_module("scripts.analysis.05_cleanlab_vision_xai")
standardize_packera_taxon = mod_xai.standardize_packera_taxon
neutralize_mounting_paper = mod_xai.neutralize_mounting_paper
classify_foliar_surface_orientation = mod_xai.classify_foliar_surface_orientation
preprocess_macro_rosette_crop = mod_xai.preprocess_macro_rosette_crop
compute_foreground_tissue_mask = mod_xai.compute_foreground_tissue_mask
foreground_masked_patch_pooling = mod_xai.foreground_masked_patch_pooling
extract_dinov2_embeddings = mod_xai.extract_dinov2_embeddings
compute_out_of_fold_probabilities = mod_xai.compute_out_of_fold_probabilities
compute_feature_space_ood_scores = mod_xai.compute_feature_space_ood_scores
run_confident_learning_audit = mod_xai.run_confident_learning_audit
compute_foreground_attention_ratio = mod_xai.compute_foreground_attention_ratio

TEST_SPECIES = [
    "Packera aurea",
    "Packera obovata",
    "Packera schweinitziana",
    "Packera cymbalaria",
]


class TestCleanlabVisionXAI(unittest.TestCase):
    """Test suite for DINOv2 macro-architecture pooling & Cleanlab XAI."""

    @classmethod
    def setUpClass(cls):
        cls.repo_root = Path(__file__).resolve().parents[2]
        cls.audit_csv = cls.repo_root / "data" / "tables" / "label_noise_audit.csv"
        cls.figure_png = cls.repo_root / "outputs" / "figures" / "GradCAM_audit_panel.png"

    def test_zero_hardcoded_species_lists_in_phase2_and_phase5(self):
        """Verify Phase 5 and Phase 2 helpers contain zero static species lists or hardcoded epithets."""
        self.assertFalse(hasattr(mod_xai, "TARGET_TAXA"))
        self.assertFalse(hasattr(mod_xai, "EASTERN_NA_PACKERA_TAXA"))
        self.assertFalse(hasattr(mod_xai, "OKABE_ITO"))
        self.assertNotIn("TARGET_TAXA", mod_xai.__all__)
        self.assertNotIn("EASTERN_NA_PACKERA_TAXA", mod_xai.__all__)

        xai_path = self.repo_root / "scripts" / "analysis" / "05_cleanlab_vision_xai.py"
        xai_src = xai_path.read_text(encoding="utf-8")
        self.assertLessEqual(len(xai_src.splitlines()), 800, "05_cleanlab_vision_xai.py must remain <= 800 lines.")

        forbidden_epithets = (
            "anonyma", "antennariifolia", "aurea", "crawfordii", "cymbalaria",
            "dubia", "glabella", "hyperborealis", "indecora", "insulanus",
            "memmingeri", "millefolium", "obovata", "paupercula", "plattensis",
            "pseudaurea", "schweinitziana", "serpenticola", "spellenbergii", "tampicana",
        )
        for ep in forbidden_epithets:
            self.assertNotIn(ep, xai_src.lower(), f"Hardcoded epithet '{ep}' must not appear in 05_cleanlab_vision_xai.py")

        step02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")
        std02_src = inspect.getsource(step02._standardize_packera_species_py).lower()
        fig02_src = inspect.getsource(step02.generate_step02_publication_figures).lower()
        for ep in forbidden_epithets:
            self.assertNotIn(ep, std02_src, f"Hardcoded epithet '{ep}' must not appear in _standardize_packera_species_py")
            self.assertNotIn(ep, fig02_src, f"Hardcoded epithet '{ep}' must not appear in generate_step02_publication_figures")

    def test_taxonomic_standardization(self):
        """Test dynamic species synonymy standardizer and taxon pattern validator across arbitrary Packera species."""
        gbif_mock_map = {
            "Senecio obovatus": {"matchType": "EXACT", "species": "Packera obovata", "usageKey": 101},
            "Senecio robbinsii": {"matchType": "EXACT", "species": "Packera schweinitziana", "usageKey": 102},
            "Senecio resedifolius": {"matchType": "EXACT", "species": "Packera cymbalaria", "usageKey": 103},
        }
        with mock_patch(
            "pygbif.species.name_backbone",
            side_effect=lambda name="", **kw: gbif_mock_map.get(str(name).strip(), {"matchType": "NONE"}),
        ), mock_patch("pygbif.species.name_usage", return_value={"results": []}):
            self.assertEqual(standardize_packera_taxon("Packera aurea (L.) Á.Löve & D.Löve"), "Packera aurea")
            self.assertEqual(standardize_packera_taxon("Senecio obovatus Muhl. ex Willd."), "Packera obovata")
            self.assertEqual(standardize_packera_taxon("Senecio robbinsii Oakes"), "Packera schweinitziana")
            self.assertEqual(standardize_packera_taxon("Senecio resedifolius Less."), "Packera cymbalaria")
            self.assertEqual(standardize_packera_taxon("Packera cana (Hook.) W.A.Weber & Á.Löve"), "Packera cana")
            self.assertEqual(standardize_packera_taxon("Packera streptanthifolia (Greene) W.A.Weber & Á.Löve"), "Packera streptanthifolia")
            self.assertEqual(standardize_packera_taxon(None), "Unknown")

        for valid_tx in (
            "Packera aurea",
            "Packera cana",
            "Packera × memmingeri",
            "Packera x hybrida",
            "Morphotype_K1",
            "Ecotype_1_fine_soil_ph_Low",
        ):
            self.assertTrue(mod_xai._is_valid_packera_taxon(valid_tx), f"Expected valid taxon pattern: {valid_tx}")
        for invalid_tx in ("", "Unknown", "Unassigned_Packera", "Senecio sp."):
            self.assertFalse(mod_xai._is_valid_packera_taxon(invalid_tx), f"Expected invalid taxon pattern: {invalid_tx}")

        # Verify lightweight regex fallback when harvester raises an exception
        with mock_patch("scripts.core.harvester.standardize_packera_taxon", side_effect=RuntimeError("offline")):
            self.assertEqual(standardize_packera_taxon("Packera cana (Hook.)"), "Packera cana")
            self.assertEqual(standardize_packera_taxon("Packera x hybrida"), "Packera × hybrida")

    def test_neutralize_mounting_paper(self):
        """Test vectorized Otsu background neutralization replaces paper with RGB(128,128,128)."""
        arr = np.full((100, 100, 3), fill_value=[225, 220, 200], dtype=np.uint8)
        y, x = np.ogrid[:100, :100]
        mask_leaf = ((x - 50) ** 2 + (y - 50) ** 2) <= 25 ** 2
        arr[mask_leaf] = [45, 95, 35]

        img = Image.fromarray(arr)
        neutralized_img = neutralize_mounting_paper(img)
        res_arr = np.array(neutralized_img)

        self.assertEqual(res_arr.shape, (100, 100, 3))
        self.assertEqual(res_arr.dtype, np.uint8)
        np.testing.assert_array_equal(res_arr[0, 0], [128, 128, 128])
        np.testing.assert_array_equal(res_arr[50, 50], [45, 95, 35])

    def test_224x224_low_pass_macro_architectural_pooling(self):
        """Verify Gaussian anti-aliasing (sigma=1.0) and strict 224x224 downsampling attenuate fine pubescence."""
        # Create a 448x448 image with a macro rosette disk and high-frequency 1-px trichome checkerboard noise
        arr = np.full((448, 448, 3), fill_value=235, dtype=np.uint8)
        y, x = np.ogrid[:448, :448]
        rosette_mask = ((x - 224) ** 2 + (y - 224) ** 2) <= 140 ** 2
        base_green = np.full((448, 448, 3), fill_value=[50, 100, 40], dtype=np.int16)
        checker = (((x + y) % 2) * 80 - 40)[:, :, None]
        noisy_rosette = np.clip(base_green + checker, 0, 255).astype(np.uint8)
        arr[rosette_mask] = noisy_rosette[rosette_mask]

        raw_center_std = float(np.std(arr[180:268, 180:268, 1].astype(np.float32)))
        macro = preprocess_macro_rosette_crop(arr, target_size=(224, 224), sigma=1.0)

        self.assertEqual(macro.shape, (224, 224, 3))
        self.assertEqual(macro.dtype, np.uint8)
        filtered_center_std = float(np.std(macro[90:134, 90:134, 1].astype(np.float32)))
        self.assertLess(
            filtered_center_std,
            raw_center_std * 0.25,
            "Gaussian anti-aliasing + 224x224 downsampling must suppress high-frequency epidermal texture.",
        )

    def test_foreground_masked_token_averaging(self):
        """Verify 16x16 foreground mask excludes background paper tokens from the 768-d L2-normalized embedding."""
        arr = np.full((224, 224, 3), fill_value=[242, 239, 230], dtype=np.uint8)
        # Plant tissue strictly in the top-left quadrant (rows 0..112, cols 0..112 -> 8x8 patch grid)
        arr[:112, :112] = [40, 95, 35]

        fg_mask_224, fg_grid_16 = compute_foreground_tissue_mask(arr, grid_size=(16, 16))
        self.assertEqual(fg_mask_224.shape, (224, 224))
        self.assertEqual(fg_grid_16.shape, (16, 16))
        self.assertGreater(float(fg_grid_16[:8, :8].mean()), 0.85)
        self.assertLess(float(fg_grid_16[8:, 8:].mean()), 0.15)

        # Construct synthetic 16x16=256 patch tokens where foreground tokens have [1, 0, ...]
        # and background paper tokens have [0, 100, ...]
        tokens = np.zeros((16, 16, 768), dtype=np.float32)
        tokens[:8, :8, 0] = 5.0
        tokens[8:, :, 1] = 100.0
        tokens[:, 8:, 1] = 100.0
        tokens_flat = tokens.reshape(256, 768)

        pooled = foreground_masked_patch_pooling(tokens_flat, fg_grid_16)
        self.assertEqual(pooled.shape, (768,))
        self.assertAlmostEqual(float(np.linalg.norm(pooled)), 1.0, places=5)
        self.assertGreater(float(pooled[0]), 0.98, "Pooled embedding must reflect foreground tissue tokens.")
        self.assertLess( abs(float(pooled[1])), 0.05, "Background paper tokens must be excluded from pooled embedding.")

    def test_classify_foliar_surface_orientation_adaxial_and_abaxial(self):
        """Verify adaxial vs. abaxial foliar orientation classification and <0.05s execution speed."""
        arr_adax = np.full((120, 120, 3), fill_value=128, dtype=np.uint8)
        y, x = np.ogrid[:120, :120]
        mask_leaf = ((x - 60) ** 2 + (y - 60) ** 2) <= 30 ** 2
        arr_adax[mask_leaf] = [45, 95, 35]
        res_adax = classify_foliar_surface_orientation(arr_adax)
        self.assertEqual(res_adax["surface_orientation"], "adaxial")

        arr_abax = np.full((120, 120, 3), fill_value=128, dtype=np.uint8)
        arr_abax[mask_leaf] = [225, 230, 220]
        res_abax = classify_foliar_surface_orientation(arr_abax)
        self.assertEqual(res_abax["surface_orientation"], "abaxial")

        t0 = time.perf_counter()
        for _ in range(25):
            _ = classify_foliar_surface_orientation(arr_adax)
        self.assertLess((time.perf_counter() - t0) / 25.0, 0.05)

    def test_dinov2_feature_extraction_synthetic(self):
        """Test DINOv2 feature extraction returns (N, 768) L2-normalized embeddings when demo=True."""
        synthetic_records = [
            {
                "catalogNumber": f"SYNTH_{i:03d}",
                "patch_path": f"/nonexistent/path/patch_{i}.jpg",
                "label_idx": i % len(TEST_SPECIES),
                "taxon": TEST_SPECIES[i % len(TEST_SPECIES)],
            }
            for i in range(12)
        ]

        with self.assertRaises(RuntimeError):
            extract_dinov2_embeddings(
                records=synthetic_records,
                model_name="dinov2_vitb14_reg",
                device="cpu",
                batch_size=4,
                demo=False,
            )

        features, labels, cat_nums = extract_dinov2_embeddings(
            records=synthetic_records,
            model_name="dinov2_vitb14_reg",
            device="cpu",
            batch_size=4,
            demo=True,
        )

        self.assertEqual(features.shape, (12, 768))
        self.assertEqual(len(labels), 12)
        self.assertEqual(len(cat_nums), 12)
        norms = np.linalg.norm(features, axis=1)
        np.testing.assert_allclose(norms, np.ones(12), atol=1e-4)

    def test_authority_anchored_institution_grouped_oof_calibration(self):
        """Verify Gold/Silver anchor training, institutionCode-grouped CV, and passive Bronze prediction."""
        rng = np.random.default_rng(42)
        n_anchors = 24
        n_bronze = 8
        n_classes = 3

        anchor_labels = np.array([i % n_classes for i in range(n_anchors)], dtype=int)
        X_anc = rng.normal(0.0, 0.15, size=(n_anchors, 768)).astype(np.float32)
        for i in range(n_anchors):
            X_anc[i, anchor_labels[i] * 50 : (anchor_labels[i] + 1) * 50] += 4.0

        # Create corrupted Tier_3_Bronze samples with inverted features that would warp decision boundary if trained on
        bronze_labels = np.array([i % n_classes for i in range(n_bronze)], dtype=int)
        X_brz = rng.normal(0.0, 0.15, size=(n_bronze, 768)).astype(np.float32)
        for i in range(n_bronze):
            true_c = (bronze_labels[i] + 1) % n_classes
            X_brz[i, true_c * 50 : (true_c + 1) * 50] += 4.0

        features = np.vstack([X_anc, X_brz])
        labels = np.concatenate([anchor_labels, bronze_labels])
        tiers = ["Tier_1_Gold" if i % 2 == 0 else "Tier_2_Silver" for i in range(n_anchors)] + ["Tier_3_Bronze"] * n_bronze
        inst_codes = [["NCU", "US", "NY", "MO"][(i // 3) % 4] for i in range(len(labels))]

        pred_probs, acc, f1 = compute_out_of_fold_probabilities(
            features=features,
            labels=labels,
            n_splits=3,
            random_state=42,
            determiner_tiers=tiers,
            institution_codes=inst_codes,
        )

        self.assertEqual(pred_probs.shape, (n_anchors + n_bronze, n_classes))
        self.assertTrue(np.allclose(pred_probs.sum(axis=1), 1.0, atol=1e-4))
        self.assertGreaterEqual(acc, 0.85, "Gold/Silver anchor OOF accuracy should remain high despite corrupted Bronze records.")

        # Passive Bronze records with swapped features should be strongly predicted as their true visual class (true_c)
        for i in range(n_bronze):
            idx = n_anchors + i
            expected_cls = (bronze_labels[i] + 1) % n_classes
            self.assertEqual(int(np.argmax(pred_probs[idx])), expected_cls)

    def test_compute_vit_patch_saliency_and_foreground_attention_ratio(self):
        """Verify 16x16 -> 224x224 ViT patch saliency and Pointing Game foreground_attention_ratio."""
        import torch
        from sklearn.linear_model import LogisticRegression

        class DummyDINOv2(torch.nn.Module):
            def forward_features(self, x: torch.Tensor):
                b = x.shape[0]
                patch_tokens = torch.zeros((b, 256, 768), dtype=torch.float32)
                patch_tokens[:, :64, 0] = 4.0
                patch_tokens[:, 64:, 0] = -2.0
                cls_token = torch.zeros((b, 768), dtype=torch.float32)
                cls_token[:, 0] = 1.0
                return {"x_norm_patchtokens": patch_tokens, "x_norm_clstoken": cls_token}

        X_train = np.random.RandomState(42).randn(16, 768).astype(np.float32)
        y_train = np.array([i % 4 for i in range(16)])
        clf = LogisticRegression(max_iter=200, random_state=42)
        clf.fit(X_train, y_train)
        clf.coef_[:] = 0.0
        clf.coef_[1, 0] = 5.0

        arr = np.full((224, 224, 3), fill_value=128, dtype=np.uint8)
        heatmap = mod_xai.compute_vit_patch_saliency(
            img_arr=arr,
            model=DummyDINOv2(),
            classifier=clf,
            target_class_idx=1,
            device="cpu",
        )
        self.assertEqual(heatmap.shape, (224, 224))
        self.assertGreater(float(heatmap[:56, :].mean()), float(heatmap[112:, :].mean()))

        # Foreground mask covering top rows where attribution is concentrated -> high ratio (>= 0.65)
        fg_top = np.zeros((224, 224), dtype=np.uint8)
        fg_top[:80, :] = 1
        ratio_high = compute_foreground_attention_ratio(heatmap, fg_top)
        self.assertGreaterEqual(ratio_high, 0.65)

        # Foreground mask covering bottom rows while attribution is on top background -> low ratio (< 0.65)
        fg_bottom = np.zeros((224, 224), dtype=np.uint8)
        fg_bottom[140:, :] = 1
        ratio_low = compute_foreground_attention_ratio(heatmap, fg_bottom)
        self.assertLess(ratio_low, 0.65)

        # Binary (2-class) LogisticRegression where coef_.shape == (1, 768) and classes_ == [0, 1]
        y_bin = np.array([i % 2 for i in range(16)])
        clf_bin = LogisticRegression(max_iter=200, random_state=42).fit(X_train, y_bin)
        self.assertEqual(clf_bin.coef_.shape, (1, 768))
        clf_bin.coef_[:] = 0.0
        clf_bin.coef_[0, 0] = 5.0
        hm_cls1 = mod_xai.compute_vit_patch_saliency(arr, model=DummyDINOv2(), classifier=clf_bin, target_class_idx=1, device="cpu")
        hm_cls0 = mod_xai.compute_vit_patch_saliency(arr, model=DummyDINOv2(), classifier=clf_bin, target_class_idx=0, device="cpu")
        self.assertEqual(hm_cls1.shape, (224, 224))
        self.assertEqual(hm_cls0.shape, (224, 224))
        self.assertGreater(float(hm_cls1[:56, :].mean()), float(hm_cls1[112:, :].mean()))
        self.assertGreater(float(hm_cls0[112:, :].mean()), float(hm_cls0[:56, :].mean()))

    def test_confident_learning_audit_schema_and_ood(self):
        """Verify Confident Learning audit exports all required Step 07 schema columns and OOD scores."""
        pred_probs = np.array([
            [0.05, 0.92, 0.02, 0.01],
            [0.05, 0.90, 0.03, 0.02],
            [0.20, 0.20, 0.60, 0.00],
            [0.88, 0.02, 0.00, 0.10],
        ])
        labels = np.array([0, 1, 2, 3])
        rng = np.random.default_rng(42)
        feats = rng.normal(0.0, 1.0, size=(4, 768)).astype(np.float32)

        records_df = pd.DataFrame({
            "catalogNumber": ["CAT001", "CAT002", "CAT003", "CAT004"],
            "patch_path": ["/p1", "/p2", "/p3", "/p4"],
            "species_raw": TEST_SPECIES,
            "label_idx": [0, 1, 2, 3],
            "determiner_tier": ["Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze", "Tier_1_Gold"],
            "foreground_attention_ratio": [0.88, 0.91, 0.79, 0.52],
        })

        audit = run_confident_learning_audit(
            pred_probs=pred_probs,
            labels=labels,
            records_df=records_df,
            class_names=TEST_SPECIES,
            error_threshold=0.85,
            features=feats,
        )

        required_step07_cols = [
            "catalogNumber",
            "plant_individual_id",
            "c_error",
            "predicted_label",
            "vision_predicted_label",
            "vision_confidence",
            "vision_second_prob",
            "vision_margin",
            "vision_ood_score",
            "vision_is_ood",
            "foreground_attention_ratio",
            "is_artifact_driven",
            "is_label_issue",
            "is_label_corrupted",
            "tomentum_taxon_discordance",
        ]
        for col in required_step07_cols:
            self.assertIn(col, audit.columns, f"Missing required Step 07 contract column: {col}")

        self.assertEqual(audit["plant_individual_id"].tolist(), [0, 0, 0, 0])
        self.assertEqual(audit["vision_is_ood"].dtype, np.dtype(bool))
        self.assertEqual(audit["is_artifact_driven"].dtype, np.dtype(bool))
        self.assertTrue(audit.loc[0, "is_label_corrupted"])
        self.assertFalse(audit.loc[1, "is_label_corrupted"])
        self.assertAlmostEqual(float(audit.loc[0, "vision_margin"]), 0.92 - 0.05, places=4)
        self.assertTrue(bool(audit.loc[3, "is_artifact_driven"]))
        self.assertIn("Downweight", audit.loc[3, "triage_action"])

    def test_rosettes_dense_ingestion_multi_crop_aggregation_and_dynamic_indexing(self):
        """Verify rosettes_dense ingestion, multi-crop individual aggregation, and dynamic 0..K-1 taxon indexing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            rosette_dir = tmp_path / "cropped_patches"
            dense_dir = rosette_dir / "rosettes_dense"
            dense_dir.mkdir(parents=True)

            # Create multiple crops for (CATANON, p0) and single crops for CATDUB, CATPAUP, CATPLAT
            for fname in (
                "CATANON_p0_rosette.jpg",
                "CATANON_p0_rosette_crop2.jpg",
                "CATDUB_p0_rosette.jpg",
                "CATPAUP_p0_rosette.jpg",
                "CATPLAT_p0_rosette.jpg",
            ):
                img = Image.new("RGB", (300, 300), color=(50, 100, 40))
                img.save(dense_dir / fname)

            # Filter voucher table to 3 taxa (excluding Packera plattensis / CATPLAT)
            vouchers_csv = tmp_path / "curated_vouchers.csv"
            pd.DataFrame({
                "catalogNumber": ["CATANON", "CATDUB", "CATPAUP"],
                "species_raw": ["Packera anonyma", "Packera dubia", "Packera paupercula"],
                "determiner_tier": ["Tier_1_Gold", "Tier_1_Gold", "Tier_2_Silver"],
                "institutionCode": ["NCU", "US", "NY"],
            }).to_csv(vouchers_csv, index=False)

            df, class_to_idx = mod_xai.load_and_link_rosette_patches(
                rosette_dir=rosette_dir,
                vouchers_csv=vouchers_csv,
                demo=False,
            )

            self.assertEqual(len(df), 3)
            self.assertNotIn("CATPLAT", set(df["catalogNumber"]))
            self.assertEqual(
                class_to_idx,
                {"Packera anonyma": 0, "Packera dubia": 1, "Packera paupercula": 2},
            )
            self.assertEqual(sorted(df["label_idx"].unique().tolist()), [0, 1, 2])
            anon_row = df[df["catalogNumber"] == "CATANON"].iloc[0]
            self.assertEqual(int(anon_row["n_rosette_patches"]), 2)

            # Verify multi-crop sheet aggregates to 1 specimen-level L2-normalized embedding per catalogNumber
            feats, lbls, cats = extract_dinov2_embeddings(df.to_dict("records"), demo=True)
            self.assertEqual(feats.shape, (3, 768))
            self.assertEqual(cats, ["CATANON", "CATDUB", "CATPAUP"])
            np.testing.assert_allclose(np.linalg.norm(feats, axis=1), np.ones(3), atol=1e-4)

            # Verify generate_gradcam_panel plots probability bars only for active_taxa present in the run
            from unittest.mock import patch as mock_patch
            captured_yticks = []
            orig_savefig = mod_xai.plt.savefig
            def _spy_savefig(*args, **kwargs):
                fig = mod_xai.plt.gcf()
                if len(fig.axes) >= 3:
                    captured_yticks.append([t.get_text() for t in fig.axes[2].get_yticklabels()])
                return orig_savefig(*args, **kwargs)
            with mock_patch.object(mod_xai.plt, "savefig", side_effect=_spy_savefig):
                mod_xai.generate_gradcam_panel(
                    df.to_dict("records"),
                    tmp_path / "gradcam_3taxa.png",
                    demo=True,
                )
            self.assertTrue(captured_yticks)
            self.assertEqual(captured_yticks[0], ["P. anonyma", "P. dubia", "P. paupercula"])

    def test_dynamic_non_dubia_species_audit_tomentum_discordance_and_figures(self):
        """Verify dynamic non-dubia taxa indexing, empirical chromatic tomentum discordance, and Step 05 figures."""
        from unittest.mock import patch as mock_patch

        non_dubia_taxa = ["Packera aurea", "Packera obovata", "Packera schweinitziana"]
        palette = mod_xai._get_dynamic_okabe_ito_map(non_dubia_taxa)
        self.assertEqual(len(set(palette.values())), 3)
        self.assertEqual(palette["Packera aurea"], "#0072B2")
        self.assertEqual(palette["Packera obovata"], "#D55E00")
        self.assertEqual(palette["Packera schweinitziana"], "#009E73")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            rosette_dir = tmp_path / "cropped_patches"
            dense_dir = rosette_dir / "rosettes_dense"
            dense_dir.mkdir(parents=True)

            voucher_rows = []
            # 3 species x 4 vouchers = 12 sheets (zero Packera dubia records)
            for idx_sp, sp_name in enumerate(non_dubia_taxa):
                for j in range(4):
                    cat = f"NONDUB_{idx_sp}_{j}"
                    is_abax = (idx_sp == 0 and j == 3)
                    # Bright desaturated RGB for abaxial sheet; darker green for adaxial sheets;
                    # slightly brighter adaxial baseline for P. schweinitziana (idx_sp == 2)
                    rgb_color = (228, 232, 224) if is_abax else ((155, 165, 150) if idx_sp == 2 else (45, 95, 35))
                    Image.new("RGB", (240, 240), color=rgb_color).save(dense_dir / f"{cat}_p0_rosette.jpg")
                    voucher_rows.append({
                        "catalogNumber": cat,
                        "species_raw": sp_name,
                        "determiner_tier": "Tier_1_Gold" if j < 3 else "Tier_2_Silver",
                        "institutionCode": ["NCU", "NY", "US", "MO"][j],
                    })

            vouchers_csv = tmp_path / "curated_vouchers.csv"
            pd.DataFrame(voucher_rows).to_csv(vouchers_csv, index=False)

            df, class_to_idx = mod_xai.load_and_link_rosette_patches(
                rosette_dir=rosette_dir,
                vouchers_csv=vouchers_csv,
                demo=False,
            )
            self.assertEqual(
                class_to_idx,
                {"Packera aurea": 0, "Packera obovata": 1, "Packera schweinitziana": 2},
            )
            self.assertNotIn("Packera dubia", set(df["taxon"]))
            self.assertNotIn("Packera dubia", class_to_idx)

            # Construct probabilities where the abaxial P. aurea sheet (NONDUB_0_3, label 0) is strongly
            # mispredicted as higher-luminance P. schweinitziana (class 2) with c_error = 0.94 > 0.85
            labels = df["label_idx"].to_numpy(dtype=int)
            pred_probs = np.zeros((len(df), 3), dtype=np.float64)
            for i, lbl in enumerate(labels):
                if df.loc[i, "catalogNumber"] == "NONDUB_0_3":
                    pred_probs[i] = [0.04, 0.02, 0.94]
                else:
                    pred_probs[i] = [0.04, 0.04, 0.04]
                    pred_probs[i, lbl] = 0.92

            active_names = ["Packera aurea", "Packera obovata", "Packera schweinitziana"]
            audit = run_confident_learning_audit(
                pred_probs=pred_probs,
                labels=labels,
                records_df=df,
                class_names=active_names,
                error_threshold=0.85,
            )

            for col_chk in ("species_standardized", "given_label", "predicted_label", "vision_predicted_label"):
                self.assertNotIn("Packera dubia", set(audit[col_chk]))

            abax_row = audit[audit["catalogNumber"] == "NONDUB_0_3"].iloc[0]
            self.assertEqual(abax_row["surface_orientation"], "abaxial")
            self.assertTrue(bool(abax_row["tomentum_taxon_discordance"]))
            adax_rows = audit[audit["catalogNumber"] != "NONDUB_0_3"]
            self.assertFalse(bool(adax_rows["tomentum_taxon_discordance"].any()))

            # Also verify graceful constant-luminance (0.0) fallback on synthetic unit-test fixtures
            synth_df = pd.DataFrame({
                "catalogNumber": ["S1", "S2"],
                "species_raw": ["Packera anonyma", "Packera dubia"],
                "surface_orientation": ["abaxial", "adaxial"],
                "mean_luminance": [0.0, 0.0],
                "mean_saturation": [0.0, 0.0],
            })
            synth_audit = run_confident_learning_audit(
                pred_probs=np.array([[0.05, 0.95], [0.05, 0.95]]),
                labels=np.array([0, 1]),
                records_df=synth_df,
                class_names=["Packera anonyma", "Packera dubia"],
                error_threshold=0.85,
            )
            self.assertTrue(bool(synth_audit.loc[0, "tomentum_taxon_discordance"]))
            self.assertFalse(bool(synth_audit.loc[1, "tomentum_taxon_discordance"]))

            # Verify generate_step05_publication_figures renders without injecting Packera dubia
            captured_fig05b_labels = []
            orig_fig_savefig = mod_xai.plt.Figure.savefig

            def _spy_fig_savefig(fig_self, *args, **kwargs):
                if len(fig_self.axes) >= 2:
                    xt = [t.get_text() for t in fig_self.axes[0].get_xticklabels()]
                    _, leg_lbls = fig_self.axes[1].get_legend_handles_labels()
                    captured_fig05b_labels.append((xt, leg_lbls))
                return orig_fig_savefig(fig_self, *args, **kwargs)

            with mock_patch.dict("os.environ", {"PACKERA_RUN_DIR": ""}, clear=False):
                with mock_patch.object(mod_xai.plt.Figure, "savefig", _spy_fig_savefig):
                    fig_paths = mod_xai.generate_step05_publication_figures(
                        audit, tmp_path / "figures" / "GradCAM_audit_panel.png"
                    )
            self.assertTrue(fig_paths["pdf"].exists())
            self.assertTrue(fig_paths["png"].exists())
            self.assertTrue(captured_fig05b_labels)
            xt_labels, leg_labels = captured_fig05b_labels[0]
            self.assertEqual(xt_labels, ["P. aurea", "P. obovata", "P. schweinitziana"])
            self.assertEqual(leg_labels, ["P. aurea", "P. obovata", "P. schweinitziana"])
            self.assertNotIn("P. dubia", xt_labels)

    def test_step02_label_blind_dissection_and_dynamic_taxon_figures(self):
        """Verify Phase 2 label-blind Dissected/Bipinnatisect routing and non-dubia dynamic Step 02 figures."""
        import cv2
        from unittest.mock import patch as mock_patch
        from scripts.vision.lm2_geometry_utils import classify_and_route_leaf

        step02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")

        # 1. _standardize_packera_species_py preserves valid binomials/morphotypes without forcing Packera dubia
        self.assertEqual(step02._standardize_packera_species_py(None), "Unassigned_Packera")
        self.assertEqual(step02._standardize_packera_species_py("Unknown"), "Unassigned_Packera")
        self.assertEqual(step02._standardize_packera_species_py("Morphotype_K1"), "Morphotype_K1")
        self.assertEqual(step02._standardize_packera_species_py("Packera aurea (L.) Á.Löve & D.Löve"), "Packera aurea")
        with mock_patch(
            "pygbif.species.name_backbone",
            side_effect=lambda name="", **kw: (
                {"matchType": "EXACT", "species": "Packera schweinitziana", "usageKey": 102}
                if "robbinsii" in str(name)
                else {"matchType": "NONE"}
            ),
        ), mock_patch("pygbif.species.name_usage", return_value={"results": []}):
            self.assertEqual(step02._standardize_packera_species_py("Senecio robbinsii Oakes"), "Packera schweinitziana")

        # 2. Label-blind pinnatifid (0.48 <= Solidity < 0.72 -> Tier 1 (Dissected)) across any species name
        pinn_mask = np.zeros((400, 400), dtype=np.uint8)
        cx, cy = 200, 200
        cv2.ellipse(pinn_mask, (cx, cy - 70), (35, 45), 0, 0, 360, 255, -1)
        cv2.line(pinn_mask, (cx, cy - 100), (cx, cy + 100), 255, thickness=16)
        for y_off in (-20, 20, 60):
            cv2.ellipse(pinn_mask, (cx - 35, cy + y_off), (35, 12), -15, 0, 360, 255, -1)
            cv2.ellipse(pinn_mask, (cx + 35, cy + y_off), (35, 12), 15, 0, 360, 255, -1)
        res_aurea = classify_and_route_leaf(pinn_mask, taxon="Packera aurea")
        res_obovata = classify_and_route_leaf(pinn_mask, taxon="Packera obovata", taxa_with_lyrate_tendency=[])
        self.assertEqual(res_aurea.tier, "Tier 1 (Dissected)")
        self.assertEqual(res_obovata.tier, "Tier 1 (Dissected)")
        self.assertEqual(res_aurea.canonical_tier, "Tier_1_Direct")

        # 3. Label-blind bipinnatisect (0.18 <= Solidity < 0.48 with >= 4 sinuses -> Tier 1 (Bipinnatisect))
        bipinn_mask = np.zeros((400, 400), dtype=np.uint8)
        cv2.line(bipinn_mask, (cx, 80), (cx, 320), 255, thickness=8)
        for y_off in (-90, -45, 0, 45, 90):
            cv2.line(bipinn_mask, (cx - 58, cy + y_off), (cx + 58, cy + y_off), 255, thickness=7)
        res_mille = classify_and_route_leaf(bipinn_mask, taxon="Packera millefolium")
        self.assertGreaterEqual(res_mille.solidity, 0.18)
        self.assertLess(res_mille.solidity, 0.48)
        self.assertEqual(res_mille.tier, "Tier 1 (Bipinnatisect)")
        self.assertEqual(res_mille.canonical_tier, "Tier_1_Direct")

        # 4. Step 02 figures dynamically title & style non-dubia species without injecting Packera dubia
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_fig_dir = Path(tmpdir) / "outputs" / "figures"
            df_multi = pd.DataFrame({
                "catalogNumber": [f"V{i}" for i in range(6)],
                "leaf_id": [1] * 6,
                "assigned_tier": ["Tier_1_Direct"] * 6,
                "species_standardized": [
                    "Packera aurea", "Packera aurea",
                    "Packera obovata", "Packera obovata",
                    "Packera schweinitziana", "Packera schweinitziana",
                ],
                "solidity": [0.85, 0.88, 0.90, 0.87, 0.62, np.nan],
                "ucs_score": [0.92] * 6,
                "aspect_ratio": [1.45, np.nan, 1.80, 1.90, 2.50, np.nan],
            })
            captured_title_a = []
            orig_save02 = step02._atomic_save_matplotlib_fig

            def _spy_save02(fig, out_base_path, dpi=300):
                if Path(out_base_path).name == "fig02a_segmentation_qc_distributions":
                    captured_title_a.append(fig.axes[0].get_title(loc="left") or fig.axes[0].get_title())
                return orig_save02(fig, out_base_path, dpi=dpi)

            with mock_patch.dict("os.environ", {"PACKERA_FOCAL_TAXON": ""}, clear=False):
                with mock_patch.object(step02, "_atomic_save_matplotlib_fig", side_effect=_spy_save02):
                    saved_multi = step02.generate_step02_publication_figures(
                        df_multi, pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), tmp_fig_dir
                    )
            self.assertTrue(saved_multi["fig02a_png"].exists())
            self.assertEqual(
                captured_title_a[-1],
                "A. Foliar Solidity vs. Blade Elongation ($N = 3$ Taxa)",
            )
            self.assertNotIn("dubia", captured_title_a[-1])

            with mock_patch.dict("os.environ", {"PACKERA_FOCAL_TAXON": "Packera aurea"}, clear=False):
                with mock_patch.object(step02, "_atomic_save_matplotlib_fig", side_effect=_spy_save02):
                    step02.generate_step02_publication_figures(
                        df_multi, pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), tmp_fig_dir
                    )
            self.assertEqual(
                captured_title_a[-1],
                r"A. Foliar Solidity vs. Blade Elongation ($\it{Packera\ aurea}$)",
            )
            self.assertNotIn("dubia", captured_title_a[-1])

            df_single = df_multi[df_multi["species_standardized"] == "Packera obovata"].copy()
            with mock_patch.dict("os.environ", {"PACKERA_FOCAL_TAXON": ""}, clear=False):
                with mock_patch.object(step02, "_atomic_save_matplotlib_fig", side_effect=_spy_save02):
                    step02.generate_step02_publication_figures(
                        df_single, pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), tmp_fig_dir
                    )
            self.assertEqual(
                captured_title_a[-1],
                r"A. Foliar Solidity vs. Blade Elongation ($\it{Packera\ obovata}$)",
            )
            self.assertNotIn("dubia", captured_title_a[-1])

    def test_step02_post_qc_rosette_masking_excludes_roots_stems_and_cauline_leaves(self):
        """Verify Step 02 runs crop_dense_rosette_tier4 post-QC, masks non-leaf pixels to [242,239,230], and saves companion mask."""
        import cv2
        step02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            vouchers_csv = tmp_path / "curated_vouchers.csv"
            pd.DataFrame({"catalogNumber": ["NCU_MASK1"], "species_raw": ["Packera aurea"]}).to_csv(vouchers_csv, index=False)
            weights_p = tmp_path / "dummy.pth"
            weights_p.write_bytes(b"0")
            sheet_p = tmp_path / "NCU_MASK1.jpg"

            # Synthetic 600x500 BGR sheet with dark roots at bottom and dark stem/cauline bract at top
            basal_mask = np.zeros((600, 500), dtype=np.uint8)
            cv2.ellipse(basal_mask, (250, 370), (65, 85), 0, 0, 360, 255, -1)

            cauline_mask = np.zeros((600, 500), dtype=np.uint8)
            pts_caul = np.array([[250, 35], [262, 100], [250, 175], [238, 100]], dtype=np.int32)
            cv2.fillPoly(cauline_mask, [pts_caul], 255)

            sheet_bgr = np.full((600, 500, 3), 215, dtype=np.uint8)
            sheet_bgr[464:580, 180:320] = (20, 30, 45)    # dark fibrous root mass right below basal leaf (inside 15% crop pad)
            sheet_bgr[cauline_mask > 0] = (25, 45, 30)    # upper stem/cauline bract
            sheet_bgr[basal_mask > 0] = (40, 105, 50)     # green basal leaf lamina
            cv2.imwrite(str(sheet_p), sheet_bgr)

            pipeline = step02.SegmentAndExtractPipeline(
                vouchers_csv=vouchers_csv,
                model_weights=weights_p,
                stage1_weights=weights_p,
                output_dir=tmp_path / "data",
                device="cpu",
                enable_geometric_qc=True,
            )
            preds = [
                ((280, 180, 460, 320), basal_mask, 0.96, 0),
                ((30, 235, 180, 265), cauline_mask, 0.88, 0),
            ]
            with mock_patch.object(pipeline.model_engine, "predict", return_value=preds):
                records, failure = pipeline.process_voucher("NCU_MASK1", sheet_p, taxon="Packera aurea")

            self.assertIsNone(failure)
            self.assertEqual(len(records), 1)
            rosette_jpg = pipeline.rosettes_dense_dir / "NCU_MASK1_p0_rosette.jpg"
            rosette_mask_png = pipeline.rosettes_dense_dir / "NCU_MASK1_p0_rosette_mask.png"
            self.assertTrue(rosette_jpg.exists())
            self.assertTrue(rosette_mask_png.exists())

            with Image.open(rosette_jpg) as _r_im:
                rosette_rgb = np.array(_r_im.convert("RGB"))
            with Image.open(rosette_mask_png) as _m_im:
                rmask = np.array(_m_im.convert("L"))
            self.assertEqual(rosette_rgb.shape[:2], rmask.shape[:2])
            # Crop height must cover only the basal leaf + 15% pad (~230 px), NOT the rejected upper cauline bract at y=35
            self.assertLess(rosette_rgb.shape[0], 300)
            # Corner / exterior root region outside dilated mask must be neutral herbarium paper [242, 239, 230]
            np.testing.assert_allclose(rosette_rgb[-2, 2], [242, 239, 230], atol=4)

    def test_tier1_mask_linking_gradcam_leaf_saliency_empirical_fg_ratio_and_oof_probs(self):
        """Verify Tier_1_Direct basal leaf mask linking, root/stem masking, empirical FG ratios, and true OOF posteriors."""
        import torch

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            data_dir = tmp_path / "data"
            basal_dir = data_dir / "cropped_patches" / "basal_leaves_raw"
            dense_dir = data_dir / "cropped_patches" / "rosettes_dense"
            tier1_mask_dir = data_dir / "masks" / "tier1_intact"
            tables_dir = data_dir / "tables"
            for d in (basal_dir, dense_dir, tier1_mask_dir, tables_dir):
                d.mkdir(parents=True, exist_ok=True)

            # Create 224x224 crops where top-left quadrant is leaf blade and bottom-right has dark brown roots
            for idx_c, (cat, sp) in enumerate([
                ("CAT_A", "Packera aurea"),
                ("CAT_B", "Packera obovata"),
                ("CAT_C", "Packera cymbalaria"),
            ]):
                rgb = np.full((224, 224, 3), fill_value=[235, 232, 222], dtype=np.uint8)
                rgb[14:110, 14:110] = [45, 110, 40]     # verified basal leaf lamina (top-left)
                rgb[126:215, 126:215] = [22, 18, 15]    # dark root/stem contaminant (bottom-right)
                Image.fromarray(rgb).save(basal_dir / f"{cat}_leaf1.jpg")

                exact_m = np.zeros((224, 224), dtype=np.uint8)
                exact_m[14:110, 14:110] = 255
                Image.fromarray(exact_m).save(tier1_mask_dir / f"{cat}_leaf1.png")

                # Also drop a companion rosette_mask.png in rosettes_dense to verify *_mask.png is never ingested as RGB
                Image.fromarray(exact_m).save(dense_dir / f"{cat}_p0_rosette_mask.png")

            # Verify 35%-padded tier1_intact mask window aligns cleanly to a tight 120x100 basal_leaves_raw RGB crop
            tight_rgb = np.full((120, 100, 3), fill_value=[30, 20, 15], dtype=np.uint8)
            padded_win_m = np.zeros((204, 170), dtype=np.uint8)
            padded_win_m[42:162, 35:135] = 255
            padded_mp = tmp_path / "padded_tier1.png"
            Image.fromarray(padded_win_m).save(padded_mp)
            masked_tight = mod_xai.apply_tier1_leaf_mask_to_rgb(tight_rgb, padded_mp, pad_px=0)
            self.assertEqual(masked_tight.shape, (120, 100, 3))
            np.testing.assert_array_equal(masked_tight[0, 0], [30, 20, 15])
            np.testing.assert_array_equal(masked_tight[119, 99], [30, 20, 15])

            vouchers_csv = tables_dir / "curated_vouchers.csv"
            pd.DataFrame({
                "catalogNumber": ["CAT_A", "CAT_B", "CAT_C"],
                "species_raw": ["Packera aurea", "Packera obovata", "Packera cymbalaria"],
                "determiner_tier": ["Tier_1_Gold", "Tier_1_Gold", "Tier_2_Silver"],
                "institutionCode": ["NCU", "NY", "US"],
            }).to_csv(vouchers_csv, index=False)

            manifest_csv = tables_dir / "extracted_leaf_manifest.csv"
            pd.DataFrame({
                "catalogNumber": ["CAT_A", "CAT_B", "CAT_C"],
                "leaf_id": [1, 1, 1],
                "assigned_tier": ["Tier_1_Direct", "Tier_1_Direct", "Tier_1_Direct"],
                "solidity": [0.89, 0.91, 0.87],
                "ucs_score": [0.94, 0.93, 0.92],
            }).to_csv(manifest_csv, index=False)

            df, _ = mod_xai.load_and_link_rosette_patches(
                rosette_dir=data_dir / "cropped_patches",
                vouchers_csv=vouchers_csv,
                demo=False,
            )
            self.assertEqual(len(df), 3)
            for _, row in df.iterrows():
                self.assertIn("basal_leaves_raw", str(row["patch_path"]))
                self.assertTrue(Path(str(row["mask_path"])).is_file())
                self.assertNotIn("_mask.png", str(row["patch_path"]))

            # Verify RosettePatchDataset masks out dark roots in bottom-right to [242,239,230] and zeroes fg_grid there
            ds = mod_xai.RosettePatchDataset(df.to_dict("records"))
            _tensor0, fg_grid0, _, _ = ds[0]
            self.assertGreater(float(fg_grid0[1:7, 1:7].mean()), 0.85)
            self.assertEqual(float(fg_grid0[10:15, 10:15].sum()), 0.0, "Dark root quadrant outside Tier_1_Direct mask must have 0 fg_grid weight.")

            # Verify empirical foreground_attention_ratio in extract_dinov2_embeddings across ALL vouchers (not hardcoded 0.85)
            class MockEnergyViT(torch.nn.Module):
                def forward_features(self, x: torch.Tensor):
                    b = x.shape[0]
                    pt = torch.zeros((b, 256, 768), dtype=torch.float32)
                    grid = pt.view(b, 16, 16, 768)
                    grid[:, 2:6, 2:6, 0] = 9.0    # energy 9*16 = 144 inside top-left leaf mask
                    grid[:, 11:13, 11:13, 0] = 3.0  # energy 3*4 = 12 outside leaf mask -> ratio = 144 / 156 = 0.9231
                    return {"x_norm_patchtokens": pt, "x_norm_clstoken": pt.mean(dim=1)}

            recs = df.to_dict("records")
            with mock_patch("torch.hub.load", return_value=MockEnergyViT()):
                feats, lbls, cats, fg_ratios = extract_dinov2_embeddings(recs, device="cpu", return_fg_ratios=True)
            self.assertEqual(len(fg_ratios), 3)
            for fg_v, r_d in zip(fg_ratios, recs):
                self.assertAlmostEqual(float(fg_v), 144.0 / 156.0, places=3)
                self.assertAlmostEqual(float(r_d["foreground_attention_ratio"]), 144.0 / 156.0, places=3)
                self.assertNotAlmostEqual(float(fg_v), 0.85, places=2)

            # Verify generate_gradcam_panel masks roots to [242,239,230], zeroes background saliency, and plots true oof_probs
            true_oof = [0.17, 0.59, 0.24]
            recs[0]["oof_probs"] = true_oof
            recs[0]["predicted_label"] = "Packera obovata"
            recs[0]["c_error"] = 0.83

            captured_panel = {}
            orig_savefig = mod_xai.plt.savefig

            def _spy_gradcam(*args, **kwargs):
                fig = mod_xai.plt.gcf()
                if len(fig.axes) >= 3:
                    im0 = np.asarray(fig.axes[0].images[0].get_array())[:, :, :3]
                    im1 = np.asarray(fig.axes[1].images[0].get_array())[:, :, :3]
                    bar_widths = [float(p.get_width()) for p in fig.axes[2].patches]
                    captured_panel["im0"] = im0
                    captured_panel["im1"] = im1
                    captured_panel["bars"] = bar_widths
                return orig_savefig(*args, **kwargs)

            active_tx = ["Packera aurea", "Packera obovata", "Packera cymbalaria"]
            with mock_patch.object(mod_xai.plt, "savefig", side_effect=_spy_gradcam):
                mod_xai.generate_gradcam_panel(
                    [recs[0]],
                    tmp_path / "gradcam_masked.png",
                    model=MockEnergyViT(),
                    active_taxa=active_tx,
                    class_names=active_tx,
                    device="cpu",
                )

            self.assertIn("im0", captured_panel)
            # Bottom-right root region (e.g. [180, 180]) must be masked to neutral paper [242, 239, 230] in Panel [1] and Panel [2]
            np.testing.assert_allclose(captured_panel["im0"][180, 180], [242, 239, 230], atol=3)
            np.testing.assert_array_equal(captured_panel["im1"][180, 180], captured_panel["im0"][180, 180])
            self.assertAlmostEqual(float(recs[0]["foreground_attention_ratio"]), 144.0 / 156.0, places=2)
            # Panel [3] bars must match true calibrated oof_probs * 100 = [17.0, 59.0, 24.0] (not synthetic equal-split bars!)
            np.testing.assert_allclose(captured_panel["bars"], [17.0, 59.0, 24.0], atol=0.2)

    def test_intraspecific_mode_morphotype_cluster_remapping_in_main(self):
        """Verify --intraspecific-mode and --morphometrics-csv audit Morphotype_K1..Kn clusters in single-species runs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            dense_dir = tmp_path / "data" / "cropped_patches" / "rosettes_dense"
            tables_dir = tmp_path / "data" / "tables"
            dense_dir.mkdir(parents=True)
            tables_dir.mkdir(parents=True)

            cats = [f"INTRA_{i}" for i in range(6)]
            for idx_i, c in enumerate(cats):
                col = (45, 105, 40) if idx_i < 3 else (90, 130, 55)
                Image.new("RGB", (224, 224), color=col).save(dense_dir / f"{c}_p0_rosette.jpg")

            vouchers_csv = tables_dir / "curated_vouchers.csv"
            pd.DataFrame({
                "catalogNumber": cats,
                "species_raw": ["Packera aurea"] * 6,
                "determiner_tier": ["Tier_1_Gold"] * 6,
                "institutionCode": ["NCU", "NY", "US", "NCU", "NY", "US"],
            }).to_csv(vouchers_csv, index=False)

            morph_csv = tables_dir / "morphometric_flags.csv"
            pd.DataFrame({
                "catalogNumber": cats,
                "plant_individual_id": [0] * 6,
                "gmm_cluster": [1, 1, 1, 2, 2, 2],
            }).to_csv(morph_csv, index=False)

            out_csv = tables_dir / "label_noise_audit.csv"
            out_fig = tmp_path / "outputs" / "figures" / "GradCAM_audit_panel.png"

            with mock_patch("torch.hub.load", return_value=mod_xai._DemoViTFeatureExtractor()):
                with mock_patch("sys.argv", [
                    "05_cleanlab_vision_xai.py",
                    "--rosette-dir", str(tmp_path / "data" / "cropped_patches"),
                    "--vouchers-csv", str(vouchers_csv),
                    "--morphometrics-csv", str(morph_csv),
                    "--intraspecific-mode",
                    "--output-csv", str(out_csv),
                    "--output-figure", str(out_fig),
                    "--no-export-figures",
                ]):
                    mod_xai.main()

            self.assertTrue(out_csv.exists())
            exported = pd.read_csv(out_csv)
            self.assertEqual(len(exported), 6)
            self.assertEqual(set(exported["given_label"]), {"Morphotype_K1", "Morphotype_K2"})
            self.assertTrue(set(exported["predicted_label"]).issubset({"Morphotype_K1", "Morphotype_K2"}))
            self.assertNotIn("oof_probs", exported.columns)

    def test_multi_individual_sheet_composite_key_preservation_and_manifest_matching(self):
        """Verify (catalogNumber, plant_individual_id) preservation across multi-plant sheets in Step 05."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            data_dir = tmp_path / "data"
            basal_dir = data_dir / "cropped_patches" / "basal_leaves_raw"
            dense_dir = data_dir / "cropped_patches" / "rosettes_dense"
            tier1_mask_dir = data_dir / "masks" / "tier1_intact"
            tables_dir = data_dir / "tables"
            for d in (basal_dir, dense_dir, tier1_mask_dir, tables_dir):
                d.mkdir(parents=True, exist_ok=True)

            # Multi-individual sheet CAT_MULTI (p0 has 2 Tier_1_Direct leaves; p1 has 1 Tier_1_Direct leaf)
            for stem, rgb_col in (
                ("CAT_MULTI_p0_leaf1", (45, 110, 42)),
                ("CAT_MULTI_p0_leaf2", (48, 112, 44)),
                ("CAT_MULTI_p1_leaf1", (95, 130, 52)),
            ):
                Image.new("RGB", (224, 224), color=rgb_col).save(basal_dir / f"{stem}.jpg")
                m_arr = np.zeros((224, 224), dtype=np.uint8)
                m_arr[20:190, 20:190] = 255
                Image.fromarray(m_arr).save(tier1_mask_dir / f"{stem}.png")

            for r_stem, r_col in (
                ("CAT_MULTI_p0_rosette", (45, 110, 42)),
                ("CAT_MULTI_p1_rosette", (95, 130, 52)),
                ("CAT_S1_p0_rosette", (50, 105, 45)),
                ("CAT_S2_p0_rosette", (92, 128, 50)),
            ):
                Image.new("RGB", (224, 224), color=r_col).save(dense_dir / f"{r_stem}.jpg")

            # curated_vouchers.csv intentionally omits plant_individual_id (sheet-level table)
            vouchers_csv = tables_dir / "curated_vouchers.csv"
            pd.DataFrame({
                "catalogNumber": ["CAT_MULTI", "CAT_S1", "CAT_S2"],
                "species_raw": ["Packera aurea", "Packera aurea", "Packera obovata"],
                "determiner_tier": ["Tier_1_Gold", "Tier_1_Gold", "Tier_1_Gold"],
                "institutionCode": ["NCU", "NY", "US"],
            }).to_csv(vouchers_csv, index=False)

            manifest_csv = tables_dir / "extracted_leaf_manifest.csv"
            pd.DataFrame({
                "catalogNumber": ["CAT_MULTI", "CAT_MULTI", "CAT_MULTI"],
                "plant_individual_id": [0, 0, 1],
                "leaf_id": [1, 2, 1],
                "assigned_tier": ["Tier_1_Direct", "Tier_1_Direct", "Tier_1_Direct"],
                "solidity": [0.90, 0.88, 0.91],
                "ucs_score": [0.95, 0.92, 0.94],
            }).to_csv(manifest_csv, index=False)

            df, _ = mod_xai.load_and_link_rosette_patches(
                rosette_dir=data_dir / "cropped_patches",
                vouchers_csv=vouchers_csv,
                demo=False,
                manifest_csv=manifest_csv,
            )
            self.assertEqual(len(df), 4)
            self.assertTrue(pd.api.types.is_integer_dtype(df["plant_individual_id"]))
            observed_keys = list(zip(df["catalogNumber"], df["plant_individual_id"]))
            self.assertEqual(
                observed_keys,
                [("CAT_MULTI", 0), ("CAT_MULTI", 1), ("CAT_S1", 0), ("CAT_S2", 0)],
            )

            # Verify Tier_1_Direct leaves matched strictly per (catalogNumber, plant_individual_id)
            row_p0 = df[(df["catalogNumber"] == "CAT_MULTI") & (df["plant_individual_id"] == 0)].iloc[0]
            row_p1 = df[(df["catalogNumber"] == "CAT_MULTI") & (df["plant_individual_id"] == 1)].iloc[0]
            self.assertEqual(len(row_p0["tier1_leaf_patches"]), 2)
            self.assertTrue(all("CAT_MULTI_p0_leaf" in p for p in row_p0["tier1_leaf_patches"]))
            self.assertEqual(len(row_p1["tier1_leaf_patches"]), 1)
            self.assertIn("CAT_MULTI_p1_leaf1", row_p1["tier1_leaf_patches"][0])

            # Verify multi-patch embedding aggregation preserves both individuals on CAT_MULTI
            with mock_patch("torch.hub.load", return_value=mod_xai._DemoViTFeatureExtractor()):
                feats, lbls, key_ord, fg_ratios = extract_dinov2_embeddings(
                    df.to_dict("records"), device="cpu", return_fg_ratios=True
                )
            self.assertEqual(feats.shape, (4, 768))
            self.assertEqual(len(lbls), 4)
            self.assertEqual(len(fg_ratios), 4)
            self.assertEqual(
                key_ord,
                [("CAT_MULTI", 0), ("CAT_MULTI", 1), ("CAT_S1", 0), ("CAT_S2", 0)],
            )

            # Verify end-to-end main() exports label_noise_audit.csv preserving plant_individual_id and vision_is_ood
            out_csv = tables_dir / "label_noise_audit.csv"
            out_fig = tmp_path / "outputs" / "figures" / "GradCAM_audit_panel.png"
            with mock_patch("torch.hub.load", return_value=mod_xai._DemoViTFeatureExtractor()):
                with mock_patch("sys.argv", [
                    "05_cleanlab_vision_xai.py",
                    "--rosette-dir", str(data_dir / "cropped_patches"),
                    "--vouchers-csv", str(vouchers_csv),
                    "--manifest-csv", str(manifest_csv),
                    "--output-csv", str(out_csv),
                    "--output-figure", str(out_fig),
                    "--no-export-figures",
                ]):
                    mod_xai.main()

            exported = pd.read_csv(out_csv)
            self.assertEqual(len(exported), 4)
            self.assertEqual(
                list(zip(exported["catalogNumber"], exported["plant_individual_id"])),
                [("CAT_MULTI", 0), ("CAT_MULTI", 1), ("CAT_S1", 0), ("CAT_S2", 0)],
            )
            for req_col in ("plant_individual_id", "vision_ood_score", "vision_is_ood", "foreground_attention_ratio", "is_artifact_driven"):
                self.assertIn(req_col, exported.columns)
            self.assertEqual(exported["vision_is_ood"].dtype, np.dtype(bool))
            self.assertEqual(exported["is_artifact_driven"].dtype, np.dtype(bool))

            # Verify intraspecific-mode joins morphometrics_csv on ["catalogNumber", "plant_individual_id"]
            vouchers_intra = tables_dir / "vouchers_intra.csv"
            pd.DataFrame({
                "catalogNumber": ["CAT_MULTI", "CAT_S1", "CAT_S2"],
                "species_raw": ["Packera aurea", "Packera aurea", "Packera aurea"],
                "determiner_tier": ["Tier_1_Gold", "Tier_1_Gold", "Tier_1_Gold"],
                "institutionCode": ["NCU", "NY", "US"],
            }).to_csv(vouchers_intra, index=False)
            morph_intra = tables_dir / "morphometric_flags.csv"
            pd.DataFrame({
                "catalogNumber": ["CAT_MULTI", "CAT_MULTI", "CAT_S1", "CAT_S2"],
                "plant_individual_id": [0, 1, 0, 0],
                "gmm_cluster": [1, 2, 1, 2],
            }).to_csv(morph_intra, index=False)
            out_intra_csv = tables_dir / "label_noise_audit_intra.csv"
            with mock_patch("torch.hub.load", return_value=mod_xai._DemoViTFeatureExtractor()):
                with mock_patch("sys.argv", [
                    "05_cleanlab_vision_xai.py",
                    "--rosette-dir", str(data_dir / "cropped_patches"),
                    "--vouchers-csv", str(vouchers_intra),
                    "--manifest-csv", str(manifest_csv),
                    "--morphometrics-csv", str(morph_intra),
                    "--intraspecific-mode",
                    "--output-csv", str(out_intra_csv),
                    "--output-figure", str(out_fig),
                    "--no-export-figures",
                ]):
                    mod_xai.main()
            exp_intra = pd.read_csv(out_intra_csv)
            self.assertEqual(len(exp_intra), 4)
            lbl_by_key = dict(zip(zip(exp_intra["catalogNumber"], exp_intra["plant_individual_id"]), exp_intra["given_label"]))
            self.assertEqual(lbl_by_key[("CAT_MULTI", 0)], "Morphotype_K1")
            self.assertEqual(lbl_by_key[("CAT_MULTI", 1)], "Morphotype_K2")

    def test_vision_is_ood_and_artifact_boolean_thresholds_and_backward_compat(self):
        """Verify run_confident_learning_audit exports boolean vision_is_ood (>= 0.65) and preserves plant_individual_id."""
        pred_probs = np.array([
            [0.90, 0.10],
            [0.85, 0.15],
            [0.15, 0.85],
            [0.10, 0.90],
        ])
        labels = np.array([0, 0, 1, 1])

        # Multi-individual sheet CAT_SHEET1 (pid 0 and 1) + boundary tests at 0.65
        records_df = pd.DataFrame({
            "catalogNumber": ["CAT_SHEET1", "CAT_SHEET1", "CAT_SHEET2", "CAT_SHEET3"],
            "plant_individual_id": [0, 1, 0, 0],
            "species_raw": ["Packera aurea", "Packera aurea", "Packera obovata", "Packera obovata"],
            "vision_ood_score": [0.20, 0.6499, 0.65, 0.82],
            "foreground_attention_ratio": [0.92, 0.65, 0.6499, 0.42],
        })
        audit = run_confident_learning_audit(
            pred_probs=pred_probs,
            labels=labels,
            records_df=records_df,
            class_names=["Packera aurea", "Packera obovata"],
        )
        self.assertEqual(len(audit), 4)
        self.assertTrue(pd.api.types.is_integer_dtype(audit["plant_individual_id"]))
        self.assertEqual(
            list(zip(audit["catalogNumber"], audit["plant_individual_id"])),
            [("CAT_SHEET1", 0), ("CAT_SHEET1", 1), ("CAT_SHEET2", 0), ("CAT_SHEET3", 0)],
        )
        self.assertEqual(audit["vision_is_ood"].dtype, np.dtype(bool))
        self.assertEqual(audit["vision_is_ood"].tolist(), [False, False, True, True])
        self.assertEqual(audit["is_artifact_driven"].dtype, np.dtype(bool))
        self.assertEqual(audit["is_artifact_driven"].tolist(), [False, False, True, True])

        # Backward compatibility when plant_individual_id is omitted from records_df
        legacy_df = pd.DataFrame({
            "catalogNumber": ["LEG_1", "LEG_2"],
            "species_raw": ["Packera aurea", "Packera obovata"],
        })
        legacy_audit = run_confident_learning_audit(
            pred_probs=np.array([[0.9, 0.1], [0.1, 0.9]]),
            labels=np.array([0, 1]),
            records_df=legacy_df,
            class_names=["Packera aurea", "Packera obovata"],
        )
        self.assertIn("plant_individual_id", legacy_audit.columns)
        self.assertTrue(pd.api.types.is_integer_dtype(legacy_audit["plant_individual_id"]))
        self.assertEqual(legacy_audit["plant_individual_id"].tolist(), [0, 0])
        self.assertIn("vision_is_ood", legacy_audit.columns)
        self.assertEqual(legacy_audit["vision_is_ood"].dtype, np.dtype(bool))

    def test_promoted_bronze_anchor_alignment_and_conditional_abaxial_gate(self):
        """Verify promoted Bronze anchors (<8 Gold/Silver) enter OOF/OOD and abaxial filtering preserves structural reassignments without chromatic discordance."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            dense_dir = tmp_p / "cropped_patches" / "rosettes_dense"
            dense_dir.mkdir(parents=True)
            for fname in ("CAT_MULTI_p0_rosette.jpg", "CAT_MULTI_p1_rosette.jpg", "CAT_RARE_p0_rosette.jpg"):
                Image.new("RGB", (240, 240), color=(45, 95, 35)).save(dense_dir / fname)

            vouchers_csv = tmp_p / "curated_vouchers.csv"
            pd.DataFrame({
                "catalogNumber": ["CAT_MULTI", "CAT_RARE"],
                "species_raw": ["Packera aurea", "Packera obovata"],
                "determiner_tier": ["Tier_3_Bronze", "Tier_3_Bronze"],
                "institutionCode": ["NCU", "NY"],
            }).to_csv(vouchers_csv, index=False)

            morph_csv = tmp_p / "morphometric_flags.csv"
            pd.DataFrame({
                "catalogNumber": ["CAT_MULTI", "CAT_MULTI", "CAT_RARE"],
                "plant_individual_id": [0, 1, 0],
                "is_promoted_bronze": [True, False, False],
                "anchor_tier_source": ["Tier_3_Bronze_Promoted", "Tier_3_Bronze_Passive", "Tier_3_Bronze_Promoted"],
            }).to_csv(morph_csv, index=False)

            linked_df, _ = mod_xai.load_and_link_rosette_patches(
                rosette_dir=tmp_p / "cropped_patches",
                vouchers_csv=vouchers_csv,
                morphometrics_csv=morph_csv,
            )
            promo_by_key = dict(
                zip(zip(linked_df["catalogNumber"], linked_df["plant_individual_id"]), linked_df["is_promoted_bronze"])
            )
            self.assertEqual(
                promo_by_key,
                {("CAT_MULTI", 0): True, ("CAT_MULTI", 1): False, ("CAT_RARE", 0): True},
            )

        # Taxon 0 has 8 Gold anchors (sufficient -> Bronze NOT promoted); Taxon 1 has only 2 Gold + 4 Promoted Bronze (<8 GS -> promoted Bronze included at weight 1.0)
        rng = np.random.default_rng(2026)
        labels = np.array([0] * 10 + [1] * 8, dtype=int)
        tiers = (
            ["Tier_1_Gold"] * 8 + ["Tier_3_Bronze", "Tier_3_Bronze"]
            + ["Tier_1_Gold"] * 2 + ["Tier_3_Bronze"] * 6
        )
        promo_flags = np.array(
            [False] * 8 + [True, False]
            + [False, False, True, True, True, True, False, False],
            dtype=bool,
        )
        anc_mask, anc_weights = mod_xai._resolve_adaptive_anchors(labels, tiers, promo_flags)
        # Taxon 0 already has 8 Gold anchors -> index 8 (promoted flag True) stays passive
        self.assertFalse(bool(anc_mask[8]))
        # Taxon 1 has 2 Gold (<8) -> indices 12..15 (promoted Bronze) become active anchors with weight 1.0, while 16..17 remain passive
        self.assertEqual(anc_mask[10:18].tolist(), [True, True, True, True, True, True, False, False])
        np.testing.assert_allclose(anc_weights[12:16], np.ones(4), atol=1e-6)

        X = rng.normal(0.0, 0.10, size=(len(labels), 64)).astype(np.float32)
        X[labels == 0, :16] += 4.0
        X[labels == 1, 16:32] += 4.0
        # Corrupt the 2 unpromoted passive Bronze rows in Taxon 1 (indices 16, 17) to look like Taxon 0
        X[16:18, 16:32] = 0.0
        X[16:18, :16] += 4.0

        oof_probs, oof_acc, _ = compute_out_of_fold_probabilities(
            features=X,
            labels=labels,
            n_splits=3,
            determiner_tiers=tiers,
            is_promoted_bronze=promo_flags,
        )
        self.assertGreaterEqual(oof_acc, 0.90)
        self.assertEqual(int(np.argmax(oof_probs[16])), 0)
        self.assertEqual(int(np.argmax(oof_probs[17])), 0)

        ood_with_promo = compute_feature_space_ood_scores(
            features=X,
            labels=labels,
            pred_probs=oof_probs,
            determiner_tiers=tiers,
            is_promoted_bronze=promo_flags,
        )
        self.assertEqual(len(ood_with_promo), len(labels))
        self.assertTrue(np.all((ood_with_promo >= 0.0) & (ood_with_promo <= 1.0)))

        # Conditional abaxial gate test in run_confident_learning_audit
        abax_records = pd.DataFrame({
            "catalogNumber": ["ABAX_TOMENTUM", "ABAX_LOW_CONF", "ABAX_LOW_FG", "ABAX_TRUE_MISID"],
            "plant_individual_id": [0, 0, 0, 0],
            "species_raw": ["Packera aurea"] * 4,
            "surface_orientation": ["abaxial"] * 4,
            "tomentum_taxon_discordance": [True, False, False, False],
            "foreground_attention_ratio": [0.88, 0.88, 0.45, 0.86],
        })
        abax_probs = np.array([
            [0.04, 0.96],  # high c_error=0.96, conf=0.96, BUT tomentum_taxon_discordance=True -> suppressed
            [0.20, 0.80],  # conf=0.80 < 0.85, c_error=0.80 <= 0.88 -> suppressed
            [0.05, 0.95],  # foreground_attention_ratio=0.45 < 0.65 -> suppressed
            [0.05, 0.95],  # tomentum_taxon_discordance=False, conf=0.95 >= 0.85, c_error=0.95 > 0.88, fg=0.86 >= 0.65 -> PRESERVED!
        ])
        abax_audit = run_confident_learning_audit(
            pred_probs=abax_probs,
            labels=np.array([0, 0, 0, 0]),
            records_df=abax_records,
            class_names=["Packera aurea", "Packera obovata"],
            error_threshold=0.85,
            filter_abaxial_artifacts=True,
        ).set_index("catalogNumber")

        self.assertFalse(bool(abax_audit.loc["ABAX_TOMENTUM", "is_label_corrupted"]))
        self.assertFalse(bool(abax_audit.loc["ABAX_TOMENTUM", "is_label_issue"]))
        self.assertFalse(bool(abax_audit.loc["ABAX_LOW_CONF", "is_label_corrupted"]))
        self.assertFalse(bool(abax_audit.loc["ABAX_LOW_FG", "is_label_corrupted"]))
        self.assertTrue(bool(abax_audit.loc["ABAX_TRUE_MISID", "is_label_corrupted"]))
        self.assertTrue(bool(abax_audit.loc["ABAX_TRUE_MISID", "is_label_issue"]))

    def test_basal_manifest_filter_dissection_neutral_weights_and_narrow_pinna_pooling(self):
        """Verify basal manifest filtering + lowest-lamina fallback, dissection-neutral leaf weights, and narrow-pinna ViT token pooling."""
        from scripts.analysis.vision_dissection_utils import (
            compute_dissection_neutral_leaf_weight,
            filter_basal_manifest_rows_py,
        )

        # 1. Dissection-neutral quality weighting: dissected and entire leaves with equal UCS get identical weights
        w_diss, sol_diss = compute_dissection_neutral_leaf_weight({"solidity": 0.36, "ucs_score": 0.92})
        w_ent, sol_ent = compute_dissection_neutral_leaf_weight({"solidity": 0.94, "ucs_score": 0.92})
        self.assertAlmostEqual(w_diss, w_ent, places=6)
        self.assertAlmostEqual(sol_diss, 0.36, places=4)
        self.assertAlmostEqual(sol_ent, 0.94, places=4)

        # 2. Positional basal gate + detached-rosette fallback + disk-glob cauline leak prevention
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            rosette_dir = tmp_path / "cropped_patches"
            leaves_dir = rosette_dir / "basal_leaves_raw"
            masks_dir = tmp_path / "masks" / "tier1_intact"
            leaves_dir.mkdir(parents=True)
            masks_dir.mkdir(parents=True)

            # Create leaf crops on disk: CAT_BASAL has leaf1 (cauline, cy=0.15) and leaf2 (basal, cy=0.78)
            # CAT_FALLBACK has leaf1 (cy=0.24) and leaf2 (cy=0.32) -> rescued by max - 0.15 fallback
            for cat, lid in [("CAT_BASAL", 1), ("CAT_BASAL", 2), ("CAT_FALLBACK", 1), ("CAT_FALLBACK", 2)]:
                img_p = leaves_dir / f"{cat}_p0_leaf{lid}.png"
                msk_p = masks_dir / f"{cat}_p0_leaf{lid}.png"
                Image.new("RGB", (224, 224), color=(60, 110, 50)).save(img_p)
                Image.new("L", (224, 224), color=255).save(msk_p)

            manifest_df = pd.DataFrame([
                {
                    "catalogNumber": "CAT_BASAL", "plant_individual_id": 0, "leaf_id": 1,
                    "assigned_tier": "Tier_1_Direct", "solidity": 0.92, "circularity": 0.65,
                    "ucs_score": 0.90, "cy_rel_plant_bbox": 0.15, "blade_area_cm2": 1.8,
                    "effective_min_blade_area_cm2": 0.4,
                    "leaf_crop_path": str(leaves_dir / "CAT_BASAL_p0_leaf1.png"),
                    "mask_path": str(masks_dir / "CAT_BASAL_p0_leaf1.png"),
                },
                {
                    "catalogNumber": "CAT_BASAL", "plant_individual_id": 0, "leaf_id": 2,
                    "assigned_tier": "Tier_1_Direct", "solidity": 0.42, "circularity": 0.14,
                    "ucs_score": 0.91, "cy_rel_plant_bbox": 0.78, "blade_area_cm2": 2.1,
                    "effective_min_blade_area_cm2": 0.4,
                    "leaf_crop_path": str(leaves_dir / "CAT_BASAL_p0_leaf2.png"),
                    "mask_path": str(masks_dir / "CAT_BASAL_p0_leaf2.png"),
                },
                {
                    "catalogNumber": "CAT_FALLBACK", "plant_individual_id": 0, "leaf_id": 1,
                    "assigned_tier": "Tier_1_Direct", "solidity": 0.48, "circularity": 0.16,
                    "ucs_score": 0.89, "cy_rel_plant_bbox": 0.24, "blade_area_cm2": 1.6,
                    "effective_min_blade_area_cm2": 0.4,
                    "leaf_crop_path": str(leaves_dir / "CAT_FALLBACK_p0_leaf1.png"),
                    "mask_path": str(masks_dir / "CAT_FALLBACK_p0_leaf1.png"),
                },
                {
                    "catalogNumber": "CAT_FALLBACK", "plant_individual_id": 0, "leaf_id": 2,
                    "assigned_tier": "Tier_1_Direct", "solidity": 0.88, "circularity": 0.58,
                    "ucs_score": 0.92, "cy_rel_plant_bbox": 0.32, "blade_area_cm2": 1.9,
                    "effective_min_blade_area_cm2": 0.4,
                    "leaf_crop_path": str(leaves_dir / "CAT_FALLBACK_p0_leaf2.png"),
                    "mask_path": str(masks_dir / "CAT_FALLBACK_p0_leaf2.png"),
                },
            ])
            filt_df = filter_basal_manifest_rows_py(manifest_df, tier_col="assigned_tier")
            self.assertEqual(len(filt_df), 3)
            self.assertEqual(
                set(zip(filt_df["catalogNumber"].tolist(), filt_df["leaf_id"].astype(int).tolist())),
                {("CAT_BASAL", 2), ("CAT_FALLBACK", 1), ("CAT_FALLBACK", 2)},
            )

            manifest_csv = tmp_path / "extracted_leaf_manifest.csv"
            vouchers_csv = tmp_path / "curated_vouchers.csv"
            manifest_df.to_csv(manifest_csv, index=False)
            pd.DataFrame({
                "catalogNumber": ["CAT_BASAL", "CAT_FALLBACK"],
                "species_raw": [TEST_SPECIES[0], TEST_SPECIES[1]],
                "determiner_tier": ["Tier_1_Gold", "Tier_1_Gold"],
                "institutionCode": ["NCU", "US"],
            }).to_csv(vouchers_csv, index=False)

            linked_df, _ = mod_xai.load_and_link_rosette_patches(
                rosette_dir=rosette_dir,
                vouchers_csv=vouchers_csv,
                demo=False,
                manifest_csv=manifest_csv,
            )
            basal_row = linked_df.loc[linked_df["catalogNumber"] == "CAT_BASAL"].iloc[0]
            # Excluded cauline leaf1 on disk must NOT be re-ingested by disk glob fallback
            self.assertEqual(int(basal_row["n_rosette_patches"]), 1)
            self.assertEqual(len(basal_row["tier1_leaf_patches"]), 1)
            self.assertIn("leaf2.png", basal_row["tier1_leaf_patches"][0])
            fb_row = linked_df.loc[linked_df["catalogNumber"] == "CAT_FALLBACK"].iloc[0]
            self.assertEqual(int(fb_row["n_rosette_patches"]), 2)
            self.assertEqual(len(fb_row["tier1_leaf_patches"]), 2)

            # 3. Continuous area-weighted ViT token pooling for narrow dissected pinnae (0.05 <= patch_frac < 0.15)
            import torch

            class PinnaDINOv2(torch.nn.Module):
                def forward_features(self, x: torch.Tensor):
                    b = x.shape[0]
                    pt = torch.zeros((b, 256, 768), dtype=torch.float32)
                    # Patch (0, 0) has +10 in dim 0; all background patches have +10 in dim 1
                    pt[:, 0, 0] = 10.0
                    pt[:, 1:, 1] = 10.0
                    return {"x_norm_patchtokens": pt, "x_norm_clstoken": torch.zeros((b, 768))}

            narrow_img = tmp_path / "narrow_pinna.png"
            narrow_msk = tmp_path / "narrow_pinna_mask.png"
            Image.new("RGB", (224, 224), color=(80, 120, 70)).save(narrow_img)
            # Put 20 foreground pixels inside top-left 14x14 patch -> patch_frac = 20/196 = 0.102 (< 0.15 old threshold, >= 0.05 new)
            m_arr = np.zeros((224, 224), dtype=np.uint8)
            m_arr[:2, :10] = 255
            Image.fromarray(m_arr, mode="L").save(narrow_msk)

            rec = [{
                "catalogNumber": "CAT_PINNA",
                "plant_individual_id": 0,
                "label_idx": 0,
                "patch_path": str(narrow_img),
                "mask_path": str(narrow_msk),
                "tier1_leaf_patches": [str(narrow_img)],
                "tier1_leaf_masks": [str(narrow_msk)],
                "tier1_leaf_weights": [1.0],
            }]
            with patch("torch.hub.load", return_value=PinnaDINOv2()):
                emb, _, _ = extract_dinov2_embeddings(rec, demo=False, db_path=tmp_path / "ephem.sqlite")
            # Foreground patch (0,0) must dominate dim 0 rather than falling back to background global average (dim 1)
            self.assertGreater(float(emb[0, 0]), 0.95)
            self.assertAlmostEqual(float(emb[0, 1]), 0.0, places=5)


if __name__ == "__main__":
    unittest.main()
