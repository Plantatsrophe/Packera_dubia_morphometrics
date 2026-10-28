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
import tempfile
import time
import unittest
from pathlib import Path
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
TARGET_TAXA = mod_xai.TARGET_TAXA


class TestCleanlabVisionXAI(unittest.TestCase):
    """Test suite for DINOv2 macro-architecture pooling & Cleanlab XAI."""

    @classmethod
    def setUpClass(cls):
        cls.repo_root = Path(__file__).resolve().parents[2]
        cls.audit_csv = cls.repo_root / "data" / "tables" / "label_noise_audit.csv"
        cls.figure_png = cls.repo_root / "outputs" / "figures" / "GradCAM_audit_panel.png"

    def test_taxonomic_standardization(self):
        """Test species synonymy standardizer for the 4 core complex taxa."""
        self.assertEqual(standardize_packera_taxon("Senecio smallii Britton"), "Packera anonyma")
        self.assertEqual(standardize_packera_taxon("Packera anonyma (Alph.Wood) W.A.Weber"), "Packera anonyma")
        self.assertEqual(standardize_packera_taxon("Senecio tomentosus Michx."), "Packera dubia")
        self.assertEqual(standardize_packera_taxon("Packera tomentosa C.Jeffrey"), "Packera dubia")
        self.assertEqual(standardize_packera_taxon("Packera dubia (Spreng.) Trock & Mabb."), "Packera dubia")
        self.assertEqual(standardize_packera_taxon("Packera paupercula var. savannarum"), "Packera paupercula")
        self.assertEqual(standardize_packera_taxon("Packera paupercula var. pseudotomentosa"), "Packera paupercula")
        self.assertEqual(standardize_packera_taxon("Senecio plattensis Nutt."), "Packera plattensis")
        self.assertEqual(standardize_packera_taxon(None), "Unknown")

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
                "label_idx": i % 4,
                "taxon": TARGET_TAXA[i % 4],
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
            "species_raw": ["Packera anonyma", "Packera dubia", "Packera paupercula", "Packera plattensis"],
            "label_idx": [0, 1, 2, 3],
            "determiner_tier": ["Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze", "Tier_1_Gold"],
            "foreground_attention_ratio": [0.88, 0.91, 0.79, 0.52],
        })

        audit = run_confident_learning_audit(
            pred_probs=pred_probs,
            labels=labels,
            records_df=records_df,
            class_names=TARGET_TAXA,
            error_threshold=0.85,
            features=feats,
        )

        required_step07_cols = [
            "catalogNumber",
            "c_error",
            "predicted_label",
            "vision_predicted_label",
            "vision_confidence",
            "vision_second_prob",
            "vision_margin",
            "vision_ood_score",
            "foreground_attention_ratio",
            "is_label_issue",
            "is_label_corrupted",
            "tomentum_taxon_discordance",
        ]
        for col in required_step07_cols:
            self.assertIn(col, audit.columns, f"Missing required Step 07 contract column: {col}")

        self.assertTrue(audit.loc[0, "is_label_corrupted"])
        self.assertFalse(audit.loc[1, "is_label_corrupted"])
        self.assertAlmostEqual(float(audit.loc[0, "vision_margin"]), 0.92 - 0.05, places=4)
        self.assertTrue(bool(audit.loc[3, "is_artifact_driven"]))
        self.assertIn("Downweight", audit.loc[3, "triage_action"])

    def test_rosettes_dense_ingestion_multi_crop_aggregation_and_dynamic_indexing(self):
        """Verify rosettes_dense ingestion, multi-crop sheet aggregation, and dynamic 0..K-1 taxon indexing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            rosette_dir = tmp_path / "cropped_patches"
            dense_dir = rosette_dir / "rosettes_dense"
            dense_dir.mkdir(parents=True)

            # Create multiple crops for CATANON and single crops for CATDUB, CATPAUP, CATPLAT
            for fname in (
                "CATANON_p0_rosette.jpg",
                "CATANON_p1_rosette.jpg",
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


if __name__ == "__main__":
    unittest.main()
