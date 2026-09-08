import os
import shutil
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

# Assume PROJECT_ROOT is up two levels
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

class TestPhase2Phase3Handshake(unittest.TestCase):
    def setUp(self):
        self.test_dir = PROJECT_ROOT / "data" / "_test_handshake"
        self.test_dir.mkdir(parents=True, exist_ok=True)
        
        self.vouchers_path = self.test_dir / "curated_vouchers.csv"
        df = pd.DataFrame([
            {"catalogNumber": "TEST001", "species_raw": "Packera dubia"}
        ])
        df.to_csv(self.vouchers_path, index=False)

    def tearDown(self):
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)

    @patch("shutil.which", return_value="/usr/bin/Rscript")
    @patch("main.verify_dir_has_files")
    @patch("main.verify_file_exists")
    @patch("subprocess.run")
    def test_end_to_end_handshake(self, mock_run, mock_verify_file, mock_verify_dir, mock_which):
        # This test ensures that main.py parses arguments and calls the
        # underlying scripts with the correct default paths for the handshake.
        from main import build_parser, run_segment, run_morphometrics
        
        parser = build_parser()
        
        # Test Phase 2
        args_segment = parser.parse_args([
            "segment", 
            "--vouchers", str(self.vouchers_path),
            "--output-dir", str(self.test_dir),
            "--weights", str(self.test_dir / "dummy.pth")
        ])
        
        # We just need to check if run_segment correctly resolves paths.
        # But run_segment does subprocess.run which we mocked.
        mock_run.return_value.returncode = 0
        
        try:
            run_segment(args_segment, cfg={"paths": {"workspace_root": str(self.test_dir), "contours_dir": str(self.test_dir / "contours")}, "segmentation": {"device": "cpu", "min_solidity": 0.7, "min_ucs": 0.7, "score_thresh": 0.5}})
        except SystemExit:
            pass # ignore sys.exit(1) due to file not founds inside verify checks
            
        # Check that the command passed to subprocess.run for segment has the right paths
        segment_cmd = mock_run.call_args_list[0][0][0]
        self.assertIn(str(self.test_dir), segment_cmd)

        # Test Phase 3
        args_morph = parser.parse_args([
            "morphometrics",
            "--vouchers", str(self.vouchers_path),
            "--input", str(self.test_dir / "contours"),
            "--manifest", str(self.test_dir / "tables" / "extracted_leaf_manifest.csv"),
            "--harmonics-out", str(self.test_dir / "leaf_efa_harmonics.csv"),
            "--flags-out", str(self.test_dir / "flags.csv"),
            "--plot-out", str(self.test_dir / "plot.pdf"),
            "--report-out", str(self.test_dir / "report.csv"),
        ])
        
        try:
            run_morphometrics(args_morph, cfg={"paths": {"contours_dir": str(self.test_dir / "contours")}, "morphometrics": {"harmonics": 12, "num_pcs": 5, "max_k": 5}})
        except SystemExit:
            pass # ignore verify_dir_has_files checks
            
        morph_cmd = mock_run.call_args_list[1][0][0]
        # Verify the --input and --manifest flags were correctly forwarded to the R script
        self.assertIn("--input", morph_cmd)
        self.assertIn(str(self.test_dir / "contours"), morph_cmd)
        self.assertIn("--manifest", morph_cmd)
        self.assertIn(str(self.test_dir / "tables" / "extracted_leaf_manifest.csv"), morph_cmd)

if __name__ == "__main__":
    unittest.main()
