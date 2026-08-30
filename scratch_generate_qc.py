import cv2
import os
from scripts.core.artifact_harvester import ArtifactHarvester
from scripts.core.botanical_annotations import extract_botanical_annotations
from scripts.core.config import CLASS_COLORS_BGR, CLASS_NAMES

def generate_qc(img_name, out_path):
    img_path_in = f"/home/brandon/Packera_dubia_morphometrics/data/raw_vouchers/{img_name}.jpg"
    image_bgr = cv2.imread(img_path_in)
    
    harvester = ArtifactHarvester()
    artifact_anns = harvester.detect_and_extract_sheet_artifacts(image_bgr, img_name)
    bot_anns = extract_botanical_annotations(image_bgr, artifact_anns)
    
    vis = image_bgr.copy()
    for ann in artifact_anns + bot_anns:
        pts = ann.polygon.astype(int)
        cname = CLASS_NAMES[ann.class_id]
        color = CLASS_COLORS_BGR[ann.class_id]
        cv2.polylines(vis, [pts], True, color, 12)
        bx, by, _, _ = ann.bbox
        cv2.putText(vis, cname, (int(bx), int(by) - 25), cv2.FONT_HERSHEY_SIMPLEX, 3.5, color, 8)
        
    vis_resized = cv2.resize(vis, (0,0), fx=0.15, fy=0.15)
    cv2.imwrite(out_path, vis_resized)
    print(f"Saved {out_path}")

generate_qc("NCU00438248", "/home/brandon/.gemini/antigravity-ide/brain/b5c908d7-1acb-4bcd-8f2c-7ed5ffd87a3c/scratch/qc_test_NCU00438248_v4.jpg")
generate_qc("LSU00220854", "/home/brandon/.gemini/antigravity-ide/brain/b5c908d7-1acb-4bcd-8f2c-7ed5ffd87a3c/scratch/qc_test_LSU00220854_v4.jpg")
