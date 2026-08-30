import cv2
import numpy as np

def analyze_thin(img_name):
    print(f"\n--- Analyzing {img_name} ---")
    img_path = f"/home/brandon/Packera_dubia_morphometrics/data/raw_vouchers/{img_name}.jpg"
    image_bgr = cv2.imread(img_path)
    h, w = image_bgr.shape[:2]
    
    lab = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2LAB)
    l_chan = lab[:, :, 0]
    _, plant_thresh = cv2.threshold(l_chan, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    
    kernel_small = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    clean_plant = cv2.morphologyEx(plant_thresh, cv2.MORPH_OPEN, kernel_small)
    
    kernel_large = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (25, 25))
    thick_plant = cv2.morphologyEx(clean_plant, cv2.MORPH_OPEN, kernel_large)
    thin_plant = cv2.subtract(clean_plant, thick_plant)
    thin_plant = cv2.morphologyEx(thin_plant, cv2.MORPH_OPEN, kernel_small)
    
    contours, _ = cv2.findContours(thin_plant, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    for cnt in contours:
        area = cv2.contourArea(cnt)
        if area > (w * h * 0.0001):
            hull = cv2.convexHull(cnt)
            solidity = area / max(1.0, cv2.contourArea(hull))
            bx, by, bw, bh = cv2.boundingRect(cnt)
            aspect = bw / max(1.0, bh)
            perimeter = cv2.arcLength(cnt, True)
            compactness = (perimeter ** 2) / area if area > 0 else 0
            
            print(f"[{by}:{by+bh}, {bx}:{bx+bw}] Area={area:.0f} | Aspect={aspect:.2f} | Solidity={solidity:.2f} | Compactness={compactness:.0f}")

analyze_thin("NCU00438248")
analyze_thin("LSU00220854")
