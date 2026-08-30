import cv2
import numpy as np

def extract_botanical_masks(image_path):
    image_bgr = cv2.imread(image_path)
    h, w = image_bgr.shape[:2]

    lab = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2LAB)
    l_chan = lab[:, :, 0]
    _, plant_thresh = cv2.threshold(l_chan, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    
    kernel_small = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    clean_plant = cv2.morphologyEx(plant_thresh, cv2.MORPH_OPEN, kernel_small)
    
    # Extract thick biomass (Rosettes, wide leaves, capitulums)
    kernel_large = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (25, 25))
    thick_plant = cv2.morphologyEx(clean_plant, cv2.MORPH_OPEN, kernel_large)
    
    # Extract thin biomass (Roots, peduncles, petioles)
    thin_plant = cv2.subtract(clean_plant, thick_plant)
    
    # Clean up the thin mask to remove 1px edge halos
    thin_plant = cv2.morphologyEx(thin_plant, cv2.MORPH_OPEN, kernel_small)
    
    cv2.imwrite("scratch/debug_thick.jpg", thick_plant)
    cv2.imwrite("scratch/debug_thin.jpg", thin_plant)

extract_botanical_masks("scratch/qc_test_LSU00220854_v3.jpg")
