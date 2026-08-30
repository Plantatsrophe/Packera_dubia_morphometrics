import cv2
import numpy as np

# Load test image
image_bgr = cv2.imread("scratch/qc_test_LSU00220854.jpg")
h, w = image_bgr.shape[:2]

lab = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2LAB)
l_chan = lab[:, :, 0]
_, plant_thresh = cv2.threshold(l_chan, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
clean_plant = cv2.morphologyEx(plant_thresh, cv2.MORPH_OPEN, kernel)

# To separate roots (thin) from rosettes (thick), we can use morphological opening with a large kernel
# A 25x25 kernel will erase anything thinner than 25 pixels (most roots)
large_kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (25, 25))
rosettes_only = cv2.morphologyEx(clean_plant, cv2.MORPH_OPEN, large_kernel)

# The roots are the difference between the original plant mask and the rosettes_only mask!
# (Well, roots plus peduncles and thin leaves)
thin_structures = cv2.subtract(clean_plant, rosettes_only)

cv2.imwrite("scratch/debug_rosettes.jpg", rosettes_only)
cv2.imwrite("scratch/debug_thin.jpg", thin_structures)
