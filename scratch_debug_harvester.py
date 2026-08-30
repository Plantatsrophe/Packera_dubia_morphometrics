import cv2
import sys
import numpy as np
from pathlib import Path

# Add project root to path
sys.path.insert(0, "/home/brandon/Packera_dubia_morphometrics")

image_path = "/home/brandon/Packera_dubia_morphometrics/data/raw_vouchers/LSU00220854.jpg"
image_bgr = cv2.imread(image_path)
if image_bgr is None:
    print("Could not read image")
    sys.exit(1)

h, w = image_bgr.shape[:2]
gray = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2GRAY)

print(f"Image dimensions: {w}x{h}")

# 1. Herbarium Label
label_y1, label_y2 = int(h * 0.60), int(h * 0.99)
label_x1, label_x2 = 0, w
label_roi = gray[label_y1:label_y2, label_x1:label_x2]
blurred = cv2.GaussianBlur(label_roi, (25, 25), 0)
_, label_thresh = cv2.threshold(blurred, 190, 255, cv2.THRESH_BINARY)
contours, _ = cv2.findContours(label_thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

print(f"Found {len(contours)} contours for herbarium_label (Thresholding)")
for cnt in contours:
    area = cv2.contourArea(cnt)
    if area > (w * h * 0.015):
        bx, by, bw, bh = cv2.boundingRect(cnt)
        rect_ratio = area / max(1.0, float(bw * bh))
        aspect = bw / max(1, bh)
        print(f"  - Large contour: area={area:.0f}, rect_ratio={rect_ratio:.2f}, aspect={aspect:.2f}, bw={bw}, bh={bh}")

# 2. Color Calibration Chart
# Look at both left and right margins
left_margin_w = int(w * 0.25)
right_margin_x = int(w * 0.75)
sat = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2HSV)[:, :, 1]
_, sat_thresh = cv2.threshold(sat, 80, 255, cv2.THRESH_BINARY)

# Mask out the middle to only look at margins
mask = np.zeros_like(sat_thresh)
mask[:, :left_margin_w] = 255
mask[:, right_margin_x:] = 255
sat_thresh_margins = cv2.bitwise_and(sat_thresh, sat_thresh, mask=mask)

kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (45, 45))
closed_sat = cv2.morphologyEx(sat_thresh_margins, cv2.MORPH_CLOSE, kernel)
cv2.imwrite("scratch/debug_color_chart_mask.jpg", closed_sat)
cv2.imwrite("scratch/debug_label_thresh.jpg", label_thresh)

chart_cnts, _ = cv2.findContours(closed_sat, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

print(f"\nFound {len(chart_cnts)} contours for color_chart (global sat in margins)")
for cnt in chart_cnts:
    area = cv2.contourArea(cnt)
    if (w * h * 0.003) < area < (w * h * 0.08):
        bx, by, bw, bh = cv2.boundingRect(cnt)
        rect_ratio = area / max(1.0, float(bw * bh))
        print(f"  - Chart contour: area={area:.0f}, rect_ratio={rect_ratio:.2f}, bx={bx}, by={by}, bw={bw}, bh={bh}")
