import cv2
import numpy as np

def test_label(img_path):
    image_bgr = cv2.imread(img_path)
    h, w = image_bgr.shape[:2]
    gray = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2GRAY)
    
    # Bottom half
    label_roi = gray[int(h*0.5):, :]
    
    # 1. Canny Edges
    edges = cv2.Canny(label_roi, 50, 150)
    
    # 2. Morphological close to join text into a solid block
    # Text is usually grouped, so a 35x35 kernel works well to bridge text lines
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (35, 35))
    closed = cv2.morphologyEx(edges, cv2.MORPH_CLOSE, kernel)
    
    # 3. Find contours
    contours, _ = cv2.findContours(closed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    output = np.zeros_like(label_roi)
    for cnt in contours:
        area = cv2.contourArea(cnt)
        if area > (w * h * 0.005): # At least 0.5% of image
            bx, by, bw, bh = cv2.boundingRect(cnt)
            # Check bounding box density instead of contour area to allow for irregular edge borders
            bbox_area = bw * bh
            if bbox_area > (w * h * 0.015):
                rect_ratio = area / bbox_area
                aspect = bw / max(1.0, bh)
                if 0.5 < aspect < 3.0 and rect_ratio > 0.4:
                    cv2.rectangle(output, (bx, by), (bx+bw, by+bh), 255, -1)
                    
    cv2.imwrite("scratch/debug_label_edges.jpg", edges)
    cv2.imwrite("scratch/debug_label_closed.jpg", closed)
    cv2.imwrite("scratch/debug_label_detected.jpg", output)

test_label("scratch/qc_test_LSU00220854_v3.jpg")
