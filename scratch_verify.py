import cv2
import numpy as np
import os

def test_texture_variance():
    image_path = "data/raw_vouchers/000331814.jpg"
    print(f"Loading image {image_path}...")
    
    img = cv2.imread(image_path)
    if img is None:
        print("Failed to load image.")
        return
        
    img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    h, w = img.shape[:2]
    print(f"Image shape: {w}x{h}")
    
    tile_size = 1024
    
    # Save patches to see them
    os.makedirs('scratch', exist_ok=True)
    
    regions = {
        "Top-Left (Pure Paper)": (0, 0),
        "Top-Right (Pure Paper/Edge)": (0, w - tile_size),
        "Center (Contains Plant)": (h // 2 - tile_size // 2, w // 2 - tile_size // 2),
        "Bottom-Right (Label/Text likely)": (h - tile_size, w - tile_size),
        "Mid-Left (Possibly just paper or small stem)": (h // 2, 0)
    }
    
    print("\nTesting Variance Metrics:")
    print(f"{'Region':<40} | {'STD (Gray)':>10} | {'Laplacian Var':>13} | {'Canny Edges':>12}")
    print("-" * 85)
    
    for name, (y, x) in regions.items():
        y = max(0, min(y, h - tile_size))
        x = max(0, min(x, w - tile_size))
        
        tile_gray = img_gray[y:y+tile_size, x:x+tile_size]
        
        std_gray = np.std(tile_gray)
        laplacian_var = cv2.Laplacian(tile_gray, cv2.CV_64F).var()
        edges = cv2.Canny(tile_gray, 50, 150)
        edge_density = np.sum(edges > 0) / (tile_size * tile_size) * 100
        
        print(f"{name:<40} | {std_gray:>10.2f} | {laplacian_var:>13.2f} | {edge_density:>11.3f}%")

if __name__ == '__main__':
    test_texture_variance()
