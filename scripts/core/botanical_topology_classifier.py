"""
Botanical Organ Topology Classifier Module
==========================================
Differentiates narrow elongated botanical organs (root_rhizome vs. cauline_stem vs. leaf_petiole)
using medial axis skeletonization, 3x3 neighbor-counting convolution kernels, and vertical
herbarium sheet priors.

Ontological Classification:
- leaf_petiole: narrow stalk attaching leaf blade to caudex/stem, or organ directly connected to a laminar blade.
- cauline_stem: continuous vertical flowering stalk / scape (elevated norm_y, low branch density <= 2).
- root_rhizome: subterranean anchoring and storage organs (basal norm_y > 0.65, highly branched >= 3).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np
from skimage.morphology import skeletonize


# 3x3 Neighbor-counting convolution kernel
# Center weight is 10, each 8-connected neighbor weight is 1.
# - Center pixel present: convolution base is 10.
# - Number of 8-connected neighbors = convolution_value - 10.
#   * Isolated point: 10 (0 neighbors)
#   * Endpoint: 11 (1 neighbor)
#   * Regular segment point: 12 (2 neighbors)
#   * Branch / Junction point: >= 13 (>= 3 neighbors)
NEIGHBOR_KERNEL: np.ndarray = np.array([
    [1,  1, 1],
    [1, 10, 1],
    [1,  1, 1]
], dtype=np.int32)

# Canonical Botanical Class Constants
ROOT_RHIZOME = "root_rhizome"
CAULINE_STEM = "cauline_stem"
LEAF_PETIOLE = "leaf_petiole"

# Default Spatial Prior Thresholds
DEFAULT_ROOT_NORM_Y_THRESHOLD: float = 0.65
DEFAULT_ROOT_MIN_BRANCH_POINTS: int = 3
DEFAULT_STEM_NORM_Y_THRESHOLD: float = 0.70
DEFAULT_STEM_MAX_BRANCH_POINTS: int = 2
DEFAULT_FALLBACK_PETIOLE_NORM_Y: float = 0.50


@dataclass
class TopologicalMetrics:
    """
    Morphometric topological summary parameters extracted from a candidate organ skeleton.
    """
    predicted_class: str
    num_branch_points: int
    num_end_points: int
    num_junction_clusters: int
    skeleton_pixel_count: int
    path_length_px: float
    chord_length_px: float
    tortuosity: float
    norm_y: float
    y_centroid_px: float
    x_centroid_px: float
    mask_area_px: int
    has_connected_blade: bool
    pixel_size_mm: Optional[float] = None
    path_length_mm: Optional[float] = None
    chord_length_mm: Optional[float] = None
    mask_area_mm2: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert metrics to a JSON/table serializable dictionary."""
        return asdict(self)


def extract_medial_axis_skeleton(mask: np.ndarray) -> np.ndarray:
    """
    Performs medial axis skeletonization on a 2D binary organ mask.

    Args:
        mask: 2D numpy array (binary, boolean, or uint8 mask).

    Returns:
        2D uint8 numpy array where 1 represents skeleton pixels and 0 is background.
    """
    if mask is None or mask.size == 0 or not np.any(mask > 0):
        if mask is not None and mask.ndim == 2:
            return np.zeros_like(mask, dtype=np.uint8)
        return np.zeros((0, 0), dtype=np.uint8)

    binary_mask = (mask > 0).astype(bool)
    skel = skeletonize(binary_mask)
    return skel.astype(np.uint8)


def analyze_skeleton_topology(
    skeleton: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, int, int]:
    """
    Convolves the 3x3 neighbor-counting kernel over the skeleton to identify branch/junction
    points (conv >= 13) and end points (conv == 11).

    Args:
        skeleton: 2D uint8 array with binary skeleton (1 = foreground, 0 = background).

    Returns:
        Tuple of:
            - conv: 2D int32 array containing convolution values.
            - branch_mask: 2D boolean mask of branch points (conv >= 13).
            - end_mask: 2D boolean mask of end points (conv == 11).
            - num_branch_points: Total count of branch point pixels.
            - num_end_points: Total count of end point pixels.
            - num_junction_clusters: Count of connected component clusters for branch junctions.
    """
    if skeleton is None or skeleton.size == 0 or not np.any(skeleton > 0):
        h, w = skeleton.shape if (skeleton is not None and skeleton.ndim == 2) else (0, 0)
        conv = np.zeros((h, w), dtype=np.int32)
        branch_mask = np.zeros((h, w), dtype=bool)
        end_mask = np.zeros((h, w), dtype=bool)
        return conv, branch_mask, end_mask, 0, 0, 0

    skel_f32 = (skeleton > 0).astype(np.float32)
    kernel_f32 = NEIGHBOR_KERNEL.astype(np.float32)
    conv_f32 = cv2.filter2D(skel_f32, cv2.CV_32F, kernel_f32, borderType=cv2.BORDER_CONSTANT)
    conv = np.round(conv_f32).astype(np.int32)

    branch_mask = (skel_f32 > 0) & (conv >= 13)
    end_mask = (skel_f32 > 0) & (conv == 11)

    num_branch_points = int(np.count_nonzero(branch_mask))
    num_end_points = int(np.count_nonzero(end_mask))

    # Calculate connected junction clusters (resolves multi-pixel junction nodes)
    if num_branch_points > 0:
        num_labels, _ = cv2.connectedComponents(branch_mask.astype(np.uint8), connectivity=8)
        num_junction_clusters = int(num_labels - 1)
    else:
        num_junction_clusters = 0

    return conv, branch_mask, end_mask, num_branch_points, num_end_points, num_junction_clusters


def calculate_mask_centroid(mask: np.ndarray) -> Tuple[float, float]:
    """
    Computes (x_centroid, y_centroid) for a 2D mask in pixel coordinates.
    """
    if mask is None or mask.size == 0 or not np.any(mask > 0):
        return 0.0, 0.0

    y_indices, x_indices = np.where(mask > 0)
    x_c = float(np.mean(x_indices))
    y_c = float(np.mean(y_indices))
    return x_c, y_c


def calculate_skeleton_path_and_chord_length(
    skeleton: np.ndarray,
    end_mask: np.ndarray
) -> Tuple[float, float, float]:
    """
    Computes the metric Euclidean path length, maximum endpoint chord length, and tortuosity.

    Path Length:
        Sum of orthogonal steps (1.0) and diagonal steps (sqrt(2)) along the skeleton graph edges.
    Chord Length:
        Maximum Euclidean distance between any pair of skeleton endpoints.
    Tortuosity:
        Ratio of path length to chord length (tau = L / C >= 1.0).

    Returns:
        Tuple of (path_length_px, chord_length_px, tortuosity).
    """
    skel_pts = np.argwhere(skeleton > 0)
    num_skel_pts = len(skel_pts)

    if num_skel_pts == 0:
        return 0.0, 0.0, 1.0
    if num_skel_pts == 1:
        return 1.0, 1.0, 1.0

    # Metric path length using 8-connected neighbor step weights
    # Orthogonal step weight = 1.0, Diagonal step weight = sqrt(2)
    # Each undirected edge is visited twice, so divide total by 2.
    ortho_kernel = np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ], dtype=np.uint8)

    diag_kernel = np.array([
        [1, 0, 1],
        [0, 0, 0],
        [1, 0, 1]
    ], dtype=np.uint8)

    skel_u8 = (skeleton > 0).astype(np.uint8)
    ortho_neighbors = cv2.filter2D(skel_u8, cv2.CV_32F, ortho_kernel, borderType=cv2.BORDER_CONSTANT)
    diag_neighbors = cv2.filter2D(skel_u8, cv2.CV_32F, diag_kernel, borderType=cv2.BORDER_CONSTANT)

    total_ortho_edges = float(np.sum(ortho_neighbors[skel_u8 > 0])) / 2.0
    total_diag_edges = float(np.sum(diag_neighbors[skel_u8 > 0])) / 2.0

    path_length = total_ortho_edges * 1.0 + total_diag_edges * np.sqrt(2.0)
    if path_length < 1.0:
        path_length = float(num_skel_pts)

    # Chord length: max pairwise Euclidean distance between endpoints
    end_pts = np.argwhere(end_mask > 0)
    if len(end_pts) >= 2:
        # Pairwise euclidean distances between all endpoints
        diff = end_pts[:, np.newaxis, :] - end_pts[np.newaxis, :, :]
        dist_matrix = np.sqrt(np.sum(diff ** 2, axis=-1))
        chord_length = float(np.max(dist_matrix))
    elif len(end_pts) == 1:
        # Distance from endpoint to furthest skeleton pixel
        diff = skel_pts - end_pts[0]
        dist_matrix = np.sqrt(np.sum(diff ** 2, axis=-1))
        chord_length = float(np.max(dist_matrix))
    else:
        # Loop or cycle without endpoints: use maximum span of skeleton
        diff = skel_pts[:, np.newaxis, :] - skel_pts[np.newaxis, :, :]
        dist_matrix = np.sqrt(np.sum(diff ** 2, axis=-1))
        chord_length = float(np.max(dist_matrix))

    if chord_length <= 1e-6:
        chord_length = 1.0

    tortuosity = max(1.0, float(path_length / chord_length))
    return path_length, chord_length, tortuosity


def classify_elongated_botanical_organ(
    mask: np.ndarray,
    full_sheet_height: float,
    y_centroid: Optional[float] = None,
    has_connected_blade: bool = False
) -> str:
    """
    Classifies a narrow elongated botanical organ into root_rhizome, cauline_stem, or leaf_petiole.

    Rules:
    1. If `has_connected_blade` is True, classify immediately as `leaf_petiole`.
    2. Perform medial axis skeletonization via `skimage.morphology.skeletonize(mask > 0)`.
    3. Convolve 3x3 neighbor-counting kernel `[[1,1,1],[1,10,1],[1,1,1]]` across skeleton
       to identify branch points (pixel value >= 13) and end points (pixel value == 11).
    4. Calculate normalized vertical position on the full sheet: norm_y = y_centroid / H_sheet.
    5. Classification logic:
       - If norm_y > 0.65 and branch points >= 3: return `root_rhizome` (subterranean, highly branched).
       - If norm_y < 0.70 and branch points <= 2: return `cauline_stem` (continuous vertical axis).
       - Otherwise, fallback to `leaf_petiole` if norm_y > 0.50 else `cauline_stem`.

    Args:
        mask: 2D numpy array representing the candidate organ instance mask.
        full_sheet_height: Total height H_sheet of the voucher image in pixels.
        y_centroid: Absolute vertical centroid position of the organ instance. If None,
                    it is calculated directly from the mask.
        has_connected_blade: Boolean flag indicating whether the organ directly attaches
                             to a segmented lamina/leaf blade.

    Returns:
        One of 'leaf_petiole', 'root_rhizome', or 'cauline_stem'.
    """
    # Rule 1: Immediate classification for connected leaf blades
    if has_connected_blade:
        return LEAF_PETIOLE

    # Rule 2: Medial axis skeletonization
    skel = extract_medial_axis_skeleton(mask)

    # Rule 3: 3x3 Neighbor-counting convolution kernel
    _, _, _, num_branch_points, _, _ = analyze_skeleton_topology(skel)

    # Determine y_centroid if not provided
    if y_centroid is None:
        _, y_centroid = calculate_mask_centroid(mask)

    # Rule 4: Normalized vertical sheet coordinate
    h_sheet = float(full_sheet_height) if full_sheet_height and full_sheet_height > 0 else 1.0
    norm_y = float(y_centroid) / h_sheet

    # Rule 5: Botanical classification decision boundary
    if norm_y > DEFAULT_ROOT_NORM_Y_THRESHOLD and num_branch_points >= DEFAULT_ROOT_MIN_BRANCH_POINTS:
        return ROOT_RHIZOME
    if norm_y < DEFAULT_STEM_NORM_Y_THRESHOLD and num_branch_points <= DEFAULT_STEM_MAX_BRANCH_POINTS:
        return CAULINE_STEM

    # Fallback boundary
    if norm_y > DEFAULT_FALLBACK_PETIOLE_NORM_Y:
        return LEAF_PETIOLE
    return CAULINE_STEM


def compute_topological_summary(
    mask: np.ndarray,
    full_sheet_height: float,
    y_centroid: Optional[float] = None,
    pixel_size_mm: Optional[float] = None,
    has_connected_blade: bool = False
) -> TopologicalMetrics:
    """
    Extracts a comprehensive suite of topological and morphometric summary metrics
    for export into phenotypic tabular datasets.

    Args:
        mask: 2D numpy array representing the candidate organ mask.
        full_sheet_height: Total height H_sheet of the voucher image in pixels.
        y_centroid: Absolute vertical centroid position. Calculated from mask if None.
        pixel_size_mm: Resolution scaling factor (mm/pixel). If provided, metric lengths
                       and areas are calculated.
        has_connected_blade: Whether the organ attaches to a leaf blade.

    Returns:
        TopologicalMetrics dataclass with extracted features and classification.
    """
    skel = extract_medial_axis_skeleton(mask)
    _, branch_mask, end_mask, num_branch_pts, num_end_pts, num_junction_clusters = analyze_skeleton_topology(skel)

    x_c, computed_y_c = calculate_mask_centroid(mask)
    if y_centroid is None:
        y_centroid = computed_y_c

    h_sheet = float(full_sheet_height) if full_sheet_height and full_sheet_height > 0 else 1.0
    norm_y = float(y_centroid) / h_sheet

    # Classification call
    pred_class = classify_elongated_botanical_organ(
        mask=mask,
        full_sheet_height=full_sheet_height,
        y_centroid=y_centroid,
        has_connected_blade=has_connected_blade
    )

    path_len_px, chord_len_px, tortuosity = calculate_skeleton_path_and_chord_length(skel, end_mask)
    skel_pixel_count = int(np.count_nonzero(skel > 0))
    mask_area_px = int(np.count_nonzero(mask > 0)) if mask is not None else 0

    path_len_mm = (path_len_px * pixel_size_mm) if pixel_size_mm is not None else None
    chord_len_mm = (chord_len_px * pixel_size_mm) if pixel_size_mm is not None else None
    mask_area_mm2 = (mask_area_px * (pixel_size_mm ** 2)) if pixel_size_mm is not None else None

    return TopologicalMetrics(
        predicted_class=pred_class,
        num_branch_points=num_branch_pts,
        num_end_points=num_end_pts,
        num_junction_clusters=num_junction_clusters,
        skeleton_pixel_count=skel_pixel_count,
        path_length_px=path_len_px,
        chord_length_px=chord_len_px,
        tortuosity=tortuosity,
        norm_y=norm_y,
        y_centroid_px=float(y_centroid),
        x_centroid_px=float(x_c),
        mask_area_px=mask_area_px,
        has_connected_blade=has_connected_blade,
        pixel_size_mm=pixel_size_mm,
        path_length_mm=path_len_mm,
        chord_length_mm=chord_len_mm,
        mask_area_mm2=mask_area_mm2,
    )


# ---------------------------------------------------------------------------
# Synthetic Morphological Generators for Unit Testing & Pipeline Prototyping
# ---------------------------------------------------------------------------

def generate_synthetic_stem_mask(
    shape: Tuple[int, int] = (1000, 1000),
    start_pt: Tuple[int, int] = (500, 100),
    end_pt: Tuple[int, int] = (500, 600),
    thickness: int = 6,
    waviness_amplitude: float = 0.0,
    waviness_frequency: float = 0.02
) -> np.ndarray:
    """
    Generates a synthetic binary mask of a linear or gently curved cauline stem.
    """
    mask = np.zeros(shape, dtype=np.uint8)
    y_vals = np.arange(start_pt[1], end_pt[1] + 1)
    if len(y_vals) < 2:
        return mask

    x_base = np.linspace(start_pt[0], end_pt[0], len(y_vals))
    x_vals = x_base + waviness_amplitude * np.sin(waviness_frequency * (y_vals - start_pt[1]))
    x_vals = np.clip(np.round(x_vals).astype(np.int32), 0, shape[1] - 1)

    pts = np.column_stack((x_vals, y_vals)).astype(np.int32)
    for i in range(len(pts) - 1):
        cv2.line(mask, tuple(pts[i]), tuple(pts[i + 1]), 255, thickness)
    return mask


def generate_synthetic_root_mask(
    shape: Tuple[int, int] = (1000, 1000),
    origin: Tuple[int, int] = (500, 750),
    num_branches: int = 5,
    max_depth: int = 3,
    branch_length: int = 60,
    thickness: int = 4
) -> np.ndarray:
    """
    Generates a synthetic binary mask of a highly branched subterranean root/rhizome system.
    """
    mask = np.zeros(shape, dtype=np.uint8)

    def draw_branch(x: int, y: int, length: int, angle_rad: float, depth: int):
        if depth > max_depth or length < 10:
            return
        x_end = int(x + length * np.sin(angle_rad))
        y_end = int(y + length * np.cos(angle_rad))

        x_end = np.clip(x_end, 0, shape[1] - 1)
        y_end = np.clip(y_end, 0, shape[0] - 1)

        t = max(1, thickness - depth)
        cv2.line(mask, (x, y), (x_end, y_end), 255, t)

        # Spawn sub-branches
        spread = 0.45
        draw_branch(x_end, y_end, int(length * 0.75), angle_rad - spread, depth + 1)
        draw_branch(x_end, y_end, int(length * 0.75), angle_rad + spread, depth + 1)
        if depth == 1:
            draw_branch(x_end, y_end, int(length * 0.85), angle_rad, depth + 1)

    # Initial main root axis
    main_len = 100
    cv2.line(mask, origin, (origin[0], min(shape[0] - 1, origin[1] + main_len)), 255, thickness)

    # Spawn primary lateral roots
    for i in range(num_branches):
        y_node = origin[1] + int((i + 1) * (main_len / (num_branches + 1)))
        draw_branch(origin[0], y_node, branch_length, -0.6, depth=1)
        draw_branch(origin[0], y_node, branch_length, 0.6, depth=1)

    return mask


def generate_synthetic_petiole_with_blade_mask(
    shape: Tuple[int, int] = (1000, 1000),
    stalk_start: Tuple[int, int] = (500, 600),
    stalk_end: Tuple[int, int] = (420, 520),
    blade_radius: int = 35,
    stalk_thickness: int = 4
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generates a synthetic petiole stalk mask and an attached leaf blade mask.

    Returns:
        Tuple of (petiole_mask, blade_mask).
    """
    petiole_mask = np.zeros(shape, dtype=np.uint8)
    blade_mask = np.zeros(shape, dtype=np.uint8)

    cv2.line(petiole_mask, stalk_start, stalk_end, 255, stalk_thickness)
    cv2.circle(blade_mask, stalk_end, blade_radius, 255, -1)

    return petiole_mask, blade_mask


if __name__ == "__main__":
    print("=" * 70)
    print("BOTANICAL ORGAN TOPOLOGY CLASSIFIER - VERIFICATION RUN")
    print("=" * 70)

    # 1. Stem
    stem = generate_synthetic_stem_mask((1000, 1000), (500, 100), (500, 500), thickness=6)
    stem_metrics = compute_topological_summary(stem, full_sheet_height=1000, pixel_size_mm=0.05)
    print(f"[Stem]     Class: {stem_metrics.predicted_class:15s} | norm_y: {stem_metrics.norm_y:.3f} | Branches: {stem_metrics.num_branch_points:2d} | Tortuosity: {stem_metrics.tortuosity:.3f}")

    # 2. Root
    root = generate_synthetic_root_mask((1000, 1000), (500, 720), num_branches=4, thickness=5)
    root_metrics = compute_topological_summary(root, full_sheet_height=1000, pixel_size_mm=0.05)
    print(f"[Root]     Class: {root_metrics.predicted_class:15s} | norm_y: {root_metrics.norm_y:.3f} | Branches: {root_metrics.num_branch_points:2d} | Tortuosity: {root_metrics.tortuosity:.3f}")

    # 3. Petiole with Blade
    petiole, _ = generate_synthetic_petiole_with_blade_mask((1000, 1000))
    petiole_metrics = compute_topological_summary(petiole, full_sheet_height=1000, pixel_size_mm=0.05, has_connected_blade=True)
    print(f"[Petiole]  Class: {petiole_metrics.predicted_class:15s} | norm_y: {petiole_metrics.norm_y:.3f} | Branches: {petiole_metrics.num_branch_points:2d} | Tortuosity: {petiole_metrics.tortuosity:.3f}")

    print("=" * 70)
    print("All synthetic organ demonstrations completed successfully.")

