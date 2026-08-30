#!/usr/bin/env python3
"""
scripts/vision/spatial_clustering.py
===================================
Adaptive DBSCAN spatial clustering of detected botanical organs on herbarium sheets.
Groups components into plant_individual_id clusters to prevent trait averaging across specimens.
"""

from __future__ import annotations

import logging
from typing import List

import numpy as np
from sklearn.cluster import DBSCAN

try:
    from scripts.vision.lm2_data_loader import LeafCandidate
except ImportError:
    from lm2_data_loader import LeafCandidate

logger = logging.getLogger("LM2_PostProcessing")


def cluster_voucher_plants_dbscan(
    candidates: List[LeafCandidate],
    sheet_width: int = 3000,
    sheet_height: int = 4000
) -> List[LeafCandidate]:
    """
    Applies adaptive DBSCAN spatial clustering (eps ≈ 0.15 * sheet_width, min_samples=2)
    to group detected organs into distinct plant_individual_id clusters.
    """
    if not candidates:
        return candidates

    # Extract centroids
    centroids = []
    for cand in candidates:
        ymin, xmin, ymax, xmax = cand.bbox
        if ymax > ymin and xmax > xmin:
            cx = (xmin + xmax) / 2.0
            cy = (ymin + ymax) / 2.0
        else:
            cx = sheet_width / 2.0
            cy = sheet_height / 2.0
        centroids.append([cx, cy])

    coords = np.array(centroids, dtype=np.float32)

    # Single leaf sheet
    if len(candidates) <= 1:
        candidates[0].plant_individual_id = 0
        return candidates

    # Adaptive eps = 0.15 * sheet_width
    adaptive_eps = 0.15 * float(sheet_width)
    db = DBSCAN(eps=adaptive_eps, min_samples=2).fit(coords)
    labels = db.labels_

    # Map cluster IDs: multi-organ clusters get 0, 1, 2, ...
    # Noise outliers (label == -1) get distinct individual IDs so no leaf is lost
    unique_clusters = sorted([l for l in set(labels) if l >= 0])
    cluster_mapping = {c: i for i, c in enumerate(unique_clusters)}

    next_id = len(cluster_mapping)
    for i, cand in enumerate(candidates):
        raw_label = labels[i]
        if raw_label >= 0:
            cand.plant_individual_id = cluster_mapping[raw_label]
        else:
            cand.plant_individual_id = next_id
            next_id += 1

    return candidates
