#!/usr/bin/env python3
"""
===============================================================================
Script: annotate_with_sam2.py
Project: Packera dubia Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Interactive Botanical Instance Segmentation Annotator powered by Segment
    Anything Model 2 (SAM 2) with ultra-high performance native X11 GUI rendering,
    viewport-first compositing, cached base layer rendering, custom polygon bounding box
    ROI selection, point prompts, knife slicing, and lag-free 60+ FPS zoom/pan navigation.
===============================================================================
"""

from __future__ import annotations

import argparse
import ctypes
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np
import torch


def get_project_root() -> Path:
    """Dynamically resolves the project root directory."""
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "data").exists() or (parent / "models").exists() or (parent / ".git").exists():
            return parent
    return current.parents[1] if len(current.parents) > 1 else current.parents[0]


PROJECT_ROOT = get_project_root()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Fallback path inclusion for SAM 2 if not installed in current site-packages
SAM2_ARCHIVE_PATH = PROJECT_ROOT / "scripts" / "_archive" / "root_artifacts" / "segment-anything-2"
if SAM2_ARCHIVE_PATH.exists() and str(SAM2_ARCHIVE_PATH) not in sys.path:
    sys.path.insert(0, str(SAM2_ARCHIVE_PATH))

try:
    from sam2.build_sam import build_sam2
    from sam2.sam2_image_predictor import SAM2ImagePredictor
except ImportError as err:
    build_sam2 = None
    SAM2ImagePredictor = None

from scripts.data_prep.sam2_geometry import (
    clip_box_to_image,
    mask_to_normalized_polygon,
    mask_to_yolo_bbox,
    polygon_interior_point,
    polygon_to_bounding_box,
    rasterize_lasso_polygon,
    split_mask_with_knife_line,
)
from scripts.data_prep.sam2_rendering import (
    CLASS_COLORS,
    CLASS_NAMES,
    apply_viewport_transform,
    compose_mask_overlay,
    overlay_candidate_mask_on_viewport,
    render_hud_overlay,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("SAM2Annotator")


class X11GUIWindow:
    """
    High-performance native X11 window for interactive botanical annotation,
    completely independent of Qt HighGUI to prevent phantom canvas panning.
    """

    def __init__(self, title: str = "SAM 2 Precision Botanical Annotator", width: int = 1280, height: int = 800):
        self.width = width
        self.height = height
        self.title = title
        self._x11 = ctypes.CDLL("libX11.so.6")
        self.disp = None
        self.win = None
        self.gc = None
        self._init_x11()

    def _init_x11(self) -> None:
        self._x11.XOpenDisplay.restype = ctypes.c_void_p
        self._x11.XOpenDisplay.argtypes = [ctypes.c_char_p]
        self._x11.XDefaultScreen.restype = ctypes.c_int
        self._x11.XDefaultScreen.argtypes = [ctypes.c_void_p]
        self._x11.XDefaultRootWindow.restype = ctypes.c_ulong
        self._x11.XDefaultRootWindow.argtypes = [ctypes.c_void_p]
        self._x11.XDefaultVisual.restype = ctypes.c_void_p
        self._x11.XDefaultVisual.argtypes = [ctypes.c_void_p, ctypes.c_int]
        self._x11.XDefaultDepth.restype = ctypes.c_int
        self._x11.XDefaultDepth.argtypes = [ctypes.c_void_p, ctypes.c_int]
        self._x11.XCreateSimpleWindow.restype = ctypes.c_ulong
        self._x11.XCreateSimpleWindow.argtypes = [
            ctypes.c_void_p, ctypes.c_ulong, ctypes.c_int, ctypes.c_int,
            ctypes.c_uint, ctypes.c_uint, ctypes.c_uint, ctypes.c_ulong, ctypes.c_ulong
        ]
        self._x11.XSelectInput.restype = ctypes.c_int
        self._x11.XSelectInput.argtypes = [ctypes.c_void_p, ctypes.c_ulong, ctypes.c_long]
        self._x11.XMapWindow.restype = ctypes.c_int
        self._x11.XMapWindow.argtypes = [ctypes.c_void_p, ctypes.c_ulong]
        self._x11.XCreateGC.restype = ctypes.c_void_p
        self._x11.XCreateGC.argtypes = [ctypes.c_void_p, ctypes.c_ulong, ctypes.c_ulong, ctypes.c_void_p]
        self._x11.XFreeGC.restype = ctypes.c_int
        self._x11.XFreeGC.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
        self._x11.XCreateImage.restype = ctypes.c_void_p
        self._x11.XCreateImage.argtypes = [
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_uint, ctypes.c_int,
            ctypes.c_int, ctypes.c_char_p, ctypes.c_uint, ctypes.c_uint,
            ctypes.c_int, ctypes.c_int
        ]
        self._x11.XPutImage.restype = ctypes.c_int
        self._x11.XPutImage.argtypes = [
            ctypes.c_void_p, ctypes.c_ulong, ctypes.c_void_p, ctypes.c_void_p,
            ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
            ctypes.c_uint, ctypes.c_uint
        ]
        self._x11.XFlush.restype = ctypes.c_int
        self._x11.XFlush.argtypes = [ctypes.c_void_p]
        self._x11.XEventsQueued.restype = ctypes.c_int
        self._x11.XEventsQueued.argtypes = [ctypes.c_void_p, ctypes.c_int]
        self._x11.XNextEvent.restype = ctypes.c_int
        self._x11.XNextEvent.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
        self._x11.XKeycodeToKeysym.restype = ctypes.c_ulong
        self._x11.XKeycodeToKeysym.argtypes = [ctypes.c_void_p, ctypes.c_ubyte, ctypes.c_int]
        self._x11.XDestroyWindow.restype = ctypes.c_int
        self._x11.XDestroyWindow.argtypes = [ctypes.c_void_p, ctypes.c_ulong]
        self._x11.XCloseDisplay.restype = ctypes.c_int
        self._x11.XCloseDisplay.argtypes = [ctypes.c_void_p]
        self._x11.XStoreName.restype = ctypes.c_int
        self._x11.XStoreName.argtypes = [ctypes.c_void_p, ctypes.c_ulong, ctypes.c_char_p]
        self._x11.XInternAtom.restype = ctypes.c_ulong
        self._x11.XInternAtom.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
        self._x11.XSetWMProtocols.restype = ctypes.c_int
        self._x11.XSetWMProtocols.argtypes = [ctypes.c_void_p, ctypes.c_ulong, ctypes.POINTER(ctypes.c_ulong), ctypes.c_int]

        self.disp = self._x11.XOpenDisplay(None)
        if not self.disp:
            raise RuntimeError("Failed to open X11 Display. Please verify DISPLAY environment variable.")
        self.screen = self._x11.XDefaultScreen(self.disp)
        self.root = self._x11.XDefaultRootWindow(self.disp)
        self.visual = self._x11.XDefaultVisual(self.disp, self.screen)
        self.depth = self._x11.XDefaultDepth(self.disp, self.screen)

        self.win = self._x11.XCreateSimpleWindow(self.disp, self.root, 50, 50, self.width, self.height, 1, 0, 0)
        self._x11.XStoreName(self.disp, self.win, self.title.encode("utf-8"))

        self.wm_delete = self._x11.XInternAtom(self.disp, b"WM_DELETE_WINDOW", False)
        atom_arr = (ctypes.c_ulong * 1)(self.wm_delete)
        self._x11.XSetWMProtocols(self.disp, self.win, atom_arr, 1)

        # ExposureMask | KeyPressMask | KeyReleaseMask | ButtonPressMask | ButtonReleaseMask | PointerMotionMask | StructureNotifyMask
        event_mask = (1 << 15) | (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 6) | (1 << 17)
        self._x11.XSelectInput(self.disp, self.win, event_mask)
        self._x11.XMapWindow(self.disp, self.win)

        self.gc = self._x11.XCreateGC(self.disp, self.win, 0, None)

        self.bgra_buffer = np.zeros((self.height, self.width, 4), dtype=np.uint8)
        self.ximage = self._x11.XCreateImage(
            self.disp, self.visual, self.depth, 2,  # ZPixmap
            0, self.bgra_buffer.ctypes.data_as(ctypes.c_char_p),
            self.width, self.height, 32, 0
        )

    def imshow(self, bgr_img: np.ndarray) -> None:
        """Transfers BGR image directly to native X11 window buffer and flushes to display."""
        h, w = bgr_img.shape[:2]
        if h != self.height or w != self.width:
            bgr_img = cv2.resize(bgr_img, (self.width, self.height))
        self.bgra_buffer[:, :, :3] = bgr_img
        self.bgra_buffer[:, :, 3] = 255
        self._x11.XPutImage(self.disp, self.win, self.gc, self.ximage, 0, 0, 0, 0, self.width, self.height)
        self._x11.XFlush(self.disp)

    def poll_events(self) -> List[Tuple[str, Any]]:
        """
        Polls queued X11 events, coalescing intermediate mouse motion to guarantee
        sub-millisecond responsiveness during mouse dragging and canvas panning.
        """
        raw_events = []
        if not hasattr(self, "disp") or not self.disp:
            return raw_events

        class XEvent(ctypes.Structure):
            _fields_ = [
                ('type', ctypes.c_int),
                ('serial', ctypes.c_ulong),
                ('send_event', ctypes.c_int),
                ('display', ctypes.c_void_p),
                ('window', ctypes.c_ulong),
                ('root', ctypes.c_ulong),
                ('subwindow', ctypes.c_ulong),
                ('time', ctypes.c_ulong),
                ('x', ctypes.c_int),
                ('y', ctypes.c_int),
                ('x_root', ctypes.c_int),
                ('y_root', ctypes.c_int),
                ('state', ctypes.c_uint),
                ('button_or_keycode', ctypes.c_uint),
                ('same_screen', ctypes.c_int),
                ('pad', ctypes.c_ulong * 10)
            ]

        evt = XEvent()
        QueuedAfterReading = 1
        latest_motion = None

        while self._x11.XEventsQueued(self.disp, QueuedAfterReading) > 0:
            self._x11.XNextEvent(self.disp, ctypes.byref(evt))
            if evt.type == 2:  # KeyPress
                keysym = self._x11.XKeycodeToKeysym(self.disp, evt.button_or_keycode, 0)
                raw_events.append(('key_press', keysym))
            elif evt.type == 4:  # ButtonPress
                btn = evt.button_or_keycode
                raw_events.append(('button_press', evt.x, evt.y, btn, evt.state))
            elif evt.type == 5:  # ButtonRelease
                btn = evt.button_or_keycode
                raw_events.append(('button_release', evt.x, evt.y, btn, evt.state))
            elif evt.type == 6:  # MotionNotify (coalesce)
                latest_motion = ('motion', evt.x, evt.y, evt.state)
            elif evt.type == 33:  # ClientMessage (WM_DELETE_WINDOW)
                raw_events.append(('close', None))

        if latest_motion is not None:
            raw_events.append(latest_motion)

        return raw_events

    def close(self) -> None:
        """Safely cleans up X11 window and display connection without double-freeing buffer."""
        if hasattr(self, "disp") and self.disp:
            try:
                if hasattr(self, "gc") and self.gc:
                    self._x11.XFreeGC(self.disp, self.gc)
                    self.gc = None
                if hasattr(self, "win") and self.win:
                    self._x11.XDestroyWindow(self.disp, self.win)
                    self.win = None
                self._x11.XCloseDisplay(self.disp)
            except Exception:
                pass
            self.disp = None


class PrecisionSAM2Annotator:
    """
    Interactive botanical annotator using SAM 2 to target Packera specimens,
    accepting point, rectangular bounding box, and custom polygonal bounding box ROI prompts.
    Outputs binary pixel masks tagged explicitly with taxonomic labels.
    """

    def __init__(
        self,
        images_dir: Union[str, Path] = "data/raw_vouchers",
        output_dir: Union[str, Path] = "data/raw_annotations",
        single_image: Optional[Union[str, Path]] = None,
        checkpoint_path: Union[str, Path] = "models/checkpoints/sam2_hiera_large.pt",
        config_path: Union[str, Path] = "sam2_hiera_l.yaml",
        window_w: int = 1280,
        window_h: int = 800,
        resume_unannotated: bool = True
    ):
        self.project_root = get_project_root()
        self.images_dir = Path(images_dir)
        self.output_dir = Path(output_dir)
        self.masks_dir = self.output_dir / "masks"
        self.labels_dir = self.output_dir / "labels"
        self.masks_dir.mkdir(parents=True, exist_ok=True)
        self.labels_dir.mkdir(parents=True, exist_ok=True)

        self.window_w = window_w
        self.window_h = window_h
        self.resume_unannotated = resume_unannotated

        if single_image:
            self.image_files = [Path(single_image)]
        else:
            exts = {".jpg", ".jpeg", ".png", ".tif", ".tiff"}
            self.image_files = sorted([
                p for p in self.images_dir.glob("*.*")
                if p.suffix.lower() in exts and not p.name.startswith(".")
            ])

        self.checkpoint_path = Path(checkpoint_path)
        self.config_path = str(config_path)

        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.predictor = None
        self._init_model()

        self.current_idx = 0
        self.active_image: Optional[np.ndarray] = None
        self.cached_base_layer: Optional[np.ndarray] = None
        self.orig_h = 1000
        self.orig_w = 1000

        self.point_coords: List[List[float]] = []
        self.point_labels: List[int] = []
        self.box_prompt: Optional[List[float]] = None
        self.candidate_mask: Optional[np.ndarray] = None
        self.saved_instances: List[Dict[str, Any]] = []

        self.zoom_level = 1.0
        self.pan_offset = [0, 0]
        self.mode = "SELECT"

        # Standard interaction states
        self.lbutton_down = False
        self.drag_start_screen = (0, 0)
        self.drag_start_img = (0, 0)
        self.is_box_dragging = False
        self.is_knife_dragging = False
        self.knife_start = (0, 0)
        self.knife_current = (0, 0)
        self.is_pan_dragging = False
        self.pan_drag_start = (0, 0)
        self.pan_offset_start = [0, 0]

        # Custom polygon bounding box selection states
        self.polygon_points: List[Tuple[int, int]] = []
        self.hover_img_pos: Tuple[int, int] = (0, 0)
        self.hover_screen_pos: Tuple[int, int] = (0, 0)

    def _init_model(self) -> None:
        """Initializes the SAM 2 model weights and image predictor."""
        if build_sam2 is None or SAM2ImagePredictor is None:
            return

        if self.checkpoint_path.exists():
            try:
                sam2_model = build_sam2(self.config_path, str(self.checkpoint_path), device=self.device)
                self.predictor = SAM2ImagePredictor(sam2_model)
                logger.info(f"Loaded SAM 2 model ({self.config_path}) onto {self.device}")
            except Exception as e:
                logger.error(f"Failed to load SAM 2 weights: {e}")

    def _update_cached_base_layer(self) -> None:
        """
        Pre-composites all saved instances onto the base image once to allow instant
        sub-millisecond viewport cropping without re-blending full resolution masks.
        """
        if self.active_image is None:
            self.cached_base_layer = None
            return
        if not self.saved_instances:
            self.cached_base_layer = self.active_image.copy()
            return

        alpha = 0.45
        overlay = self.active_image.copy()
        for inst in self.saved_instances:
            mask = inst["mask"]
            class_id = inst["class_id"]
            color = CLASS_COLORS.get(class_id, (0, 255, 0))

            locs = mask > 0
            if np.any(locs):
                overlay[locs] = (overlay[locs] * (1.0 - alpha) + np.array(color, dtype=np.uint8) * alpha).astype(np.uint8)
                contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                cv2.drawContours(overlay, contours, -1, color, 2)

        self.cached_base_layer = overlay

    def run_inference(self) -> None:
        """Executes SAM 2 inference on current prompts."""
        if self.predictor is None or self.active_image is None:
            return

        pts = np.array(self.point_coords, dtype=np.float32) if self.point_coords else None
        lbls = np.array(self.point_labels, dtype=np.int32) if self.point_labels else None
        box = np.array(self.box_prompt, dtype=np.float32) if self.box_prompt else None

        if pts is None and box is None:
            self.candidate_mask = None
            return

        try:
            masks, scores, _ = self.predictor.predict(
                point_coords=pts,
                point_labels=lbls,
                box=box,
                multimask_output=False
            )
            if masks is not None and len(masks) > 0:
                self.candidate_mask = (masks[0] > 0.0).astype(np.uint8) * 255
        except Exception as e:
            logger.error(f"SAM 2 prediction error: {e}")

    def finalize_polygon_selection(self) -> None:
        """
        Converts marked custom polygon selection points into a precise SAM 2 bounding box
        prompt and constrains the candidate mask within the custom polygon boundary.
        """
        if len(self.polygon_points) < 3:
            return

        poly_box = polygon_to_bounding_box(self.polygon_points)
        if not poly_box:
            return

        bx0, by0, bx1, by1 = poly_box
        if abs(bx1 - bx0) < 3 or abs(by1 - by0) < 3:
            return

        self.box_prompt = [float(bx0), float(by0), float(bx1), float(by1)]

        # Determine interior guide point (pole of inaccessibility)
        interior_pt = polygon_interior_point(self.polygon_points, self.orig_h, self.orig_w)
        if interior_pt:
            self.point_coords = [[interior_pt[0], interior_pt[1]]]
            self.point_labels = [1]
        else:
            self.point_coords = []
            self.point_labels = []

        # Run SAM 2 segmentation inference
        self.run_inference()

        # Constrain resulting mask to user's custom polygon
        poly_mask = rasterize_lasso_polygon(self.polygon_points, self.orig_h, self.orig_w)
        if self.candidate_mask is not None and np.count_nonzero(self.candidate_mask) > 0:
            constrained = cv2.bitwise_and(self.candidate_mask, poly_mask)
            if np.count_nonzero(constrained) > 0:
                self.candidate_mask = constrained
            else:
                self.candidate_mask = poly_mask
        else:
            self.candidate_mask = poly_mask

        logger.info(f"Finalized custom polygon bounding box with {len(self.polygon_points)} vertices")
        self.polygon_points = []

    def save_current_sheet(self) -> None:
        """Saves current sheet's instances to YOLO polygon text file and binary PNG masks."""
        if not hasattr(self, "image_files") or not self.image_files:
            return
        if self.current_idx >= len(self.image_files):
            return

        current_file = self.image_files[self.current_idx]
        voucher_id = current_file.stem
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.masks_dir.mkdir(parents=True, exist_ok=True)

        txt_file = self.output_dir / f"{voucher_id}.txt"
        with open(txt_file, "w", encoding="utf-8") as f:
            for inst in self.saved_instances:
                poly = inst.get("polygon", [])
                c_id = inst.get("class_id", 0)
                poly_str = " ".join([str(round(v, 6)) for v in poly])
                f.write(f"{c_id} {poly_str}\n")

        for idx, inst in enumerate(self.saved_instances):
            label = inst.get("label", CLASS_NAMES[inst.get("class_id", 0)])
            mask_dest = self.masks_dir / f"{voucher_id}_inst{idx:02d}_{label}.png"
            b_mask = inst.get("binary_mask")
            if b_mask is not None:
                uint8_mask = (b_mask.astype(np.uint8)) * 255 if b_mask.dtype == bool else b_mask.astype(np.uint8)
                cv2.imwrite(str(mask_dest), uint8_mask)

        logger.info(f"Saved {len(self.saved_instances)} instances for voucher {voucher_id}")

    def load_active_image(self) -> bool:
        """Loads current voucher sheet into memory and feeds it to SAM 2 image predictor."""
        if not hasattr(self, "image_files") or not self.image_files:
            return False
        if self.current_idx >= len(self.image_files):
            return False

        img_path = self.image_files[self.current_idx]
        self.active_image = cv2.imread(str(img_path))
        if self.active_image is None:
            logger.error(f"Could not load image from {img_path}")
            return False

        self.orig_h, self.orig_w = self.active_image.shape[:2]
        self.point_coords = []
        self.point_labels = []
        self.box_prompt = None
        self.candidate_mask = None
        self.saved_instances = []
        self.zoom_level = 1.0
        self.pan_offset = [0, 0]
        self.polygon_points = []

        self._update_cached_base_layer()

        if self.predictor is not None:
            img_rgb = cv2.cvtColor(self.active_image, cv2.COLOR_BGR2RGB)
            self.predictor.set_image(img_rgb)
            logger.info(f"Loaded voucher [{self.current_idx + 1}/{len(self.image_files)}]: {img_path.name}")
        return True

    def _viewport_to_image(
        self,
        vx: int,
        vy: int,
        scale_x: float,
        scale_y: float,
        crop_x0: int,
        crop_y0: int
    ) -> Tuple[int, int]:
        """Converts viewport display window coordinates back to native resolution image space."""
        ix = int(crop_x0 + vx / max(scale_x, 1e-6))
        iy = int(crop_y0 + vy / max(scale_y, 1e-6))
        ix = max(0, min(self.orig_w - 1, ix))
        iy = max(0, min(self.orig_h - 1, iy))
        return ix, iy

    def _image_to_viewport(
        self,
        ix: Union[int, float],
        iy: Union[int, float],
        scale_x: float,
        scale_y: float,
        crop_x0: int,
        crop_y0: int
    ) -> Tuple[int, int]:
        """Converts native resolution image coordinates to viewport display window space."""
        vx = int((ix - crop_x0) * scale_x)
        vy = int((iy - crop_y0) * scale_y)
        return vx, vy

    def run(self) -> None:
        """Main interactive GUI annotation event loop."""
        if not hasattr(self, "image_files") or not self.image_files:
            logger.warning(f"No valid voucher images found in {self.images_dir}")
            return

        window_title = "SAM 2 Precision Botanical Annotator"
        win = X11GUIWindow(title=window_title, width=self.window_w, height=self.window_h)

        if not self.load_active_image():
            logger.error("Failed to load initial voucher image.")
            win.close()
            return

        transform: Tuple[float, float, int, int] = (1.0, 1.0, 0, 0)

        try:
            while True:
                if self.active_image is None:
                    break

                current_voucher = self.image_files[self.current_idx].stem

                # 1. Fast sub-millisecond viewport crop from pre-composited base layer
                base_source = self.cached_base_layer if self.cached_base_layer is not None else self.active_image
                viewport_img, transform = apply_viewport_transform(
                    base_source,
                    self.zoom_level,
                    tuple(self.pan_offset),
                    self.window_w,
                    self.window_h
                )
                scale_x, scale_y, crop_x0, crop_y0 = transform

                # 2. Fast viewport-space candidate mask overlay
                if self.candidate_mask is not None:
                    viewport_img = overlay_candidate_mask_on_viewport(
                        viewport_img,
                        self.candidate_mask,
                        transform,
                        alpha=0.55
                    )

                # 3. Draw vector graphics directly in viewport space (0.1 ms)
                # Prompt points
                for pt, lbl in zip(self.point_coords, self.point_labels):
                    vx, vy = self._image_to_viewport(pt[0], pt[1], scale_x, scale_y, crop_x0, crop_y0)
                    if 0 <= vx < self.window_w and 0 <= vy < self.window_h:
                        color = (0, 255, 0) if lbl == 1 else (0, 0, 255)
                        cv2.circle(viewport_img, (vx, vy), 6, color, -1)
                        cv2.circle(viewport_img, (vx, vy), 8, (255, 255, 255), 1)

                # Bounding box prompt
                if self.box_prompt:
                    bx0, by0 = self._image_to_viewport(self.box_prompt[0], self.box_prompt[1], scale_x, scale_y, crop_x0, crop_y0)
                    bx1, by1 = self._image_to_viewport(self.box_prompt[2], self.box_prompt[3], scale_x, scale_y, crop_x0, crop_y0)
                    cv2.rectangle(viewport_img, (bx0, by0), (bx1, by1), (255, 200, 0), 2)

                # Knife line
                if getattr(self, "is_knife_dragging", False) and hasattr(self, "knife_start") and hasattr(self, "knife_current"):
                    k0 = self._image_to_viewport(self.knife_start[0], self.knife_start[1], scale_x, scale_y, crop_x0, crop_y0)
                    k1 = self._image_to_viewport(self.knife_current[0], self.knife_current[1], scale_x, scale_y, crop_x0, crop_y0)
                    cv2.line(viewport_img, k0, k1, (0, 0, 255), 3)

                # Custom polygon vertices and preview
                if self.polygon_points:
                    v_pts = [self._image_to_viewport(p[0], p[1], scale_x, scale_y, crop_x0, crop_y0) for p in self.polygon_points]

                    # Shaded preview fill
                    if len(v_pts) >= 3:
                        poly_ov = viewport_img.copy()
                        pts_arr = np.array(v_pts, dtype=np.int32).reshape((-1, 1, 2))
                        cv2.fillPoly(poly_ov, [pts_arr], (0, 220, 255))
                        cv2.addWeighted(poly_ov, 0.25, viewport_img, 0.75, 0, viewport_img)

                    # Connecting lines
                    for i in range(len(v_pts) - 1):
                        cv2.line(viewport_img, v_pts[i], v_pts[i + 1], (0, 240, 255), 2, cv2.LINE_AA)

                    # Dynamic live preview line to cursor position
                    if self.mode == "POLYGON" and hasattr(self, "hover_screen_pos"):
                        cv2.line(viewport_img, v_pts[-1], self.hover_screen_pos, (0, 200, 255), 1, cv2.LINE_AA)

                    # Vertex markers
                    for idx, (vx, vy) in enumerate(v_pts):
                        if 0 <= vx < self.window_w and 0 <= vy < self.window_h:
                            if idx == 0:
                                cv2.circle(viewport_img, (vx, vy), 7, (0, 255, 0), -1)
                                cv2.circle(viewport_img, (vx, vy), 9, (255, 255, 255), 2)
                            else:
                                cv2.circle(viewport_img, (vx, vy), 5, (0, 220, 255), -1)
                                cv2.circle(viewport_img, (vx, vy), 7, (255, 255, 255), 1)

                # 4. Render HUD status bar
                hud_display = render_hud_overlay(
                    viewport_img,
                    current_voucher,
                    self.current_idx,
                    len(self.image_files),
                    self.saved_instances,
                    self.mode,
                    self.zoom_level,
                    tuple(self.pan_offset),
                )

                # If in polygon mode with >= 3 vertices, draw start-point closure target hint
                if self.mode == "POLYGON" and len(self.polygon_points) >= 3:
                    start_vx, start_vy = self._image_to_viewport(
                        self.polygon_points[0][0], self.polygon_points[0][1],
                        scale_x, scale_y, crop_x0, crop_y0
                    )
                    if 0 <= start_vx < self.window_w and 70 <= start_vy < self.window_h:
                        dist_to_start = max(abs(self.hover_screen_pos[0] - start_vx), abs(self.hover_screen_pos[1] - start_vy))
                        ring_color = (0, 255, 0) if dist_to_start <= 15 else (200, 200, 200)
                        cv2.circle(hud_display, (start_vx, start_vy), 14, ring_color, 2, cv2.LINE_AA)
                        if dist_to_start <= 15:
                            cv2.putText(hud_display, "Click to Close", (start_vx + 16, start_vy + 5),
                                        cv2.FONT_HERSHEY_SIMPLEX, 0.45, (0, 255, 0), 1, cv2.LINE_AA)

                win.imshow(hud_display)

                # 5. Process all native X11 GUI events
                should_exit = False
                events = win.poll_events()

                for ev in events:
                    ev_type = ev[0]
                    if ev_type == "button_press":
                        _, vx, vy, btn, state = ev
                        ix, iy = self._viewport_to_image(vx, vy, scale_x, scale_y, crop_x0, crop_y0)

                        if btn == 1:  # Left Button Down
                            self.is_pan_dragging = False

                            if self.mode == "POLYGON":
                                if len(self.polygon_points) >= 3:
                                    start_vx, start_vy = self._image_to_viewport(
                                        self.polygon_points[0][0], self.polygon_points[0][1],
                                        scale_x, scale_y, crop_x0, crop_y0
                                    )
                                    dist_to_start = max(abs(vx - start_vx), abs(vy - start_vy))
                                    if dist_to_start <= 15:
                                        self.finalize_polygon_selection()
                                    else:
                                        self.polygon_points.append((ix, iy))
                                else:
                                    self.polygon_points.append((ix, iy))

                            elif self.mode == "KNIFE" or (state & 0x0004):
                                self.is_knife_dragging = True
                                self.knife_start = (ix, iy)
                                self.knife_current = (ix, iy)

                            else:
                                self.lbutton_down = True
                                self.drag_start_screen = (vx, vy)
                                self.drag_start_img = (ix, iy)
                                self.is_box_dragging = False

                        elif btn == 3:  # Right Button Down
                            self.is_pan_dragging = False
                            if self.mode == "POLYGON" and len(self.polygon_points) >= 3:
                                self.finalize_polygon_selection()
                            else:
                                self.point_coords.append([float(ix), float(iy)])
                                self.point_labels.append(0)
                                self.run_inference()

                        elif btn == 2:  # Middle Button Down -> Pan Canvas
                            if not getattr(self, "lbutton_down", False):
                                self.is_pan_dragging = True
                                self.pan_drag_start = (vx, vy)
                                self.pan_offset_start = list(self.pan_offset)

                        elif btn == 4:  # Wheel Up -> Zoom In
                            self.zoom_level = min(10.0, round(self.zoom_level * 1.25, 2))

                        elif btn == 5:  # Wheel Down -> Zoom Out
                            self.zoom_level = max(1.0, round(self.zoom_level / 1.25, 2))
                            if self.zoom_level <= 1.0:
                                self.pan_offset = [0, 0]

                    elif ev_type == "motion":
                        _, vx, vy, state = ev
                        ix, iy = self._viewport_to_image(vx, vy, scale_x, scale_y, crop_x0, crop_y0)
                        self.hover_img_pos = (ix, iy)
                        self.hover_screen_pos = (vx, vy)

                        if getattr(self, "lbutton_down", False) and self.mode != "POLYGON":
                            dist_screen = max(abs(vx - self.drag_start_screen[0]), abs(vy - self.drag_start_screen[1]))
                            if dist_screen > 5:
                                self.is_box_dragging = True
                                x0, y0 = self.drag_start_img
                                self.box_prompt = [float(min(x0, ix)), float(min(y0, iy)), float(max(x0, ix)), float(max(y0, iy))]

                        elif getattr(self, "is_knife_dragging", False):
                            self.knife_current = (ix, iy)

                        elif getattr(self, "is_pan_dragging", False):
                            dx = vx - self.pan_drag_start[0]
                            dy = vy - self.pan_drag_start[1]
                            img_dx = int(dx / max(scale_x, 1e-6))
                            img_dy = int(dy / max(scale_y, 1e-6))
                            crop_w = int(self.orig_w / max(self.zoom_level, 1.0))
                            crop_h = int(self.orig_h / max(self.zoom_level, 1.0))
                            if self.zoom_level > 1.0:
                                min_x = -int(crop_w * 0.85)
                                max_x = int(self.orig_w - crop_w * 0.15)
                                min_y = -int(crop_h * 0.85)
                                max_y = int(self.orig_h - crop_h * 0.15)
                                self.pan_offset[0] = max(min_x, min(max_x, self.pan_offset_start[0] - img_dx))
                                self.pan_offset[1] = max(min_y, min(max_y, self.pan_offset_start[1] - img_dy))
                            else:
                                self.pan_offset = [0, 0]

                    elif ev_type == "button_release":
                        _, vx, vy, btn, state = ev
                        ix, iy = self._viewport_to_image(vx, vy, scale_x, scale_y, crop_x0, crop_y0)

                        if btn == 1 and self.mode != "POLYGON":
                            if self.is_knife_dragging:
                                self.is_knife_dragging = False
                                if self.candidate_mask is not None:
                                    self.candidate_mask = split_mask_with_knife_line(
                                        self.candidate_mask, self.knife_start, (ix, iy), line_thickness=5
                                    )
                            elif getattr(self, "lbutton_down", False):
                                self.lbutton_down = False
                                if self.is_box_dragging:
                                    self.is_box_dragging = False
                                    x0, y0 = self.drag_start_img
                                    if abs(ix - x0) > 5 and abs(iy - y0) > 5:
                                        self.box_prompt = [float(min(x0, ix)), float(min(y0, iy)), float(max(x0, ix)), float(max(y0, iy))]
                                        self.run_inference()
                                    else:
                                        self.box_prompt = None
                                else:
                                    self.point_coords.append([float(ix), float(iy)])
                                    self.point_labels.append(1)
                                    self.run_inference()

                        elif btn == 2:
                            self.is_pan_dragging = False

                    elif ev_type == "key_press":
                        _, sym = ev

                        # Class assignment '0' through '6'
                        if ord('0') <= sym <= ord('6'):
                            class_id = sym - ord('0')
                            if class_id < len(CLASS_NAMES) and self.candidate_mask is not None and np.count_nonzero(self.candidate_mask) > 0:
                                label = CLASS_NAMES[class_id]
                                contours, _ = cv2.findContours(self.candidate_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                                if contours:
                                    largest = max(contours, key=cv2.contourArea)
                                    poly = []
                                    for p in largest.reshape(-1, 2):
                                        poly.append(float(p[0]) / self.orig_w)
                                        poly.append(float(p[1]) / self.orig_h)
                                    self.saved_instances.append({
                                        "class_id": class_id,
                                        "label": label,
                                        "polygon": poly,
                                        "binary_mask": self.candidate_mask > 0,
                                        "mask": self.candidate_mask.copy()
                                    })
                                    self.candidate_mask = None
                                    self.point_coords = []
                                    self.point_labels = []
                                    self.box_prompt = None
                                    self.polygon_points = []
                                    self._update_cached_base_layer()
                                    logger.info(f"Assigned instance #{len(self.saved_instances)} -> '{label}' (Class {class_id})")

                        # 'P' or 'p': Toggle Custom Polygon Bounding Box Mode
                        elif sym in (ord('p'), ord('P')):
                            self.mode = "POLYGON" if self.mode != "POLYGON" else "SELECT"
                            if self.mode != "POLYGON":
                                self.polygon_points = []
                            logger.info(f"Mode switched to {self.mode}")

                        # 'K' or 'k': Toggle Knife Slicing Mode
                        elif sym in (ord('k'), ord('K')):
                            self.mode = "KNIFE" if self.mode != "KNIFE" else "SELECT"
                            logger.info(f"Mode switched to {self.mode}")

                        # 'Enter' / 'Return' / 'Space'
                        elif sym in (0xff0d, 0xff8d, 0x20, ord('f'), ord('F')):
                            if self.mode == "POLYGON" and len(self.polygon_points) >= 3:
                                self.finalize_polygon_selection()
                            elif self.mode != "POLYGON" and sym in (0xff0d, 0xff8d):
                                self.save_current_sheet()
                                self.current_idx += 1
                                if self.current_idx < len(self.image_files):
                                    self.load_active_image()
                                else:
                                    logger.info("Finished annotating all vouchers in directory.")
                                    should_exit = True
                                    break

                        # 'Backspace' / 'Delete': Undo last vertex in polygon mode
                        elif sym in (0xff08, 0xffff):
                            if self.mode == "POLYGON" and self.polygon_points:
                                self.polygon_points.pop()
                                logger.info(f"Removed vertex. Remaining polygon vertices: {len(self.polygon_points)}")

                        # 'C' or 'c': Clear current candidate prompts or polygon points
                        elif sym in (ord('c'), ord('C')):
                            if self.polygon_points:
                                self.polygon_points = []
                                logger.info("Cleared custom polygon selection points.")
                            else:
                                self.point_coords = []
                                self.point_labels = []
                                self.box_prompt = None
                                self.candidate_mask = None
                                logger.info("Cleared active prompts and candidate mask.")

                        # 'U' or 'u': Undo last saved instance
                        elif sym in (ord('u'), ord('U')):
                            if self.saved_instances:
                                popped = self.saved_instances.pop()
                                self._update_cached_base_layer()
                                logger.info(f"Removed instance: {popped.get('label')}")

                        # Zoom in / out
                        elif sym in (ord('z'), ord('Z')):
                            self.zoom_level = min(10.0, round(self.zoom_level * 1.3, 2))

                        elif sym in (ord('x'), ord('X')):
                            self.zoom_level = max(1.0, round(self.zoom_level / 1.3, 2))
                            if self.zoom_level <= 1.0:
                                self.pan_offset = [0, 0]

                        # Panning: W / A / S / D & Arrow Keys
                        elif sym in (ord('w'), ord('W'), 0xff52):  # Up
                            if self.zoom_level > 1.0:
                                crop_h = int(self.orig_h / self.zoom_level)
                                min_y = -int(crop_h * 0.85)
                                step = max(20, int(150 / self.zoom_level))
                                self.pan_offset[1] = max(min_y, self.pan_offset[1] - step)

                        elif sym in (ord('s'), ord('S'), 0xff54):  # Down
                            if self.zoom_level > 1.0:
                                crop_h = int(self.orig_h / self.zoom_level)
                                max_y = int(self.orig_h - crop_h * 0.15)
                                step = max(20, int(150 / self.zoom_level))
                                self.pan_offset[1] = min(max_y, self.pan_offset[1] + step)

                        elif sym in (ord('a'), ord('A'), 0xff51):  # Left
                            if self.zoom_level > 1.0:
                                crop_w = int(self.orig_w / self.zoom_level)
                                min_x = -int(crop_w * 0.85)
                                step = max(20, int(150 / self.zoom_level))
                                self.pan_offset[0] = max(min_x, self.pan_offset[0] - step)

                        elif sym in (ord('d'), ord('D'), 0xff53):  # Right
                            if self.zoom_level > 1.0:
                                crop_w = int(self.orig_w / self.zoom_level)
                                max_x = int(self.orig_w - crop_w * 0.15)
                                step = max(20, int(150 / self.zoom_level))
                                self.pan_offset[0] = min(max_x, self.pan_offset[0] + step)

                        # Voucher navigation
                        elif sym in (ord('n'), ord('N'), ord('v'), ord('V')):  # Next voucher
                            self.save_current_sheet()
                            self.current_idx += 1
                            if self.current_idx < len(self.image_files):
                                self.load_active_image()
                            else:
                                logger.info("Finished annotating all vouchers in directory.")
                                should_exit = True
                                break

                        elif sym in (ord('b'), ord('B')):  # Previous voucher
                            self.save_current_sheet()
                            if self.current_idx > 0:
                                self.current_idx -= 1
                                self.load_active_image()

                        # Quit / Esc
                        elif sym in (ord('q'), ord('Q'), 0xff1b):
                            if self.mode == "POLYGON" and self.polygon_points:
                                self.polygon_points = []
                            elif self.mode == "POLYGON":
                                self.mode = "SELECT"
                            else:
                                self.save_current_sheet()
                                logger.info("Exiting annotator safely.")
                                should_exit = True
                                break

                    elif ev_type == "close":
                        self.save_current_sheet()
                        logger.info("Exiting annotator safely.")
                        should_exit = True
                        break

                if should_exit:
                    break

                time.sleep(0.005)  # Responsive loop tick
        finally:
            win.close()


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments for the SAM 2 botanical annotator."""
    parser = argparse.ArgumentParser(description="SAM 2 Interactive Botanical Annotator for Packera")
    parser.add_argument("--images-dir", type=str, default="data/raw_vouchers", help="Input vouchers directory")
    parser.add_argument("--output-dir", type=str, default="data/raw_annotations", help="Annotations output directory")
    parser.add_argument("--single-image", type=str, default=None, help="Target a specific image file")
    parser.add_argument("--checkpoint", type=str, default="models/checkpoints/sam2_hiera_large.pt", help="SAM2 weights")
    parser.add_argument("--config", type=str, default="sam2_hiera_l.yaml", help="SAM2 model configuration")
    return parser.parse_args()


def main() -> None:
    """Main execution entrypoint."""
    args = parse_args()
    annotator = PrecisionSAM2Annotator(
        images_dir=args.images_dir,
        output_dir=args.output_dir,
        single_image=args.single_image,
        checkpoint_path=args.checkpoint,
        config_path=args.config
    )
    logger.info("Launching Precision SAM 2 Annotator GUI...")
    annotator.run()


if __name__ == "__main__":
    main()
