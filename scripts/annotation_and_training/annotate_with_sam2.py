#!/usr/bin/env python3
"""
===============================================================================
Script: annotate_with_sam2.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Interactive Botanical Instance Segmentation Annotator powered by Segment
    Anything Model 2 (SAM 2) with high-performance native X11 GUI rendering,
    smooth cursor-centered navigation, multimask proposal cycling, high-contrast
    inspection modes, two-click knife slicing, single-pixel morphological margin
    tuning, instant class commitment, and incremental COCO JSON export.
===============================================================================
"""

from __future__ import annotations

import argparse
import ctypes
import json
import logging
import os
import re
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

SAM2_ARCHIVE_PATH = PROJECT_ROOT / "scripts" / "_archive" / "root_artifacts" / "segment-anything-2"
if SAM2_ARCHIVE_PATH.exists() and str(SAM2_ARCHIVE_PATH) not in sys.path:
    sys.path.insert(0, str(SAM2_ARCHIVE_PATH))

try:
    from sam2.build_sam import build_sam2
    from sam2.sam2_image_predictor import SAM2ImagePredictor
except ImportError:
    build_sam2 = None
    SAM2ImagePredictor = None

from scripts.annotation_and_training.sam2_annotator_utils import (
    CLASS_COLORS,
    CLASS_NAMES,
    KEYSYM_TO_CLASS,
    PCD_CLASS_MAPPING,
    PCD_COCO_CATEGORIES,
    apply_knife_cut,
    apply_mask_dilation,
    apply_mask_erosion,
    apply_morphological_tuning,
    apply_viewport_transform,
    calculate_cursor_centered_zoom,
    clamp_viewport_pan,
    clip_box_to_image,
    compose_mask_overlay,
    convert_masks_to_coco_dataset,
    export_coco_annotations,
    image_to_viewport_coords,
    mask_to_normalized_polygon,
    mask_to_polygons,
    mask_to_yolo_bbox,
    overlay_candidate_mask_on_viewport,
    parse_mask_filename,
    polygon_interior_point,
    polygon_to_bounding_box,
    rasterize_lasso_polygon,
    render_hud_overlay,
    get_undo_button_rect,
    get_prev_button_rect,
    get_next_button_rect,
    get_save_button_rect,
    save_coco_json,
    split_mask_with_knife_line,
    viewport_to_image_coords,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("SAM2Annotator")


class XKeyEvent(ctypes.Structure):
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
        ('keycode', ctypes.c_uint),
        ('same_screen', ctypes.c_int),
    ]


class XButtonEvent(ctypes.Structure):
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
        ('button', ctypes.c_uint),
        ('same_screen', ctypes.c_int),
    ]


class XMotionEvent(ctypes.Structure):
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
        ('is_hint', ctypes.c_char),
        ('same_screen', ctypes.c_int),
    ]


class XClientMessageEvent(ctypes.Structure):
    _fields_ = [
        ('type', ctypes.c_int),
        ('serial', ctypes.c_ulong),
        ('send_event', ctypes.c_int),
        ('display', ctypes.c_void_p),
        ('window', ctypes.c_ulong),
        ('message_type', ctypes.c_ulong),
        ('format', ctypes.c_int),
        ('data_l', ctypes.c_long * 5),
    ]


class XConfigureEvent(ctypes.Structure):
    _fields_ = [
        ('type', ctypes.c_int),
        ('serial', ctypes.c_ulong),
        ('send_event', ctypes.c_int),
        ('display', ctypes.c_void_p),
        ('event', ctypes.c_ulong),
        ('window', ctypes.c_ulong),
        ('x', ctypes.c_int),
        ('y', ctypes.c_int),
        ('width', ctypes.c_int),
        ('height', ctypes.c_int),
        ('border_width', ctypes.c_int),
        ('above', ctypes.c_ulong),
        ('override_redirect', ctypes.c_int),
    ]


class XEvent(ctypes.Union):
    _fields_ = [
        ('type', ctypes.c_int),
        ('xkey', XKeyEvent),
        ('xbutton', XButtonEvent),
        ('xmotion', XMotionEvent),
        ('xclient', XClientMessageEvent),
        ('xconfigure', XConfigureEvent),
        ('pad', ctypes.c_long * 24),
    ]


class X11GUIWindow:
    """
    High-performance native X11 window for interactive botanical annotation,
    independent of Qt HighGUI to prevent phantom canvas panning.
    """

    def __init__(self, title: str = "SAM 2 Precision Botanical Annotator", width: int = 1280, height: int = 800):
        self.width = width
        self.height = height
        self.title = title
        self._x11 = ctypes.CDLL("libX11.so.6")
        self.disp = None
        self.win = None
        self.gc = None
        self.ximage = None
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
        self._x11.XPending.restype = ctypes.c_int
        self._x11.XPending.argtypes = [ctypes.c_void_p]
        self._x11.XEventsQueued.restype = ctypes.c_int
        self._x11.XEventsQueued.argtypes = [ctypes.c_void_p, ctypes.c_int]
        self._x11.XNextEvent.restype = ctypes.c_int
        self._x11.XNextEvent.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
        self._x11.XPeekEvent.restype = ctypes.c_int
        self._x11.XPeekEvent.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
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

        # Include ButtonMotionMask, Button1-3MotionMasks so drag motion is reported
        event_mask = (
            (1 << 0)   # KeyPressMask
            | (1 << 1) # KeyReleaseMask
            | (1 << 2) # ButtonPressMask
            | (1 << 3) # ButtonReleaseMask
            | (1 << 6) # PointerMotionMask
            | (1 << 8) # Button1MotionMask
            | (1 << 9) # Button2MotionMask
            | (1 << 10) # Button3MotionMask
            | (1 << 13) # ButtonMotionMask
            | (1 << 15) # ExposureMask
            | (1 << 17) # StructureNotifyMask
        )
        self._x11.XSelectInput(self.disp, self.win, event_mask)
        self._x11.XMapWindow(self.disp, self.win)
        self._x11.XFlush(self.disp)

        self.gc = self._x11.XCreateGC(self.disp, self.win, 0, None)
        self.bgra_buffer = np.zeros((self.height, self.width, 4), dtype=np.uint8)
        self.ximage = self._x11.XCreateImage(
            self.disp, self.visual, self.depth, 2,
            0, self.bgra_buffer.ctypes.data_as(ctypes.c_char_p),
            self.width, self.height, 32, 0
        )

    def _resize_buffer(self, new_width: int, new_height: int) -> None:
        """Reallocates backbuffer and XImage on window resize."""
        if new_width <= 0 or new_height <= 0 or (new_width == self.width and new_height == self.height):
            return
        self.width = new_width
        self.height = new_height
        self.bgra_buffer = np.zeros((self.height, self.width, 4), dtype=np.uint8)
        self.ximage = self._x11.XCreateImage(
            self.disp, self.visual, self.depth, 2,
            0, self.bgra_buffer.ctypes.data_as(ctypes.c_char_p),
            self.width, self.height, 32, 0
        )

    def imshow(self, bgr_img: np.ndarray) -> None:
        """Transfers BGR image directly to native X11 window buffer."""
        h, w = bgr_img.shape[:2]
        if h != self.height or w != self.width:
            bgr_img = cv2.resize(bgr_img, (self.width, self.height))
        self.bgra_buffer[:, :, :3] = bgr_img
        self.bgra_buffer[:, :, 3] = 255
        self._x11.XPutImage(self.disp, self.win, self.gc, self.ximage, 0, 0, 0, 0, self.width, self.height)
        self._x11.XFlush(self.disp)

    def poll_events(self) -> List[Tuple[str, Any]]:
        """Polls queued X11 events with 192-byte safe union, auto-repeat filter, and resize/close tracking."""
        raw_events = []
        if not hasattr(self, "disp") or not self.disp:
            return raw_events

        evt = XEvent()

        while self._x11.XPending(self.disp) > 0:
            self._x11.XNextEvent(self.disp, ctypes.byref(evt))
            if evt.type == 2:  # KeyPress
                keysym = self._x11.XKeycodeToKeysym(self.disp, evt.xkey.keycode, 0)
                raw_events.append(('key_press', keysym, evt.xkey.state))
            elif evt.type == 3:  # KeyRelease
                # Detect and swallow fake X11 keyboard auto-repeat releases
                is_repeat = False
                if self._x11.XPending(self.disp) > 0:
                    next_evt = XEvent()
                    self._x11.XPeekEvent(self.disp, ctypes.byref(next_evt))
                    if (next_evt.type == 2 and
                        next_evt.xkey.keycode == evt.xkey.keycode and
                        next_evt.xkey.time == evt.xkey.time):
                        self._x11.XNextEvent(self.disp, ctypes.byref(next_evt))
                        is_repeat = True
                if not is_repeat:
                    keysym = self._x11.XKeycodeToKeysym(self.disp, evt.xkey.keycode, 0)
                    raw_events.append(('key_release', keysym, evt.xkey.state))
            elif evt.type == 4:  # ButtonPress
                btn = evt.xbutton.button
                raw_events.append(('button_press', evt.xbutton.x, evt.xbutton.y, btn, evt.xbutton.state))
            elif evt.type == 5:  # ButtonRelease
                btn = evt.xbutton.button
                raw_events.append(('button_release', evt.xbutton.x, evt.xbutton.y, btn, evt.xbutton.state))
            elif evt.type == 6:  # MotionNotify
                raw_events.append(('motion', evt.xmotion.x, evt.xmotion.y, evt.xmotion.state))
            elif evt.type == 17:  # DestroyNotify (Window destroyed by WM)
                raw_events.append(('close', None))
            elif evt.type == 22:  # ConfigureNotify (Window Resize)
                cfg = evt.xconfigure
                if cfg.width > 0 and cfg.height > 0 and (cfg.width != self.width or cfg.height != self.height):
                    self._resize_buffer(cfg.width, cfg.height)
                    raw_events.append(('resize', cfg.width, cfg.height))
            elif evt.type == 33:  # ClientMessage (WM_DELETE_WINDOW)
                if evt.xclient.data_l[0] == self.wm_delete or evt.xclient.message_type == self.wm_delete:
                    raw_events.append(('close', None))

        return raw_events

    def close(self) -> None:
        """Safely cleans up X11 window and display connection."""
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
    Interactive botanical annotator using SAM 2 on Packera herbarium vouchers.
    Features cursor-centered smooth zooming, multimask proposal cycling,
    quality inspection modes, botanical knife slicing, morphological tuning,
    instant classification commit, and incremental COCO JSON dataset export.
    """

    def __init__(
        self,
        images_dir: Union[str, Path] = "data/raw_vouchers",
        output_dir: Union[str, Path] = "data/raw_annotations",
        coco_output: Union[str, Path] = "data/annotations/annotations_packera_train.json",
        single_image: Optional[Union[str, Path]] = None,
        checkpoint_path: Union[str, Path] = "models/checkpoints/sam2_hiera_large.pt",
        config_path: Union[str, Path] = "sam2_hiera_l.yaml",
        window_w: int = 1280,
        window_h: int = 800,
        resume_last: bool = False,
        resume_unannotated: bool = False,
        target_voucher: Optional[str] = None,
        target_index: Optional[int] = None,
        tier: str = "1",
        vouchers_csv: Union[str, Path] = "data/tables/curated_vouchers.csv",
    ):
        self.project_root = get_project_root()
        self.images_dir = Path(images_dir)
        self.output_dir = Path(output_dir)
        self.coco_output = Path(coco_output)
        self.masks_dir = self.output_dir / "masks"
        self.labels_dir = self.output_dir / "labels"

        self.masks_dir.mkdir(parents=True, exist_ok=True)
        self.labels_dir.mkdir(parents=True, exist_ok=True)
        self.coco_output.parent.mkdir(parents=True, exist_ok=True)

        self.window_w = window_w
        self.window_h = window_h
        self.resume_last = resume_last
        self.resume_unannotated = resume_unannotated
        self.target_voucher = str(target_voucher).strip() if target_voucher else None
        self.target_index = target_index
        self.tier = str(tier).strip().lower()
        self.vouchers_csv = Path(vouchers_csv)
        self.voucher_tier_map: Dict[str, str] = {}

        if self.vouchers_csv.exists():
            try:
                import pandas as pd
                df_tiers = pd.read_csv(self.vouchers_csv)
                if "catalogNumber" in df_tiers.columns and "determiner_tier" in df_tiers.columns:
                    for _, row in df_tiers.iterrows():
                        cat = str(row["catalogNumber"]).strip()
                        dt = str(row["determiner_tier"]).strip()
                        self.voucher_tier_map[cat] = dt
                        if "image_path" in row and pd.notna(row["image_path"]):
                            stem = Path(str(row["image_path"])).stem
                            self.voucher_tier_map[stem] = dt
            except Exception as e:
                logger.warning(f"Could not load voucher tiers from {self.vouchers_csv}: {e}")

        self.tier_label: Optional[str] = None
        if single_image:
            self.image_files = [Path(single_image)]
        else:
            exts = {".jpg", ".jpeg", ".png", ".tif", ".tiff"}
            all_images = sorted([
                p for p in self.images_dir.glob("*.*")
                if p.suffix.lower() in exts and not p.name.startswith(".")
            ])
            if self.tier in ("1", "tier1", "tier_1", "gold", "tier_1_gold") and self.voucher_tier_map:
                self.tier_label = "Tier 1 Gold"
                filtered = [
                    p for p in all_images
                    if self.voucher_tier_map.get(p.stem) == "Tier_1_Gold"
                ]
                if filtered:
                    self.image_files = filtered
                    logger.info(f"Filtered to {len(self.image_files)} Tier 1 (Gold monograph authority) vouchers.")
                else:
                    logger.warning("No Tier 1 vouchers matched in metadata table; using all images.")
                    self.image_files = all_images
            elif self.tier in ("2", "tier2", "tier_2", "silver") and self.voucher_tier_map:
                self.tier_label = "Tier 2 Silver"
                filtered = [p for p in all_images if self.voucher_tier_map.get(p.stem) == "Tier_2_Silver"]
                self.image_files = filtered if filtered else all_images
            elif self.tier in ("3", "tier3", "tier_3", "bronze") and self.voucher_tier_map:
                self.tier_label = "Tier 3 Bronze"
                filtered = [p for p in all_images if self.voucher_tier_map.get(p.stem) == "Tier_3_Bronze"]
                self.image_files = filtered if filtered else all_images
            else:
                self.tier_label = "All Tiers" if self.tier in ("all", "all_tiers", "0") else None
                self.image_files = all_images

        # If a specific target voucher was requested outside the active tier filter, include it
        if self.target_voucher and not any(p.stem == self.target_voucher or self.target_voucher in p.name for p in self.image_files):
            for p in self.images_dir.glob("*.*"):
                if p.stem == self.target_voucher or self.target_voucher in p.name:
                    self.image_files.insert(0, p)
                    break

        self.checkpoint_path = Path(checkpoint_path)
        self.config_path = str(config_path)
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.predictor = None
        self._init_model()

        self.is_dirty: bool = False
        self.current_idx = 0
        if not single_image and self.image_files:
            if self.target_voucher:
                for idx, p in enumerate(self.image_files):
                    if p.stem == self.target_voucher or self.target_voucher in p.name:
                        self.current_idx = idx
                        logger.info(f"Jumped to specified voucher [{self.current_idx + 1}/{len(self.image_files)}]: {p.name}")
                        break
                else:
                    logger.warning(f"Specified voucher '{self.target_voucher}' not found in {self.images_dir}; starting at index 0.")
            elif self.target_index is not None:
                self.current_idx = max(0, min(len(self.image_files) - 1, int(self.target_index) - 1))
                logger.info(f"Jumped to voucher at index [{self.current_idx + 1}/{len(self.image_files)}]: {self.image_files[self.current_idx].name}")
            elif self.resume_last:
                target_stem = None
                state_file = self.output_dir / ".session_state.json"
                if state_file.exists():
                    try:
                        with open(state_file, "r", encoding="utf-8") as f:
                            state = json.load(f)
                            target_stem = state.get("last_annotated_voucher") or state.get("current_voucher")
                    except Exception as e:
                        logger.warning(f"Failed to read session state: {e}")

                if not target_stem and self.masks_dir.exists():
                    # Fallback: identify voucher with the most recently modified mask file
                    mask_files = list(self.masks_dir.glob("*_inst*.png"))
                    if mask_files:
                        latest_mask = max(mask_files, key=lambda p: p.stat().st_mtime)
                        target_stem = latest_mask.stem.split("_inst")[0]

                if not target_stem:
                    annotated_stems = {
                        p.stem for p in self.output_dir.glob("*.txt")
                        if p.stat().st_size > 0
                    }
                    if annotated_stems:
                        for idx in range(len(self.image_files) - 1, -1, -1):
                            if self.image_files[idx].stem in annotated_stems:
                                target_stem = self.image_files[idx].stem
                                break

                found_idx = None
                if target_stem:
                    for idx, p in enumerate(self.image_files):
                        if p.stem == target_stem:
                            found_idx = idx
                            break

                if found_idx is None:
                    annotated_stems = {
                        p.stem for p in self.output_dir.glob("*.txt")
                        if p.stat().st_size > 0
                    }
                    if self.masks_dir.exists():
                        for m in self.masks_dir.glob("*_inst*.png"):
                            annotated_stems.add(m.stem.split("_inst")[0])
                    for idx in range(len(self.image_files) - 1, -1, -1):
                        if self.image_files[idx].stem in annotated_stems:
                            found_idx = idx
                            break

                if found_idx is not None:
                    self.current_idx = found_idx
                    logger.info(
                        f"Resumed to voucher [{self.current_idx + 1}/{len(self.image_files)}]: {self.image_files[self.current_idx].name}"
                    )
                else:
                    logger.info("No prior annotations found in active tier to resume; starting at index 0.")
            elif self.resume_unannotated:
                for idx, p in enumerate(self.image_files):
                    txt_p = self.output_dir / f"{p.stem}.txt"
                    if not txt_p.exists() or txt_p.stat().st_size == 0:
                        self.current_idx = idx
                        logger.info(
                            f"Resumed to first unannotated voucher [{self.current_idx + 1}/{len(self.image_files)}]: {p.name}"
                        )
                        break
                else:
                    logger.info("All vouchers have annotations; starting at index 0.")

        self.active_image: Optional[np.ndarray] = None
        self.cached_base_layer: Optional[np.ndarray] = None
        self.orig_h = 1000
        self.orig_w = 1000

        # Prompt buffers
        self.point_coords: List[List[float]] = []
        self.point_labels: List[int] = []
        self.box_prompt: Optional[List[float]] = None
        self.polygon_points: List[Tuple[int, int]] = []

        # SAM 2 Multimask candidate state
        self.candidate_masks: List[np.ndarray] = []
        self.candidate_scores: List[float] = []
        self.active_mask_idx: int = 0
        self.candidate_mask: Optional[np.ndarray] = None
        self.saved_instances: List[Dict[str, Any]] = []

        # Viewport Navigation
        self.zoom_level: float = 1.0
        self.pan_offset: List[int] = [0, 0]
        self.mode: str = "SELECT"

        # Quality Inspection Modes
        self.view_mode: str = "FILL"        # "FILL" or "CONTOUR"
        self.is_peeking: bool = False       # Hold-to-peek 'v' key
        self.overlay_alpha: float = 0.55    # Adjust via '[' and ']'

        # Mouse & Key Interaction tracking
        self.lbutton_down: bool = False
        self.rbutton_down: bool = False
        self.space_down: bool = False
        self.drag_start_screen: Tuple[int, int] = (0, 0)
        self.drag_start_img: Tuple[int, int] = (0, 0)
        self.is_box_dragging: bool = False
        self.is_pan_dragging: bool = False
        self.pan_drag_start: Tuple[int, int] = (0, 0)
        self.pan_offset_start: List[int] = [0, 0]

        # Two-Click Knife Tool
        self.knife_pt_a: Optional[Tuple[int, int]] = None
        self.hover_img_pos: Tuple[int, int] = (0, 0)
        self.hover_screen_pos: Tuple[int, int] = (0, 0)

    def _init_model(self) -> None:
        """Initializes SAM 2 model weights and image predictor."""
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
        """Pre-composites saved instances onto base image."""
        if self.active_image is None:
            self.cached_base_layer = None
            return
        if not self.saved_instances:
            self.cached_base_layer = self.active_image.copy()
            return
        self.cached_base_layer = compose_mask_overlay(self.active_image, self.saved_instances, alpha=self.overlay_alpha)

    def run_inference(self) -> None:
        """Executes SAM 2 inference with multimask_output=True to retrieve 3 proposal granularities."""
        pts = np.array(self.point_coords, dtype=np.float32) if self.point_coords else None
        lbls = np.array(self.point_labels, dtype=np.int32) if self.point_labels else None
        box = np.array(self.box_prompt, dtype=np.float32) if self.box_prompt else None

        if pts is None and box is None:
            self.candidate_masks = []
            self.candidate_scores = []
            self.candidate_mask = None
            return

        if self.predictor is None or self.active_image is None:
            return

        try:
            masks, scores, _ = self.predictor.predict(
                point_coords=pts,
                point_labels=lbls,
                box=box,
                multimask_output=True
            )
            if masks is not None and len(masks) > 0:
                self.candidate_masks = [(m > 0.0).astype(np.uint8) * 255 for m in masks]
                self.candidate_scores = [float(s) for s in scores]
                if self.active_mask_idx >= len(self.candidate_masks):
                    self.active_mask_idx = min(1, len(self.candidate_masks) - 1)
                self.candidate_mask = self.candidate_masks[self.active_mask_idx]
        except Exception as e:
            logger.error(f"SAM 2 prediction error: {e}")

    def cycle_multimask_proposal(self) -> None:
        """Cycles through candidate granularities (1: sub-lobe, 2: blade, 3: clump)."""
        if not self.candidate_masks:
            return
        self.active_mask_idx = (self.active_mask_idx + 1) % len(self.candidate_masks)
        self.candidate_mask = self.candidate_masks[self.active_mask_idx]
        iou_val = self.candidate_scores[self.active_mask_idx] if self.active_mask_idx < len(self.candidate_scores) else None
        iou_str = f"{iou_val:.2f}" if iou_val is not None else "N/A"
        logger.info(f"Cycled mask proposal -> #{self.active_mask_idx + 1}/3 (IoU: {iou_str})")

    def finalize_polygon_selection(self) -> None:
        """Converts marked polygon vertices into bounding box and interior point prompt."""
        if len(self.polygon_points) < 3:
            return
        poly_box = polygon_to_bounding_box(self.polygon_points)
        if not poly_box:
            return

        bx0, by0, bx1, by1 = poly_box
        if abs(bx1 - bx0) < 3 or abs(by1 - by0) < 3:
            return

        self.box_prompt = [float(bx0), float(by0), float(bx1), float(by1)]
        interior_pt = polygon_interior_point(self.polygon_points, self.orig_h, self.orig_w)
        if interior_pt:
            self.point_coords = [[interior_pt[0], interior_pt[1]]]
            self.point_labels = [1]
        else:
            self.point_coords = []
            self.point_labels = []

        self.run_inference()
        poly_mask = rasterize_lasso_polygon(self.polygon_points, self.orig_h, self.orig_w)
        if self.candidate_mask is not None and np.count_nonzero(self.candidate_mask) > 0:
            constrained = cv2.bitwise_and(self.candidate_mask, poly_mask)
            self.candidate_mask = constrained if np.count_nonzero(constrained) > 0 else poly_mask
        else:
            self.candidate_mask = poly_mask

        if self.candidate_masks and self.active_mask_idx < len(self.candidate_masks):
            self.candidate_masks[self.active_mask_idx] = self.candidate_mask.copy()

        logger.info(f"Finalized custom polygon bounding box with {len(self.polygon_points)} vertices")
        self.polygon_points = []

    def undo_last_point(self) -> bool:
        """
        Removes the most recently placed inclusion or exclusion point prompt,
        or cancels uncommitted knife / polygon / bounding box constraints,
        and re-executes SAM 2 inference in real time.
        """
        if self.mode == "KNIFE" and self.knife_pt_a is not None:
            self.knife_pt_a = None
            logger.info("Undid Knife Point A.")
            return True

        if self.mode == "POLYGON" and self.polygon_points:
            popped_poly = self.polygon_points.pop()
            logger.info(f"Undid last polygon vertex at {popped_poly}. Remaining vertices: {len(self.polygon_points)}")
            return True

        if self.point_coords:
            popped_pt = self.point_coords.pop()
            popped_lbl = self.point_labels.pop() if self.point_labels else 1
            lbl_name = "inclusion (+)" if popped_lbl == 1 else "exclusion (-)"
            logger.info(
                f"Undid last prompt point: {lbl_name} at ({popped_pt[0]:.1f}, {popped_pt[1]:.1f}). "
                f"Remaining points: {len(self.point_coords)}"
            )
            self.run_inference()
            return True

        if self.box_prompt is not None:
            self.box_prompt = None
            logger.info("Undid bounding box prompt.")
            self.run_inference()
            return True

        logger.info("No active prompt points or constraints to undo.")
        return False

    def commit_active_instance(self, class_id: int) -> None:
        """Instant Class Commit: assigns class, commits mask, resets prompts, readies next leaf."""
        if class_id < 0 or class_id >= len(CLASS_NAMES):
            return
        if self.candidate_mask is None or np.count_nonzero(self.candidate_mask) == 0:
            return

        label = CLASS_NAMES[class_id]
        poly_str = mask_to_normalized_polygon(self.candidate_mask, class_id)
        if poly_str:
            poly = [float(val) for val in poly_str.split()[1:]]
            self.saved_instances.append({
                "class_id": class_id,
                "label": label,
                "polygon": poly,
                "binary_mask": self.candidate_mask > 0,
                "mask": self.candidate_mask.copy()
            })
            self.candidate_mask = None
            self.candidate_masks = []
            self.candidate_scores = []
            self.active_mask_idx = 0
            self.point_coords = []
            self.point_labels = []
            self.box_prompt = None
            self.polygon_points = []
            self.knife_pt_a = None
            self.is_dirty = True
            self._update_cached_base_layer()
            logger.info(f"Committed instance #{len(self.saved_instances)} -> '{label}' (Class {class_id})")

    def _save_session_state(self) -> None:
        """Saves current navigation position to session state file."""
        if not hasattr(self, "image_files") or not self.image_files or self.current_idx >= len(self.image_files):
            return
        state_file = self.output_dir / ".session_state.json"
        try:
            current_voucher = self.image_files[self.current_idx].stem
            state = {}
            if state_file.exists():
                with open(state_file, "r", encoding="utf-8") as f:
                    state = json.load(f)
            state["current_voucher"] = current_voucher
            state["current_index"] = self.current_idx
            state["last_active_time"] = time.time()
            with open(state_file, "w", encoding="utf-8") as f:
                json.dump(state, f, indent=2)
        except Exception as e:
            logger.debug(f"Could not write session state: {e}")

    def _save_annotated_state(self, voucher_id: str) -> None:
        """Records the voucher that was just annotated/saved."""
        state_file = self.output_dir / ".session_state.json"
        try:
            state = {}
            if state_file.exists():
                with open(state_file, "r", encoding="utf-8") as f:
                    state = json.load(f)
            state["last_annotated_voucher"] = voucher_id
            state["last_annotated_index"] = self.current_idx
            state["current_voucher"] = voucher_id
            state["current_index"] = self.current_idx
            state["last_annotated_time"] = time.time()
            with open(state_file, "w", encoding="utf-8") as f:
                json.dump(state, f, indent=2)
        except Exception as e:
            logger.debug(f"Could not record annotated state: {e}")

    def prev_voucher(self) -> None:
        """Navigates to the previous voucher sheet WITHOUT saving or modifying annotations."""
        if self.current_idx > 0:
            self.current_idx -= 1
            self.load_active_image()
            self._save_session_state()
            logger.info(f"Navigated to previous voucher [{self.current_idx + 1}/{len(self.image_files)}]: {self.image_files[self.current_idx].name}")
        else:
            logger.info("Already at first voucher sheet.")

    def next_voucher(self) -> None:
        """Navigates to the next voucher sheet WITHOUT saving or modifying annotations."""
        if self.current_idx + 1 < len(self.image_files):
            self.current_idx += 1
            self.load_active_image()
            self._save_session_state()
            logger.info(f"Navigated to next voucher [{self.current_idx + 1}/{len(self.image_files)}]: {self.image_files[self.current_idx].name}")
        else:
            logger.info("Already at last voucher sheet.")

    def save_current_sheet(self, sync_coco: bool = True) -> None:
        """Autosaves active sheet's annotations to YOLO txt, PNG masks, and optionally syncs COCO JSON."""
        if not hasattr(self, "image_files") or not self.image_files or self.current_idx >= len(self.image_files):
            return

        current_file = self.image_files[self.current_idx]
        voucher_id = current_file.stem
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.masks_dir.mkdir(parents=True, exist_ok=True)

        txt_file = self.output_dir / f"{voucher_id}.txt"
        if not self.saved_instances:
            if txt_file.exists() and txt_file.stat().st_size > 0:
                logger.info(f"No new instances committed for voucher {voucher_id}; preserving existing annotations.")
                return

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

        # Remove lingering mask files if instances were deleted/undone
        for old_mask in self.masks_dir.glob(f"{voucher_id}_inst*.png"):
            m = re.match(rf"^{re.escape(voucher_id)}_inst(\d+)_.+\.png$", old_mask.name)
            if m and int(m.group(1)) >= len(self.saved_instances):
                old_mask.unlink(missing_ok=True)

        self.is_dirty = False
        self._save_annotated_state(voucher_id)
        logger.info(f"Saved {len(self.saved_instances)} instances for voucher {voucher_id}")

        if sync_coco:
            try:
                export_coco_annotations(
                    masks_dir=self.masks_dir,
                    output_coco_path=self.coco_output,
                    images_dir=self.images_dir,
                )
                logger.info(f"Synchronized COCO dataset to {self.coco_output}")
            except Exception as err:
                logger.warning(f"COCO synchronization warning: {err}")

    def _load_existing_instances(self, voucher_id: str) -> List[Dict[str, Any]]:
        """
        Loads pre-existing annotations for the given voucher sheet from raw_annotations (.txt and masks).
        Populates self.saved_instances so prior annotations render directly in the GUI.
        """
        instances: List[Dict[str, Any]] = []
        txt_path = self.output_dir / f"{voucher_id}.txt"

        # 1. Primary path: load from YOLO polygon text file if present
        if txt_path.exists() and txt_path.stat().st_size > 0:
            try:
                with open(txt_path, "r", encoding="utf-8") as f:
                    for idx, line in enumerate(f):
                        parts = line.strip().split()
                        if len(parts) < 7:
                            continue
                        try:
                            class_id = int(parts[0])
                            poly = [float(v) for v in parts[1:]]
                        except ValueError:
                            continue

                        label = CLASS_NAMES[class_id] if 0 <= class_id < len(CLASS_NAMES) else f"class_{class_id}"

                        mask = None
                        expected_mask = self.masks_dir / f"{voucher_id}_inst{idx:02d}_{label}.png"
                        if expected_mask.exists():
                            mask = cv2.imread(str(expected_mask), cv2.IMREAD_GRAYSCALE)
                        else:
                            matches = list(self.masks_dir.glob(f"{voucher_id}_inst{idx:02d}_*.png"))
                            if matches:
                                mask = cv2.imread(str(matches[0]), cv2.IMREAD_GRAYSCALE)

                        if mask is not None:
                            if mask.shape[:2] != (self.orig_h, self.orig_w):
                                mask = cv2.resize(mask, (self.orig_w, self.orig_h), interpolation=cv2.INTER_NEAREST)
                        else:
                            mask = np.zeros((self.orig_h, self.orig_w), dtype=np.uint8)
                            pts = np.array(poly, dtype=np.float32).reshape(-1, 2)
                            pts[:, 0] *= self.orig_w
                            pts[:, 1] *= self.orig_h
                            cv2.fillPoly(mask, [np.int32(pts)], 255)

                        instances.append({
                            "class_id": class_id,
                            "label": label,
                            "polygon": poly,
                            "binary_mask": mask > 0,
                            "mask": mask,
                        })
            except Exception as e:
                logger.warning(f"Failed to load annotations from {txt_path}: {e}")

        # 2. Fallback path: if no txt file (or empty), check if masks exist on disk
        if not instances and self.masks_dir.exists():
            mask_files = sorted(self.masks_dir.glob(f"{voucher_id}_inst*.png"))
            for idx, mask_path in enumerate(mask_files):
                _, label, _ = parse_mask_filename(mask_path.name)
                class_id = CLASS_NAMES.index(label) if label in CLASS_NAMES else 0
                mask = cv2.imread(str(mask_path), cv2.IMREAD_GRAYSCALE)
                if mask is None:
                    continue
                if mask.shape[:2] != (self.orig_h, self.orig_w):
                    mask = cv2.resize(mask, (self.orig_w, self.orig_h), interpolation=cv2.INTER_NEAREST)

                poly_str = mask_to_normalized_polygon(mask, class_id)
                poly = [float(val) for val in poly_str.split()[1:]] if poly_str else []

                instances.append({
                    "class_id": class_id,
                    "label": label,
                    "polygon": poly,
                    "binary_mask": mask > 0,
                    "mask": mask,
                })

        return instances

    def load_active_image(self) -> bool:
        """Loads current voucher sheet into memory and initializes SAM 2 image embeddings."""
        if not hasattr(self, "image_files") or not self.image_files or self.current_idx >= len(self.image_files):
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
        self.candidate_masks = []
        self.candidate_scores = []
        self.active_mask_idx = 0
        self.candidate_mask = None
        self.zoom_level = 1.0
        self.pan_offset = [0, 0]
        self.polygon_points = []
        self.knife_pt_a = None
        self.is_dirty = False

        # Repopulate previously completed annotations
        voucher_id = img_path.stem
        self.saved_instances = self._load_existing_instances(voucher_id)
        if self.saved_instances:
            logger.info(f"Loaded {len(self.saved_instances)} existing instance annotations for voucher {voucher_id}")

        self._update_cached_base_layer()

        if self.predictor is not None:
            img_rgb = cv2.cvtColor(self.active_image, cv2.COLOR_BGR2RGB)
            self.predictor.set_image(img_rgb)
            logger.info(f"Loaded voucher [{self.current_idx + 1}/{len(self.image_files)}]: {img_path.name}")
        return True

    def _viewport_to_image(
        self, vx: int, vy: int, scale_x: float, scale_y: float, crop_x0: int, crop_y0: int
    ) -> Tuple[int, int]:
        ix = int(crop_x0 + vx / max(scale_x, 1e-6))
        iy = int(crop_y0 + vy / max(scale_y, 1e-6))
        return max(0, min(self.orig_w - 1, ix)), max(0, min(self.orig_h - 1, iy))

    def _image_to_viewport(
        self, ix: Union[int, float], iy: Union[int, float], scale_x: float, scale_y: float, crop_x0: int, crop_y0: int
    ) -> Tuple[int, int]:
        return int((ix - crop_x0) * scale_x), int((iy - crop_y0) * scale_y)

    def run(self) -> None:
        """Main interactive GUI annotation event loop."""
        if not hasattr(self, "image_files") or not self.image_files:
            logger.warning(f"No valid voucher images found in {self.images_dir}")
            return

        win = X11GUIWindow(title="SAM 2 Precision Botanical Annotator", width=self.window_w, height=self.window_h)
        if not self.load_active_image():
            logger.error("Failed to load initial voucher image.")
            win.close()
            return

        try:
            while True:
                if self.cached_base_layer is None or self.active_image is None:
                    break

                current_voucher = self.image_files[self.current_idx].stem
                render_base = self.active_image if self.is_peeking else self.cached_base_layer

                viewport_img, transform = apply_viewport_transform(
                    render_base,
                    self.zoom_level,
                    tuple(self.pan_offset),
                    self.window_w,
                    self.window_h,
                )
                scale_x, scale_y, crop_x0, crop_y0 = transform

                # Candidate mask overlay (suppressed while holding 'v' to peek)
                if not self.is_peeking and self.candidate_mask is not None:
                    viewport_img = overlay_candidate_mask_on_viewport(
                        viewport_img,
                        self.candidate_mask,
                        transform,
                        alpha=self.overlay_alpha,
                        view_mode=self.view_mode,
                    )

                # Overlays & Prompts (suppressed while holding 'v')
                if not self.is_peeking:
                    # Point prompts
                    for pt, lbl in zip(self.point_coords, self.point_labels):
                        vx, vy = self._image_to_viewport(pt[0], pt[1], scale_x, scale_y, crop_x0, crop_y0)
                        if 0 <= vx < self.window_w and 0 <= vy < self.window_h:
                            c = (0, 255, 0) if lbl == 1 else (0, 0, 255)
                            cv2.circle(viewport_img, (vx, vy), 5, c, -1)
                            cv2.circle(viewport_img, (vx, vy), 7, (255, 255, 255), 1)

                    # Bounding box prompt
                    if self.box_prompt is not None:
                        bx0, by0, bx1, by1 = self.box_prompt
                        b_vx0, b_vy0 = self._image_to_viewport(bx0, by0, scale_x, scale_y, crop_x0, crop_y0)
                        b_vx1, b_vy1 = self._image_to_viewport(bx1, by1, scale_x, scale_y, crop_x0, crop_y0)
                        cv2.rectangle(viewport_img, (b_vx0, b_vy0), (b_vx1, b_vy1), (0, 255, 255), 2)

                    # Two-Click Knife line preview
                    if self.mode == "KNIFE" and self.knife_pt_a is not None:
                        k0 = self._image_to_viewport(self.knife_pt_a[0], self.knife_pt_a[1], scale_x, scale_y, crop_x0, crop_y0)
                        k1 = self._image_to_viewport(self.hover_img_pos[0], self.hover_img_pos[1], scale_x, scale_y, crop_x0, crop_y0)
                        cv2.line(viewport_img, k0, k1, (0, 0, 255), 2, cv2.LINE_AA)
                        cv2.circle(viewport_img, k0, 4, (0, 255, 255), -1)

                    # Polygon lasso preview
                    if self.polygon_points:
                        v_pts = [self._image_to_viewport(p[0], p[1], scale_x, scale_y, crop_x0, crop_y0) for p in self.polygon_points]
                        if len(v_pts) >= 3:
                            poly_ov = viewport_img.copy()
                            pts_arr = np.array(v_pts, dtype=np.int32).reshape((-1, 1, 2))
                            cv2.fillPoly(poly_ov, [pts_arr], (0, 220, 255))
                            cv2.addWeighted(poly_ov, 0.25, viewport_img, 0.75, 0, viewport_img)

                        for i in range(len(v_pts) - 1):
                            cv2.line(viewport_img, v_pts[i], v_pts[i + 1], (0, 240, 255), 2, cv2.LINE_AA)

                        if self.mode == "POLYGON" and hasattr(self, "hover_screen_pos"):
                            cv2.line(viewport_img, v_pts[-1], self.hover_screen_pos, (0, 200, 255), 1, cv2.LINE_AA)

                        for idx, (vx, vy) in enumerate(v_pts):
                            if 0 <= vx < self.window_w and 0 <= vy < self.window_h:
                                c = (0, 255, 0) if idx == 0 else (0, 220, 255)
                                cv2.circle(viewport_img, (vx, vy), 6 if idx == 0 else 5, c, -1)
                                cv2.circle(viewport_img, (vx, vy), 8 if idx == 0 else 7, (255, 255, 255), 1)

                # Render Top Status HUD banner
                active_iou = (
                    self.candidate_scores[self.active_mask_idx]
                    if (self.candidate_mask is not None and self.candidate_scores and self.active_mask_idx < len(self.candidate_scores))
                    else None
                )
                hud_display = render_hud_overlay(
                    display_img=viewport_img,
                    voucher_name=current_voucher,
                    voucher_idx=self.current_idx,
                    total_vouchers=len(self.image_files),
                    saved_instances=self.saved_instances,
                    mode=self.mode,
                    zoom_level=self.zoom_level,
                    pan_offset=tuple(self.pan_offset),
                    candidate_idx=self.active_mask_idx if self.candidate_mask is not None else None,
                    candidate_total=len(self.candidate_masks) if self.candidate_masks else 3,
                    candidate_iou=active_iou,
                    view_mode=self.view_mode,
                    alpha=self.overlay_alpha,
                    is_dirty=self.is_dirty,
                    tier_label=self.tier_label,
                )
                win.imshow(hud_display)

                should_exit = False
                events = win.poll_events()

                for ev in events:
                    ev_type = ev[0]

                    if ev_type == "button_press":
                        _, vx, vy, btn, state = ev
                        ix, iy = self._viewport_to_image(vx, vy, scale_x, scale_y, crop_x0, crop_y0)
                        self.press_mode = self.mode

                        if btn == 1:  # Left Button
                            if vy < 70:
                                # 1. Check Undo Point
                                bx0, by0, bx1, by1 = get_undo_button_rect(self.window_w)
                                if bx0 <= vx <= bx1 and by0 <= vy <= by1:
                                    self.clicked_undo_button = True
                                    self.undo_last_point()
                                    continue

                                # 2. Check Previous Voucher (< Prev)
                                px0, py0, px1, py1 = get_prev_button_rect(self.window_w)
                                if px0 <= vx <= px1 and py0 <= vy <= py1:
                                    self.prev_voucher()
                                    continue

                                # 3. Check Next Voucher (Next >)
                                nx0, ny0, nx1, ny1 = get_next_button_rect(self.window_w)
                                if nx0 <= vx <= nx1 and ny0 <= vy <= ny1:
                                    self.next_voucher()
                                    continue

                                # 4. Check Save Voucher (Save)
                                sx0, sy0, sx1, sy1 = get_save_button_rect(self.window_w)
                                if sx0 <= vx <= sx1 and sy0 <= vy <= sy1:
                                    self.save_current_sheet(sync_coco=False)
                                    continue
                            self.lbutton_down = True
                            self.drag_start_screen = (vx, vy)
                            self.drag_start_img = (ix, iy)
                            if self.space_down:
                                self.is_pan_dragging = True
                                self.pan_drag_start = (vx, vy)
                                self.pan_offset_start = list(self.pan_offset)
                            elif self.mode == "POLYGON":
                                if len(self.polygon_points) >= 3:
                                    s_vx, s_vy = self._image_to_viewport(
                                        self.polygon_points[0][0], self.polygon_points[0][1],
                                        scale_x, scale_y, crop_x0, crop_y0
                                    )
                                    if max(abs(vx - s_vx), abs(vy - s_vy)) <= 15:
                                        self.finalize_polygon_selection()
                                    else:
                                        self.polygon_points.append((ix, iy))
                                else:
                                    self.polygon_points.append((ix, iy))
                            elif self.mode == "KNIFE":
                                self.knife_click_handled = True
                                self.lbutton_down = False
                                if self.knife_pt_a is None:
                                    self.knife_pt_a = (ix, iy)
                                    logger.info(f"Knife Point A set at {(ix, iy)}. Click Point B to sever mask.")
                                else:
                                    if self.candidate_mask is not None:
                                        # Gather foreground prompts to automatically retain target component
                                        keep_pts = [
                                            (float(pt[0]), float(pt[1]))
                                            for pt, lbl in zip(self.point_coords, self.point_labels)
                                            if lbl == 1
                                        ]
                                        bbox = tuple(self.box_prompt) if self.box_prompt is not None else None

                                        # Fallback to interior point if polygon lasso was used
                                        if not keep_pts and self.polygon_points:
                                            poly_pt = polygon_interior_point(self.polygon_points, self.orig_h, self.orig_w)
                                            if poly_pt is not None:
                                                keep_pts = [poly_pt]

                                        self.candidate_mask = split_mask_with_knife_line(
                                            self.candidate_mask,
                                            self.knife_pt_a,
                                            (ix, iy),
                                            line_thickness=2,
                                            dilation_px=2,
                                            keep_points=keep_pts if (keep_pts or bbox is not None) else None,
                                            bounding_box=bbox,
                                        )
                                        if self.candidate_masks and self.active_mask_idx < len(self.candidate_masks):
                                            self.candidate_masks[self.active_mask_idx] = self.candidate_mask.copy()
                                        logger.info(
                                            f"Applied 2-click knife cut from {self.knife_pt_a} to {(ix, iy)} "
                                            f"(retained target prompt component, pruned severed fragment)."
                                        )
                                    self.knife_pt_a = None
                                    self.mode = "SELECT"
                                    logger.info("Knife cut complete -> switched back to SELECT mode.")
                            else:
                                self.is_box_dragging = False

                        elif btn == 2:  # Middle Button Drag Pan
                            self.is_pan_dragging = True
                            self.pan_drag_start = (vx, vy)
                            self.pan_offset_start = list(self.pan_offset)

                        elif btn == 3:  # Right Button (Drag to pan or click for background prompt)
                            self.rbutton_down = True
                            self.drag_start_screen = (vx, vy)
                            self.drag_start_img = (ix, iy)
                            self.pan_drag_start = (vx, vy)
                            self.pan_offset_start = list(self.pan_offset)

                        elif btn == 4:  # Wheel Up -> Cursor-Centered Zoom In
                            self.zoom_level, self.pan_offset = calculate_cursor_centered_zoom(
                                self.zoom_level, tuple(self.pan_offset), vx, vy, zoom_in=True,
                                target_w=self.window_w, target_h=self.window_h,
                                orig_w=self.orig_w, orig_h=self.orig_h, step=1.15, min_zoom=1.0, max_zoom=16.0
                            )
                            crop_w = max(10, int(self.orig_w / max(self.zoom_level, 1.0)))
                            crop_h = max(10, int(self.orig_h / max(self.zoom_level, 1.0)))
                            scale_x = self.window_w / max(crop_w, 1)
                            scale_y = self.window_h / max(crop_h, 1)
                            crop_x0 = self.pan_offset[0]
                            crop_y0 = self.pan_offset[1]

                        elif btn == 5:  # Wheel Down -> Cursor-Centered Zoom Out
                            self.zoom_level, self.pan_offset = calculate_cursor_centered_zoom(
                                self.zoom_level, tuple(self.pan_offset), vx, vy, zoom_in=False,
                                target_w=self.window_w, target_h=self.window_h,
                                orig_w=self.orig_w, orig_h=self.orig_h, step=1.15, min_zoom=1.0, max_zoom=16.0
                            )
                            crop_w = max(10, int(self.orig_w / max(self.zoom_level, 1.0)))
                            crop_h = max(10, int(self.orig_h / max(self.zoom_level, 1.0)))
                            scale_x = self.window_w / max(crop_w, 1)
                            scale_y = self.window_h / max(crop_h, 1)
                            crop_x0 = self.pan_offset[0]
                            crop_y0 = self.pan_offset[1]

                    elif ev_type == "motion":
                        _, vx, vy, state = ev
                        ix, iy = self._viewport_to_image(vx, vy, scale_x, scale_y, crop_x0, crop_y0)
                        self.hover_img_pos = (ix, iy)
                        self.hover_screen_pos = (vx, vy)

                        if getattr(self, "rbutton_down", False) and not self.is_pan_dragging:
                            if max(abs(vx - self.drag_start_screen[0]), abs(vy - self.drag_start_screen[1])) > 5:
                                self.is_pan_dragging = True

                        if self.is_pan_dragging:
                            dx = vx - self.pan_drag_start[0]
                            dy = vy - self.pan_drag_start[1]
                            img_dx = int(dx / max(scale_x, 1e-6))
                            img_dy = int(dy / max(scale_y, 1e-6))
                            crop_w = int(self.orig_w / max(self.zoom_level, 1.0))
                            crop_h = int(self.orig_h / max(self.zoom_level, 1.0))
                            min_x = -int(crop_w * 0.85)
                            max_x = int(self.orig_w - crop_w * 0.15)
                            min_y = -int(crop_h * 0.85)
                            max_y = int(self.orig_h - crop_h * 0.15)
                            self.pan_offset[0] = max(min_x, min(max_x, self.pan_offset_start[0] - img_dx))
                            self.pan_offset[1] = max(min_y, min(max_y, self.pan_offset_start[1] - img_dy))
                        elif getattr(self, "lbutton_down", False) and self.mode == "SELECT" and not self.space_down:
                            if max(abs(vx - self.drag_start_screen[0]), abs(vy - self.drag_start_screen[1])) > 5:
                                self.is_box_dragging = True
                                x0, y0 = self.drag_start_img
                                self.box_prompt = [float(min(x0, ix)), float(min(y0, iy)), float(max(x0, ix)), float(max(y0, iy))]

                    elif ev_type == "button_release":
                        _, vx, vy, btn, state = ev
                        ix, iy = self._viewport_to_image(vx, vy, scale_x, scale_y, crop_x0, crop_y0)

                        if btn == 1:
                            if getattr(self, "clicked_undo_button", False):
                                self.clicked_undo_button = False
                                continue
                            if getattr(self, "knife_click_handled", False):
                                self.knife_click_handled = False
                                self.lbutton_down = False
                                continue
                            self.lbutton_down = False
                            if self.is_pan_dragging:
                                self.is_pan_dragging = False
                            elif self.is_box_dragging:
                                self.is_box_dragging = False
                                self.run_inference()
                            elif self.mode == "SELECT" and getattr(self, "press_mode", "SELECT") == "SELECT":
                                self.point_coords.append([float(ix), float(iy)])
                                self.point_labels.append(1)
                                self.run_inference()
                        elif btn == 2:
                            self.is_pan_dragging = False
                        elif btn == 3:
                            was_panning = self.is_pan_dragging
                            self.rbutton_down = False
                            self.is_pan_dragging = False
                            if not was_panning:
                                if self.mode == "POLYGON" and len(self.polygon_points) >= 3:
                                    self.finalize_polygon_selection()
                                elif self.mode == "KNIFE":
                                    self.knife_pt_a = None
                                else:
                                    self.point_coords.append([float(ix), float(iy)])
                                    self.point_labels.append(0)
                                    self.run_inference()

                    elif ev_type == "resize":
                        _, new_w, new_h = ev
                        self.window_w = new_w
                        self.window_h = new_h

                    elif ev_type == "key_press":
                        sym = ev[1]
                        state = ev[2] if len(ev) > 2 else 0

                        # Undo prompt point / constraint (Ctrl+Z, 'z'/'Z', Backspace, Delete, Keypad Del/Decimal)
                        if sym in (ord('z'), ord('Z'), 0xff08, 0xffff, 0xffae, 0xff9f):
                            self.undo_last_point()

                        # Instant Class Commit (0 through 6, top-row or numpad)
                        elif sym in KEYSYM_TO_CLASS:
                            class_id = KEYSYM_TO_CLASS[sym]
                            if 0 <= class_id < len(CLASS_NAMES):
                                self.commit_active_instance(class_id)

                        # Multimask proposal cycling (Tab key or Keypad /)
                        elif sym in (0xff09, 0xffaf):
                            self.cycle_multimask_proposal()

                        # View mode toggle (Fill vs Contour) ('o' / 'O' or Keypad *)
                        elif sym in (ord('o'), ord('O'), 0xffaa):
                            self.view_mode = "CONTOUR" if self.view_mode == "FILL" else "FILL"
                            logger.info(f"Toggled view mode -> {self.view_mode}")

                        # Hold-to-peek ('v' / 'V')
                        elif sym in (ord('v'), ord('V')):
                            self.is_peeking = True

                        # Viewport Fit to Window ('f' / 'F')
                        elif sym in (ord('f'), ord('F')):
                            self.zoom_level = 1.0
                            self.pan_offset = [0, 0]
                            logger.info("Reset viewport to 1.0x (Fit to Window)")

                        # Overlay alpha adjustment ('[' and ']')
                        elif sym == ord('['):
                            self.overlay_alpha = max(0.1, round(self.overlay_alpha - 0.1, 2))
                            self._update_cached_base_layer()
                            logger.info(f"Adjusted overlay alpha -> {self.overlay_alpha:.2f}")
                        elif sym == ord(']'):
                            self.overlay_alpha = min(0.9, round(self.overlay_alpha + 0.1, 2))
                            self._update_cached_base_layer()
                            logger.info(f"Adjusted overlay alpha -> {self.overlay_alpha:.2f}")

                        # Morphological margin tuning (+ / = and - / _)
                        elif sym in (ord('+'), ord('='), 0xffab):
                            if self.candidate_mask is not None:
                                self.candidate_mask = apply_morphological_tuning(self.candidate_mask, "dilate", 3)
                                if self.candidate_masks and self.active_mask_idx < len(self.candidate_masks):
                                    self.candidate_masks[self.active_mask_idx] = self.candidate_mask.copy()
                                logger.info("Applied 3x3 binary dilation (+ margin)")
                        elif sym in (ord('-'), ord('_'), 0xffad):
                            if self.candidate_mask is not None:
                                self.candidate_mask = apply_morphological_tuning(self.candidate_mask, "erode", 3)
                                if self.candidate_masks and self.active_mask_idx < len(self.candidate_masks):
                                    self.candidate_masks[self.active_mask_idx] = self.candidate_mask.copy()
                                logger.info("Applied 3x3 binary erosion (- margin)")

                        # Two-click knife tool ('k' / 'K')
                        elif sym in (ord('k'), ord('K')):
                            self.mode = "KNIFE" if self.mode != "KNIFE" else "SELECT"
                            self.knife_pt_a = None
                            logger.info(f"Switched tool mode -> {self.mode}")

                        # Polygon lasso tool ('p' / 'P')
                        elif sym in (ord('p'), ord('P')):
                            self.mode = "POLYGON" if self.mode != "POLYGON" else "SELECT"
                            if self.mode != "POLYGON":
                                self.polygon_points = []
                            logger.info(f"Switched tool mode -> {self.mode}")

                        # Space bar down for panning
                        elif sym == 0x0020:
                            self.space_down = True

                        # Enter to finalize polygon OR save current sheet
                        elif sym in (0xff0d, 0xff8d):
                            if self.mode == "POLYGON" and len(self.polygon_points) >= 3:
                                self.finalize_polygon_selection()
                            else:
                                self.save_current_sheet(sync_coco=False)

                        # Explicit Save key ('s' / 'S')
                        elif sym in (ord('s'), ord('S')):
                            self.save_current_sheet(sync_coco=False)

                        # Clear current prompts ('c' / 'C')
                        elif sym in (ord('c'), ord('C')):
                            self.point_coords = []
                            self.point_labels = []
                            self.box_prompt = None
                            self.candidate_mask = None
                            self.candidate_masks = []
                            self.candidate_scores = []
                            self.polygon_points = []
                            self.knife_pt_a = None

                        # Undo last instance ('u' / 'U')
                        elif sym in (ord('u'), ord('U')):
                            if self.saved_instances:
                                popped = self.saved_instances.pop()
                                self.is_dirty = True
                                self._update_cached_base_layer()
                                logger.info(f"Removed instance: {popped.get('label')}")

                        # Next sheet without saving ('n' / 'N' or Right Arrow)
                        elif sym in (ord('n'), ord('N'), 0xff53):
                            self.next_voucher()

                        # Previous sheet without saving ('b' / 'B' or Left Arrow)
                        elif sym in (ord('b'), ord('B'), 0xff51):
                            self.prev_voucher()

                        # Quit ('q' / 'Q' / Escape)
                        elif sym in (ord('q'), ord('Q'), 0xff1b):
                            if self.is_dirty:
                                self.save_current_sheet(sync_coco=True)
                            should_exit = True
                            break

                    elif ev_type == "key_release":
                        sym = ev[1]
                        if sym in (ord('v'), ord('V')):
                            self.is_peeking = False
                        elif sym == 0x0020:
                            self.space_down = False
                            if self.is_pan_dragging:
                                self.is_pan_dragging = False

                    elif ev_type == "close":
                        if self.is_dirty:
                            logger.info("Window close requested (GUI 'X' button) -> saving modified sheet before exit...")
                            self.save_current_sheet(sync_coco=True)
                        else:
                            logger.info("Window close requested (GUI 'X' button) -> closing cleanly without modifying annotations.")
                        should_exit = True
                        break

                if should_exit:
                    break

                time.sleep(0.005)
        finally:
            win.close()


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments for the SAM 2 botanical annotator."""
    parser = argparse.ArgumentParser(description="SAM 2 Interactive Botanical Annotator for Packera")
    parser.add_argument("--images-dir", type=str, default="data/raw_vouchers", help="Input raw vouchers directory")
    parser.add_argument("--output-dir", type=str, default="data/raw_annotations", help="Annotations output directory")
    parser.add_argument("--output-coco", type=str, default="data/annotations/annotations_packera_train.json", help="Path for standardized COCO JSON")
    parser.add_argument("--single-image", type=str, default=None, help="Target a specific image file")
    parser.add_argument("--checkpoint", type=str, default="models/checkpoints/sam2_hiera_large.pt", help="SAM 2 weights path")
    parser.add_argument("--config", type=str, default="sam2_hiera_l.yaml", help="SAM 2 model configuration")
    parser.add_argument("--export-coco", action="store_true", help="Compile and export COCO dataset from existing masks in headless mode")
    parser.add_argument("--masks-dir", type=str, default="data/raw_annotations/masks", help="Masks directory for --export-coco")
    parser.add_argument("--val-split", type=float, default=0.0, help="Validation partition fraction (e.g. 0.20)")
    parser.add_argument("--resume-last", action="store_true", help="Automatically resume from the most recently annotated voucher")
    parser.add_argument("--resume-unannotated", action="store_true", help="Automatically jump to the first unannotated voucher in sequence")
    parser.add_argument("--voucher", type=str, default=None, help="Jump directly to a specific voucher ID (e.g. 1048)")
    parser.add_argument("--index", type=int, default=None, help="Jump directly to a 1-based voucher index (e.g. 5)")
    parser.add_argument("--tier", type=str, default="1", choices=["1", "2", "3", "all"], help="Filter vouchers by determiner authority tier (default: '1' for Tier 1 Gold monograph authorities; 'all' for all 5347 vouchers)")
    parser.add_argument("--vouchers-csv", type=str, default="data/tables/curated_vouchers.csv", help="Path to curated vouchers metadata CSV")
    return parser.parse_args()


def main() -> None:
    """Main execution entrypoint."""
    args = parse_args()

    if args.export_coco:
        logger.info(f"Compiling COCO dataset from {args.masks_dir} -> {args.output_coco}...")
        coco_doc = export_coco_annotations(
            masks_dir=args.masks_dir,
            output_coco_path=args.output_coco,
            images_dir=args.images_dir,
            val_split=args.val_split,
        )
        logger.info(
            f"Successfully compiled COCO annotations: {len(coco_doc['images'])} images, "
            f"{len(coco_doc['annotations'])} annotations."
        )
        return

    annotator = PrecisionSAM2Annotator(
        images_dir=args.images_dir,
        output_dir=args.output_dir,
        coco_output=args.output_coco,
        single_image=args.single_image,
        checkpoint_path=args.checkpoint,
        config_path=args.config,
        resume_last=args.resume_last,
        resume_unannotated=args.resume_unannotated,
        target_voucher=args.voucher,
        target_index=args.index,
        tier=args.tier,
        vouchers_csv=args.vouchers_csv,
    )
    logger.info("Launching Precision SAM 2 Annotator GUI...")
    annotator.run()


if __name__ == "__main__":
    main()
