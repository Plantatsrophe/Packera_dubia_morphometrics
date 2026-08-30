#!/usr/bin/env python3
"""
Convenience wrapper redirecting to scripts/vision/configure_leafmachine2.py
"""
import sys
from pathlib import Path

VISION_SCRIPT = Path(__file__).resolve().parents[1] / "vision" / "configure_leafmachine2.py"
if str(VISION_SCRIPT.parent.parent) not in sys.path:
    sys.path.insert(0, str(VISION_SCRIPT.parent.parent))

from scripts.vision.configure_leafmachine2 import (
    main,
    generate_high_performance_config,
    load_config_yaml,
    save_config_yaml,
    run_leafmachine2,
)

if __name__ == "__main__":
    main()
