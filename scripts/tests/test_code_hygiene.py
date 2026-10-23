#!/usr/bin/env python3
"""
===============================================================================
Script: test_code_hygiene.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Static AST and text/grep hygiene test suite verifying that all active .py
    and .R files in the repository (excluding archived script folders) are 100%
    isolated from legacy archived code and data directories.
===============================================================================
"""

from __future__ import annotations

import ast
import sys
import unittest
from pathlib import Path
from typing import List

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

ARCHIVE_MARKER = "_" + "archive"
SELF_PATH = Path(__file__).resolve()


def iter_active_source_files() -> List[Path]:
    """Return all active .py and .R files under scripts/ and root, excluding archive directories."""
    files: List[Path] = []
    for root_py in PROJECT_ROOT.glob("*.py"):
        files.append(root_py.resolve())

    scripts_dir = PROJECT_ROOT / "scripts"
    for ext in ("*.py", "*.R"):
        for path in sorted(scripts_dir.rglob(ext)):
            resolved = path.resolve()
            rel_parts = resolved.relative_to(scripts_dir.resolve()).parts
            if ARCHIVE_MARKER in rel_parts:
                continue
            files.append(resolved)
    return files


class TestCodeHygiene(unittest.TestCase):
    """Static AST and grep hygiene checks for archive isolation across active .py and .R files."""

    def test_no_active_files_reference_legacy_folders_via_grep(self):
        """Scan all active .py and .R files and fail if any references the archive folder."""
        violations: List[str] = []
        for src_path in iter_active_source_files():
            content = src_path.read_text(encoding="utf-8", errors="ignore")
            for lineno, line in enumerate(content.splitlines(), start=1):
                if ARCHIVE_MARKER in line:
                    rel = src_path.relative_to(PROJECT_ROOT)
                    violations.append(f"{rel}:{lineno}: {line.strip()}")

        self.assertEqual(
            violations,
            [],
            f"Active source files must not reference '{ARCHIVE_MARKER}':\n" + "\n".join(violations),
        )

    def test_no_active_python_ast_imports_or_constants_reference_legacy_folders(self):
        """Parse AST of all active .py files and verify zero imports or string literals reference archive."""
        ast_violations: List[str] = []
        for src_path in iter_active_source_files():
            if src_path.suffix != ".py":
                continue
            source = src_path.read_text(encoding="utf-8")
            tree = ast.parse(source, filename=str(src_path))
            rel = src_path.relative_to(PROJECT_ROOT)

            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        if ARCHIVE_MARKER in alias.name:
                            ast_violations.append(f"{rel}:{node.lineno} import {alias.name}")
                elif isinstance(node, ast.ImportFrom):
                    mod_name = node.module or ""
                    if ARCHIVE_MARKER in mod_name:
                        ast_violations.append(f"{rel}:{node.lineno} from {mod_name} import ...")
                elif isinstance(node, ast.Constant) and isinstance(node.value, str):
                    if ARCHIVE_MARKER in node.value:
                        ast_violations.append(f"{rel}:{node.lineno} string constant '{node.value}'")

        self.assertEqual(
            ast_violations,
            [],
            f"AST inspection detected '{ARCHIVE_MARKER}' references in active Python files:\n"
            + "\n".join(ast_violations),
        )


if __name__ == "__main__":
    unittest.main()
