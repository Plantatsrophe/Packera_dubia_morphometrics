"""Stdlib-only flat-path, atomic file I/O, text decoding, and HTTP GET helpers."""

from __future__ import annotations

import csv
import io
import os
import tempfile
import urllib.request
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Union


def ensure_flat_filename(path_or_name: Union[str, Path]) -> str:
    """Validate that `path_or_name` is a flat single filename with no directory traversal."""
    raw = str(path_or_name).strip()
    if not raw or "/" in raw or "\\" in raw or raw in {".", ".."} or Path(raw).name != raw:
        raise ValueError(f"Unsafe or non-flat filename: {path_or_name!r}")
    return raw


def to_text(payload: Union[bytes, bytearray, str]) -> str:
    """Decode raw bytes/bytearray payload as UTF-8 text with replacement, or return str as is."""
    if isinstance(payload, str):
        return payload
    if isinstance(payload, (bytes, bytearray)):
        return bytes(payload).decode("utf-8", errors="replace")
    return str(payload)


def atomic_write_bytes(path: Union[str, Path], data: bytes) -> Path:
    """Atomically write `data` to `path` via a sibling temporary file, fsync, and replace."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        dir=str(target.parent),
        prefix=f".{target.stem}_",
        suffix=f"{target.suffix}.tmp",
    )
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp_path, target)
    except Exception:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass
        raise
    return target


def atomic_write_text(path: Union[str, Path], text: str, encoding: str = "utf-8") -> Path:
    """Atomically write string `text` to `path` encoded as `encoding`."""
    return atomic_write_bytes(path, (text or "").encode(encoding))


def atomic_write_csv_rows(
    path: Union[str, Path],
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
) -> Path:
    """Atomically serialize row dictionaries with fixed `fieldnames` to UTF-8 CSV at `path`."""
    cols = list(fieldnames)
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    writer.writerows({k: ("" if row.get(k) is None else row.get(k)) for k in cols} for row in rows)
    return atomic_write_text(path, buf.getvalue())


def default_http_get(
    url: str,
    timeout: float = 30.0,
    user_agent: str = "Plantmatics/2.0 (UNC Herbarium NCU)",
    *,
    headers: Optional[Mapping[str, str]] = None,
    data: Optional[bytes] = None,
) -> bytes:
    """Perform a simple HTTP GET (or POST when `data` is given) returning raw response bytes."""
    req_headers: dict[str, str] = {"User-Agent": user_agent}
    if headers:
        req_headers.update({str(k): str(v) for k, v in headers.items()})
    req = urllib.request.Request(
        url,
        data=data,
        headers=req_headers,
        method="POST" if data is not None else "GET",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()

