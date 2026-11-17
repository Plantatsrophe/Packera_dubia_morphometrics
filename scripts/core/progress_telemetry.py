"""Atomic single-line JSONL progress telemetry and hardware-aware EWMA/SQLite ETA tracking."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sqlite3
import threading
import time
from typing import Any, Callable, Dict, Mapping, Optional, Sequence, Tuple, Union

_JSONL_WRITE_LOCK = threading.Lock()
_ACCEL_CACHE_LOCK = threading.Lock()
_CACHED_ACCEL_PROVIDER: Optional[str] = None
_INLINE_TRACKERS: Dict[Tuple[int, str, str], "StepETATracker"] = {}
_BENCHMARK_DDL = (
    "CREATE TABLE IF NOT EXISTS step_telemetry_benchmarks ("
    "step_id INTEGER NOT NULL, substage TEXT NOT NULL DEFAULT 'all', "
    "accelerator_provider TEXT NOT NULL, cpu_threads INTEGER NOT NULL, "
    "mean_items_per_sec REAL NOT NULL, n_runs INTEGER NOT NULL DEFAULT 1, "
    "updated_at TEXT NOT NULL, PRIMARY KEY (step_id, substage, accelerator_provider, cpu_threads));"
)


def detect_accelerator_provider() -> str:
    """Detect active ONNX Runtime execution provider with safe CPU fallback."""
    env_override = os.environ.get("PLANTMATICS_ACCELERATOR_PROVIDER", "").strip()
    if env_override:
        return env_override
    global _CACHED_ACCEL_PROVIDER
    if _CACHED_ACCEL_PROVIDER is not None:
        return _CACHED_ACCEL_PROVIDER
    with _ACCEL_CACHE_LOCK:
        if _CACHED_ACCEL_PROVIDER is not None:
            return _CACHED_ACCEL_PROVIDER
        provider = "CPUExecutionProvider"
        try:
            import onnxruntime as ort  # type: ignore

            avail = set(ort.get_available_providers() or ())
            for cand in ("CUDAExecutionProvider", "DmlExecutionProvider", "CoreMLExecutionProvider", "CPUExecutionProvider"):
                if cand in avail:
                    provider = cand
                    break
        except Exception:
            provider = "CPUExecutionProvider"
        _CACHED_ACCEL_PROVIDER = provider
        return provider


def detect_cpu_threads() -> int:
    """Resolve effective CPU worker thread count from config/env or host cores."""
    try:
        from scripts.core.config import resolve_cpu_threads

        val = resolve_cpu_threads()
        if val is not None and int(val) >= 1:
            return int(val)
    except Exception:
        pass
    for key in ("PLANTMATICS_MAX_WORKERS", "PACKERA_MAX_WORKERS"):
        raw = os.environ.get(key, "").strip()
        if raw.isdigit() and int(raw) >= 1:
            return int(raw)
    return max(1, int(os.cpu_count() or 1))


def _resolve_run_dir(run_dir: Optional[Union[str, Path]] = None) -> Optional[Path]:
    """Resolve active run directory from explicit argument, environment, or config."""
    if run_dir is not None:
        return Path(run_dir).resolve()
    for key in ("PLANTMATICS_RUN_DIR", "PACKERA_RUN_DIR"):
        raw = os.environ.get(key, "").strip()
        if raw:
            return Path(raw).resolve()
    try:
        from scripts.core.config import get_active_run_dir

        active = get_active_run_dir()
        if active is not None:
            return Path(active).resolve()
    except Exception:
        pass
    return None


def _resolve_telemetry_db_path(explicit_db: Optional[Union[str, Path]] = None) -> Optional[Path]:
    """Resolve SQLite path for telemetry benchmarks without polluting canonical DBs in tests."""
    if explicit_db is not None:
        raw = str(explicit_db).strip()
        return None if (not raw or raw == ":memory:") else Path(raw).resolve()
    if os.environ.get("PLANTMATICS_DEMO_MODE") == "1" or os.environ.get("PACKERA_DEMO_MODE") == "1":
        return None
    env_db = os.environ.get("PLANTMATICS_DB_PATH", "").strip() or os.environ.get("PACKERA_DB_PATH", "").strip()
    if env_db:
        if env_db == ":memory:":
            return None
        p = Path(env_db).resolve()
        return p if p.is_file() else None
    if "PYTEST_CURRENT_TEST" in os.environ:
        return None
    try:
        from scripts.core.db import resolve_db_path

        cand = resolve_db_path()
        if cand and str(cand) != ":memory:" and Path(cand).is_file():
            return Path(cand).resolve()
    except Exception:
        pass
    return None


def _write_jsonl_record(run_dir: Path, payload: Mapping[str, Any]) -> None:
    """Append one JSON record atomically to `<run_dir>/logs/progress.jsonl`."""
    log_path = Path(run_dir) / "logs" / "progress.jsonl"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    line_bytes = (json.dumps(payload, separators=(",", ":"), ensure_ascii=False) + "\n").encode("utf-8")
    with _JSONL_WRITE_LOCK:
        fd = os.open(str(log_path), os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o666)
        try:
            os.write(fd, line_bytes)
        finally:
            os.close(fd)


class StepETATracker:
    """Track step throughput via EWMA (alpha=0.25) with SQLite benchmark fallback."""

    def __init__(
        self, step: int, substage: str = "all", total: int = 0, cached_skipped: int = 0,
        unit: str = "items", alpha: float = 0.25, substage_weight: float = 1.0,
        substage_offset: float = 0.0, run_dir: Optional[Union[str, Path]] = None,
        db_path: Optional[Union[str, Path]] = None, accelerator_provider: Optional[str] = None,
        cpu_threads: Optional[int] = None, clock: Optional[Callable[[], float]] = None,
    ) -> None:
        self.step, self.substage = int(step), str(substage or "all")
        self.total, self.cached_skipped = max(0, int(total)), max(0, int(cached_skipped))
        self.n_active = max(0, self.total - self.cached_skipped)
        self.done, self.unit, self.alpha = 0, str(unit or "items"), float(alpha)
        self.substage_weight, self.substage_offset = float(substage_weight), float(substage_offset)
        self.run_dir = Path(run_dir).resolve() if run_dir is not None else None
        self.db_path = db_path
        self.accelerator_provider = str(accelerator_provider or detect_accelerator_provider())
        self.cpu_threads = max(1, int(cpu_threads if cpu_threads is not None else detect_cpu_threads()))
        self._clock = clock or time.monotonic
        self._start_t = self._last_t = float(self._clock())
        self._total_active_dt, self._ewma_updates = 0.0, 0
        self.rate_ewma: Optional[float] = None
        self._benchmark_recorded = False

    def _open_db(self) -> Optional[sqlite3.Connection]:
        db_file = _resolve_telemetry_db_path(self.db_path)
        if db_file is None or (self.db_path is None and not db_file.is_file()):
            return None
        try:
            db_file.parent.mkdir(parents=True, exist_ok=True)
            conn = sqlite3.connect(str(db_file), timeout=5.0)
            conn.execute("PRAGMA busy_timeout = 5000;")
            conn.execute(_BENCHMARK_DDL)
            return conn
        except sqlite3.Error:
            return None

    def query_benchmark_rate(self, step: Optional[int] = None, substage: Optional[str] = None) -> Optional[float]:
        """Query mean historical items/sec from SQLite `step_telemetry_benchmarks`."""
        eff_step, eff_sub = int(self.step if step is None else step), str(self.substage if substage is None else substage)
        conn = self._open_db()
        if conn is None:
            return None
        try:
            queries = (
                ("SELECT mean_items_per_sec FROM step_telemetry_benchmarks WHERE step_id = ? AND substage = ? AND accelerator_provider = ? AND cpu_threads = ? LIMIT 1;", (eff_step, eff_sub, self.accelerator_provider, self.cpu_threads)),
                ("SELECT mean_items_per_sec FROM step_telemetry_benchmarks WHERE step_id = ? AND substage IN (?, 'all') AND accelerator_provider = ? ORDER BY CASE WHEN substage = ? THEN 0 ELSE 1 END, n_runs DESC, updated_at DESC LIMIT 1;", (eff_step, eff_sub, self.accelerator_provider, eff_sub)),
                ("SELECT mean_items_per_sec FROM step_telemetry_benchmarks WHERE step_id = ? ORDER BY CASE WHEN substage = ? THEN 0 ELSE 1 END, n_runs DESC, updated_at DESC LIMIT 1;", (eff_step, eff_sub)),
            )
            for sql, params in queries:
                row = conn.execute(sql, params).fetchone()
                if row and row[0] is not None and float(row[0]) > 0.0:
                    return float(row[0])
            return None
        except sqlite3.Error:
            return None
        finally:
            conn.close()

    def record_benchmark(self, mean_items_per_sec: float) -> bool:
        """Upsert running mean items/sec into SQLite `step_telemetry_benchmarks`."""
        rate_val = float(mean_items_per_sec)
        conn = self._open_db() if rate_val > 0.0 else None
        if conn is None:
            return False
        sql = (
            "INSERT INTO step_telemetry_benchmarks (step_id, substage, accelerator_provider, cpu_threads, mean_items_per_sec, n_runs, updated_at) "
            "VALUES (?, ?, ?, ?, ?, 1, ?) ON CONFLICT(step_id, substage, accelerator_provider, cpu_threads) DO UPDATE SET "
            "mean_items_per_sec = ((step_telemetry_benchmarks.mean_items_per_sec * step_telemetry_benchmarks.n_runs) + excluded.mean_items_per_sec) / (step_telemetry_benchmarks.n_runs + 1), "
            "n_runs = step_telemetry_benchmarks.n_runs + 1, updated_at = excluded.updated_at;"
        )
        try:
            with conn:
                conn.execute(sql, (self.step, self.substage, self.accelerator_provider, self.cpu_threads, rate_val, datetime.now(timezone.utc).isoformat()))
            self._benchmark_recorded = True
            return True
        except sqlite3.Error:
            return False
        finally:
            conn.close()

    def estimate_step_eta(self, step: int, n_items: int, substage: str = "all") -> Optional[float]:
        """Estimate ETA in seconds for an unstarted downstream step from SQLite benchmarks."""
        if int(n_items) <= 0:
            return 0.0
        rate = self.query_benchmark_rate(step=int(step), substage=str(substage))
        return None if (rate is None or rate <= 0.0) else round(float(n_items) / float(rate), 4)

    def estimate_downstream_etas(self, steps: Sequence[int] = (1, 2, 3, 4, 5, 6, 7), step_totals: Optional[Mapping[int, int]] = None) -> Dict[int, Optional[float]]:
        """Estimate ETAs across downstream pipeline steps 1..7 from SQLite benchmarks."""
        totals = step_totals or {}
        return {int(s): self.estimate_step_eta(int(s), int(totals.get(int(s), 1))) for s in steps}

    def update(self, done: int, total: Optional[int] = None, cached_skipped: Optional[int] = None, status: str = "running", emit: bool = True, dt: Optional[float] = None) -> Dict[str, Any]:
        """Update active progress, EWMA rate, and ETA, optionally emitting a JSONL event."""
        if total is not None:
            self.total = max(0, int(total))
        if cached_skipped is not None:
            self.cached_skipped = max(0, int(cached_skipped))
        self.n_active = max(0, self.total - self.cached_skipped)
        done_i, now = max(0, int(done)), float(self._clock())
        delta_t = float(dt) if dt is not None else max(0.0, now - self._last_t)
        delta_done = done_i - self.done
        if delta_done > 0 and delta_t > 0.0:
            inst_rate = float(delta_done) / float(delta_t)
            self.rate_ewma = inst_rate if self.rate_ewma is None else (self.alpha * inst_rate + (1.0 - self.alpha) * self.rate_ewma)
            self._ewma_updates += 1
            self._total_active_dt += delta_t
            self._last_t = now
        elif delta_done > 0 and dt is None:
            self._last_t = now
        self.done = done_i

        eff_status = str(status or "running")
        if eff_status == "skipped_cached" or self.n_active == 0 or done_i >= self.n_active or eff_status == "completed":
            eta_sec: Optional[float] = 0.0
            eff_rate = float(self.rate_ewma or 0.0)
        elif self._ewma_updates >= 3 and self.rate_ewma is not None and self.rate_ewma > 0.0:
            eff_rate = float(self.rate_ewma)
            eta_sec = max(0.0, float(self.n_active - done_i)) / eff_rate
        else:
            bench_rate = self.query_benchmark_rate()
            if bench_rate is not None and bench_rate > 0.0:
                eff_rate = float(self.rate_ewma) if (self.rate_ewma is not None and self.rate_ewma > 0.0) else bench_rate
                eta_sec = max(0.0, float(self.n_active - done_i)) / bench_rate
            else:
                eff_rate, eta_sec = float(self.rate_ewma or 0.0), None

        if eff_status == "completed" and self.n_active > 0 and done_i > 0 and not self._benchmark_recorded:
            eff_dt = self._total_active_dt if self._total_active_dt > 0.0 else max(0.0, now - self._start_t)
            mean_rate = (float(done_i) / eff_dt) if eff_dt > 0.0 else float(self.rate_ewma or 0.0)
            if mean_rate > 0.0:
                self.record_benchmark(mean_rate)

        if emit:
            return emit_progress(
                step=self.step, substage=self.substage, done=done_i, total=self.total, unit=self.unit,
                cached_skipped=self.cached_skipped, substage_weight=self.substage_weight, substage_offset=self.substage_offset,
                status=eff_status, run_dir=self.run_dir, rate_items_per_sec=round(eff_rate, 4),
                eta_sec=None if eta_sec is None else round(eta_sec, 4),
                accelerator_provider=self.accelerator_provider, cpu_threads=self.cpu_threads, _from_tracker=True,
            )
        return {
            "step": self.step, "substage": self.substage, "done": done_i, "total": self.total,
            "cached_skipped": self.cached_skipped, "n_active": self.n_active,
            "rate_items_per_sec": round(eff_rate, 4), "eta_sec": None if eta_sec is None else round(eta_sec, 4),
        }

    def finish(self, done: Optional[int] = None, status: str = "completed") -> Dict[str, Any]:
        """Finalize step tracking, persist SQLite benchmark, and emit terminal telemetry event."""
        return self.update(done=(self.n_active if done is None else max(0, int(done))), status=status, emit=True)


def emit_progress(
    step: int, substage: str, done: int, total: int, unit: str = "items",
    cached_skipped: int = 0, substage_weight: float = 1.0, substage_offset: float = 0.0,
    status: str = "running", run_dir: Optional[Path] = None,
    *, rate_items_per_sec: Optional[float] = None, eta_sec: Optional[float] = None,
    accelerator_provider: Optional[str] = None, cpu_threads: Optional[int] = None,
    db_path: Optional[Union[str, Path]] = None, _from_tracker: bool = False,
) -> Dict[str, Any]:
    """Build and optionally append a single-line JSON progress telemetry record."""
    resolved_dir = _resolve_run_dir(run_dir)
    total_i, cached_i, done_i = max(0, int(total)), max(0, int(cached_skipped)), max(0, int(done))
    n_active, eff_status = max(0, total_i - cached_i), str(status or "running")
    provider = str(accelerator_provider or detect_accelerator_provider())
    threads_i = max(1, int(cpu_threads if cpu_threads is not None else detect_cpu_threads()))

    if not _from_tracker and rate_items_per_sec is None and eta_sec is None:
        t_key = (int(step), str(substage or "all"), str(resolved_dir or ""))
        trk = _INLINE_TRACKERS.get(t_key)
        if trk is None or done_i == 0 or eff_status == "started" or trk.total != total_i or trk.cached_skipped != cached_i or done_i < trk.done:
            trk = StepETATracker(
                int(step), str(substage or "all"), total_i, cached_i, unit,
                substage_weight=substage_weight, substage_offset=substage_offset,
                run_dir=resolved_dir, db_path=db_path, accelerator_provider=provider, cpu_threads=threads_i,
            )
            _INLINE_TRACKERS[t_key] = trk
        snap = trk.update(done=done_i, status=eff_status, emit=False)
        eff_rate, eff_eta = float(snap["rate_items_per_sec"]), snap["eta_sec"]
    else:
        eff_rate = float(rate_items_per_sec or 0.0)
        eff_eta = 0.0 if (eff_status in {"skipped_cached", "completed"} or n_active == 0 or done_i >= n_active) else eta_sec

    if eff_status in {"completed", "skipped_cached"}:
        local_frac = 1.0
    elif n_active > 0:
        local_frac = min(1.0, max(0.0, float(done_i) / float(n_active)))
    else:
        local_frac = 1.0 if (total_i > 0 and cached_i >= total_i) else 0.0
    progress_frac = round(min(1.0, max(0.0, float(substage_offset) + float(substage_weight) * local_frac)), 6)

    record: Dict[str, Any] = {
        "timestamp": datetime.now(timezone.utc).isoformat(), "step": int(step), "substage": str(substage),
        "done": done_i, "total": total_i, "cached_skipped": cached_i, "n_active": n_active,
        "unit": str(unit or "items"), "substage_weight": float(substage_weight), "substage_offset": float(substage_offset),
        "progress_frac": progress_frac, "status": eff_status, "rate_items_per_sec": round(eff_rate, 4),
        "eta_sec": None if eff_eta is None else round(float(eff_eta), 4),
        "accelerator_provider": provider, "cpu_threads": threads_i,
    }
    if resolved_dir is not None:
        try:
            _write_jsonl_record(resolved_dir, record)
        except OSError:
            pass
    return record
