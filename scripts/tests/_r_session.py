#!/usr/bin/env python3
"""
scripts/tests/_r_session.py
===========================
Persistent per-stage warm Rscript worker session for the Packera morphometrics
test suite. Boots at most one R interpreter per pipeline stage (Steps 03 to 06),
preloads heavy CRAN/GitHub namespaces once in RAM, snapshots and restores
.GlobalEnv and environment variables between test invocations, and returns
standard subprocess.CompletedProcess objects.
"""

from __future__ import annotations

import atexit
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Union

PROJECT_ROOT = Path(__file__).resolve().parents[2]

_MANAGED_ENV_KEYS = (
    "PACKERA_DB_PATH",
    "PACKERA_RUN_DIR",
    "PACKERA_SKIP_FIGURES",
)

_WORKER_R_CODE = r"""
local({
  options(rgl.useNULL = TRUE)
  boot_wd <- getwd()
  preload_script <- Sys.getenv("PACKERA_PRELOAD_R_SCRIPT", unset = "")
  Sys.unsetenv("PACKERA_PRELOAD_R_SCRIPT")

  if (nzchar(preload_script)) {
    boot_err <- NULL
    tryCatch(
      suppressMessages(source(preload_script, local = .GlobalEnv)),
      error = function(e) { boot_err <<- conditionMessage(e) }
    )
    if (!is.null(boot_err)) {
      cat("<<<PACKERA_BOOT_ERROR>>>:", boot_err, "\n", sep = "")
      flush(stdout())
      quit(save = "no", status = 2L)
    }
  }

  while (grDevices::dev.cur() > 1L) {
    tryCatch(grDevices::dev.off(), error = function(e) break)
  }

  boot_snapshot <- as.list.environment(.GlobalEnv, all.names = TRUE)
  boot_names <- names(boot_snapshot)

  restore_global_state <- function() {
    while (grDevices::dev.cur() > 1L) {
      tryCatch(grDevices::dev.off(), error = function(e) break)
    }
    tryCatch(setwd(boot_wd), error = function(e) NULL)
    cur_names <- ls(envir = .GlobalEnv, all.names = TRUE)
    extra_names <- setdiff(cur_names, boot_names)
    if (length(extra_names) > 0L) {
      rm(list = extra_names, envir = .GlobalEnv)
    }
    for (nm in boot_names) {
      if (!exists(nm, envir = .GlobalEnv, inherits = FALSE) ||
          !identical(get(nm, envir = .GlobalEnv, inherits = FALSE), boot_snapshot[[nm]])) {
        assign(nm, boot_snapshot[[nm]], envir = .GlobalEnv)
      }
    }
    if (exists(".packera_taxon_syn_env", envir = .GlobalEnv, inherits = FALSE)) {
      syn_e <- get(".packera_taxon_syn_env", envir = .GlobalEnv, inherits = FALSE)
      if (is.environment(syn_e)) {
        syn_e$db_sig <- ""
        syn_e$syn_map <- character(0)
        syn_e$py_memo <- character(0)
      }
    }
    if (exists(".AUTH_WEIGHT_CFG_ENV", envir = .GlobalEnv, inherits = FALSE)) {
      aw_e <- get(".AUTH_WEIGHT_CFG_ENV", envir = .GlobalEnv, inherits = FALSE)
      if (is.environment(aw_e)) {
        rm(list = ls(envir = aw_e, all.names = TRUE), envir = aw_e)
      }
    }
  }

  cat("<<<PACKERA_BOOT_READY>>>\n")
  flush(stdout())

  con_in <- file("stdin", open = "r")
  on.exit(tryCatch(close(con_in), error = function(e) NULL), add = TRUE)

  repeat {
    spec_line <- readLines(con_in, n = 1L, warn = FALSE)
    if (length(spec_line) == 0L) break
    spec_line <- trimws(spec_line)
    if (!nzchar(spec_line)) next
    if (identical(spec_line, "__QUIT__")) break

    hdr <- readLines(spec_line, warn = FALSE, encoding = "UTF-8")
    mode <- hdr[1L]
    entry_fn <- hdr[2L]
    req_cwd <- hdr[3L]
    script_path <- hdr[4L]
    args_file <- hdr[5L]
    env_keys_file <- hdr[6L]
    env_vals_file <- hdr[7L]
    out_file <- hdr[8L]
    err_file <- hdr[9L]
    status_file <- hdr[10L]
    done_token <- hdr[11L]

    cli_args <- if (nzchar(args_file) && file.exists(args_file)) {
      readLines(args_file, warn = FALSE, encoding = "UTF-8")
    } else {
      character(0)
    }
    env_keys <- if (nzchar(env_keys_file) && file.exists(env_keys_file)) {
      readLines(env_keys_file, warn = FALSE, encoding = "UTF-8")
    } else {
      character(0)
    }
    env_vals <- if (nzchar(env_vals_file) && file.exists(env_vals_file)) {
      readLines(env_vals_file, warn = FALSE, encoding = "UTF-8")
    } else {
      character(0)
    }

    prev_env <- stats::setNames(vector("list", length(env_keys)), env_keys)
    for (idx in seq_along(env_keys)) {
      k <- env_keys[idx]
      prev_env[[k]] <- Sys.getenv(k, unset = NA_character_)
      v <- if (idx <= length(env_vals)) env_vals[idx] else "__UNSET__"
      if (identical(v, "__UNSET__")) {
        Sys.unsetenv(k)
      } else {
        do.call(Sys.setenv, stats::setNames(list(v), k))
      }
    }

    err_chunks <- character(0)
    job_ok <- TRUE

    tryCatch(
      withCallingHandlers(
        utils::capture.output({
          if (nzchar(req_cwd)) setwd(req_cwd)
          if (identical(mode, "cli")) {
            if (identical(entry_fn, "main")) {
              main(cli_args)
            } else {
              fn_obj <- get(entry_fn, envir = .GlobalEnv, mode = "function")
              opts_obj <- parse_args_robust(cli_args)
              fn_obj(opts_obj)
            }
          } else {
            source(script_path, local = .GlobalEnv)
          }
          invisible(NULL)
        }, file = out_file),
        message = function(m) {
          err_chunks <<- c(err_chunks, conditionMessage(m))
          invokeRestart("muffleMessage")
        },
        warning = function(w) {
          err_chunks <<- c(err_chunks, paste0("Warning message: ", conditionMessage(w), "\n"))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        job_ok <<- FALSE
        err_chunks <<- c(err_chunks, paste0("Error: ", conditionMessage(e), "\n"))
      }
    )

    for (k in names(prev_env)) {
      pv <- prev_env[[k]]
      if (is.na(pv)) {
        Sys.unsetenv(k)
      } else {
        do.call(Sys.setenv, stats::setNames(list(pv), k))
      }
    }
    restore_global_state()

    writeLines(paste(err_chunks, collapse = ""), con = err_file, sep = "")
    writeLines(if (job_ok) "0" else "1", con = status_file)
    cat(done_token, "\n", sep = "")
    flush(stdout())
  }
})
"""

_STAGE_CONFIG: Dict[str, tuple[str, str]] = {
    "step03": ("scripts/morphometrics/03_fourier_extractor.R", "run_fourier_extraction"),
    "step04": ("scripts/morphometrics/04_gmm_morphotools.R", "run_gmm_morphotools_pipeline"),
    "step05": ("scripts/analysis/05_multimodal_spatial_rf.R", "main"),
    "step06": ("scripts/analysis/06_triage_dashboard_synthesis.R", "main"),
}

_SCRIPT_TO_STAGE: Dict[str, str] = {
    "03_fourier_extractor.R": "step03",
    "04_gmm_morphotools.R": "step04",
    "05_multimodal_spatial_rf.R": "step05",
    "06_triage_dashboard_synthesis.R": "step06",
}


def rscript_available() -> bool:
    """Return True when Rscript is installed and fast Python-only audit mode is off."""
    if os.environ.get("PACKERA_FAST_AUDIT", "").strip() == "1":
        return False
    return shutil.which("Rscript") is not None


class WarmRSession:
    """Persistent warm Rscript process preloaded with a stage's R environment."""

    def __init__(self, stage_key: str) -> None:
        if stage_key not in _STAGE_CONFIG:
            raise ValueError(f"Unknown R stage key: {stage_key}")
        self.stage_key = stage_key
        self.preload_rel, self.entry_fn = _STAGE_CONFIG[stage_key]
        self._tmp_dir: Optional[tempfile.TemporaryDirectory[str]] = None
        self._proc: Optional[subprocess.Popen[str]] = None
        self._job_seq = 0

    def _ensure_started(self) -> None:
        if self._proc is not None and self._proc.poll() is None:
            return
        self.close()
        rscript_bin = shutil.which("Rscript")
        if not rscript_bin:
            raise RuntimeError("Rscript binary not found on PATH.")
        self._tmp_dir = tempfile.TemporaryDirectory(prefix=f"packera_r_{self.stage_key}_")
        tmp_path = Path(self._tmp_dir.name)
        worker_file = tmp_path / "warm_worker.R"
        worker_file.write_text(_WORKER_R_CODE, encoding="utf-8")

        boot_env = os.environ.copy()
        boot_env["PACKERA_PRELOAD_R_SCRIPT"] = str(PROJECT_ROOT / self.preload_rel)
        for k in _MANAGED_ENV_KEYS:
            boot_env.pop(k, None)

        self._proc = subprocess.Popen(
            [rscript_bin, "--vanilla", str(worker_file)],
            cwd=str(PROJECT_ROOT),
            env=boot_env,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert self._proc.stdout is not None
        boot_lines: List[str] = []
        while True:
            line = self._proc.stdout.readline()
            if not line:
                rc = self._proc.wait()
                self.close()
                raise RuntimeError(
                    f"Warm R worker ({self.stage_key}) exited during boot (code {rc}):\n"
                    + "".join(boot_lines)
                )
            if "<<<PACKERA_BOOT_READY>>>" in line:
                break
            if line.startswith("<<<PACKERA_BOOT_ERROR>>>:"):
                rc = self._proc.wait()
                self.close()
                raise RuntimeError(
                    f"Warm R worker ({self.stage_key}) failed to preload {self.preload_rel} (code {rc}):\n"
                    + "".join(boot_lines)
                    + line
                )
            boot_lines.append(line)

    def _dispatch(
        self,
        *,
        mode: str,
        cli_args: Sequence[str] = (),
        script_path: Optional[Path] = None,
        cwd: Optional[Union[str, Path]] = None,
        env: Optional[Mapping[str, str]] = None,
        skip_figures: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        self._ensure_started()
        assert self._proc is not None and self._proc.stdin is not None and self._proc.stdout is not None
        assert self._tmp_dir is not None

        self._job_seq += 1
        jid = f"{self._job_seq:04d}"
        jdir = Path(self._tmp_dir.name)

        args_file = jdir / f"job_{jid}.args"
        env_keys_file = jdir / f"job_{jid}.envk"
        env_vals_file = jdir / f"job_{jid}.envv"
        out_file = jdir / f"job_{jid}.out"
        err_file = jdir / f"job_{jid}.err"
        status_file = jdir / f"job_{jid}.status"
        spec_file = jdir / f"job_{jid}.spec"
        done_token = f"<<<PACKERA_RPC_DONE_{jid}>>>"

        args_file.write_text(
            "\n".join(str(a) for a in cli_args) + ("\n" if cli_args else ""),
            encoding="utf-8",
        )

        env_map: Dict[str, str] = {}
        for k in _MANAGED_ENV_KEYS:
            if env is not None and k in env:
                env_map[k] = str(env[k])
            elif env is None and k in os.environ:
                env_map[k] = os.environ[k]
            else:
                env_map[k] = "__UNSET__"
        if skip_figures:
            env_map["PACKERA_SKIP_FIGURES"] = "1"
        if env is not None:
            for k, v in env.items():
                if k not in os.environ or os.environ.get(k) != v:
                    if "\n" not in k and "\n" not in str(v):
                        env_map[k] = str(v)

        keys_list = list(env_map.keys())
        vals_list = [env_map[k] for k in keys_list]
        env_keys_file.write_text("\n".join(keys_list) + "\n", encoding="utf-8")
        env_vals_file.write_text("\n".join(vals_list) + "\n", encoding="utf-8")
        out_file.write_text("", encoding="utf-8")
        err_file.write_text("", encoding="utf-8")

        req_cwd = str(Path(cwd).resolve()) if cwd is not None else str(PROJECT_ROOT)
        s_path = str(script_path.resolve()) if script_path is not None else ""

        spec_lines = [
            mode,
            self.entry_fn if mode == "cli" else "",
            req_cwd,
            s_path,
            str(args_file),
            str(env_keys_file),
            str(env_vals_file),
            str(out_file),
            str(err_file),
            str(status_file),
            done_token,
        ]
        spec_file.write_text("\n".join(spec_lines) + "\n", encoding="utf-8")

        self._proc.stdin.write(str(spec_file) + "\n")
        self._proc.stdin.flush()

        passthrough_lines: List[str] = []
        while True:
            line = self._proc.stdout.readline()
            if not line:
                rc = self._proc.wait() or 1
                out_txt = out_file.read_text(encoding="utf-8", errors="replace") if out_file.exists() else ""
                err_txt = err_file.read_text(encoding="utf-8", errors="replace") if err_file.exists() else ""
                err_txt = err_txt + "".join(passthrough_lines)
                self.close()
                return subprocess.CompletedProcess(
                    args=["Rscript", self.preload_rel, *cli_args],
                    returncode=rc,
                    stdout=out_txt,
                    stderr=err_txt,
                )
            if done_token in line:
                break
            passthrough_lines.append(line)

        rc_txt = status_file.read_text(encoding="utf-8").strip() if status_file.exists() else "1"
        rc = 0 if rc_txt == "0" else 1
        out_txt = out_file.read_text(encoding="utf-8", errors="replace") if out_file.exists() else ""
        err_txt = err_file.read_text(encoding="utf-8", errors="replace") if err_file.exists() else ""
        if passthrough_lines:
            err_txt = "".join(passthrough_lines) + err_txt

        return subprocess.CompletedProcess(
            args=["Rscript", self.preload_rel if mode == "cli" else s_path, *cli_args],
            returncode=rc,
            stdout=out_txt,
            stderr=err_txt,
        )

    def run_cli(
        self,
        cli_args: Sequence[Union[str, Path]],
        *,
        cwd: Optional[Union[str, Path]] = None,
        env: Optional[Mapping[str, str]] = None,
        skip_figures: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        return self._dispatch(
            mode="cli",
            cli_args=[str(a) for a in cli_args],
            cwd=cwd,
            env=env,
            skip_figures=skip_figures,
        )

    def run_script(
        self,
        script_or_code: Union[Path, str],
        *,
        cwd: Optional[Union[str, Path]] = None,
        env: Optional[Mapping[str, str]] = None,
        skip_figures: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        self._ensure_started()
        assert self._tmp_dir is not None
        if isinstance(script_or_code, Path):
            s_path = script_or_code
        else:
            s_path = Path(self._tmp_dir.name) / f"inline_{self._job_seq + 1:04d}.R"
            s_path.write_text(script_or_code, encoding="utf-8")
        return self._dispatch(
            mode="script",
            script_path=s_path,
            cwd=cwd,
            env=env,
            skip_figures=skip_figures,
        )

    def close(self) -> None:
        if self._proc is not None:
            try:
                if self._proc.poll() is None and self._proc.stdin is not None:
                    try:
                        self._proc.stdin.write("__QUIT__\n")
                        self._proc.stdin.flush()
                        self._proc.stdin.close()
                    except Exception:
                        pass
                    try:
                        self._proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        self._proc.kill()
                        self._proc.wait(timeout=2)
            finally:
                if self._proc.stdout is not None:
                    try:
                        self._proc.stdout.close()
                    except Exception:
                        pass
                self._proc = None
        if self._tmp_dir is not None:
            try:
                self._tmp_dir.cleanup()
            except Exception:
                pass
            self._tmp_dir = None


_SESSIONS: Dict[str, WarmRSession] = {}


def get_warm_r_session(stage_key: str) -> WarmRSession:
    session = _SESSIONS.get(stage_key)
    if session is None:
        session = WarmRSession(stage_key)
        _SESSIONS[stage_key] = session
    return session


def run_warm_r_cmd(
    cmd: Sequence[Union[str, Path]],
    *,
    stage_key: Optional[str] = None,
    cwd: Optional[Union[str, Path]] = None,
    env: Optional[Mapping[str, str]] = None,
    skip_figures: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Drop-in warm replacement for subprocess.run(['Rscript', ...])."""
    parts = [str(c) for c in cmd]
    if len(parts) >= 3 and parts[1] == "-e":
        sk = stage_key or "step03"
        return get_warm_r_session(sk).run_script(parts[2], cwd=cwd, env=env, skip_figures=skip_figures)
    if len(parts) >= 2:
        script_name = Path(parts[1]).name
        resolved_stage = _SCRIPT_TO_STAGE.get(script_name)
        if resolved_stage is not None:
            return get_warm_r_session(resolved_stage).run_cli(
                parts[2:], cwd=cwd, env=env, skip_figures=skip_figures
            )
        sk = stage_key or "step04"
        return get_warm_r_session(sk).run_script(Path(parts[1]), cwd=cwd, env=env, skip_figures=skip_figures)
    raise ValueError(f"Unsupported Rscript command sequence: {parts}")


def close_all_warm_r_sessions() -> None:
    for sess in list(_SESSIONS.values()):
        sess.close()
    _SESSIONS.clear()


atexit.register(close_all_warm_r_sessions)
