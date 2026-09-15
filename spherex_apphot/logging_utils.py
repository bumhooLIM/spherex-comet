"""
Run logging.

The primitive pipeline printed progress to stdout and, on failure, appended a
bare target name to ``failed_targets.log``; the exception itself was swallowed
by a bare ``except`` (review item R3).  This module gives every run a real log
file that records the configuration, per-exposure warnings *with the filename
that caused them*, and a traceback for anything that fails.
"""

from __future__ import annotations

import logging
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

__all__ = ["setup_logging", "get_logger", "log_config", "utcnow_iso"]

_LOGGER_NAME = "spherex_apphot"
_FMT = "%(asctime)s %(levelname)-7s %(name)s: %(message)s"
_DATEFMT = "%Y-%m-%d %H:%M:%S"


def utcnow_iso() -> str:
    """Current UTC time as a second-resolution ISO-8601 string."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Return the package logger, or a child of it."""
    base = logging.getLogger(_LOGGER_NAME)
    return base if name is None else base.getChild(name)


def setup_logging(
    log_dir: Optional[Path] = None,
    *,
    run_name: str = "run",
    level: int = logging.INFO,
    console: bool = True,
    console_level: Optional[int] = None,
) -> Path | None:
    """
    Configure the package logger with a console handler and a file handler.

    Parameters
    ----------
    log_dir : Path or None
        Directory for the log file.  ``None`` disables the file handler, which
        is what the notebook wants when it only needs console output.
    run_name : str
        Stem of the log file; a UTC timestamp is appended.
    level : int
        Threshold for the file handler.
    console : bool
        Whether to attach a stdout handler.
    console_level : int, optional
        Threshold for the console handler (defaults to ``level``).

    Returns
    -------
    Path or None
        The log file that was opened, or ``None``.

    Notes
    -----
    Existing handlers are removed first, so calling this twice from a notebook
    cell does not duplicate every message.
    """
    logger = get_logger()
    logger.setLevel(min(level, console_level or level))
    for h in list(logger.handlers):
        logger.removeHandler(h)
        h.close()
    logger.propagate = False

    fmt = logging.Formatter(_FMT, datefmt=_DATEFMT)

    if console:
        ch = logging.StreamHandler(stream=sys.stdout)
        ch.setLevel(console_level or level)
        ch.setFormatter(fmt)
        logger.addHandler(ch)

    if log_dir is None:
        return None

    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    path = log_dir / f"{run_name}_{stamp}.log"
    fh = logging.FileHandler(path, mode="w", encoding="utf-8")
    fh.setLevel(level)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    logger.info("log file: %s", path)
    logger.info("python %s | pid %d | cwd %s", sys.version.split()[0], os.getpid(), Path.cwd())
    return path


def log_config(cfg, logger: Optional[logging.Logger] = None) -> None:
    """Write the full configuration to the log, one line per field."""
    log = logger or get_logger()
    log.info("config hash %s", cfg.hash)
    for line in cfg.to_json().splitlines():
        log.debug("  cfg %s", line.rstrip())
