"""Run logging: console plus a per-run file under ``results/comspec/logs``."""

from __future__ import annotations

import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

__all__ = ["get_logger", "setup_logging", "utcnow_iso"]

_NAME = "spherex_comspec"
_FMT = "%(asctime)s %(levelname)-7s %(name)s: %(message)s"


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def get_logger(name: Optional[str] = None) -> logging.Logger:
    base = logging.getLogger(_NAME)
    return base if name is None else base.getChild(name)


def setup_logging(log_dir: Optional[Path] = None, *, run_name: str = "comspec",
                  level: int = logging.INFO, console: bool = True) -> Optional[Path]:
    """Attach a console handler and, when ``log_dir`` is given, a file handler."""
    log = get_logger()
    log.setLevel(level)
    for h in list(log.handlers):
        log.removeHandler(h)
        h.close()
    log.propagate = False
    fmt = logging.Formatter(_FMT, datefmt="%Y-%m-%d %H:%M:%S")
    if console:
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(fmt)
        log.addHandler(ch)
    if log_dir is None:
        return None
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    path = log_dir / f"{run_name}_{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}.log"
    fh = logging.FileHandler(path, mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    log.addHandler(fh)
    log.info("log file: %s", path)
    return path
