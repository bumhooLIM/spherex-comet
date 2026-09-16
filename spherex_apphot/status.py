"""
Machine-readable run status table.

``results/apphot/status.csv`` holds one row per (target, run) with enough information
to answer "what has been processed, what failed, and why" without grepping
logs, and it is what the batch driver consults to decide whether a target can
be skipped.

That resume check is also the fix for review item R1: the primitive orchestrator
tested for ``<objdesig>.csv`` while the worker wrote ``<objdesig without
spaces>.csv``, so every target whose designation contains a space -- most comets
-- was reprocessed on every run.  Both sides now go through :func:`slugify`.
"""

from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

from .logging_utils import get_logger, utcnow_iso

__all__ = ["slugify", "StatusLog", "STATUS_COLUMNS"]

log = get_logger("status")

STATUS_COLUMNS = [
    "objdesig", "slug", "status", "timestamp", "elapsed_s",
    "n_exposures", "n_exposures_ok", "n_missing_fits", "n_read_errors",
    "n_epochs", "n_rows_out", "n_ap_skipped", "n_badphot", "n_gaia",
    "config_hash", "output_csv", "message",
]

_SLUG_RE = re.compile(r"[^A-Za-z0-9_.+-]")


def slugify(objdesig: str) -> str:
    """
    Canonical filesystem name for a target designation.

    ``"2021 G2"`` -> ``"2021G2"``.  This is the *only* place the mapping is
    defined; the writer and the resume check both call it.
    """
    return _SLUG_RE.sub("", str(objdesig).strip())


class StatusLog:
    """
    Append-or-update table of per-target pipeline outcomes.

    The table is keyed on ``slug``: re-running a target replaces its row rather
    than appending a duplicate.  Writes go to a temporary file in the same
    directory and are then ``os.replace``-d into place, so an interrupted run
    cannot leave a half-written status table behind.
    """

    def __init__(self, path: Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    # -- reading ---------------------------------------------------------
    def read(self) -> pd.DataFrame:
        """Return the current table (empty, correctly-typed, if absent)."""
        if not self.path.exists():
            return pd.DataFrame(columns=STATUS_COLUMNS)
        try:
            return pd.read_csv(self.path, dtype={"objdesig": str, "slug": str})
        except (pd.errors.EmptyDataError, pd.errors.ParserError) as exc:
            log.warning("status table %s unreadable (%s); starting a new one", self.path, exc)
            return pd.DataFrame(columns=STATUS_COLUMNS)

    def is_done(self, objdesig: str) -> bool:
        """True if the target has a recorded ``status == 'ok'``."""
        df = self.read()
        if df.empty:
            return False
        hit = df[df["slug"] == slugify(objdesig)]
        return bool(len(hit) and (hit["status"].iloc[-1] == "ok"))

    # -- writing ---------------------------------------------------------
    def update(self, objdesig: str, status: str, **fields: Any) -> Dict[str, Any]:
        """
        Insert or replace the row for ``objdesig`` and flush to disk.

        Parameters
        ----------
        objdesig : str
            Target designation as it appears in the database.
        status : {'ok', 'failed', 'skipped', 'running', 'empty'}
            Outcome.  ``'running'`` is written before the work starts so that an
            interrupted run is visible as such rather than simply absent.
        **fields
            Any of :data:`STATUS_COLUMNS`; unknown keys are kept as extra
            columns rather than dropped.

        Returns
        -------
        dict
            The row that was written.
        """
        row: Dict[str, Any] = {c: None for c in STATUS_COLUMNS}
        row.update(objdesig=str(objdesig), slug=slugify(objdesig),
                   status=status, timestamp=utcnow_iso())
        row.update(fields)

        # Rebuilt from records rather than via pd.concat: concatenating a
        # single row whose optional fields are all None triggers pandas'
        # empty/all-NA dtype deprecation, and the table is small enough that
        # this is both clearer and free.
        existing = self.read()
        records = ([r for r in existing.to_dict("records")
                    if str(r.get("slug")) != row["slug"]] if len(existing) else [])
        records.append(row)
        extra = [c for r in records for c in r if c not in STATUS_COLUMNS]
        cols = STATUS_COLUMNS + list(dict.fromkeys(extra))
        df = pd.DataFrame(records, columns=cols).sort_values("slug", kind="stable")
        self._atomic_write(df)
        return row

    def _atomic_write(self, df: pd.DataFrame) -> None:
        fd, tmp = tempfile.mkstemp(dir=str(self.path.parent), suffix=".csv.tmp")
        try:
            with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
                df.to_csv(fh, index=False)
            os.replace(tmp, self.path)
        except BaseException:
            Path(tmp).unlink(missing_ok=True)
            raise

    def summary(self) -> str:
        """One-line-per-status summary of the table."""
        df = self.read()
        if df.empty:
            return "status table is empty"
        counts = df["status"].value_counts().to_dict()
        parts = ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
        return f"{len(df)} targets: {parts}"
