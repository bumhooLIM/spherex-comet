"""
Reading a target list from a spreadsheet or text file.

The project's working list of comets lives in
``data/reference/sx_comet_list_ver2607.xlsx``, so the pipeline needs to be able to take
"the targets in that file" as an argument rather than 68 designations typed on a
command line.

The reader is deliberately forgiving about layout -- these lists are maintained
by hand -- but strict about one thing: it drops summary rows.  The v2607 sheet
ends with a ``Total`` row whose ``desig`` cell holds ``68``, the *count* of
targets, which would otherwise be processed as a comet called "68".
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence

from .logging_utils import get_logger

__all__ = ["read_target_list", "DESIG_COLUMNS"]

log = get_logger("targetlist")

#: Column names accepted as the designation column, in priority order.
DESIG_COLUMNS = ("desig", "objdesig", "designation", "target", "name", "comet")

#: Values in the first column that mark a summary row rather than a target.
_SUMMARY_MARKERS = {"total", "sum", "count", "n", "subtotal"}


def read_target_list(
    path: Path,
    column: Optional[str] = None,
    sheet: Optional[str] = None,
) -> List[str]:
    """
    Read target designations from ``.xlsx``, ``.csv`` or a plain text file.

    Parameters
    ----------
    path : Path
        Spreadsheet or text file.  A text file is read one designation per line,
        with ``#`` comments and blank lines ignored.
    column : str, optional
        Designation column.  Auto-detected from :data:`DESIG_COLUMNS` otherwise.
    sheet : str, optional
        Worksheet name; the first sheet is used otherwise.

    Returns
    -------
    list of str
        Designations in file order, whitespace-trimmed, duplicates removed.

    Raises
    ------
    FileNotFoundError, ValueError
        If the file is missing, or no designation column can be identified.

    Notes
    -----
    Rows whose first column reads ``Total`` (or similar) are dropped: a summary
    row's designation cell typically holds a *count*, which would otherwise be
    processed as a target.
    """
    import pandas as pd

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"target list not found: {path}")

    if path.suffix.lower() in (".txt", ".lis", ".list", ""):
        names = [ln.split("#")[0].strip() for ln in path.read_text().splitlines()]
        out = [n for n in names if n]
        log.info("read %d target(s) from %s", len(out), path.name)
        return list(dict.fromkeys(out))

    if path.suffix.lower() in (".xlsx", ".xls", ".xlsm"):
        df = pd.read_excel(path, sheet_name=sheet or 0)
    else:
        df = pd.read_csv(path)

    col = column
    if col is None:
        lower = {str(c).strip().lower(): c for c in df.columns}
        col = next((lower[c] for c in DESIG_COLUMNS if c in lower), None)
    if col is None:
        raise ValueError(
            f"{path.name}: no designation column found "
            f"(looked for {DESIG_COLUMNS}; columns are {list(df.columns)})")

    work = df
    if len(df.columns):
        first = df[df.columns[0]].astype(str).str.strip().str.lower()
        summary = first.isin(_SUMMARY_MARKERS)
        if summary.any():
            log.info("dropping %d summary row(s) from %s", int(summary.sum()), path.name)
            work = df[~summary]

    names = (work[col].dropna().astype(str).str.strip())
    out = [n for n in names if n and n.lower() != "nan"]
    out = list(dict.fromkeys(out))
    log.info("read %d target(s) from %s (column %r)", len(out), path.name, col)
    return out
