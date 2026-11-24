"""Read detector traces exported as CSV (e.g. by tools/thermo_channels)."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from glyverity.io_agilent import Trace


def read_trace_csv(path: str | Path, column: str | None = None,
                   detector: str = "FLD") -> Trace:
    """Trace from a CSV whose first column is time in minutes.

    ``column`` names the signal column; by default the first column after time is used.
    Rows with a missing or non-numeric value are skipped.
    """
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    head, body = rows[0], rows[1:]
    k = 1 if column is None else head.index(column)
    t, v = [], []
    for r in body:
        try:
            ti, vi = float(r[0]), float(r[k])
        except (IndexError, ValueError):
            continue
        t.append(ti)
        v.append(vi)
    order = np.argsort(t)
    return Trace(head[k], detector, np.asarray(t)[order], np.asarray(v)[order],
                 {"source": str(path)})
