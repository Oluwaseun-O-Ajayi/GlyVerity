"""Extracted-ion chromatograms (XICs) from centroided MS1 scans."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np


def xic(centroid_file, keep: np.ndarray | None, targets: Sequence[float], ppm: float = 10.0,
        t_lo: float | None = None, t_hi: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Summed intensity within ``ppm`` of each target m/z, per kept scan.

    Returns
    -------
    times, intensities
        ``times`` has one entry per kept scan; ``intensities`` has shape
        (n_scans, n_targets).
    """
    times = np.asarray(centroid_file.xlabels, float)
    sel = np.ones(times.size, bool) if keep is None else np.asarray(keep, bool).copy()
    if t_lo is not None:
        sel &= times >= t_lo
    if t_hi is not None:
        sel &= times <= t_hi
    idx = np.flatnonzero(sel)
    tg = np.asarray(targets, float)
    lo, hi = tg * (1 - ppm * 1e-6), tg * (1 + ppm * 1e-6)
    out = np.zeros((idx.size, tg.size))
    for row, i in enumerate(idx):
        mz, it = centroid_file.scan(int(i))
        mz = np.asarray(mz, float)
        it = np.asarray(it, float)
        if mz.size == 0:
            continue
        order = np.argsort(mz)
        mz, it = mz[order], it[order]
        cum = np.concatenate([[0.0], np.cumsum(it)])
        a = np.searchsorted(mz, lo, side="left")
        b = np.searchsorted(mz, hi, side="right")
        out[row] = cum[b] - cum[a]
    return times[idx], out


@dataclass(frozen=True)
class PeakSummary:
    """Apex and area of one XIC."""

    apex_time: float
    apex_intensity: float
    area: float


def summarise(times: np.ndarray, trace: np.ndarray) -> PeakSummary:
    """Apex time, apex intensity and trapezoidal area (intensity x min)."""
    if trace.size == 0 or not np.any(trace > 0):
        return PeakSummary(float("nan"), 0.0, 0.0)
    k = int(np.argmax(trace))
    return PeakSummary(float(times[k]), float(trace[k]), float(np.trapezoid(trace, times)))
