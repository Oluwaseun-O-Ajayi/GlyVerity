"""Estimate the retention-time offset between two detectors (for example FLD and MS)."""

from __future__ import annotations

import numpy as np


def _resample(times: np.ndarray, values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    y = np.interp(grid, times, values)
    y = y - np.mean(y)
    s = np.std(y)
    return y / s if s > 0 else y


def estimate_lag(
    t_ref: np.ndarray,
    y_ref: np.ndarray,
    t_other: np.ndarray,
    y_other: np.ndarray,
    window: tuple[float, float],
    max_lag: float = 1.0,
    step: float = 0.001,
) -> tuple[float, float]:
    """Return (lag, correlation) that best aligns ``other`` to ``ref``.

    A positive lag means features appear later in ``other`` than in ``ref``;
    subtract it from ``other`` times to align. Both traces are resampled on a
    common grid of spacing ``step`` (min) inside ``window``, centred on their
    mean and scaled to unit standard deviation, so the returned value is the
    Pearson correlation (between -1 and 1) at the best lag.
    """
    lo, hi = window
    if hi <= lo:
        raise ValueError("window end must be greater than start")
    grid = np.arange(lo, hi, step)
    a = _resample(np.asarray(t_ref, float), np.asarray(y_ref, float), grid)
    n_lag = round(max_lag / step)
    best_lag, best_r = 0.0, -np.inf
    for k in range(-n_lag, n_lag + 1):
        shifted = np.asarray(t_other, float) - k * step
        b = _resample(shifted, np.asarray(y_other, float), grid)
        r = float(np.dot(a, b) / grid.size)
        if r > best_r:
            best_lag, best_r = k * step, r
    return best_lag, best_r
