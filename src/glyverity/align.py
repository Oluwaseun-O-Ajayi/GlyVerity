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


def targeted_ion_current(spectra, target_mz: np.ndarray, ppm: float = 10.0) -> tuple[np.ndarray,
                                                                                    np.ndarray]:
    """Per-spectrum summed intensity of peaks within ``ppm`` of any target m/z.

    ``spectra`` yields (time in minutes, m/z array, intensity array). Background ions dominate a
    total ion current in many glycan runs; summing only candidate glycan m/z values gives a
    trace that follows the glycans and can be aligned with a fluorescence trace.
    """
    tg = np.sort(np.asarray(target_mz, float))
    times, vals = [], []
    for t, mz, it in spectra:
        mz = np.asarray(mz, float)
        it = np.asarray(it, float)
        if mz.size == 0 or tg.size == 0:
            times.append(t)
            vals.append(0.0)
            continue
        i = np.clip(np.searchsorted(tg, mz), 1, tg.size - 1) if tg.size > 1 else np.zeros(mz.size, int)
        near = np.minimum(np.abs(mz - tg[np.maximum(i - 1, 0)]), np.abs(mz - tg[i]))
        times.append(t)
        vals.append(float(it[near <= mz * ppm * 1e-6].sum()))
    return np.asarray(times), np.asarray(vals)
