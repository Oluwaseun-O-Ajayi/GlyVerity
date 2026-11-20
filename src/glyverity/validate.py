"""Validation of MS glycan abundances against fluorescence (FLD) peak areas.

A 2-AB label carries one fluorophore per glycan, so FLD peak area is proportional to the molar
amount of each labelled glycan, independent of ionization. MS abundance is not: sialylated
glycans spread their signal over Na/K salts, isotope-error picks and dimers that a protonated-
ion measurement misses. If the cleaning step is right, summing every ion assigned to a species
should bring MS shares closer to FLD shares than the protonated ions alone, and should flatten
any trend of MS/FLD with the number of sialic acids.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np


@dataclass
class FLDWindow:
    """One FLD peak between two valleys, with its baseline-corrected area."""

    start: float
    end: float
    apex: float
    area: float
    members: list[int] = field(default_factory=list)  # indices of the species inside


def _smooth(y: np.ndarray, n: int) -> np.ndarray:
    return np.convolve(y, np.ones(n) / n, mode="same") if n > 1 else y


def fld_window(times: np.ndarray, values: np.ndarray, t: float,
               smooth: int = 5, max_half: float = 1.0) -> tuple[float, float]:
    """Valley-to-valley bounds of the FLD peak containing time ``t``.

    From the trace point nearest ``t`` (moved uphill to the local apex), each bound is the first
    local minimum of the smoothed trace, or ``max_half`` min from the apex, whichever is closer.
    """
    x = np.asarray(times, float)
    y = _smooth(np.asarray(values, float), smooth)
    i = int(np.clip(np.searchsorted(x, t), 0, x.size - 1))
    while 0 < i < x.size - 1 and max(y[i - 1], y[i + 1]) > y[i]:  # climb to the apex
        i = i - 1 if y[i - 1] > y[i + 1] else i + 1
    lo = i
    while lo > 0 and y[lo - 1] <= y[lo] and x[i] - x[lo - 1] <= max_half:
        lo -= 1
    hi = i
    while hi < x.size - 1 and y[hi + 1] <= y[hi] and x[hi + 1] - x[i] <= max_half:
        hi += 1
    return float(x[lo]), float(x[hi])


def fld_area(times: np.ndarray, values: np.ndarray, start: float, end: float) -> float:
    """Area above the straight baseline joining the trace at ``start`` and ``end``."""
    x = np.asarray(times, float)
    y = np.asarray(values, float)
    m = (x >= start) & (x <= end)
    if m.sum() < 2:
        return 0.0
    xs, ys = x[m], y[m]
    base = np.interp(xs, [xs[0], xs[-1]], [ys[0], ys[-1]])
    return float(np.trapezoid(np.clip(ys - base, 0.0, None), xs))


def windows_for(times: np.ndarray, values: np.ndarray, apexes: Sequence[float],
                smooth: int = 5, max_half: float = 1.0) -> list[FLDWindow]:
    """FLD windows for species apexes (FLD clock); species in the same peak share a window."""
    out: list[FLDWindow] = []
    for k, t in enumerate(apexes):
        lo, hi = fld_window(times, values, t, smooth, max_half)
        for w in out:
            if abs(w.start - lo) < 1e-9 and abs(w.end - hi) < 1e-9:
                w.members.append(k)
                break
        else:
            x = np.asarray(times, float)
            m = (x >= lo) & (x <= hi)
            apex = float(x[m][np.argmax(np.asarray(values, float)[m])]) if m.any() else t
            out.append(FLDWindow(lo, hi, apex, fld_area(times, values, lo, hi), [k]))
    return out


def shares(v: Sequence[float]) -> np.ndarray:
    a = np.asarray(v, float)
    s = a.sum()
    return a / s if s > 0 else np.full(a.shape, np.nan)


def total_variation(a: Sequence[float], b: Sequence[float]) -> float:
    """Half the summed absolute difference of two share vectors (0 = identical, 1 = disjoint)."""
    return float(0.5 * np.sum(np.abs(shares(a) - shares(b))))


def log_ratio_slope(ms: Sequence[float], fld: Sequence[float], x: Sequence[float]) -> float:
    """Least-squares slope of log10(MS share / FLD share) against ``x`` (e.g. sialic acids)."""
    r = np.log10(shares(ms) / shares(fld))
    xx = np.asarray(x, float)
    ok = np.isfinite(r)
    if ok.sum() < 2 or np.ptp(xx[ok]) == 0:
        return float("nan")
    return float(np.polyfit(xx[ok], r[ok], 1)[0])
