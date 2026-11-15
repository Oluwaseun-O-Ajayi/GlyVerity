"""Summed MS1 spectra over a retention-time window and composition matching.

This is a first-pass annotator for inspecting features by hand. It is not the
FDR-controlled assignment step; it lists every candidate within tolerance so
that ambiguity stays visible.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from glyverity import chemistry as ch

ISOTOPE_SPACING = 1.0033548378  # 13C minus 12C (u)


def merge_peaks(mz: np.ndarray, intensity: np.ndarray, ppm: float = 10.0) -> tuple[np.ndarray, np.ndarray]:
    """Merge centroids from many scans into one peak list.

    Peaks are sorted by m/z and grouped while each step to the next peak is
    within ``ppm``. Each group gives its intensity-weighted mean m/z and its
    summed intensity.
    """
    mz = np.asarray(mz, float)
    intensity = np.asarray(intensity, float)
    if mz.size == 0:
        return mz, intensity
    order = np.argsort(mz)
    mz, intensity = mz[order], intensity[order]
    breaks = np.flatnonzero(np.diff(mz) > mz[:-1] * ppm * 1e-6) + 1
    out_mz, out_int = [], []
    for seg in np.split(np.arange(mz.size), breaks):
        w = intensity[seg]
        total = float(w.sum())
        out_mz.append(float(np.dot(mz[seg], w) / total) if total > 0 else float(mz[seg].mean()))
        out_int.append(total)
    return np.array(out_mz), np.array(out_int)


def window_spectrum(centroid_file, keep: np.ndarray | None, t_lo: float, t_hi: float,
                    ppm: float = 10.0) -> tuple[np.ndarray, np.ndarray, int]:
    """Sum the centroided scans with ``t_lo <= time <= t_hi`` (and ``keep`` True).

    Returns merged m/z, merged intensity and the number of scans used.
    """
    times = np.asarray(centroid_file.xlabels, float)
    sel = (times >= t_lo) & (times <= t_hi)
    if keep is not None:
        sel &= keep
    idx = np.flatnonzero(sel)
    mzs, ints = [], []
    for i in idx:
        m, y = centroid_file.scan(int(i))
        mzs.append(np.asarray(m, float))
        ints.append(np.asarray(y, float))
    if not mzs:
        return np.array([]), np.array([]), 0
    m, y = merge_peaks(np.concatenate(mzs), np.concatenate(ints), ppm)
    return m, y, int(idx.size)


def enumerate_compositions(ranges: Mapping[str, tuple[int, int]]) -> list[dict[str, int]]:
    """All compositions within inclusive count ranges, excluding the empty one."""
    names = list(ranges)
    grids = [range(lo, hi + 1) for lo, hi in ranges.values()]
    out = []
    for counts in itertools.product(*grids):
        comp = {n: c for n, c in zip(names, counts) if c}
        if comp:
            out.append(comp)
    return out


def composition_string(comp: Mapping[str, int]) -> str:
    """Compact text form, for example ``Hex5HexNAc4NeuAc2``."""
    order = ["Hex", "HexNAc", "dHex", "NeuAc", "NeuGc"]
    return "".join(f"{k}{comp[k]}" for k in order if comp.get(k))


DEFAULT_CARRIERS: Sequence[Mapping[str, int]] = (
    {"H+": 1}, {"Na+": 1}, {"K+": 1}, {"NH4+": 1},
)


def carrier_sets(charge: int) -> list[dict[str, int]]:
    """Protons only, or protons plus one Na+, K+ or NH4+, for a given charge."""
    sets = [{"H+": charge}]
    for other in ("Na+", "K+", "NH4+"):
        s = {other: 1}
        if charge > 1:
            s["H+"] = charge - 1
        sets.append(s)
    return sets


@dataclass(frozen=True)
class Match:
    """One candidate explanation of an observed m/z."""

    composition: str
    label: str
    charge: int
    carriers: str
    theoretical_mz: float
    ppm_error: float


def match_mz(observed: float, compositions: Iterable[Mapping[str, int]], labels: Sequence[str],
             charges: Sequence[int] = (1, 2, 3), ppm: float = 20.0) -> list[Match]:
    """Every composition, label, charge and carrier combination within ``ppm`` of ``observed``."""
    out = []
    for comp in compositions:
        for label in labels:
            neutral = ch.glycan_neutral_mass(comp, label)
            for z in charges:
                for car in carrier_sets(z):
                    theo = ch.mz(neutral, z, car)
                    err = (observed - theo) / theo * 1e6
                    if abs(err) <= ppm:
                        cstr = "+".join(f"{n}{k}" for k, n in car.items())
                        out.append(Match(composition_string(comp), label, z, cstr, theo, err))
    return sorted(out, key=lambda m: abs(m.ppm_error))


def isotope_parent(mz: np.ndarray, intensity: np.ndarray, i: int, charges: Sequence[int] = (1, 2, 3),
                   ppm: float = 10.0) -> tuple[int, int] | None:
    """If peak ``i`` sits one 13C spacing above a more intense peak, return (parent index, charge)."""
    for z in charges:
        target = mz[i] - ISOTOPE_SPACING / z
        j = int(np.argmin(np.abs(mz - target)))
        if abs(mz[j] - target) <= target * ppm * 1e-6 and intensity[j] > intensity[i]:
            return j, z
    return None
