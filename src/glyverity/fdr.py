"""Target-decoy q-values for composition assignment of species neutral masses.

Decoys follow the octet-rule-violation design of Wang et al. (J. Proteome Res. 2018, 17,
2328-2334, doi:10.1021/acs.jproteome.8b00019): each target mass plus or minus an odd number of
hydrogen atoms. These formulas cannot exist as closed-shell molecules, but they keep the mass
defect of real CHNO molecules, so a decoy library samples the same region of the mass axis as
the targets. Decoy masses that fall within the search tolerance of any target are removed.

FDR is estimated by separate search: for a score cutoff t (ppm error), the mean number of
species whose best decoy match is within t, over the number whose best target match is within
t, across the decoy libraries. q-values are the running minimum of FDR from the largest cutoff
down. Exchangeability (decoys behave like random target matches) is checked in an error band
where true matches are rare: there the target and decoy match counts should be equal.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from glyverity import chemistry as ch

H = ch.ATOMIC_MASS["H"]


def odd_h_shifts(counts: Sequence[int] = (3, 5, 7, 9), both_signs: bool = True) -> list[float]:
    """Decoy mass shifts of odd numbers of hydrogen atoms (Da)."""
    out = [k * H for k in counts]
    if both_signs:
        out += [-k * H for k in counts]
    return out


def decoy_libraries(targets: np.ndarray, shifts: Sequence[float],
                    ppm: float) -> tuple[list[np.ndarray], list[int]]:
    """One sorted decoy library per shift, with decoys within ``ppm`` of any target removed.

    Returns the libraries and the number removed from each.
    """
    t = np.sort(np.asarray(targets, float))
    libs, removed = [], []
    for d in shifts:
        dec = t + d
        err = best_errors(dec, t)
        keep = err > ppm
        libs.append(np.sort(dec[keep]))
        removed.append(int((~keep).sum()))
    return libs, removed


def best_errors(query: np.ndarray, library: np.ndarray) -> np.ndarray:
    """Absolute ppm error from each query mass to the nearest mass in a sorted library."""
    q = np.asarray(query, float)
    lib = np.asarray(library, float)
    if lib.size == 0:
        return np.full(q.shape, np.inf)
    i = np.clip(np.searchsorted(lib, q), 1, lib.size - 1) if lib.size > 1 else np.zeros(q.shape, int)
    lo = lib[np.maximum(i - 1, 0)]
    hi = lib[i]
    d = np.minimum(np.abs(q - lo), np.abs(q - hi))
    return d / q * 1e6


def qvalues(target_err: np.ndarray, decoy_errs: Sequence[np.ndarray],
            max_ppm: float) -> np.ndarray:
    """q-value for each query whose best target error is within ``max_ppm`` (nan otherwise)."""
    te = np.asarray(target_err, float)
    q = np.full(te.shape, np.nan)
    hit = te <= max_ppm
    if not hit.any():
        return q
    cut = np.sort(te[hit])
    n_t = np.arange(1, cut.size + 1)
    n_d = np.mean([np.searchsorted(np.sort(d), cut, side="right") for d in decoy_errs], axis=0)
    fdr = np.minimum(n_d / n_t, 1.0)
    qcum = np.minimum.accumulate(fdr[::-1])[::-1]
    idx = np.searchsorted(cut, te[hit], side="right") - 1
    q[hit] = qcum[idx]
    return q


def exchangeability(target_err: np.ndarray, decoy_errs: Sequence[np.ndarray],
                    band: tuple[float, float]) -> tuple[int, float, float]:
    """Target and mean decoy match counts with best error inside ``band`` (ppm), and their ratio.

    Only queries with no target match below the band's lower edge are counted, so true matches
    do not enter. A ratio near 1 supports exchangeability; above 1 means decoys under-sample
    random matches and q-values are too optimistic.
    """
    lo, hi = band
    te = np.asarray(target_err, float)
    free = te > lo
    n_t = int(np.sum(free & (te <= hi)))
    n_d = float(np.mean([np.sum(free & (np.asarray(d) > lo) & (np.asarray(d) <= hi))
                         for d in decoy_errs]))
    return n_t, n_d, (n_t / n_d if n_d else np.inf)
