"""Target-decoy q-values for composition assignment of species neutral masses.

Decoys follow the octet-rule-violation design of Wang et al. (J. Proteome Res. 2018, 17,
2328-2334, doi:10.1021/acs.jproteome.8b00019): each target mass plus or minus an odd number of
hydrogen atoms. These formulas cannot exist as closed-shell molecules, but they keep the mass
defect of real CHNO molecules, so a decoy library samples the same region of the mass axis as
the targets. Decoy masses that fall within the search tolerance of any target are removed.

Decoy counts are scaled by target library size over decoy library size, so removing decoys
does not make the estimate optimistic. FDR is estimated by separate search: for a score cutoff t (ppm error), the mean number of
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


def signed_errors(query: np.ndarray, library: np.ndarray) -> np.ndarray:
    """Signed ppm error (query minus nearest library mass) for each query mass."""
    q = np.asarray(query, float)
    lib = np.asarray(library, float)
    i = np.clip(np.searchsorted(lib, q), 1, lib.size - 1)
    lo, hi = lib[i - 1], lib[i]
    near = np.where(np.abs(q - lo) <= np.abs(q - hi), lo, hi)
    return (q - near) / near * 1e6


def recalibrate(masses: np.ndarray, errors: np.ndarray, anchors: np.ndarray,
                min_anchors: int = 5) -> tuple[np.ndarray, float, float, int]:
    """Remove a systematic mass error, linear in mass, fitted on anchor species.

    Theil-Sen fit (median of pairwise slopes, then median intercept) of signed ppm error against
    mass over the anchors; with fewer than 5 anchors only the median offset is removed. Returns
    corrected masses, the slope (ppm per kDa), the intercept at 2 kDa (ppm) and the anchor count.
    With fewer than ``min_anchors`` anchors nothing is corrected: one or two anchors would set
    the correction from the very species being scored, and that species would then match its
    own correction exactly.
    The same correction is applied to every query mass, so targets and decoys are scored alike.
    """
    m = np.asarray(masses, float)
    e = np.asarray(errors, float)
    a = np.asarray(anchors, bool)
    n = int(a.sum())
    if n < min_anchors:
        return m.copy(), 0.0, 0.0, n
    x = m[a] / 1000.0 - 2.0
    y = e[a]
    if n >= 5:
        i, j = np.triu_indices(n, 1)
        dx = x[j] - x[i]
        ok = np.abs(dx) > 0.05
        slope = float(np.median((y[j] - y[i])[ok] / dx[ok])) if ok.any() else 0.0
    else:
        slope = 0.0
    icpt = float(np.median(y - slope * x))
    corr = icpt + slope * (m / 1000.0 - 2.0)
    return m * (1 - corr * 1e-6), slope, icpt, n


def qvalues(target_err: np.ndarray, decoy_errs: Sequence[np.ndarray],
            max_ppm: float, scale: Sequence[float] | None = None) -> np.ndarray:
    """q-value for each query whose best target error is within ``max_ppm`` (nan otherwise).

    ``scale`` multiplies each library's decoy count, normally target library size over decoy
    library size, so libraries thinned by :func:`decoy_libraries` are not under-counted.
    """
    te = np.asarray(target_err, float)
    q = np.full(te.shape, np.nan)
    hit = te <= max_ppm
    if not hit.any():
        return q
    cut = np.sort(te[hit])
    n_t = np.arange(1, cut.size + 1)
    sc = np.ones(len(decoy_errs)) if scale is None else np.asarray(scale, float)
    n_d = np.mean([k * np.searchsorted(np.sort(d), cut, side="right")
                   for k, d in zip(sc, decoy_errs)], axis=0)
    fdr = np.minimum(n_d / n_t, 1.0)
    qcum = np.minimum.accumulate(fdr[::-1])[::-1]
    idx = np.searchsorted(cut, te[hit], side="right") - 1
    q[hit] = qcum[idx]
    return q


def exchangeability(target_err: np.ndarray, decoy_errs: Sequence[np.ndarray],
                    band: tuple[float, float],
                    scale: Sequence[float] | None = None) -> tuple[int, float, float]:
    """Target and mean decoy match counts with best error inside ``band`` (ppm), and their ratio.

    Only queries with no target match below the band's lower edge are counted, so true matches
    do not enter. A ratio near 1 supports exchangeability; above 1 means decoys under-sample
    random matches and q-values are too optimistic.
    """
    lo, hi = band
    te = np.asarray(target_err, float)
    free = te > lo
    n_t = int(np.sum(free & (te <= hi)))
    sc = np.ones(len(decoy_errs)) if scale is None else np.asarray(scale, float)
    n_d = float(np.mean([k * np.sum(free & (np.asarray(d) > lo) & (np.asarray(d) <= hi))
                         for k, d in zip(sc, decoy_errs)]))
    return n_t, n_d, (n_t / n_d if n_d else np.inf)
