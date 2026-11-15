"""Evidence layers that explain candidates as artifacts of other candidates.

AdductIsotopeIsobar
    A candidate whose neutral mass equals another co-eluting, more intense
    candidate's mass plus k adduct exchanges (Na+, K+ or NH4+ for H+) plus n
    13C spacings is flagged as explained by that candidate.

LabelSideProductSeries
    For a reference reducing-end group (for example 2-AB) and a test group,
    the per-composition area ratio, retention shift and fluorescence are
    compared. A consistent ratio and a same-sign shift across compositions
    marks the test group as a side-product series of the reference label.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from glyverity import chemistry as ch
from glyverity.annotate import ISOTOPE_SPACING
from glyverity.candidates import Candidate
from glyverity.evidence import Evidence, EvidenceLayer, register

ADDUCT_EXCHANGE = {
    "Na": ch.CHARGE_CARRIERS["Na+"] - ch.PROTON_MASS,
    "K": ch.CHARGE_CARRIERS["K+"] - ch.PROTON_MASS,
    "NH4": ch.CHARGE_CARRIERS["NH4+"] - ch.PROTON_MASS,
}


def exchange_combos(max_total: int = 4) -> list[tuple[dict[str, int], float]]:
    """All mixtures of Na+, K+ and NH4+ replacing H+ with 1..``max_total`` exchanges in total.

    Exchanges are not limited by charge: sialic acid carboxylates form neutral Na and K salts,
    so a 2+ ion can carry more than two non-proton cations.
    """
    out = []
    names = list(ADDUCT_EXCHANGE)
    for a in range(max_total + 1):
        for b in range(max_total + 1 - a):
            for c in range(max_total + 1 - a - b):
                if a + b + c == 0:
                    continue
                counts = dict(zip(names, (a, b, c)))
                out.append(({k: v for k, v in counts.items() if v},
                            sum(n * ADDUCT_EXCHANGE[k] for k, n in counts.items())))
    return out


def _exchange_label(counts: dict[str, int]) -> str:
    return "+".join(f"{n}x{k}" for k, n in counts.items())


def relation_offsets(max_adducts: int = 4, max_isotope: int = 5,
                     min_isotope: int = -1) -> list[tuple[str, float]]:
    """Mass offsets (mixed Na/K/NH4 exchanges) + n*13C, excluding zero, with readable labels.

    Ordered from simplest (fewest exchanges and isotope steps) to most complex. Negative ``n``
    covers a weaker, lighter candidate whose own M+1 channel picks up the stronger one.
    """
    rows = []
    for n in range(min_isotope, max_isotope + 1):
        if n != 0:
            rows.append((abs(n), f"M{n:+d}", n * ISOTOPE_SPACING))
        for counts, off in exchange_combos(max_adducts):
            iso = f"M{n:+d}" if n else "M+0"
            rows.append((sum(counts.values()) + abs(n), f"{_exchange_label(counts)} {iso}",
                         off + n * ISOTOPE_SPACING))
    rows.sort(key=lambda r: r[0])
    return [(label, off) for _, label, off in rows]


@register
class AdductIsotopeIsobar(EvidenceLayer):
    """Flag candidates explained by an adduct or isotope peak of a co-eluting stronger candidate."""

    name = "adduct_isotope_isobar"
    version = "1"

    def __init__(self, ppm: float = 5.0, rt_tol: float = 0.15, max_adducts: int = 4,
                 max_isotope: int = 5, min_isotope: int = -1):
        self.ppm = ppm
        self.rt_tol = rt_tol
        self.offsets = relation_offsets(max_adducts, max_isotope, min_isotope)

    def evaluate(self, features: Sequence[Candidate]) -> list[Evidence]:
        live = [c for c in features if c.area > 0 and np.isfinite(c.apex_time)]
        live.sort(key=lambda c: -c.area)
        explained: set[str] = set()
        out = []
        for i, c in enumerate(live):
            for d in live[:i]:  # only stronger candidates can explain c
                if d.key in explained or abs(c.apex_time - d.apex_time) > self.rt_tol:
                    continue
                diff = c.neutral_mass - d.neutral_mass
                for label, off in self.offsets:
                    err = (diff - off) / c.neutral_mass * 1e6
                    if abs(err) <= self.ppm:
                        explained.add(c.key)
                        out.append(Evidence(c.key, self.name, "explained_by", float(err),
                                            {"by": d.key, "relation": label,
                                             "rt_diff": c.apex_time - d.apex_time}))
                        break
                if c.key in explained:
                    break
        return out


R = ch.RESIDUE_MASS
NEUTRAL_LOSSES = {
    "NeuAc": R["NeuAc"],
    "NeuGc": R["NeuGc"],
    "2NeuAc": 2 * R["NeuAc"],
    "Hex+NeuAc": R["Hex"] + R["NeuAc"],
    "Hex": R["Hex"],
    "dHex": R["dHex"],
    "HexNAc": R["HexNAc"],
    "Hex+HexNAc": R["Hex"] + R["HexNAc"],
    "H2O": ch.WATER,
    "NeuAc+H2O": R["NeuAc"] + ch.WATER,
}


@register
class InSourceFragment(EvidenceLayer):
    """Flag a weaker candidate that co-elutes with a stronger one and is lighter by a common
    in-source neutral loss (sialic acid, terminal Hex+NeuAc, hexose, fucose, HexNAc, water).

    Both candidates must carry the same reducing-end group, since in-source fragmentation
    removes residues from the non-reducing end. Co-elution within ``rt_tol`` is required.
    """

    name = "in_source_fragment"
    version = "1"

    def __init__(self, ppm: float = 10.0, rt_tol: float = 0.15,
                 losses: dict[str, float] | None = None):
        self.ppm = ppm
        self.rt_tol = rt_tol
        self.losses = dict(NEUTRAL_LOSSES if losses is None else losses)

    def evaluate(self, features: Sequence[Candidate]) -> list[Evidence]:
        live = [c for c in features if c.area > 0 and np.isfinite(c.apex_time) and not c.flags]
        live.sort(key=lambda c: -c.area)
        out, done = [], set()
        for i, c in enumerate(live):
            for d in live[:i]:
                if (d.key in done or d.end_group != c.end_group
                        or abs(c.apex_time - d.apex_time) > self.rt_tol):
                    continue
                loss = d.neutral_mass - c.neutral_mass
                for name, m in self.losses.items():
                    err = (loss - m) / c.neutral_mass * 1e6
                    if abs(err) <= self.ppm:
                        done.add(c.key)
                        out.append(Evidence(c.key, self.name, "fragment_of", float(err),
                                            {"by": d.key, "loss": name,
                                             "rt_diff": c.apex_time - d.apex_time}))
                        break
                if c.key in done:
                    break
        return out


@dataclass(frozen=True)
class SeriesSummary:
    """Per-group statistics of a candidate side-product series."""

    reference: str
    test: str
    n: int
    ratio_median: float
    ratio_spread: float
    rt_shift_median: float
    same_sign_fraction: float
    fluorescence_fraction: float
    n_fluorescence_pairs: int
    is_series: bool


def fld_peaks(times: np.ndarray, values: np.ndarray, min_prominence: float = 0.01,
              smooth: int = 5) -> tuple[np.ndarray, np.ndarray]:
    """Apex times and heights of FLD peaks, each counted once.

    Local maxima of the moving-average trace whose prominence (height above the higher of the
    two minima separating it from taller peaks or the trace ends) is at least
    ``min_prominence`` times the trace maximum. Heights are read from the unsmoothed trace.
    """
    t = np.asarray(times, float)
    y = np.asarray(values, float)
    ys = np.convolve(y, np.ones(smooth) / smooth, mode="same") if smooth > 1 else y
    top = float(np.max(ys)) if ys.size else 0.0
    if top <= 0:
        return np.array([]), np.array([])
    idx = [i for i in range(1, len(ys) - 1) if ys[i] >= ys[i - 1] and ys[i] > ys[i + 1]]
    keep = []
    for i in idx:
        left = ys[:i][::-1]
        hl = np.argmax(left > ys[i]) if np.any(left > ys[i]) else len(left)
        right = ys[i + 1:]
        hr = np.argmax(right > ys[i]) if np.any(right > ys[i]) else len(right)
        base = max(left[:hl].min() if hl else ys[i], right[:hr].min() if hr else ys[i])
        if ys[i] - base >= min_prominence * top:
            keep.append(i)
    keep = np.array(keep, int)
    return t[keep], y[keep]


def fld_tail_profile(times: np.ndarray, values: np.ndarray, apexes: Sequence[float],
                     heights: Sequence[float], window: float = 2.5, n_grid: int = 101,
                     min_peaks: int = 2) -> dict[int, tuple[np.ndarray, np.ndarray, int]]:
    """Median normalised FLD peak shape on each side of an apex, from isolated peaks.

    ``apexes`` and ``heights`` are FLD peaks (see :func:`fld_peaks`). A side of a peak is used
    when no other peak lies within ``window`` min on that side. Returns
    {+1: (dt, profile, n_peaks), -1: (...)} for the trailing and leading side; a side with
    fewer than ``min_peaks`` isolated peaks is left out. Profiles are made non-increasing.
    """
    dt = np.linspace(0.0, window, n_grid)
    ap = np.asarray(apexes, float)
    out: dict[int, tuple[np.ndarray, np.ndarray, int]] = {}
    for side in (1, -1):
        rows = []
        for a, h in zip(ap, heights):
            if not (h > 0 and np.isfinite(a)):
                continue
            ahead = side * (ap - a)
            if np.any((ahead > 0) & (ahead <= window)):
                continue
            prof = np.interp(a + side * dt, times, values) / h
            rows.append(np.minimum.accumulate(np.clip(prof, 0.0, None)))
        if len(rows) >= min_peaks:
            out[side] = (dt, np.median(np.vstack(rows), axis=0), len(rows))
    return out


def predicted_tail(t: float, apexes: Sequence[float], heights: Sequence[float],
                   profile: dict[int, tuple[np.ndarray, np.ndarray, int]],
                   own: float = 0.1) -> float:
    """FLD expected at time ``t`` from the tails of FLD peaks more than ``own`` min away."""
    total = 0.0
    for a, h in zip(apexes, heights):
        if not (h > 0 and np.isfinite(a)) or abs(t - a) <= own:
            continue
        side = 1 if t > a else -1
        if side not in profile:
            continue
        dt, prof, _ = profile[side]
        d = abs(t - a)
        if d <= dt[-1]:
            total += h * float(np.interp(d, dt, prof))
    return total


@register
class LabelSideProductSeries(EvidenceLayer):
    """Detect a reducing-end side-product series relative to a fluorescent reference label.

    Pairs are formed per elemental composition when both forms have signal, the test apex is
    within ``max_shift`` min of the reference apex, and the reference area is at least
    ``min_ref_fraction`` of the strongest reference candidate. Consistency is judged with
    robust statistics: the spread of the area ratio is the median absolute deviation of
    log10(ratio) (0.1 is about a 26% spread), and the retention shifts must share a sign in
    at least ``min_same_sign`` of pairs. Fluorescence fraction is the median of observed FLD at
    the test apex over the FLD expected if the test form fluoresced like the reference.
    """

    name = "label_side_product_series"
    version = "2"

    def __init__(self, reference: str, test: str, min_n: int = 3, max_spread: float = 0.1,
                 max_shift: float = 2.0, min_ref_fraction: float = 0.01,
                 min_same_sign: float = 0.8, max_fluorescence_fraction: float = 0.5,
                 fld_clearance: float = 0.3,
                 fld_trace: tuple[np.ndarray, np.ndarray] | None = None,
                 fld_own: float = 0.1):
        self.reference = reference
        self.test = test
        self.min_n = min_n
        self.max_spread = max_spread
        self.max_shift = max_shift
        self.min_ref_fraction = min_ref_fraction
        self.min_same_sign = min_same_sign
        self.max_ff = max_fluorescence_fraction
        self.fld_clearance = fld_clearance
        self.fld_trace = fld_trace
        self.fld_own = fld_own
        self.profile: dict[int, tuple[np.ndarray, np.ndarray, int]] = {}
        self.n_fld_peaks = 0
        self.summary: SeriesSummary | None = None
        self.pairs: list[dict] = []

    def evaluate(self, features: Sequence[Candidate]) -> list[Evidence]:
        def strongest(group):
            best: dict[str, Candidate] = {}
            for c in features:
                if c.end_group == group and c.area > 0 and (
                        c.name not in best or c.area > best[c.name].area):
                    best[c.name] = c
            return best

        ref = strongest(self.reference)
        tst = strongest(self.test)
        top = max((c.area for c in ref.values()), default=0.0)
        pairs = [(ref[k], tst[k]) for k in sorted(ref.keys() & tst.keys())
                 if ref[k].area >= self.min_ref_fraction * top
                 and abs(tst[k].apex_time - ref[k].apex_time) <= self.max_shift]
        n = len(pairs)
        if n == 0:
            self.summary = SeriesSummary(self.reference, self.test, 0, np.nan, np.nan, np.nan,
                                         np.nan, np.nan, 0, False)
            return []
        ratios = np.array([t.area / r.area for r, t in pairs])
        shifts = np.array([t.apex_time - r.apex_time for r, t in pairs])
        logr = np.log10(ratios)
        spread = float(np.median(np.abs(logr - np.median(logr))))
        med_shift = float(np.median(shifts))
        same = float(np.mean(np.sign(shifts) == np.sign(med_shift))) if med_shift != 0 else 0.0
        refs = [c for c in features
                if c.end_group == self.reference and c.area > 0 and np.isfinite(c.apex_time)]
        ref_apexes = [c.apex_time for c in refs]
        self.profile = {}
        if self.fld_trace is not None:
            pk_t, pk_h = fld_peaks(self.fld_trace[0], self.fld_trace[1])
            self.n_fld_peaks = len(pk_t)
            self.profile = fld_tail_profile(self.fld_trace[0], self.fld_trace[1], pk_t, pk_h)
        ff = []
        self.pairs = []
        for r, t in pairs:
            row = {"name": r.name, "ref_rt": r.apex_time, "test_rt": t.apex_time,
                   "ratio": t.area / r.area, "fld_ref": r.fld_at_apex, "fld_test": t.fld_at_apex,
                   "fld_expected": r.fld_at_apex * t.area / r.area, "fld_tail": np.nan,
                   "used": False, "model_fail": False}
            self.pairs.append(row)
            if not (r.fld_at_apex > 0 and np.isfinite(t.fld_at_apex)):
                continue
            if self.profile:
                tail = predicted_tail(t.apex_time, pk_t, pk_h, self.profile, self.fld_own)
                row["fld_tail"] = tail
                if tail > 1.1 * t.fld_at_apex:
                    row["model_fail"] = True
                    continue
                if tail > row["fld_expected"]:
                    continue
                row["used"] = True
                ff.append(max(t.fld_at_apex - tail, 0.0) / row["fld_expected"])
            else:
                if any(abs(t.apex_time - a) < self.fld_clearance for a in ref_apexes):
                    continue
                row["used"] = True
                ff.append(t.fld_at_apex / row["fld_expected"])
        ffm = float(np.median(ff)) if ff else np.nan
        consistent = n >= self.min_n and spread <= self.max_spread and same >= self.min_same_sign
        fluorescence_ok = (not ff) or ffm <= self.max_ff
        is_series = bool(consistent and fluorescence_ok)
        self.summary = SeriesSummary(self.reference, self.test, n, float(np.median(ratios)),
                                     spread, med_shift, same, ffm, len(ff), is_series)
        if not is_series:
            return []
        return [Evidence(t.key, self.name, "side_product_of", float(t.area / r.area),
                         {"reference": r.key, "rt_shift": t.apex_time - r.apex_time})
                for r, t in pairs]


@dataclass(frozen=True)
class RecurringDelta:
    """A neutral mass difference seen between several co-eluting candidate pairs."""

    delta: float
    count: int
    pairs: tuple[tuple[str, str], ...]


def recurring_deltas(candidates: Sequence[Candidate], rt_tol: float = 0.15, tol_da: float = 0.005,
                     min_count: int = 2, min_abs_delta: float = 0.5) -> list[RecurringDelta]:
    """Cluster neutral mass differences (weaker minus stronger) between co-eluting candidates.

    Only unflagged candidates with signal are used. Differences closer than ``tol_da`` join a
    cluster; clusters with at least ``min_count`` pairs are returned, most frequent first. This
    reports relations for inspection and flags nothing.
    """
    live = [c for c in candidates if c.area > 0 and np.isfinite(c.apex_time) and not c.flags]
    live.sort(key=lambda c: -c.area)
    diffs = []
    for i, c in enumerate(live):
        for d in live[:i]:
            if abs(c.apex_time - d.apex_time) <= rt_tol:
                delta = c.neutral_mass - d.neutral_mass
                if abs(delta) >= min_abs_delta:
                    diffs.append((delta, (c.key, d.key)))
    diffs.sort()
    clusters: list[list[tuple[float, tuple[str, str]]]] = []
    for item in diffs:
        if clusters and item[0] - clusters[-1][-1][0] <= tol_da:
            clusters[-1].append(item)
        else:
            clusters.append([item])
    out = [RecurringDelta(float(np.mean([x[0] for x in cl])), len(cl), tuple(x[1] for x in cl))
           for cl in clusters if len(cl) >= min_count]
    return sorted(out, key=lambda r: (-r.count, abs(r.delta)))
