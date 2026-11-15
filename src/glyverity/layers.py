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


def relation_offsets(max_adducts: int = 2, max_isotope: int = 5) -> list[tuple[str, float]]:
    """Mass offsets k*(adduct - H) + n*13C, excluding the zero offset, with readable labels."""
    out = []
    for name, d in ADDUCT_EXCHANGE.items():
        for k in range(max_adducts + 1):
            for n in range(max_isotope + 1):
                if k == 0 and n == 0:
                    continue
                if k == 0 and name != "Na":
                    continue  # pure isotope offsets listed once
                label = f"{k}x{name} M+{n}" if k else f"M+{n}"
                out.append((label, k * d + n * ISOTOPE_SPACING))
    return out


@register
class AdductIsotopeIsobar(EvidenceLayer):
    """Flag candidates explained by an adduct or isotope peak of a co-eluting stronger candidate."""

    name = "adduct_isotope_isobar"
    version = "1"

    def __init__(self, ppm: float = 5.0, rt_tol: float = 0.15, max_adducts: int = 2,
                 max_isotope: int = 5):
        self.ppm = ppm
        self.rt_tol = rt_tol
        self.offsets = relation_offsets(max_adducts, max_isotope)

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
    is_series: bool


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
                 min_same_sign: float = 0.8, max_fluorescence_fraction: float = 0.5):
        self.reference = reference
        self.test = test
        self.min_n = min_n
        self.max_spread = max_spread
        self.max_shift = max_shift
        self.min_ref_fraction = min_ref_fraction
        self.min_same_sign = min_same_sign
        self.max_ff = max_fluorescence_fraction
        self.summary: SeriesSummary | None = None

    def evaluate(self, features: Sequence[Candidate]) -> list[Evidence]:
        ref = {c.name: c for c in features if c.end_group == self.reference and c.area > 0}
        tst = {c.name: c for c in features if c.end_group == self.test and c.area > 0}
        top = max((c.area for c in ref.values()), default=0.0)
        pairs = [(ref[k], tst[k]) for k in sorted(ref.keys() & tst.keys())
                 if ref[k].area >= self.min_ref_fraction * top
                 and abs(tst[k].apex_time - ref[k].apex_time) <= self.max_shift]
        n = len(pairs)
        if n == 0:
            self.summary = SeriesSummary(self.reference, self.test, 0, np.nan, np.nan, np.nan,
                                         np.nan, np.nan, False)
            return []
        ratios = np.array([t.area / r.area for r, t in pairs])
        shifts = np.array([t.apex_time - r.apex_time for r, t in pairs])
        logr = np.log10(ratios)
        spread = float(np.median(np.abs(logr - np.median(logr))))
        med_shift = float(np.median(shifts))
        same = float(np.mean(np.sign(shifts) == np.sign(med_shift))) if med_shift != 0 else 0.0
        ff = np.array([
            t.fld_at_apex / (r.fld_at_apex * (t.area / r.area))
            if r.fld_at_apex > 0 and np.isfinite(t.fld_at_apex) else np.nan
            for r, t in pairs
        ])
        ffm = float(np.nanmedian(ff)) if np.any(np.isfinite(ff)) else np.nan
        is_series = (n >= self.min_n and spread <= self.max_spread and same >= self.min_same_sign
                     and np.isfinite(ffm) and ffm <= self.max_ff)
        self.summary = SeriesSummary(self.reference, self.test, n, float(np.median(ratios)),
                                     spread, med_shift, same, ffm, is_series)
        if not is_series:
            return []
        return [Evidence(t.key, self.name, "side_product_of", float(t.area / r.area),
                         {"reference": r.key, "rt_shift": t.apex_time - r.apex_time})
                for r, t in pairs]
