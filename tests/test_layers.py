"""Tests for glyverity.layers, using the mass relations seen in the fetuin run."""

import numpy as np
import pytest

from glyverity import chemistry as ch
from glyverity.candidates import Candidate
from glyverity.layers import AdductIsotopeIsobar, LabelSideProductSeries, relation_offsets

TAG = ch.formula_mass({"C": 3, "H": 9, "N": 1}) - ch.ATOMIC_MASS["O"]


def _cand(comp, group, rt, area, fld=np.nan):
    added = {"2-AB": ch.LABELS["2-AB"].added_mass, "free": 0.0, "tag": TAG}[group]
    return Candidate(comp, group, ch.composition_mass(comp) + ch.WATER + added, rt, area, area, fld)


def test_relation_offsets_contains_na_m4():
    offs = dict(relation_offsets())
    assert offs["1xNa M+4"] == pytest.approx(25.995, abs=1e-3)


def test_free_fucosylated_explained_by_na_isotope_of_2ab():
    base = {"Hex": 6, "HexNAc": 5, "NeuAc": 2}
    real = _cand(base, "2-AB", 31.70, 2.0e5)
    fake = _cand(dict(base, dHex=1), "free", 31.66, 1.2e4)
    ev = AdductIsotopeIsobar(ppm=5, rt_tol=0.15).evaluate([real, fake])
    assert len(ev) == 1
    assert ev[0].feature_id == fake.key
    assert ev[0].detail["by"] == real.key
    assert ev[0].detail["relation"] == "1xNa M+4"


def test_isobar_needs_coelution():
    base = {"Hex": 6, "HexNAc": 5, "NeuAc": 2}
    real = _cand(base, "2-AB", 31.70, 2.0e5)
    other = _cand(dict(base, dHex=1), "free", 34.0, 1.2e4)
    assert AdductIsotopeIsobar(rt_tol=0.15).evaluate([real, other]) == []


def test_side_product_series_detected():
    comps = [{"Hex": 6, "HexNAc": 5, "NeuAc": s} for s in range(1, 5)]
    cands = []
    for i, c in enumerate(comps):
        cands.append(_cand(c, "2-AB", 26 + 3 * i, 1e5 * (i + 1), fld=2.0))
        cands.append(_cand(c, "tag", 25 + 3 * i, 0.43e5 * (i + 1), fld=0.1))
    layer = LabelSideProductSeries("2-AB", "tag")
    ev = layer.evaluate(cands)
    assert layer.summary.is_series
    assert layer.summary.ratio_spread == pytest.approx(0.0, abs=1e-9)
    assert layer.summary.ratio_median == pytest.approx(0.43)
    assert layer.summary.rt_shift_median == pytest.approx(-1.0)
    assert len(ev) == 4


def test_side_product_rejected_when_fluorescent():
    comps = [{"Hex": 6, "HexNAc": 5, "NeuAc": s} for s in range(1, 5)]
    cands = []
    for i, c in enumerate(comps):
        cands.append(_cand(c, "2-AB", 26 + 3 * i, 1e5, fld=2.0))
        cands.append(_cand(c, "tag", 25 + 3 * i, 0.4e5, fld=0.8))  # fluorescence as expected
    layer = LabelSideProductSeries("2-AB", "tag")
    assert layer.evaluate(cands) == []
    assert not layer.summary.is_series


class _Scans:
    """Two scans; the 2-AB A3S2 [M+2H]2+ ion is present in the second only."""

    def __init__(self):
        from glyverity.candidates import target_mzs
        m = ch.composition_mass({"Hex": 6, "HexNAc": 5, "NeuAc": 2}) + ch.WATER + ch.LABELS["2-AB"].added_mass
        self.mz2 = target_mzs(m, (2,), 1)[0]
        self.xlabels = np.array([30.0, 31.0])

    def scan(self, i):
        if i == 0:
            return np.array([500.0]), np.array([5.0])
        return np.array([500.0, self.mz2]), np.array([5.0, 100.0])


def test_extract_candidates_single_pass():
    from glyverity.candidates import extract_candidates
    comps = [{"Hex": 6, "HexNAc": 5, "NeuAc": 2}, {"Hex": 5, "HexNAc": 4}]
    groups = {"2-AB": ch.LABELS["2-AB"].added_mass, "free": 0.0}
    cands = extract_candidates(_Scans(), None, comps, groups, (29.0, 32.0))
    found = {c.key: c for c in cands}
    assert found["Hex6HexNAc5NeuAc2|2-AB"].apex_time == 31.0
    assert found["Hex6HexNAc5NeuAc2|2-AB"].apex_intensity == pytest.approx(100.0)
    assert found["Hex5HexNAc4|free"].area == 0.0
    assert len(cands) == 4


def test_isomeric_compositions_merged():
    from glyverity.candidates import group_isomeric
    a = {"Hex": 6, "HexNAc": 5, "NeuAc": 3}
    b = {"Hex": 5, "HexNAc": 5, "dHex": 1, "NeuAc": 2, "NeuGc": 1}
    groups = group_isomeric([b, a, {"Hex": 5, "HexNAc": 4}])
    assert len(groups) == 2
    merged = next(g for g in groups if len(g) == 2)
    assert merged[0] == a  # display choice: fewer NeuGc


def test_series_ignores_far_and_weak_pairs():
    comps = [{"Hex": 6, "HexNAc": 5, "NeuAc": s} for s in range(1, 5)]
    cands = []
    for i, c in enumerate(comps):
        cands.append(_cand(c, "2-AB", 26 + 3 * i, 1e5, fld=2.0))
        cands.append(_cand(c, "tag", 25 + 3 * i, 0.43e5, fld=0.1))
    # noise: a test form 6 min away and a negligible reference
    noise = {"Hex": 3, "HexNAc": 6}
    cands.append(_cand(noise, "2-AB", 20.0, 50.0, fld=0.1))
    cands.append(_cand(noise, "tag", 26.0, 5e4, fld=0.0))
    layer = LabelSideProductSeries("2-AB", "tag")
    layer.evaluate(cands)
    assert layer.summary.n == 4
    assert layer.summary.is_series
