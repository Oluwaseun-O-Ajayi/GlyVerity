"""Tests for glyverity.features on a synthetic run with known glycan envelopes."""

import numpy as np
import pytest

pytest.importorskip("pyopenms")

from glyverity import chemistry as ch
from glyverity.features import find_features, to_msexperiment


def _envelope(formula, z):
    import pyopenms as oms
    ef = oms.EmpiricalFormula(ch.formula_string(formula))
    iso = ef.getIsotopeDistribution(oms.CoarseIsotopePatternGenerator(6))
    neutral = ef.getMonoWeight()
    mz = np.array([(neutral + k * 1.0033548378 + z * ch.PROTON_MASS) / z
                   for k in range(len(iso.getContainer()))])
    ab = np.array([p.getIntensity() for p in iso.getContainer()])
    return mz, ab / ab.max()


def _synthetic():
    a3s2 = ch.glycan_formula({"Hex": 6, "HexNAc": 5, "NeuAc": 2}, ch.label_net_formula("2-AB"))
    a2s2 = ch.glycan_formula({"Hex": 5, "HexNAc": 4, "NeuAc": 2}, ch.label_net_formula("2-AB"))
    species = [(a3s2, 2, 31.7, 1e6), (a3s2, 3, 31.7, 6e5), (a2s2, 2, 23.6, 4e5)]
    rng = np.random.default_rng(0)
    scans = []
    for t in np.arange(20.0, 34.0, 0.017):  # about 1 s per scan
        mzs, its = [rng.uniform(400, 2000, 30)], [rng.uniform(50, 300, 30)]
        for f, z, rt, amp in species:
            mz, ab = _envelope(f, z)
            h = amp * np.exp(-((t - rt) / 0.08) ** 2)
            if h > 100:
                mzs.append(mz)
                its.append(ab * h)
        scans.append((t, np.concatenate(mzs), np.concatenate(its)))
    return scans


def test_finds_known_envelopes():
    feats = find_features(to_msexperiment(_synthetic()), ppm=10, noise=20)
    m_a3s2 = ch.glycan_neutral_mass({"Hex": 6, "HexNAc": 5, "NeuAc": 2}, "2-AB")
    m_a2s2 = ch.glycan_neutral_mass({"Hex": 5, "HexNAc": 4, "NeuAc": 2}, "2-AB")

    def hit(mass, z, rt):
        return [f for f in feats if f.charge == z and abs(f.neutral_mass - mass) / mass < 5e-6
                and abs(f.rt - rt) < 0.1]

    assert hit(m_a3s2, 2, 31.7)
    assert hit(m_a3s2, 3, 31.7)
    assert hit(m_a2s2, 2, 23.6)


def test_label_net_formula_matches_mass():
    net = ch.label_net_formula("2-AB")
    assert ch.formula_mass(net) == pytest.approx(ch.LABELS["2-AB"].added_mass, abs=1e-9)
    f = ch.glycan_formula({"Hex": 5, "HexNAc": 2})
    assert ch.formula_mass(f) == pytest.approx(ch.glycan_neutral_mass({"Hex": 5, "HexNAc": 2}))
