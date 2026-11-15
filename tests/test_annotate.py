"""Tests for glyverity.annotate."""

import numpy as np
import pytest

from glyverity import annotate as an
from glyverity import chemistry as ch


def test_merge_peaks_groups_within_ppm():
    mz = np.array([1000.000, 1000.005, 1200.0])
    inten = np.array([1.0, 3.0, 2.0])
    m, y = an.merge_peaks(mz, inten, ppm=10)
    assert m.size == 2
    assert m[0] == pytest.approx(1000.00375)
    np.testing.assert_allclose(y, [4.0, 2.0])


def test_enumerate_and_string():
    comps = an.enumerate_compositions({"Hex": (0, 1), "HexNAc": (0, 1)})
    assert len(comps) == 3
    assert an.composition_string({"HexNAc": 4, "Hex": 5, "NeuAc": 2}) == "Hex5HexNAc4NeuAc2"


def test_match_finds_labelled_sialylated_biantennary():
    comp = {"Hex": 5, "HexNAc": 4, "NeuAc": 2}
    obs = ch.mz(ch.glycan_neutral_mass(comp, "2-AB"), 2)
    hits = an.match_mz(obs, [comp, {"Hex": 5, "HexNAc": 2}], ["free", "2-AB"], ppm=5)
    assert hits[0].composition == "Hex5HexNAc4NeuAc2"
    assert hits[0].label == "2-AB"
    assert hits[0].charge == 2
    assert hits[0].ppm_error == pytest.approx(0.0, abs=1e-6)


def test_isotope_parent():
    mz = np.array([1000.0, 1000.0 + an.ISOTOPE_SPACING / 2])
    inten = np.array([10.0, 6.0])
    assert an.isotope_parent(mz, inten, 1) == (0, 2)
    assert an.isotope_parent(mz, inten, 0) is None


class _C:
    def __init__(self):
        self.xlabels = np.array([1.0, 2.0, 3.0])
        self._s = [(np.array([500.0]), np.array([1.0])),
                   (np.array([500.001]), np.array([3.0])),
                   (np.array([700.0]), np.array([5.0]))]

    def scan(self, i):
        return self._s[i]


def test_window_spectrum_respects_window_and_mask():
    m, y, n = an.window_spectrum(_C(), np.array([True, True, False]), 0.5, 3.5)
    assert n == 2
    assert m.size == 1
    assert y[0] == pytest.approx(4.0)
