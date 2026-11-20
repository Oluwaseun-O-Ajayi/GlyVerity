"""Tests for glyverity.validate on synthetic FLD traces."""

import numpy as np
import pytest

from glyverity.validate import fld_area, fld_window, log_ratio_slope, total_variation, windows_for


def _trace(peaks, sigma=0.08):
    t = np.linspace(10, 40, 9001)
    y = np.full_like(t, 0.05)
    for a, h in peaks:
        y += h * np.exp(-0.5 * ((t - a) / sigma) ** 2)
    return t, y


def test_window_and_area_of_isolated_peak():
    t, y = _trace([(20.0, 2.0)])
    lo, hi = fld_window(t, y, 20.03)
    assert lo < 19.7 and hi > 20.3
    assert fld_area(t, y, lo, hi) == pytest.approx(2.0 * 0.08 * np.sqrt(2 * np.pi), rel=0.03)


def test_coeluting_species_share_a_window_and_neighbours_split():
    t, y = _trace([(20.0, 2.0), (21.0, 1.0)])
    ws = windows_for(t, y, [20.0, 20.02, 21.0])
    assert len(ws) == 2
    assert ws[0].members == [0, 1]
    assert ws[1].area == pytest.approx(0.5 * ws[0].area, rel=0.05)


def test_total_variation_and_slope():
    assert total_variation([1, 1], [1, 1]) == 0.0
    assert total_variation([1, 0], [0, 1]) == 1.0
    fld = [1.0, 1.0, 1.0]
    ms = [1.0, 0.5, 0.25]  # each extra sialic acid halves the protonated MS response
    assert log_ratio_slope(ms, fld, [1, 2, 3]) == pytest.approx(np.log10(0.5))
