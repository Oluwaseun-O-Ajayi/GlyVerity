"""Tests for glyverity.xic."""

import numpy as np
import pytest

from glyverity.xic import summarise, xic


class _C:
    def __init__(self):
        self.xlabels = np.array([1.0, 2.0, 3.0, 4.0])
        self._s = [
            (np.array([500.000, 800.0]), np.array([1.0, 9.0])),
            (np.array([500.004, 600.0]), np.array([4.0, 9.0])),
            (np.array([499.990]), np.array([2.0])),
            (np.array([500.001]), np.array([7.0])),
        ]

    def scan(self, i):
        return self._s[i]


def test_xic_ppm_window_and_mask():
    t, y = xic(_C(), np.array([True, True, True, False]), [500.0], ppm=10)
    np.testing.assert_allclose(t, [1.0, 2.0, 3.0])
    # 500.004 is +8 ppm (in); 499.990 is -20 ppm (out)
    np.testing.assert_allclose(y[:, 0], [1.0, 4.0, 0.0])


def test_xic_time_window():
    t, _ = xic(_C(), None, [500.0], t_lo=2.0, t_hi=3.0)
    np.testing.assert_allclose(t, [2.0, 3.0])


def test_summarise():
    s = summarise(np.array([0.0, 1.0, 2.0]), np.array([0.0, 2.0, 0.0]))
    assert s.apex_time == 1.0
    assert s.apex_intensity == 2.0
    assert s.area == pytest.approx(2.0)
    assert summarise(np.array([0.0]), np.array([0.0])).area == 0.0
