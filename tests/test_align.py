"""Tests for glyverity.align."""

import numpy as np
import pytest

from glyverity.align import estimate_lag


def _peaks(t, centres, width=0.1):
    return sum(np.exp(-((t - c) / width) ** 2) for c in centres)


def test_recovers_known_lag():
    t1 = np.linspace(10, 50, 8000)
    t2 = np.linspace(10, 50, 2000)
    centres = [19.4, 23.6, 26.8, 31.7, 36.8]
    lag, r = estimate_lag(t1, _peaks(t1, centres), t2, _peaks(t2, [c + 0.12 for c in centres]),
                          window=(15, 45), max_lag=0.5, step=0.002)
    assert lag == pytest.approx(0.12, abs=0.004)
    assert 0.9 < r <= 1.0 + 1e-9


def test_bad_window():
    t = np.linspace(0, 1, 10)
    with pytest.raises(ValueError):
        estimate_lag(t, t, t, t, window=(5, 5))
