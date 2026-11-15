"""Tests for glyverity.io_agilent using stand-in objects shaped like rainbow-api data files."""

import numpy as np
import pytest

from glyverity import io_agilent as ia


class _Trace:
    def __init__(self, name, detector, times, values, metadata=None):
        self.name = name
        self.detector = detector
        self.xlabels = np.asarray(times)
        self.data = np.asarray(values).reshape(-1, 1)
        self.metadata = metadata or {}


class _Centroid:
    def __init__(self, times, scans):
        self.name = "MSPeak.bin"
        self.detector = "MS"
        self.xlabels = np.asarray(times)
        self._scans = scans
        self.metadata = {}

    def scan(self, i):
        return self._scans[i]


def test_fld_selected_and_flattened():
    files = [
        _Trace("FLD1A.cg", "FLD", [0.0, 0.1, 0.2], [1.0, 5.0, 2.0], {"unit": "EU"}),
        _Trace("DAD1A.cg", "UV", [0.0, 0.1], [0.3, 0.4]),
    ]
    out = ia.fld_traces(files)
    assert len(out) == 1
    assert out[0].name == "FLD1A.cg"
    np.testing.assert_allclose(out[0].values, [1.0, 5.0, 2.0])
    assert out[0].metadata["unit"] == "EU"


def test_tic_sums_each_scan():
    c = _Centroid([1.0, 2.0], [
        (np.array([500.0, 600.0]), np.array([10.0, 30.0])),
        (np.array([700.0]), np.array([5.0])),
    ])
    tic = ia.ms_tic(c)
    np.testing.assert_allclose(tic.times, [1.0, 2.0])
    np.testing.assert_allclose(tic.values, [40.0, 5.0])


def test_inventory_reports_ranges():
    items = ia.inventory([_Trace("FLD1A.cg", "FLD", [0.5, 1.0, 30.0], [0, 0, 0])])
    assert items[0].n_points == 3
    assert items[0].t_start == pytest.approx(0.5)
    assert items[0].t_end == pytest.approx(30.0)
