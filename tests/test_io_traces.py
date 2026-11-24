"""Tests for glyverity.io_traces."""

import numpy as np

from glyverity.io_traces import read_trace_csv


def test_read_trace_csv(tmp_path):
    p = tmp_path / "x.csv"
    p.write_text("time_min,Emission_1,other\n0.2,5,1\n0.1,4,1\n0.3,,1\n0.4,6,1\n")
    tr = read_trace_csv(p)
    assert tr.name == "Emission_1"
    assert np.allclose(tr.times, [0.1, 0.2, 0.4])
    assert np.allclose(tr.values, [4, 5, 6])
    assert read_trace_csv(p, "other").name == "other"
