"""Tests for glyverity.fdr on simulated mass lists."""

import numpy as np
import pytest

from glyverity.fdr import best_errors, decoy_libraries, exchangeability, odd_h_shifts, qvalues


def _sim(seed=0, n_true=200, n_null=800, ppm_sd=1.5):
    rng = np.random.default_rng(seed)
    targets = np.sort(rng.uniform(1000, 4000, 3000))
    true = rng.choice(targets, n_true) * (1 + rng.normal(0, ppm_sd, n_true) * 1e-6)
    null = rng.uniform(1000, 4000, n_null)
    return targets, np.concatenate([true, null]), n_true


def test_best_errors_nearest():
    lib = np.array([1000.0, 2000.0])
    err = best_errors(np.array([1000.001, 1999.99]), lib)
    assert err == pytest.approx([1.0, 5.0], rel=1e-3)


def test_decoys_removed_near_targets():
    t = np.array([1000.0, 1000.0 + 3 * 1.00782503207])
    libs, removed = decoy_libraries(t, [3 * 1.00782503207], ppm=10)
    assert removed[0] == 1
    assert libs[0].size == 1


def test_qvalues_track_true_fdr_and_exchangeability_near_one():
    targets, query, n_true = _sim()
    libs, _ = decoy_libraries(targets, odd_h_shifts(), ppm=30)
    te = best_errors(query, targets)
    des = [best_errors(query, lib) for lib in libs]
    q = qvalues(te, des, max_ppm=10)
    accepted = q <= 0.05
    is_null = np.arange(query.size) >= n_true
    obs_fdr = (accepted & is_null).sum() / accepted.sum()
    assert obs_fdr < 0.1
    assert accepted[:n_true].mean() > 0.8
    _, _, ratio = exchangeability(te, des, (15, 30))
    assert 0.8 < ratio < 1.25


def test_scale_restores_thinned_library_count():
    targets, query, _ = _sim(seed=1)
    full = np.sort(targets + 5 * 1.00782503207)
    thin = full[::2]
    te = best_errors(query, targets)
    q_full = qvalues(te, [best_errors(query, full)], max_ppm=10)
    q_thin = qvalues(te, [best_errors(query, thin)], max_ppm=10, scale=[full.size / thin.size])
    ok = ~np.isnan(q_full)
    assert np.nanmedian(np.abs(q_full[ok] - q_thin[ok])) < 0.02


def test_recalibrate_removes_linear_drift():
    from glyverity.fdr import recalibrate, signed_errors
    targets, _, _ = _sim(seed=2)
    rng = np.random.default_rng(2)
    true = rng.choice(targets, 40)
    drift = 4.0 + 2.0 * (true / 1000 - 2.0)  # ppm
    meas = true * (1 + (drift + rng.normal(0, 0.5, 40)) * 1e-6)
    e = signed_errors(meas, targets)
    fixed, slope, icpt, n = recalibrate(meas, e, np.ones(40, bool))
    assert n == 40
    assert slope == pytest.approx(2.0, abs=0.5)
    assert icpt == pytest.approx(4.0, abs=0.5)
    assert np.median(np.abs(signed_errors(fixed, targets))) < 0.7


def test_recalibrate_needs_enough_anchors():
    from glyverity.fdr import recalibrate
    m = np.array([2000.0, 2500.0, 3000.0])
    e = np.array([9.0, 1.0, -2.0])
    fixed, _, icpt, n = recalibrate(m, e, np.array([True, False, False]))
    assert n == 1 and icpt == 0.0 and np.array_equal(fixed, m)
