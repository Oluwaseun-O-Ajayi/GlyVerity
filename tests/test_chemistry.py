"""Tests for glyverity.chemistry against published reference masses."""

import pytest

from glyverity import chemistry as ch


def test_water_mass():
    assert ch.WATER == pytest.approx(18.010565, abs=1e-6)


@pytest.mark.parametrize(
    "residue, expected",
    [
        ("Hex", 162.052823),
        ("HexNAc", 203.079373),
        ("dHex", 146.057909),
        ("NeuAc", 291.095417),
        ("NeuGc", 307.090331),
    ],
)
def test_residue_masses(residue, expected):
    assert ch.RESIDUE_MASS[residue] == pytest.approx(expected, abs=1e-6)


def test_man5_free_neutral_mass():
    # Man5 (Hex5HexNAc2), free reducing end: 1234.4334
    assert ch.glycan_neutral_mass({"Hex": 5, "HexNAc": 2}) == pytest.approx(1234.4334, abs=1e-4)


def test_label_added_masses():
    # 2-AB adds 120.069 Da and procainamide 219.173 Da (Waters application note 720004212).
    assert ch.LABELS["2-AB"].added_mass == pytest.approx(120.0687, abs=1e-4)
    assert ch.LABELS["ProA"].added_mass == pytest.approx(219.1735, abs=1e-4)


def test_procainamide_acetonitrile_series():
    # Procainamide plus acetonitrile gives the +260.200 Da side series seen in GPST000788.
    total = ch.LABELS["ProA"].added_mass + ch.ARTIFACT_SHIFTS["acetonitrile"]
    assert total == pytest.approx(260.2001, abs=1e-4)


def test_mz_protonated():
    m = ch.glycan_neutral_mass({"Hex": 5, "HexNAc": 2})
    assert ch.mz(m, 1) == pytest.approx(m + ch.PROTON_MASS)
    assert ch.mz(m, 2) == pytest.approx((m + 2 * ch.PROTON_MASS) / 2)


def test_mz_mixed_carriers():
    m = 1000.0
    expected = (m + ch.PROTON_MASS + ch.CHARGE_CARRIERS["Na+"]) / 2
    assert ch.mz(m, 2, {"H+": 1, "Na+": 1}) == pytest.approx(expected)


def test_mz_rejects_bad_input():
    with pytest.raises(ValueError):
        ch.mz(1000.0, 0)
    with pytest.raises(ValueError):
        ch.mz(1000.0, 2, {"H+": 1})


def test_negative_count_rejected():
    with pytest.raises(ValueError):
        ch.composition_mass({"Hex": -1})
