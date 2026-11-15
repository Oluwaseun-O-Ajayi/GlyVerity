"""Tests for glyverity.census on hand-built features."""

import pytest

from glyverity import chemistry as ch
from glyverity.census import assign_species, classify, group_species, neutral_for
from glyverity.features import MSFeature

TAG = ch.formula_mass({"C": 3, "H": 9, "N": 1}) - ch.ATOMIC_MASS["O"]
GROUPS = {"2-AB": ch.LABELS["2-AB"].added_mass, "free": 0.0, "tag": TAG}
A3S2 = {"Hex": 6, "HexNAc": 5, "NeuAc": 2}


def _f(fid, neutral, z, rt, inten, carriers=None):
    return MSFeature(fid, ch.mz(neutral, z, carriers), z, rt, inten, 4)


def test_charge_states_and_adduct_grouped():
    m = ch.glycan_neutral_mass(A3S2, "2-AB")
    feats = [_f(0, m, 2, 31.70, 5e6), _f(1, m, 3, 31.71, 3e6),
             _f(2, m, 2, 31.70, 8e5, {"H+": 1, "Na+": 1}), _f(3, 1500.0, 2, 20.0, 1e5)]
    sp = group_species(feats)
    assert len(sp) == 2
    main = max(sp, key=lambda s: s.intensity)
    assert sorted(form for _, form in main.members) == ["1Na", "H", "H"]
    assert main.neutral_mass == pytest.approx(m, rel=1e-9)


def test_sodium_adduct_with_two_charge_states_absorbed():
    m = ch.glycan_neutral_mass(A3S2, "2-AB")
    feats = [_f(0, m, 2, 31.75, 1.7e6), _f(1, m, 3, 31.70, 1.7e6),
             _f(2, m, 2, 31.70, 9.7e5, {"H+": 1, "Na+": 1}),
             _f(3, m, 3, 31.70, 3e5, {"H+": 2, "Na+": 1})]
    sp = group_species(feats)
    assert len(sp) == 1
    assert sorted(form for _, form in sp[0].members) == ["1Na", "1Na", "H", "H"]


def test_charge_states_0p14_min_apart_grouped():
    m = ch.glycan_neutral_mass({"Hex": 6, "HexNAc": 5, "NeuAc": 3}, "2-AB")
    sp = group_species([_f(0, m, 3, 36.837, 2.2e6), _f(1, m, 2, 36.695, 1.1e6)])
    assert len(sp) == 1


def test_neutral_for_adduct():
    m = 2000.0
    f = _f(0, m, 3, 10.0, 1.0, {"H+": 2, "Na+": 1})
    assert neutral_for(f, "Na+", 1) == pytest.approx(m)


def test_assign_and_classify():
    m_ab = ch.glycan_neutral_mass(A3S2, "2-AB")
    m_tag = ch.composition_mass(A3S2) + ch.WATER + TAG
    feats = [_f(0, m_ab, 2, 31.70, 5e6), _f(1, m_ab, 2, 31.70, 8e5, {"H+": 1, "Na+": 1}),
             _f(2, m_tag, 2, 30.65, 2e6), _f(3, 1234.5678, 2, 25.0, 1e5)]
    asg = assign_species(group_species(feats), [A3S2, {"Hex": 5, "HexNAc": 4}], GROUPS)
    cls = classify(asg, "2-AB")
    assert cls[0] == "labelled glycan"
    assert cls[1] == "adduct ion of a species"
    assert cls[2] == "other end group"  # before the series layer flags it
    assert cls[3] == "unassigned"


def test_end_to_end_synthetic():
    np = pytest.importorskip("numpy")
    oms = pytest.importorskip("pyopenms")
    from glyverity.candidates import Candidate  # noqa: F401
    from glyverity.features import find_features, to_msexperiment
    from glyverity.layers import InSourceFragment

    def env(formula, z, carriers=None):
        ef = oms.EmpiricalFormula(ch.formula_string(formula))
        iso = ef.getIsotopeDistribution(oms.CoarseIsotopePatternGenerator(6)).getContainer()
        m0 = ch.mz(ef.getMonoWeight(), z, carriers)
        ab = np.array([p.getIntensity() for p in iso])
        return np.array([m0 + k * 1.0033548378 / z for k in range(len(ab))]), ab / ab.max()

    ab_net = ch.label_net_formula("2-AB")
    tag_net = {"C": 3, "H": 9, "N": 1, "O": -1}
    a3s2 = ch.glycan_formula(A3S2, ab_net)
    frag = ch.glycan_formula({"Hex": 5, "HexNAc": 5, "NeuAc": 1}, ab_net)
    tag = ch.glycan_formula(A3S2, tag_net)
    planted = [(a3s2, 2, None, 31.7, 1e6), (a3s2, 3, None, 31.7, 6e5),
               (a3s2, 2, {"H+": 1, "Na+": 1}, 31.7, 1.5e5), (frag, 2, None, 31.7, 4e4),
               (tag, 2, None, 30.65, 4e5)]
    rng = np.random.default_rng(1)
    scans = []
    for t in np.arange(29.5, 33.0, 0.017):
        mzs, its = [rng.uniform(400, 2000, 30)], [rng.uniform(50, 300, 30)]
        for f, z, car, rt, amp in planted:
            h = amp * np.exp(-((t - rt) / 0.08) ** 2)
            if h > 100:
                mz, a = env(f, z, car)
                mzs.append(mz)
                its.append(a * h)
        scans.append((t, np.concatenate(mzs), np.concatenate(its)))

    feats = find_features(to_msexperiment(scans), noise=20)
    asg = assign_species(group_species(feats), [A3S2, {"Hex": 5, "HexNAc": 5, "NeuAc": 1}],
                         GROUPS)
    cands = [c for a in asg for c in a.matches]
    by_key = {c.key: c for c in cands}
    for e in InSourceFragment(rt_tol=0.1).evaluate(cands):
        by_key[e.feature_id].flags.append("isf")
    cls = classify(asg, "2-AB")
    got = {}
    for f in feats:
        got.setdefault(cls[f.feature_id], []).append(round(f.mz, 2))
    assert len(got["labelled glycan"]) == 2          # 2+ and 3+ of A3S2 2-AB
    assert len(got["adduct ion of a species"]) == 1  # [M+H+Na]2+
    assert len(got["in-source fragment"]) == 1
    assert len(got["other end group"]) == 1          # tag form (series needs >= 3 compositions)
    assert "unassigned" not in got


def test_repair_split_envelope():
    from glyverity.annotate import ISOTOPE_SPACING
    from glyverity.census import repair_charges
    mono = 1519.5226
    feats = [MSFeature(0, mono, 1, 36.742, 1.8e5, 1),
             MSFeature(1, mono + ISOTOPE_SPACING / 2, 1, 36.742, 2.3e5, 1),
             MSFeature(2, 700.0, 1, 36.742, 1e4, 1)]
    out, n = repair_charges(feats)
    assert n == 1
    z2 = [f for f in out if f.charge == 2]
    assert len(z2) == 1 and z2[0].mz == mono and z2[0].intensity == pytest.approx(4.1e5)


def test_mixed_salt_adduct_absorbed():
    m = ch.glycan_neutral_mass({"Hex": 6, "HexNAc": 5, "NeuAc": 2}, "2-AB")
    na, k = ch.CHARGE_CARRIERS["Na+"], ch.CHARGE_CARRIERS["K+"]
    p = ch.PROTON_MASS
    salt = m + 2 * (na - p) + (k - p)  # 2 Na + 1 K replacing 3 H (one as a neutral salt)
    feats = [_f(0, m, 2, 31.70, 3e6), _f(1, salt, 2, 31.80, 2e5)]
    sp = group_species(feats)
    assert len(sp) == 1
    assert any(form.startswith("2Na+1K") for _, form in sp[0].members)


def test_oxonium_list_contains_neuac_hex_hexnac():
    from glyverity.census import oxonium_ions
    ions = dict(oxonium_ions())
    assert ions["Hex1HexNAc1NeuAc1"] == pytest.approx(657.2349, abs=1e-4)
    assert ions["HexNAc1"] == pytest.approx(204.0867, abs=1e-4)


def test_isf_with_isotope_error_is_fragment_not_adduct():
    a3s3 = {"Hex": 6, "HexNAc": 5, "NeuAc": 3}
    m = ch.glycan_neutral_mass(a3s3, "2-AB")
    frag = m - ch.RESIDUE_MASS["NeuAc"] + 1.0033548378  # -NeuAc, monoisotopic picked +1
    feats = [_f(0, m, 3, 36.84, 2.2e6), _f(1, frag, 2, 36.70, 1.4e5)]
    asg = assign_species(group_species(feats), [a3s3], GROUPS)
    cls = classify(asg, "2-AB")
    assert cls[0] == "labelled glycan"
    assert cls[1] == "in-source fragment"


def test_repair_merges_peaks_0p26_min_apart():
    from glyverity.annotate import ISOTOPE_SPACING
    from glyverity.census import repair_charges
    feats = [MSFeature(0, 1194.4162, 1, 23.797, 2.7e5, 1),
             MSFeature(1, 1194.4162 + ISOTOPE_SPACING / 2, 1, 23.541, 2.9e5, 1)]
    out, n = repair_charges(feats)
    assert n == 1 and out[0].charge == 2


def test_lone_1plus_feature_read_as_2plus_salt():
    # A3S3 with one Na and one K replacing protons; OpenMS saw one peak and called it 1+
    m = ch.glycan_neutral_mass({"Hex": 6, "HexNAc": 5, "NeuAc": 3}, "2-AB")
    salt_mz = ch.mz(m, 2, {"Na+": 1, "K+": 1})  # true 2+ ion
    parent = [_f(0, m, 2, 36.70, 2e6), _f(1, m, 3, 36.72, 3e6)]
    lone = MSFeature(2, salt_mz, 1, 36.80, 2e5, 1)
    sp = group_species([*parent, lone])
    assert len(sp) == 1
    forms = [form for _, form in sp[0].members]
    assert any("read as 2+" in f for f in forms)
