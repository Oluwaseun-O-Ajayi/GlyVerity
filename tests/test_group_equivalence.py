"""The indexed, vectorised group_species gives the same species as the loop reference."""

import sys
from pathlib import Path

import numpy as np
import pytest

from glyverity import chemistry as ch
from glyverity.census import group_species
from glyverity.features import MSFeature

sys.path.insert(0, str(Path(__file__).parent))
from _ref.group_reference import group_species_reference

COMPS = [{"Hex": h, "HexNAc": n, "NeuAc": s} for h, n in ((5, 4), (6, 5)) for s in range(4)]


def _features(seed):
    rng = np.random.default_rng(seed)
    feats, fid = [], 0
    H, NA, K = ch.ATOMIC_MASS["H"], ch.ATOMIC_MASS["Na"], ch.ATOMIC_MASS["K"]
    for c in COMPS:
        m = ch.glycan_neutral_mass(c, "2-AB")
        rt = rng.uniform(15, 45)
        big = rng.uniform(1e5, 1e7)
        for z in (2, 3):  # charge states
            feats.append(MSFeature(fid, ch.mz(m, z), z, rt + rng.normal(0, 0.02), big / z, 3))
            fid += 1
        for nna, nk, iso in ((1, 0, 0), (2, 0, 0), (1, 1, 0), (0, 0, 1), (1, 0, 2)):  # salts
            mm = m + nna * (NA - H) + nk * (K - H) + iso * 1.0033548378
            z = int(rng.integers(2, 4))
            feats.append(MSFeature(fid, (mm + z * ch.PROTON_MASS) / z, z,
                                   rt + rng.normal(0, 0.05), big * rng.uniform(0.02, 0.3),
                                   int(rng.integers(1, 4))))
            fid += 1
        lone = m + 2 * (NA - H)  # salt seen as a lone 1+ peak
        feats.append(MSFeature(fid, (lone + 2 * ch.PROTON_MASS) / 2, 1, rt + 0.03,
                               big * 0.05, 1))
        fid += 1
    for _ in range(150):  # unrelated features
        z = int(rng.integers(1, 4))
        feats.append(MSFeature(fid, rng.uniform(600, 2000), z, rng.uniform(15, 45),
                               rng.uniform(1e3, 3e6), int(rng.integers(1, 4))))
        fid += 1
    return feats


def _summary(species):
    return sorted((round(s.neutral_mass, 6), round(s.rt, 6),
                   tuple(sorted((f.feature_id, form) for f, form in s.members)))
                  for s in species)


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("kw", [{}, {"decoy_shift": 0.41},
                                {"decoy_shift": 0.41, "decoy_families": ["adduct"]},
                                {"family_cut": {"adduct": 0.5, "reread": 0.3, "isf_sia": 0.0}}])
def test_same_species_and_records(seed, kw):
    feats = _features(seed)
    r_new, r_ref = [], []
    new = group_species(feats, record=r_new, **kw)
    ref = group_species_reference(feats, record=r_ref, **kw)
    assert _summary(new) == _summary(ref)
    assert sorted(r_new) == pytest.approx(sorted(r_ref))
