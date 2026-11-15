"""Feature-first census: group features into species, assign compositions, classify intensity.

Each observed MS1 feature is counted once. Features are grouped into species (one neutral
molecule seen at several charges and as adduct ions), species are matched to compositions with
each reducing-end group, the artifact layers run on species, and every feature's intensity is
assigned to exactly one class.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

from glyverity import chemistry as ch
from glyverity.candidates import Candidate, group_isomeric
from glyverity.features import MSFeature

ADDUCT_FORMS = {  # name: (number of non-proton carriers, carrier)
    "Na": "Na+",
    "K": "K+",
    "NH4": "NH4+",
}


def neutral_for(feature: MSFeature, carrier: str | None = None, k: int = 0) -> float:
    """Neutral mass of a feature if it carries ``k`` ``carrier`` ions and protons for the rest."""
    z = feature.charge
    added = (z - k) * ch.PROTON_MASS + (k * ch.CHARGE_CARRIERS[carrier] if carrier else 0.0)
    return feature.mz * z - added


@dataclass
class Species:
    """One neutral molecule and the features that represent it."""

    species_id: int
    neutral_mass: float
    rt: float
    members: list[tuple[MSFeature, str]] = field(default_factory=list)  # (feature, ion form)

    @property
    def intensity(self) -> float:
        return float(sum(f.intensity for f, _ in self.members))

    @property
    def protonated_intensity(self) -> float:
        return float(sum(f.intensity for f, form in self.members if form == "H"))


def group_species(features: Sequence[MSFeature], ppm: float = 10.0, rt_tol: float = 0.2,
                  max_adducts: int = 2) -> list[Species]:
    """Group co-eluting features into species.

    Pass 1: features whose protonated neutral masses agree within ``ppm`` and whose apexes are
    within ``rt_tol`` min are charge states of one species. Pass 2: a whole species is absorbed
    as adduct ions of a stronger co-eluting species when its (protonated) neutral mass equals
    the parent's plus k x (adduct - H), k = 1..``max_adducts``, for Na+, K+ or NH4+. This keeps
    the 2+ and 3+ sodium adduct ions together with their parent.
    """
    feats = sorted(features, key=lambda f: -f.intensity)
    species: list[Species] = []
    used: set[int] = set()
    for f in feats:  # pass 1: charge states
        if f.feature_id in used:
            continue
        m = neutral_for(f)
        s = Species(len(species), m, f.rt, [(f, "H")])
        used.add(f.feature_id)
        for g in feats:
            if g.feature_id in used or g.charge == f.charge or abs(g.rt - f.rt) > rt_tol:
                continue
            if abs(neutral_for(g) - m) / m * 1e6 <= ppm:
                s.members.append((g, "H"))
                used.add(g.feature_id)
        _refresh(s)
        species.append(s)

    # pass 2: absorb adduct species into stronger co-eluting parents
    exchange = {name: ch.CHARGE_CARRIERS[carrier] - ch.PROTON_MASS
                for name, carrier in ADDUCT_FORMS.items()}
    parents = sorted(species, key=lambda s: -s.intensity)
    absorbed: set[int] = set()
    for s in sorted(species, key=lambda s: s.intensity):
        if s.species_id in absorbed:
            continue
        for parent in parents:
            if (parent.species_id == s.species_id or parent.species_id in absorbed
                    or parent.intensity <= s.intensity or abs(parent.rt - s.rt) > rt_tol):
                continue
            hit = None
            for name, d in exchange.items():
                for k in range(1, max_adducts + 1):
                    if abs(s.neutral_mass - (parent.neutral_mass + k * d)) / parent.neutral_mass \
                            * 1e6 <= ppm:
                        hit = f"{k}{name}"
                        break
                if hit:
                    break
            if hit:
                parent.members.extend((f, hit) for f, _ in s.members)
                absorbed.add(s.species_id)
                break
    return [s for s in species if s.species_id not in absorbed]


def _refresh(s: Species) -> None:
    """Recompute a species' neutral mass and time from its protonated members."""
    prot = [f for f, form in s.members if form == "H"]
    w = np.array([f.intensity for f in prot])
    s.neutral_mass = float(np.dot([neutral_for(f) for f in prot], w) / w.sum())
    s.rt = float(np.dot([f.rt for f in prot], w) / w.sum())


@dataclass
class Assignment:
    """All composition and end-group matches of one species."""

    species: Species
    matches: list[Candidate]


def assign_species(species: Sequence[Species], compositions: Sequence[Mapping[str, int]],
                   end_groups: Mapping[str, float], ppm: float = 10.0) -> list[Assignment]:
    """Match each species' neutral mass to every composition x end group within ``ppm``.

    Isomeric compositions are merged first, so one formula gives one match per end group.
    """
    table = []
    for iso in group_isomeric(compositions):
        base = ch.composition_mass(iso[0]) + ch.WATER
        for group, added in end_groups.items():
            table.append((base + added, iso, group))
    masses = np.array([t[0] for t in table])
    order = np.argsort(masses)
    masses = masses[order]
    table = [table[i] for i in order]
    out = []
    for s in species:
        lo = np.searchsorted(masses, s.neutral_mass * (1 - ppm * 1e-6))
        hi = np.searchsorted(masses, s.neutral_mass * (1 + ppm * 1e-6), side="right")
        matches = []
        for j in range(lo, hi):
            mass, iso, group = table[j]
            c = Candidate(dict(iso[0]), group, mass, s.rt, s.protonated_intensity, s.intensity,
                          alternatives=[dict(x) for x in iso])
            c.uid = f"{c.name}|{group}|s{s.species_id}"
            matches.append(c)
        out.append(Assignment(s, matches))
    return out


CLASSES = ("labelled glycan", "other end group", "side-product series", "in-source fragment",
           "adduct/isotope", "adduct ion of a species", "ambiguous", "unassigned")


def classify(assignments: Sequence[Assignment], label: str) -> dict[int, str]:
    """Class of each feature (by feature_id) after the layers have set candidate flags.

    A species with one match takes that match's class; with several non-flagged matches it is
    ambiguous; with none it is unassigned. Adduct-ion members of a species are counted as
    'adduct ion of a species'.
    """
    out: dict[int, str] = {}
    for a in assignments:
        live = [c for c in a.matches if not c.flags]
        if not a.matches:
            cls = "unassigned"
        elif not live:
            f0 = a.matches[0].flags[0]
            cls = ("in-source fragment" if f0.startswith("isf") else
                   "side-product series" if f0.startswith("side_product") else "adduct/isotope")
        elif len(live) > 1:
            cls = "ambiguous"
        else:
            cls = "labelled glycan" if live[0].end_group == label else "other end group"
        for f, form in a.species.members:
            out[f.feature_id] = cls if form == "H" else "adduct ion of a species"
    return out


def unassigned_report(assignments: Sequence[Assignment], classes: Mapping[int, str],
                      top: int = 20, rt_tol: float = 0.2) -> list[tuple[MSFeature, float, str]]:
    """Strongest unassigned features, each with the neutral mass difference to the strongest
    co-eluting assigned species (feature neutral, protonated, minus species neutral)."""
    assigned = [a.species for a in assignments
                if a.matches and any(classes.get(f.feature_id) not in (None, "unassigned")
                                     for f, _ in a.species.members)]
    label = {a.species.species_id: (a.matches[0].name + "(" + a.matches[0].end_group + ")")
             for a in assignments if a.matches}
    rows = []
    feats = [f for a in assignments for f, _ in a.species.members
             if classes.get(f.feature_id) == "unassigned"]
    for f in sorted(feats, key=lambda f: -f.intensity)[:top]:
        near = [s for s in assigned if abs(s.rt - f.rt) <= rt_tol]
        if near:
            s = max(near, key=lambda s: s.intensity)
            rows.append((f, neutral_for(f) - s.neutral_mass, label[s.species_id]))
        else:
            rows.append((f, float("nan"), ""))
    return rows
