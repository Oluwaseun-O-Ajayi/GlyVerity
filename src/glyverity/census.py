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
from glyverity.annotate import ISOTOPE_SPACING, enumerate_compositions
from glyverity.candidates import Candidate, group_isomeric
from glyverity.features import MSFeature
from glyverity.layers import NEUTRAL_LOSSES, exchange_combos

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


def repair_charges(features: Sequence[MSFeature], ppm: float = 10.0,
                   rt_tol: float = 0.3) -> tuple[list[MSFeature], int]:
    """Merge co-eluting 1+ features spaced by 13C/z (z = 2, 3, 4) into one z-charged envelope.

    The feature finder sometimes splits one multiply charged envelope into several 1+
    features. The lowest-m/z feature of a chain becomes the monoisotopic peak; intensities are
    summed. Returns the repaired list and the number of merged envelopes.
    """
    ones = sorted([f for f in features if f.charge == 1], key=lambda f: f.mz)
    others = [f for f in features if f.charge != 1]
    taken: set[int] = set()
    merged: list[MSFeature] = []
    for f in ones:
        if f.feature_id in taken:
            continue
        best: list[MSFeature] = [f]
        best_z = 1
        for z in (4, 3, 2):
            chain = [f]
            while True:
                target = chain[-1].mz + ISOTOPE_SPACING / z
                nxt = [g for g in ones if g.feature_id not in taken and g not in chain
                       and abs(g.mz - target) / target * 1e6 <= ppm and abs(g.rt - f.rt) <= rt_tol]
                if not nxt:
                    break
                chain.append(max(nxt, key=lambda g: g.intensity))
            if len(chain) > len(best):
                best, best_z = chain, z
        if best_z > 1:
            taken.update(g.feature_id for g in best)
            w = np.array([g.intensity for g in best])
            rt = float(np.dot([g.rt for g in best], w) / w.sum())
            merged.append(MSFeature(f.feature_id, f.mz, best_z, rt, float(w.sum()), len(best)))
        else:
            taken.add(f.feature_id)
            merged.append(f)
    n_merged = sum(1 for m in merged if m.charge > 1)
    return others + merged, n_merged


def group_species(features: Sequence[MSFeature], ppm: float = 10.0, rt_tol: float = 0.25,
                  max_adducts: int = 4, decoy_shift: float = 0.0) -> list[Species]:
    """Group co-eluting features into species.

    Pass 1: features whose protonated neutral masses agree within ``ppm`` and whose apexes are
    within ``rt_tol`` min are charge states of one species. Pass 2: a whole species is absorbed
    as ions of a stronger co-eluting species when its (protonated) neutral mass equals the
    parent's plus any mixture of up to ``max_adducts`` Na+/K+/NH4+ exchanges (sialic acid salts
    allow more exchanges than the charge) plus 0 to 4 13C steps (monoisotopic peak picked one
    or more isotopes too high), or minus a common in-source neutral loss (with 0 to 2 13C
    steps), or twice the parent (a non-covalent [2M + ...] dimer, with the same exchanges and
    isotope errors; tried last). Absorbed members keep the relation as their ion form, so in-source fragments
    are classified separately from adduct and isotope ions.

    ``decoy_shift`` (Da) is added to every pass-2 relation offset. With a shift that matches no
    chemistry, every absorption is a chance match, so decoy runs estimate how many target
    absorptions are coincidences (see :func:`absorption_stats`).
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

    # pass 2: absorb adduct and isotope-error species into stronger co-eluting parents
    relations = []  # (rank, label, offset, multiple of the parent mass)
    for n in range(5):
        if n:
            relations.append((n, f"isotope+{n}", n * ISOTOPE_SPACING, 1))
        for counts, off in exchange_combos(max_adducts):
            label = "+".join(f"{v}{k}" for k, v in counts.items()) + (f" M+{n}" if n else "")
            relations.append((sum(counts.values()) + n, label, off + n * ISOTOPE_SPACING, 1))
    for loss, m in NEUTRAL_LOSSES.items():  # in-source fragments, with isotope error 0..2
        for n in range(3):
            relations.append((2 + n, f"isf:-{loss}" + (f" M+{n}" if n else ""),
                              -m + n * ISOTOPE_SPACING, 1))
    # non-covalent dimers [2M + ...], tried only after every single-molecule relation
    dimer_ex = [({}, 0.0), *exchange_combos(max_adducts)]
    for n in range(5):
        for counts, off in dimer_ex:
            label = "+".join(f"{v}{k}" for k, v in counts.items())
            label = "dimer:2M" + (f"+{label}" if label else "") + (f" M+{n}" if n else "")
            relations.append((100 + sum(counts.values()) + n, label, off + n * ISOTOPE_SPACING, 2))
    relations = [(rank, label, off + decoy_shift, k) for rank, label, off, k in relations]
    relations.sort(key=lambda r: r[0])
    parents = sorted(species, key=lambda s: -s.intensity)
    absorbed: set[int] = set()
    for s in sorted(species, key=lambda s: s.intensity):
        if s.species_id in absorbed:
            continue
        for parent in parents:
            if (parent.species_id == s.species_id or parent.species_id in absorbed
                    or parent.intensity <= s.intensity or abs(parent.rt - s.rt) > rt_tol):
                continue
            # a lone 1+ feature may be a multiply charged ion whose partner peaks were missed
            trial = [(s.neutral_mass, "")]
            if len(s.members) == 1 and s.members[0][0].charge == 1:
                f1 = s.members[0][0]
                trial += [(f1.mz * z - z * ch.PROTON_MASS, f" (read as {z}+)") for z in (2, 3)]
            hit = None
            for mass, note in trial:
                for _, label, off, k in relations:
                    target = k * parent.neutral_mass + off
                    if abs(mass - target) / target * 1e6 <= ppm:
                        hit = label + note
                        break
                if hit:
                    break
            if hit:
                parent.members.extend((f, hit) for f, _ in s.members)
                absorbed.add(s.species_id)
                break
    return [s for s in species if s.species_id not in absorbed]


def absorption_stats(species: Sequence[Species]) -> dict[str, float]:
    """Count features absorbed in pass 2 (ion form other than "H"), in total and for the
    lone 1+ features re-read at a higher charge."""
    out = {"n": 0, "intensity": 0.0, "n_reread": 0, "intensity_reread": 0.0}
    for s in species:
        for f, form in s.members:
            if form == "H":
                continue
            out["n"] += 1
            out["intensity"] += f.intensity
            if "(read as" in form:
                out["n_reread"] += 1
                out["intensity_reread"] += f.intensity
    return out


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


def oxonium_ions() -> list[tuple[str, float]]:
    """Singly charged non-reducing-end oxonium ions (B ions) of up to 5 residues, and their
    water-loss forms when they contain a sialic acid."""
    ranges = {"Hex": (0, 2), "HexNAc": (0, 2), "dHex": (0, 1), "NeuAc": (0, 2), "NeuGc": (0, 1)}
    out = []
    from glyverity.annotate import composition_string
    for comp in enumerate_compositions(ranges):
        if sum(comp.values()) > 5:
            continue
        m = ch.composition_mass(comp) + ch.PROTON_MASS
        name = composition_string(comp)
        out.append((name, m))
        if comp.get("NeuAc") or comp.get("NeuGc"):
            out.append((name + "-H2O", m - ch.WATER))
    return out


def reclassify_oxonium(assignments: Sequence[Assignment], classes: dict[int, str],
                       ppm: float = 10.0, rt_tol: float = 0.3) -> int:
    """Mark unassigned 1+ features that match an oxonium ion and co-elute with an assigned
    glycan species as in-source fragments. Returns the number reclassified."""
    ions = oxonium_ions()
    glycan_rts = [a.species.rt for a in assignments
                  if any(classes.get(f.feature_id) == "labelled glycan" for f, _ in a.species.members)
                  or any(classes.get(f.feature_id) == "side-product series" for f, _ in a.species.members)]
    n = 0
    for a in assignments:
        for f, _ in a.species.members:
            if classes.get(f.feature_id) != "unassigned" or f.charge != 1:
                continue
            if not any(abs(f.rt - r) <= rt_tol for r in glycan_rts):
                continue
            if any(abs(f.mz - m) / m * 1e6 <= ppm for _, m in ions):
                classes[f.feature_id] = "in-source fragment"
                n += 1
    return n


CLASSES = ("labelled glycan", "other end group", "side-product series", "in-source fragment",
           "adduct/isotope", "adduct ion of a species", "dimer ion", "ambiguous", "unassigned")
ADDUCT_ION = "adduct ion of a species"  # adduct, salt or isotope-error ion of a stronger species


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
            if form == "H":
                out[f.feature_id] = cls
            elif form.startswith("isf:"):
                out[f.feature_id] = "in-source fragment"
            elif form.startswith("dimer:"):
                out[f.feature_id] = "dimer ion"
            else:
                out[f.feature_id] = ADDUCT_ION
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
