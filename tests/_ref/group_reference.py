"""Frozen loop implementation of census.group_species (2026-11-23), for equivalence tests."""
# ruff: noqa
from collections.abc import Mapping, Sequence

from glyverity import chemistry as ch
from glyverity.census import (FAMILIES, ISOTOPE_SPACING, NEUTRAL_LOSSES, Species, _refresh,
                              exchange_combos, neutral_for, relation_family)
from glyverity.features import MSFeature


def group_species_reference(features: Sequence[MSFeature], ppm: float = 10.0, rt_tol: float = 0.25,
                  max_adducts: int = 4, decoy_shift: float = 0.0,
                  decoy_families: Sequence[str] | None = None,
                  family_cut: Mapping[str, float] | None = None,
                  record: list[tuple[str, float]] | None = None) -> list[Species]:
    """Group co-eluting features into species.

    Pass 1: features whose protonated neutral masses agree within ``ppm`` and whose apexes are
    within ``rt_tol`` min are charge states of one species. Pass 2: a whole species is absorbed
    as ions of a stronger co-eluting species when its (protonated) neutral mass equals the
    parent's plus any mixture of up to ``max_adducts`` Na+/K+/NH4+ exchanges (sialic acid salts
    allow more exchanges than the charge) plus 0 to 4 13C steps (monoisotopic peak picked one
    or more isotopes too high), or minus a common in-source neutral loss (with 0 to 2 13C
    steps), or twice the parent (a non-covalent [2M + ...] dimer, with the same exchanges and
    isotope errors; tried last, and only when the features of both the candidate and the parent
    have a measured charge, i.e. two or more isotopes, because [2M+2zH]2z+ has the m/z of
    [M+zH]z+). Ions already absorbed into a species keep their own relation when that species
    is absorbed in turn (form "<own> via <species relation>"). Absorbed members keep the relation as their ion form, so in-source fragments
    are classified separately from adduct and isotope ions.

    ``decoy_shift`` (Da) is added to the pass-2 relation offsets of the families in
    ``decoy_families`` (all families when None). Shifting one family while the others stay
    real keeps the decoy's pool of candidate features the same as the target's; shifting all
    at once leaves features free that real adducts would claim, which inflates the decoy count
    of families tested later. With a shift that matches no
    chemistry, every absorption is a chance match, so decoy runs estimate how many target
    absorptions are coincidences (see :func:`absorption_stats`).

    Each candidate relation gets a score max(ppm error / ``ppm``, |retention difference| /
    ``rt_tol``) in [0, 1]: adducts, isotope errors and in-source fragments form in the source,
    so they co-elute exactly with their parent, while chance matches spread over both the mass
    and the time window. ``family_cut`` sets a score cutoff per relation family (see
    :func:`relation_family` and :func:`family_cutoffs`); ``record`` collects (family, score)
    for every absorption, for target and decoy runs alike.
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
    shifted = set(FAMILIES if decoy_families is None else decoy_families)
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
            # a lone 1+ feature may be a multiply charged ion whose partner peaks were missed;
            # "lone" counts pass-1 members only, so ions already absorbed into s do not block it
            own = [f for f, form in s.members if form == "H"]
            trial = [(s.neutral_mass, "")]
            if len(own) == 1 and own[0].charge == 1:
                f1 = own[0]
                trial += [(f1.mz * z - z * ch.PROTON_MASS, f" (read as {z}+)") for z in (2, 3)]
            hit = None
            measured = all(f.n_isotopes >= 2 for f in own) and all(
                f.n_isotopes >= 2 for f, form in parent.members if form == "H")
            for mass, note in trial:
                for _, label, off, k in relations:
                    # a dimer needs a measured charge: [2M+2zH]2z+ sits at the m/z of [M+zH]z+
                    if k == 2 and (note or not measured):
                        continue
                    fam = relation_family(label + note)
                    target = k * parent.neutral_mass + off + (decoy_shift if fam in shifted else 0.0)
                    err = abs(mass - target) / target * 1e6
                    if err > ppm:
                        continue
                    score = max(err / ppm, abs(parent.rt - s.rt) / rt_tol)
                    limit = 1.0 if family_cut is None else family_cut.get(fam, 1.0)
                    if score <= limit:
                        hit = label + note
                        if record is not None:
                            record.append((fam, score))
                        break
                if hit:
                    break
            if hit:
                # members s had absorbed keep their own relation, chained to s's relation
                parent.members.extend((f, hit if form == "H" else f"{form} via {hit}")
                                      for f, form in s.members)
                absorbed.add(s.species_id)
                break
    return [s for s in species if s.species_id not in absorbed]

