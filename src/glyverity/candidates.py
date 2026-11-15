"""Targeted candidates: composition and reducing-end group hypotheses with extracted MS1 evidence."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

from glyverity import annotate as an
from glyverity import chemistry as ch
from glyverity.xic import summarise, xic


def composition_formula(comp: Mapping[str, int]) -> tuple[tuple[str, int], ...]:
    """Elemental formula of the residues in a composition, as sorted (element, count) pairs."""
    total: dict[str, int] = {}
    for residue, n in comp.items():
        for el, k in ch.RESIDUE_FORMULA[residue].items():
            total[el] = total.get(el, 0) + k * n
    return tuple(sorted(total.items()))


def _display_order(comp: Mapping[str, int]) -> tuple[int, int, int]:
    """Display preference among isomeric compositions: fewer NeuGc, then fewer dHex.

    This only chooses which name to print first. MS1 cannot distinguish them, and
    all alternatives are kept.
    """
    return (comp.get("NeuGc", 0), comp.get("dHex", 0), -comp.get("NeuAc", 0))


def group_isomeric(compositions: Sequence[Mapping[str, int]]) -> list[list[dict[str, int]]]:
    """Group compositions with identical elemental formulas; each group sorted for display."""
    groups: dict[tuple, list[dict[str, int]]] = {}
    for c in compositions:
        groups.setdefault(composition_formula(c), []).append(dict(c))
    return [sorted(g, key=_display_order) for g in groups.values()]


@dataclass
class Candidate:
    """One elemental composition with one reducing-end group, and its extracted MS1 signal.

    ``composition`` is the display choice; ``alternatives`` lists every composition with the
    same elemental formula (including ``composition``), which MS1 alone cannot separate.
    """

    composition: dict[str, int]
    end_group: str
    neutral_mass: float
    apex_time: float = float("nan")
    apex_intensity: float = 0.0
    area: float = 0.0
    fld_at_apex: float = float("nan")
    flags: list[str] = field(default_factory=list)
    alternatives: list[dict[str, int]] = field(default_factory=list)

    @property
    def name(self) -> str:
        return an.composition_string(self.composition)

    @property
    def alternative_names(self) -> list[str]:
        return [an.composition_string(c) for c in self.alternatives if c != self.composition]

    @property
    def key(self) -> str:
        return f"{self.name}|{self.end_group}"


def target_mzs(neutral: float, charges: Sequence[int] = (2, 3), isotopes: int = 2) -> list[float]:
    """Protonated m/z for each charge, monoisotopic plus ``isotopes - 1`` 13C peaks."""
    out = []
    for z in charges:
        mono = ch.mz(neutral, z)
        out += [mono + k * an.ISOTOPE_SPACING / z for k in range(isotopes)]
    return out


def extract_candidates(
    centroid_file,
    keep: np.ndarray | None,
    compositions: Sequence[Mapping[str, int]],
    end_groups: Mapping[str, float],
    window: tuple[float, float],
    ppm: float = 10.0,
    fld: tuple[np.ndarray, np.ndarray] | None = None,
    lag: float = 0.0,
    charges: Sequence[int] = (2, 3),
) -> list[Candidate]:
    """Build one candidate per elemental composition and end group (isomeric compositions are
    merged), with XIC apex, area and FLD at apex.

    ``end_groups`` maps a name to the mass it adds to the free glycan.
    ``fld`` is (times, values); its value is read at ``apex - lag``.
    """
    lo, hi = window
    specs, targets = [], []
    for iso in group_isomeric(compositions):
        comp = iso[0]
        for group, added in end_groups.items():
            neutral = ch.composition_mass(comp) + ch.WATER + added
            tg = target_mzs(neutral, charges)
            specs.append((comp, iso, group, neutral, len(targets), len(targets) + len(tg)))
            targets += tg
    times, traces = xic(centroid_file, keep, targets, ppm, lo, hi)  # one pass over the scans
    out = []
    for comp, iso, group, neutral, a, b in specs:
        s = summarise(times, traces[:, a:b].sum(axis=1))
        c = Candidate(dict(comp), group, neutral, s.apex_time, s.apex_intensity, s.area,
                      alternatives=[dict(x) for x in iso])
        if fld is not None and np.isfinite(s.apex_time):
            c.fld_at_apex = float(np.interp(s.apex_time - lag, fld[0], fld[1]))
        out.append(c)
    return out
