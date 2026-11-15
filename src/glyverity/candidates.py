"""Targeted candidates: composition and reducing-end group hypotheses with extracted MS1 evidence."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

from glyverity import annotate as an
from glyverity import chemistry as ch
from glyverity.xic import summarise, xic


@dataclass
class Candidate:
    """One composition with one reducing-end group, and its extracted MS1 signal."""

    composition: dict[str, int]
    end_group: str
    neutral_mass: float
    apex_time: float = float("nan")
    apex_intensity: float = 0.0
    area: float = 0.0
    fld_at_apex: float = float("nan")
    flags: list[str] = field(default_factory=list)

    @property
    def name(self) -> str:
        return an.composition_string(self.composition)

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
    """Build one candidate per composition and end group, with XIC apex, area and FLD at apex.

    ``end_groups`` maps a name to the mass it adds to the free glycan.
    ``fld`` is (times, values); its value is read at ``apex - lag``.
    """
    lo, hi = window
    specs, targets = [], []
    for comp in compositions:
        for group, added in end_groups.items():
            neutral = ch.composition_mass(comp) + ch.WATER + added
            tg = target_mzs(neutral, charges)
            specs.append((comp, group, neutral, len(targets), len(targets) + len(tg)))
            targets += tg
    times, traces = xic(centroid_file, keep, targets, ppm, lo, hi)  # one pass over the scans
    out = []
    for comp, group, neutral, a, b in specs:
        s = summarise(times, traces[:, a:b].sum(axis=1))
        c = Candidate(dict(comp), group, neutral, s.apex_time, s.apex_intensity, s.area)
        if fld is not None and np.isfinite(s.apex_time):
            c.fld_at_apex = float(np.interp(s.apex_time - lag, fld[0], fld[1]))
        out.append(c)
    return out
