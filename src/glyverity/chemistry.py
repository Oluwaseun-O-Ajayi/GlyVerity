"""Mass calculations for released glycans, reducing-end labels and adducts.

All masses are monoisotopic and computed from elemental formulas, so every
value traces back to the atomic masses defined in this module.

Conventions
-----------
* A glycan composition is a mapping of residue name to count, for example
  ``{"Hex": 5, "HexNAc": 2}`` for Man5.
* Residue masses are dehydrated (glycosidic) residues. The neutral mass of a
  free released glycan is the sum of residue masses plus one water.
* A reductive-amination label adds ``label - H2O + H2`` to the free glycan.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

# Monoisotopic atomic masses (u).
ATOMIC_MASS: Mapping[str, float] = MappingProxyType(
    {
        "C": 12.0,
        "H": 1.00782503207,
        "N": 14.0030740048,
        "O": 15.99491461956,
        "Na": 22.9897692809,
        "K": 38.96370668,
    }
)

ELECTRON_MASS = 0.00054857990946
PROTON_MASS = ATOMIC_MASS["H"] - ELECTRON_MASS


def formula_mass(formula: Mapping[str, int]) -> float:
    """Return the monoisotopic mass of an elemental formula.

    Parameters
    ----------
    formula:
        Mapping of element symbol to atom count, for example ``{"C": 6, "H": 12, "O": 6}``.

    Raises
    ------
    KeyError
        If an element is not defined in :data:`ATOMIC_MASS`.
    """
    return sum(ATOMIC_MASS[element] * count for element, count in formula.items())


WATER = formula_mass({"H": 2, "O": 1})
H2 = formula_mass({"H": 2})

# Dehydrated residue formulas.
RESIDUE_FORMULA: Mapping[str, Mapping[str, int]] = MappingProxyType(
    {
        "Hex": {"C": 6, "H": 10, "O": 5},
        "HexNAc": {"C": 8, "H": 13, "N": 1, "O": 5},
        "dHex": {"C": 6, "H": 10, "O": 4},
        "NeuAc": {"C": 11, "H": 17, "N": 1, "O": 8},
        "NeuGc": {"C": 11, "H": 17, "N": 1, "O": 9},
    }
)

RESIDUE_MASS: Mapping[str, float] = MappingProxyType(
    {name: formula_mass(f) for name, f in RESIDUE_FORMULA.items()}
)


@dataclass(frozen=True)
class ReducingEndLabel:
    """A reducing-end modification and the mass it adds to a free glycan."""

    name: str
    added_mass: float
    reagent_formula: Mapping[str, int] | None = None


def reductive_amination_label(name: str, reagent_formula: Mapping[str, int]) -> ReducingEndLabel:
    """Build a label added by reductive amination (Schiff base, then reduction).

    The added mass is ``reagent - H2O + H2``.
    """
    added = formula_mass(reagent_formula) - WATER + H2
    return ReducingEndLabel(name=name, added_mass=added, reagent_formula=dict(reagent_formula))


LABELS: Mapping[str, ReducingEndLabel] = MappingProxyType(
    {
        "free": ReducingEndLabel("free", 0.0),
        "reduced": ReducingEndLabel("reduced", H2),
        # 2-aminobenzamide, C7H8N2O
        "2-AB": reductive_amination_label("2-AB", {"C": 7, "H": 8, "N": 2, "O": 1}),
        # procainamide, 4-amino-N-[2-(diethylamino)ethyl]benzamide, C13H21N3O
        "ProA": reductive_amination_label("ProA", {"C": 13, "H": 21, "N": 3, "O": 1}),
    }
)

# Neutral mass shifts observed as artifact series (added to the labelled glycan).
ARTIFACT_SHIFTS: Mapping[str, float] = MappingProxyType(
    {
        "acetonitrile": formula_mass({"C": 2, "H": 3, "N": 1}),
        "water_loss": -WATER,
    }
)

# Charge carriers: mass added per charge, relative to the neutral molecule.
CHARGE_CARRIERS: Mapping[str, float] = MappingProxyType(
    {
        "H+": PROTON_MASS,
        "Na+": ATOMIC_MASS["Na"] - ELECTRON_MASS,
        "K+": ATOMIC_MASS["K"] - ELECTRON_MASS,
        "NH4+": formula_mass({"N": 1, "H": 4}) - ELECTRON_MASS,
    }
)


def composition_mass(composition: Mapping[str, int]) -> float:
    """Return the summed residue mass of a composition (no reducing-end water).

    Raises
    ------
    ValueError
        If a count is negative.
    KeyError
        If a residue name is unknown.
    """
    total = 0.0
    for residue, count in composition.items():
        if count < 0:
            raise ValueError(f"negative count for {residue}: {count}")
        total += RESIDUE_MASS[residue] * count
    return total


def glycan_neutral_mass(composition: Mapping[str, int], label: str = "free") -> float:
    """Neutral monoisotopic mass of a released glycan with a reducing-end label."""
    return composition_mass(composition) + WATER + LABELS[label].added_mass


def mz(neutral_mass: float, charge: int, carriers: Mapping[str, int] | None = None) -> float:
    """Return m/z for a positive ion.

    Parameters
    ----------
    neutral_mass:
        Neutral monoisotopic mass of the molecule.
    charge:
        Positive charge state.
    carriers:
        Mapping of carrier name to count, summing to ``charge``. Defaults to
        all protons.
    """
    if charge <= 0:
        raise ValueError("charge must be a positive integer")
    if carriers is None:
        carriers = {"H+": charge}
    if sum(carriers.values()) != charge:
        raise ValueError("carrier counts must sum to the charge state")
    added = sum(CHARGE_CARRIERS[name] * n for name, n in carriers.items())
    return (neutral_mass + added) / charge


def label_net_formula(label: str) -> dict[str, int]:
    """Net elemental change a reducing-end label makes to the free glycan.

    Reductive amination adds ``reagent - O`` (reagent - H2O + H2). ``free`` adds nothing and
    ``reduced`` adds H2.
    """
    if label == "free":
        return {}
    if label == "reduced":
        return {"H": 2}
    reagent = LABELS[label].reagent_formula
    if reagent is None:
        raise ValueError(f"no reagent formula recorded for label {label}")
    net = dict(reagent)
    net["O"] = net.get("O", 0) - 1
    return {k: v for k, v in net.items() if v}


def glycan_formula(composition: Mapping[str, int], extra: Mapping[str, int] | None = None) -> dict[str, int]:
    """Elemental formula of a free released glycan (residues + H2O) plus ``extra`` atoms."""
    total: dict[str, int] = {"H": 2, "O": 1}
    for residue, n in composition.items():
        for el, k in RESIDUE_FORMULA[residue].items():
            total[el] = total.get(el, 0) + k * n
    for el, k in (extra or {}).items():
        total[el] = total.get(el, 0) + k
    return {k: v for k, v in total.items() if v}


def formula_string(formula: Mapping[str, int]) -> str:
    """Hill-order formula string accepted by OpenMS (negative counts allowed, e.g. ``O-1``)."""
    order = ["C", "H"] + sorted(k for k in formula if k not in ("C", "H"))
    return "".join(f"{el}{formula[el]}" for el in order if formula.get(el))
