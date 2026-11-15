"""Evidence-layer plugin interface.

Every source of evidence about a feature (isotope pattern, adduct grouping,
label side-product series, in-source fragment, retention, ion mobility CCS,
biosynthetic feasibility) is an evidence layer. Layers are registered by name
and combined by the engine, so new layers can be added without changing the
core.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class Feature:
    """One LC-MS feature: a monoisotopic m/z, charge, retention time and intensity."""

    feature_id: str
    mz: float
    charge: int
    rt: float
    intensity: float
    extra: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Evidence:
    """A statement from one layer about one feature."""

    feature_id: str
    layer: str
    kind: str
    score: float
    detail: Mapping[str, Any] = field(default_factory=dict)


class EvidenceLayer(ABC):
    """Base class for evidence layers.

    Subclasses set :attr:`name` and :attr:`version` and implement
    :meth:`evaluate`. A layer must be deterministic for a given input and
    configuration so that results are reproducible.
    """

    name: str = ""
    version: str = "0"

    @abstractmethod
    def evaluate(self, features: Sequence[Any]) -> list[Evidence]:
        """Return evidence for the given features or candidates."""


_REGISTRY: dict[str, type[EvidenceLayer]] = {}


def register(cls: type[EvidenceLayer]) -> type[EvidenceLayer]:
    """Class decorator that registers an evidence layer under ``cls.name``."""
    if not cls.name:
        raise ValueError("evidence layer must define a non-empty name")
    if cls.name in _REGISTRY and _REGISTRY[cls.name] is not cls:
        raise ValueError(f"evidence layer '{cls.name}' is already registered")
    _REGISTRY[cls.name] = cls
    return cls


def available_layers() -> list[str]:
    """Names of all registered evidence layers, sorted."""
    return sorted(_REGISTRY)


def get_layer(name: str) -> type[EvidenceLayer]:
    """Return the registered layer class for ``name``."""
    try:
        return _REGISTRY[name]
    except KeyError as exc:
        raise KeyError(f"no evidence layer named '{name}'") from exc


def _unregister(name: str) -> None:
    """Remove a layer (used by tests)."""
    _REGISTRY.pop(name, None)


LayerFactory = Callable[[], EvidenceLayer]
