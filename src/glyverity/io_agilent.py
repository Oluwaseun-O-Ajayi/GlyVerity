"""Read fluorescence (FLD) traces and centroided MS data from Agilent MassHunter .d folders.

Parsing of the vendor binaries is delegated to ``rainbow-api`` (LGPL-3.0), an
optional dependency installed with ``pip install glyverity[agilent]``. This
module only selects the traces GlyVerity needs and reduces them to plain
NumPy arrays.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass(frozen=True)
class Trace:
    """A one-dimensional chromatographic trace."""

    name: str
    detector: str
    times: np.ndarray
    values: np.ndarray
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class InventoryItem:
    """One data file found in a .d folder, as reported by the reader."""

    name: str
    detector: str
    n_points: int
    t_start: float
    t_end: float


def _times(datafile: Any) -> np.ndarray:
    return np.asarray(datafile.xlabels, dtype=float)


def inventory(datafiles: Iterable[Any]) -> list[InventoryItem]:
    """Summarise every data file the reader returned."""
    items = []
    for df in datafiles:
        t = _times(df)
        items.append(
            InventoryItem(
                name=str(df.name),
                detector=str(df.detector),
                n_points=int(t.size),
                t_start=float(t[0]) if t.size else float("nan"),
                t_end=float(t[-1]) if t.size else float("nan"),
            )
        )
    return items


def fld_traces(datafiles: Iterable[Any]) -> list[Trace]:
    """Return every fluorescence trace (detector ``FLD``)."""
    out = []
    for df in datafiles:
        if str(df.detector) != "FLD":
            continue
        values = np.asarray(df.data, dtype=float).reshape(-1)
        out.append(Trace(str(df.name), "FLD", _times(df), values, dict(df.metadata or {})))
    return out


def ms_tic(centroid_file: Any) -> Trace:
    """Total ion current from a centroided MS data file (sum of intensities per scan)."""
    times = _times(centroid_file)
    tic = np.empty(times.size, dtype=float)
    for i in range(times.size):
        _, intensities = centroid_file.scan(i)
        tic[i] = float(np.sum(intensities))
    return Trace(str(centroid_file.name), "MS", times, tic, dict(centroid_file.metadata or {}))


def read_d(path: str) -> tuple[list[InventoryItem], list[Trace], Trace | None]:
    """Read an Agilent .d folder.

    Returns
    -------
    inventory, fld, tic
        Every data file found, the fluorescence traces, and the MS TIC from
        ``MSPeak.bin`` (``None`` if no centroided MS file is present).
    """
    try:
        import rainbow
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError("install the Agilent reader with: pip install 'glyverity[agilent]'") from exc

    directory = rainbow.read(path, centroid=True)
    files = list(directory.datafiles)
    ms = [df for df in files if str(df.detector) == "MS" and hasattr(df, "scan")]
    tic = ms_tic(ms[0]) if ms else None
    return inventory(files), fld_traces(files), tic
