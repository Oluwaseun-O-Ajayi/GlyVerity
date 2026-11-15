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


def ms_tic(centroid_file: Any, keep: np.ndarray | None = None) -> Trace:
    """Total ion current from a centroided MS data file (sum of intensities per scan).

    Parameters
    ----------
    keep:
        Optional boolean mask, one entry per scan. Only scans where it is True
        are used (for example MS1 scans only).
    """
    times = _times(centroid_file)
    idx = np.arange(times.size) if keep is None else np.flatnonzero(keep)
    tic = np.empty(idx.size, dtype=float)
    for j, i in enumerate(idx):
        _, intensities = centroid_file.scan(int(i))
        tic[j] = float(np.sum(intensities))
    return Trace(str(centroid_file.name), "MS", times[idx], tic, dict(centroid_file.metadata or {}))


def scan_records(acqdata_path: str) -> list[dict]:
    """Per-scan records from MSScan.bin (scan time, MS level and other vendor fields).

    Uses rainbow-api's MassHunter record reader. These are internal rainbow
    functions, so the rainbow version is pinned in the optional dependency.
    """
    import os

    from rainbow.agilent import masshunter as mh

    types = mh.parse_scan_xsd(os.path.join(acqdata_path, "MSScan.xsd"))
    return list(mh.read_scan_records(
        os.path.join(acqdata_path, "MSScan.bin"), types, mh.count_scans(acqdata_path)))


def ms1_mask(records: list[dict]) -> np.ndarray | None:
    """Boolean mask of MS1 scans, or None if the records carry no ``MSLevel`` field."""
    if not records or "MSLevel" not in records[0]:
        return None
    return np.array([int(r["MSLevel"]) == 1 for r in records])


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

    import os

    directory = rainbow.read(path, centroid=True)
    files = list(directory.datafiles)
    ms = [df for df in files if str(df.detector) == "MS" and hasattr(df, "scan")]
    tic = None
    if ms:
        records = scan_records(os.path.join(path, "AcqData"))
        mask = ms1_mask(records)
        if mask is not None and mask.size != len(ms[0].xlabels):
            mask = None
        tic = ms_tic(ms[0], keep=mask)
        tic.metadata["scan_fields"] = sorted(records[0]) if records else []
        tic.metadata["ms1_filtered"] = mask is not None
        tic.metadata["n_scans_total"] = len(ms[0].xlabels)
        if records and "MSLevel" in records[0]:
            levels = [int(r["MSLevel"]) for r in records]
            tic.metadata["ms_level_counts"] = {lv: levels.count(lv) for lv in sorted(set(levels))}
    return inventory(files), fld_traces(files), tic
