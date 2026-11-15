"""Compare MS1 signal of the same glycan compositions under different reducing-end groups.

For each composition and reducing-end group, XICs are summed over charge states
(2+ and 3+, protonated) and the monoisotopic plus first 13C isotope peaks.
Each series is reported with apex time, apex intensity, area, and the
fluorescence signal at the apex (shifted by the measured FLD-to-MS lag).

Reducing-end groups:
  2-AB         2-aminobenzamide label (+120.0687)
  free         unlabelled reducing end
  C3H9N-RA     net +C3H9N-O (+43.0786), hypothesis from the fetuin MS2 data

MS signal is not molar abundance: ionisation efficiency differs between groups.

Usage
-----
    python scripts/label_series.py RUN.d OUT_DIR --window 15 50 --lag 0.037
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np

from glyverity import annotate as an
from glyverity import chemistry as ch
from glyverity.io_agilent import fld_traces, ms1_mask, scan_records
from glyverity.xic import summarise, xic

TAG = ch.ReducingEndLabel("C3H9N-RA", ch.formula_mass({"C": 3, "H": 9, "N": 1}) - ch.ATOMIC_MASS["O"])
GROUPS = {"2-AB": ch.LABELS["2-AB"].added_mass, "free": 0.0, "C3H9N-RA": TAG.added_mass}

# Complex N-glycans expected on bovine fetuin and close relatives.
COMPOSITIONS = [
    {"Hex": h, "HexNAc": n, "NeuAc": s, **({"dHex": f} if f else {})}
    for (h, n) in [(5, 4), (6, 5)]
    for s in range(5 if n == 5 else 3)
    for f in (0, 1)
]


def target_mzs(neutral: float) -> list[float]:
    out = []
    for z in (2, 3):
        mono = ch.mz(neutral, z)
        out += [mono, mono + an.ISOTOPE_SPACING / z]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--window", nargs=2, type=float, default=(15.0, 50.0))
    ap.add_argument("--lag", type=float, default=0.0, help="MS minus FLD time (min)")
    ap.add_argument("--ppm", type=float, default=10.0)
    args = ap.parse_args(argv)

    import rainbow

    args.out.mkdir(parents=True, exist_ok=True)
    directory = rainbow.read(str(args.run), centroid=True)
    files = list(directory.datafiles)
    ms = next(df for df in files if str(df.detector) == "MS" and hasattr(df, "scan"))
    fld = fld_traces(files)
    keep = ms1_mask(scan_records(os.path.join(str(args.run), "AcqData")))
    lo, hi = args.window

    rows, totals = [], {g: 0.0 for g in GROUPS}
    for comp in COMPOSITIONS:
        name = an.composition_string(comp)
        for group, added in GROUPS.items():
            neutral = ch.composition_mass(comp) + ch.WATER + added
            times, traces = xic(ms, keep, target_mzs(neutral), args.ppm, lo, hi)
            s = summarise(times, traces.sum(axis=1))
            fld_at = float("nan")
            if fld and np.isfinite(s.apex_time):
                fld_at = float(np.interp(s.apex_time - args.lag, fld[0].times, fld[0].values))
            rows.append([name, group, f"{neutral:.4f}", f"{s.apex_time:.3f}",
                         f"{s.apex_intensity:.0f}", f"{s.area:.0f}", f"{fld_at:.3f}"])
            totals[group] += s.area

    path = args.out / "label_series.tsv"
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["composition", "reducing_end", "neutral_mass", "apex_min", "apex_intensity",
                    "area", "fld_at_apex"])
        w.writerows(rows)

    print(f"{'composition':22s} {'end':9s} {'apex':>7s} {'apex int':>10s} {'area':>12s} {'FLD':>7s}")
    for r in rows:
        if float(r[5]) > 0:
            print(f"{r[0]:22s} {r[1]:9s} {r[3]:>7s} {r[4]:>10s} {r[5]:>12s} {r[6]:>7s}")
    grand = sum(totals.values())
    print("\nMS1 area by reducing-end group (sum over listed compositions):")
    for g, v in totals.items():
        share = 100 * v / grand if grand else 0.0
        print(f"  {g:9s} {v:14.0f}  {share:5.1f}%")
    print(f"table: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
