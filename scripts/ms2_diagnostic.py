"""Sum MS/MS scans for given precursors and report diagnostic fragment ions.

For each precursor m/z, MS2 scans whose ``MzOfInterest`` lies within
``--tol`` of it are summed (peaks merged within 20 ppm). The script reports
the intensity of each diagnostic ion (within ``--ppm``) and the strongest
fragments.

Diagnostic ions (singly charged):
  HexNAc oxonium                    204.0866
  HexHexNAc oxonium                 366.1395
  NeuAc oxonium                     292.1027
  NeuAc - H2O                       274.0921
  Y1, free reducing GlcNAc          222.0972
  Y1, GlcNAc + 2-AB                 342.1660
  Y1, GlcNAc + tag (+43.0786 net)   265.1758   (hypothesis: C3H9N reductive amination)

Usage
-----
    python scripts/ms2_diagnostic.py RUN.d --precursors 1316.51 878.01 1355.00 903.67 --top 12
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np

from glyverity import annotate as an
from glyverity import chemistry as ch
from glyverity.io_agilent import scan_records

P = ch.PROTON_MASS
GLCNAC_FREE = ch.RESIDUE_MASS["HexNAc"] + ch.WATER
DIAGNOSTIC = {
    "HexNAc oxonium": ch.RESIDUE_MASS["HexNAc"] + P,
    "HexHexNAc oxonium": ch.RESIDUE_MASS["Hex"] + ch.RESIDUE_MASS["HexNAc"] + P,
    "NeuAc oxonium": ch.RESIDUE_MASS["NeuAc"] + P,
    "NeuAc - H2O": ch.RESIDUE_MASS["NeuAc"] - ch.WATER + P,
    "Y1 free GlcNAc": GLCNAC_FREE + P,
    "Y1 GlcNAc+2-AB": GLCNAC_FREE + ch.LABELS["2-AB"].added_mass + P,
    "Y1 GlcNAc+tag(+43.0786)": GLCNAC_FREE + 43.0786 + P,
}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("--precursors", nargs="+", type=float, required=True)
    ap.add_argument("--tol", type=float, default=0.05, help="precursor tolerance (m/z)")
    ap.add_argument("--ppm", type=float, default=20.0, help="fragment tolerance")
    ap.add_argument("--top", type=int, default=12)
    args = ap.parse_args(argv)

    import rainbow

    directory = rainbow.read(str(args.run), centroid=True)
    ms = next(df for df in directory.datafiles if str(df.detector) == "MS" and hasattr(df, "scan"))
    recs = scan_records(os.path.join(str(args.run), "AcqData"))
    levels = np.array([int(r["MSLevel"]) for r in recs])
    targets = np.array([float(r.get("MzOfInterest") or 0.0) for r in recs])
    times = np.asarray(ms.xlabels, float)

    print("diagnostic ions:", ", ".join(f"{k} {v:.4f}" for k, v in DIAGNOSTIC.items()))
    for prec in args.precursors:
        idx = np.flatnonzero((levels == 2) & (np.abs(targets - prec) <= args.tol))
        print(f"\n=== precursor {prec:.4f}: {idx.size} MS2 scans", end="")
        if idx.size == 0:
            print(" ===")
            continue
        print(f", RT {times[idx].min():.2f} to {times[idx].max():.2f} min ===")
        mz = np.concatenate([np.asarray(ms.scan(int(i))[0], float) for i in idx])
        it = np.concatenate([np.asarray(ms.scan(int(i))[1], float) for i in idx])
        mz, it = an.merge_peaks(mz, it, ppm=20.0)
        base = float(it.max()) if it.size else 1.0
        for name, theo in DIAGNOSTIC.items():
            j = int(np.argmin(np.abs(mz - theo))) if mz.size else -1
            if j >= 0 and abs(mz[j] - theo) / theo * 1e6 <= args.ppm:
                print(f"  {name:26s} {theo:9.4f}  found {mz[j]:9.4f} "
                      f"({(mz[j] - theo) / theo * 1e6:+.1f} ppm)  {100 * it[j] / base:5.1f}% of base")
            else:
                print(f"  {name:26s} {theo:9.4f}  not found")
        order = np.argsort(-it)[: args.top]
        print("  strongest fragments: " + ", ".join(f"{mz[k]:.4f} ({100 * it[k] / base:.0f}%)"
                                                     for k in order))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
