"""Released glycans or glycopeptides? Diagnostic-ion counts in the MS2 spectra of an mzML file.

Every glycan-containing MS2 spectrum (released glycan or glycopeptide) tends to show the HexNAc
oxonium ion (m/z 204.0867). Only a released, 2-AB-labelled N-glycan gives the Y1 ion
[GlcNAc + 2-AB + H]+ at m/z 342.1660, because its reducing-end GlcNAc carries the label; a
glycopeptide's reducing-end GlcNAc is attached to asparagine instead. The share of
oxonium-positive spectra that also show 342.1660 therefore separates the two cases, provided
the low-mass region was recorded. Precursor charges are reported too.

Usage: python scripts/ms2_mzml_check.py run.mzML [--ppm 20] [--da 0.3]
"""

from __future__ import annotations

import argparse
from collections import Counter

import numpy as np

IONS = {
    "HexNAc oxonium 204.0867": 204.0867,
    "HexHexNAc oxonium 366.1395": 366.1395,
    "NeuAc oxonium 292.1027": 292.1027,
    "Y1 GlcNAc+2-AB 342.1660": 342.1660,
}


def has(mz: np.ndarray, it: np.ndarray, target: float, tol_da: float, min_rel: float) -> bool:
    k = np.abs(mz - target) <= tol_da
    return bool(k.any() and it[k].max() >= min_rel * it.max())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mzml")
    ap.add_argument("--ppm", type=float, default=20.0, help="tolerance for high-resolution MS2")
    ap.add_argument("--da", type=float, default=0.3, help="tolerance for ion-trap MS2 (Da)")
    ap.add_argument("--min-rel", type=float, default=0.01,
                    help="ion must reach this fraction of the spectrum's base peak")
    args = ap.parse_args()

    import pyopenms as oms

    exp = oms.MSExperiment()
    oms.MzMLFile().load(args.mzml, exp)
    ms2 = [s for s in exp.getSpectra() if s.getMSLevel() == 2 and s.size()]
    low = [s for s in ms2 if s.get_peaks()[0].min() < 200.0]
    print(f"MS2 spectra: {len(ms2)}; with peaks below m/z 200 (low mass recorded): {len(low)}")

    for label, tol in (("high resolution", None), ("ion trap", args.da)):
        counts = Counter()
        z_ox = Counter()
        z_all = Counter()
        for s in low:
            mz, it = s.get_peaks()
            z = s.getPrecursors()[0].getCharge() if s.getPrecursors() else 0
            z_all[z] += 1
            hit = {n: has(mz, it, t, tol if tol else t * args.ppm * 1e-6, args.min_rel)
                   for n, t in IONS.items()}
            for n, h in hit.items():
                counts[n] += h
            if hit["HexNAc oxonium 204.0867"]:
                z_ox[z] += 1
                counts["oxonium and Y1 342"] += hit["Y1 GlcNAc+2-AB 342.1660"]
        tol_txt = f"{args.ppm:g} ppm" if tol is None else f"{tol:g} Da"
        print(f"\n{label} tolerance ({tol_txt}), ion >= {args.min_rel:.0%} of base peak:")
        for n in IONS:
            print(f"  {n:<28} {counts[n]:6d} spectra")
        n_ox = counts["HexNAc oxonium 204.0867"]
        if n_ox:
            print(f"  oxonium-positive spectra that also show Y1 342.166: "
                  f"{counts['oxonium and Y1 342']} of {n_ox} "
                  f"({counts['oxonium and Y1 342'] / n_ox:.1%})")
        print(f"  precursor charge, oxonium-positive: {dict(sorted(z_ox.items()))}")
        print(f"  precursor charge, all low-mass MS2: {dict(sorted(z_all.items()))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
