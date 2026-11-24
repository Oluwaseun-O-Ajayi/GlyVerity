"""Diagnostic ions in the MS2 spectra of chosen precursors (mzML input).

For each neutral mass given, MS2 spectra are selected whose precursor m/z matches
[M + zH]z+ (z = 1 to 4) within --ppm, and the listed fragment ions are reported with their
intensity relative to the spectrum's base peak (0 if absent within --frag-ppm).

Usage:
  python scripts/ms2_precursor_ions.py run.mzML --mass 1787.6822 2078.7786 2369.8735
      [--ions 204.0867 223.1288 265.1758 342.1660 369.1867] [--ppm 10] [--frag-ppm 20]
"""

from __future__ import annotations

import argparse

import numpy as np

PROTON = 1.007276466812
DEFAULT_IONS = {
    204.0867: "HexNAc oxonium",
    366.1395: "HexHexNAc oxonium",
    223.1288: "Y1 amino-alditol",
    369.1867: "Y1 amino-alditol + Fuc",
    265.1758: "Y1 C3H9N tag",
    342.1660: "Y1 2-AB",
    488.2239: "Y1 2-AB + Fuc",
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mzml")
    ap.add_argument("--mass", type=float, nargs="+", required=True, help="neutral masses (Da)")
    ap.add_argument("--ions", type=float, nargs="*", default=sorted(DEFAULT_IONS))
    ap.add_argument("--ppm", type=float, default=10.0, help="precursor tolerance")
    ap.add_argument("--frag-ppm", type=float, default=20.0, help="fragment tolerance")
    args = ap.parse_args()

    import pyopenms as oms

    exp = oms.MSExperiment()
    oms.MzMLFile().load(args.mzml, exp)
    ms2 = [s for s in exp.getSpectra() if s.getMSLevel() == 2 and s.getPrecursors()]
    print(f"{len(ms2)} MS2 spectra with a precursor")
    head = "  ".join(f"{x:9.4f}" for x in args.ions)
    for m in args.mass:
        print(f"\nneutral mass {m:.4f}")
        found = 0
        for s in ms2:
            p = s.getPrecursors()[0]
            pmz = p.getMZ()
            z = next((z for z in (1, 2, 3, 4)
                      if abs(pmz - (m + z * PROTON) / z) <= pmz * args.ppm * 1e-6), None)
            if z is None:
                continue
            mz, it = s.get_peaks()
            if it.size == 0:
                continue
            if found == 0:
                print(f"  RT min   precursor  z  base      {head}")
            rel = []
            for ion in args.ions:
                k = np.abs(mz - ion) <= ion * args.frag_ppm * 1e-6
                rel.append(it[k].max() / it.max() if k.any() else 0.0)
            print(f"  {s.getRT() / 60:6.2f}  {pmz:10.4f}  {z}  {it.max():8.3g}  "
                  + "  ".join(f"{r:9.3f}" for r in rel))
            found += 1
        if not found:
            print("  no MS2 spectrum selected this precursor")
    print("\nions: " + "; ".join(f"{x:.4f} {DEFAULT_IONS.get(x, '')}".strip() for x in args.ions))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
