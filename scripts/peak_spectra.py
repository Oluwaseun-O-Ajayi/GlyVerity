"""List the strongest MS1 peaks in retention-time windows and match them to glycan compositions.

For each window the MS1 scans are summed, peaks within 10 ppm are merged, and
the strongest peaks are listed. Peaks that sit one 13C spacing above a more
intense peak are marked as isotopes. Every other peak is matched against
glycan compositions as free (unlabelled) and 2-AB labelled ions, with H+, Na+,
K+ or NH4+ charge carriers, charges 1 to 3.

Output: one TSV per window in OUT_DIR, plus a printed summary.

Usage
-----
    python scripts/peak_spectra.py RUN.d OUT_DIR --windows 30.4-31.0 35.5-36.1 --top 15
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

from glyverity import annotate as an
from glyverity.io_agilent import ms1_mask, scan_records

RANGES = {"Hex": (0, 9), "HexNAc": (0, 7), "dHex": (0, 2), "NeuAc": (0, 4), "NeuGc": (0, 2)}


def parse_window(text: str) -> tuple[float, float]:
    lo, hi = (float(x) for x in text.split("-"))
    if hi <= lo:
        raise argparse.ArgumentTypeError(f"window end must exceed start: {text}")
    return lo, hi


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--windows", nargs="+", type=parse_window, required=True,
                    help="retention-time windows in MS time, e.g. 30.4-31.0")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--ppm", type=float, default=20.0, help="matching tolerance")
    args = ap.parse_args(argv)

    import rainbow

    args.out.mkdir(parents=True, exist_ok=True)
    directory = rainbow.read(str(args.run), centroid=True)
    ms = [df for df in directory.datafiles if str(df.detector) == "MS" and hasattr(df, "scan")]
    if not ms:
        print("no centroided MS data in this run")
        return 1
    records = scan_records(os.path.join(str(args.run), "AcqData"))
    keep = ms1_mask(records)
    comps = [c for c in an.enumerate_compositions(RANGES) if c.get("HexNAc", 0) >= 1]
    labels = ["free", "2-AB"]

    for lo, hi in args.windows:
        mz, inten, n_scans = an.window_spectrum(ms[0], keep, lo, hi)
        order = mz.size and sorted(range(mz.size), key=lambda k: -inten[k])[: args.top]
        path = args.out / f"window_{lo:.2f}-{hi:.2f}.tsv"
        print(f"\n=== {lo:.2f} to {hi:.2f} min: {n_scans} MS1 scans, {mz.size} merged peaks ===")
        with path.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["rank", "mz", "intensity", "isotope_of_mz", "composition", "label",
                        "charge", "carriers", "theoretical_mz", "ppm_error"])
            for rank, i in enumerate(order or [], start=1):
                iso = an.isotope_parent(mz, inten, i)
                if iso is not None:
                    j, z = iso
                    w.writerow([rank, f"{mz[i]:.4f}", f"{inten[i]:.0f}", f"{mz[j]:.4f} (z={z})",
                                "", "", "", "", "", ""])
                    print(f"{rank:>3} {mz[i]:10.4f} {inten[i]:12.0f}  isotope of {mz[j]:.4f} (z={z})")
                    continue
                matches = an.match_mz(mz[i], comps, labels, ppm=args.ppm)
                if not matches:
                    w.writerow([rank, f"{mz[i]:.4f}", f"{inten[i]:.0f}", "", "no match", "", "", "",
                                "", ""])
                    print(f"{rank:>3} {mz[i]:10.4f} {inten[i]:12.0f}  no match")
                    continue
                for m in matches:
                    w.writerow([rank, f"{mz[i]:.4f}", f"{inten[i]:.0f}", "", m.composition, m.label,
                                m.charge, m.carriers, f"{m.theoretical_mz:.4f}",
                                f"{m.ppm_error:+.1f}"])
                best = matches[0]
                more = f" (+{len(matches) - 1} more)" if len(matches) > 1 else ""
                print(f"{rank:>3} {mz[i]:10.4f} {inten[i]:12.0f}  {best.composition} {best.label} "
                      f"[M+{best.carriers}]{best.charge}+ {best.ppm_error:+.1f} ppm{more}")
        print(f"table: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
