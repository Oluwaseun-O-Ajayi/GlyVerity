"""Extract targeted candidates from one run and apply the artifact evidence layers.

Candidates: every composition in the search space x every reducing-end group.
Layers:
  AdductIsotopeIsobar        candidates explained by Na+/K+/NH4+ adducts or 13C peaks of a
                             co-eluting stronger candidate
  LabelSideProductSeries     each non-reference end group tested as a side-product series of
                             the reference label

Outputs candidates.tsv (one row per candidate with signal, with flags) and a printed summary.

Usage
-----
    python scripts/run_layers.py RUN.d OUT_DIR --label 2-AB --window 15 50 --lag 0.037
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

from glyverity import annotate as an
from glyverity import chemistry as ch
from glyverity.candidates import extract_candidates
from glyverity.io_agilent import fld_traces, ms1_mask, scan_records
from glyverity.layers import AdductIsotopeIsobar, LabelSideProductSeries, recurring_deltas

TAG_C3H9N = ch.formula_mass({"C": 3, "H": 9, "N": 1}) - ch.ATOMIC_MASS["O"]
SEARCH = {"Hex": (3, 9), "HexNAc": (2, 6), "dHex": (0, 1), "NeuAc": (0, 4), "NeuGc": (0, 1)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--label", default="2-AB", choices=["2-AB", "ProA"])
    ap.add_argument("--window", nargs=2, type=float, default=(15.0, 50.0))
    ap.add_argument("--lag", type=float, default=0.0)
    ap.add_argument("--ppm", type=float, default=10.0, help="XIC tolerance")
    ap.add_argument("--min-area", type=float, default=1000.0)
    ap.add_argument("--isobar-ppm", type=float, default=None,
                    help="adduct/isotope relation tolerance (default: same as --ppm)")
    args = ap.parse_args(argv)

    import rainbow

    args.out.mkdir(parents=True, exist_ok=True)
    files = list(rainbow.read(str(args.run), centroid=True).datafiles)
    ms = next(df for df in files if str(df.detector) == "MS" and hasattr(df, "scan"))
    fl = fld_traces(files)
    fld = (fl[0].times, fl[0].values) if fl else None
    keep = ms1_mask(scan_records(os.path.join(str(args.run), "AcqData")))

    groups = {args.label: ch.LABELS[args.label].added_mass, "free": 0.0, "C3H9N-RA": TAG_C3H9N}
    comps = [c for c in an.enumerate_compositions(SEARCH)
             if c.get("HexNAc", 0) >= 2 and c.get("Hex", 0) >= 3
             and c.get("NeuAc", 0) + c.get("NeuGc", 0) <= max(c["HexNAc"] - 2, 0)]
    from glyverity.candidates import group_isomeric
    print(f"{len(comps)} compositions ({len(group_isomeric(comps))} distinct formulas) x "
          f"{len(groups)} end groups")
    cands = extract_candidates(ms, keep, comps, groups, tuple(args.window), args.ppm, fld, args.lag)
    live = [c for c in cands if c.area >= args.min_area]
    print(f"{len(live)} candidates with area >= {args.min_area:g}")

    isobar = AdductIsotopeIsobar(ppm=args.isobar_ppm or args.ppm, rt_tol=0.15)
    by_key = {c.key: c for c in live}
    for e in isobar.evaluate(live):
        by_key[e.feature_id].flags.append(
            f"isobar:{e.detail['relation']} of {e.detail['by']} ({e.score:+.1f} ppm)")

    for test in [g for g in groups if g != args.label]:
        layer = LabelSideProductSeries(args.label, test)
        ev = layer.evaluate([c for c in live if not c.flags])
        s = layer.summary
        print(f"\nseries test {test} vs {args.label}: n={s.n}, ratio median {s.ratio_median:.3f}, "
              f"spread (MAD log10) {s.ratio_spread:.3f}, RT shift median {s.rt_shift_median:+.2f} "
              f"min, same sign {s.same_sign_fraction:.2f}, fluorescence fraction "
              f"{s.fluorescence_fraction:.2f} -> {'SERIES' if s.is_series else 'not a series'}")
        for e in ev:
            by_key[e.feature_id].flags.append(f"side_product:{test} of {e.detail['reference']}")

    path = args.out / "candidates.tsv"
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["composition", "end_group", "neutral_mass", "apex_min", "area", "fld_at_apex",
                    "flags", "isomeric_alternatives"])
        for c in sorted(live, key=lambda c: -c.area):
            w.writerow([c.name, c.end_group, f"{c.neutral_mass:.4f}", f"{c.apex_time:.3f}",
                        f"{c.area:.0f}", f"{c.fld_at_apex:.3f}", "; ".join(c.flags),
                        ", ".join(c.alternative_names)])

    print("\nsummary by end group (candidates, MS1 area):")
    tot_all = sum(c.area for c in live)
    for g in groups:
        sub = [c for c in live if c.end_group == g]
        a = sum(c.area for c in sub)
        iso = sum(c.area for c in sub if any(f.startswith("isobar") for f in c.flags))
        side = sum(c.area for c in sub if any(f.startswith("side_product") for f in c.flags)
                   and not any(f.startswith("isobar") for f in c.flags))
        print(f"  {g:9s} {len(sub):4d} candidates | area {a:12.0f} | adduct/isotope {iso:10.0f} | "
              f"side-product series {side:10.0f}")
    iso_all = sum(c.area for c in live if any(f.startswith("isobar") for f in c.flags))
    side_all = sum(c.area for c in live if any(f.startswith("side_product") for f in c.flags)
                   and not any(f.startswith("isobar") for f in c.flags))
    if tot_all:
        print(f"\nof {tot_all:.0f} candidate MS1 area: adduct/isotope {100 * iso_all / tot_all:.1f}%, "
              f"side-product series {100 * side_all / tot_all:.1f}%, "
              f"unflagged {100 * (tot_all - iso_all - side_all) / tot_all:.1f}%")

    print("\nrecurring mass differences between co-eluting unflagged candidates (top 10):")
    for r in recurring_deltas(live)[:10]:
        ex = ", ".join(f"{a.split('|')[0]}({a.split('|')[1]}) vs {b.split('|')[0]}({b.split('|')[1]})"
                       for a, b in r.pairs[:3])
        print(f"  {r.delta:+10.4f} Da  x{r.count}:  {ex}")
    print(f"table: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
