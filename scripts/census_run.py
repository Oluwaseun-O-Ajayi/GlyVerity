"""Feature-first census of one Agilent run: features, species, assignments, layers, shares.

Steps
  1. MS1 scans in the window go into an OpenMS experiment (optionally written as mzML).
  2. OpenMS FeatureFindingMetabo finds isotope-envelope features.
  3. Features are grouped into species (charge states, Na+/K+/NH4+ adduct ions).
  4. Species are matched to compositions x reducing-end groups.
  5. Layers run on the matches: adduct/isotope isobars, in-source fragments, side-product series.
  6. Every feature's intensity is assigned to one class; shares are printed.

Usage
-----
    python scripts/census_run.py RUN.d OUT_DIR --label 2-AB --window 15 50 [--noise 100] [--mzml]
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np

from glyverity import annotate as an
from glyverity import chemistry as ch
from glyverity.census import (
    CLASSES,
    assign_species,
    classify,
    group_species,
    reclassify_oxonium,
    repair_charges,
    unassigned_report,
)
from glyverity.features import find_features, scans_from_centroid, to_msexperiment, write_mzml
from glyverity.io_agilent import fld_traces, ms1_mask, scan_records
from glyverity.layers import AdductIsotopeIsobar, InSourceFragment, LabelSideProductSeries

TAG_C3H9N = ch.formula_mass({"C": 3, "H": 9, "N": 1}) - ch.ATOMIC_MASS["O"]
NH3_RA = ch.formula_mass({"N": 1, "H": 3}) - ch.ATOMIC_MASS["O"]  # amino-alditol (hypothesis)
SEARCH = {"Hex": (3, 9), "HexNAc": (2, 6), "dHex": (0, 1), "NeuAc": (0, 4), "NeuGc": (0, 1)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--label", default="2-AB", choices=["2-AB", "ProA"])
    ap.add_argument("--window", nargs=2, type=float, default=(15.0, 50.0))
    ap.add_argument("--ppm", type=float, default=10.0)
    ap.add_argument("--noise", type=float, default=100.0)
    ap.add_argument("--lag", type=float, default=0.0, help="MS minus FLD time (min)")
    ap.add_argument("--rt-tol", type=float, default=0.2, help="co-elution tolerance (min)")
    ap.add_argument("--max-shift", type=float, default=6.0,
                    help="largest retention shift between a label and a side-product form (min)")
    ap.add_argument("--mzml", action="store_true", help="also write the MS1 window as mzML")
    args = ap.parse_args(argv)

    import rainbow

    args.out.mkdir(parents=True, exist_ok=True)
    files = list(rainbow.read(str(args.run), centroid=True).datafiles)
    ms = next(df for df in files if str(df.detector) == "MS" and hasattr(df, "scan"))
    keep = ms1_mask(scan_records(os.path.join(str(args.run), "AcqData")))
    lo, hi = args.window
    exp = to_msexperiment(list(scans_from_centroid(ms, keep, lo, hi)))
    if args.mzml:
        path = args.out / "ms1_window.mzML"
        write_mzml(exp, str(path))
        print(f"mzML: {path}")

    feats = find_features(exp, ppm=args.ppm, noise=args.noise)
    n_raw = len(feats)
    feats, n_fix = repair_charges(feats, ppm=args.ppm)
    species = group_species(feats, ppm=args.ppm, rt_tol=args.rt_tol)
    print(f"{exp.getNrSpectra()} MS1 scans -> {n_raw} features ({n_fix} split envelopes "
          f"re-merged, {len(feats)} after repair) -> {len(species)} species")

    groups = {args.label: ch.LABELS[args.label].added_mass, "free": 0.0, "C3H9N-RA": TAG_C3H9N,
              "NH3-RA": NH3_RA}
    comps = [c for c in an.enumerate_compositions(SEARCH)
             if c.get("HexNAc", 0) >= 2 and c.get("Hex", 0) >= 3
             and c.get("NeuAc", 0) + c.get("NeuGc", 0) <= c["HexNAc"] - 1]
    asg = assign_species(species, comps, groups, ppm=args.ppm)
    cands = [c for a in asg for c in a.matches]
    by_key = {c.key: c for c in cands}
    fl = fld_traces(files)
    if fl:
        for c in cands:
            c.fld_at_apex = float(np.interp(c.apex_time - args.lag, fl[0].times, fl[0].values))

    for e in AdductIsotopeIsobar(ppm=args.ppm, rt_tol=args.rt_tol).evaluate(cands):
        by_key[e.feature_id].flags.append(f"isobar:{e.detail['relation']} of {e.detail['by']}")
    for e in InSourceFragment(ppm=args.ppm, rt_tol=args.rt_tol).evaluate(cands):
        by_key[e.feature_id].flags.append(f"isf:-{e.detail['loss']} of {e.detail['by']}")
    for test in [g for g in groups if g != args.label]:
        layer = LabelSideProductSeries(args.label, test, max_shift=args.max_shift)
        ev = layer.evaluate([c for c in cands if not c.flags])
        s = layer.summary
        print(f"series {test} vs {args.label}: n={s.n}, ratio {s.ratio_median:.3f}, "
              f"spread {s.ratio_spread:.3f}, shift {s.rt_shift_median:+.2f} min, "
              f"fluorescence fraction {s.fluorescence_fraction:.2f} -> "
              f"{'SERIES' if s.is_series else 'not a series'}")
        for e in ev:
            by_key[e.feature_id].flags.append(f"side_product:{test}")

    cls = classify(asg, args.label)
    n_ox = reclassify_oxonium(asg, cls, ppm=args.ppm, rt_tol=args.rt_tol)
    print(f"oxonium ions reclassified as in-source fragments: {n_ox}")
    total = sum(f.intensity for f in feats)
    print(f"\nfeature intensity by class (n features, share of {total:.4g}):")
    for k in CLASSES:
        sub = [f for f in feats if cls.get(f.feature_id) == k]
        share = 100 * sum(f.intensity for f in sub) / total if total else 0.0
        print(f"  {k:24s} {len(sub):5d}  {share:5.1f}%")

    print("\nstrongest unassigned features (delta = neutral minus co-eluting assigned species):")
    for f, delta, near in unassigned_report(asg, cls):
        d = f"{delta:+10.4f}" if np.isfinite(delta) else "         -"
        print(f"  {f.mz:10.4f} z={f.charge} {f.rt:7.3f} min {f.intensity:11.0f}  {d}  {near}")

    sp_of = {f.feature_id: (a.species.species_id, form) for a in asg for f, form in a.species.members}
    best = {a.species.species_id: "; ".join(f"{c.name}({c.end_group}){'[' + ','.join(c.flags) + ']' if c.flags else ''}"
                                            for c in a.matches) for a in asg}
    with (args.out / "census_features.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["feature_id", "mz", "charge", "rt_min", "intensity", "species", "ion_form",
                    "class", "matches"])
        for f in sorted(feats, key=lambda f: -f.intensity):
            sid, form = sp_of[f.feature_id]
            w.writerow([f.feature_id, f"{f.mz:.4f}", f.charge, f"{f.rt:.3f}", f"{f.intensity:.0f}",
                        sid, form, cls[f.feature_id], best[sid]])
    print(f"table: {args.out / 'census_features.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
