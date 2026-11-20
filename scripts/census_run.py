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
from glyverity.annotate import composition_string
from glyverity.census import (
    CLASSES,
    FAMILIES,
    absorption_stats,
    assign_species,
    classify,
    family_cutoffs,
    group_species,
    reclassify_oxonium,
    repair_charges,
    target_table,
    unassigned_report,
)
from glyverity.fdr import best_errors, decoy_libraries, exchangeability, odd_h_shifts, qvalues
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
    ap.add_argument("--rt-tol", type=float, default=0.25, help="co-elution tolerance (min)")
    ap.add_argument("--max-shift", type=float, default=6.0,
                    help="largest retention shift between a label and a side-product form (min)")
    ap.add_argument("--mzml", action="store_true", help="also write the MS1 window as mzML")
    ap.add_argument("--decoy-shifts", default="0.23,0.41,0.59,0.77",
                    help="comma-separated Da shifts for decoy grouping runs ('' to skip)")
    ap.add_argument("--relation-fdr", type=float, default=0.05,
                    help="largest chance share allowed per relation family")
    ap.add_argument("--series-detail", action="store_true",
                    help="print every label/side-product pair with its fluorescence values")
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
    rec: list = []
    species = group_species(feats, ppm=args.ppm, rt_tol=args.rt_tol, record=rec)
    print(f"{exp.getNrSpectra()} MS1 scans -> {n_raw} features ({n_fix} split envelopes "
          f"re-merged, {len(feats)} after repair) -> {len(species)} species at {args.ppm:g} ppm")
    total = sum(f.intensity for f in feats)
    shifts = [float(x) for x in args.decoy_shifts.split(",") if x.strip()]
    if shifts:
        tgt = absorption_stats(species)
        # one decoy family at a time, the other families real (same candidate pool as target)
        drecs = [[] for _ in shifts]
        for fam in FAMILIES:
            for i, d in enumerate(shifts):
                r: list = []
                group_species(feats, ppm=args.ppm, rt_tol=args.rt_tol, decoy_shift=d,
                              decoy_families=[fam], record=r)
                drecs[i].extend(x for x in r if x[0] == fam)
        cut = family_cutoffs(rec, drecs, fdr=args.relation_fdr)
        print(f"relation families at {args.ppm:g} ppm / {args.rt_tol:g} min (target / decoy mean / "
              f"chance share; decoy shifts one family at a time) and score cutoff for <= "
              f"{args.relation_fdr:.0%} chance (score = max(ppm/{args.ppm:g}, dRT/{args.rt_tol:g})):")
        for fam in FAMILIES:
            tn = tgt[f"n:{fam}"]
            dm = float(np.mean([sum(1 for f, _ in r if f == fam) for r in drecs]))
            share = f"{dm / tn:6.1%}" if tn else "     -"
            print(f"      {fam:<9} {tn:5d} {dm:7.1f}  {share}   -> {cut[fam]:.3f} "
                  f"({cut[fam] * args.ppm:.1f} ppm, {cut[fam] * args.rt_tol:.3f} min)")
        rec = []
        species = group_species(feats, ppm=args.ppm, rt_tol=args.rt_tol, family_cut=cut,
                                record=rec)
        nd_fam = 0.0
        for fam in FAMILIES:
            counts = []
            for d in shifts:
                r = []
                group_species(feats, ppm=args.ppm, rt_tol=args.rt_tol, decoy_shift=d,
                              decoy_families=[fam], family_cut=cut, record=r)
                counts.append(sum(1 for f, _ in r if f == fam))
            nd_fam += float(np.mean(counts))
        nt = len(rec)
        tgt = absorption_stats(species)
        print(f"after family tolerances: {len(species)} species; {nt} absorptions "
              f"({100 * tgt['intensity'] / total:.1f}% of intensity), decoy sum over families "
              f"{nd_fam:.1f} -> estimated chance share {nd_fam / max(nt, 1):.1%}")

    groups = {args.label: ch.LABELS[args.label].added_mass, "free": 0.0, "C3H9N-RA": TAG_C3H9N,
              "NH3-RA": NH3_RA}
    comps = [c for c in an.enumerate_compositions(SEARCH)
             if c.get("HexNAc", 0) >= 2 and c.get("Hex", 0) >= 3
             and c.get("NeuAc", 0) + c.get("NeuGc", 0) <= c["HexNAc"] - 1]
    asg = assign_species(species, comps, groups, ppm=args.ppm)

    # composition q-values: odd-hydrogen decoys (octet-rule violations), separate search
    tmass = np.array([m for m, _, _ in target_table(comps, groups)])
    shifts_h = odd_h_shifts()
    libs, removed = decoy_libraries(tmass, shifts_h, ppm=args.ppm)
    smass = np.array([s.neutral_mass for s in species])
    te = best_errors(smass, np.sort(tmass))
    des = [best_errors(smass, lib) for lib in libs]
    scale = [len(tmass) / max(len(lib), 1) for lib in libs]
    qv = qvalues(te, des, max_ppm=args.ppm, scale=scale)
    q_of = {s.species_id: (te[i], qv[i]) for i, s in enumerate(species)}
    n_t, n_d, ratio = exchangeability(te, des, (15.0, 30.0), scale=scale)
    print(f"\ncomposition FDR: {len(tmass)} target masses, {len(libs)} odd-H decoy libraries "
          f"(+/-3,5,7,9 H; {sum(removed)} decoys within {args.ppm:g} ppm of a target removed)")
    print(f"  exchangeability, best error 15-30 ppm: target {n_t}, decoy mean {n_d:.1f}, "
          f"ratio {ratio:.2f} (1 expected; decoy counts scaled for removed decoys)")
    tsorted = sorted(target_table(comps, groups), key=lambda r: r[0])
    tm = np.array([r[0] for r in tsorted])
    band = [i for i in np.argsort(-np.array([s.intensity for s in species]))
            if 15.0 < te[i] <= 30.0]
    if band:
        print("  target matches in the 15-30 ppm band (strongest first):")
        for i in band[:12]:
            s = species[i]
            j = int(np.argmin(np.abs(tm - s.neutral_mass)))
            m, iso, grp = tsorted[j]
            zs = ",".join(sorted({str(f.charge) for f, form in s.members if form == "H"}))
            nm = composition_string(iso[0])
            print(f"    {s.neutral_mass:10.4f} {s.rt:7.3f} min z={zs:<5} {s.intensity:10.0f} "
                  f"{(s.neutral_mass - m) / m * 1e6:+6.1f} ppm  {nm}({grp})")
    sp_int = np.array([s.intensity for s in species])
    tot_int = sum(f.intensity for f in feats)
    for label, mask in (("target match", te <= args.ppm), ("q <= 0.05", qv <= 0.05),
                        ("q <= 0.01", qv <= 0.01)):
        print(f"  species with {label:<12}: {int(mask.sum()):4d} of {len(species)}, "
              f"{100 * sp_int[mask].sum() / tot_int:5.1f}% of feature intensity")
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
        trace = (np.asarray(fl[0].times) + args.lag, np.asarray(fl[0].values)) if fl else None
        layer = LabelSideProductSeries(args.label, test, max_shift=args.max_shift,
                                       fld_trace=trace)
        ev = layer.evaluate([c for c in cands if not c.flags])
        s = layer.summary
        fl_txt = (f"fluorescence fraction {s.fluorescence_fraction:.2f} "
              f"({s.n_fluorescence_pairs} clear pairs)" if s.n_fluorescence_pairs
              else "fluorescence untested")
        print(f"series {test} vs {args.label}: n={s.n}, ratio {s.ratio_median:.3f}, "
              f"spread {s.ratio_spread:.3f}, shift {s.rt_shift_median:+.2f} min, {fl_txt} -> "
              f"{'SERIES' if s.is_series else 'not a series'}")
        if args.series_detail and layer.pairs:
            print(f"    FLD peaks detected: {layer.n_fld_peaks}; tail shape "
                  f"{'measured' if layer.profile else 'not measurable (no isolated peaks), fixed clearance used'}")
            if layer.profile:
                sides = ", ".join(f"{'trailing' if k > 0 else 'leading'} from {v[2]} peaks"
                                  for k, v in sorted(layer.profile.items()))
                print(f"    tail correction: {layer.n_fld_peaks} FLD peaks; shape {sides}")
            print("    composition                     ref_rt  test_rt  ratio   FLD_ref  "
                  "FLD_test  FLD_expected  FLD_tail  used")
            for r in layer.pairs:
                print(f"    {r['name']:<30} {r['ref_rt']:7.3f} {r['test_rt']:8.3f} {r['ratio']:6.3f} "
                      f"{r['fld_ref']:9.3g} {r['fld_test']:9.3g} {r['fld_expected']:13.3g} "
                      f"{r['fld_tail']:9.3g}  "
                      f"{'yes' if r['used'] else 'no'}"
                      f"{'  (model fail)' if r['model_fail'] else ''}")
        for e in ev:
            by_key[e.feature_id].flags.append(f"side_product:{test}")

    cls = classify(asg, args.label)
    n_ox = reclassify_oxonium(asg, cls, ppm=args.ppm)
    print(f"oxonium ions reclassified as in-source fragments: {n_ox}")
    total = sum(f.intensity for f in feats)
    print(f"\nfeature intensity by class (n features, share of {total:.4g}):")
    for k in CLASSES:
        sub = [f for f in feats if cls.get(f.feature_id) == k]
        share = 100 * sum(f.intensity for f in sub) / total if total else 0.0
        print(f"  {k:24s} {len(sub):5d}  {share:5.1f}%")

    sp_cls = {a.species.species_id: cls[a.species.members[0][0].feature_id] for a in asg}
    print("\nspecies intensity by class of its own ions, split at q = 0.05:")
    for k in CLASSES:
        ids = [sid for sid, c in sp_cls.items() if c == k]
        if not ids:
            continue
        sub = [s for s in species if s.species_id in ids]
        ok = sum(s.intensity for s in sub if q_of[s.species_id][1] <= 0.05)
        rest = sum(s.intensity for s in sub) - ok
        print(f"  {k:24s} q<=0.05 {100 * ok / tot_int:5.1f}%   other {100 * rest / tot_int:5.1f}%")

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
                    "class", "matches", "species_target_ppm", "species_q"])
        for f in sorted(feats, key=lambda f: -f.intensity):
            sid, form = sp_of[f.feature_id]
            w.writerow([f.feature_id, f"{f.mz:.4f}", f.charge, f"{f.rt:.3f}", f"{f.intensity:.0f}",
                        sid, form, cls[f.feature_id], best[sid], f"{q_of[sid][0]:.2f}",
                        "" if np.isnan(q_of[sid][1]) else f"{q_of[sid][1]:.4f}"])
    print(f"table: {args.out / 'census_features.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
