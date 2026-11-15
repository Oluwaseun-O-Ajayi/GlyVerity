"""Inventory an Agilent .d folder and plot fluorescence against the MS TIC.

Writes, into the output folder:
  inventory.tsv   every data file found (name, detector, points, time range)
  fld_<name>.csv  each fluorescence trace (time, signal)
  ms_tic.csv      MS total ion current from MSPeak.bin (time, TIC)
  fld_vs_tic.pdf  two-panel figure on a shared time axis (also .png at 600 dpi)

Usage
-----
    python scripts/inspect_agilent_d.py RUN.d OUT_DIR
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from glyverity.io_agilent import read_d


def write_inventory(items, path: Path) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["name", "detector", "n_points", "t_start", "t_end"])
        for it in items:
            w.writerow([it.name, it.detector, it.n_points, f"{it.t_start:.4f}", f"{it.t_end:.4f}"])


def write_trace(trace, path: Path, value_name: str) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time", value_name])
        for t, v in zip(trace.times, trace.values):
            w.writerow([f"{t:.5f}", f"{v:.6g}"])


def _window(trace, lo, hi):
    if lo is None or hi is None:
        return trace.times, trace.values
    m = (trace.times >= lo) & (trace.times <= hi)
    return trace.times[m], trace.values[m]


def plot(fld, tic, out_base: Path, lo=None, hi=None) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(3.5, 3.2), constrained_layout=True)
    for tr in fld:
        unit = tr.metadata.get("unit", "")
        x, y = _window(tr, lo, hi)
        axes[0].plot(x, y, color="black", linewidth=0.7)
        axes[0].set_ylabel(f"Fluorescence ({unit})" if unit else "Fluorescence")
    x, y = _window(tic, lo, hi)
    axes[1].plot(x, y, color="black", linewidth=0.7)
    axes[1].set_ylabel("MS TIC (counts)")
    axes[1].set_xlabel("Retention time (min)")
    for ax, letter in zip(axes, "AB"):
        ax.text(-0.30, 1.05, letter, transform=ax.transAxes, fontweight="bold", va="top")
        ax.spines[["top", "right"]].set_visible(False)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    fig.savefig(out_base.with_suffix(".pdf"))
    fig.savefig(out_base.with_suffix(".png"), dpi=600)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run", type=Path, help="Agilent .d folder")
    ap.add_argument("out", type=Path, help="output folder")
    ap.add_argument("--window", nargs=2, type=float, metavar=("START", "END"),
                    help="plot only this retention-time window (min)")
    args = ap.parse_args(argv)
    lo, hi = args.window if args.window else (None, None)

    if not args.run.is_dir():
        ap.error(f"not a folder: {args.run}")
    args.out.mkdir(parents=True, exist_ok=True)

    items, fld, tic = read_d(str(args.run))
    write_inventory(items, args.out / "inventory.tsv")
    print(f"{len(items)} data files:")
    for it in items:
        print(f"  {it.name:<20} {it.detector:<6} {it.n_points:>7} points  "
              f"{it.t_start:8.3f} to {it.t_end:8.3f}")

    for tr in fld:
        write_trace(tr, args.out / f"fld_{Path(tr.name).stem}.csv", "fluorescence")
    if tic is not None:
        write_trace(tic, args.out / "ms_tic.csv", "tic")

    print(f"FLD traces: {len(fld)}; MS TIC: {'yes' if tic is not None else 'no'}")
    if tic is not None:
        md = tic.metadata
        print(f"scan record fields: {md.get('scan_fields')}")
        print(f"MS level counts: {md.get('ms_level_counts', 'no MSLevel field')}")
        print(f"TIC from MS1 scans only: {md.get('ms1_filtered')} "
              f"({tic.times.size} of {md.get('n_scans_total')} scans)")
    if fld and tic is not None:
        plot(fld, tic, args.out / "fld_vs_tic", lo, hi)
        print(f"figure: {args.out / 'fld_vs_tic.pdf'}")
    else:
        print("figure skipped: needs at least one FLD trace and the MS TIC")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
