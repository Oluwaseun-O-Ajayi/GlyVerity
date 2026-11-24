"""Figure: MS shares against FLD shares before and after cleaning, four runs.

Panels A-D: one run each. Each point is one FLD peak; x = FLD share, y = MS share, open
orange = protonated ions only, filled blue = all ions assigned to the species. The dashed line
is y = x. Axes use a square-root scale so small and large peaks are both readable.
Panel E: TV(protonated) - TV(all assigned) with 95% bootstrap intervals over FLD peaks
(2000 resamples, seed 0, as in census_run.py); positive = cleaning moves MS toward FLD.
The shuffle-control p-values need raw intensities, which tables written before 2026-11-23
lack; they are reported from the census logs, not recomputed here.

Inputs are the fld_validation.tsv tables written by census_run.py.
Usage: python scripts/fig_fld_validation.py <out_prefix>
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from glyverity.validate import bootstrap_tv_difference, total_variation

HERE = Path(__file__).resolve().parents[1] / "data" / "results" / "fld_validation"
RUNS = [
    ("fetuin_2ab.tsv", "Fetuin, 2-AB", "Agilent Q-TOF, lab 1"),
    ("ifnb_2ab.tsv", "IFN-β, 2-AB", "Orbitrap Fusion, lab 2"),
    ("ifnb_proa_a.tsv", "IFN-β, ProA (a)", "Orbitrap Fusion, lab 2"),
    ("ifnb_proa_b.tsv", "IFN-β, ProA (b)", "Orbitrap Fusion, lab 2"),
]
BLUE, ORANGE = "#2a78d6", "#eb6834"  # validated pair (dataviz validator, light surface)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def load(name: str) -> dict[str, np.ndarray]:
    with open(HERE / name, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {k: np.array([float(r[k]) for r in rows])
            for k in ("fld_share", "ms_protonated_share", "ms_all_share")}


def main() -> int:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "fig_fld_validation")
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
        "ytick.color": INK2, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(7.0, 4.4))  # two-column width
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.72], hspace=0.42, wspace=0.42)
    stats = []
    for i, (fn, title, sub) in enumerate(RUNS):
        d = load(fn)
        ax = fig.add_subplot(gs[0, i])
        lim = 0.36
        ax.plot([0, lim], [0, lim], ls=(0, (3, 2)), lw=0.7, color=INK2, zorder=1)
        ax.scatter(d["fld_share"], d["ms_protonated_share"], s=14, facecolors="none",
                   edgecolors=ORANGE, linewidths=0.8, zorder=3, label="Protonated ions")
        ax.scatter(d["fld_share"], d["ms_all_share"], s=14, color=BLUE, edgecolors="#fcfcfb",
                   linewidths=0.5, zorder=4, label="All assigned ions")
        for a in (ax.set_xscale, ax.set_yscale):
            a("function", functions=(np.sqrt, np.square))
        ticks = [0, 0.01, 0.05, 0.1, 0.2, 0.3]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels(["0", ".01", ".05", ".1", ".2", ".3"])
        ax.set_yticklabels(["0", ".01", ".05", ".1", ".2", ".3"])
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect("equal")
        ax.grid(color=GRID, lw=0.5, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        n = len(d["fld_share"])
        tv_p = total_variation(d["ms_protonated_share"], d["fld_share"])
        tv_a = total_variation(d["ms_all_share"], d["fld_share"])
        ax.set_title(f"{title}\n{sub}\n{n} FLD peaks\nTV {tv_p:.3f} \u2192 {tv_a:.3f}",
                     fontsize=6.5, color=INK, loc="left", pad=4, linespacing=1.3)
        ax.set_xlabel("FLD share")
        if i == 0:
            ax.set_ylabel("MS share")
        n = len(d["fld_share"])
        tv_p = total_variation(d["ms_protonated_share"], d["fld_share"])
        tv_a = total_variation(d["ms_all_share"], d["fld_share"])
        ax.text(-0.02, 1.47, "ABCD"[i], transform=ax.transAxes, fontsize=9,
                fontweight="bold", color=INK, va="bottom")
        diff, lo, hi = bootstrap_tv_difference(d["ms_protonated_share"], d["ms_all_share"],
                                               d["fld_share"])
        stats.append((title, n, diff, lo, hi))
        if i == 0:
            handles, labels = ax.get_legend_handles_labels()

    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=2,
               frameon=False, fontsize=7, handletextpad=0.3, columnspacing=1.5)
    ax = fig.add_subplot(gs[1, :4])
    y = np.arange(len(stats))[::-1]
    for yy, (title, n, diff, lo, hi) in zip(y, stats):
        ax.plot([lo, hi], [yy, yy], color=BLUE, lw=2, solid_capstyle="round")
        ax.scatter([diff], [yy], s=28, color=BLUE, edgecolors="#fcfcfb", linewidths=0.8,
                   zorder=3)
        ax.text(0.215, yy, f"{diff:+.3f}  [{lo:+.3f}, {hi:+.3f}]",
                va="center", fontsize=6, color=INK2)
    ax.axvline(0, color=INK2, lw=0.7, ls=(0, (3, 2)))
    ax.set_yticks(y)
    ax.set_yticklabels([f"{t} (n = {n})" for t, n, *_ in stats])
    ax.set_xlim(-0.15, 0.32)
    ax.set_xlabel("TV(protonated) − TV(all assigned), 95% bootstrap interval")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.text(-0.02, 1.02, "E", transform=ax.transAxes, fontsize=9, fontweight="bold",
            color=INK, ha="right")
    top = len(stats) - 0.35
    ax.text(-0.006, top, "\u2190 away from FLD", fontsize=6, color=INK2, ha="right", va="bottom")
    ax.text(0.006, top, "toward FLD \u2192", fontsize=6, color=INK2, ha="left", va="bottom")
    ax.set_ylim(-0.5, len(stats) - 0.05)

    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{out}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight",
                    facecolor="white")
    for title, n, diff, lo, hi in stats:
        print(f"{title:<22} n={n:3d}  diff {diff:+.3f}  [{lo:+.3f}, {hi:+.3f}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
