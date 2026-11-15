"""Screen PRIDE harvest records for released-glycan LC-MS datasets.

Reads the raw JSON files written by ``harvest_pride.py``, removes duplicates,
and flags each project from its title, description and protocols using
transparent keyword rules. The output is a triage list for manual review,
not a final classification.

Usage
-----
    python scripts/screen_pride.py data/catalog/harvest data/catalog/pride_screened.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

TEXT_FIELDS = ["title", "projectDescription", "sampleProcessingProtocol", "dataProcessingProtocol"]

# Each rule: flag name -> regex (case-insensitive) applied to the joined text.
RULES = {
    "release": r"pngase|endo-?h\b|glycan release|released (n-|o-)?glycans?|released glycans?"
    r"|free (n-|o-)?glycans?|beta-elimination|hydrazinolysis",
    "glycopeptide": r"glycopeptide|glyco-?proteom|pglyco|byonic|msfragger-glyco|glycoform site"
    r"|site-specific glycosylation|intact glyco",
    "label_2AB": r"\b2-?ab\b|2-aminobenzamide",
    "label_2AA": r"\b2-?aa\b|2-aminobenzoic|anthranilic acid",
    "label_ProA": r"procainamide",
    "label_RFMS": r"rapifluor",
    "label_InstantPC": r"instantpc|instant pc",
    "permethylated": r"permethyl",
    "reduced": r"\breduced (n-|o-)?glycans?|alditol|sodium borohydride|nabh4",
    "hilic": r"hilic",
    "pgc": r"porous graphiti[sz]ed carbon|\bpgc\b|hypercarb",
    "rp": r"\bc18\b|reversed[- ]phase",
    "ce": r"capillary electrophoresis|\bce-ms\b|cesi",
    "maldi": r"maldi",
    "fluorescence": r"fluorescen|\bflr\b|\bfld\b|uplc-flr",
}
COMPILED = {k: re.compile(v, re.IGNORECASE) for k, v in RULES.items()}
RAW_EXT = (".raw", ".d", ".wiff", ".wiff2", ".mzml", ".mzxml", ".d.zip", ".baf", ".tdf")


def load_records(folder: Path) -> dict[str, dict]:
    """Load and de-duplicate all PRIDE JSON records in ``folder``."""
    records: dict[str, dict] = {}
    for path in sorted(folder.glob("pride_*.json")):
        for rec in json.loads(path.read_text()):
            acc = rec.get("accession")
            if acc:
                records.setdefault(acc, rec)
    return records


def flags_for(rec: dict) -> dict[str, bool]:
    """Apply the keyword rules to one record."""
    text = " ".join(str(rec.get(f) or "") for f in TEXT_FIELDS)
    return {k: bool(rx.search(text)) for k, rx in COMPILED.items()}


def raw_file_count(rec: dict) -> int:
    """Count files whose names end in a recognised raw or open MS format."""
    names = rec.get("projectFileNames") or []
    return sum(1 for n in names if str(n).lower().endswith(RAW_EXT))


def category(fl: dict[str, bool]) -> str:
    """Triage category from flags."""
    labelled = any(fl[k] for k in ("label_2AB", "label_2AA", "label_ProA", "label_RFMS",
                                   "label_InstantPC", "permethylated", "reduced"))
    released = fl["release"] or labelled
    if released and not fl["glycopeptide"]:
        return "released_only"
    if released and fl["glycopeptide"]:
        return "mixed"
    if fl["glycopeptide"]:
        return "glycopeptide_only"
    return "unclear"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("harvest_dir", type=Path)
    ap.add_argument("out_csv", type=Path)
    args = ap.parse_args(argv)

    recs = load_records(args.harvest_dir)
    cols = ["accession", "category", "year", "instruments", "lab_pi", "n_raw_files",
            *RULES.keys(), "title"]
    rows = []
    for acc, rec in sorted(recs.items()):
        fl = flags_for(rec)
        rows.append({
            "accession": acc,
            "category": category(fl),
            "year": str(rec.get("submissionDate") or "")[:4],
            "instruments": "; ".join(rec.get("instruments") or []),
            "lab_pi": "; ".join(rec.get("labPIs") or []),
            "n_raw_files": raw_file_count(rec),
            **{k: int(v) for k, v in fl.items()},
            "title": rec.get("title", ""),
        })
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["category"]] = counts.get(r["category"], 0) + 1
    print(f"{len(rows)} projects -> {args.out_csv}")
    for k in sorted(counts):
        print(f"  {k}: {counts[k]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
