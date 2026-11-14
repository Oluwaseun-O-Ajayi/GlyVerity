"""Harvest candidate released-glycan datasets from the PRIDE Archive search API.

Run on the cluster (the workspace that wrote this script has no access to EBI).
Writes the raw API responses and a flat CSV for manual screening.

Usage
-----
    python scripts/harvest_pride.py --out data/catalog/harvest

The endpoint and response fields are checked at run time. If the API shape
differs from what this script expects, the raw JSON is still saved so nothing
is lost, and the script reports which fields were missing.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

PRIDE_SEARCH = "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects"
KEYWORDS = [
    "glycomics",
    "N-glycan",
    "released glycans",
    "glycan profiling",
    "procainamide",
    "2-aminobenzamide",
    "RapiFluor",
    "permethylated glycans",
]
FIELDS = ["accession", "title", "submissionDate", "publicationDate", "instruments", "keywords"]


def fetch(url: str, retries: int = 3, wait: float = 2.0) -> object:
    """GET a URL and parse JSON, retrying on network errors."""
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return json.load(resp)
        except (OSError, ValueError) as exc:  # URLError, timeouts, bad JSON
            last = exc
            time.sleep(wait * (attempt + 1))
    raise RuntimeError(f"failed after {retries} attempts: {url}") from last


def search(keyword: str, page_size: int, max_pages: int) -> list[dict]:
    """Return all project records for one keyword."""
    records: list[dict] = []
    for page in range(max_pages):
        query = urllib.parse.urlencode({"keyword": keyword, "pageSize": page_size, "page": page})
        data = fetch(f"{PRIDE_SEARCH}?{query}")
        batch = data if isinstance(data, list) else data.get("_embedded", {}).get("projects", [])
        if not batch:
            break
        records.extend(batch)
        if len(batch) < page_size:
            break
    return records


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--page-size", type=int, default=100)
    parser.add_argument("--max-pages", type=int, default=20)
    args = parser.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    seen: dict[str, dict] = {}
    for kw in KEYWORDS:
        recs = search(kw, args.page_size, args.max_pages)
        (args.out / f"pride_{kw.replace(' ', '_')}.json").write_text(json.dumps(recs, indent=1))
        print(f"{kw}: {len(recs)} records")
        for r in recs:
            acc = r.get("accession")
            if acc:
                seen.setdefault(acc, r).setdefault("_matched", set()).add(kw)

    missing = sorted({f for r in seen.values() for f in FIELDS if f not in r})
    if missing:
        print(f"note: fields absent from some records: {missing}", file=sys.stderr)

    out_csv = args.out / "pride_candidates.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(FIELDS + ["matched_keywords"])
        for acc in sorted(seen):
            r = seen[acc]
            row = [json.dumps(r.get(f)) if isinstance(r.get(f), (list, dict)) else r.get(f, "")
                   for f in FIELDS]
            w.writerow(row + [";".join(sorted(r["_matched"]))])
    print(f"{len(seen)} unique projects written to {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
