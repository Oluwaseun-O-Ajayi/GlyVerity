# Public glycomics dataset catalog

`datasets.csv` lists candidate datasets for the census and for GlyVerity Bench.
It is a seed list. The systematic harvest runs on the cluster with
`scripts/harvest_pride.py` and by browsing GlycoPOST and MassIVE (see below).

## Columns

| Column | Meaning |
|---|---|
| dataset_id | Repository accession or DOI |
| repository | GlycoPOST, MassIVE, PRIDE, Zenodo, other |
| year | Year of deposition or publication |
| group_or_first_author | Lab or first author |
| instrument | Mass spectrometer (and IM device if any) |
| lc_mode | HILIC, PGC, RP, CE, nanoLC |
| label | none (native), 2-AB, 2-AA, ProA, RapiFluor-MS, InstantPC, permethylated, reduced |
| glycan_class | N-glycans, O-glycans, HMOs and so on |
| fluorescence_paired | yes, no, or consensus values; whether an independent quantitative reference exists |
| raw_format | vendor raw, mzML, mzXML |
| licence | Licence stated by the repository record. Record it before any reuse |
| in_scope | yes if released-glycan LC-MS (PhD-stage scope) |
| status | unverified, listed_in_paper, listed_in_report, verified_downloadable, embargoed, raw_not_found |
| source_checked | Where the row's facts came from |
| notes | Anything that affects use |

## Rules

* A row moves to `verified_downloadable` only after a file has been downloaded and opened.
* Never fill a cell from memory. Use `to_check` until confirmed.
* Record the licence for every dataset before it enters the census.

## Manual harvest (browser)

1. GlycoPOST: https://glycopost.glycosmos.org/ , browse public projects, filter for LC-MS.
2. MassIVE: https://massive.ucsd.edu/ , search keywords "glycan", "glycomics", "N-glycan".
3. For each hit, add a row with the columns above.
