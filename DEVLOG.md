# GlyVerity development log

## 2026-11-14

**Goal:** Start Phase 1: repo scaffold, prior-art check, claim verification, seed dataset catalog.

**Done**
* Package skeleton (`src/glyverity`), MIT, Python 3.10+, CI on 3.10 to 3.13 (ruff, pytest).
* `chemistry.py`: residue, label and charge-carrier masses from elemental formulas.
  17 tests pass, including Man5 = 1234.4334 Da, 2-AB = +120.0687 Da, ProA = +219.1735 Da.
* `evidence.py`: plugin interface and registry for evidence layers.
* `docs/phase1/prior_art.md`: 11 released-glycan tools, 4 adjacent precedents, 4 libraries.
* `docs/phase1/claims_verification.md`: 9 claims checked.
* `docs/phase1/go_no_go.md`: interim conditional go, refined gap statement.
* `data/catalog/datasets.csv`: 9 seed rows; `scripts/harvest_pride.py` for the cluster.

**Findings**
* MS1 composition decoys already exist (GlycReSoft 2012). The novelty is artifact
  pre-annotation (incl. label side-products) plus calibrated q-values, not FDR itself.
* GlycoGenius (Nat. Commun. 2025) is the closest active tool: MS1, adducts, isotopes,
  ProA and 2-AB, but threshold filters only, no decoys described, no ISF or label side-products.
* The khipu preprint was misstated in the brief: in metabolomics ISFs are under 10% of features.
* The FDA September 2025 guidance is "Development of Therapeutic Protein Biosimilars:
  Comparative Analytical Assessment and Other Quality-Related Considerations" (FDA-2019-D-2102).
* The +260.200 Da series in GPST000788 equals ProA plus acetonitrile (41.0265 Da).
* Repositories (Zenodo, GlycoPOST, PRIDE, MassIVE) are not reachable from the build workspace;
  harvesting runs on the cluster.

**Evening update**
* PRIDE harvest run on the cluster and pushed (474 unique projects).
* `scripts/screen_pride.py` added: keyword triage, then manual check of protocols and files.
* 17 in-scope PRIDE projects added to the catalog (2016 to 2026). Two run fluorescence and MS
  in the same LC run: PXD077477 (ProA and 2-AB) and PXD078897 (RapiFluor-MS).
* Most public data are PGC/reduced or C18/permethylated, not labelled HILIC. Decision: labelled
  HILIC is the validation core; the census includes all labels.

**Next**
* Download one glycan raw file each from PXD077477 and PXD078897; check for a fluorescence channel.
* Browse GlycoPOST and MassIVE; add rows to the catalog.
* Count fluorescence-paired and standards-based ground truth sets (go/no-go condition 1).
* Finish prior-art items marked "pending" and "not yet checked".

## 2026-11-15

**Goal:** Close the Phase 1 verification items and stage the first ground-truth data.

**Done**
* FDA September 2025 guidance wording checked against the PDF (pp. 11, 12, 14); claim reworded.
* Scope decision recorded: labelled HILIC is the validation core; the census includes all labels.
* Zenodo 4705546 (Smith 2021, CC BY 4.0): three runs downloaded to the cluster at
  `~/glyverity_data/raw/zenodo_4705546`; MD5 sums match `data/catalog/checksums/`.
* Both Agilent runs (fetuin, IgG) contain `FLD1.cg` (fluorescence) and `MSPeak.bin` (MS).

* Added `glyverity.io_agilent` (reads FLD traces and MSPeak.bin centroids via rainbow-api 1.5.5,
  LGPL-3.0, optional extra `agilent`) and `scripts/inspect_agilent_d.py` (inventory, CSV export,
  FLD vs TIC figure). Tested on stand-in data only; not yet run on the real files.

* Decision: no patent. GlyVerity stays open (MIT) for unrestricted adoption. Added CITATION.cff.

* First real run (fetuin Agilent): FLD and MS both read. FLD dominated by an early saturated
  peak (about 1 to 2 min, likely excess 2-AB); glycans elute about 18 to 46 min. MS TIC
  alternates scan to scan, consistent with interleaved MS1 and MS/MS scans.
* Reader now reads per-scan records, keeps MS1 scans when an MSLevel field exists, and reports
  the fields found. Script takes `--window START END`. rainbow-api pinned to 1.5.5.

* Fetuin Agilent run, window 15 to 50 min: 2050 MS1 of 9346 scans (7296 MS2). FLD and MS1 TIC
  share apexes near 19.4, 23.6, 26.8, 29.2, 31.7, 34.2, 36.8 and 39.6 min (read from the
  figure, not yet computed). MS1 TIC peaks near 30.7 and 35.8 min (about 0.9e6) and smaller ones
  near 33.3, 41.7 and 44.5 min have no FLD counterpart. Hypothesis to test: MS-only species
  (unlabelled, in-source fragments, adducts, contaminants).
* IgG Agilent run: MS not readable; the deposit lacks MSScan.bin, MSScan.xsd and MSTS.xml.
* Added `glyverity.align.estimate_lag` (cross-correlation; tested on a known 0.12 min shift)
  and lag reporting in the inspect script. pyproject licence moved to SPDX form.

**Next**
* Measure the FLD to MS lag on the fetuin run.
* Extract the m/z values under the MS-only peaks and check them against labelled and
  unlabelled glycan compositions, adducts and in-source fragments.
* Download one glycan raw file each from PXD077477 and PXD078897; check for a fluorescence channel.
* Browse GlycoPOST and MassIVE; add rows to the catalog.
