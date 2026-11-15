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

* Lag re-run: +0.037 min, Pearson r 0.903 (window 15 to 50 min).
* Hand check of windowed spectra (2+ and 3+ ions agree within 2 ppm in every case):
  * 31.7 min (FLD and MS): Hex6HexNAc5NeuAc2 2-AB, +2.9 ppm. 36.8 min: Hex6HexNAc5NeuAc3 2-AB, +1.6 ppm.
  * 30.7 min (MS only): neutral 2630.9995 Da. 35.8 min (MS only): neutral 2922.0940 Da.
    They differ by one NeuAc (291.0945 Da). Each is 76.99 Da lighter than the 2-AB form of the
    glycan eluting about 1 min later. No composition (Hex, HexNAc, dHex, NeuAc, NeuGc; free,
    2-AB, reduced; acetyl, sulfate, phosphate, water loss) fits within 5 ppm.
  * Hypothesis: free glycan + 43.0786 Da (reductive amination with a C3H9N amine) after a
    +2 ppm calibration offset. Test: Y1 ion at m/z 265.1758 in MS2.
* First-pass matcher lessons for Phase 2: rank by more than ppm; constrain charge by isotope
  spacing; do not assume the monoisotopic peak is the most intense; the M+1 of a 2+ ion
  (1355.5039) matched Man5-2AB 1+ at -4 ppm.
* Added `scripts/ms2_diagnostic.py`.
* MS2 result: the four MS-only precursors (1316.51, 878.01, 1462.05, 975.04) show glycan oxonium
  ions (204.087, 366.140, 292.103, 274.092) and Y1 at 265.17 (GlcNAc + 43.08), and no 2-AB Y1
  (342.166). The 2-AB controls 1355.00, 903.67 and 1000.70 show 342.165 and no 265.17
  (1500.55 had 7 scans and showed neither). Conclusion: triantennary di- and trisialylated
  glycans carrying a non-fluorescent reducing-end group of net +C3H9N-O (+43.0786 Da).
  Predicted m/z for that group: 1316.5042, 878.0052, 1462.0519, 975.0370 (observed +0.9 to
  +3.2 ppm, in line with the +2 ppm offset seen on the 2-AB glycans).
* Possible origin (hypothesis): acetone was used to remove unreacted 2-AB after labelling with
  cyanoborohydride (Zenodo record description); reductive amination with acetone can add an
  isopropyl group. Not testable without new experiments; report masses and fragments only.
* Added `glyverity.xic` and `scripts/label_series.py` (MS1 area per reducing-end group).

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

* Fetuin run: MS lags FLD by +0.037 min (2.2 s). The printed correlation (1.021) was wrong:
  the estimator centred on the median, not the mean. Fixed to a Pearson correlation; re-run due.
* IgG run confirmed: no MSScan files in AcqData.
* Added `glyverity.annotate` and `scripts/peak_spectra.py`: summed MS1 spectra per window,
  isotope flags, and composition matching (free and 2-AB; H+, Na+, K+, NH4+; z 1 to 3).

* Lag re-run: +0.037 min, Pearson r 0.903 (window 15 to 50 min).
* Hand check of windowed spectra (2+ and 3+ ions agree within 2 ppm in every case):
  * 31.7 min (FLD and MS): Hex6HexNAc5NeuAc2 2-AB, +2.9 ppm. 36.8 min: Hex6HexNAc5NeuAc3 2-AB, +1.6 ppm.
  * 30.7 min (MS only): neutral 2630.9995 Da. 35.8 min (MS only): neutral 2922.0940 Da.
    They differ by one NeuAc (291.0945 Da). Each is 76.99 Da lighter than the 2-AB form of the
    glycan eluting about 1 min later. No composition (Hex, HexNAc, dHex, NeuAc, NeuGc; free,
    2-AB, reduced; acetyl, sulfate, phosphate, water loss) fits within 5 ppm.
  * Hypothesis: free glycan + 43.0786 Da (reductive amination with a C3H9N amine) after a
    +2 ppm calibration offset. Test: Y1 ion at m/z 265.1758 in MS2.
* First-pass matcher lessons for Phase 2: rank by more than ppm; constrain charge by isotope
  spacing; do not assume the monoisotopic peak is the most intense; the M+1 of a 2+ ion
  (1355.5039) matched Man5-2AB 1+ at -4 ppm.
* Added `scripts/ms2_diagnostic.py`.
* MS2 result: the four MS-only precursors (1316.51, 878.01, 1462.05, 975.04) show glycan oxonium
  ions (204.087, 366.140, 292.103, 274.092) and Y1 at 265.17 (GlcNAc + 43.08), and no 2-AB Y1
  (342.166). The 2-AB controls 1355.00, 903.67 and 1000.70 show 342.165 and no 265.17
  (1500.55 had 7 scans and showed neither). Conclusion: triantennary di- and trisialylated
  glycans carrying a non-fluorescent reducing-end group of net +C3H9N-O (+43.0786 Da).
  Predicted m/z for that group: 1316.5042, 878.0052, 1462.0519, 975.0370 (observed +0.9 to
  +3.2 ppm, in line with the +2 ppm offset seen on the 2-AB glycans).
* Possible origin (hypothesis): acetone was used to remove unreacted 2-AB after labelling with
  cyanoborohydride (Zenodo record description); reductive amination with acetone can add an
  isopropyl group. Not testable without new experiments; report masses and fragments only.
* Added `glyverity.xic` and `scripts/label_series.py` (MS1 area per reducing-end group).

**Next**
* Extract the m/z values under the MS-only peaks and check them against labelled and
  unlabelled glycan compositions, adducts and in-source fragments.
* Download one glycan raw file each from PXD077477 and PXD078897; check for a fluorescence channel.
* Browse GlycoPOST and MassIVE; add rows to the catalog.
