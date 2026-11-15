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
* label_series result: tag/2-AB area 0.350 to 0.461 over 8 compositions (pooled 0.436, 30.4%
  of tag plus 2-AB signal); tag elutes 0.48 to 1.15 min earlier. All six "free + dHex" hits are
  the M+4 isotope of the Na adduct of the co-eluting 2-AB glycan (1.9 to 3.0 ppm). The script's
  own summary line mixes in that interference and is not to be quoted.
* Full record: `docs/findings/2026-11-15_zenodo4705546_fetuin.md`.

* Added `glyverity.candidates` (targeted candidates, single-pass XIC extraction) and
  `glyverity.layers` with two registered evidence layers: `AdductIsotopeIsobar` (k Na/K/NH4
  exchanges plus n 13C, co-elution, stronger partner) and `LabelSideProductSeries` (ratio CV,
  same-sign RT shift, fluorescence fraction). Tests reproduce the fetuin relations.
  `scripts/run_layers.py` applies both over 350 compositions x 3 end groups.

* First `run_layers.py` run (fetuin): main 2-AB glycans unflagged; Hex8HexNAc6NeuAc1 2-AB
  correctly flagged as the Na adduct of the C3H9N-RA form of Hex6HexNAc5NeuAc3 (+0.8 ppm).
  Two flaws: (1) compositions with identical elemental formulas (dHex + NeuGc = Hex + NeuAc)
  were separate rows with identical signal; (2) the series test paired noise rows minutes
  apart (CV 1.3 to 1.5, mixed shift signs). Fixed: candidates merged by formula
  (350 compositions -> 290 formulas, 60 with isomeric alternatives listed); series test now
  pairs within 2 min, reference >= 1% of top, robust spread (MAD of log10 ratio) and
  same-sign fraction >= 0.8.
* Lead: two improbable 2-AB hits (Hex3HexNAc6NeuAc2/3NeuGc1) co-elute with the two main
  glycans at +24.011 Da; no adduct-isotope relation within 5 ppm.

* Re-run: C3H9N-RA vs 2-AB is a SERIES (n=24, ratio median 0.381, MAD log10 0.069, shift
  median -0.93 min, same sign 1.00, fluorescence fraction 0.43); free vs 2-AB n=0. Main 2-AB
  glycans unflagged. Hex8HexNAc6 2-AB flagged as Na adduct of the C3H9N-RA form of
  Hex6HexNAc5NeuAc2 (+0.9 ppm).
* Isotope shadow found: Hex6HexNAc6dHex1NeuAc1 and Hex5HexNAc6NeuAc1NeuGc1 (27.817 min, same
  area) differ by dHex + Hex - NeuGc = 1.0204 Da, 6 ppm from one 13C step (1.0034 Da), inside
  the 10 ppm XIC window. Isobar tolerance now defaults to the XIC tolerance.
  Hypothesis to test on GPST000788: part of the +1.020 Da series may be this composition
  coincidence rather than a label side-product.
* Summary now split by artifact type; added `recurring_deltas` (report-only).

* Re-run (10 ppm isobar): adduct/isotope 22.0%, side-product series 12.1%, unflagged 65.8%
  of candidate MS1 area. Series stable (n=20, ratio 0.381, MAD 0.077, shift -0.90 min).
  Recurring deltas: -1.0204 x4 (lighter weaker candidate's M+1 channel on the stronger one's
  monoisotopic peak; fixed by allowing n = -1 in isobar offsets); +15.9949 and +31.9898
  (NeuGc vs NeuAc, Hex vs dHex), -41.9841 (2-AB vs Hex), +25.0316 (-2 Hex + HexNAc + dHex):
  composition arithmetic between hypotheses, mostly improbable compositions.
* Conclusion: "unflagged" is not "real". A wide composition space yields coincidental
  matches; the decoy FDR step and a biosynthetic plausibility prior are needed next.

* Confirmed: -1.0204 pairs flagged; 7 main 2-AB glycans unflagged. Error found: the search
  filter (sialic acids <= HexNAc - 2) excluded Hex6HexNAc5NeuAc4, a real fluorescent glycan
  (39.68 min). Filter relaxed to <= HexNAc - 1 (476 compositions, 392 formulas).
* Delta mining surfaced -453.1482 Da x2 (Hex + NeuAc) between co-eluting 2-AB candidates and
  the two main glycans: in-source loss of a sialylated galactose. Added `InSourceFragment`
  layer (same end group, co-elution, weaker and lighter by a listed neutral loss).
* Summary now reports four classes: adduct/isotope, in-source fragment, side-product series,
  unflagged.

* Re-run (476 compositions): all 8 main 2-AB glycans present and unflagged, including
  Hex6HexNAc5NeuAc4. In-source fragments flagged: Hex5HexNAc5NeuAc1/2/3 2-AB as -Hex+NeuAc of
  A3S2/A3S3/A3S4, co-eluting. Series unchanged (n=20, ratio 0.381, MAD 0.060, -0.91 min).
* Ambiguous: Hex5HexNAc4NeuAc3 2-AB (36.33 min, FLD 1.89) flagged as Na adduct of
  Hex9HexNAc5 2-AB (hypotheses 5.2 ppm apart). "Stronger explains weaker" cannot decide; needs
  MS2 of m/z 1317.98 or a plausibility prior. Recorded as ambiguous, not artifact.
* Methods notes: (1) relation errors are between theoretical masses, so residue-loss relations
  are exact; evidence is co-elution plus the 10 ppm XIC window. (2) Area shares summed over
  candidates depend on the search space (adduct/isotope 22% -> 32% when the space grew), because
  one ion can fall in several candidates' windows. Not a census metric.

* Feature-first pipeline built: `glyverity.features` (MS1 scans -> OpenMS MSExperiment;
  MassTraceDetection, ElutionPeakDetection, FeatureFindingMetabo; mzML export; works with old
  and new pyOpenMS call styles), `glyverity.census` (species grouping over charge states and
  Na/K/NH4 adduct ions, composition assignment, one class per feature), and
  `scripts/census_run.py`. pyOpenMS 3.6.0 finding: with the noise threshold above all
  background, MassTraceDetection returns no traces (min_sample_rate); keep noise below
  background. Synthetic end-to-end test classifies every planted feature correctly.

* First census run (fetuin, pyOpenMS 3.6.0 on cluster): 668 MS1 scans -> 1139 features ->
  987 species. Top features correct (A3S3, A3S2, A2S2 2-AB, labelled glycan). Shares:
  labelled 23.1%, other end group 10.7%, ISF 0.7%, adduct/isotope 2.0%, adduct ion 11.6%,
  ambiguous 0.4%, unassigned 51.5%. Three faults: (1) A3S3 2+ and 3+ apexes 0.14 min apart
  split into two species (rt_tol 0.1 too tight; now 0.2); (2) [M+H+Na]2+ of A3S2 (1365.99,
  6th strongest) unassigned because its 3+ partner formed a species first and pass 2 only
  re-homed single features (now whole species are compared); (3) series test lacked FLD
  (now read, with --lag). Added an unassigned-feature report with deltas to co-eluting species.
  These shares are superseded by the next run.

* Census run 2: 1139 features -> 892 species. C3H9N-RA SERIES (n=8, ratio 0.403, spread 0.059,
  shift -0.89 min, fluorescence fraction 0.34). Shares: labelled 21.8%, other end group 2.9%,
  side-product 7.7%, ISF 0.9%, adduct/isotope 0.3%, adduct ion 25.9%, ambiguous 0.4%,
  unassigned 40.0%. Unassigned patterns (checked by arithmetic):
  * mixed Na/K salts beyond the charge: 1530.5130 2+ = A3S3 + Na + K - 2H; 1395.9567 2+ =
    A3S2 + 2Na + K - 3H (2.5 ppm);
  * envelopes split into 1+ features (pairs 0.5017 apart) and monoisotopic picked +4 (1502.5533);
  * oxonium ion 657.235 (NeuAc-Hex-HexNAc, 657.2349) co-eluting with main glycans;
  * 2880.047 Da at 41.6 min (2+ 1441.0292, 3+ 961.0240), no FLD = free A3S3 + NH3 - O
    (amino-alditol, 2880.0423, 1.6 ppm). Hypothesis; MS2 test: Y1 at 223.1288.
* Added: mixed-cation exchange combos (up to 4, any mix of Na/K/NH4) in species grouping and
  isobar layer; isotope-error relations (+1 to +4) in species grouping; `repair_charges`
  (merge split envelopes); oxonium-ion reclassification; NH3-RA end group; amino-alditol Y1 in
  `ms2_diagnostic.py`.

* Census run 3: 6 split envelopes merged; 630 species. C3H9N-RA SERIES (n=8, ratio 0.421,
  spread 0.033). Shares: labelled 21.7%, other end group 7.2%, side-product 7.7%, ISF 2.7%,
  adduct/isotope 0.1%, adduct ion 38.8%, ambiguous 0.2%, unassigned 21.6%.
* MS2 of the 2880 Da species: 961.02 (3+) shows Y1 223.1278 (amino-alditol, -4.9 ppm, 4.9%)
  and neither 342.166 nor 265.176; 1441.03 (2+) shows no Y1 of any kind. Controls clean
  (1000.70: 342 only; 975.04: 265 only). Amino-alditol supported, not proven.
* Fixes: split-envelope repair window 0.05 -> 0.3 min (real pairs 0.1 to 0.26 min apart);
  series max shift 2 -> 6 min (NH3-RA elutes about 4.8 min after 2-AB); in-source losses
  with 0 to 2 isotope errors added to species grouping and classified as in-source
  fragments (e.g. 1355.5004 2+ = A3S3 - NeuAc, M+1 picked).

**Next (Phase 2 design correction)**
* Feature-first census: detect each MS1 feature once (isotope envelope, charge, apex), then
  classify it (labelled glycan, side product, adduct/isotope, in-source fragment, unexplained).
  Shares of observed feature intensity are the census metric.
* Common input format: write Agilent runs (via rainbow) to mzML with pyOpenMS; Thermo runs via
  MSConvert; feature finding with OpenMS FeatureFinderMetabo or an envelope finder on centroids.
* Then decoys and q-values at the feature-assignment level.
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
* label_series result: tag/2-AB area 0.350 to 0.461 over 8 compositions (pooled 0.436, 30.4%
  of tag plus 2-AB signal); tag elutes 0.48 to 1.15 min earlier. All six "free + dHex" hits are
  the M+4 isotope of the Na adduct of the co-eluting 2-AB glycan (1.9 to 3.0 ppm). The script's
  own summary line mixes in that interference and is not to be quoted.
* Full record: `docs/findings/2026-11-15_zenodo4705546_fetuin.md`.

* Added `glyverity.candidates` (targeted candidates, single-pass XIC extraction) and
  `glyverity.layers` with two registered evidence layers: `AdductIsotopeIsobar` (k Na/K/NH4
  exchanges plus n 13C, co-elution, stronger partner) and `LabelSideProductSeries` (ratio CV,
  same-sign RT shift, fluorescence fraction). Tests reproduce the fetuin relations.
  `scripts/run_layers.py` applies both over 350 compositions x 3 end groups.

* First `run_layers.py` run (fetuin): main 2-AB glycans unflagged; Hex8HexNAc6NeuAc1 2-AB
  correctly flagged as the Na adduct of the C3H9N-RA form of Hex6HexNAc5NeuAc3 (+0.8 ppm).
  Two flaws: (1) compositions with identical elemental formulas (dHex + NeuGc = Hex + NeuAc)
  were separate rows with identical signal; (2) the series test paired noise rows minutes
  apart (CV 1.3 to 1.5, mixed shift signs). Fixed: candidates merged by formula
  (350 compositions -> 290 formulas, 60 with isomeric alternatives listed); series test now
  pairs within 2 min, reference >= 1% of top, robust spread (MAD of log10 ratio) and
  same-sign fraction >= 0.8.
* Lead: two improbable 2-AB hits (Hex3HexNAc6NeuAc2/3NeuGc1) co-elute with the two main
  glycans at +24.011 Da; no adduct-isotope relation within 5 ppm.

* Re-run: C3H9N-RA vs 2-AB is a SERIES (n=24, ratio median 0.381, MAD log10 0.069, shift
  median -0.93 min, same sign 1.00, fluorescence fraction 0.43); free vs 2-AB n=0. Main 2-AB
  glycans unflagged. Hex8HexNAc6 2-AB flagged as Na adduct of the C3H9N-RA form of
  Hex6HexNAc5NeuAc2 (+0.9 ppm).
* Isotope shadow found: Hex6HexNAc6dHex1NeuAc1 and Hex5HexNAc6NeuAc1NeuGc1 (27.817 min, same
  area) differ by dHex + Hex - NeuGc = 1.0204 Da, 6 ppm from one 13C step (1.0034 Da), inside
  the 10 ppm XIC window. Isobar tolerance now defaults to the XIC tolerance.
  Hypothesis to test on GPST000788: part of the +1.020 Da series may be this composition
  coincidence rather than a label side-product.
* Summary now split by artifact type; added `recurring_deltas` (report-only).

* Re-run (10 ppm isobar): adduct/isotope 22.0%, side-product series 12.1%, unflagged 65.8%
  of candidate MS1 area. Series stable (n=20, ratio 0.381, MAD 0.077, shift -0.90 min).
  Recurring deltas: -1.0204 x4 (lighter weaker candidate's M+1 channel on the stronger one's
  monoisotopic peak; fixed by allowing n = -1 in isobar offsets); +15.9949 and +31.9898
  (NeuGc vs NeuAc, Hex vs dHex), -41.9841 (2-AB vs Hex), +25.0316 (-2 Hex + HexNAc + dHex):
  composition arithmetic between hypotheses, mostly improbable compositions.
* Conclusion: "unflagged" is not "real". A wide composition space yields coincidental
  matches; the decoy FDR step and a biosynthetic plausibility prior are needed next.

* Confirmed: -1.0204 pairs flagged; 7 main 2-AB glycans unflagged. Error found: the search
  filter (sialic acids <= HexNAc - 2) excluded Hex6HexNAc5NeuAc4, a real fluorescent glycan
  (39.68 min). Filter relaxed to <= HexNAc - 1 (476 compositions, 392 formulas).
* Delta mining surfaced -453.1482 Da x2 (Hex + NeuAc) between co-eluting 2-AB candidates and
  the two main glycans: in-source loss of a sialylated galactose. Added `InSourceFragment`
  layer (same end group, co-elution, weaker and lighter by a listed neutral loss).
* Summary now reports four classes: adduct/isotope, in-source fragment, side-product series,
  unflagged.

* Re-run (476 compositions): all 8 main 2-AB glycans present and unflagged, including
  Hex6HexNAc5NeuAc4. In-source fragments flagged: Hex5HexNAc5NeuAc1/2/3 2-AB as -Hex+NeuAc of
  A3S2/A3S3/A3S4, co-eluting. Series unchanged (n=20, ratio 0.381, MAD 0.060, -0.91 min).
* Ambiguous: Hex5HexNAc4NeuAc3 2-AB (36.33 min, FLD 1.89) flagged as Na adduct of
  Hex9HexNAc5 2-AB (hypotheses 5.2 ppm apart). "Stronger explains weaker" cannot decide; needs
  MS2 of m/z 1317.98 or a plausibility prior. Recorded as ambiguous, not artifact.
* Methods notes: (1) relation errors are between theoretical masses, so residue-loss relations
  are exact; evidence is co-elution plus the 10 ppm XIC window. (2) Area shares summed over
  candidates depend on the search space (adduct/isotope 22% -> 32% when the space grew), because
  one ion can fall in several candidates' windows. Not a census metric.

* Feature-first pipeline built: `glyverity.features` (MS1 scans -> OpenMS MSExperiment;
  MassTraceDetection, ElutionPeakDetection, FeatureFindingMetabo; mzML export; works with old
  and new pyOpenMS call styles), `glyverity.census` (species grouping over charge states and
  Na/K/NH4 adduct ions, composition assignment, one class per feature), and
  `scripts/census_run.py`. pyOpenMS 3.6.0 finding: with the noise threshold above all
  background, MassTraceDetection returns no traces (min_sample_rate); keep noise below
  background. Synthetic end-to-end test classifies every planted feature correctly.

* First census run (fetuin, pyOpenMS 3.6.0 on cluster): 668 MS1 scans -> 1139 features ->
  987 species. Top features correct (A3S3, A3S2, A2S2 2-AB, labelled glycan). Shares:
  labelled 23.1%, other end group 10.7%, ISF 0.7%, adduct/isotope 2.0%, adduct ion 11.6%,
  ambiguous 0.4%, unassigned 51.5%. Three faults: (1) A3S3 2+ and 3+ apexes 0.14 min apart
  split into two species (rt_tol 0.1 too tight; now 0.2); (2) [M+H+Na]2+ of A3S2 (1365.99,
  6th strongest) unassigned because its 3+ partner formed a species first and pass 2 only
  re-homed single features (now whole species are compared); (3) series test lacked FLD
  (now read, with --lag). Added an unassigned-feature report with deltas to co-eluting species.
  These shares are superseded by the next run.

* Census run 2: 1139 features -> 892 species. C3H9N-RA SERIES (n=8, ratio 0.403, spread 0.059,
  shift -0.89 min, fluorescence fraction 0.34). Shares: labelled 21.8%, other end group 2.9%,
  side-product 7.7%, ISF 0.9%, adduct/isotope 0.3%, adduct ion 25.9%, ambiguous 0.4%,
  unassigned 40.0%. Unassigned patterns (checked by arithmetic):
  * mixed Na/K salts beyond the charge: 1530.5130 2+ = A3S3 + Na + K - 2H; 1395.9567 2+ =
    A3S2 + 2Na + K - 3H (2.5 ppm);
  * envelopes split into 1+ features (pairs 0.5017 apart) and monoisotopic picked +4 (1502.5533);
  * oxonium ion 657.235 (NeuAc-Hex-HexNAc, 657.2349) co-eluting with main glycans;
  * 2880.047 Da at 41.6 min (2+ 1441.0292, 3+ 961.0240), no FLD = free A3S3 + NH3 - O
    (amino-alditol, 2880.0423, 1.6 ppm). Hypothesis; MS2 test: Y1 at 223.1288.
* Added: mixed-cation exchange combos (up to 4, any mix of Na/K/NH4) in species grouping and
  isobar layer; isotope-error relations (+1 to +4) in species grouping; `repair_charges`
  (merge split envelopes); oxonium-ion reclassification; NH3-RA end group; amino-alditol Y1 in
  `ms2_diagnostic.py`.

* Census run 3: 6 split envelopes merged; 630 species. C3H9N-RA SERIES (n=8, ratio 0.421,
  spread 0.033). Shares: labelled 21.7%, other end group 7.2%, side-product 7.7%, ISF 2.7%,
  adduct/isotope 0.1%, adduct ion 38.8%, ambiguous 0.2%, unassigned 21.6%.
* MS2 of the 2880 Da species: 961.02 (3+) shows Y1 223.1278 (amino-alditol, -4.9 ppm, 4.9%)
  and neither 342.166 nor 265.176; 1441.03 (2+) shows no Y1 of any kind. Controls clean
  (1000.70: 342 only; 975.04: 265 only). Amino-alditol supported, not proven.
* Fixes: split-envelope repair window 0.05 -> 0.3 min (real pairs 0.1 to 0.26 min apart);
  series max shift 2 -> 6 min (NH3-RA elutes about 4.8 min after 2-AB); in-source losses
  with 0 to 2 isotope errors added to species grouping and classified as in-source
  fragments (e.g. 1355.5004 2+ = A3S3 - NeuAc, M+1 picked).

**Next (Phase 2 design correction)**
* Feature-first census: detect each MS1 feature once (isotope envelope, charge, apex), then
  classify it (labelled glycan, side product, adduct/isotope, in-source fragment, unexplained).
  Shares of observed feature intensity are the census metric.
* Common input format: write Agilent runs (via rainbow) to mzML with pyOpenMS; Thermo runs via
  MSConvert; feature finding with OpenMS FeatureFinderMetabo or an envelope finder on centroids.
* Then decoys and q-values at the feature-assignment level.
* Extract the m/z values under the MS-only peaks and check them against labelled and
  unlabelled glycan compositions, adducts and in-source fragments.
* Download one glycan raw file each from PXD077477 and PXD078897; check for a fluorescence channel.
* Browse GlycoPOST and MassIVE; add rows to the catalog.

## 2026-11-15 (late): series fluorescence check, lone 1+ features

**Goal:** fix the two faults found in census run 4.

* Run 4: 19 split envelopes merged; 585 species. C3H9N-RA SERIES (n=8, ratio 0.411, spread
  0.024, shift -0.89 min, FF 0.34). NH3-RA n=8, ratio 0.188, spread 0.068, shift +4.57 min,
  FF 3.11 -> rejected. Shares: labelled 21.5%, other end group 7.0%, side-product 7.7%,
  ISF 3.9%, adduct/isotope 0.1%, adduct ion 41.4%, ambiguous 0.2%, unassigned 18.3%.
* NH3-RA rejection traced to the test, not the data: forms eluting 4.6 min after their 2-AB
  partner sit on other 2-AB glycans' FLD peaks, and a small ratio makes the expected FLD tiny.
  Fix: pairs whose test apex lies within 0.3 min of any 2-AB apex are left out of the
  fluorescence fraction; with no clear pair the result is "untested" and the decision rests on
  ratio and shift consistency. This also lets the series test run on datasets with no FLD.
* Free-isotope alternative for NH3-RA (+1.0034 vs +1.0316 Da) ruled out: the free
  monoisotopic species is absent (free series n=0).
* Split-envelope repair now uses intensity-weighted RT (1194.4162, A2S2 disodium salt, had
  taken the first peak's time, 0.21 min from its parent).
* Co-elution tolerance 0.2 -> 0.25 min; oxonium window 0.3 min (657.235 is 0.205 min from A3S4).
* A lone 1+ feature that matches no relation is retried as 2+ and 3+ (1530.5130 = A3S3 Na+K
  salt as 2+). Labelled "(read as z+)" so these can be counted separately.
* Tests: 61 passed (3 new).

**Next:** re-run census on the cluster; then decoys and q-values.

## 2026-11-15 (late, 2): run 5 and a decoy check on grouping

* Run 5: 428 species (585 in run 4). NH3-RA n=8, ratio 0.180, spread 0.034, shift +4.57 min,
  fluorescence fraction 2.46 over 6 clear pairs -> still rejected. C3H9N-RA SERIES (ratio
  0.424, spread 0.033, FF 0.19 over 5 clear pairs). Shares: labelled 21.4%, other end group
  6.9%, side-product 7.7%, ISF 4.7%, adduct/isotope 0.1%, adduct ion 48.3%, ambiguous 0.1%,
  unassigned 10.7%.
* Concern: the unassigned drop (18.3 -> 10.7%) went almost entirely into "adduct ion"
  (41.4 -> 48.3%). Pass 2 tests about 200 offsets, three charge readings for lone 1+ features
  and several parents per window at 10 ppm, so chance matches are possible. Added a decoy
  estimate: every relation offset is shifted by a non-chemical amount (0.23, 0.41, 0.59,
  0.77 Da) and the absorptions counted. Not trusted until that is run.
* A formula search for the recurring +217.878 Da delta (seen on two 2-AB parents) returns many
  salt/acid combinations within 2 mDa. Not identified; a combinatorial search alone cannot
  decide it. Needs the isotope pattern (37Cl, 41K) of those features.
* NH3-RA: the fluorescence test is not being adjusted again to change its answer. Added
  `--series-detail` to print each pair's FLD values for inspection.

## 2026-11-15 (late, 3): decoy result, dimers, FLD tail correction

* Decoy grouping (offset shifts 0.23/0.41/0.59/0.77 Da): 38.0 of 665 absorbed features
  (5.7%), 3.1% of absorbed intensity, are chance matches; re-read 1+ features 10.5 of 117
  (9.0%). Absorbed intensity is 51.0% of total, of which about 1.6 points are chance. The
  grouping step holds; the re-read rule is the weakest part and stays reported separately.
* 2008.3909 (3+, 36.837 min) = [2M+Na+2H]3+ of Hex6HexNAc5NeuAc3 2-AB picked at M+2,
  0.3 ppm, co-eluting with the monomer (36.790 min). Added non-covalent dimer relations
  (2M plus the same exchanges and isotope errors, tried after all monomer relations) and a
  "dimer ion" class. Covered by the decoy shift.
* NH3-RA pair table: high fluorescence fractions sit on tails of strong 2-AB peaks
  (31.269 min is 0.46 min after A3S2 2-AB, FLD 6.36; reads 1.76 vs 0.38 expected). The two
  pairs about 2 min or more from any paired 2-AB apex (41.586, 44.492 min) read 7% and 13% of
  expected. The fixed 0.3 min clearance is narrower than the detector's tailing.
* Replaced the fixed clearance with a tail correction: median normalised FLD peak shape from
  isolated 2-AB peaks (each side separately), neighbours' predicted tails subtracted at the
  test apex; a pair is used only when the predicted tail is no larger than the FLD expected
  from the test form. Rule set from 2-AB peaks only; applies to both series. Tests cover a
  non-fluorescent form on a tail (accepted) and a fluorescent one (rejected).
* Tests: 65 passed.
