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

## 2026-11-15 (late, 4): run 6 exposed two faults

* Run 6: 378 species; chance absorptions rose to 10.3% of absorbed features (5.7% before
  dimers), decoy mean 73.5 vs 38.0. About 35 extra chance absorptions against 115 target
  dimer features: the dimer relation as built had a chance rate near 30%. Cause: about 175
  more offsets at twice the parent mass, and dimers allowed for features without a measured
  charge ([2M+2zH]2z+ sits at the m/z of [M+zH]z+, so only isotope spacing tells them apart).
  Fix: dimers only for species whose features all have 2 or more isotopes, never for re-read
  1+ features. Decoy counts now reported per relation family.
* Tail correction failed a sanity check: predicted tails exceeded observed FLD (6.14 vs 0.51
  at 30.713 min; 15 vs 4.65 at 36.542 min). Cause: neighbours were taken from MS label
  candidates, so co-eluting candidates counted one FLD peak several times, and weak
  candidates carried a neighbour's FLD as their own height. Fix: FLD peaks detected in the
  trace itself (each counted once, 1% prominence); a pair is dropped as a model failure if
  the predicted tail exceeds the observed FLD by more than 10%.
* Run-6 verdicts not used: C3H9N-RA "untested" and NH3-RA "SERIES" both came from the faulty
  correction. The raw NH3-RA values at the two pairs far from strong peaks stand on their
  own: 41.586 min FLD 0.092 vs 1.29 expected; 44.492 min 0.032 vs 0.236.
* Tests: 69 passed.

## 2026-11-15 (late, 5): per-family tolerances from decoys

* Run 7 (dimer fix, trace-based FLD peaks): 390 species. Decoy chance share by relation
  family at 10 ppm: adduct 497 / 30.5 (6.1%), isotope 47 / 0.5 (1.1%), in-source fragment
  56 / 24.8 (44.2%), dimer 103 / 10.2 (10.0%). In-source fragment grouping at 10 ppm is close
  to half chance.
* FLD tail shape could not be measured: no FLD peak in this dense trace is isolated within
  2.5 min on either side, so the fixed 0.3 min clearance was used (FLD_tail nan). NH3-RA
  verdict returns to run 5 (FF 2.47, 6 pairs). Decision: stop changing the fluorescence test
  on this run. NH3-RA stays classed as "other end group" (conservative). Status: MS2 Y1 223.128
  and ratio/shift consistency support it; fluorescence inconclusive (the two pairs far from
  strong peaks read 7% and 13% of expected; the rest sit on tails).
* Added decoy-calibrated tolerances: each absorption records (family, ppm error) in target
  and decoy runs; per family, the largest ppm cutoff with decoy/target <= 5% is used for the
  final grouping, and the chance share is re-estimated after. Lone 1+ re-reads are their own
  family. This is target-decoy FDR control on the cleaning step, ahead of composition q-values.
* Tests: 71 passed.

## 2026-11-15 (late, 6): relation score uses mass and retention time

* Run 8 (ppm-only family tolerances): adduct 5.81 ppm, isotope 7.53, in-source fragment 0
  (off), dimer 7.32, re-read 1.25. Final chance share 7.5% (target 5%); unassigned 13.2%.
* Real ions were cut: 1524.0330 as 2+ = A3S3 2-AB 2Na M+3 (0.6 ppm); 864.3259 as 3+ = A3S2
  NH3-RA M+1 (2.2 ppm, same RT as 1002.7067); 1355.5004 (2+) and 904.0043 (3+) = A3S3 - NeuAc
  M+1 (2.4 and 0.2 ppm). Mass error alone cannot separate them from chance among about 200
  relations.
* Score is now max(ppm/10, |dRT|/0.25 min). Source-formed ions co-elute exactly with their
  parent; chance matches spread over the time window. Cutoffs per family are set on this
  score by the decoys as before.
* Open: cutoffs are chosen on the first grouping pass, so the final chance share can exceed
  the target (7.5% vs 5% in run 8). Report the achieved value; iterate if it stays high.
* Tests: 71 passed.

## 2026-11-15 (late, 7): run 9, chance share under target

* Run 9 (score = max(ppm/10, dRT/0.25)): cutoffs adduct 0.822 (8.2 ppm, 0.206 min), isotope
  0.869, dimer 0.755, re-read 0.600 (6.0 ppm, 0.150 min), in-source fragment 0 (off).
  Final: 470 species, 593 absorptions (48.1% of intensity), decoy mean 27.5, estimated
  chance share 4.6% (target 5%). Shares: labelled 21.7%, other end group 7.2%, side-product
  7.7%, ISF 3.1% (oxonium and candidate-level ISF only), adduct/isotope 0.2%, adduct ion
  45.5%, dimer 2.6%, ambiguous 0.2%, unassigned 11.9%.
* Series verdicts unchanged: C3H9N-RA SERIES (ratio 0.416, spread 0.035, FF 0.19 over 5
  pairs); NH3-RA not a series on fluorescence (FF 2.57), classed as other end group.

**Next session**
* In-source fragments are still off: the family pools 10 losses x 3 isotope errors, so
  decoy matches are dense. Real ones are in the unassigned list (1355.5004 2+ and 904.0043
  3+ = A3S3 2-AB - NeuAc M+1; dRT 0.095 and 0.047 min). Try sub-families (sialic acid losses
  vs Hex/HexNAc/dHex/H2O, which also match biosynthetic neighbours) with their own cutoffs.
* 1524.0330 (1+, 36.742 min) = A3S3 2-AB 2Na M+3 as 2+ (0.6 ppm, dRT 0.048 min) yet still
  unassigned. Hypothesis: not a lone feature (grouped in pass 1 with a partner), so it is never
  re-read. Checked: feature 814, species 38, form H, class unassigned; no 762.5 partner. Still
  to check: other members of species 38 (e.g. a 3+ at m/z 508.7), and the RT of the A3S3
  2-AB parent species, since the re-read cutoff includes dRT (0.150 min).
* Recurring +217.88 / +212.87 / +202.92 / +195.90 deltas on 2-AB parents: inspect isotope
  patterns (37Cl, 41K) before naming anything.
* Then composition decoys and q-values (go/no-go condition 3).

## 2026-11-15 (late, 8): two grouping bugs found from the 1524.0330 check

* Species 38 = 1524.0330 (1+) plus feature 891 (1590.9753, 2+) labelled
  "dimer:2M+3K+1NH4 M+3" of that 1+ ion, a likely chance match on an unidentified parent.
  Because 891 was absorbed first, 1524.0330 no longer counted as lone and was never re-read
  as the 2Na M+3 salt of A3S3 2-AB (0.6 ppm, dRT 0.05 min).
* Second bug in the same table: when a species was absorbed, every member took the new
  relation as its form, so ions it had absorbed lost their own relation (e.g. 2015.3848 1+
  shown as an A3S3 dimer, which is impossible at that mass).
* Fixes: "lone" counts pass-1 members only; dimers need measured charges on both the
  candidate and the parent's own features; absorbed members keep their own relation as
  "<own> via <species relation>". Test added that fails on the old code.
* Tests: 72 passed. Re-run needed: dimer counts and the family cutoffs will change.

## 2026-11-15 (late, 9): run 10 after the grouping fixes

* 382 species at 10 ppm. Families (target / decoy / chance): adduct 359 / 21.8 / 6.1%,
  isotope 39 / 0.2 / 0.6%, isf 52 / 23.5 / 45.2% (off), dimer 101 / 8.0 / 7.9%, re-read
  160 / 12.8 / 8.0%. Cutoffs: adduct 0.822, isotope 0.869, dimer 0.755, re-read 0.619.
  Final: 460 species, 603 absorptions (48.9% of intensity), chance share 4.6%.
* 1524.0330 now absorbed (re-reads rose from 101 to 160 once absorbed members stopped
  blocking the lone check). Dimer class barely changed (100 features, 2.7%): check its
  members by charge next session before trusting it.
* Shares: labelled 21.7%, other end group 7.2%, side-product 7.7%, ISF 2.9%, adduct/isotope
  0.4%, adduct ion 46.1%, dimer 2.7%, ambiguous 0.2%, unassigned 11.1%.
* Next: sialic-acid loss sub-family for in-source fragments (1355.5004, 904.0043 still
  unassigned); dimer member check; +217.88 delta isotope patterns; then composition decoys
  and q-values.

## 2026-11-20: dimer check, in-source fragment sub-families

* Dimer ions from run 10 by charge: z=1 2, z=2 7, z=3 66, z=4 25. Of the first 39 listed,
  33 are dimers of the two most intense species (species 0 = A3S3 2-AB, 2-AB at 36.6 to
  37.0 min; species 1 = A3S2 2-AB at 31.6 to 31.9 min) and 3 are dimers of unassigned
  species with no composition. Dimer formation scales with the square of concentration, so
  concentration in the two strongest glycans is the expected pattern; chance matches would
  scatter across species. 3+ dimers cannot coincide with any monomer charge state; 4+ dimers
  carrying an odd number of Na/K cannot either. Class kept.
* Chained forms ("... via 1Na") reach 5 to 6 exchanges in total; covered by the decoy shift,
  but report the count of chained dimers separately in the paper.
* In-source fragment family split: sialic acid losses (NeuAc, NeuGc, 2NeuAc, Hex+NeuAc,
  NeuAc+H2O) vs other losses (Hex, HexNAc, dHex, Hex+HexNAc, H2O), each with its own decoy
  cutoff. Tests: 73 passed.

## 2026-11-20 (2): decoy design fault, fixed

* Run 11 (ISF split): isf_sia 32 target / 22.5 decoy (70.3%), isf_other 20 / 1.0 (5.0%,
  cutoff 0.915). Final chance 4.5%, unassigned 10.9%; 1355.5004 and 904.0043 still unassigned.
* 22.5 decoys from 15 sialic-loss offsets, against 21.8 from about 200 adduct offsets, is not
  a chance rate. Cause: every family was shifted at once, so no feature was absorbed as an
  adduct in decoy runs; features that real adducts claim in the target run stayed free and
  were matched by families tested later. Late families (ISF, re-read) had inflated decoy
  counts; earlier per-family cutoffs are biased low for them. The pooled estimate was
  conservative.
* Fix: one family shifted at a time, the others real, so the decoy sees the target's
  candidate pool. Final chance share = sum of per-family decoy counts over target absorptions.
  Grouping now runs 1 + 24 + 1 + 24 times. Tests: 74 passed.

## 2026-11-20 (3): run 12, per-family decoys

* Families (target / decoy / chance -> cutoff): adduct 359 / 11.0 / 3.1% -> 0.998
  (10.0 ppm, 0.250 min); isotope 39 / 0.2 / 0.6% -> 0.869; isf_sia 32 / 9.0 / 28.1% -> 0.217
  (2.2 ppm, 0.054 min); isf_other 20 / 0.2 / 1.2% -> 0.915; dimer 101 / 6.5 / 6.4% -> 0.857;
  re-read 160 / 10.8 / 6.7% -> 0.619 (6.2 ppm, 0.155 min).
* With one family shifted at a time, the adduct decoy count halved (21.8 -> 11.0) and isf_sia
  fell from 22.5 to 9.0, confirming the pool bias.
* Final: 424 species, 639 absorptions (51.2% of intensity), per-family decoy sum 23.5,
  estimated chance share 3.7%. Shares: labelled 21.4%, other end group 6.9%, side-product
  7.7%, ISF 4.6%, adduct/isotope 0.1%, adduct ion 46.1%, dimer 2.9%, ambiguous 0.1%,
  unassigned 10.1%. 1355.5004 and 904.0043 now absorbed as A3S3 - NeuAc.
* Grouping step settled for this run: every tolerance set by decoys, chance share 3.7%.
* Next: composition decoys and q-values (go/no-go condition 3); +217.88 delta isotope
  patterns; second dataset.

## 2026-11-20 (4): composition q-values

* Literature check on decoys for MS1 formula assignment: JUMPm (Wang et al., J. Proteome Res.
  2018, 17, 2328-2334, doi:10.1021/acs.jproteome.8b00019) adds odd numbers of H (octet-rule
  violation; +1H, and +3H to +9H tested); pySM (Palmer et al., Nat. Methods 2017,
  doi:10.1038/nmeth.4072) uses implausible adducts. Odd-H decoys keep the mass defect of CHNO
  molecules, which fractional shifts do not; glycans occupy a narrow mass-defect band, so a
  fractional-shift decoy would under-sample chance matches and make q-values optimistic.
  (Same question applies to the grouping decoys; to check.)
* New module fdr.py: decoy libraries target +/- 3, 5, 7, 9 H (decoys within the search
  tolerance of any target removed), separate-search FDR, monotone q-values, and an
  exchangeability check (target vs decoy match counts at 15-30 ppm, where true matches are
  rare; ratio 1 expected).
* census_run.py: q-value per species after grouping; counts and intensity shares at q <= 0.05
  and 0.01; class shares split at q = 0.05; species_target_ppm and species_q columns in the TSV.
* Tests: 77 passed (simulated data: q <= 0.05 gives observed FDR < 10% and > 80% recall;
  exchangeability ratio within 0.8 to 1.25).

## 2026-11-20 (5): first composition q-values on fetuin

* 1568 target masses; 8 odd-H libraries, 1400 of 12544 decoys removed (within 10 ppm of a
  target). Exchangeability at 15-30 ppm: target 18, decoy mean 10.0, ratio 1.80.
* 97 of 424 species match a target (85.1% of feature intensity); q <= 0.05: 73 species
  (82.9%); q <= 0.01: 3 species (0.3%).
* Class split at q = 0.05 (species intensity incl. absorbed ions): labelled glycan 50.5% vs
  0.2%; other end group 13.7 vs 1.9; side-product 17.9 vs 0.0; ISF 0.8 vs 2.3; unassigned
  0 vs 12.6.
* Faults and reading:
  - Removing decoys thinned the libraries by 11% without correction (optimistic). Fixed:
    decoy counts scaled by target size / decoy library size.
  - Ratio 1.80 on 18 vs 10 counts. Either the decoys under-sample chance matches, or true
    glycans with poor mass accuracy fall in 15-30 ppm. Added a listing of the band species
    to tell these apart.
  - q <= 0.01 is nearly unreachable: with about 0.5 target masses per Da and ~330 non-glycan
    species, chance matches within 1-2 ppm occur in every decoy library. Accurate mass alone
    has an FDR floor of a few percent at this search-space size. Other evidence (isotope fit,
    retention, series, MS2) is needed for 1%. Paper point, not a bug.
* Tests: 78 passed.

## 2026-11-20 (6): scaled decoys, band species, recalibration

* With decoy counts scaled for removed decoys: exchangeability 18 vs 11.4 (ratio 1.58);
  q <= 0.05: 68 species, 70.0% of intensity; labelled glycan class 37.6% at q <= 0.05 vs
  13.0% above.
* Band (15-30 ppm) species are not fetuin glycans with poor accuracy: assignments such as
  Hex8HexNAc6NeuGc1(C3H9N-RA) and Hex9HexNAc6dHex1(C3H9N-RA); masses near 2883 to 2884 Da
  recur at 32.7, 44.3 and 44.6 min (1442.5236 2+ is among the strongest unassigned). A
  structured unknown family repeated at several retention times raises variance of the band
  count; 18 vs 11.4 does not establish a decoy bias on one run. Re-test across datasets.
* Losing 13% of intensity from labelled-glycan species above q = 0.05 points to the score:
  ppm error alone ranks an intense glycan at 6 ppm below a weak unknown at 1 ppm. First fix:
  mass recalibration (Theil-Sen, linear in mass) on anchor species seen in 2+ charge states
  within 10 ppm; composition assignment and q-values on recalibrated masses; raw q kept for
  comparison; strongest matched species printed with raw and recalibrated errors.
* Tests: 79 passed.

## 2026-11-20 (7): recalibration result

* 21 anchors; offset +1.68 ppm at 2 kDa, slope -0.05 ppm/kDa; anchor |error| median
  1.68 -> 0.32 ppm. Major fetuin glycans now within 1 ppm.
* q <= 0.05: 56 species, 81.3% of intensity (raw: 68, 70.0%); q <= 0.01: 7 species, 14.5%.
  Labelled-glycan species 50.0% at q <= 0.05 vs 0.6% above (13.0% before recalibration).
* q of the top glycans clusters at 0.012: one scaled decoy hit (0.14) over about 12 targets
  at sub-ppm errors. q <= 0.01 is limited by decoy count granularity on one run.
* Exchangeability 17 vs 12.0 (1.42). go_no_go.md updated: condition 3 design done,
  exchangeability not yet shown.

## 2026-11-20 (8): FLD validation built

* New module validate.py: FLD windows valley to valley around each species apex (FLD clock =
  MS time - lag), straight-line baseline, area above it; species in one FLD peak share it.
  Agreement measures: total variation distance between MS and FLD share vectors, and slope of
  log10(MS share / FLD share) against sialic acid count.
* census_run.py: for 2-AB labelled-glycan species with q <= 0.05 and one live match, compares
  three MS measures with FLD: strongest single ion, protonated charge states, all ions
  assigned to the species (salts, isotope picks, fragments, dimers). fld_validation.tsv.
* Prediction: if cleaning is right, "all assigned ions" has the smallest TV distance and the
  flattest sialic acid slope. Known limits: MS response still differs by structure; dimer ions
  are counted once though they carry two molecules; weak side-product fluorescence can sit in
  a window.
* Tests: 82 passed.

## 2026-11-20 (9): first FLD validation result, controls added

* 21 confident 2-AB glycans (q <= 0.05) in 13 FLD peaks. TV distance to FLD shares:
  strongest ion 0.140, protonated ions 0.092, all assigned ions 0.055. Slope of
  log10(MS/FLD) per sialic acid: -0.026, +0.006, +0.052.
* Supports cleaning on the distance measure (all assigned ions 40% closer than protonated).
  The slope prediction failed: protonated ions were already flat and all assigned ions
  over-weight sialylated glycans (+0.052 per sialic acid, about 13%). Sialic acid salts add
  more ions to highly sialylated species; whether this is over-absorption or a real
  ionization difference is open.
* Shared FLD peaks to inspect: 36.73 min (A3S3 with a Hex5HexNAc4NeuAc2 assignment) and
  23.60 min (A2S2 with Hex8HexNAc4NeuAc2; MS share 0.124-0.128 vs FLD 0.081).
* Controls added: bootstrap 95% interval over FLD peaks for TV(protonated) - TV(all), and a
  decoy-grouping control (species built with every relation shifted, so absorbed ions are
  chance) that must do worse than real grouping. Tests: 83 passed.

## 2026-11-20 (10): FLD controls, run 1

* TV(protonated) - TV(all assigned) = +0.037, 95% bootstrap over 13 peaks +0.002 to +0.086.
  Decoy-grouped ions TV 0.091 (protonated 0.092): chance relations absorb almost nothing, so
  this control only repeats the protonated comparison; it does not test whether adding the
  same signal to the wrong peaks would help as much.
* Added a shuffle control: the extra signal per peak (all assigned minus protonated) is moved
  to random peaks 5000 times; one-sided p = fraction of shuffles at least as close to FLD.
* Window note: several FLD windows run the full +/-1 min from the apex (e.g. 35.741-37.728),
  so valleys were not found within 1 min and neighbouring peaks may share a window.
* Status: support is positive but marginal on one run (lower bound +0.002). Needs more runs
  with FLD before it can be claimed. Tests: 84 passed.

## 2026-11-20 (11): shuffle control result

* Extra signal (all assigned minus protonated) moved to random peaks: median TV 0.389 vs
  0.055 as assigned; p = 0.0002 (1/5001, no shuffle as close; the smallest p this test can
  give with 5000 shuffles).
* Reading: two nulls bracket the result. Random placement of the extra signal (shuffle) is an
  easy null, because the extra signal is large and misplacing it distorts the profile badly
  (0.389, worse than protonated ions alone). Placement proportional to each peak's
  protonated signal leaves shares unchanged, i.e. TV = protonated 0.092; against that null the
  gain is +0.037 (95% bootstrap +0.002 to +0.086). Claim for one run: the assigned ions go to
  the right peaks (strong evidence), and they improve agreement with FLD beyond what scaling
  protonated signal gives (marginal evidence).

## 2026-11-23: second run, fetuin on Thermo Orbitrap

* IgG Agilent run ruled out: in MassHunter data, MSPeak.bin holds (m/z, intensity) pairs
  only; scan offsets, retention times and MS level are in MSScan.bin, absent from the deposit
  (checked in rainbow 1.5.5 masshunter.py). Retention times cannot be recovered.
* Fetuin Thermo .raw (Zenodo 4705546, 601,751,870 bytes) converted on the cluster with
  ThermoRawFileParser 1.4.5 (bioconda, mono 6.12.0.199 from conda-forge; Miniforge3/24.11.3-0
  module): indexed mzML, peak picking on, metadata JSON; 44,946 scans, 0 errors, 0 warnings.
* census_run.py reads .mzML input (MS1 in the time window; no FLD, so FLD steps are skipped).
  Tests: 85 passed.

## 2026-11-23 (2): grouping speed, stage log, batch job

* Orbitrap fetuin census (window 10-70 min, noise 3000) still running after 40 min on the
  login node at full CPU, 0.95 GB. Cause: group_species compared every feature with every
  other in Python loops (pass 1) and tested ~380 relations one at a time per species pair
  (pass 2), and the decoy calibration repeats grouping about 50 times.
* Rewrote group_species: candidates found through a retention-time index (window widened by
  1e-9 min, exact test kept), relations tested as arrays, intensities kept in a running
  table. The loop version is frozen in tests/_ref and an equivalence test compares species,
  ion forms and absorption records on 24 synthetic cases (6 seeds x 4 settings: plain,
  all-family decoy, single-family decoy, family cutoffs). 2,568 synthetic features: 5.19 s ->
  0.45 s, identical species.
* census_run.py prints stage times to stderr and flushes output lines; scripts/slurm/
  census.sbatch runs a census as a Sapelo2 batch job (partition batch, 8 GB, 12 h).
* Tests: 109 passed.

## 2026-11-23 (3): Orbitrap fetuin census, job 50515340

* 1880 MS1 scans (10-70 min) -> 8103 features -> 2778 species; run time 6.6 min on a compute
  node. Relation decoys: adduct 14.3%, isotope 28.2%, isf_sia 98.6%, isf_other 16.8%, dimer
  54.1%, re-read 39.9% chance at full width; after cutoffs 43.2% chance share overall.
  Unassigned 92.6% of intensity; one species at q <= 0.05.
* Strongest features (654.3154 3+, 635.3317 2+, 840.4513 3+, 865.6773 4+ with 1153.8997 3+)
  do not look like released glycans. With the method name (AB_22Apr20_GYPeptideID_100min),
  MS2 precursor charges 2-6 and A3S3 2-AB 2+ at about 1e-4 of the top features, the file may
  be a fetuin glycopeptide run. Hypothesis only; the Zenodo title says released glycans on both
  instruments. Test: scripts/ms2_mzml_check.py (share of HexNAc-oxonium MS2 that also show the
  2-AB Y1 ion 342.1660).
* Reading if confirmed: the decoys refused to certify relations on non-glycan-dominated data
  (a negative control for the cleaning step), and the run is out of scope for the released-
  glycan census.
* Bug: recalibration ran on 1 anchor (+9.01 ppm) and that anchor then scored q = 0 against its
  own correction. Fixed: at least 5 anchors, otherwise no correction. Tests: 110 passed.

## 2026-11-23 (4): Orbitrap file is not released glycans

* ms2_mzml_check.py: 38,256 MS2, 22,922 with low mass recorded. At 20 ppm: HexNAc oxonium in
  8,244; HexHexNAc 6,529; NeuAc 5,462; 2-AB Y1 342.1660 in 182; oxonium-positive with Y1: 49 of
  8,244 (0.6%). Oxonium-positive precursor charges: 2+ 1,396; 3+ 2,353; 4+ 2,006; 5+ 1,752;
  6+ 737. At 0.3 Da (looser): 953 of 12,149 (7.8%).
* Reading: glycan-containing spectra are plentiful, but almost none carry the labelled reducing
  end, and most precursors are 3+ to 6+. Consistent with glycopeptides (the GYPeptideID method),
  not released 2-AB glycans. Excluded from the census; catalog note updated. A positive control
  (share of 342.166 in oxonium-positive MS2 of the Agilent released-glycan run) would make the
  test quantitative; to do.
* Second FLD run still needed: next candidates PXD077477 and PXD078897 (Thermo; FLD in the
  raw files unverified), then GlycoPOST and MassIVE searches.

## 2026-11-23 (5): Thermo detector lister works

* tools/thermo_channels built on the cluster (mcs from the trfp environment, netstandard facade
  referenced). Fetuin Orbitrap .raw: 1 MS device; 4 analog devices, all from a Dionex
  NCS-3500RS nano pump (pressure and flow channels); 1 MS analog device with channel "UV 214"
  (volts); no UV, PDA or fluorescence device. A nano-LC run with 214 nm UV fits the
  glycopeptide reading; no FLD ground truth in this file.

## 2026-11-23 (6): second FLD run found (PXD077477)

* PXD077477 (Moon, BI-NEX; interferon-beta R27T, PNGase F): three released-glycan .raw files
  (2-AB 40C 600-2500; procainamide 40C 600-2500; procainamide 40C); the rest are peptide maps.
  Downloaded 20231212_R27T_2AB_40C_600-2500.raw over HTTPS, 331,967,983 bytes (matches PRIDE).
* Detectors (tools/thermo_channels): MS; Dionex DGP-3600RS pump pressures; Dionex FLD-3400RS
  with UV-type channel Emission_1 (28,501 points, 0-95 min) and analog FLD_FlowCell; TCC-3200
  oven and cooler temperatures. SDRF: R27T_IFNB_PNGaseF_2AB_glycan, technical replicate 1.
* census_run.py: --fld-csv/--fld-column read an exported trace; --lag auto estimates the
  FLD-to-MS delay by cross-correlation with the MS1 TIC (align.estimate_lag, +/-2 min).
  New io_traces.read_trace_csv (a row with an empty value used to desynchronise time and
  value lists; caught by its test). Tests: 111 passed.

## 2026-11-23 (7): PXD077477 2-AB run inspected

* Traces exported (FLD Emission_1 28,501 points, 0-95 min). Lister crashed on the empty
  "Other" device (scans -1); fixed to skip devices without scans. mzML: 9,356 scans
  (9,201 MS1, 155 MS2), 0 errors.
* FLD: injection front 1.33e7 at 0-5 min, baseline 2.6e3-6.6e3, glycan peaks 25-75 min
  (max 3.7e4-1.6e5 per 5-min bin). MS1 peak intensities low (percentiles 1/10/50/90/99:
  722/855/1200/1752/2705) and the MS1 TIC nearly flat 0-65 min (2.2e8-3.4e8 per 5 min), so
  background dominates the TIC.
* --lag auto now aligns FLD with a glycan-targeted ion current (summed MS1 intensity within
  10 ppm of candidate labelled-glycan m/z, 1+ to 3+), and prints the TIC correlation for
  comparison. Tests: 112 passed.

## 2026-11-23 (8): PXD077477 2-AB census, job 50516575

* 5,429 MS1 scans (22-78 min), noise 1000 -> 96 features -> 47 species; relation decoys all
  0.0 (too few features to calibrate); recalibration on 8 anchors (+0.57 ppm, anchors 0.48 ->
  0.16 ppm). 33 of 47 species at q <= 0.05 (79.8% of intensity); 25 at q <= 0.01 (71.8%).
  Lag from glycan ion current -0.010 min (r 0.560; TIC r 0.299).
* FLD validation, 6 glycans in 6 peaks: TV strongest ion 0.234, protonated 0.330, all ions
  0.266; difference +0.064 (95% -0.081 to +0.184); shuffle p = 0.27. Inconclusive at n = 6.
  Sialic acid slope: protonated -0.233, all ions -0.127 (cleaning flattens it here; in the
  fetuin run it overshot to +0.052).
* Finding: the C3H9N-RA and NH3-RA end groups give identical formulas through
  NeuGc + C3H9N - O = HexNAc + dHex + NH3 - O (C14H26N2O8). The strongest species,
  1787.6822 (1+, 38.17 and 40.02 min), 2078.7786 and 2369.8735 differ by NeuAc
  (291.0964, 291.0949) and read as the NH3-RA forms of FA2, A2G2S1F and A2G2S2F, the 2-AB
  series of this protein. A second lab's data (Korean lab, IFN-beta) showing the +1.0316 Da
  series supports the amino-alditol hypothesis. Proposed mechanism: PNGase F glycosylamine
  reduced by NaCNBH3 during labelling. Precedent for the product class: Figl and Altmann,
  Proteomics 2018 (doi:10.1002/pmic.201700330), under reductive alkaline release.
* New: resolve_by_label_support flags isomeric readings whose composition has no labelled
  partner in the run. Tests: 113 passed.
* Feature count too low (96): MS1 intensities sit near the floor (median 1200); rerun with
  noise 500.

## 2026-11-23 (9): reruns at commit 516dc67

* PXD077477 2-AB, noise 500 (job 50516741): 185 features -> 92 species; relation chance
  1.9%; 12 species resolved by labelled-glycan support; 34 of 92 species at q <= 0.05 (60.8%).
  Feature intensity: labelled glycan 7.0%, other end group 26.6%, adduct ions 31.3%, dimer
  3.0%, ambiguous 4.2%, unassigned 26.2%. NH3-RA forms outweigh 2-AB forms in MS signal; not
  a molar statement (primary amine vs 2-AB amide ionize differently; no FLD for NH3-RA).
* Retention offsets 2-AB -> NH3-RA: FA2 26.22 -> 38.17 (+12.0); A2G2S1F 32.04 -> 42.41-44.68;
  A2G2S2F 37.62 -> 47.01 / 48.51. Series test n = 0 because max_shift was 6 min (set on the
  fetuin gradient). Rerun with --max-shift 15.
* FA2 NH3-RA form shows two peaks (38.17, 40.01 min) though FA2 has no linkage isomers.
  Test with MS2: Y1 amino-alditol + Fuc 369.1867 (core-fucosylated reading) vs Y1 C3H9N tag
  265.1758 (other reading); scripts/ms2_precursor_ions.py.
* FLD validation, 8 glycans in 8 peaks: TV strongest 0.262, protonated 0.352, all ions 0.275
  (difference +0.078, 95% -0.108 to +0.190; shuffle p = 0.19). Inconclusive.
* Fetuin rerun (job 50516755): unchanged except 3 species resolved by support; ambiguous
  0.1% -> 0.0%.

## 2026-11-23 (10): MS2 supports the amino-alditol reading in PXD077477

* PXD077477 MS2 are HCD with a low-resolution readout (e.g. HexNAc oxonium at 203.86-203.93,
  HexHexNAc at 365.90-366.02). At 20 ppm nothing matched; at +/-1000 ppm (about +/-0.3 Da):
  2-AB FA2 control shows Y1 2-AB 342 (0.148) and 488 (0.071), no 223/369; 2078.78 and
  2369.87 (8 spectra) show 223.13 (0.076-0.106) and 369.19 (0.049-0.130), 265 at most 0.012,
  342 at most 0.004. Written up in docs/findings/2026-11-23_amino_alditol_two_labs.md.

## 2026-11-23 (11): NH3-RA series passes in PXD077477 (job 50517166)

* With --max-shift 15: NH3-RA vs 2-AB n = 5, median ratio 1.817 (pairs 0.910-8.616), MAD
  0.074, shift +9.40 min, fluorescence fraction 0.09 (5 clear pairs; tail shape from 3
  isolated FLD peaks) -> SERIES. 8 species resolved by labelled-glycan support. Findings doc
  updated. Evidence in this run: mass (<0.5 ppm), retention offset, no fluorescence, MS2 Y1.
* Results now copied to the user's PC with scp via xfer.gacrc.uga.edu.

## 2026-11-23 (12): literature, full texts

* Wang et al. J Proteomics 2018 (doi:10.1016/j.jprot.2018.06.002): glycosylamine -> open-ring
  imine -> NaBH3CN -> 1-amino-alditol; also ammonium acetate + NaBH3CN at 70 C.
* Figl and Altmann, Proteomics 2018: fetuin trisialo 1-amino-1-deoxy-glycitol at m/z 961.02
  (3+), the same ion as the fetuin census species.
* Cosenza et al. ARKIVOC 2011: Amadori route to epimeric aminoalditols in acetic acid; a
  hypothesis for the paired NH3-RA peaks.
* Moon et al. IJMS 2026 (PXD077477): glycans dried after EB10 clean-up and labelled with Ludger
  kits; no hydrolysis step stated.
* The IgG .d zip uploaded again is the deposit already examined (no MSScan files); not used.

## 2026-11-23 (13): PXD077477 procainamide runs (jobs 50517803, 50517804)

* 20231212_R27T_PROC_40C_600-2500: 444 features -> 198 species; relation chance 1.1%;
  recalibration on 12 anchors (0.41 -> 0.43 ppm). 69 species at q <= 0.05 (82.6% of
  intensity). Lag -0.008 min (r 0.611 glycan current; TIC r 0.736 here). Labelled glycan
  55.4% of feature intensity, adduct ions 27.1%, unassigned 10.0%.
  FLD validation: 53 glycans in 38 peaks. TV strongest 0.346, protonated 0.361, all ions
  0.295; difference +0.066 (95% +0.027 to +0.111); shuffle p = 0.0028. Sialic slope
  -0.400 -> -0.340.
* 20231120_R27T_PROC_40C: 585 features -> 221 species; relation chance 5.1% (isf_sia off,
  55.6% chance); recalibration on 17 anchors (0.67 -> 0.28 ppm); 91 species at q <= 0.05
  (76.6%). FLD validation: 53 glycans in 35 peaks. TV 0.319 / 0.355 / 0.299; difference
  +0.056 (95% +0.025 to +0.101); shuffle p = 0.018. Sialic slope -0.279 -> -0.211.
* FLD validation across runs, TV(protonated) - TV(all assigned):
  fetuin 2-AB +0.037 (+0.002, +0.085), p 0.0002, 13 peaks;
  IFN-beta 2-AB +0.078 (-0.108, +0.190), p 0.19, 8 peaks;
  IFN-beta ProA a +0.066 (+0.027, +0.111), p 0.0028, 38 peaks;
  IFN-beta ProA b +0.056 (+0.025, +0.101), p 0.018, 35 peaks.
  Same sign in 4 of 4; interval above zero in 3 of 4. Residual MS-vs-FLD distance after
  cleaning stays large for ProA (about 0.30): cleaning reduces, does not remove, MS bias.
  Sialic acid trend: fetuin +0.052 after cleaning; IFN-beta ProA -0.34 and -0.21
  (sialylated glycans under-represented in MS); direction depends on run.
* NH3-RA vs ProA: run a n = 7, ratio 0.072, MAD 0.102, shift +8.49 min, FLD fraction 0.29
  (1 clear pair); run b n = 15, ratio 0.173, MAD 0.108, shift +8.22 min, FLD fraction 0.37
  (7 clear pairs). Both "not a series" only because MAD exceeds the preset 0.10. Threshold
  not changed after seeing data. Form present under procainamide with a consistent shift, so
  it does not depend on the 2-AB reagent. Ratio to the label form is much lower with ProA
  (0.07-0.17) than with 2-AB (1.8): ProA carries a tertiary amine and ionizes far better.
* Exchangeability 3.13 and 2.95 in ProA runs (24 and 23 target vs 7.7 and 7.8 decoy in
  15-30 ppm). The band holds an intense series 2955.15 / 3101.22 / 3392.32 / 3683.41
  (NeuAc steps) about +17 to +19 ppm from C3H9N-RA readings, at several retention times; a
  real unexplained species family, not chance. Open.
