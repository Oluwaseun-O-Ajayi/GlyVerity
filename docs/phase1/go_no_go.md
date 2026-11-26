# Phase 1: go/no-go (interim, 2026-11-14)

## Recommendation: CONDITIONAL GO

Proceed to Phase 2, with the gap statement narrowed and three conditions checked before
Dec 1.

## Evidence for going ahead

1. The community names the problem. Dagstuhl 25351 (doi:10.4230/DagRep.15.8.46) states that
   glycan identification "often lacks robust FDR control" and asks for composition-level FDR,
   GlyTouCan IDs and explicit structural specificity.
2. Reporting standards now ask for what GlyVerity Report produces. MIRAGE MS v1.0 (2025;
   doi:10.1016/j.mcpro.2025.101473) asks for GlyTouCan IDs and search-space details.
3. The newest MS1 composition tool, GlycoGenius (Nat. Commun. 2025), filters on fixed score
   thresholds with no decoys or q-values described, and does not address in-source fragments or
   label side-products.
4. No tool checked models label side-product or label-solvent adduct series. Your GPST000788
   data contain such a series (ProA plus acetonitrile, +260.2001 Da, which matches
   219.1735 + 41.0265 to four decimals in `tests/test_chemistry.py`).

## Evidence that weakens the original case

1. Decoys for glycan compositions are not new. GlycReSoft (2012) used randomized elemental
   compositions to estimate FDR at MS1; GlycoNote, GlycanDIA Finder and GlycomicsAT use decoys at
   MS2. "First FDR for glycomics" cannot be claimed.
2. The khipu result was misstated in the brief. In metabolomics, in-source fragments are under
   10% of features; most abundant features are isotopologue and adduct patterns of real
   compounds. The census framing still works, but the expected answer for glycomics is open and
   may be "fewer artifacts than feared".
3. Fluorescence-paired public raw data may be scarce. Only one candidate (GPST000246) has been
   seen so far besides your Zenodo set, and its MS mode is unconfirmed. If fewer than about five
   independent paired sets exist, the "decade census validated against fluorescence" design
   must change (see conditions).
4. Your own dataset is embargoed to 2027-11-14, after the planned submission.

## Refined gap statement

"Released-glycan LC-MS tools assign compositions after limited or no grouping of related ions,
and where error rates are estimated, the decoys have not been tested for calibration against an
independent reference. GlyVerity first explains features that arise from isotopes, adducts,
in-source fragments and reducing-end label side-products, then assigns the remaining features
with composition-level q-values whose calibration is checked against fluorescence-quantified
profiles and characterized standards."

## Conditions to meet before Dec 1

1. Ground truth count. Harvest PRIDE (script ready), browse GlycoPOST and MassIVE, and count
   in-scope datasets with an independent reference (fluorescence, characterized standards such
   as GPST000578, or spike-ins). Go if there are at least 3 independent references across at
   least 2 labels; otherwise reframe validation around standards and simulation.
2. Census size. Count in-scope downloadable LC-MS datasets by year, lab, instrument and label.
   The "decade" claim stands only if the years in the catalog span it.
3. Decoy design. Decide how MS1 decoys will be generated so that they share the property
   distributions of targets (mass defect, isotope pattern, retention behaviour). The FragPipe
   MBG design (+11 Da, shared RT and IM) is the nearest precedent; test whether +11 Da decoys
   are trivially rejected by isotope-pattern scoring before adopting it. Burger, Ge and Li
   (Nat. Biotechnol. 2026) name decoy exchangeability as the assumption target-decoy methods
   must justify; GlyVerity will report that check with every result. If decoys cannot be made
   exchangeable, the fallback is an empirical-null estimate from standards or entrapment
   compositions.

## Update 2026-11-14 (evening): PRIDE harvest screened

Source: `data/catalog/harvest/` (PRIDE search, 8 keywords, 474 unique projects), screened by
`scripts/screen_pride.py` into `data/catalog/pride_screened.csv`, then checked by hand against
each record's protocols and file list. Counts below are from that check.

* 17 PRIDE projects are in scope (released-glycan LC-MS with raw files listed), spanning
  2016 to 2026 with 2019 missing. The decade span is met from PRIDE alone, before GlycoPOST and
  MassIVE are added.
* Labels in those 17 (plus two with label not yet confirmed): permethylated 4, reduced 4,
  native 3, 2-AB 3, reduced and permethylated 2, ProA and 2-AB 1, RapiFluor-MS 1,
  isobaric tag 1.
* Condition 1 (ground truth), progress: two PRIDE projects run fluorescence and MS on the same
  LC, PXD077477 (ProA and 2-AB, Orbitrap Fusion, Ex/Em 310/370 nm) and PXD078897 (RapiFluor-MS,
  HILIC-FLD-MS, Q Exactive HF-X). Together with the Zenodo 2-AB set and the GPST000578
  standards, that covers 3 labels. Not yet met: whether the fluorescence channel is stored in
  the deposited raw files has to be checked by opening one file from each.

### New finding that changes scope

Labelled HILIC, the biopharma QC workflow the brief targets, is a small share of public
repository data (3 of 17 in-scope PRIDE projects: PXD003498 native HILIC, PXD077477,
PXD078897). Most public released-glycan LC-MS data are research workflows: PGC with reduced
glycans, or C18 with permethylated glycans. Two options:

1. Keep labelled HILIC as the validation core (where fluorescence ground truth exists) and
   make the census label-agnostic, so permethylated, reduced and native data count too.
   The chemistry module already separates label from composition, so this costs little.
2. Restrict the census to labelled HILIC, accepting a much smaller census.

Decision (2026-11-14): option 1. Labelled HILIC is the validation core; the census includes all labels.

## Update 2026-11-15: first real-data evidence

In one public fetuin 2-AB run (Zenodo 4705546), two artifact types appear that fluorescence
cannot see and composition search misreports: a non-fluorescent reducing-end series
(+43.0786 Da, 30.4% of tag plus 2-AB MS1 signal) and six false "free fucosylated" glycans that
are sodium-adduct isotope peaks. Details: `docs/findings/2026-11-15_zenodo4705546_fetuin.md`.
This supports the gap statement on real data. It is one run; the census must show how common
these are.

## Things to do outside the code

* Decide whether to shorten the GPST000788 embargo or arrange reviewer access.

## Update 2026-11-20: condition 3 (decoy design), first result

* Composition decoys: target formula masses plus or minus 3, 5, 7 or 9 H (octet-rule
  violation; JUMPm, Wang et al., J. Proteome Res. 2018, doi:10.1021/acs.jproteome.8b00019).
  They keep the CHNO mass defect, and their isotope patterns are those of the target
  formulas, so isotope-pattern scoring cannot reject them trivially (the concern raised for
  +11 Da decoys above). Separate-search FDR, decoy counts scaled for decoys removed near
  targets.
* Relation (cleaning) decoys: per-family offset shifts, one family at a time; tolerances set at
  <= 5% chance per family; final chance share 3.7% on the fetuin run.
* Fetuin run after recalibration on 21 multi-charge anchors (+1.68 ppm offset; anchor error
  1.68 -> 0.32 ppm): 56 of 424 species at q <= 0.05 carry 81.3% of feature intensity;
  labelled-glycan species 50.0% at q <= 0.05 vs 0.6% above.
* Exchangeability check (15-30 ppm band): 17 target vs 12.0 decoy (ratio 1.42). Not settled
  on one run; band species include a recurring unknown near 2883 to 2884 Da. Test across
  datasets; if the ratio stays above 1, add entrapment compositions (the stated fallback).
* Status: design decided and implemented; exchangeability not yet shown.

## Update 2026-11-25: status of the three conditions before Dec 1

### Condition 1, ground truth: met as written, narrow

Fluorescence traces are stored in the deposited files and are read by GlyVerity (Agilent .d
directly; Thermo .raw through tools/thermo_channels). Four runs with FLD and MS on the same LC:

| Run | Lab | Instrument | Label | FLD peaks | TV(protonated) - TV(all assigned), 95% CI | Shuffle p |
|---|---|---|---|---|---|---|
| Fetuin, Zenodo 4705546 | Charles River | Agilent Q-TOF | 2-AB | 13 | +0.037 (+0.002, +0.085) | 0.0002 |
| IFN-beta, PXD077477 | BI-NEX | Orbitrap Fusion | 2-AB | 8 | +0.078 (-0.108, +0.190) | 0.19 |
| IFN-beta, PXD077477 (a) | BI-NEX | Orbitrap Fusion | ProA | 38 | +0.066 (+0.027, +0.111) | 0.0028 |
| IFN-beta, PXD077477 (b) | BI-NEX | Orbitrap Fusion | ProA | 35 | +0.056 (+0.025, +0.101) | 0.018 |

Three independent references (fetuin 2-AB; IFN-beta 2-AB; IFN-beta ProA) across two labels:
the stated bar. Only two laboratories, though. PXD078897 (RapiFluor-MS, Q Exactive HF-X, a
third lab and label) is next. Figure: docs/figures/fig_fld_validation.pdf.

Result: summing the ions GlyVerity assigns to each species moves MS shares toward FLD shares
in 4 of 4 runs; the interval excludes zero in 3 of 4, and moving the same ions to random
peaks makes agreement worse. Limits: after cleaning, MS still differs from FLD by TV about
0.30 in the ProA runs, and the MS bias with sialic acid count changes sign between runs, so
cleaning reduces MS bias but does not remove it.

### Condition 2, census size: catalog met, census not yet run

17 in-scope PRIDE projects, 2016 to 2026 (2019 missing); GlycoPOST and MassIVE not yet
searched. One file was out of scope on inspection (Zenodo 4705546 Orbitrap file:
glycopeptide-like MS2, nano-LC with UV 214, no FLD). Thermo .raw files are now handled on the
cluster (ThermoRawFileParser 1.4.5, mzML input, batch script), which most of the catalog needs.

### Condition 3, decoy design: design done, exchangeability not yet shown

* Relation (cleaning) decoys: offsets shifted one family at a time, tolerance per family set
  at <= 5% chance on a score combining ppm error and retention difference. Achieved chance
  share: 3.7%, 1.9%, 1.1%, 5.1% (four FLD runs). A run dominated by non-glycan features
  (glycopeptide file) gave 43%, i.e. the method declined to certify relations there.
* Composition decoys: odd-hydrogen (JUMPm design), counts scaled for removed decoys, mass
  recalibration on >= 5 multi-charge anchors. q <= 0.05 species carry 81.3%, 60.8%, 82.6%,
  76.6% of feature intensity.
* Exchangeability (15-30 ppm band, ratio 1 expected): 1.42, 0.31, 3.13, 2.95. In the ProA runs
  the excess comes from one intense family (2955.15, 3101.22, 3392.32, 3683.41 Da; NeuAc steps;
  17-19 ppm from the nearest reading), real chemistry outside the search space, not chance.
  The check must exclude explained-but-unlisted families or the search space must cover them.
  Open; the stated fallback (entrapment compositions) remains available.

### Finding that strengthens the case

An unlabelled 1-amino-1-deoxy-glycitol series (+NH3 - O, +1.0316 Da on the free glycan) appears
in both labs' data and under both reductive-amination labels. Evidence: mass (< 0.5 ppm),
consistent retention offset (+4.6 min fetuin; +8 to +9.4 min IFN-beta), FLD at 9% of
expectation where tails could be corrected (IFN-beta 2-AB), and MS2 Y1 ions 223.13 and
369.19 absent from the in-run 2-AB control. In the IFN-beta 2-AB run these forms carry more MS
signal than the labelled glycans (26.6% vs 7.0% of feature intensity) and are absent from the
published profile. The product class and its NaBH3CN chemistry are known (Figl and Altmann
2018; Wang et al. 2018); its presence as an unrecognized artifact of routine PNGase F +
reductive-amination profiles is not reported in the papers checked so far. Details:
docs/findings/2026-11-23_amino_alditol_two_labs.md.

Testable prediction: RapiFluor-MS labels the glycosylamine directly with no reductant, so the
+1.0316 Da series should be absent or far weaker in PXD078897. A positive result there would
argue against the reduction mechanism.

### Recommendation: GO for Phase 2

Conditions 1 and 3 have working implementations with real-data results; the remaining gaps
(a third lab, exchangeability, census breadth) are Phase 2 and Phase 3 work, not reasons to
stop. Phase 2 priorities in order: PXD078897 (RFMS, third lab); exchangeability check that
handles unexplained real families; then the census pipeline over the catalog.
