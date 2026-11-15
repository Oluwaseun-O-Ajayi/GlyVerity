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

## Things to do outside the code

* Decide whether to shorten the GPST000788 embargo or arrange reviewer access.
