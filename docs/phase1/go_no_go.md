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
   are trivially rejected by isotope-pattern scoring before adopting it.

## Things to do outside the code

* Read the FDA September 2025 PDF and confirm the "reproducible and reliable" wording and page.
* Confirm the Zenodo 4705546 record contents in a browser.
* Decide whether to shorten the GPST000788 embargo or arrange reviewer access.
