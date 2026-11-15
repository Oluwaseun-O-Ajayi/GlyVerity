# Phase 1: verification of the claims in the project brief

Checked 2026-11-14. Each claim is marked VERIFIED, CORRECTED or UNVERIFIED.

## 1. Dagstuhl Seminar 25351: VERIFIED

* Report: Gundry, Palmblad, Wilhelm (eds.), "Computational Proteomics (Dagstuhl Seminar 25351)",
  Dagstuhl Reports 15(8):46, doi:10.4230/DagRep.15.8.46. Seminar held 24 to 29 August 2025.
* The seminar is titled Computational Proteomics; glycomics is one working group (Section 3.3).
* Section 3.3.3 states: "Glycan/glycopeptide identification often lacks robust FDR control"
  and proposes multi-level FDR (PSM, composition, topology).
* Section 3.3.2 recommends GlyTouCan IDs at minimum at composition level and explicit
  declaration of structural specificity.
* Section 3.3.4 recommends benchmarking with known perturbation datasets for normalization and
  differential analysis. The report does not use the term "ground truth" in this section.

## 2. khipu "dark matter" preprint: CORRECTED

* Chi, Mitchell, Zheng, (Thapa in a later version), Li. bioRxiv doi:10.1101/2025.02.04.636472.
  First version titled "Systematic pre-annotation explains the 'dark matter' in LC-MS
  metabolomics"; a later version is titled "Assessing the metabolomics 'dark matter' by a
  detectable khipu model".
* What it found (61 public LC-MS datasets of human plasma/serum): in-source fragments are under
  10% of features; most abundant features have identifiable ion patterns (isotopologues,
  adducts); most features come from real compounds, and the number of compounds is much
  smaller than the number of features.
* The brief said most high-abundance features are "isotopes/adducts/fragments/contaminants".
  That overstates it. Accurate wording: "In LC-MS metabolomics, most abundant features belong to
  ion patterns (isotopologues and adducts) of a smaller set of compounds, and in-source fragments
  account for under 10% of features (Chi et al.)."
* "Nobody has asked this of glycomics": no such study found in searches to date. Keep as
  "we found no equivalent analysis for glycomics", not as a certainty.

## 3. FDA September 2025 guidance: VERIFIED (title and wording corrected)

* Document: "Development of Therapeutic Protein Biosimilars: Comparative Analytical Assessment
  and Other Quality-Related Considerations; Guidance for Industry", final, September 2025,
  docket FDA-2019-D-2102, Federal Register notice 2025-17226 (9 September 2025). It replaces the
  2019 draft and the 2015 "Quality Considerations" guidance.
* Checked in the FDA PDF (page numbers are PDF pages of 26, checked by OOA 2026-11-15):
  * p. 14: characterization methods "do not necessarily need to be validated", but should be
    scientifically sound, fit for intended use, and give results that are "reproducible and
    reliable". This applies to characterization methods in general, not to glycosylation
    profiling specifically.
  * p. 11: sponsors should understand product heterogeneity, citing "the nature, location, and
    levels of glycosylation" as an example.
  * p. 12: functional assays should be read alongside glycosylation and other attributes.
* The brief's wording ("stresses glycosylation profiling that is reproducible and reliable")
  joins two separate statements. Accurate wording: "FDA's 2025 biosimilar guidance asks that
  characterization methods give reproducible and reliable results and that sponsors understand
  the nature, location and levels of glycosylation."
* Do not confuse with the separate 29 October 2025 draft guidance on comparative efficacy studies.
* FDA guidances have no DOI. Cite by title, date and docket number.

## 4. MIRAGE MS guidelines: VERIFIED with date detail

* Guideline document: "MIRAGE Guidelines for Reporting Mass Spectrometry Experiments for Glycan
  and Glycoproteomics Data", Version 1.0, dated 10 October 2025 (Beilstein-Institut PDF posted
  26 November 2025). Replaces the 2013 MIRAGE MS guidelines (doi:10.3762/mirage.2).
* Paper: Lisacek et al., "Update and New Implementation of the MIRAGE Reporting Guidelines for
  Mass Spectrometry Experiments in Glycoscience", Mol. Cell. Proteomics 2025,
  doi:10.1016/j.mcpro.2025.101473 (published 24 November 2025).
* The guideline asks for GlyTouCan IDs for reported glycans, which matches GlyVerity Report.

## 5. GlycoPOST: VERIFIED

* GlycoPOST paper: Nucleic Acids Res. 2021, 49(D1):D1523 to D1528, PMID 33174597 (DOI pending).
* Size: 357 projects and 20,864 files as of 25 August 2024 (Anal. Bioanal. Chem. 2024,
  doi:10.1007/s00216-024-05673-3).

## 6. HUPO Human Glycoproteomics Initiative: CORRECTED (scope)

* First HGI study: Kawahara et al., Nat. Methods 2021, doi:10.1038/s41592-021-01309-x.
  It evaluated glycopeptide search software, not released glycans. Cite it as a model of
  community benchmarking, not as evidence about glycomics. A second HGI study is under way
  per the HUPO site.

## 7. NIST interlaboratory mAb study: VERIFIED, raw data UNVERIFIED

* De Leoz et al., Mol. Cell. Proteomics 2020, 19(1):11 to 30, doi:10.1074/mcp.RA119.001677.
* 76 laboratories, 103 reports, 116 compositions reported in total, 57 with consensus abundance
  values; per-laboratory counts ranged from 4 to 48 compositions.
* Raw data availability: not found yet. Next step is to read the data availability statement
  and supplementary files.

## 8. Zenodo 2-AB set (doi:10.5281/zenodo.4705546): VERIFIED

* Smith R., "Released N-glycans from IgG and Fetuin on Agilent QTOF and Orbitrap", Zenodo,
  version v1, published 20 April 2021; concept DOI 10.5281/zenodo.4705545. Data prepared by the
  Biologics MS group, Charles River Labs. Licence CC BY 4.0.
* Samples: 2-AB labelled N-glycans from human serum IgG and from fetuin (fetal bovine serum),
  released with PNGase F. Three files (998.5 MB): fetuin Agilent run
  `MASOP00220_AB_CRL_07JAN2021_Fetuin_001.d.zip` (384.9 MB), IgG Agilent run (11.9 MB) and
  fetuin Orbitrap run `AB_23Apr20_0074_test_Fetuin.raw` (601.8 MB). MD5 sums are in
  `data/catalog/checksums/zenodo_4705546.md5`.
* The IgG Agilent run (`MASOP00011_AB_0332_16DEC2020_hIgG_NoDE_40ul_001.d.zip`) contains a
  fluorescence channel (`FLD1.cg`) and MS data (`MSPeak.bin`). The fetuin run's file list is
  still to be checked for an FLD channel.

## 9. Brief's problem statement: CORRECTED

* "Reports glycan compositions with no error rate" is true of routine practice in the tools
  reviewed, but decoy methods do exist (see `prior_art.md`). The paper must cite GlycReSoft 2012,
  GlycoNote, GlycanDIA Finder and GlycoGenius and state the gap precisely.
