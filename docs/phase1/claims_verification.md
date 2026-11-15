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

## 3. FDA September 2025 guidance: CORRECTED (title)

* Document: "Development of Therapeutic Protein Biosimilars: Comparative Analytical Assessment
  and Other Quality-Related Considerations; Guidance for Industry", final, September 2025,
  docket FDA-2019-D-2102, Federal Register notice 2025-17226 (9 September 2025). It replaces the
  2019 draft and the 2015 "Quality Considerations" guidance.
* A secondary source (Foley and Lardner, December 2025) reports that pages 11 to 12 say
  characterization methods should give results that are "reproducible and reliable", and that
  page 11 covers glycosylation profiling. The FDA PDF itself has not yet been read here.
  Read the PDF and confirm page and wording before quoting.
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

## 8. Zenodo 2-AB set (doi:10.5281/zenodo.4705546): CORRECTED (sample)

* Licence: CC BY 4.0. One Agilent run, `MASOP00011_AB_0332_16DEC2020_hIgG_NoDE_40ul_001.d.zip`,
  acquired with method `AB_18DEC2020_HILIC_2AB`.
* The sample is human IgG, not fetuin as the brief stated.
* The .d folder holds a fluorescence channel (`FLD1.cg`), a DAD channel and MS data
  (`MSPeak.bin`). No profile-mode MS file is listed, so the MS data are likely centroided.
* Title and authors still to be recorded from the record page.

## 9. Brief's problem statement: CORRECTED

* "Reports glycan compositions with no error rate" is true of routine practice in the tools
  reviewed, but decoy methods do exist (see `prior_art.md`). The paper must cite GlycReSoft 2012,
  GlycoNote, GlycanDIA Finder and GlycoGenius and state the gap precisely.
