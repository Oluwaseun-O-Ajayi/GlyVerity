# Phase 1: prior-art check

Status as of 2026-11-14. Each entry was checked against the publisher page, PubMed Central,
the preprint server or the project repository. "Pending" means the identifier has not yet been
confirmed and must not be cited until it is.

## Question asked of each tool

1. Does it assign glycan compositions from MS1 (LC-MS) data?
2. Does it estimate an error rate (decoys, FDR, q-values), and at what level?
3. Which artifacts does it model: isotopes, adducts, in-source fragments (ISF), label side-products?
4. Which reducing-end labels does it support?

## Released-glycan LC-MS tools (closest to GlyVerity)

| Tool | Reference | Level | Error rate | Artifacts modelled | Notes |
|---|---|---|---|---|---|
| GlycReSoft (2012) | Maxwell et al., PLoS ONE 2012, doi:10.1371/journal.pone.0045474 | MS1, LC-MS | Yes. False positives and FDR estimated with randomized elemental compositions as decoys | Charge states and adducts summarized per composition; input deconvoluted by Decon2LS/DeconTools | Earliest MS1 composition decoy approach found. Must be cited as direct precedent |
| GlycReSoft network smoothing | Klein, Carvalho, Zaia, Bioinformatics 2018, doi:10.1093/bioinformatics/bty397 | MS1, LC-MS | Likelihood scoring with network (biosynthetic) regularization | Isotopic fit, peak shape, mass accuracy, ISF correlation cited from earlier work | Reads mzML/mzXML. Composition network idea overlaps GlyVerity's planned biosynthetic prior |
| GlyQ-IQ | Kronewitter et al., Anal. Chem. 2014, doi:10.1021/ac501492f | MS1 targeted, LC-MS | No FDR found; 99% of assignments passed manual validation | Isotope profiles, glycan family relationships, in-source fragmentation | Human serum native N-glycans |
| MultiGlycan | Yu et al., Bioinformatics 2013, doi:10.1093/bioinformatics/btt190 | MS1 (LC-MS, MALDI) | No FDR found | Overlapping isotope envelopes; ESI adducts | Windows C# |
| MultiGlycan-ESI | Hu et al., Rapid Commun. Mass Spectrom. 2015, doi:10.1002/rcm.7093 | MS1, LC-ESI-MS | No FDR found | Merged adducts and charge states | Permethylated glycans |
| GlycoGenius | Loponte et al., Nat. Commun. 2025, 16:10335, doi:10.1038/s41467-025-65265-2 | MS1 plus MS2 annotation | No FDR, q-value or decoy described in the preprint text (doi:10.1101/2025.03.10.642485). Uses fixed thresholds: isotopic fit 0.8, curve fit 0.5, ppm plus or minus 20, S/N 3 | Isotopes and user-set adducts. ISF not addressed. Label side-products not mentioned | Supports 2-AB, 2-AA, ProA, GirP, permethylation, reduced end. mzML input. Closest active competitor; benchmark against it |
| GlycanDIA Finder | Xie et al., Nat. Commun. 2025, 16:7075, doi:10.1038/s41467-025-61473-y | MS1-centric DIA with MS2 | Yes. Random mass shift applied to theoretical fragment ions, 100 decoy sets, 2% FDR | Isotope peak check; adducts and ISF not mentioned | Native glycans, PGC, Orbitrap. Data MassIVE MSV000093677. Decoys are MS2-level |
| GlycoNote | Anal. Chem. 2023 (described in the GlycanDIA paper as earlier work by the same team) | MS/MS | Target-decoy with iterative decoy searching; random mass shift | Not checked | DOI pending |
| GlyCombo | Kelly and Ashwood, J. Am. Soc. Mass Spectrom. 2024; preprint doi:10.26434/chemrxiv-2024-3j008 | Precursor composition enumeration | No FDR | Multiple adducts; off-by-one isotope errors | Apache-2.0. Labels incl. 2-AB, ProA, RapiFluor-MS, InstantPC. Journal DOI pending |
| GlycoAnnotateR | J. Am. Soc. Mass Spectrom. 2025, doi:10.1021/jasms.5c00093 | De novo composition annotation, MS1 | Reports annotations at an FDR in one use case (MALDI-MSI); method not yet read | Not checked | R package |
| GlycomicsAT | bioRxiv 2025, doi:10.1101/2025.07.02.662857 | MS/MS | Decoys from scrambled monosaccharide masses, total mass kept | Not checked | Transformer-based |

## Adjacent precedents (not released glycans)

| Item | Reference | Why it matters |
|---|---|---|
| FragPipe Match-Between-Glycans | bioRxiv 2026, doi:10.64898/2026.02.18.706650 | MS1-level decoys for glycopeptides: decoy shares RT and IM with target and adds 11 Da; LDA scoring; entrapment check gave 0.6% FDR for implausible glycans. Closest MS1 decoy design to adapt |
| pGlyco | Zeng et al. (PMC4853738) | Spectrum-based glycan decoys for glycopeptides |
| GlycReSoft glycopeptide engine | Klein, Carvalho, Zaia, Nat. Commun. 2024, doi:10.1038/s41467-024-50338-5 | Joint peptide and glycan FDR; glycan network smoothing |
| HUPO HGI first study | Kawahara et al., Nat. Methods 2021, doi:10.1038/s41592-021-01309-x | Community benchmark design, glycoproteomics only |

## General FDR guidance

| Item | Reference | Why it matters |
|---|---|---|
| FDR taxonomy for omics | Burger, Ge, Li, "A unified and practical taxonomy for false discovery rate control in omics research", Nat. Biotechnol. Comment, 13 Nov 2026 (DOI pending) | Groups FDR methods as traditional, empirical-null, or target-decoy/knockoff. GlyVerity's composition FDR is in the target-decoy family, whose key assumption is decoy exchangeability. Cite it to justify the decoy checks in go/no-go condition 3. Not glycomics-specific; builds no tool |

## Libraries GlyVerity can build on

| Library | Reference | Licence | Use |
|---|---|---|---|
| glypy | Klein and Zaia, J. Proteome Res. 2019, doi:10.1021/acs.jproteome.9b00367 | Apache-2.0 | Composition and structure objects; GlyTouCan access |
| glycowork | Thomes, Burkholz, Bojar, Glycobiology 2021, doi:10.1093/glycob/cwab067 | Check | Biosynthetic feasibility priors (later layer) |
| CandyCrunch | Urban et al., Nat. Methods 2024, doi:10.1038/s41592-024-02314-6 | Check | MS/MS structure prediction; topology evidence (later layer) |
| khipu | Li and Zheng, Anal. Chem. 2023 (DOI pending); github.com/shuzhao-li/khipu | Check | Pre-annotation by mass-difference trees; supports derivatization links |

## Not yet checked

GlycoWorkbench, GlycoMod, Skyline glycan workflows, Byos, SimGlycan, LaCyTools/MassyTools,
nQuant, GlyGen and GlyCosmos analysis tools. These go on the Phase 1 task list before Dec 1.

## What this means for the gap

* "Glycomics has no FDR" is false as a blanket claim. Decoy-based estimates exist at MS1
  (GlycReSoft 2012) and MS2 (GlycoNote, GlycanDIA Finder, GlycomicsAT).
* Not found in any tool checked so far: a pre-annotation stage that explains label
  side-products and label-solvent adduct series before assignment, combined with composition
  q-values whose calibration is tested against independent ground truth. This is the defensible
  gap, stated in `go_no_go.md`.
* Absence in our search is not proof of absence. The statement is limited to the tools listed.

## Update 2026-11-23: amino-alditol (+NH3 - O) side product

* Background, verified: PNGase F releases N-glycans as glycosylamines (e.g. US 12,092,644;
  Tarentino et al. cited there). Reductive amination labels such as 2-AB use NaCNBH3.
* Closest precedent found: Figl R, Altmann F. Reductive alkaline release of N-glycans
  generates a variety of unexpected, useful products. Proteomics 2018.
  doi:10.1002/pmic.201700330. Under reductive alkaline release (1 M NaBH4), about one third of
  N-glycans appeared as 1-amino-1-deoxy-glycitols, and glycosylamines were seen when acidic
  conditions were avoided.
* Not found yet: a report of 1-amino-1-deoxy-glycitols formed during PNGase F release
  followed by 2-AB reductive amination (glycosylamine reduced by NaCNBH3 before hydrolysis).
  Quick search only; search further before any novelty claim.
