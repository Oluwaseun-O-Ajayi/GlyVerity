# Amino-alditol (+NH3 - O) forms of released N-glycans in two independent 2-AB datasets

Status: working finding, 2026-11-23. Every number below traces to a run logged in DEVLOG.md.

## Observation

Released, 2-AB-labelled N-glycan runs from two unrelated laboratories contain a series of
species whose neutral masses equal the free glycan + 1.0316 Da (net +NH3 - O), i.e. the
1-amino-1-deoxy-glycitol of the same compositions as the 2-AB glycans.

| Dataset | Instrument | Glycoprotein | Evidence |
|---|---|---|---|
| Zenodo 4705546 (Charles River) | Agilent Q-TOF, HILIC-FLD-MS | bovine fetuin | 8 pairs with the 2-AB series; area ratio 0.17-0.18, spread 0.03-0.04, shift +4.57 min; MS2 of 961.02 (3+): Y1 223.1278, no 342.166 |
| PXD077477 (Moon, BI-NEX) | Orbitrap Fusion, HILIC-FLD-MS | interferon-beta R27T | FA2, A2G2S1F, A2G2S2F forms (1787.6822, 2078.7786, 2369.8735; NeuAc steps 291.0964, 291.0949); retention +9 to +12 min after the 2-AB forms; MS2 below |

## Series and fluorescence test (PXD077477, job 50517166, --max-shift 15)

LabelSideProductSeries NH3-RA vs 2-AB: n = 5, median ratio 1.817 (per pair 0.910 to 8.616;
MAD of log10 0.074), median shift +9.40 min, fluorescence fraction 0.09 over 5 clear pairs
(tail shape measured from 3 isolated FLD peaks) -> SERIES.

| Composition | 2-AB RT | NH3-RA RT | MS ratio | FLD at test apex | FLD expected if fluorescent |
|---|---|---|---|---|---|
| Hex5HexNAc4dHex1 | 26.219 | 38.166 | 1.817 | 9.95e3 | 2.08e5 |
| Hex5HexNAc4dHex1NeuAc1 | 32.042 | 44.683 | 0.910 | 1.13e4 | 1.27e5 |
| Hex5HexNAc4dHex1NeuAc2 | 37.609 | 47.014 | 2.084 | 4.05e4 | 3.01e5 |
| Hex7HexNAc6dHex1NeuAc2 | 50.230 | 59.631 | 8.616 | 2.93e4 | 5.99e5 |
| Hex7HexNAc6dHex1NeuAc3 | 54.086 | 62.858 | 1.533 | 2.87e4 | 1.29e5 |

The forms carry about 9% of the fluorescence a 2-AB-like form would give: consistent with no
fluorophore. In the fetuin run the same test was inconclusive (no isolated FLD peaks; the
fixed-clearance fallback gave 2.64 on peaks that sit on 2-AB tails).

## MS2 test (PXD077477, HCD, low-resolution readout; fragment tolerance about +/-0.3 Da)

Fragment ions relative to base peak. Readings: Y1 amino-alditol 223.1288 and + Fuc 369.1867
(core-fucosylated amino-alditol); Y1 with C3H9N tag 265.1758 (the isomeric alternative,
Hex5HexNAc3NeuGc1-type with the tag); Y1 2-AB 342.1660 and + Fuc 488.2239.

| Precursor (neutral) | RT min | 223.13 | 369.19 | 265.18 | 342.17 | 488.22 |
|---|---|---|---|---|---|---|
| 2078.7786 (2+) | 42.34 | 0.104 | 0.115 | 0.000 | 0.000 | 0.000 |
| 2078.7786 (2+) | 42.92 | 0.099 | 0.103 | 0.000 | 0.004 | 0.000 |
| 2078.7786 (2+) | 43.04 | 0.097 | 0.077 | 0.009 | 0.000 | 0.000 |
| 2078.7786 (2+) | 44.76 | 0.076 | 0.082 | 0.005 | 0.000 | 0.000 |
| 2369.8735 (2+) | 46.91 | 0.106 | 0.130 | 0.012 | 0.000 | 0.000 |
| 2369.8735 (2+) | 47.04 | 0.094 | 0.104 | 0.000 | 0.000 | 0.000 |
| 2369.8735 (2+) | 48.44 | 0.105 | 0.050 | 0.000 | 0.000 | 0.000 |
| 2369.8735 (2+) | 48.56 | 0.100 | 0.049 | 0.000 | 0.000 | 0.000 |
| 1906.7191 (2-AB FA2, control) | 26.18 | 0.000 | 0.000 | 0.000 | 0.148 | 0.071 |

The control shows the 2-AB Y1 ions and none of the amino-alditol ions; the test species show
the amino-alditol Y1 and Y1 + Fuc and essentially no tag or 2-AB Y1. A Na-adduct precursor
of 2369.87 (1196.9332, 2+, 47.05 min) also shows 223.01 and 368.99 among its top peaks.
Both retention peaks of each composition carry the same Y1 pattern. The 1787.68 species was
not selected for MS2.

## Proposed origin

PNGase F releases N-glycans as glycosylamines. If 2-AB reductive amination (NaCNBH3) starts
before the glycosylamine hydrolyses, the reductant can reduce it to the 1-amino-1-deoxy-glycitol
(free glycan - O + N + 3H = +1.0316 Da). Precedent for the product class under a different
release chemistry: Figl R, Altmann F, Proteomics 2018, doi:10.1002/pmic.201700330 (about one
third of N-glycans as 1-amino-1-deoxy-glycitols after reductive alkaline release).

## Literature check, round 2 (2026-11-23)

* Wang et al. (corresponding author wangzhf@nwu.edu.cn), "Reductive chemical release of
  N-glycans as 1-amino-alditols and subsequent 9-fluorenylmethyloxycarbonyl labeling for MS and
  LC/MS analysis", J Proteomics 2018, 187: 47-58 (PMID 29885470; DOI to confirm). Aqueous
  ammonia with NaBH3CN releases N-glycans as 1-amino-alditols on purpose. So NaBH3CN with
  ammonia converts N-glycans to exactly this product.
* This opens a second route besides direct glycosylamine reduction: reductive amination of
  the free glycan by ammonia or ammonium (from glycosylamine hydrolysis, or from ammonium
  buffers often used with PNGase F) competing with 2-AB during labelling. Both routes give
  the same 1-amino-1-deoxy-glycitol (+1.0316 Da), so MS cannot separate them; sample
  preparation records (buffers, clean-up before labelling) could.
* Still not found: a report of this product as an artifact in PNGase F + 2-AB (or other
  reductive amination) released-glycan profiles. Two search rounds so far.
* Also noted: an imidazoline by-product of fluorescent labelling of PNGase F-released IgG
  N-glycans was reported in work on alpha-picoline borane for reductive amination (seen as a
  ResearchGate record only; full reference to obtain before citing).

## PXD077477 sample preparation (from the paper)

Moon D, Kim G, Park M, Park B, Kim NY, Son W, Shin YK, Song K. Relative quantitative analysis
of site-specific N-linked glycosylation in hyperglycosylated interferon-beta via mass
spectrometry. Int J Mol Sci 2026, 27, 7828. doi:10.3390/ijms27177828 (section 4.2):

* R27T produced in CHO cells. 0.2 mg protein, 1 uL 0.5 M DTT and 10 uL 1% RapiGest SF,
  65 C for 30 min; then 2 uL PNGase F (NEB P0705L). No release buffer, incubation time or
  temperature is stated.
* Products purified on a Ludger EB10 cartridge, eluted glycans dried in a centrifugal
  evaporator, then labelled with a Ludger kit (2-AB LT-KAB-A2 or procainamide LT-KPROC-24).
  No glycosylamine hydrolysis step (for example a mild acid step) is described between
  release and labelling. Reductant in these kits: to confirm from the kit sheets.
* HILIC: Waters BEH Glycan Amide 2.1 x 150 mm; A = 100 mM ammonium formate in water,
  B = acetonitrile; 25 to 40% A over 60 min; 0.5 mL/min; 60 C for 2-AB, 40 C for
  procainamide. FLD Ex/Em 330/420 (2-AB) and 310/370 (procainamide). Orbitrap Fusion,
  "resolution 120,000, scan range 600-2500 m/z"; top 15, precursor charge 1-3, HCD.
* The paper reports no unlabelled, +1 Da or amino-alditol species; quantification used FLD and
  XICs of labelled glycans.

Reading: released glycosylamines were dried and labelled with no stated hydrolysis step, which
fits direct reduction of the glycosylamine. No ammonium buffer is stated for the release, so
the ammonia route is not supported by the methods text (ammonium formate is only in the
LC mobile phase, after labelling). The amino-alditol forms carry no fluorophore, so the
paper's FLD-based quantification is unaffected; MS-based (XIC) shares of labelled glycans
would not include them either.

## Limits

* MS intensity is not molar amount; the amine ionizes differently from the 2-AB amide and the
  form does not fluoresce, so its molar share is not known from these data.
* MS2 readout at low resolution; assignments rest on +/-0.3 Da windows, with the in-run 2-AB
  control as the specificity check.
* A literature search beyond one round is needed before any claim that this is unreported for
  PNGase F + reductive amination workflows.
* Exact isomer: NeuGc + C3H9N - O = HexNAc + dHex + NH3 - O (C14H26N2O8). The MS2 above favours
  the core-fucosylated amino-alditol reading for these species.
