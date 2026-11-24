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

## Limits

* MS intensity is not molar amount; the amine ionizes differently from the 2-AB amide and the
  form does not fluoresce, so its molar share is not known from these data.
* MS2 readout at low resolution; assignments rest on +/-0.3 Da windows, with the in-run 2-AB
  control as the specificity check.
* A literature search beyond one round is needed before any claim that this is unreported for
  PNGase F + reductive amination workflows.
* Exact isomer: NeuGc + C3H9N - O = HexNAc + dHex + NH3 - O (C14H26N2O8). The MS2 above favours
  the core-fucosylated amino-alditol reading for these species.
