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
* Most public data are PGC/reduced or C18/permethylated, not labelled HILIC. Scope decision
  needed (see `docs/phase1/go_no_go.md`).

**Next**
* Download one glycan raw file each from PXD077477 and PXD078897; check for a fluorescence channel.
* Browse GlycoPOST and MassIVE; add rows to the catalog.
* Count fluorescence-paired and standards-based ground truth sets (go/no-go condition 1).
* Finish prior-art items marked "pending" and "not yet checked".
