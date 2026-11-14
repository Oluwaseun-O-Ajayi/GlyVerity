# GlyVerity

**Artifact-aware, FDR-controlled glycan assignment for released-glycan LC-MS.**

Glycomics profiling reports glycan compositions without an error rate, and many detected
features are isotopes, adducts, label side-products or in-source fragments rather than
distinct glycans. GlyVerity aims to:

- **Clean:** group and annotate features explained by isotopes, adducts, label
  side-products and in-source fragments
- **Assign:** match compositions with target-decoy false discovery rate (q-values)
- **Report:** state the specificity of each assignment (composition, topology, linkage)
  with GlyTouCan identifiers and MIRAGE-ready output

GlyVerity reads mzML, builds on [pyOpenMS](https://pyopenms.readthedocs.io), and extends
[LC-MS Data Processor](https://doi.org/10.5281/zenodo.21652568).

## Status

In development (Phase 1). Not yet ready for use. Progress is logged in [DEVLOG.md](DEVLOG.md).

## Install (development)

```bash
git clone https://github.com/Oluwaseun-O-Ajayi/GlyVerity.git
cd GlyVerity
python -m pip install -e ".[dev]"      # add ",ms" for pyOpenMS
pytest -q
```

## Documentation

* [Architecture](docs/architecture.md)
* [Phase 1: prior art](docs/phase1/prior_art.md)
* [Phase 1: claim verification](docs/phase1/claims_verification.md)
* [Phase 1: go/no-go](docs/phase1/go_no_go.md)
* [Dataset catalog](data/catalog/README.md)

## Author

Oluwaseun O. Ajayi, Department of Chemistry, University of Georgia
ORCID: [0000-0003-0040-7217](https://orcid.org/0000-0003-0040-7217)

## License

MIT
