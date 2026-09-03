# Multi-cohort neonatal lung endothelial functional genomics

This repository is the versioned analysis companion to the manuscript:

> **Cross-cohort single-cell genomics of neonatal hyperoxia reveals reproducible endothelial responses and cohort-dependent model rankings**

The submission version is frozen as tag **`v1.0.0-jgg`**.

Archival DOI: **https://doi.org/10.5281/zenodo.22263929**

The repository separates three claims that are often conflated in single-cell computational studies: within-species biological replication, held-out perturbation prediction, and cross-species transfer.

## Study design

1. **Mouse discovery and internal replication** — GSE216046 and GSE151974 define a fixed endothelial response using animals, not cells, as statistical units.
2. **External mouse replication** — GSE243129, GSE209664 and GSE230672 test the locked genes and pathways without redefining them.
3. **Perturbation-prediction benchmark** — seven methods share fold-specific 1,800-gene feature sets, animal-grouped validation, ten fixed training seeds and joint animal bootstrap resampling in GSE151974 and the independently trained GSE243129 cohort.
4. **Human transfer test** — a versioned LungMAP object provides the 21-donor primary test; GSE275938 is retained only as a technically confounded four-donor sensitivity analysis.

The seven reportable methods are gene-space shift, PCA shift, VAE, conditional VAE, protocol-adapted scGen, protocol-adapted CPA and entropy-regularized Sinkhorn optimal transport. The scGen-style and CPA-style entries preserve the defining perturbation operators under a common split and tuning budget; they are not executions of the official legacy packages. The Sinkhorn entry is an explicit optimal-transport baseline in fold-trained PCA space and is not CellOT.

## Main conclusion

The mouse endothelial response reproduces across independent cohorts, but perturbation-prediction rankings depend on cohort and endpoint. The independent human atlas does not support a large, broad directional mouse-to-human transfer. Equivalence testing excludes effects of at least 0.5 standard deviations for the primary strict signature but does not exclude effects near 0.3 standard deviations.

## Repository layout

- `scripts/` — benchmark, sensitivity, restoration, validation and reporting code
- `config/release_restore_manifest.tsv` — expected project-relative paths and SHA-256 values
- `results/` — compact, machine-readable benchmark and sensitivity summaries
- `source_data/` — numerical source tables included with the code snapshot
- `figures/jgg_submission/` — final JGG figures in PDF and PNG
- `manuscript/jgg/` — final JGG manuscript and supplementary-information sources
- `reports/r10/` — executed audit notebook and technical report for the frozen scientific analysis
- `docs/JGG_CLEAN_REPRODUCTION.md` — clean-room protocol and verified test summary
- `docs/JGG_CLEAN_REPRODUCTION_LOG.txt` — automatically generated, path-sanitized public run log

Large public inputs are not redistributed. Accession identifiers and the versioned LungMAP object are documented in the manuscript and release manifest.

## Lightweight reproduction

Install the versions in `requirements.txt`, place the companion archive `Supplementary_Data_JGG.zip` beside the repository, and run:

```bash
python run_jgg_release.py --skip-heavy --skip-render --data-archive Supplementary_Data_JGG.zip
```

The command verifies and restores 118 released assets, recreates the benchmark figures and source tables, regenerates the audit documents, and runs the two-cohort and exact-animal sensitivity validations. Omit `--skip-render` only when LibreOffice is available. The lightweight route consumes frozen model outputs; it does not claim to refit neural networks or reprocess the largest public matrices.

Full refitting requires the public raw inputs identified in the manuscript and substantially more compute. The feature locks are preserved across entry points: external mouse and human outcomes cannot redefine the 198-gene response, 33-gene multi-subtype signature, strict 25-gene cross-mouse subset or 162 mouse pathway set.

## Reproducibility safeguards

- every model matrix contains exactly 1,800 genes selected inside the outer training fold;
- no internally selected disease-set gene is forced into a primary model matrix;
- fitting and validation animals do not overlap;
- stochastic training uses ten fixed seeds per fold;
- seed variation is reported separately from biological bootstrap variation;
- all models in a cohort share the same folds and animal-resampling ledger;
- GSE243129 exact-subset and leave-one-animal analyses quantify influence without claiming additional biological replication;
- all analysis paths are project-relative and restored with SHA-256 verification.

## Software and licensing

- Python 3.12 or later; versions are listed in `requirements.txt`
- R 4.5 or later; versions are listed in `environment_R_packages.tsv`
- source code and executable workflow files: MIT License
- author-generated figures, source tables, derived result tables and documentation: CC BY 4.0

See `LICENSE` and `LICENSES.md` for details.

The fixed release is archived at https://doi.org/10.5281/zenodo.22263929. The GitHub tag `v1.0.0-jgg` and the Zenodo record identify the same code, source data and derived results.
