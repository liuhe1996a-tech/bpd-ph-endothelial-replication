# Neonatal lung endothelial virtual-cell benchmark

Reproducible code for the manuscript **Simple baselines match variational autoencoders in a cross-dataset neonatal lung endothelial benchmark**.

The repository contains two linked workflows:

1. cross-dataset mouse endothelial discovery and replication with an exploratory human projection;
2. held-out virtual-cell benchmarking, independent capillary calibration, regulon analysis, p53 perturbation calibration and animal-balanced ligand-receptor analysis.

## Data

All inputs are public. The main accessions are GSE216046, GSE151974 and GSE275938. Large expression matrices are downloaded locally and are not committed to Git. Fixed URLs and SHA-256 checksums are recorded in `config/download_manifest.tsv` and the accompanying manifests.

## Environment

- Python 3.12: `pip install -r requirements.txt` and `pip install -r AI_upgrade/requirements_ai.txt`
- R 4.5.1: install versions listed in `environment/R_packages.tsv`

## Reproduce

Copy `AI_upgrade/workflow_config.example.json` to a local configuration file and update project-relative paths. Inspect both command plans with:

```text
python run_bib_workflow.py --ai-config AI_upgrade/workflow_config.example.json
```

Execute the complete workflow with:

```text
python run_bib_workflow.py --execute --rscript /path/to/Rscript --ai-config /path/to/workflow_config.json
```

The base workflow can be skipped when its frozen outputs already exist:

```text
python run_bib_workflow.py --execute --skip-base --ai-config /path/to/workflow_config.json
```

The virtual-cell models use seed 20260807. Other fixed seeds are documented in the workflow README files. Predictions are evaluated at the animal level, and the VAE results are reported even where they do not outperform simpler baselines.

## Release

The code corresponding to the Briefings in Bioinformatics submission is tagged `v1.1.0-bib-submission`. A permanent archival DOI can be added through Zenodo after the authors enable repository archiving.
