# Benchmark-guided AI upgrade

This directory reproduces the additional analyses used in the manuscript. The
workflow was designed to test whether a virtual-cell model adds predictive
value, not to assume that a deep model is superior.

## Scientific scope

1. Build a sparse mouse lung endothelial matrix from GSE151974.
2. Hold out each age-by-hyperoxia condition and compare linear, PCA, VAE and
   conditional-VAE predictions at the animal level.
3. Infer transcription-factor activities with a curated CollecTRI network and
   calibrate p53-related predictions against an independent genetic study.
4. Evaluate ligand-receptor expression routes using animal-balanced
   pseudobulks and external CellChatDB annotation.
5. Compare virtual-cell models by paired bootstrap and assemble figures.

The VAE and conditional VAE did not consistently outperform the simple
baselines. This negative benchmark is retained in the manuscript to prevent
overclaiming. Network perturbations and ligand-receptor routes are predictive
or expression-supported results; they are not evidence of physical binding or
causality.

## Run

Copy `workflow_config.example.json`, adjust only project-relative paths, then
run:

```text
python scripts/00_run_ai_upgrade.py --config workflow_config.json --dry-run
python scripts/00_run_ai_upgrade.py --config workflow_config.json
```

To resume a partial run:

```text
python scripts/00_run_ai_upgrade.py --config workflow_config.json --from-step regulon
```

The random seed is 20260807. The large public expression matrices are not
redistributed in this software archive; accession identifiers and input
checksums are documented in the manuscript and data manifests. The original
mouse and human differential-expression workflow is provided in the parent
archive and should be run before this extension when starting from downloaded
GEO files.

## Environment

Python and R package versions are listed in `requirements_ai.txt` and in the
Nature Portfolio Reporting Summary. The workflow accepts all paths through the
JSON configuration and contains no author-specific filesystem paths.
