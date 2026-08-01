# BPD-PH endothelial replication and cross-species projection

Reproducible code for the manuscript **“Replicated capillary endothelial
responses to neonatal hyperoxia with exploratory projection to human
BPD-associated pulmonary hypertension.”**

The workflow integrates GSE216046 mouse bulk RNA-seq, GSE151974 animal-level
endothelial pseudobulk replication, and an explicitly exploratory GSE275938
human endothelial projection. Biological replication is defined at the animal
or donor level. All public inputs have fixed URLs and SHA-256 checksums.

## Reproduce

1. Install Python 3.12 packages with `pip install -r requirements.txt`.
2. Install R 4.5.1 packages listed in `environment/R_packages.tsv`.
3. Run `python run_all.py` for a dry-run command plan.
4. Run `python run_all.py --execute --rscript /path/to/Rscript`.

Outputs are created under `04_processed/`, `07_results/`, `08_figures/`, and
`09_plotting_data/`. Raw inputs are downloaded under `02_raw/` and are not
tracked by Git.

## Reproducibility status

The frozen R5 workflow was executed from a clean output directory on
2026-07-31. Nine core result tables matched the manuscript-freeze results,
and all 26 figure QA checks passed. See `docs/R5_clean_reproduction_audit.json`.

## Citation and release

Citation metadata are provided in `CITATION.cff`. The analysis code corresponding
to the submitted manuscript is publicly available in this repository. A permanent
archival DOI will be added upon publication.
