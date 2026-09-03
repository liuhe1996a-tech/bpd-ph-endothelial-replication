# JGG clean-room reproduction

Date: 3 September 2026  
Platform tested: Windows 11 (AMD64)  
Python: 3.12.13

The software archive was extracted into a newly created directory that did not contain the authors' working project. Dependencies were supplied from an isolated project environment according to `requirements.txt`.

Command:

```bash
python run_jgg_release.py --skip-heavy --skip-render --data-archive Supplementary_Data_JGG.zip
```

Result: **PASS**.

- 118 release assets were verified through `release_restore_manifest.tsv` and restored to documented project-relative paths.
- Two-cohort benchmark validation passed 20 of 20 checks.
- Exact GSE243129 animal-sensitivity validation passed 12 of 12 checks.
- Both cohorts contained exactly 1,800 outer-training-fold HVGs per fold.
- All stochastic methods used ten fixed seeds per fold.
- Fitting and validation animals had zero overlap.
- Each cohort contained 500 shared biological-bootstrap identifiers.
- Benchmark figures, numerical source tables, the executed notebook and the technical report were regenerated.

The lightweight test consumes frozen released model outputs and does not claim to refit neural models. Cross-platform reproducibility is defined by successful checksum restoration and numerical agreement within the released validation tolerances; line endings, archive metadata and floating-point terminal digits can differ across operating systems.

A public log is generated automatically at the end of the command by `scripts/59_write_jgg_public_reproduction_log.py`. It reports only the relative command, software versions and PASS summaries; author-specific absolute paths are not written.

The lightweight route does not invoke R. The R version reported in the public log describes the frozen upstream assets, whereas the execution-platform line describes the system used for the lightweight replay.
