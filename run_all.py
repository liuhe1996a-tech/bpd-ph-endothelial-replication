"""Unified entry point for the frozen discovery-replication-projection workflow."""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true", help="Run commands; default is dry-run.")
    parser.add_argument("--rscript", default=os.environ.get("RSCRIPT", "Rscript"))
    parser.add_argument("--from-stage", default="download")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    py = sys.executable
    r = args.rscript
    scripts = root / "scripts"
    mouse = scripts / "04_mouse_endothelium"
    cross = scripts / "05_cross_species"
    pathways = scripts / "06_pathways"
    reporting = scripts / "06_reporting"
    exact = root / "07_results" / "mouse_endothelium_exact"
    processed_mouse = root / "04_processed" / "GSE151974"
    h5_root = root / "02_raw" / "GSE275938" / "GSE275938_RAW"
    steps = [
        ("download", [py, str(root / "download_required_data.py"), "--project-root", str(root)]),
        ("msigdb", [py, str(scripts / "00_setup" / "download_mouse_msigdb.py")]),
        ("mouse_bulk", [r, str(mouse / "GSE216046_DESeq2.R"), str(root / "02_raw/GSE216046/GSE216046_gene_count.csv.gz"), str(exact)]),
        ("mouse_sc_build", [py, str(mouse / "build_GSE151974_endothelial_pseudobulk.py")]),
        ("mouse_sc_model", [r, str(mouse / "GSE151974_edgeR_pseudobulk.R"), str(processed_mouse / "GSE151974_endothelial_animal_pseudobulk_counts.tsv.gz"), str(processed_mouse / "GSE151974_endothelial_animal_pseudobulk_metadata.tsv"), str(exact)]),
        ("composition", [py, str(mouse / "GSE151974_composition_permutation_sensitivity.py")]),
        ("mouse_integrate", [py, str(mouse / "integrate_exact_mouse_results.py")]),
        ("pathways", [r, str(pathways / "mouse_pathway_fgsea.R")]),
        ("pathway_postprocess", [py, str(pathways / "postprocess_mouse_pathway_replication.py")]),
        ("pathway_revision", [py, str(pathways / "postprocess_mouse_pathway_revision.py")]),
        ("orthology", [py, str(cross / "build_ncbi_mouse_human_ortholog_mapping.py")]),
        ("human_build", [py, str(cross / "build_GSE275938_endothelial_subtype_pseudobulk.py")]),
        ("human_tmm", [r, str(cross / "GSE275938_endothelial_subtype_TMM.R")]),
        ("human_sensitivity", [py, str(cross / "GSE275938_human_subtype_sensitivity.py")]),
        ("human_technical_audit", [py, str(cross / "GSE275938_donor_technical_audit.py")]),
        ("matched_null", [py, str(cross / "GSE275938_expression_matched_null.py")]),
        ("human_postprocess", [py, str(cross / "postprocess_GSE275938_human_results.py")]),
        ("composition_restricted", [r, str(cross / "GSE275938_composition_restricted_TMM.R")]),
        ("signature_sensitivity", [py, str(cross / "signature_and_composition_sensitivity.py")]),
        ("h5_projection_audit", [py, str(cross / "build_GSE275938_projection_audit.py"), "--project-root", str(root), "--h5-root", str(h5_root)]),
        ("human_figure", [r, str(reporting / "plot_GSE275938_human_subtype_sensitivity.R")]),
        ("main_figures", [r, str(reporting / "build_manuscript_main_figures.R")]),
        ("validate_figures", [r, str(reporting / "validate_manuscript_main_figures.R")]),
    ]
    stage_names = [name for name, _ in steps]
    if args.from_stage not in stage_names:
        raise SystemExit(f"Unknown stage {args.from_stage!r}; choose from {stage_names}")
    environment = os.environ.copy()
    environment["BPD_PH_PROJECT_ROOT"] = str(root)
    environment["GSE151974_RAW_ROOT"] = str(root / "02_raw" / "GSE151974")
    environment["GSE275938_RAW_ROOT"] = str(root / "02_raw" / "GSE275938")
    environment["GSE275938_H5_ROOT"] = str(h5_root)
    start = stage_names.index(args.from_stage)
    for name, command in steps[start:]:
        print(f"[{name}] {' '.join(command)}")
        if args.execute:
            if command[0] == r and shutil.which(r) is None and not Path(r).exists():
                raise RuntimeError("Rscript was not found; pass --rscript with its path.")
            subprocess.run(command, cwd=root, env=environment, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
