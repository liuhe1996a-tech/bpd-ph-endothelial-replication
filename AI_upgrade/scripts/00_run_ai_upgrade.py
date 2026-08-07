"""Run the benchmark-guided AI upgrade from project-relative inputs.

The workflow deliberately benchmarks virtual-cell models against simple
baselines and treats all perturbation outputs as predictions rather than as
experimental evidence.  Paths are read from a JSON configuration file and
resolved relative to ``project_root``.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


STEP_ORDER = [
    "matrix",
    "virtual_cell",
    "regulon",
    "ligand_receptor",
    "external_calibration",
    "paired_comparison",
    "figures",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--from-step", choices=STEP_ORDER, default=STEP_ORDER[0])
    parser.add_argument("--to-step", choices=STEP_ORDER, default=STEP_ORDER[-1])
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def command_path(value: str, project_root: Path) -> str:
    path = Path(value)
    if not path.is_absolute():
        path = project_root / path
    return str(path.resolve())


def main() -> int:
    args = parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    root_value = Path(config.get("project_root", "."))
    project_root = (
        root_value.resolve()
        if root_value.is_absolute()
        else (config_path.parent / root_value).resolve()
    )
    paths = config["paths"]
    parameters = config.get("parameters", {})
    scripts = Path(__file__).resolve().parent

    p = lambda key: command_path(paths[key], project_root)
    py = sys.executable
    virtual_args = parameters.get("virtual_cell", {})
    lr_args = parameters.get("ligand_receptor", {})

    commands: dict[str, list[str]] = {
        "matrix": [
            py, str(scripts / "01_build_mouse_endothelial_matrix.py"),
            "--matrix", p("raw_count_matrix"),
            "--metadata", p("raw_cell_metadata"),
            "--output-dir", p("endothelial_matrix_dir"),
        ],
        "virtual_cell": [
            py, str(scripts / "02_virtual_cell_age_holdout.py"),
            "--matrix", str(Path(p("endothelial_matrix_dir")) / "GSE151974_endothelial_counts_cells_by_genes.npz"),
            "--metadata", str(Path(p("endothelial_matrix_dir")) / "GSE151974_endothelial_cell_metadata.tsv.gz"),
            "--genes", str(Path(p("endothelial_matrix_dir")) / "GSE151974_endothelial_genes.tsv"),
            "--signature", p("mouse_signature"),
            "--output-dir", p("virtual_cell_dir"),
            "--heldout-ages", *[str(x) for x in virtual_args.get("heldout_ages", ["P3", "P7", "P14"])],
            "--n-hvg", str(virtual_args.get("n_hvg", 1800)),
            "--latent-dim", str(virtual_args.get("latent_dim", 32)),
            "--max-epochs", str(virtual_args.get("max_epochs", 80)),
            "--batch-size", str(virtual_args.get("batch_size", 256)),
            "--beta", str(virtual_args.get("beta", 0.0001)),
            "--bootstrap", str(virtual_args.get("bootstrap", 500)),
            "--seed", str(virtual_args.get("seed", 20260807)),
        ],
        "regulon": [
            py, str(scripts / "03_regulon_virtual_perturbation.py"),
            "--matrix", str(Path(p("endothelial_matrix_dir")) / "GSE151974_endothelial_counts_cells_by_genes.npz"),
            "--metadata", str(Path(p("endothelial_matrix_dir")) / "GSE151974_endothelial_cell_metadata.tsv.gz"),
            "--genes", str(Path(p("endothelial_matrix_dir")) / "GSE151974_endothelial_genes.tsv"),
            "--network", p("collectri_network"),
            "--signature", p("mouse_signature"),
            "--jci-hyperoxia", p("jci_hyperoxia"),
            "--jci-p53-null", p("jci_p53_null"),
            "--jci-p53-ec", p("jci_p53_endothelial"),
            "--output-dir", p("regulon_dir"),
        ],
        "ligand_receptor": [
            py, str(scripts / "04_animal_balanced_ligand_receptor.py"),
            "--matrix", p("raw_count_matrix"),
            "--metadata", p("raw_cell_metadata"),
            "--lr-resource", p("omnipath_lr_resource"),
            "--external-cellchat", p("external_cellchat_results"),
            "--output-dir", p("ligand_receptor_dir"),
            "--chunk-size", str(lr_args.get("chunk_size", 100)),
            "--min-cells", str(lr_args.get("min_cells", 20)),
            "--min-logcpm", str(lr_args.get("min_logcpm", 0.6931471805599453)),
        ],
        "external_calibration": [
            py, str(scripts / "05_external_virtual_cell_calibration.py"),
            "--predictions", str(Path(p("virtual_cell_dir")) / "virtual_cell_all_gene_predictions.tsv.gz"),
            "--external-hyperoxia", p("jci_hyperoxia"),
            "--internal-metrics", str(Path(p("virtual_cell_dir")) / "virtual_cell_benchmark_metrics.tsv"),
            "--output-dir", p("external_calibration_dir"),
        ],
        "paired_comparison": [
            py, str(scripts / "06_virtual_cell_paired_comparison.py"),
            "--bootstrap", str(Path(p("virtual_cell_dir")) / "virtual_cell_benchmark_bootstrap.tsv.gz"),
            "--output-dir", p("paired_comparison_dir"),
        ],
        "figures": [
            py, str(scripts / "07_build_ai_upgrade_figures.py"),
            "--results-root", p("results_root"),
            "--output-dir", p("figures_dir"),
        ],
    }

    start = STEP_ORDER.index(args.from_step)
    stop = STEP_ORDER.index(args.to_step)
    if start > stop:
        raise SystemExit("--from-step must precede --to-step")

    for step in STEP_ORDER[start : stop + 1]:
        command = commands[step]
        print(f"[{step}]", " ".join(command), flush=True)
        if not args.dry_run:
            subprocess.run(command, check=True, cwd=project_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
