#!/usr/bin/env python3
"""Run the R10 two-cohort benchmark and publication-reporting workflow.

The lightweight path consumes released, checksum-verified outputs and rebuilds
the benchmark figures, manuscript files, and validation reports.  It does not
claim to refit neural models.  The full path refits all models when raw inputs
and R are available.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROJECT_PYTHON = ROOT / "python_packages"
SEEDS = [str(seed) for seed in range(20260817, 20260827)]


def existing(candidates) -> str | None:
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return str(candidate)
    return None


def resolve_rscript() -> str:
    candidate = existing(
        [
            os.environ.get("RSCRIPT"),
            shutil.which("Rscript"),
            Path(os.environ.get("ProgramFiles", "")) / "R" / "R-4.5.1" / "bin" / "Rscript.exe",
            Path(os.environ.get("ProgramFiles", "")) / "R" / "R-4.4.3" / "bin" / "Rscript.exe",
        ]
    )
    if candidate:
        return candidate
    raise FileNotFoundError("Rscript was not found; set RSCRIPT to the installed executable.")


def registry_libreoffice_candidates() -> list[str]:
    if os.name != "nt":
        return []
    candidates: list[str] = []
    try:
        import winreg

        keys = (
            (winreg.HKEY_LOCAL_MACHINE, r"SOFTWARE\LibreOffice\UNO\InstallPath"),
            (winreg.HKEY_LOCAL_MACHINE, r"SOFTWARE\WOW6432Node\LibreOffice\UNO\InstallPath"),
            (winreg.HKEY_CURRENT_USER, r"SOFTWARE\LibreOffice\UNO\InstallPath"),
        )
        for hive, key_name in keys:
            try:
                with winreg.OpenKey(hive, key_name) as key:
                    value, _ = winreg.QueryValueEx(key, "")
                    install = Path(value)
                    candidates.extend([str(install / "soffice.com"), str(install / "soffice.exe")])
            except OSError:
                continue
    except ImportError:
        pass
    return candidates


def resolve_libreoffice() -> str:
    roots = [
        os.environ.get("LIBREOFFICE_HOME"),
        os.environ.get("ProgramFiles"),
        os.environ.get("ProgramFiles(x86)"),
        os.environ.get("LOCALAPPDATA"),
    ]
    candidates: list[object] = [
        os.environ.get("SOFFICE"),
        shutil.which("soffice.com"),
        shutil.which("soffice.exe"),
    ]
    candidates.extend(registry_libreoffice_candidates())
    for root in roots:
        if not root:
            continue
        base = Path(root)
        candidates.extend(
            [
                base / "LibreOffice" / "program" / "soffice.com",
                base / "LibreOffice" / "program" / "soffice.exe",
                base / "Programs" / "LibreOffice" / "program" / "soffice.com",
                base / "Programs" / "LibreOffice" / "program" / "soffice.exe",
            ]
        )
    candidate = existing(candidates)
    if candidate:
        return candidate
    raise FileNotFoundError(
        "LibreOffice was not found through SOFFICE, LIBREOFFICE_HOME, PATH, "
        "the Windows registry, or common install roots. Set SOFFICE explicitly."
    )


def environment() -> dict[str, str]:
    env = dict(os.environ)
    if PROJECT_PYTHON.exists():
        env["PYTHONPATH"] = str(PROJECT_PYTHON) + os.pathsep + env.get("PYTHONPATH", "")
    return env


def run(label: str, command: list[str], dry_run: bool) -> None:
    print(f"[{label}] {' '.join(command)}", flush=True)
    if not dry_run:
        subprocess.run(command, cwd=ROOT, check=True, env=environment())


def ensure_lightweight_assets(py: str, data_archive: Path | None, dry_run: bool) -> None:
    if data_archive:
        run(
            "Restore checksum-verified release assets",
            [
                py,
                str(ROOT / "scripts" / "56_restore_release_assets.py"),
                "--root",
                str(ROOT),
                "--data-archive",
                str(data_archive.resolve()),
            ],
            dry_run,
        )
    elif (ROOT / "release_restore_manifest.tsv").exists():
        run(
            "Restore checksum-verified release assets",
            [
                py,
                str(ROOT / "scripts" / "56_restore_release_assets.py"),
                "--root",
                str(ROOT),
                "--source-root",
                str(ROOT),
            ],
            dry_run,
        )
    if dry_run:
        return
    required = [
        ROOT / "results" / "virtual_cell_benchmark_r10" / "virtual_cell_point_metrics.tsv",
        ROOT / "results" / "virtual_cell_benchmark_r10" / "GSE243129_point_metrics.tsv",
        ROOT / "results" / "virtual_cell_r10_sensitivity" / "GSE243129_exact_animal_sensitivity_audit.json",
        ROOT / "results" / "LungMAP_human_validation" / "LungMAP_signature_models.tsv",
        ROOT / "results" / "LungMAP_human_validation_corrected" / "LungMAP_primary_signature_equivalence_tests.tsv",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Lightweight release assets are missing. Extract Supplementary_Data_R10.zip "
            "or pass --data-archive. Missing: " + "; ".join(missing)
        )


def render_documents(dry_run: bool) -> None:
    source = ROOT / "manuscript" / "corrected_r10"
    destination = source / "rendered_pdf"
    if dry_run:
        print(f"[Render documents] LibreOffice -> {destination}")
        return
    soffice = resolve_libreoffice()
    temp_parent = Path(os.environ.get("R10_TEMP_DIR", tempfile.gettempdir())).resolve()
    temp_parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="r10_lo_work_", dir=temp_parent))
    output = Path(tempfile.mkdtemp(prefix="r10_lo_pdf_", dir=temp_parent))
    profile = Path(tempfile.mkdtemp(prefix="r10_lo_profile_", dir=temp_parent))
    try:
        documents = sorted(source.glob("*.docx"))
        if not documents:
            raise FileNotFoundError(f"No DOCX files found in {source}")
        for document in documents:
            shutil.copy2(document, staging / document.name)
        profile_uri = profile.as_uri()
        for document in sorted(staging.glob("*.docx")):
            libreoffice_env = dict(os.environ)
            libreoffice_env.pop("PYTHONHOME", None)
            libreoffice_env.pop("PYTHONPATH", None)
            subprocess.run(
                [
                    soffice,
                    f"-env:UserInstallation={profile_uri}",
                    "--headless",
                    "--convert-to",
                    "pdf:writer_pdf_Export",
                    "--outdir",
                    str(output),
                    str(document),
                ],
                check=True,
                cwd=staging,
                env=libreoffice_env,
            )
        destination.mkdir(parents=True, exist_ok=True)
        for pdf in output.glob("*.pdf"):
            shutil.copy2(pdf, destination / pdf.name)
        print(f"[Render documents] {len(list(destination.glob('*.pdf')))} PDFs written with {soffice}")
    finally:
        shutil.rmtree(staging, ignore_errors=True)
        shutil.rmtree(output, ignore_errors=True)
        shutil.rmtree(profile, ignore_errors=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-heavy", action="store_true")
    parser.add_argument("--skip-render", action="store_true")
    parser.add_argument("--data-archive", type=Path)
    args = parser.parse_args()
    py = sys.executable

    heavy: list[tuple[str, list[str]]] = []
    if not args.skip_heavy:
        rscript = "Rscript" if args.dry_run else resolve_rscript()
        heavy = [
            (
                "Export GSE243129 benchmark inputs",
                [rscript, str(ROOT / "scripts" / "43_export_GSE243129_virtual_cell_inputs.R"), str(ROOT)],
            ),
            (
                "Fit primary VAE and conditional-VAE baselines",
                [
                    py, str(ROOT / "scripts" / "27_virtual_cell_leakage_free_benchmark.py"),
                    "--matrix", str(ROOT / "raw/GSE151974_processed/GSE151974_endothelial_counts_cells_by_genes.npz"),
                    "--metadata", str(ROOT / "raw/GSE151974_processed/GSE151974_endothelial_cell_metadata.tsv.gz"),
                    "--genes", str(ROOT / "raw/GSE151974_processed/GSE151974_endothelial_genes.tsv"),
                    "--signature", str(ROOT / "raw/frozen_prior_results/Supplementary_Data_3_mouse_replicated_198_genes.tsv"),
                    "--output-dir", str(ROOT / "results/virtual_cell_hvg_only_r10"),
                    "--max-epochs", "60", "--bootstrap", "500", "--seeds", *SEEDS,
                ],
            ),
            (
                "Fit primary scGen-style, CPA-style and Sinkhorn models",
                [
                    py, str(ROOT / "scripts" / "44_run_expanded_virtual_cell_models.py"),
                    "--cohort", "GSE151974",
                    "--matrix", str(ROOT / "raw/GSE151974_processed/GSE151974_endothelial_counts_cells_by_genes.npz"),
                    "--metadata", str(ROOT / "raw/GSE151974_processed/GSE151974_endothelial_cell_metadata.tsv.gz"),
                    "--genes", str(ROOT / "raw/GSE151974_processed/GSE151974_endothelial_genes.tsv"),
                    "--signature", str(ROOT / "raw/frozen_prior_results/Supplementary_Data_3_mouse_replicated_198_genes.tsv"),
                    "--output-dir", str(ROOT / "results/virtual_cell_GSE151974_expanded_r10"),
                    "--heldout-ages", "P3", "P7", "P14", "--cell-types", "Cap", "Cap-a",
                    "--reuse-r8-dir", str(ROOT / "results/virtual_cell_hvg_only_r10"),
                    "--max-epochs", "60", "--bootstrap", "500", "--seeds", *SEEDS,
                ],
            ),
            (
                "Fit all seven methods in GSE243129",
                [
                    py, str(ROOT / "scripts" / "44_run_expanded_virtual_cell_models.py"),
                    "--cohort", "GSE243129",
                    "--matrix", str(ROOT / "raw/GSE243129_virtual_cell/GSE243129_WT_capillary_counts_cells_by_genes.mtx.gz"),
                    "--metadata", str(ROOT / "raw/GSE243129_virtual_cell/GSE243129_WT_capillary_cell_metadata.tsv.gz"),
                    "--genes", str(ROOT / "raw/GSE243129_virtual_cell/GSE243129_WT_capillary_genes.tsv"),
                    "--signature", str(ROOT / "raw/frozen_prior_results/Supplementary_Data_3_mouse_replicated_198_genes.tsv"),
                    "--output-dir", str(ROOT / "results/virtual_cell_GSE243129_expanded_r10"),
                    "--heldout-ages", "P7", "P14", "--cell-types", "Cap", "Cap-a",
                    "--include-existing-models", "--max-epochs", "60", "--bootstrap", "500", "--seeds", *SEEDS,
                ],
            ),
        ]
    else:
        ensure_lightweight_assets(py, args.data_archive, args.dry_run)

    integration: list[tuple[str, list[str]]] = []
    if not args.skip_heavy:
        integration = [
            (
                "Integrate independently fitted cohorts",
                [
                    py, str(ROOT / "scripts" / "45_integrate_expanded_virtual_cell_benchmark.py"),
                    "--r8-dir", str(ROOT / "results/virtual_cell_hvg_only_r10"),
                    "--primary-addition-dir", str(ROOT / "results/virtual_cell_GSE151974_expanded_r10"),
                    "--secondary-dir", str(ROOT / "results/virtual_cell_GSE243129_expanded_r10"),
                    "--output-dir", str(ROOT / "results/virtual_cell_benchmark_r10"),
                ],
            ),
            (
                "External GSE266988 calibration",
                [
                    py, str(ROOT / "scripts" / "31_external_calibration_corrected_benchmark.py"),
                    "--predictions", str(ROOT / "results/virtual_cell_benchmark_r10/virtual_cell_gene_predictions.tsv.gz"),
                    "--external-hyperoxia", str(ROOT / "raw/GSE266988/jci.insight.182880.sdt2.xlsx"),
                    "--output-dir", str(ROOT / "results/virtual_cell_external_calibration_r10"),
                ],
            ),
            (
                "Exact GSE243129 animal-subset analysis",
                [py, str(ROOT / "scripts" / "54_GSE243129_exact_animal_sensitivity.py"), "--root", str(ROOT)],
            ),
        ]

    reporting = [
        (
            "Exact GSE243129 animal-sensitivity validation",
            [py, str(ROOT / "scripts" / "57_validate_r10_sensitivity.py"), "--root", str(ROOT)],
        ),
        (
            "R10 benchmark figures and source data",
            [
                py, str(ROOT / "scripts" / "46_plot_expanded_virtual_cell_benchmark.py"),
                "--benchmark-dir", str(ROOT / "results/virtual_cell_benchmark_r10"),
                "--figure-dir", str(ROOT / "figures/r10"),
                "--source-dir", str(ROOT / "source_data/r10"),
                "--sensitivity-dir", str(ROOT / "results/virtual_cell_r10_sensitivity"),
            ],
        ),
        (
            "Two-cohort benchmark validation",
            [
                py, str(ROOT / "scripts" / "47_validate_expanded_virtual_cell_release.py"),
                "--benchmark-dir", str(ROOT / "results/virtual_cell_benchmark_r10"),
                "--output", str(ROOT / "logs/r10_expanded_virtual_cell_validation"),
            ],
        ),
        (
            "Executed R10 audit notebook and technical report",
            [py, str(ROOT / "scripts" / "55_build_r10_benchmark_notebook_and_report.py"), "--root", str(ROOT)],
        ),
    ]
    for label, command in heavy + integration + reporting:
        run(label, command, args.dry_run)
    if not args.skip_render:
        render_documents(args.dry_run)
    run(
        "Path-sanitized public reproduction log",
        [py, str(ROOT / "scripts" / "59_write_jgg_public_reproduction_log.py"), "--root", str(ROOT)],
        args.dry_run,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
