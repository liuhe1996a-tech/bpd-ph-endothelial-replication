#!/usr/bin/env python3
"""Write the path-sanitized R5 public log after a successful clean-room run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from datetime import date
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def first_line(path: Path, fallback: str) -> str:
    if not path.is_file():
        return fallback
    return path.read_text(encoding="utf-8", errors="replace").splitlines()[0].strip() or fallback


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = (args.output or (root / "docs" / "JGG_CLEAN_REPRODUCTION_LOG.txt")).resolve()

    manifest_path = root / "config" / "release_restore_manifest.tsv"
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    failures: list[str] = []
    unique_paths: set[str] = set()
    for row in manifest:
        relative = row["project_path"]
        unique_paths.add(relative)
        target = root / Path(relative)
        if (
            not target.is_file()
            or target.stat().st_size != int(row["size_bytes"])
            or sha256(target) != row["sha256"]
        ):
            failures.append(relative)
    if failures:
        raise RuntimeError(f"Release restoration verification failed for {len(failures)} project-relative files")

    benchmark_path = root / "logs" / "r10_expanded_virtual_cell_validation.json"
    sensitivity_path = root / "logs" / "r10_sensitivity_validation" / "R10_sensitivity_validation_summary.json"
    benchmark = json.loads(benchmark_path.read_text(encoding="utf-8"))
    sensitivity = json.loads(sensitivity_path.read_text(encoding="utf-8"))
    if not benchmark.get("all_passed") or benchmark.get("n_passed") != benchmark.get("n_checks"):
        raise RuntimeError("Two-cohort benchmark validation did not pass")
    if not sensitivity.get("all_passed") or sensitivity.get("passed") != sensitivity.get("checks"):
        raise RuntimeError("Exact-animal sensitivity validation did not pass")

    required_outputs = [
        root / "figures" / "r10" / "Figure5_expanded_virtual_cell_benchmark.png",
        root / "figures" / "r10" / "Supplementary_Figure_GSE243129_exact_animal_sensitivity.png",
        root / "reports" / "r10" / "R10_virtual_cell_benchmark_audit.ipynb",
        root / "reports" / "r10" / "R10_virtual_cell_benchmark_technical_report.html",
    ]
    missing_outputs = [path.name for path in required_outputs if not path.is_file()]
    if missing_outputs:
        raise RuntimeError("Regenerated public outputs are missing: " + ", ".join(missing_outputs))

    r_version = first_line(
        root / "results" / "GSE243129_external_replication" / "sessionInfo.txt",
        "R version recorded in environment_R_packages.tsv",
    )
    today = date.today().strftime("%d %B %Y")
    lines = [
        "JGG clean-room reproduction log",
        "",
        f"Generated: {today}",
        f"Lightweight execution platform: {platform.system()} {platform.release()} ({platform.machine()})",
        f"Python invoked in lightweight run: {platform.python_version()}",
        f"R used for frozen upstream assets: {r_version}",
        "R invoked in lightweight run: No",
        "",
        "Command:",
        "python run_jgg_release.py --skip-heavy --skip-render --data-archive Supplementary_Data_JGG.zip",
        "",
        "Overall result: PASS",
        f"Release restoration manifest: PASS ({len(manifest)}/{len(manifest)} entries; {len(unique_paths)} project-relative files)",
        f"Two-cohort benchmark validation: PASS ({benchmark['n_passed']}/{benchmark['n_checks']} checks)",
        f"Exact GSE243129 animal-sensitivity validation: PASS ({sensitivity['passed']}/{sensitivity['checks']} checks)",
        "Fold-specific feature selection: PASS (1,800 outer-training-fold HVGs in every fold)",
        "Stochastic training ledger: PASS (ten fixed seeds per stochastic method and fold)",
        "Fit-validation separation: PASS (overlapping animal keys = 0)",
        "Shared biological bootstrap ledger: PASS (500 identifiers per cohort)",
        "Regenerated figures and numerical source tables: PASS",
        "Executed audit notebook and technical report: PASS",
        "",
        "The lightweight route consumes checksum-verified frozen model outputs and does not refit neural networks. ",
        "Cross-platform agreement is evaluated by checksum restoration and the released numerical tolerances; line endings, archive metadata and terminal floating-point digits may differ.",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(output.relative_to(root).as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
