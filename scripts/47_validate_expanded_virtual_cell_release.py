"""Validate the R10 two-cohort, seven-method perturbation benchmark release."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


METHODS = {
    "gene_linear", "pca_latent", "vae_latent", "cvae_counterfactual",
    "scgen_adapted", "cpa_adapted", "sinkhorn_ot",
}
SEEDED = {
    "vae_latent", "cvae_counterfactual", "scgen_adapted", "cpa_adapted", "sinkhorn_ot"
}
ENDPOINTS = {
    "primary_fold_hvg",
    "secondary_replicated_198_hvg_intersection",
    "secondary_signature_33_hvg_intersection",
}


def check_cohort(
    point_path: Path,
    bootstrap_path: Path,
    feature_path: Path,
    split_path: Path,
    heldout_ages: set[str],
    expected_bootstrap: int,
) -> list[dict[str, object]]:
    checks: list[dict[str, object]] = []
    point = pd.read_csv(point_path, sep="\t")
    bootstrap = pd.read_csv(bootstrap_path, sep="\t")
    features = pd.read_csv(feature_path, sep="\t")
    splits = pd.read_csv(split_path, sep="\t")

    def add(name: str, passed: bool, detail: object) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})

    observed_methods = set(point.loc[point["cell_type"].isin(["Cap", "Cap-a"]), "method"])
    add("all seven reportable methods present", METHODS <= observed_methods, sorted(observed_methods))
    add("all three endpoints present", ENDPOINTS <= set(point["endpoint"]), sorted(point["endpoint"].unique()))
    add("expected held-out ages", set(point["heldout_age"]) == heldout_ages, sorted(point["heldout_age"].unique()))
    add(
        "exactly 1800 fold-only HVGs",
        features.groupby("heldout_age")["primary_fold_hvg"].sum().eq(1800).all(),
        features.groupby("heldout_age")["primary_fold_hvg"].sum().to_dict(),
    )
    add(
        "feature matrix contains no non-HVG entries",
        features.groupby("heldout_age").size().eq(1800).all()
        and features["primary_fold_hvg"].astype(bool).all(),
        features.groupby("heldout_age").size().to_dict(),
    )
    seed_counts = (
        point.loc[point["method"].isin(SEEDED)]
        .groupby(["heldout_age", "method"], observed=True)["seed"].nunique()
    )
    add("ten seeds for every stochastic method and fold", seed_counts.eq(10).all(), seed_counts.to_dict())
    overlap = splits.groupby(["heldout_age", "seed", "animal_id"], observed=True)["role"].nunique()
    overlapping_animal_keys = int((overlap > 1).sum())
    add(
        "no fit-validation animal overlap",
        overlapping_animal_keys == 0,
        f"overlapping animal keys={overlapping_animal_keys}",
    )
    add(
        "500 shared biological bootstrap identifiers",
        bootstrap["bootstrap"].nunique() == expected_bootstrap
        and bootstrap.groupby(["heldout_age", "cell_type", "method", "seed", "endpoint"], observed=True)["bootstrap"].nunique().eq(expected_bootstrap).all(),
        int(bootstrap["bootstrap"].nunique()),
    )
    point_key = ["heldout_age", "cell_type", "method", "seed", "endpoint"]
    boot_key = point_key + ["bootstrap"]
    add("point metric keys unique", not point.duplicated(point_key).any(), int(point.duplicated(point_key).sum()))
    add("bootstrap metric keys unique", not bootstrap.duplicated(boot_key).any(), int(bootstrap.duplicated(boot_key).sum()))
    return checks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    checks: list[dict[str, object]] = []
    checks.extend(
        {**row, "cohort": "GSE151974"}
        for row in check_cohort(
            args.benchmark_dir / "virtual_cell_point_metrics.tsv",
            args.benchmark_dir / "virtual_cell_strictly_paired_bootstrap.tsv.gz",
            args.benchmark_dir / "virtual_cell_fold_feature_manifest.tsv",
            args.benchmark_dir / "virtual_cell_grouped_validation_splits.tsv",
            {"P3", "P7", "P14"}, 500,
        )
    )
    checks.extend(
        {**row, "cohort": "GSE243129"}
        for row in check_cohort(
            args.benchmark_dir / "GSE243129_point_metrics.tsv",
            args.benchmark_dir / "GSE243129_strictly_paired_bootstrap.tsv.gz",
            args.benchmark_dir / "GSE243129_fold_feature_manifest.tsv",
            args.benchmark_dir / "GSE243129_grouped_validation_splits.tsv",
            {"P7", "P14"}, 500,
        )
    )
    frame = pd.DataFrame(checks)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output.with_suffix(".tsv"), sep="\t", index=False)
    report = {
        "analysis_version": "R10",
        "n_checks": int(len(frame)),
        "n_passed": int(frame["passed"].sum()),
        "all_passed": bool(frame["passed"].all()),
        "checks": checks,
    }
    args.output.with_suffix(".json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(frame.to_string(index=False))
    if not report["all_passed"]:
        raise SystemExit("R10 validation failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
