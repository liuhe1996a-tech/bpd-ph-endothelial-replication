"""Integrate the frozen R8 benchmark with R9 methods and an external cohort.

The script never recomputes model predictions.  It verifies that the added
GSE151974 methods used the frozen R8 features, validation animals and joint
bootstrap ledger, then combines them with the original four methods.  The
independent GSE243129 benchmark is retained as a separate cohort and summarized
with the same endpoint and paired-bootstrap definitions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


METRICS = ("pearson_effect", "spearman_effect", "rmse_effect", "direction_accuracy")
CAPILLARY = ("Cap", "Cap-a")
DETERMINISTIC = ("identity", "gene_linear", "pca_latent")
SEEDED = (
    "vae_latent",
    "cvae_counterfactual",
    "scgen_adapted",
    "cpa_adapted",
    "sinkhorn_ot",
)
ADDED = ("scgen_adapted", "cpa_adapted", "sinkhorn_ot")
BASELINES = ("gene_linear", "pca_latent")


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def stable_frame(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    return frame.sort_values(columns, kind="stable").reset_index(drop=True)


def assert_frames_equal(left: pd.DataFrame, right: pd.DataFrame, columns: list[str], label: str) -> None:
    left = stable_frame(left[columns].copy(), columns)
    right = stable_frame(right[columns].copy(), columns)
    try:
        pd.testing.assert_frame_equal(left, right, check_dtype=False, check_like=False)
    except AssertionError as exc:
        raise AssertionError(f"Frozen-ledger mismatch for {label}: {exc}") from exc


def paired_comparisons(bootstrap: pd.DataFrame) -> pd.DataFrame:
    data = bootstrap.loc[bootstrap["cell_type"].isin(CAPILLARY)].copy()
    rows: list[dict[str, object]] = []
    for endpoint, endpoint_data in data.groupby("endpoint", sort=False):
        averaged = (
            endpoint_data.groupby(["bootstrap", "method"], observed=True)[list(METRICS)]
            .mean()
            .reset_index()
        )
        available = set(averaged["method"])
        for model in SEEDED:
            if model not in available:
                continue
            for baseline in BASELINES:
                if baseline not in available:
                    continue
                for metric in METRICS:
                    model_values = averaged.loc[
                        averaged["method"].eq(model), ["bootstrap", metric]
                    ]
                    base_values = averaged.loc[
                        averaged["method"].eq(baseline), ["bootstrap", metric]
                    ]
                    paired = model_values.merge(
                        base_values,
                        on="bootstrap",
                        how="inner",
                        suffixes=("_model", "_baseline"),
                        validate="one_to_one",
                    )
                    if metric == "rmse_effect":
                        delta = paired[f"{metric}_baseline"] - paired[f"{metric}_model"]
                    else:
                        delta = paired[f"{metric}_model"] - paired[f"{metric}_baseline"]
                    values = delta.to_numpy(dtype=float)
                    rows.append(
                        {
                            "endpoint": endpoint,
                            "model": model,
                            "baseline": baseline,
                            "metric": metric,
                            "median_delta": float(np.nanmedian(values)),
                            "lower_2_5": float(np.nanquantile(values, 0.025)),
                            "upper_97_5": float(np.nanquantile(values, 0.975)),
                            "win_fraction": float(np.nanmean(values > 0)),
                            "n_paired_bootstrap_values": int(np.isfinite(values).sum()),
                        }
                    )
    return pd.DataFrame(rows)


def seed_variance(point: pd.DataFrame) -> pd.DataFrame:
    data = point.loc[
        point["cell_type"].isin(CAPILLARY) & point["method"].isin(SEEDED)
    ].copy()
    seed_level = (
        data.groupby(["endpoint", "method", "seed"], observed=True)[list(METRICS)]
        .mean()
        .reset_index()
    )
    out = seed_level.groupby(["endpoint", "method"], observed=True)[list(METRICS)].agg(
        ["mean", "std", "min", "max"]
    )
    out.columns = ["__".join(column) for column in out.columns]
    return out.reset_index()


def method_summary(point: pd.DataFrame, cohort: str) -> pd.DataFrame:
    data = point.loc[
        point["cell_type"].isin(CAPILLARY)
        & point["endpoint"].eq("primary_fold_hvg")
        & point["method"].isin(DETERMINISTIC + SEEDED)
    ].copy()
    task = (
        data.groupby(["heldout_age", "cell_type", "method"], observed=True)[list(METRICS)]
        .mean()
        .reset_index()
    )
    overall = task.groupby("method", observed=True)[list(METRICS)].mean().reset_index()
    rank_values = task.pivot_table(
        index=["heldout_age", "cell_type"], columns="method", values="spearman_effect"
    ).rank(axis=1, ascending=False, method="average")
    ranks = rank_values.mean(axis=0).rename("mean_task_rank").reset_index()
    overall = overall.merge(ranks, on="method", how="left", validate="one_to_one")
    overall.insert(0, "cohort", cohort)
    overall["n_tasks"] = task.groupby("method").size().reindex(overall["method"]).to_numpy()
    return overall.sort_values("spearman_effect", ascending=False, kind="stable")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_primary(
    r8_dir: Path,
    addition_dir: Path,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    r8_point = read_tsv(r8_dir / "virtual_cell_point_metrics.tsv")
    add_point = read_tsv(addition_dir / "expanded_point_metrics.tsv")
    r8_boot = read_tsv(r8_dir / "virtual_cell_strictly_paired_bootstrap.tsv.gz")
    add_boot = read_tsv(addition_dir / "expanded_strictly_paired_bootstrap.tsv.gz")
    r8_pred = read_tsv(r8_dir / "virtual_cell_gene_predictions.tsv.gz")
    add_pred = read_tsv(addition_dir / "expanded_gene_predictions.tsv.gz")

    if set(add_point["method"]) != set(ADDED):
        raise AssertionError(f"Unexpected primary addition methods: {sorted(add_point['method'].unique())}")

    r8_features = read_tsv(r8_dir / "virtual_cell_fold_feature_manifest.tsv")
    add_features = read_tsv(addition_dir / "expanded_fold_feature_manifest.tsv")
    feature_columns = [
        "heldout_age", "model_feature_position", "gene_index", "gene",
        "primary_fold_hvg", "secondary_replicated_198_hvg_intersection",
        "secondary_signature_33_hvg_intersection",
    ]
    assert_frames_equal(r8_features, add_features, feature_columns, "GSE151974 feature manifest")

    r8_draws = read_tsv(r8_dir / "virtual_cell_joint_capillary_bootstrap_draws.tsv.gz")
    add_draws = read_tsv(addition_dir / "expanded_joint_capillary_bootstrap_draws.tsv.gz")
    draw_columns = ["heldout_age", "bootstrap", "oxygen_group", "draw_position", "animal_id"]
    assert_frames_equal(r8_draws, add_draws, draw_columns, "GSE151974 bootstrap draw ledger")

    r8_splits = read_tsv(r8_dir / "virtual_cell_grouped_validation_splits.tsv")
    add_splits = read_tsv(addition_dir / "expanded_grouped_validation_splits.tsv")
    split_columns = ["heldout_age", "seed", "animal_id", "age", "oxygen", "role"]
    assert_frames_equal(r8_splits, add_splits, split_columns, "GSE151974 validation split ledger")

    point = pd.concat([r8_point, add_point], ignore_index=True)
    bootstrap = pd.concat([r8_boot, add_boot], ignore_index=True)
    predictions = pd.concat([r8_pred, add_pred], ignore_index=True)
    histories = pd.concat(
        [
            read_tsv(r8_dir / "virtual_cell_training_histories.tsv.gz"),
            read_tsv(addition_dir / "expanded_training_histories.tsv.gz"),
        ],
        ignore_index=True,
    )
    comparisons = paired_comparisons(bootstrap)
    variance = seed_variance(point)

    point.to_csv(output_dir / "virtual_cell_point_metrics.tsv", sep="\t", index=False)
    bootstrap.to_csv(
        output_dir / "virtual_cell_strictly_paired_bootstrap.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    predictions.to_csv(
        output_dir / "virtual_cell_gene_predictions.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    histories.to_csv(
        output_dir / "virtual_cell_training_histories.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    comparisons.to_csv(output_dir / "virtual_cell_paired_model_comparisons.tsv", sep="\t", index=False)
    variance.to_csv(output_dir / "virtual_cell_model_training_variance.tsv", sep="\t", index=False)
    shutil.copy2(
        r8_dir / "virtual_cell_fold_feature_manifest.tsv",
        output_dir / "virtual_cell_fold_feature_manifest.tsv",
    )
    shutil.copy2(
        r8_dir / "virtual_cell_grouped_validation_splits.tsv",
        output_dir / "virtual_cell_grouped_validation_splits.tsv",
    )
    shutil.copy2(
        r8_dir / "virtual_cell_joint_capillary_bootstrap_draws.tsv.gz",
        output_dir / "virtual_cell_joint_capillary_bootstrap_draws.tsv.gz",
    )
    return point, bootstrap


def write_secondary(source_dir: Path, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = {
        "expanded_point_metrics.tsv": "GSE243129_point_metrics.tsv",
        "expanded_strictly_paired_bootstrap.tsv.gz": "GSE243129_strictly_paired_bootstrap.tsv.gz",
        "expanded_gene_predictions.tsv.gz": "GSE243129_gene_predictions.tsv.gz",
        "expanded_fold_feature_manifest.tsv": "GSE243129_fold_feature_manifest.tsv",
        "expanded_grouped_validation_splits.tsv": "GSE243129_grouped_validation_splits.tsv",
        "expanded_training_histories.tsv.gz": "GSE243129_training_histories.tsv.gz",
        "expanded_joint_capillary_bootstrap_draws.tsv.gz": "GSE243129_joint_capillary_bootstrap_draws.tsv.gz",
        "expanded_capillary_animal_arrays.npz": "GSE243129_capillary_animal_arrays.npz",
        "expanded_benchmark_audit.json": "GSE243129_benchmark_audit.json",
    }
    for source, target in mapping.items():
        shutil.copy2(source_dir / source, output_dir / target)
    point = read_tsv(source_dir / "expanded_point_metrics.tsv")
    bootstrap = read_tsv(source_dir / "expanded_strictly_paired_bootstrap.tsv.gz")
    paired_comparisons(bootstrap).to_csv(
        output_dir / "GSE243129_paired_model_comparisons.tsv", sep="\t", index=False
    )
    seed_variance(point).to_csv(
        output_dir / "GSE243129_model_training_variance.tsv", sep="\t", index=False
    )
    return point, bootstrap


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r8-dir", type=Path, required=True)
    parser.add_argument("--primary-addition-dir", type=Path, required=True)
    parser.add_argument("--secondary-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    primary_point, primary_boot = write_primary(
        args.r8_dir, args.primary_addition_dir, args.output_dir
    )
    secondary_point, secondary_boot = write_secondary(args.secondary_dir, args.output_dir)
    summary = pd.concat(
        [
            method_summary(primary_point, "GSE151974"),
            method_summary(secondary_point, "GSE243129"),
        ],
        ignore_index=True,
    )
    summary.to_csv(args.output_dir / "cross_cohort_method_summary.tsv", sep="\t", index=False)

    file_hashes = {}
    for path in sorted(args.output_dir.iterdir()):
        if path.is_file() and path.name != "R9_expanded_benchmark_audit.json":
            file_hashes[path.name] = sha256(path)
    audit = {
        "analysis_version": "R9 expanded virtual-cell benchmark",
        "primary_cohort": "GSE151974",
        "independent_cohort": "GSE243129",
        "methods": list(DETERMINISTIC + SEEDED),
        "protocol": (
            "outer-age holdout; exactly 1,800 outer-training HVGs; animal-grouped "
            "validation; ten shared logical seeds for stochastic methods; joint animal "
            "bootstrap shared by capillary subtypes and all methods"
        ),
        "primary_rows": {"point": len(primary_point), "bootstrap": len(primary_boot)},
        "secondary_rows": {"point": len(secondary_point), "bootstrap": len(secondary_boot)},
        "frozen_primary_ledgers_verified": True,
        "output_sha256": file_hashes,
    }
    (args.output_dir / "R9_expanded_benchmark_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(summary[["cohort", "method", "spearman_effect", "mean_task_rank"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
