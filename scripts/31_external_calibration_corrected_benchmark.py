"""Externally calibrate the outer-fold-locked P14 virtual-cell predictions.

All methods are evaluated against the independent GSE266988 capillary
hyperoxia contrast without model selection or refitting.  Deep-model metrics
are retained for each of ten training seeds and summarized across seeds.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


GENE_SETS = {
    "primary_fold_hvg": "primary_fold_hvg",
    "secondary_replicated_198_hvg_intersection": "secondary_replicated_198_hvg_intersection",
    "secondary_signature_33_hvg_intersection": "secondary_signature_33_hvg_intersection",
}


def safe_corr(x: np.ndarray, y: np.ndarray, kind: str) -> float:
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    if kind == "pearson":
        return float(pearsonr(x, y).statistic)
    return float(spearmanr(x, y).statistic)


def read_external(path: Path) -> pd.DataFrame:
    frame = pd.read_excel(path)
    if "gene" not in frame.columns:
        frame = frame.rename(columns={frame.columns[0]: "gene"})
    if "avg_log2FC" not in frame.columns:
        raise ValueError("Expected avg_log2FC in independent capillary table")
    frame["external_effect"] = pd.to_numeric(frame["avg_log2FC"], errors="coerce")
    frame["gene"] = frame["gene"].astype(str)
    return (
        frame.dropna(subset=["gene", "external_effect"])
        .groupby("gene", as_index=False)["external_effect"]
        .mean()
    )


def metric_record(
    subset: pd.DataFrame,
    cell_type: str,
    method: str,
    seed: int,
    endpoint: str,
    effect_column: str,
) -> dict[str, object]:
    x = subset[effect_column].to_numpy(dtype=float)
    y = subset["external_effect"].to_numpy(dtype=float)
    return {
        "cell_type": cell_type,
        "method": method,
        "seed": seed,
        "endpoint": endpoint,
        "n_genes": int(len(subset)),
        "pearson": safe_corr(x, y, "pearson"),
        "spearman": safe_corr(x, y, "spearman"),
        "direction_accuracy": float(np.mean(np.sign(x) == np.sign(y))),
        "rmse": float(np.sqrt(np.mean((x - y) ** 2))),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--external-hyperoxia", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    predictions = pd.read_csv(args.predictions, sep="\t")
    predictions = predictions.loc[
        predictions["heldout_age"].eq("P14")
        & predictions["cell_type"].isin(["Cap", "Cap-a"])
    ].copy()
    predictions["predicted_effect"] = (
        predictions["predicted_hyperoxia_mean"] - predictions["normoxia_mean"]
    )
    predictions["observed_effect"] = (
        predictions["observed_hyperoxia_mean"] - predictions["normoxia_mean"]
    )
    external = read_external(args.external_hyperoxia)
    merged = predictions.merge(external, on="gene", how="inner", validate="many_to_one")

    rows: list[dict[str, object]] = []
    for (cell_type, method, seed), group in merged.groupby(
        ["cell_type", "method", "seed"], observed=True
    ):
        for endpoint, flag in GENE_SETS.items():
            subset = group.loc[group[flag].astype(bool)].copy()
            rows.append(
                metric_record(
                    subset, str(cell_type), str(method), int(seed), endpoint,
                    "predicted_effect",
                )
            )
    for cell_type, group in merged.groupby("cell_type", observed=True):
        observed = group.drop_duplicates("gene")
        for endpoint, flag in GENE_SETS.items():
            subset = observed.loc[observed[flag].astype(bool)].copy()
            rows.append(
                metric_record(
                    subset, str(cell_type), "observed", -1, endpoint,
                    "observed_effect",
                )
            )
    metrics = pd.DataFrame(rows).drop_duplicates(
        ["cell_type", "method", "seed", "endpoint"]
    )
    metrics.to_csv(
        args.output_dir / "external_calibration_seed_level_metrics.tsv",
        sep="\t", index=False,
    )
    merged.to_csv(
        args.output_dir / "external_calibration_gene_level.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )

    summary = (
        metrics.groupby(["cell_type", "method", "endpoint"], observed=True)
        .agg(
            n_training_seeds=("seed", "size"),
            n_genes=("n_genes", "first"),
            pearson_mean=("pearson", "mean"),
            pearson_min=("pearson", "min"),
            pearson_max=("pearson", "max"),
            spearman_mean=("spearman", "mean"),
            spearman_min=("spearman", "min"),
            spearman_max=("spearman", "max"),
            direction_accuracy_mean=("direction_accuracy", "mean"),
            direction_accuracy_min=("direction_accuracy", "min"),
            direction_accuracy_max=("direction_accuracy", "max"),
            rmse_mean=("rmse", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(
        args.output_dir / "external_calibration_summary.tsv",
        sep="\t", index=False,
    )
    audit = {
        "heldout_fold": "P14 with all P14 cells excluded from model fitting and feature selection",
        "external_dataset": "GSE266988 independent capillary hyperoxia contrast",
        "deep_model_seed_policy": "all ten fixed analysis seeds retained and summarized",
        "model_feature_policy": "exactly 1,800 genes selected from outer-training cells; no disease-set genes forced into the model matrix",
        "models_selected_using_external_outcomes": False,
        "interpretive_boundary": "External predictive calibration, not causal validation.",
        "external_unique_genes": int(external["gene"].nunique()),
        "merged_unique_genes": int(merged["gene"].nunique()),
    }
    (args.output_dir / "external_calibration_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(
        summary.loc[
            summary["endpoint"].isin(
                ["primary_fold_hvg", "secondary_replicated_198_hvg_intersection"]
            )
        ].to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
