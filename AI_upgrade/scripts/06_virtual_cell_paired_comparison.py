"""Paired bootstrap comparison of virtual-cell models and simple baselines."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


LABELS = {
    "gene_linear": "Gene-space linear",
    "pca_latent": "PCA latent shift",
    "vae_latent": "VAE latent shift",
    "cvae_counterfactual": "Conditional VAE",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(args.bootstrap, sep="\t")
    data = data[
        data["gene_set"].eq("replicated_198")
        & data["cell_type"].isin(["Cap", "Cap-a"])
        & data["method"].isin(LABELS)
    ].copy()
    metrics = ["pearson_effect", "spearman_effect", "direction_accuracy", "rmse_effect"]
    task_avg = (
        data.groupby(["bootstrap", "method"], observed=True)[metrics]
        .mean()
        .reset_index()
    )
    wide = task_avg.pivot(index="bootstrap", columns="method", values=metrics)
    rows: list[dict[str, object]] = []
    comparisons = [
        ("vae_latent", "gene_linear"),
        ("vae_latent", "pca_latent"),
        ("cvae_counterfactual", "gene_linear"),
        ("cvae_counterfactual", "pca_latent"),
    ]
    for deep, baseline in comparisons:
        for metric in metrics:
            if metric == "rmse_effect":
                delta = wide[(metric, baseline)] - wide[(metric, deep)]
                label = "RMSE improvement"
            else:
                delta = wide[(metric, deep)] - wide[(metric, baseline)]
                label = metric.replace("_effect", "").replace("_", " ").title()
            rows.append({
                "deep_model": deep,
                "baseline": baseline,
                "metric": metric,
                "metric_label": label,
                "median_delta": float(delta.median()),
                "lower_2_5": float(delta.quantile(0.025)),
                "upper_97_5": float(delta.quantile(0.975)),
                "bootstrap_win_fraction": float((delta > 0).mean()),
                "n_bootstrap": int(delta.notna().sum()),
            })
    summary = pd.DataFrame(rows)
    summary["comparison"] = (
        summary["deep_model"].map(LABELS)
        + " vs "
        + summary["baseline"].map(LABELS)
    )
    summary.to_csv(args.output_dir / "virtual_cell_paired_bootstrap_comparisons.tsv", sep="\t", index=False)

    plot = summary[summary["metric"].isin(["pearson_effect", "direction_accuracy", "rmse_effect"])].copy()
    plot["y"] = np.arange(len(plot))
    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    colors = {"Pearson": "#4C78A8", "Direction Accuracy": "#59A14F", "RMSE improvement": "#E15759"}
    for metric_label, group in plot.groupby("metric_label", observed=True):
        ax.errorbar(
            group["median_delta"],
            group["y"],
            xerr=[group["median_delta"] - group["lower_2_5"], group["upper_97_5"] - group["median_delta"]],
            fmt="o",
            color=colors[metric_label],
            capsize=3,
            label=metric_label,
        )
    ax.axvline(0, color="#444444", linewidth=0.9)
    ax.set_yticks(plot["y"], plot["comparison"] + " | " + plot["metric_label"])
    ax.set_xlabel("Paired bootstrap performance difference (>0 favors deep model)")
    ax.set_ylabel("")
    ax.set_title("Deep virtual-cell models do not consistently exceed simple baselines")
    ax.legend(frameon=False, loc="lower right")
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(args.output_dir / "virtual_cell_paired_bootstrap.png", dpi=300, bbox_inches="tight")
    fig.savefig(args.output_dir / "virtual_cell_paired_bootstrap.pdf", bbox_inches="tight")
    plt.close(fig)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
