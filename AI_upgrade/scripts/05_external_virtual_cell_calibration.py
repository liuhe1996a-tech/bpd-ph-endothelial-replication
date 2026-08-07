"""Externally calibrate age-held-out virtual-cell predictions.

P14 predictions generated without P14 hyperoxia cells are compared with a
published independent capillary endothelial hyperoxia contrast (JCI Insight
2025, Supplemental Table 2). Results are reported for all shared model genes,
the frozen 198-gene replicated set, and the frozen 33-gene signature. The
analysis compares every model rather than selecting a model on the external
dataset.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr


METHOD_LABELS = {
    "identity": "No-change baseline",
    "gene_linear": "Gene-space linear",
    "pca_latent": "PCA latent shift",
    "vae_latent": "VAE latent shift",
    "cvae_counterfactual": "Conditional VAE",
    "observed": "Observed P14 effect",
}


def safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    if method == "pearson":
        return float(pearsonr(x, y).statistic)
    return float(spearmanr(x, y).statistic)


def read_external(path: Path) -> pd.DataFrame:
    frame = pd.read_excel(path)
    first = frame.columns[0]
    if first != "gene":
        frame = frame.rename(columns={first: "gene"})
    frame["external_effect"] = pd.to_numeric(frame["avg_log2FC"], errors="coerce")
    frame = frame.dropna(subset=["gene", "external_effect"])
    frame["gene"] = frame["gene"].astype(str)
    return frame.groupby("gene", as_index=False)["external_effect"].mean()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--external-hyperoxia", type=Path, required=True)
    parser.add_argument("--internal-metrics", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pred = pd.read_csv(args.predictions, sep="\t")
    pred = pred[pred["heldout_age"].eq("P14") & pred["cell_type"].isin(["Cap", "Cap-a"])].copy()
    external = read_external(args.external_hyperoxia)
    merged = pred.merge(external, on="gene", how="inner", validate="many_to_one")

    rows: list[dict[str, object]] = []
    gene_sets = {
        "all_shared_model_genes": np.ones(len(merged), dtype=bool),
        "replicated_198": merged["gene_set_198"].astype(bool).to_numpy(),
        "signature_33": merged["gene_set_33"].astype(bool).to_numpy(),
    }
    for (cell_type, method), group in merged.groupby(["cell_type", "method"], observed=True):
        for gene_set in gene_sets:
            if gene_set == "all_shared_model_genes":
                subset = group
            elif gene_set == "replicated_198":
                subset = group[group["gene_set_198"].astype(bool)]
            else:
                subset = group[group["gene_set_33"].astype(bool)]
            x = subset["predicted_effect"].to_numpy(float)
            y = subset["external_effect"].to_numpy(float)
            rows.append({
                "cell_type": cell_type,
                "method": method,
                "gene_set": gene_set,
                "n_genes": len(subset),
                "pearson": safe_corr(x, y, "pearson"),
                "spearman": safe_corr(x, y, "spearman"),
                "direction_accuracy": float(np.mean(np.sign(x) == np.sign(y))),
            })

        observed = group.drop_duplicates("gene")
        for gene_set in gene_sets:
            if gene_set == "all_shared_model_genes":
                subset = observed
            elif gene_set == "replicated_198":
                subset = observed[observed["gene_set_198"].astype(bool)]
            else:
                subset = observed[observed["gene_set_33"].astype(bool)]
            x = subset["observed_effect"].to_numpy(float)
            y = subset["external_effect"].to_numpy(float)
            rows.append({
                "cell_type": cell_type,
                "method": "observed",
                "gene_set": gene_set,
                "n_genes": len(subset),
                "pearson": safe_corr(x, y, "pearson"),
                "spearman": safe_corr(x, y, "spearman"),
                "direction_accuracy": float(np.mean(np.sign(x) == np.sign(y))),
            })
    metrics = pd.DataFrame(rows).drop_duplicates(["cell_type", "method", "gene_set"])
    metrics["method_label"] = metrics["method"].map(METHOD_LABELS)
    metrics.to_csv(args.output_dir / "external_virtual_cell_calibration_metrics.tsv", sep="\t", index=False)
    merged.to_csv(args.output_dir / "external_virtual_cell_gene_level.tsv.gz", sep="\t", index=False, compression="gzip")

    internal = pd.read_csv(args.internal_metrics, sep="\t")
    internal = internal[
        internal["gene_set"].eq("replicated_198")
        & internal["cell_type"].isin(["Cap", "Cap-a"])
        & ~internal["method"].eq("identity")
    ]
    internal_summary = (
        internal.groupby("method", observed=True)
        .agg(
            median_pearson_effect=("pearson_effect", "median"),
            median_spearman_effect=("spearman_effect", "median"),
            median_rmse_effect=("rmse_effect", "median"),
            median_direction_accuracy=("direction_accuracy", "median"),
        )
        .reset_index()
    )
    internal_summary["method_label"] = internal_summary["method"].map(METHOD_LABELS)
    internal_summary.to_csv(args.output_dir / "internal_capillary_benchmark_summary.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.4))
    plot_internal = internal_summary.sort_values("median_pearson_effect", ascending=False)
    sns.barplot(data=plot_internal, y="method_label", x="median_pearson_effect", color="#4C78A8", ax=axes[0])
    axes[0].set_title("Age-held-out internal benchmark")
    axes[0].set_xlabel("Median Pearson r, 198 genes")
    axes[0].set_ylabel("")
    axes[0].set_xlim(0, 1)

    external_plot = metrics[
        metrics["gene_set"].eq("replicated_198")
        & ~metrics["method"].eq("identity")
    ].copy()
    sns.barplot(data=external_plot, y="method_label", x="spearman", hue="cell_type", palette={"Cap": "#59A14F", "Cap-a": "#E15759"}, ax=axes[1])
    axes[1].axvline(0, color="#555555", linewidth=0.8)
    axes[1].set_title("Independent capillary hyperoxia")
    axes[1].set_xlabel("Spearman r, 198 genes")
    axes[1].set_ylabel("")
    axes[1].legend(title="Held-out subtype", frameon=False)

    sig = metrics[
        metrics["gene_set"].eq("signature_33")
        & ~metrics["method"].eq("identity")
    ].copy()
    sns.barplot(data=sig, y="method_label", x="direction_accuracy", hue="cell_type", palette={"Cap": "#59A14F", "Cap-a": "#E15759"}, ax=axes[2])
    axes[2].axvline(0.5, color="#555555", linewidth=0.8, linestyle="--")
    axes[2].set_title("External signature direction")
    axes[2].set_xlabel("Direction agreement, 33 genes")
    axes[2].set_ylabel("")
    axes[2].set_xlim(0, 1)
    axes[2].legend(title="Held-out subtype", frameon=False)
    sns.despine(fig=fig)
    fig.tight_layout()
    fig.savefig(args.output_dir / "virtual_cell_external_calibration.png", dpi=300, bbox_inches="tight")
    fig.savefig(args.output_dir / "virtual_cell_external_calibration.pdf", bbox_inches="tight")
    plt.close(fig)

    audit = {
        "external_shared_gene_rows": int(len(merged)),
        "external_unique_genes": int(merged["gene"].nunique()),
        "models_compared": sorted(merged["method"].unique().tolist()),
        "selection_rule": "All models compared; no model selected using external outcomes.",
        "interpretive_boundary": "Predictive calibration, not experimental validation of a causal counterfactual.",
    }
    (args.output_dir / "external_virtual_cell_calibration_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(metrics.sort_values(["gene_set", "cell_type", "spearman"], ascending=[True, True, False]).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
