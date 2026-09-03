"""Plot the seven-method, two-cohort leakage-free perturbation benchmark."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


METHODS = [
    "gene_linear",
    "pca_latent",
    "vae_latent",
    "cvae_counterfactual",
    "scgen_adapted",
    "cpa_adapted",
    "sinkhorn_ot",
]
LABELS = {
    "gene_linear": "Gene shift",
    "pca_latent": "PCA shift",
    "vae_latent": "VAE",
    "cvae_counterfactual": "Conditional VAE",
    "scgen_adapted": "scGen-style",
    "cpa_adapted": "CPA-style",
    "sinkhorn_ot": "Sinkhorn OT",
}
COLORS = {
    "gene_linear": "#4776B4",
    "pca_latent": "#2A9D8F",
    "vae_latent": "#F4A261",
    "cvae_counterfactual": "#D65A6F",
    "scgen_adapted": "#7B6FD0",
    "cpa_adapted": "#A66B3F",
    "sinkhorn_ot": "#5B8C5A",
}


def panel(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="top")


def primary_heat(point: pd.DataFrame, ages: list[str]) -> pd.DataFrame:
    selected = point.loc[
        point["endpoint"].eq("primary_fold_hvg")
        & point["cell_type"].isin(["Cap", "Cap-a"])
        & point["method"].isin(METHODS)
    ].copy()
    selected["task"] = selected["heldout_age"] + " " + selected["cell_type"]
    order = [f"{age} {cell}" for age in ages for cell in ["Cap", "Cap-a"]]
    return (
        selected.groupby(["task", "method"], observed=True)["spearman_effect"]
        .mean().unstack("method").reindex(index=order, columns=METHODS)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--sensitivity-dir", type=Path)
    args = parser.parse_args()
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    args.source_dir.mkdir(parents=True, exist_ok=True)

    primary = pd.read_csv(args.benchmark_dir / "virtual_cell_point_metrics.tsv", sep="\t")
    secondary = pd.read_csv(args.benchmark_dir / "GSE243129_point_metrics.tsv", sep="\t")
    comparisons = pd.read_csv(
        args.benchmark_dir / "virtual_cell_paired_model_comparisons.tsv", sep="\t"
    )
    cross = pd.read_csv(args.benchmark_dir / "cross_cohort_method_summary.tsv", sep="\t")

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.0)
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 9.0))

    # A: primary disease cohort.
    heat_primary = primary_heat(primary, ["P3", "P7", "P14"])
    sns.heatmap(
        heat_primary, ax=axes[0, 0], cmap="YlGnBu", vmin=0, vmax=0.8,
        annot=True, fmt=".2f", linewidths=0.5,
        cbar_kws={"label": "Spearman correlation"},
    )
    axes[0, 0].set_title("GSE151974 held-out-age benchmark")
    axes[0, 0].set_xlabel("")
    axes[0, 0].set_ylabel("Held-out age and subtype")
    axes[0, 0].set_xticklabels([LABELS[m] for m in METHODS], rotation=28, ha="right")
    axes[0, 0].tick_params(axis="y", rotation=0)
    panel(axes[0, 0], "A")
    heat_primary.reset_index().to_csv(
        args.source_dir / "Figure5A_GSE151974_task_metrics.tsv", sep="\t", index=False
    )

    # B: biological bootstrap, with every stochastic model compared with PCA.
    comp = comparisons.loc[
        comparisons["endpoint"].eq("primary_fold_hvg")
        & comparisons["metric"].eq("spearman_effect")
        & comparisons["baseline"].eq("pca_latent")
    ].copy()
    comp["label"] = comp["model"].map(LABELS)
    comp = comp.sort_values("median_delta", kind="stable")
    y = np.arange(len(comp))
    axes[0, 1].errorbar(
        comp["median_delta"], y,
        xerr=np.vstack([
            comp["median_delta"] - comp["lower_2_5"],
            comp["upper_97_5"] - comp["median_delta"],
        ]),
        fmt="o", color="#314A67", ecolor="#7C90A6", capsize=3,
    )
    axes[0, 1].axvline(0, color="#555555", linestyle="--", linewidth=0.9)
    axes[0, 1].set_yticks(y, comp["label"])
    axes[0, 1].set_xlabel("Paired Spearman difference (model - PCA shift)")
    axes[0, 1].set_title("GSE151974 joint-animal bootstrap")
    panel(axes[0, 1], "B")
    comp.to_csv(
        args.source_dir / "Figure5B_GSE151974_paired_vs_PCA.tsv", sep="\t", index=False
    )

    # C: independent cohort, trained and evaluated only within that cohort.
    heat_secondary = primary_heat(secondary, ["P7", "P14"])
    sns.heatmap(
        heat_secondary, ax=axes[1, 0], cmap="YlGnBu", vmin=0, vmax=0.8,
        annot=True, fmt=".2f", linewidths=0.5,
        cbar_kws={"label": "Spearman correlation"},
    )
    axes[1, 0].set_title("Small GSE243129 benchmark (8 animals)")
    axes[1, 0].set_xlabel("")
    axes[1, 0].set_ylabel("Held-out age and subtype")
    axes[1, 0].set_xticklabels([LABELS[m] for m in METHODS], rotation=28, ha="right")
    axes[1, 0].tick_params(axis="y", rotation=0)
    panel(axes[1, 0], "C")
    heat_secondary.reset_index().to_csv(
        args.source_dir / "Figure5C_GSE243129_task_metrics.tsv", sep="\t", index=False
    )

    # D: average method performance in the two independently trained cohorts.
    pivot = cross.loc[cross["method"].isin(METHODS)].pivot(
        index="method", columns="cohort", values="spearman_effect"
    ).reindex(METHODS)
    for method, row in pivot.iterrows():
        axes[1, 1].scatter(
            row["GSE151974"], row["GSE243129"], s=65,
            color=COLORS[method], edgecolor="white", linewidth=0.7, zorder=3,
            label=LABELS[method],
        )
    lim = [min(-0.02, float(np.nanmin(pivot.to_numpy())) - 0.03),
           max(0.62, float(np.nanmax(pivot.to_numpy())) + 0.03)]
    axes[1, 1].plot(lim, lim, color="#999999", linestyle="--", linewidth=0.9)
    axes[1, 1].set_xlim(lim)
    axes[1, 1].set_ylim(lim)
    axes[1, 1].set_xlabel("Mean Spearman, GSE151974 (six tasks)")
    axes[1, 1].set_ylabel("Mean Spearman, GSE243129 (four tasks)")
    axes[1, 1].set_title("Point-estimate ordering differs across cohorts")
    axes[1, 1].legend(frameon=False, fontsize=7.5, loc="upper left", ncol=2)
    panel(axes[1, 1], "D")
    pivot.reset_index().to_csv(
        args.source_dir / "Figure5D_cross_cohort_method_summary.tsv", sep="\t", index=False
    )

    fig.tight_layout(h_pad=2.7, w_pad=2.8)
    stem = args.figure_dir / "Figure5_expanded_virtual_cell_benchmark"
    fig.savefig(stem.with_suffix(".png"), dpi=450, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {stem.with_suffix('.png')}")

    if args.sensitivity_dir:
        exact = pd.read_csv(
            args.sensitivity_dir / "GSE243129_all_nonempty_subset_method_ranks.tsv",
            sep="\t",
        )
        frequency_files = {
            "81 non-empty subsets":
                "GSE243129_all_nonempty_subset_rank_frequency.tsv",
            "16 one-pair subsets":
                "GSE243129_one_pair_rank_frequency.tsv",
            "8 leave-one-animal subsets":
                "GSE243129_leave_one_animal_rank_frequency.tsv",
        }
        sensitivity_fig, sensitivity_axes = plt.subplots(1, 2, figsize=(12.2, 4.5))
        plot_exact = exact.copy()
        plot_exact["Method"] = plot_exact["method"].map(LABELS)
        method_labels = [LABELS[method] for method in METHODS]
        sns.boxplot(
            data=plot_exact,
            x="Method",
            y="mean_spearman_across_four_tasks",
            order=method_labels,
            color="#BBD7E8",
            fliersize=2,
            linewidth=0.9,
            ax=sensitivity_axes[0],
        )
        sns.stripplot(
            data=plot_exact,
            x="Method",
            y="mean_spearman_across_four_tasks",
            order=method_labels,
            color="#345B73",
            alpha=0.34,
            size=2.2,
            jitter=0.18,
            ax=sensitivity_axes[0],
        )
        sensitivity_axes[0].set_xlabel("")
        sensitivity_axes[0].set_ylabel("Mean Spearman across four tasks")
        sensitivity_axes[0].set_title("Exact non-empty animal-subset analysis")
        sensitivity_axes[0].tick_params(axis="x", rotation=32)
        for label in sensitivity_axes[0].get_xticklabels():
            label.set_ha("right")
        panel(sensitivity_axes[0], "A")

        frequency_frames = []
        for label, filename in frequency_files.items():
            frame = pd.read_csv(args.sensitivity_dir / filename, sep="\t")
            frame["Sensitivity analysis"] = label
            frame["Method"] = frame["method"].map(LABELS)
            frequency_frames.append(frame)
        frequencies = pd.concat(frequency_frames, ignore_index=True)
        sns.barplot(
            data=frequencies,
            x="Method",
            y="fraction_ranked_first",
            hue="Sensitivity analysis",
            order=method_labels,
            palette=["#345B73", "#6A9C89", "#D28B5C"],
            ax=sensitivity_axes[1],
        )
        sensitivity_axes[1].set_xlabel("")
        sensitivity_axes[1].set_ylabel("Fraction ranked first")
        sensitivity_axes[1].set_ylim(0, 1.05)
        sensitivity_axes[1].set_title("First-place frequency under animal removal")
        sensitivity_axes[1].tick_params(axis="x", rotation=32)
        for label in sensitivity_axes[1].get_xticklabels():
            label.set_ha("right")
        legend = sensitivity_axes[1].get_legend()
        if legend is not None:
            legend.remove()
        sensitivity_axes[1].text(
            0.98,
            0.98,
            "Blue: 81 subsets   Green: 16 one-pair   Orange: 8 leave-one-out",
            transform=sensitivity_axes[1].transAxes,
            ha="right",
            va="top",
            fontsize=7.2,
        )
        panel(sensitivity_axes[1], "B")
        sensitivity_fig.tight_layout(w_pad=2.2)
        sensitivity_stem = (
            args.figure_dir / "Supplementary_Figure_GSE243129_exact_animal_sensitivity"
        )
        sensitivity_fig.savefig(
            sensitivity_stem.with_suffix(".png"), dpi=450, bbox_inches="tight"
        )
        sensitivity_fig.savefig(
            sensitivity_stem.with_suffix(".pdf"), bbox_inches="tight"
        )
        plt.close(sensitivity_fig)
        plot_exact.to_csv(
            args.source_dir / "Supplementary_Figure_GSE243129_exact_subset_metrics.tsv",
            sep="\t",
            index=False,
        )
        frequencies.to_csv(
            args.source_dir / "Supplementary_Figure_GSE243129_rank_frequencies.tsv",
            sep="\t",
            index=False,
        )
        print(f"Wrote {sensitivity_stem.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
