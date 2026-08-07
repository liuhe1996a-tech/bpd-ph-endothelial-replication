"""Build the upgraded manuscript figures from frozen row-level results."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns


INK = "#28323C"
BLUE = "#4C78A8"
ORANGE = "#F58518"
GREEN = "#59A14F"
RED = "#E15759"
PURPLE = "#B279A2"
GOLD = "#D4A72C"
GREY = "#8D99A6"

METHOD_LABELS = {
    "gene_linear": "Gene-space linear",
    "pca_latent": "PCA latent shift",
    "vae_latent": "VAE latent shift",
    "cvae_counterfactual": "Conditional VAE",
    "observed": "Observed P14 effect",
}


def save(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.png", dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_design(output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(16, 8.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.text(
        0.035,
        0.94,
        "Figure 1. Benchmark-guided multi-cohort evidence design",
        fontsize=22,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.035,
        0.895,
        "Independent biological units, held-out prediction and genetic calibration constrain the computational claims.",
        fontsize=12,
        color="#65727E",
    )

    boxes = [
        (0.035, 0.50, 0.155, 0.28, "Mouse\ndiscovery", "GSE216046\nP14 purified endothelium\n4 air + 4 hyperoxia\nDESeq2 and ranked GSEA", "#E5EFF7", BLUE),
        (0.225, 0.50, 0.155, 0.28, "Independent\nreplication", "GSE151974\n36 animals; P3, P7, P14\nAnimal × subtype pseudobulk\n198 genes; 33-gene signature", "#E5EFF7", BLUE),
        (0.415, 0.50, 0.155, 0.28, "Virtual-cell\nbenchmark", "Age–hyperoxia held out\nLinear and PCA baselines\nVAE and conditional VAE\nExternal capillary calibration", "#F4EBD2", GOLD),
        (0.605, 0.50, 0.155, 0.28, "Network\ninterpretation", "CollecTRI TF activities\nRegulon-constrained interventions\nAnimal-balanced ligand–receptor\nNo causal edge claim", "#F2E7EE", PURPLE),
        (0.795, 0.50, 0.17, 0.28, "External projection", "JCI p53 perturbation\nGenetic calibration\nGSE275938 human donors\nExploratory direction only", "#E7F2EC", GREEN),
    ]
    for x, y, w, h, title, body, fill, edge in boxes:
        patch = FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.012,rounding_size=0.012",
            linewidth=2.2,
            edgecolor=edge,
            facecolor=fill,
        )
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h - 0.068, title, ha="center", va="center", fontsize=11.1, fontweight="bold", color=INK, linespacing=1.05)
        ax.text(x + w / 2, y + 0.105, body, ha="center", va="center", fontsize=9.4, color=INK, linespacing=1.35)
    for left, right in zip(boxes[:-1], boxes[1:]):
        x1 = left[0] + left[2]
        x2 = right[0]
        ax.add_patch(FancyArrowPatch((x1 + 0.004, 0.64), (x2 - 0.004, 0.64), arrowstyle="-|>", mutation_scale=20, linewidth=2, color="#687680"))

    guardrails = [
        (0.07, "Unit of inference", "sample / animal / donor"),
        (0.31, "Prediction audit", "held-out age–hyperoxia condition"),
        (0.55, "External calibration", "independent hyperoxia and p53 data"),
        (0.79, "Interpretive boundary", "prioritization, not causal validation"),
    ]
    for x, title, body in guardrails:
        patch = FancyBboxPatch((x, 0.20), 0.18, 0.15, boxstyle="round,pad=0.01", linewidth=1.3, edgecolor="#B5C0C8", facecolor="#F6F8FA")
        ax.add_patch(patch)
        ax.text(x + 0.09, 0.295, title, ha="center", va="center", fontsize=10.8, fontweight="bold", color=INK)
        ax.text(x + 0.09, 0.245, body, ha="center", va="center", fontsize=8.7, color="#596772")
    ax.text(0.50, 0.095, "Human data never define the mouse replicated set; deep-model superiority is not assumed.", ha="center", fontsize=11.2, color="#596772", style="italic")
    save(fig, output_dir, "Figure1_benchmark_guided_design")


def build_virtual_cell(results_root: Path, output_dir: Path) -> None:
    metrics = pd.read_csv(results_root / "virtual_cell_all_ages" / "virtual_cell_benchmark_metrics.tsv", sep="\t")
    paired = pd.read_csv(results_root / "virtual_cell_model_comparison" / "virtual_cell_paired_bootstrap_comparisons.tsv", sep="\t")
    external = pd.read_csv(results_root / "external_virtual_cell_calibration" / "external_virtual_cell_calibration_metrics.tsv", sep="\t")

    sns.set_theme(style="whitegrid", font_scale=0.92)
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, left=0.07, right=0.98, bottom=0.08, top=0.87, hspace=0.40, wspace=0.30)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    fig.suptitle("Figure 5. Virtual-cell predictions require simple baselines and external calibration", x=0.04, ha="left", fontsize=21, fontweight="bold", color=INK)
    fig.text(0.04, 0.90, "Models were trained without hyperoxia cells from the held-out age; performance was scored on animal-level effects.", fontsize=11.5, color="#65727E")

    keep_methods = ["gene_linear", "pca_latent", "vae_latent", "cvae_counterfactual"]
    heat = metrics.query("gene_set == 'replicated_198' and cell_type in ['Cap','Cap-a'] and method in @keep_methods").copy()
    heat["task"] = heat["heldout_age"] + " " + heat["cell_type"]
    heat["method_label"] = heat["method"].map(METHOD_LABELS)
    task_order = [f"{age} {ct}" for age in ["P3", "P7", "P14"] for ct in ["Cap", "Cap-a"]]
    method_order = [METHOD_LABELS[m] for m in keep_methods]
    mat = heat.pivot(index="method_label", columns="task", values="pearson_effect").reindex(index=method_order, columns=task_order)
    sns.heatmap(mat, vmin=0.5, vmax=1.0, cmap="Blues", annot=True, fmt=".2f", cbar_kws={"label": "Pearson r"}, ax=ax_a)
    ax_a.set_title("A  Held-out age–hyperoxia prediction of the 198-gene effect", loc="left", fontweight="bold", color=INK)
    ax_a.set_xlabel("")
    ax_a.set_ylabel("")
    ax_a.tick_params(axis="x", rotation=30)

    forest = paired.query("metric in ['pearson_effect','direction_accuracy']").copy()
    forest["comparison_short"] = forest.apply(lambda r: f"{METHOD_LABELS[r.deep_model]} vs {METHOD_LABELS[r.baseline]}", axis=1)
    forest["metric_label"] = forest["metric"].map({"pearson_effect": "Pearson", "direction_accuracy": "Direction"})
    comps = forest["comparison_short"].drop_duplicates().tolist()
    ypos = {c: i for i, c in enumerate(comps[::-1])}
    offsets = {"Pearson": -0.13, "Direction": 0.13}
    colors = {"Pearson": BLUE, "Direction": GREEN}
    for _, row in forest.iterrows():
        y = ypos[row["comparison_short"]] + offsets[row["metric_label"]]
        ax_b.errorbar(row["median_delta"], y, xerr=[[row["median_delta"] - row["lower_2_5"]], [row["upper_97_5"] - row["median_delta"]]], fmt="o", color=colors[row["metric_label"]], capsize=3, label=row["metric_label"])
    ax_b.axvline(0, color="#444444", linewidth=1)
    ax_b.set_yticks(range(len(comps)), comps[::-1])
    ax_b.set_xlabel("Paired bootstrap difference (>0 favors deep model)")
    ax_b.set_title("B  Deep models do not consistently outperform baselines", loc="left", fontweight="bold", color=INK)
    handles, labels = ax_b.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax_b.legend(unique.values(), unique.keys(), frameon=False, loc="lower right")

    ext_198 = external.query("gene_set == 'replicated_198' and method != 'identity'").copy()
    ext_198["method_label"] = ext_198["method"].map(METHOD_LABELS)
    ext_order = [METHOD_LABELS[m] for m in ["observed", "gene_linear", "pca_latent", "vae_latent", "cvae_counterfactual"]]
    sns.barplot(data=ext_198, x="spearman", y="method_label", hue="cell_type", order=ext_order, palette={"Cap": BLUE, "Cap-a": ORANGE}, ax=ax_c)
    ax_c.set_xlim(0, 0.8)
    ax_c.set_xlabel("Spearman r with independent capillary hyperoxia")
    ax_c.set_ylabel("")
    ax_c.set_title("C  Independent-cohort calibration, 198-gene set", loc="left", fontweight="bold", color=INK)
    ax_c.legend(title="Held-out subtype", frameon=False, loc="lower right")

    ext_33 = external.query("gene_set == 'signature_33' and method != 'identity'").copy()
    ext_33["method_label"] = ext_33["method"].map(METHOD_LABELS)
    sns.barplot(data=ext_33, x="direction_accuracy", y="method_label", hue="cell_type", order=ext_order, palette={"Cap": BLUE, "Cap-a": ORANGE}, ax=ax_d)
    ax_d.set_xlim(0, 1.05)
    ax_d.axvline(0.5, color="#555555", linestyle="--", linewidth=1)
    ax_d.set_xlabel("Direction agreement in independent cohort")
    ax_d.set_ylabel("")
    ax_d.set_title("D  Independent signature-direction recovery", loc="left", fontweight="bold", color=INK)
    ax_d.legend(title="Held-out subtype", frameon=False, loc="lower right")

    for ax in (ax_a, ax_b, ax_c, ax_d):
        ax.tick_params(colors=INK)
    save(fig, output_dir, "Figure5_virtual_cell_benchmark")


def build_network(results_root: Path, output_dir: Path) -> None:
    tf = pd.read_csv(results_root / "regulon" / "GSE151974_collectri_ulm_activities.tsv", sep="\t")
    priority = pd.read_csv(results_root / "regulon" / "regulator_priority_summary.tsv", sep="\t")
    jci = pd.read_csv(results_root / "regulon" / "JCI_p53_genetic_calibration_ulm.tsv", sep="\t")
    genes = pd.read_csv(results_root / "regulon" / "JCI_signature_gene_genetic_calibration.tsv", sep="\t")
    lr_path = results_root / "ligand_receptor" / "GSE151974_recurrent_lr_summary.tsv"
    lr = pd.read_csv(lr_path, sep="\t") if lr_path.exists() else pd.DataFrame()

    sns.set_theme(style="whitegrid", font_scale=0.90)
    fig = plt.figure(figsize=(16, 10.5))
    gs = fig.add_gridspec(2, 2, left=0.07, right=0.98, bottom=0.08, top=0.87, hspace=0.40, wspace=0.35)
    ax_a, ax_b, ax_c, ax_d = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]
    fig.suptitle("Figure 6. Network-constrained analyses prioritize context-dependent endothelial regulators", x=0.04, ha="left", fontsize=21, fontweight="bold", color=INK)
    fig.text(0.04, 0.90, "Regulon and ligand–receptor results are expression-supported hypotheses; genetic data provide the independent calibration layer.", fontsize=11.4, color="#65727E")

    top_tfs = priority.head(12)["tf"].tolist()
    local = tf.query("tf in @top_tfs and contrast.str.contains('Cap')", engine="python").copy()
    contrast_order = [f"{age}__{ct}" for age in ["P3", "P7", "P14"] for ct in ["Cap", "Cap-a"]]
    mat = local.pivot(index="tf", columns="contrast", values="ulm_score").reindex(index=top_tfs, columns=contrast_order)
    mat.columns = [c.replace("__", " ") for c in mat.columns]
    sns.heatmap(mat, cmap="Reds", center=0, annot=True, fmt=".1f", cbar_kws={"label": "ULM activity"}, ax=ax_a)
    ax_a.set_title("A  Recurrent capillary TF activities", loc="left", fontweight="bold", color=INK)
    ax_a.set_xlabel("")
    ax_a.set_ylabel("")
    ax_a.tick_params(axis="x", rotation=30)

    p53 = jci.query("tf == 'Trp53'").copy()
    labels = {
        "JCI_hyperoxia_vs_room_air": "Independent hyperoxia vs room air",
        "JCI_p53_null_vs_control_hyperoxia": "Global p53-null vs control",
        "JCI_p53_EC_deletion_vs_control_hyperoxia": "Endothelial p53 deletion vs control",
    }
    p53["label"] = p53["contrast"].map(labels)
    sns.barplot(data=p53, x="ulm_score", y="label", color=ORANGE, ax=ax_b)
    ax_b.axvline(0, color="#444444", linewidth=1)
    ax_b.set_xlabel("Trp53 regulon activity score")
    ax_b.set_ylabel("")
    ax_b.set_title("B  Genetic p53 calibration is context dependent", loc="left", fontweight="bold", color=INK)
    for i, row in p53.reset_index(drop=True).iterrows():
        ax_b.text(row["ulm_score"] + 0.12, i, f"q={row['fdr_bh']:.3g}", va="center", ha="left", fontsize=9)

    sig = genes.query("included_high_confidence_ge3 == True").dropna(subset=["JCI_hyperoxia_vs_room_air", "JCI_p53_null_vs_control_hyperoxia"]).copy()
    sig["reversed"] = sig["JCI_hyperoxia_vs_room_air"] * sig["JCI_p53_null_vs_control_hyperoxia"] < 0
    sig = sig.sort_values("JCI_hyperoxia_vs_room_air")
    colors = np.where(sig["reversed"], GREEN, GREY)
    ax_c.scatter(sig["JCI_hyperoxia_vs_room_air"], sig["JCI_p53_null_vs_control_hyperoxia"], c=colors, s=52, edgecolor="white", linewidth=0.6)
    ax_c.axhline(0, color="#555555", linewidth=0.8)
    ax_c.axvline(0, color="#555555", linewidth=0.8)
    ax_c.set_xlabel("Independent hyperoxia log2 fold-change")
    ax_c.set_ylabel("Global p53-null effect under hyperoxia")
    nrev = int(sig["reversed"].sum())
    ax_c.set_title(f"C  Global p53 loss reverses {nrev}/{len(sig)} evaluable signature genes", loc="left", fontweight="bold", color=INK)
    for _, row in sig.loc[sig["reversed"]].sort_values("JCI_p53_null_vs_control_hyperoxia").head(6).iterrows():
        ax_c.text(row["JCI_hyperoxia_vs_room_air"], row["JCI_p53_null_vs_control_hyperoxia"], " " + row["gene"], fontsize=8, va="center")

    if not lr.empty:
        top = lr.query("positive_ages >= 2 and nominal_ages >= 1").head(12).copy()
        top["route"] = top["source"] + " → " + top["target"] + ": " + top["pair"]
        top = top.sort_values("priority_score", ascending=True)
        ax_d.barh(top["route"], top["priority_score"], color=np.where(top["reported_in_external_p53_cellchat"], GOLD, PURPLE))
        ax_d.set_xlabel("Recurrence priority score")
        ax_d.set_ylabel("")
        ax_d.set_title("D  Recurrent expression-supported communication routes", loc="left", fontweight="bold", color=INK)
        ax_d.text(0.98, 0.03, "gold: also reported in independent p53 CellChat", transform=ax_d.transAxes, ha="right", fontsize=8.5, color="#65727E")
    else:
        ax_d.axis("off")
        ax_d.text(0.5, 0.5, "Ligand–receptor results pending", ha="center", va="center")

    for ax in (ax_a, ax_b, ax_c, ax_d):
        ax.tick_params(colors=INK)
    save(fig, output_dir, "Figure6_regulatory_and_communication")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    build_design(args.output_dir)
    build_virtual_cell(args.results_root, args.output_dir)
    build_network(args.results_root, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
