"""Infer recurrent TF activity and benchmark a regulon-constrained perturbation.

CollecTRI signed mouse regulons and decoupler ULM are applied to animal-
balanced hyperoxia contrasts from GSE151974.  A simple network projection asks
which TF activity component, if returned to baseline, most reduces the
replicated endothelial signature.  Published p53-null and endothelial-specific
p53-deletion contrasts (JCI Insight 2025, Supplemental Tables 3 and 4) are held
out as genetic perturbation calibration data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from scipy.stats import pearsonr, spearmanr


AGES = ("P3", "P7", "P14")
CELL_TYPES = ("Cap", "Cap-a", "Art", "Vein", "Lymph")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_log1p_10k(counts: sparse.csr_matrix) -> sparse.csr_matrix:
    library = np.asarray(counts.sum(axis=1)).ravel()
    if np.any(library <= 0):
        raise ValueError("Cells with zero library size detected.")
    norm = counts.multiply((10000.0 / library)[:, None]).tocsr().astype(np.float32)
    norm.data = np.log1p(norm.data)
    return norm


def animal_balanced_contrasts(
    norm: sparse.csr_matrix,
    metadata: pd.DataFrame,
    genes: np.ndarray,
) -> pd.DataFrame:
    rows: list[pd.Series] = []
    labels: list[str] = []
    for age in AGES:
        for cell_type in CELL_TYPES:
            group_means: dict[str, np.ndarray] = {}
            for oxygen in ("Normoxia", "Hyperoxia"):
                subset = (
                    metadata["Age"].eq(age)
                    & metadata["CellType"].eq(cell_type)
                    & metadata["Oxygen"].eq(oxygen)
                )
                submeta = metadata.loc[subset]
                values = norm[subset.to_numpy()]
                animal_rows: list[np.ndarray] = []
                for animal in sorted(submeta["animal_id"].unique()):
                    animal_mask = submeta["animal_id"].eq(animal).to_numpy()
                    animal_rows.append(np.asarray(values[animal_mask].mean(axis=0)).ravel())
                if not animal_rows:
                    raise ValueError(f"Missing group: {age} {cell_type} {oxygen}")
                group_means[oxygen] = np.vstack(animal_rows).mean(axis=0)
            labels.append(f"{age}__{cell_type}")
            rows.append(
                pd.Series(
                    group_means["Hyperoxia"] - group_means["Normoxia"],
                    index=genes,
                )
            )
    return pd.DataFrame(rows, index=labels)


def long_ulm(
    effects: pd.DataFrame, net: pd.DataFrame, cohort: str
) -> pd.DataFrame:
    scores, padj = dc.mt.ulm(effects, net, tmin=10, verbose=True)
    score_long = scores.rename_axis("contrast").reset_index().melt(
        id_vars="contrast", var_name="tf", value_name="ulm_score"
    )
    padj_long = padj.rename_axis("contrast").reset_index().melt(
        id_vars="contrast", var_name="tf", value_name="fdr_bh"
    )
    merged = score_long.merge(padj_long, on=["contrast", "tf"], validate="one_to_one")
    merged.insert(0, "cohort", cohort)
    return merged


def read_jci_table(path: Path, label: str) -> pd.DataFrame:
    frame = pd.read_excel(path)
    first = frame.columns[0]
    if first != "gene":
        frame = frame.rename(columns={first: "gene"})
    frame = frame.dropna(subset=["gene", "avg_log2FC"])
    frame["gene"] = frame["gene"].astype(str)
    values = pd.to_numeric(frame["avg_log2FC"], errors="coerce")
    frame = frame.loc[values.notna()].copy()
    frame["avg_log2FC"] = values.loc[values.notna()].astype(float)
    frame = frame.groupby("gene", as_index=False)["avg_log2FC"].mean()
    return frame.set_index("gene")[["avg_log2FC"]].T.rename(
        index={"avg_log2FC": label}
    )


def ols_beta(y: np.ndarray, target_idx: np.ndarray, weights: np.ndarray) -> float:
    n = len(y)
    sum_w = float(weights.sum())
    sum_w2 = float(np.square(weights).sum())
    sum_y = float(y.sum())
    sum_wy = float(np.dot(weights, y[target_idx]))
    denominator = sum_w2 - (sum_w * sum_w / n)
    if denominator <= 1e-12:
        return float("nan")
    return (sum_wy - sum_w * sum_y / n) / denominator


def perturbation_projection(
    effects: pd.DataFrame,
    net: pd.DataFrame,
    signature: pd.DataFrame,
) -> pd.DataFrame:
    genes = effects.columns.astype(str).to_numpy()
    gene_to_idx = {gene: i for i, gene in enumerate(genes)}
    signature_198 = np.array(
        [gene_to_idx[g] for g in signature["gene"] if g in gene_to_idx], dtype=int
    )
    signature_33 = np.array(
        [
            gene_to_idx[row.gene]
            for row in signature.itertuples()
            if row.included_high_confidence_ge3 and row.gene in gene_to_idx
        ],
        dtype=int,
    )
    tf_network: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for tf, group in net.groupby("source", sort=False):
        group = group.loc[group["target"].isin(gene_to_idx)].drop_duplicates("target")
        if len(group) < 10:
            continue
        idx = np.array([gene_to_idx[g] for g in group["target"]], dtype=int)
        weights = group["weight"].to_numpy(dtype=float)
        tf_network[str(tf)] = (idx, weights)

    rows: list[dict[str, object]] = []
    for contrast, effect_row in effects.iterrows():
        y = effect_row.to_numpy(dtype=float)
        for tf, (target_idx, weights) in tf_network.items():
            beta = ols_beta(y, target_idx, weights)
            if not np.isfinite(beta):
                continue
            contribution = np.zeros(len(y), dtype=float)
            contribution[target_idx] = beta * weights
            residual = y - contribution
            row: dict[str, object] = {
                "contrast": contrast,
                "tf": tf,
                "beta": beta,
                "n_targets": len(target_idx),
                "suggested_intervention": "inhibit" if beta > 0 else "activate",
            }
            for name, idx in (
                ("all", np.arange(len(y))),
                ("replicated_198", signature_198),
                ("signature_33", signature_33),
            ):
                before = float(np.sqrt(np.mean(np.square(y[idx]))))
                after = float(np.sqrt(np.mean(np.square(residual[idx]))))
                row[f"rmse_before_{name}"] = before
                row[f"rmse_after_{name}"] = after
                row[f"rescue_fraction_{name}"] = (
                    1.0 - after / before if before > 0 else float("nan")
                )
            rows.append(row)
    return pd.DataFrame(rows)


def build_candidate_summary(
    activities: pd.DataFrame,
    projection: pd.DataFrame,
) -> pd.DataFrame:
    primary = activities.loc[
        activities["contrast"].str.endswith(("__Cap", "__Cap-a"))
    ].copy()
    med = primary.groupby("tf")["ulm_score"].median()
    primary["median_direction"] = primary["tf"].map(np.sign(med))
    primary["direction_match"] = (
        np.sign(primary["ulm_score"]) == primary["median_direction"]
    )
    activity_summary = primary.groupby("tf").agg(
        median_ulm_score=("ulm_score", "median"),
        mean_abs_ulm_score=("ulm_score", lambda x: float(np.mean(np.abs(x)))),
        significant_contrasts=("fdr_bh", lambda x: int(np.sum(x < 0.05))),
        direction_consistency=("direction_match", "mean"),
        median_fdr=("fdr_bh", "median"),
    )
    proj_primary = projection.loc[
        projection["contrast"].str.endswith(("__Cap", "__Cap-a"))
    ]
    projection_summary = proj_primary.groupby("tf").agg(
        n_targets=("n_targets", "max"),
        median_beta=("beta", "median"),
        median_rescue_198=("rescue_fraction_replicated_198", "median"),
        median_rescue_33=("rescue_fraction_signature_33", "median"),
        positive_rescue_contrasts_33=(
            "rescue_fraction_signature_33", lambda x: int(np.sum(x > 0))
        ),
    )
    summary = activity_summary.join(projection_summary, how="inner").reset_index()
    summary["suggested_intervention"] = np.where(
        summary["median_beta"] > 0, "inhibit", "activate"
    )
    summary["priority_score"] = (
        summary["mean_abs_ulm_score"]
        * summary["direction_consistency"]
        * (1 + summary["significant_contrasts"] / 6)
        * np.maximum(summary["median_rescue_33"], 0)
    )
    return summary.sort_values(
        ["priority_score", "significant_contrasts", "mean_abs_ulm_score"],
        ascending=False,
    )


def rank_p53(external: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for contrast, group in external.groupby("contrast"):
        ordered = group.sort_values("ulm_score", ascending=True).reset_index(drop=True)
        hit = ordered.loc[ordered["tf"].eq("Trp53")]
        if hit.empty:
            continue
        index = int(hit.index[0])
        rows.append({
            "contrast": contrast,
            "trp53_ulm_score": float(hit.iloc[0]["ulm_score"]),
            "trp53_fdr_bh": float(hit.iloc[0]["fdr_bh"]),
            "rank_most_negative": index + 1,
            "n_tfs": len(ordered),
            "negative_rank_percentile": 1.0 - index / max(1, len(ordered) - 1),
        })
    return pd.DataFrame(rows)


def plot_results(
    activities: pd.DataFrame,
    external: pd.DataFrame,
    candidates: pd.DataFrame,
    output: Path,
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    eligible = candidates.loc[
        (candidates["significant_contrasts"] >= 3)
        & (candidates["direction_consistency"] >= 0.8)
    ].head(14)
    top_tfs = eligible["tf"].tolist()
    if "Trp53" not in top_tfs:
        top_tfs = ["Trp53", *top_tfs[:13]]
    primary = activities.loc[
        activities["tf"].isin(top_tfs)
        & activities["contrast"].str.endswith(("__Cap", "__Cap-a"))
    ]
    heat = primary.pivot(index="tf", columns="contrast", values="ulm_score").reindex(top_tfs)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6.2), gridspec_kw={"width_ratios": [1.65, 0.8, 1.05]})
    sns.heatmap(
        heat, cmap="vlag", center=0, ax=axes[0], cbar_kws={"label": "ULM TF activity"}
    )
    axes[0].set_title("Recurrent capillary TF activities")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("")
    axes[0].tick_params(axis="x", rotation=40, labelsize=9)

    trp53 = pd.concat([
        activities.loc[activities["tf"].eq("Trp53"), ["contrast", "ulm_score"]].assign(source="GSE151974"),
        external.loc[external["tf"].eq("Trp53"), ["contrast", "ulm_score"]].assign(source="JCI genetic calibration"),
    ], ignore_index=True)
    sns.barplot(data=trp53, y="contrast", x="ulm_score", hue="source", ax=axes[1], dodge=False)
    axes[1].axvline(0, color="#555555", lw=0.8)
    axes[1].set_title("Trp53 activity reversal")
    axes[1].set_xlabel("ULM score")
    axes[1].set_ylabel("")
    axes[1].legend(fontsize=8, loc="lower right")

    plot_candidates = candidates.loc[
        (candidates["significant_contrasts"] >= 3)
        & (candidates["direction_consistency"] >= 0.8)
    ].head(12).sort_values("median_rescue_33")
    colors = np.where(plot_candidates["suggested_intervention"].eq("inhibit"), "#C06C84", "#355C7D")
    axes[2].barh(plot_candidates["tf"], plot_candidates["median_rescue_33"], color=colors)
    axes[2].axvline(0, color="#555555", lw=0.8)
    axes[2].set_xlabel("Median projected reduction\nin 33-gene effect magnitude")
    axes[2].set_title("Regulon-constrained interventions")
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--network", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--jci-hyperoxia", type=Path, required=True)
    parser.add_argument("--jci-p53-null", type=Path, required=True)
    parser.add_argument("--jci-p53-ec", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    counts = sparse.load_npz(args.matrix).tocsr()
    metadata = pd.read_csv(args.metadata, sep="\t")
    genes = pd.read_csv(args.genes, sep="\t")["gene"].astype(str).to_numpy()
    net = pd.read_csv(args.network, sep="\t")[["source", "target", "weight"]].drop_duplicates()
    signature = pd.read_csv(args.signature, sep="\t")
    norm = normalize_log1p_10k(counts)
    effects = animal_balanced_contrasts(norm, metadata, genes)
    effects.to_csv(args.output_dir / "GSE151974_animal_balanced_age_subtype_effects.tsv.gz", sep="\t", compression="gzip")

    activities = long_ulm(effects, net, "GSE151974")
    activities.to_csv(args.output_dir / "GSE151974_collectri_ulm_activities.tsv", sep="\t", index=False)
    projection = perturbation_projection(effects, net, signature)
    projection.to_csv(args.output_dir / "GSE151974_regulon_perturbation_projection.tsv.gz", sep="\t", index=False, compression="gzip")
    candidates = build_candidate_summary(activities, projection)

    external_effects = pd.concat([
        read_jci_table(args.jci_hyperoxia, "JCI_hyperoxia_vs_room_air"),
        read_jci_table(args.jci_p53_null, "JCI_p53_null_vs_control_hyperoxia"),
        read_jci_table(args.jci_p53_ec, "JCI_p53_EC_deletion_vs_control_hyperoxia"),
    ], axis=0, join="inner", sort=True)
    external = long_ulm(external_effects, net, "JCI_Insight_2025")
    external.to_csv(args.output_dir / "JCI_p53_genetic_calibration_ulm.tsv", sep="\t", index=False)
    p53_rank = rank_p53(external)
    p53_rank.to_csv(args.output_dir / "JCI_p53_target_recovery.tsv", sep="\t", index=False)

    shared_tfs = activities.pivot(index="contrast", columns="tf", values="ulm_score")
    p14_mean = shared_tfs.loc[["P14__Cap", "P14__Cap-a"]].mean(axis=0)
    ext_pivot = external.pivot(index="contrast", columns="tf", values="ulm_score")
    external_corr: list[dict[str, object]] = []
    for contrast, row in ext_pivot.iterrows():
        shared = p14_mean.index.intersection(row.dropna().index)
        external_corr.append({
            "comparison": f"GSE151974_P14_capillary_mean_vs_{contrast}",
            "n_tfs": len(shared),
            "pearson": float(pearsonr(p14_mean[shared], row[shared]).statistic),
            "spearman": float(spearmanr(p14_mean[shared], row[shared]).statistic),
        })
    pd.DataFrame(external_corr).to_csv(
        args.output_dir / "external_tf_activity_concordance.tsv", sep="\t", index=False
    )

    external_wide = external.pivot(index="tf", columns="contrast", values="ulm_score")
    external_fdr_wide = external.pivot(index="tf", columns="contrast", values="fdr_bh")
    for contrast in external_wide.columns:
        candidates[f"external_{contrast}_ulm"] = candidates["tf"].map(
            external_wide[contrast]
        )
        candidates[f"external_{contrast}_fdr"] = candidates["tf"].map(
            external_fdr_wide[contrast]
        )
    hyper_col = "external_JCI_hyperoxia_vs_room_air_ulm"
    hyper_fdr_col = "external_JCI_hyperoxia_vs_room_air_fdr"
    candidates["external_hyperoxia_direction_match"] = (
        np.sign(candidates["median_ulm_score"])
        == np.sign(candidates[hyper_col])
    )
    candidates["network_evidence_class"] = np.select(
        [
            (candidates["n_targets"] >= 100)
            & (candidates["significant_contrasts"] >= 4)
            & (candidates["direction_consistency"] >= 0.8)
            & candidates["external_hyperoxia_direction_match"]
            & (candidates[hyper_fdr_col] < 0.05),
            (candidates["n_targets"] >= 50)
            & (candidates["significant_contrasts"] >= 3)
            & (candidates["direction_consistency"] >= 0.8),
        ],
        ["high", "moderate"],
        default="screening_only",
    )
    evidence_order = {"high": 0, "moderate": 1, "screening_only": 2}
    candidates["_evidence_order"] = candidates["network_evidence_class"].map(
        evidence_order
    )
    candidates = candidates.sort_values(
        ["_evidence_order", "priority_score"], ascending=[True, False]
    ).drop(columns="_evidence_order")
    candidates.to_csv(args.output_dir / "regulator_priority_summary.tsv", sep="\t", index=False)

    jci_gene_effects = external_effects.T.reset_index(names="gene")
    signature_gene_calibration = signature.merge(
        jci_gene_effects, on="gene", how="left"
    )
    signature_gene_calibration["p53_null_reverses_hyperoxia"] = (
        signature_gene_calibration["JCI_hyperoxia_vs_room_air"]
        * signature_gene_calibration["JCI_p53_null_vs_control_hyperoxia"] < 0
    )
    signature_gene_calibration["p53_ec_deletion_reverses_hyperoxia"] = (
        signature_gene_calibration["JCI_hyperoxia_vs_room_air"]
        * signature_gene_calibration["JCI_p53_EC_deletion_vs_control_hyperoxia"] < 0
    )
    signature_gene_calibration.to_csv(
        args.output_dir / "JCI_signature_gene_genetic_calibration.tsv",
        sep="\t",
        index=False,
    )
    plot_results(activities, external, candidates, args.output_dir / "regulon_virtual_perturbation_summary")

    audit = {
        "network": str(args.network),
        "network_sha256": sha256(args.network),
        "network_edges": int(len(net)),
        "network_tfs": int(net["source"].nunique()),
        "contrasts": effects.index.tolist(),
        "evaluation_unit": "animal-balanced cell-type means",
        "external_tables": {
            str(path): sha256(path)
            for path in (args.jci_hyperoxia, args.jci_p53_null, args.jci_p53_ec)
        },
        "decoupler_version": dc.__version__,
        "top_candidates": candidates.head(15)["tf"].tolist(),
        "high_confidence_regulators": candidates.loc[
            candidates["network_evidence_class"].eq("high"), "tf"
        ].tolist(),
        "p53_null_signature_reversal": {
            "n_evaluable": int(
                signature_gene_calibration[
                    "JCI_p53_null_vs_control_hyperoxia"
                ].notna().sum()
            ),
            "n_opposite_direction": int(
                signature_gene_calibration["p53_null_reverses_hyperoxia"].sum()
            ),
        },
    }
    (args.output_dir / "regulon_virtual_perturbation_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print("\nTop regulator candidates:")
    print(candidates.head(20).to_string(index=False))
    print("\nP53 genetic calibration:")
    print(p53_rank.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
