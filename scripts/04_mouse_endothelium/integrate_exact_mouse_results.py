"""Integrate the formal mouse bulk and animal-level pseudobulk analyses.

This submission-facing script is intentionally self-contained. It reads only
the outputs produced by ``GSE216046_DESeq2.R`` and
``GSE151974_edgeR_pseudobulk.R`` in the current workflow. Legacy Python bulk
sensitivity results and previously published differential-state workbooks are
not required.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXACT_DIR = PROJECT_ROOT / "07_results" / "mouse_endothelium_exact"
PRIMARY_MODEL = "age_by_oxygen_P14_contrast"


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def bh_adjust(values: pd.Series) -> np.ndarray:
    """Benjamini-Hochberg adjustment while preserving missing values."""
    array = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    adjusted = np.full(array.shape, np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(array))
    if not len(finite_indices):
        return adjusted
    finite = array[finite_indices]
    order = np.argsort(finite)
    ranked = finite[order]
    scaled = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    scaled = np.minimum.accumulate(scaled[::-1])[::-1]
    scaled = np.clip(scaled, 0.0, 1.0)
    restored = np.empty_like(scaled)
    restored[order] = scaled
    adjusted[finite_indices] = restored
    return adjusted


def log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
    )


def hypergeometric_upper_tail(
    overlap: int,
    population: int,
    successes_population: int,
    draws: int,
) -> float:
    upper = min(successes_population, draws)
    logs = [
        log_comb(successes_population, k)
        + log_comb(population - successes_population, draws - k)
        - log_comb(population, draws)
        for k in range(overlap, upper + 1)
    ]
    finite = [value for value in logs if math.isfinite(value)]
    if not finite:
        return float("nan")
    maximum = max(finite)
    return min(
        1.0,
        math.exp(maximum)
        * sum(math.exp(value - maximum) for value in finite),
    )


def cross_dataset_replication(
    bulk: pd.DataFrame,
    single_cell: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare bulk and subtype-level effects under the frozen screen."""
    bulk_lookup = bulk.drop_duplicates("gene").set_index("gene")
    bulk_universe = set(bulk_lookup.index.astype(str))
    overlap_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for cell_type, subtype_frame in single_cell.groupby(
        "CellType",
        sort=False,
    ):
        subtype_lookup = (
            subtype_frame.drop_duplicates("gene").set_index("gene")
        )
        universe = sorted(
            bulk_universe.intersection(subtype_lookup.index.astype(str))
        )
        bulk_shared = bulk_lookup.loc[universe]
        subtype_shared = subtype_lookup.loc[universe]
        bulk_significant = set(
            bulk_shared.index[
                bulk_shared["fdr_bh"].lt(0.05)
                & bulk_shared["log2FC_hyperoxia_vs_air"].abs().ge(1)
            ]
        )
        subtype_significant = set(
            subtype_shared.index[
                subtype_shared["p_adj.loc__Hyperoxia"].lt(0.05)
                & subtype_shared["logFC__Hyperoxia"].abs().ge(1)
            ]
        )
        overlap = sorted(bulk_significant & subtype_significant)
        concordant = {
            gene
            for gene in overlap
            if np.sign(
                bulk_shared.loc[gene, "log2FC_hyperoxia_vs_air"]
            )
            == np.sign(subtype_shared.loc[gene, "logFC__Hyperoxia"])
        }
        summary_rows.append(
            {
                "CellType": cell_type,
                "population_genes_shared": len(universe),
                "bulk_significant_genes": len(bulk_significant),
                "single_cell_significant_genes": len(
                    subtype_significant
                ),
                "significant_overlap_genes": len(overlap),
                "direction_concordant_overlap_genes": len(concordant),
                "direction_concordance_fraction": (
                    len(concordant) / len(overlap)
                    if overlap
                    else np.nan
                ),
                "hypergeometric_p_value": hypergeometric_upper_tail(
                    len(overlap),
                    len(universe),
                    len(subtype_significant),
                    len(bulk_significant),
                ),
                "spearman_all_shared_gene_effects": float(
                    bulk_shared["log2FC_hyperoxia_vs_air"]
                    .rank()
                    .corr(subtype_shared["logFC__Hyperoxia"].rank())
                ),
            }
        )
        for gene in overlap:
            overlap_rows.append(
                {
                    "gene": gene,
                    "CellType": cell_type,
                    "bulk_log2FC_hyperoxia_vs_air": bulk_shared.loc[
                        gene,
                        "log2FC_hyperoxia_vs_air",
                    ],
                    "bulk_p_value": bulk_shared.loc[gene, "p_value"],
                    "bulk_fdr_bh": bulk_shared.loc[gene, "fdr_bh"],
                    "single_cell_logFC_hyperoxia": subtype_shared.loc[
                        gene,
                        "logFC__Hyperoxia",
                    ],
                    "single_cell_p_value": subtype_shared.loc[
                        gene,
                        "p_val__Hyperoxia",
                    ],
                    "single_cell_local_fdr": subtype_shared.loc[
                        gene,
                        "p_adj.loc__Hyperoxia",
                    ],
                    "direction_concordant": gene in concordant,
                }
            )

    summary = pd.DataFrame(summary_rows)
    summary["hypergeometric_fdr_bh"] = bh_adjust(
        summary["hypergeometric_p_value"]
    )
    overlap_table = pd.DataFrame(overlap_rows)
    if len(overlap_table):
        breadth = (
            overlap_table.loc[overlap_table["direction_concordant"]]
            .groupby("gene")["CellType"]
            .nunique()
            .rename("concordant_endothelial_subtypes")
        )
        overlap_table = overlap_table.merge(
            breadth,
            on="gene",
            how="left",
        )
        overlap_table["concordant_endothelial_subtypes"] = (
            overlap_table["concordant_endothelial_subtypes"]
            .fillna(0)
            .astype(int)
        )
        overlap_table = overlap_table.sort_values(
            [
                "concordant_endothelial_subtypes",
                "bulk_fdr_bh",
                "single_cell_local_fdr",
            ],
            ascending=[False, True, True],
        )
    return summary, overlap_table


def direction_concordant_genes(
    table: pd.DataFrame,
    model: str,
    cell_types: set[str] | None = None,
    minimum_breadth: int = 1,
) -> set[str]:
    subset = table.loc[
        table["model"].eq(model) & table["direction_concordant"]
    ].copy()
    if cell_types is not None:
        subset = subset.loc[subset["CellType"].isin(cell_types)]
    breadth = subset.groupby("gene")["CellType"].nunique()
    return set(breadth.index[breadth.ge(minimum_breadth)])


def main() -> int:
    EXACT_DIR.mkdir(parents=True, exist_ok=True)
    deseq = pd.read_csv(
        EXACT_DIR / "GSE216046_DESeq2_all_genes.tsv",
        sep="\t",
    )
    exact_bulk = deseq.rename(
        columns={
            "log2FoldChange_mle": "log2FC_hyperoxia_vs_air",
            "pvalue": "p_value",
            "padj": "fdr_bh",
        }
    )[
        [
            "gene",
            "baseMean",
            "log2FC_hyperoxia_vs_air",
            "p_value",
            "fdr_bh",
            "log2FoldChange_normal_shrunk",
            "lfcSE_normal_shrunk",
            "pvalue_lfcThreshold1",
            "padj_lfcThreshold1",
        ]
    ].copy()
    exact_bulk["passes_fdr05_abs_log2fc1"] = (
        exact_bulk["fdr_bh"].lt(0.05)
        & exact_bulk["log2FC_hyperoxia_vs_air"].abs().ge(1)
    )

    raw_dsa = pd.read_csv(
        EXACT_DIR / "GSE151974_raw_pseudobulk_edgeR_all_results.tsv.gz",
        sep="\t",
    )
    raw_replication_frames: list[pd.DataFrame] = []
    raw_convergent_frames: list[pd.DataFrame] = []
    for model_name, model_frame in raw_dsa.groupby(
        "model",
        sort=False,
    ):
        standardized = model_frame.rename(
            columns={
                "logFC": "logFC__Hyperoxia",
                "PValue": "p_val__Hyperoxia",
                "FDR": "p_adj.loc__Hyperoxia",
            }
        ).copy()
        replication, convergent = cross_dataset_replication(
            exact_bulk,
            standardized,
        )
        replication.insert(0, "model", model_name)
        convergent.insert(0, "model", model_name)
        raw_replication_frames.append(replication)
        raw_convergent_frames.append(convergent)
    raw_replication = pd.concat(raw_replication_frames, ignore_index=True)
    raw_convergent = pd.concat(raw_convergent_frames, ignore_index=True)

    primary_replication = raw_replication.loc[
        raw_replication["model"].eq(PRIMARY_MODEL)
    ].copy()
    primary_convergent = raw_convergent.loc[
        raw_convergent["model"].eq(PRIMARY_MODEL)
    ].copy()

    exact_bulk_treat = exact_bulk.copy()
    exact_bulk_treat["p_value"] = exact_bulk_treat[
        "pvalue_lfcThreshold1"
    ]
    exact_bulk_treat["fdr_bh"] = exact_bulk_treat[
        "padj_lfcThreshold1"
    ]
    treat_replication_frames: list[pd.DataFrame] = []
    treat_convergent_frames: list[pd.DataFrame] = []
    for model_name, model_frame in raw_dsa.groupby(
        "model",
        sort=False,
    ):
        standardized = model_frame.rename(
            columns={
                "logFC": "logFC__Hyperoxia",
                "PValue_treat_lfc1": "p_val__Hyperoxia",
                "FDR_treat_lfc1": "p_adj.loc__Hyperoxia",
            }
        ).copy()
        replication, convergent = cross_dataset_replication(
            exact_bulk_treat,
            standardized,
        )
        replication.insert(0, "model", model_name)
        convergent.insert(0, "model", model_name)
        treat_replication_frames.append(replication)
        treat_convergent_frames.append(convergent)
    treat_replication = pd.concat(
        treat_replication_frames,
        ignore_index=True,
    )
    treat_convergent = pd.concat(
        treat_convergent_frames,
        ignore_index=True,
    )

    primary_all = direction_concordant_genes(
        raw_convergent,
        PRIMARY_MODEL,
    )
    primary_multi_subtype = direction_concordant_genes(
        raw_convergent,
        PRIMARY_MODEL,
        minimum_breadth=3,
    )
    nonrare_types = {"Cap", "Cap-a", "Art"}
    primary_nonrare = direction_concordant_genes(
        raw_convergent,
        PRIMARY_MODEL,
        cell_types=nonrare_types,
    )
    primary_nonrare_breadth2 = direction_concordant_genes(
        raw_convergent,
        PRIMARY_MODEL,
        cell_types=nonrare_types,
        minimum_breadth=2,
    )

    sensitivity_models = [
        "P14_oxygen_only",
        "age_adjusted_all_ages",
        "age_by_oxygen_P14_contrast_min10_cells",
        "P14_oxygen_only_min10_cells",
        "P14_oxygen_only_min20_cells",
    ]
    sensitivity_rows: list[dict[str, object]] = []
    for model_name in sensitivity_models:
        model_genes = direction_concordant_genes(
            raw_convergent,
            model_name,
        )
        union = primary_all | model_genes
        sensitivity_rows.append(
            {
                "comparison": model_name,
                "primary_genes": len(primary_all),
                "comparison_genes": len(model_genes),
                "overlap_genes": len(primary_all & model_genes),
                "primary_retention_fraction": (
                    len(primary_all & model_genes) / len(primary_all)
                    if primary_all
                    else np.nan
                ),
                "jaccard": (
                    len(primary_all & model_genes) / len(union)
                    if union
                    else np.nan
                ),
            }
        )
    nonrare_union = primary_all | primary_nonrare
    sensitivity_rows.append(
        {
            "comparison": "exclude_Vein_and_Lymph_any_of_Cap_Cap-a_Art",
            "primary_genes": len(primary_all),
            "comparison_genes": len(primary_nonrare),
            "overlap_genes": len(primary_all & primary_nonrare),
            "primary_retention_fraction": (
                len(primary_all & primary_nonrare) / len(primary_all)
                if primary_all
                else np.nan
            ),
            "jaccard": (
                len(primary_all & primary_nonrare) / len(nonrare_union)
                if nonrare_union
                else np.nan
            ),
        }
    )
    sensitivity_summary = pd.DataFrame(sensitivity_rows)

    membership = (
        primary_convergent.loc[
            primary_convergent["direction_concordant"]
        ]
        .groupby("gene")
        .agg(
            concordant_endothelial_subtypes=("CellType", "nunique"),
            replicated_subtype_names=(
                "CellType",
                lambda values: ";".join(sorted(set(values))),
            ),
            bulk_log2FC_hyperoxia_vs_air=(
                "bulk_log2FC_hyperoxia_vs_air",
                "first",
            ),
            bulk_fdr_bh=("bulk_fdr_bh", "first"),
        )
        .reset_index()
    )
    membership["included_high_confidence_ge3"] = membership[
        "concordant_endothelial_subtypes"
    ].ge(3)
    membership["included_after_excluding_vein_lymph_ge2_of3"] = (
        membership["gene"].isin(primary_nonrare_breadth2)
    )

    write_tsv(
        raw_replication,
        EXACT_DIR
        / "raw_rerun_mouse_endothelial_cross_dataset_replication_summary.tsv",
    )
    write_tsv(
        raw_convergent,
        EXACT_DIR / "raw_rerun_mouse_endothelial_convergent_genes.tsv",
    )
    write_tsv(
        primary_replication,
        EXACT_DIR
        / "mouse_endothelial_cross_dataset_replication_summary.tsv",
    )
    write_tsv(
        primary_convergent,
        EXACT_DIR / "mouse_endothelial_convergent_genes.tsv",
    )
    write_tsv(
        treat_replication,
        EXACT_DIR
        / "mouse_endothelial_cross_dataset_replication_glmTreat_lfc1.tsv",
    )
    write_tsv(
        treat_convergent,
        EXACT_DIR
        / "mouse_endothelial_convergent_genes_glmTreat_lfc1.tsv",
    )
    write_tsv(
        sensitivity_summary,
        EXACT_DIR / "mouse_signature_model_sensitivity.tsv",
    )
    write_tsv(
        membership,
        EXACT_DIR / "mouse_high_confidence_signature_membership.tsv",
    )

    summary = {
        "analysis_version": "portable_exact_DESeq2_edgeR_1.1",
        "genes_after_cpm_filter": int(len(exact_bulk)),
        "genes_fdr05": int(exact_bulk["fdr_bh"].lt(0.05).sum()),
        "genes_fdr05_abs_mle_log2fc1": int(
            exact_bulk["passes_fdr05_abs_log2fc1"].sum()
        ),
        "genes_lfcThreshold1_fdr05": int(
            exact_bulk["padj_lfcThreshold1"].lt(0.05).sum()
        ),
        "cross_dataset": {
            "primary_GSE151974_model": PRIMARY_MODEL,
            "replication_by_endothelial_subtype": (
                primary_replication.to_dict(orient="records")
            ),
            "convergent_rows": int(len(primary_convergent)),
            "unique_overlap_genes": int(
                primary_convergent["gene"].nunique()
            ),
            "unique_direction_concordant_genes": int(
                primary_convergent.loc[
                    primary_convergent["direction_concordant"],
                    "gene",
                ].nunique()
            ),
            "multi_subtype_genes_ge3": len(primary_multi_subtype),
            "nonrare_endothelial_genes_any_of_Cap_Cap-a_Art": len(
                primary_nonrare
            ),
            "nonrare_endothelial_genes_ge2_of3": len(
                primary_nonrare_breadth2
            ),
        },
        "method_notes": [
            "The integration reads only outputs generated by the current DESeq2 and edgeR workflow.",
            "GSE216046 was analyzed with DESeq2 using mouse as the unit of inference.",
            "GSE151974 was aggregated to animal-by-endothelial-subtype pseudobulks.",
            "The primary subtype model uses categorical age, oxygen and their interaction and tests the P14 hyperoxia-versus-normoxia contrast.",
            "FDR<0.05 plus |observed log2FC|>=1 is the descriptive replication screen; DESeq2 lfcThreshold=1 and edgeR glmTreat(lfc=1) are effect-threshold sensitivities.",
        ],
    }
    (EXACT_DIR / "exact_mouse_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "genes_after_cpm_filter": summary[
                    "genes_after_cpm_filter"
                ],
                "unique_direction_concordant_genes": summary[
                    "cross_dataset"
                ]["unique_direction_concordant_genes"],
                "multi_subtype_genes_ge3": summary["cross_dataset"][
                    "multi_subtype_genes_ge3"
                ],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
