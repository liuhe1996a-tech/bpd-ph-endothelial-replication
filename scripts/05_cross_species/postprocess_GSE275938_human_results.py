"""Create compact, manuscript-facing human endothelial result tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULT_DIR = PROJECT_ROOT / "07_results" / "human_endothelial_subtypes"
CORE_SUBTYPES = ["All Endothelial", "gCap", "aCap"]


def main() -> int:
    gene_effects = pd.read_csv(
        RESULT_DIR / "GSE275938_subtype_locked_gene_effects.tsv.gz",
        sep="\t",
        low_memory=False,
    )
    primary = gene_effects.loc[
        gene_effects["normalization"].eq("TMM")
        & gene_effects["eligible_expression"]
        & gene_effects["mouse_concordant_endothelial_subtypes"].ge(3)
        & gene_effects["subtype"].isin(CORE_SUBTYPES)
    ].copy()
    robust_gene_names = (
        primary.groupby("gene")
        .filter(
            lambda group: (
                group["subtype"].nunique() == len(CORE_SUBTYPES)
                and group[
                    "all_four_pairwise_mouse_direction_concordant"
                ].all()
            )
        )["gene"]
        .drop_duplicates()
        .sort_values()
    )
    robust = primary.loc[primary["gene"].isin(robust_gene_names)].copy()
    effects_wide = robust.pivot(
        index=[
            "gene",
            "mouse_bulk_log2FC",
            "mouse_concordant_endothelial_subtypes",
        ],
        columns="subtype",
        values="human_BPD_PH_minus_BPD_log2CPM",
    ).reset_index()
    effects_wide.columns.name = None
    effects_wide = effects_wide.rename(
        columns={
            "All Endothelial": "human_effect_all_endothelial",
            "gCap": "human_effect_gCap",
            "aCap": "human_effect_aCap",
        }
    )
    effects_wide["minimum_human_effect_across_core_subtypes"] = effects_wide[
        [
            "human_effect_all_endothelial",
            "human_effect_gCap",
            "human_effect_aCap",
        ]
    ].min(axis=1)
    effects_wide = effects_wide.sort_values(
        [
            "mouse_concordant_endothelial_subtypes",
            "minimum_human_effect_across_core_subtypes",
        ],
        ascending=[False, False],
    )
    effects_wide.to_csv(
        RESULT_DIR / "GSE275938_capillary_robust_core_genes.tsv",
        sep="\t",
        index=False,
    )
    effects_wide.to_csv(
        RESULT_DIR
        / "GSE275938_exploratory_human_concordant_subset.tsv",
        sep="\t",
        index=False,
    )

    donor_expression = pd.read_csv(
        RESULT_DIR / "GSE275938_subtype_locked_gene_expression.tsv.gz",
        sep="\t",
        low_memory=False,
    )
    subset_donor_expression = donor_expression.loc[
        donor_expression["gene"].isin(effects_wide["gene"])
        & donor_expression["subtype"].isin(CORE_SUBTYPES),
        [
            "sample",
            "condition",
            "subtype",
            "gene",
            "raw_count",
            "detected",
            "log2CPM_TMM",
            "log2CPM_library",
        ],
    ].copy()
    subset_donor_expression["subset_definition"] = (
        "post-projection exploratory annotation: one-to-one ortholog, "
        "mouse high-confidence breadth >=3, expression-eligible, and "
        "direction-concordant in every BPD+PH-versus-BPD donor pairing "
        "within All Endothelial, gCap, and aCap"
    )
    subset_donor_expression.to_csv(
        RESULT_DIR
        / "GSE275938_exploratory_subset_donor_counts_and_log2CPM.tsv",
        sep="\t",
        index=False,
    )

    sensitivity = pd.read_csv(
        RESULT_DIR / "GSE275938_subtype_locked_signature_sensitivity.tsv",
        sep="\t",
    )
    downsampling = pd.read_csv(
        RESULT_DIR / "GSE275938_equal_cell_downsampling_summary.tsv",
        sep="\t",
    )
    random_null = pd.read_csv(
        RESULT_DIR / "GSE275938_expression_matched_random_set_summary.tsv",
        sep="\t",
    )
    compact = (
        sensitivity.loc[
            sensitivity["normalization"].eq("TMM")
            & sensitivity["signature"].eq("mouse_subtypes_ge3")
            & sensitivity["minimum_cells_across_donors"].ge(10)
        ][
            [
                "subtype",
                "minimum_cells_across_donors",
                "genes_eligible",
                "BPD_PH_minus_BPD_mean_score",
                "pairwise_positive_fraction",
                "all_four_pairwise_same_direction_as_full",
                "all_leave_one_out_same_direction_as_full",
            ]
        ]
        .merge(
            downsampling.loc[
                downsampling["signature"].eq("mouse_subtypes_ge3")
            ][
                [
                    "subtype",
                    "cells_per_donor",
                    "median_effect",
                    "sampling_p2_5",
                    "sampling_p97_5",
                    "positive_effect_fraction",
                ]
            ],
            on="subtype",
            how="left",
        )
        .merge(
            random_null.loc[
                random_null["signature"].eq("mouse_subtypes_ge3")
            ][
                [
                    "subtype",
                    "observed_effect",
                    "null_median",
                    "null_p2_5",
                    "null_p97_5",
                    "empirical_p_greater",
                    "empirical_p_two_sided_abs_effect",
                    "empirical_q_greater_BH_within_signature",
                    "empirical_q_two_sided_BH_within_signature",
                ]
            ],
            on="subtype",
            how="left",
        )
        .sort_values(
            ["minimum_cells_across_donors", "subtype"],
            ascending=[False, True],
        )
    )
    compact.to_csv(
        RESULT_DIR / "GSE275938_human_subtype_compact_evidence.tsv",
        sep="\t",
        index=False,
    )
    print(effects_wide.to_string(index=False))
    print("\nCompact evidence:")
    print(compact.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
