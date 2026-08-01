"""Mouse-signature definition and human-composition sensitivity analyses.

The script does not create new donor-level replication. It quantifies how the
human directional projection changes under stricter mouse effect thresholds,
alternative endothelial-subtype breadth rules, exclusion of abCap and systemic
venous EC from the all-endothelial aggregate, and equal weighting of the five
subtypes recovered with at least ten cells in every disease donor.
"""

from __future__ import annotations

import json
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MOUSE_DIR = PROJECT_ROOT / "07_results" / "mouse_endothelium_exact"
CROSS_DIR = PROJECT_ROOT / "07_results" / "cross_species"
HUMAN_DIR = PROJECT_ROOT / "07_results" / "human_endothelial_subtypes"
OUTPUT_DIR = PROJECT_ROOT / "07_results" / "sensitivity_analyses"

BPD = ["BPD 1", "BPD 2"]
PH = ["BPD+PH 1", "BPD+PH 2"]
SAMPLES = BPD + PH
COMMON_FIVE = [
    "gCap",
    "aCap",
    "Arterial EC",
    "Pulmonary venous EC",
    "Lymphatic",
]


def oriented_score(
    expression: pd.DataFrame,
    directions: pd.Series,
) -> tuple[pd.Series, int]:
    expression = expression.loc[SAMPLES, directions.index]
    standard_deviation = expression.std(axis=0, ddof=0)
    variable = standard_deviation.index[standard_deviation.gt(0)]
    if not len(variable):
        return pd.Series(np.nan, index=SAMPLES), 0
    standardized = expression[variable].sub(
        expression[variable].mean(axis=0),
        axis=1,
    ).div(standard_deviation[variable], axis=1)
    scores = standardized.mul(directions[variable], axis=1).mean(axis=1)
    return scores, int(len(variable))


def contrast_summary(scores: pd.Series) -> dict[str, float | bool]:
    scores = scores.reindex(SAMPLES)
    full = float(scores[PH].mean() - scores[BPD].mean())
    pairwise = np.asarray(
        [
            scores[ph_sample] - scores[bpd_sample]
            for ph_sample, bpd_sample in product(PH, BPD)
        ],
        dtype=float,
    )
    leave_one_out = []
    for removed in SAMPLES:
        remaining_bpd = [sample for sample in BPD if sample != removed]
        remaining_ph = [sample for sample in PH if sample != removed]
        leave_one_out.append(
            float(
                scores[remaining_ph].mean()
                - scores[remaining_bpd].mean()
            )
        )
    leave_one_out = np.asarray(leave_one_out, dtype=float)
    return {
        "oriented_effect": full,
        "pairwise_min": float(pairwise.min()),
        "pairwise_median": float(np.median(pairwise)),
        "pairwise_max": float(pairwise.max()),
        "all_pairwise_positive": bool(np.all(pairwise > 0)),
        "leave_one_donor_out_min": float(leave_one_out.min()),
        "leave_one_donor_out_median": float(np.median(leave_one_out)),
        "leave_one_donor_out_max": float(leave_one_out.max()),
        "all_leave_one_donor_out_positive": bool(
            np.all(leave_one_out > 0)
        ),
    }


def eligible_genes(frame: pd.DataFrame) -> pd.Index:
    detection = frame.groupby("gene").agg(
        total_raw_count=("raw_count", "sum"),
        donors_detected=("detected", "sum"),
    )
    return detection.index[
        detection["total_raw_count"].ge(10)
        & detection["donors_detected"].ge(2)
    ]


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    mapping = pd.read_csv(
        CROSS_DIR / "formal_one_to_one_mouse_signature.tsv",
        sep="\t",
    )
    mapping = mapping.dropna(subset=["human_gene"]).drop_duplicates("gene")
    mapping["mouse_direction"] = np.sign(
        mapping["mouse_bulk_log2FC"]
    ).astype(int)
    mapping = mapping.set_index("gene", drop=False)

    bulk = pd.read_csv(
        MOUSE_DIR / "GSE216046_DESeq2_all_genes.tsv",
        sep="\t",
    ).set_index("gene")
    edge = pd.read_csv(
        MOUSE_DIR / "GSE151974_raw_pseudobulk_edgeR_all_results.tsv.gz",
        sep="\t",
    )
    edge = edge.loc[
        edge["model"].eq("age_by_oxygen_P14_contrast")
    ].copy()
    edge["bulk_direction"] = np.sign(
        edge["gene"].map(bulk["log2FoldChange_mle"])
    )
    edge["treat_same_direction"] = (
        edge["FDR_treat_lfc1"].lt(0.05)
        & np.sign(edge["logFC"]).eq(edge["bulk_direction"])
    )
    edge_treat = (
        edge.groupby("gene", as_index=True)
        .agg(
            treat_concordant_subtypes=(
                "treat_same_direction",
                "sum",
            ),
            treat_concordant_subtype_names=(
                "CellType",
                lambda values: ";".join(
                    sorted(
                        edge.loc[
                            values.index[
                                edge.loc[values.index, "treat_same_direction"]
                            ],
                            "CellType",
                        ].unique()
                    )
                ),
            ),
        )
    )

    mouse_audit = mapping[
        [
            "gene",
            "human_gene",
            "mouse_bulk_log2FC",
            "concordant_endothelial_subtypes",
            "replicated_subtype_names",
        ]
    ].copy()
    mouse_audit["bulk_formal_lfc1_fdr05"] = mouse_audit["gene"].map(
        bulk["padj_lfcThreshold1"].lt(0.05)
    ).fillna(False)
    mouse_audit["edgeR_treat_concordant_subtypes"] = (
        mouse_audit["gene"]
        .map(edge_treat["treat_concordant_subtypes"])
        .fillna(0)
        .astype(int)
    )
    mouse_audit["edgeR_treat_concordant_subtype_names"] = (
        mouse_audit["gene"]
        .map(edge_treat["treat_concordant_subtype_names"])
        .fillna("")
    )
    subtype_sets = mouse_audit["replicated_subtype_names"].fillna("").map(
        lambda value: set(value.split(";")) if value else set()
    )
    mouse_audit["replicated_in_both_Cap_and_Cap_a"] = subtype_sets.map(
        lambda values: {"Cap", "Cap-a"}.issubset(values)
    )
    mouse_audit["original_multi_subtype_ge3"] = mouse_audit[
        "concordant_endothelial_subtypes"
    ].ge(3)
    mouse_audit["dual_threshold_any_treat"] = (
        mouse_audit["original_multi_subtype_ge3"]
        & mouse_audit["bulk_formal_lfc1_fdr05"]
        & mouse_audit["edgeR_treat_concordant_subtypes"].ge(1)
    )
    mouse_audit["dual_threshold_treat_ge3"] = (
        mouse_audit["original_multi_subtype_ge3"]
        & mouse_audit["bulk_formal_lfc1_fdr05"]
        & mouse_audit["edgeR_treat_concordant_subtypes"].ge(3)
    )
    mouse_audit.to_csv(
        OUTPUT_DIR / "mouse_33_gene_effect_threshold_audit.tsv",
        sep="\t",
        index=False,
    )

    definitions: dict[str, set[str]] = {
        "mouse_subtypes_ge2": set(
            mouse_audit.loc[
                mouse_audit["concordant_endothelial_subtypes"].ge(2),
                "gene",
            ]
        ),
        "mouse_subtypes_ge3": set(
            mouse_audit.loc[
                mouse_audit["concordant_endothelial_subtypes"].ge(3),
                "gene",
            ]
        ),
        "mouse_subtypes_ge4": set(
            mouse_audit.loc[
                mouse_audit["concordant_endothelial_subtypes"].ge(4),
                "gene",
            ]
        ),
        "mouse_subtypes_ge5": set(
            mouse_audit.loc[
                mouse_audit["concordant_endothelial_subtypes"].ge(5),
                "gene",
            ]
        ),
        "Cap_and_Cap_a_both": set(
            mouse_audit.loc[
                mouse_audit["replicated_in_both_Cap_and_Cap_a"],
                "gene",
            ]
        ),
        "ge3_bulk_formal_lfc1": set(
            mouse_audit.loc[
                mouse_audit["original_multi_subtype_ge3"]
                & mouse_audit["bulk_formal_lfc1_fdr05"],
                "gene",
            ]
        ),
        "ge3_edgeR_treat_any": set(
            mouse_audit.loc[
                mouse_audit["original_multi_subtype_ge3"]
                & mouse_audit["edgeR_treat_concordant_subtypes"].ge(1),
                "gene",
            ]
        ),
        "ge3_edgeR_treat_ge3": set(
            mouse_audit.loc[
                mouse_audit["original_multi_subtype_ge3"]
                & mouse_audit["edgeR_treat_concordant_subtypes"].ge(3),
                "gene",
            ]
        ),
        "ge3_dual_threshold_any_treat": set(
            mouse_audit.loc[
                mouse_audit["dual_threshold_any_treat"],
                "gene",
            ]
        ),
        "ge3_dual_threshold_treat_ge3": set(
            mouse_audit.loc[
                mouse_audit["dual_threshold_treat_ge3"],
                "gene",
            ]
        ),
    }

    definition_rows = []
    for name, mouse_genes in definitions.items():
        projected = mapping.index.intersection(list(mouse_genes))
        definition_rows.append(
            {
                "signature": name,
                "mouse_genes": len(mouse_genes),
                "one_to_one_human_orthologs": len(projected),
                "definition": {
                    "mouse_subtypes_ge2": "replicated in at least 2 of 5 mouse endothelial subtypes",
                    "mouse_subtypes_ge3": "replicated in at least 3 of 5 mouse endothelial subtypes",
                    "mouse_subtypes_ge4": "replicated in at least 4 of 5 mouse endothelial subtypes",
                    "mouse_subtypes_ge5": "replicated in all 5 mouse endothelial subtypes",
                    "Cap_and_Cap_a_both": "replicated in both mouse Cap and Cap-a",
                    "ge3_bulk_formal_lfc1": "ge3 set also passing DESeq2 lfcThreshold=1 at FDR<0.05",
                    "ge3_edgeR_treat_any": "ge3 set passing edgeR glmTreat in at least 1 concordant subtype",
                    "ge3_edgeR_treat_ge3": "ge3 set passing edgeR glmTreat in at least 3 concordant subtypes",
                    "ge3_dual_threshold_any_treat": "ge3 set passing DESeq2 lfcThreshold=1 and edgeR glmTreat in at least 1 subtype",
                    "ge3_dual_threshold_treat_ge3": "ge3 set passing DESeq2 lfcThreshold=1 and edgeR glmTreat in at least 3 subtypes",
                }[name],
            }
        )
    pd.DataFrame(definition_rows).to_csv(
        OUTPUT_DIR / "mouse_signature_definition_sensitivity.tsv",
        sep="\t",
        index=False,
    )

    expression = pd.read_csv(
        HUMAN_DIR / "GSE275938_subtype_locked_gene_expression.tsv.gz",
        sep="\t",
        low_memory=False,
    )
    cell_coverage = pd.read_csv(
        HUMAN_DIR / "GSE275938_disease_endothelial_subtype_cell_coverage.tsv",
        sep="\t",
    ).set_index("subtype")
    common_expression = pd.read_csv(
        HUMAN_DIR
        / "GSE275938_common_five_subtypes_aggregate_expression.tsv.gz",
        sep="\t",
        low_memory=False,
    )
    common_meta = pd.read_csv(
        HUMAN_DIR / "GSE275938_common_five_subtypes_aggregate_metadata.tsv",
        sep="\t",
    )

    summary_rows: list[dict] = []
    donor_rows: list[dict] = []

    def score_frame(
        frame: pd.DataFrame,
        analysis_population: str,
        minimum_cells: int,
        included_subtypes: str,
        excluded_subtypes: str,
    ):
        eligible = eligible_genes(frame)
        matrix = frame.pivot(
            index="sample",
            columns="gene",
            values="log2CPM_TMM",
        ).reindex(SAMPLES)
        for signature_name, mouse_genes in definitions.items():
            projected_mouse = mapping.index.intersection(list(mouse_genes))
            projected = mapping.loc[projected_mouse]
            requested_human = pd.Index(projected["human_gene"])
            scored_human = requested_human.intersection(eligible).intersection(
                matrix.columns
            )
            direction_by_human = (
                projected.set_index("human_gene")["mouse_direction"]
                .reindex(scored_human)
                .astype(int)
            )
            scores, variable = oriented_score(
                matrix.loc[:, scored_human],
                direction_by_human,
            )
            for sample in SAMPLES:
                donor_rows.append(
                    {
                        "analysis_population": analysis_population,
                        "sample": sample,
                        "condition": "BPD+PH" if sample in PH else "BPD",
                        "signature": signature_name,
                        "minimum_cells_across_donors": minimum_cells,
                        "genes_requested": len(projected_mouse),
                        "genes_eligible": len(scored_human),
                        "genes_variable": variable,
                        "oriented_score": float(scores[sample]),
                    }
                )
            summary_rows.append(
                {
                    "analysis_population": analysis_population,
                    "signature": signature_name,
                    "minimum_cells_across_donors": minimum_cells,
                    "included_subtypes": included_subtypes,
                    "excluded_subtypes": excluded_subtypes,
                    "eligible_for_primary_subtype_reporting": (
                        minimum_cells >= 10
                    ),
                    "genes_requested": len(projected_mouse),
                    "genes_eligible": len(scored_human),
                    "genes_variable": variable,
                    **contrast_summary(scores),
                }
            )

    for subtype in expression["subtype"].drop_duplicates():
        score_frame(
            expression.loc[expression["subtype"].eq(subtype)].copy(),
            subtype,
            int(cell_coverage.loc[subtype, "minimum_cells_across_donors"]),
            subtype,
            "",
        )
    score_frame(
        common_expression,
        "All Endothelial excluding abCap and Systemic venous EC",
        int(common_meta["cells"].min()),
        ";".join(COMMON_FIVE),
        "abCap;Systemic venous EC",
    )

    donor_scores = pd.DataFrame(donor_rows)
    summary = pd.DataFrame(summary_rows)

    # Composition-balanced score: each common subtype contributes equal weight
    # after its own within-subtype standardization and direction orientation.
    balanced_rows = []
    balanced_donors = []
    for signature_name in definitions:
        selected = donor_scores.loc[
            donor_scores["analysis_population"].isin(COMMON_FIVE)
            & donor_scores["signature"].eq(signature_name)
        ].copy()
        wide = selected.pivot(
            index="sample",
            columns="analysis_population",
            values="oriented_score",
        ).reindex(index=SAMPLES, columns=COMMON_FIVE)
        if wide.isna().any().any():
            raise RuntimeError(
                f"Incomplete balanced score for {signature_name}"
            )
        scores = wide.mean(axis=1)
        eligible_counts = selected.groupby("analysis_population")[
            "genes_eligible"
        ].first()
        for sample in SAMPLES:
            balanced_donors.append(
                {
                    "analysis_population": "Composition-balanced common five subtypes",
                    "sample": sample,
                    "condition": "BPD+PH" if sample in PH else "BPD",
                    "signature": signature_name,
                    "minimum_cells_across_donors": int(
                        cell_coverage.loc[COMMON_FIVE, "minimum_cells_across_donors"].min()
                    ),
                    "genes_requested": int(
                        selected["genes_requested"].max()
                    ),
                    "genes_eligible": (
                        f"{int(eligible_counts.min())}-"
                        f"{int(eligible_counts.max())}"
                    ),
                    "genes_variable": "subtype-specific",
                    "oriented_score": float(scores[sample]),
                }
            )
        balanced_rows.append(
            {
                "analysis_population": "Composition-balanced common five subtypes",
                "signature": signature_name,
                "minimum_cells_across_donors": int(
                    cell_coverage.loc[COMMON_FIVE, "minimum_cells_across_donors"].min()
                ),
                "included_subtypes": ";".join(COMMON_FIVE),
                "excluded_subtypes": "abCap;Systemic venous EC",
                "eligible_for_primary_subtype_reporting": True,
                "genes_requested": int(selected["genes_requested"].max()),
                "genes_eligible": (
                    f"{int(eligible_counts.min())}-"
                    f"{int(eligible_counts.max())}"
                ),
                "genes_variable": "subtype-specific",
                **contrast_summary(scores),
            }
        )
    summary = pd.concat(
        [summary, pd.DataFrame(balanced_rows)],
        ignore_index=True,
    )
    donor_scores = pd.concat(
        [donor_scores, pd.DataFrame(balanced_donors)],
        ignore_index=True,
    )
    summary.to_csv(
        OUTPUT_DIR / "GSE275938_signature_definition_and_composition_sensitivity.tsv",
        sep="\t",
        index=False,
    )
    donor_scores.to_csv(
        OUTPUT_DIR / "GSE275938_signature_definition_donor_scores.tsv",
        sep="\t",
        index=False,
    )

    random_summary = pd.read_csv(
        HUMAN_DIR / "GSE275938_expression_matched_random_set_summary.tsv",
        sep="\t",
    )
    random_six = random_summary.loc[
        random_summary["signature"].eq("mouse_subtypes_ge3")
    ].copy()
    if len(random_six) != 6:
        raise AssertionError(
            f"Expected six tested populations, found {len(random_six)}"
        )
    random_six.to_csv(
        OUTPUT_DIR / "GSE275938_six_population_random_set_results.tsv",
        sep="\t",
        index=False,
    )

    high = mouse_audit.loc[
        mouse_audit["original_multi_subtype_ge3"]
    ].copy()
    key_summary = {
        "multi_subtype_mouse_signature_genes": int(len(high)),
        "pass_DESeq2_formal_lfc1_FDR05": int(
            high["bulk_formal_lfc1_fdr05"].sum()
        ),
        "pass_edgeR_glmTreat_any_concordant_subtype": int(
            high["edgeR_treat_concordant_subtypes"].ge(1).sum()
        ),
        "pass_edgeR_glmTreat_at_least_2_subtypes": int(
            high["edgeR_treat_concordant_subtypes"].ge(2).sum()
        ),
        "pass_edgeR_glmTreat_at_least_3_subtypes": int(
            high["edgeR_treat_concordant_subtypes"].ge(3).sum()
        ),
        "pass_both_DESeq2_and_edgeR_any": int(
            high["dual_threshold_any_treat"].sum()
        ),
        "pass_both_DESeq2_and_edgeR_at_least_3": int(
            high["dual_threshold_treat_ge3"].sum()
        ),
        "human_cell_threshold_for_primary_subtype_reporting": 10,
        "excluded_from_primary_subtype_reporting": [
            "abCap",
            "Systemic venous EC",
        ],
        "all_endothelial_original_includes_all_seven_subtypes": True,
        "composition_restricted_aggregate_excludes": [
            "abCap",
            "Systemic venous EC",
        ],
        "composition_balanced_method": (
            "Equal arithmetic mean of donor-level oriented scores across "
            "gCap, aCap, Arterial EC, Pulmonary venous EC and Lymphatic "
            "after subtype-specific standardization."
        ),
    }
    (
        OUTPUT_DIR / "signature_and_composition_sensitivity_summary.json"
    ).write_text(
        json.dumps(key_summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(key_summary, ensure_ascii=False, indent=2))
    print(
        summary.loc[
            summary["signature"].isin(
                [
                    "mouse_subtypes_ge3",
                    "mouse_subtypes_ge4",
                    "ge3_dual_threshold_any_treat",
                    "ge3_dual_threshold_treat_ge3",
                    "Cap_and_Cap_a_both",
                ]
            ),
            [
                "analysis_population",
                "signature",
                "genes_eligible",
                "oriented_effect",
                "leave_one_donor_out_min",
            ],
        ].to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
