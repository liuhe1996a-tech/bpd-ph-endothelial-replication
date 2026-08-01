"""Donor-level human endothelial subtype and low-cell sensitivity analysis.

The four disease donors are never expanded into cell-level replicates.
Inferential p-values are intentionally omitted because BPD and BPD+PH each
contain only two donors.
"""

from __future__ import annotations

import json
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROCESSED_DIR = PROJECT_ROOT / "04_processed" / "GSE275938"
RESULT_DIR = PROJECT_ROOT / "07_results" / "human_endothelial_subtypes"
EXPRESSION_PATH = RESULT_DIR / "GSE275938_subtype_locked_gene_expression.tsv.gz"
PSEUDOBULK_METADATA_PATH = (
    PROCESSED_DIR / "GSE275938_endothelial_subtype_pseudobulk_metadata.tsv"
)
CELL_METADATA_PATH = (
    PROCESSED_DIR / "GSE275938_disease_endothelial_cell_metadata.tsv.gz"
)
CELL_COUNTS_PATH = (
    PROCESSED_DIR / "GSE275938_disease_endothelial_locked_gene_counts.npz"
)
LOCKED_EVIDENCE_PATH = (
    PROJECT_ROOT
    / "07_results"
    / "cross_species"
    / "locked_mouse_signature_human_BPD_PH_evidence.tsv"
)

BPD_SAMPLES = ["BPD 1", "BPD 2"]
PH_SAMPLES = ["BPD+PH 1", "BPD+PH 2"]
SAMPLE_ORDER = BPD_SAMPLES + PH_SAMPLES
NORMALIZATIONS = {
    "TMM": "log2CPM_TMM",
    "library_CPM": "log2CPM_library",
}
SIGNATURE_THRESHOLDS = {
    "mouse_subtypes_ge1": 1,
    "mouse_subtypes_ge2": 2,
    "mouse_subtypes_ge3": 3,
}
CELL_THRESHOLDS = [1, 5, 10, 20, 50, 100]
ALL_ENDOTHELIAL = "All Endothelial"
RNG_SEED = 20260726
DOWNSAMPLING_REPEATS = 500
RAREFACTION_REPEATS = 500


def oriented_signature_score(
    expression: pd.DataFrame,
    mouse_directions: pd.Series,
) -> tuple[pd.Series, int]:
    expression = expression.loc[:, mouse_directions.index]
    standard_deviation = expression.std(axis=0, ddof=0)
    variable_genes = standard_deviation.index[standard_deviation > 0]
    if not len(variable_genes):
        return pd.Series(np.nan, index=expression.index), 0
    z_scores = expression[variable_genes].sub(
        expression[variable_genes].mean(axis=0),
        axis=1,
    ).div(standard_deviation[variable_genes], axis=1)
    oriented = z_scores.mul(mouse_directions[variable_genes], axis=1)
    return oriented.mean(axis=1), len(variable_genes)


def donor_contrast_summary(scores: pd.Series) -> dict[str, float | bool]:
    scores = scores.reindex(SAMPLE_ORDER)
    full_effect = float(scores[PH_SAMPLES].mean() - scores[BPD_SAMPLES].mean())
    pairwise = np.asarray(
        [
            scores[ph_sample] - scores[bpd_sample]
            for ph_sample, bpd_sample in product(PH_SAMPLES, BPD_SAMPLES)
        ],
        dtype=float,
    )
    leave_one_out = []
    for removed_sample in SAMPLE_ORDER:
        remaining_bpd = [
            sample for sample in BPD_SAMPLES if sample != removed_sample
        ]
        remaining_ph = [
            sample for sample in PH_SAMPLES if sample != removed_sample
        ]
        leave_one_out.append(
            float(scores[remaining_ph].mean() - scores[remaining_bpd].mean())
        )
    leave_one_out_values = np.asarray(leave_one_out, dtype=float)
    full_sign = np.sign(full_effect)
    return {
        "BPD_PH_minus_BPD_mean_score": full_effect,
        "pairwise_min": float(np.min(pairwise)),
        "pairwise_median": float(np.median(pairwise)),
        "pairwise_max": float(np.max(pairwise)),
        "pairwise_positive_fraction": float(np.mean(pairwise > 0)),
        "all_four_pairwise_same_direction_as_full": bool(
            full_sign != 0 and np.all(np.sign(pairwise) == full_sign)
        ),
        "leave_one_donor_out_min": float(np.min(leave_one_out_values)),
        "leave_one_donor_out_median": float(np.median(leave_one_out_values)),
        "leave_one_donor_out_max": float(np.max(leave_one_out_values)),
        "all_leave_one_out_same_direction_as_full": bool(
            full_sign != 0
            and np.all(np.sign(leave_one_out_values) == full_sign)
        ),
    }


def eligible_genes_for_subtype(
    subtype_expression: pd.DataFrame,
) -> pd.Index:
    gene_detection = subtype_expression.groupby("gene").agg(
        total_raw_count=("raw_count", "sum"),
        donors_detected=("detected", "sum"),
    )
    return gene_detection.index[
        gene_detection["total_raw_count"].ge(10)
        & gene_detection["donors_detected"].ge(2)
    ]


def spearman(x: pd.Series, y: pd.Series) -> float:
    valid = x.notna() & y.notna()
    if valid.sum() < 2:
        return float("nan")
    return float(x[valid].rank().corr(y[valid].rank()))


def main() -> int:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    expression = pd.read_csv(EXPRESSION_PATH, sep="\t", low_memory=False)
    pseudobulk_metadata = pd.read_csv(
        PSEUDOBULK_METADATA_PATH,
        sep="\t",
    )
    locked = pd.read_csv(LOCKED_EVIDENCE_PATH, sep="\t")
    locked = locked.dropna(subset=["human_gene"]).drop_duplicates("human_gene")
    locked = locked.set_index("human_gene", drop=False)
    locked["mouse_direction"] = np.sign(locked["mouse_bulk_log2FC"]).astype(int)
    membership_rows = []
    for signature_name, minimum_mouse_subtypes in SIGNATURE_THRESHOLDS.items():
        members = locked.loc[
            locked["concordant_endothelial_subtypes"].ge(
                minimum_mouse_subtypes
            )
        ]
        for gene, row in members.iterrows():
            membership_rows.append(
                {
                    "signature": signature_name,
                    "minimum_concordant_mouse_endothelial_subtypes": (
                        minimum_mouse_subtypes
                    ),
                    "human_gene": gene,
                    "mouse_gene": row["gene"],
                    "mouse_direction": int(row["mouse_direction"]),
                    "weight_before_eligibility_filter": 1.0,
                    "scoring_rule": (
                        "equal-weight mean of donor-standardized log2CPM "
                        "after multiplication by mouse direction"
                    ),
                }
            )
    pd.DataFrame(membership_rows).to_csv(
        RESULT_DIR / "GSE275938_signature_membership_and_weights.tsv",
        sep="\t",
        index=False,
    )

    disease_metadata = pseudobulk_metadata.loc[
        pseudobulk_metadata["sample"].isin(SAMPLE_ORDER)
    ].copy()
    cell_coverage = (
        disease_metadata.pivot(
            index="subtype",
            columns="sample",
            values="cells",
        )
        .reindex(columns=SAMPLE_ORDER)
        .reset_index()
    )
    cell_coverage["minimum_cells_across_donors"] = cell_coverage[
        SAMPLE_ORDER
    ].min(axis=1)
    cell_coverage["maximum_cells_across_donors"] = cell_coverage[
        SAMPLE_ORDER
    ].max(axis=1)
    cell_coverage["max_to_min_cell_ratio"] = (
        cell_coverage["maximum_cells_across_donors"]
        / cell_coverage["minimum_cells_across_donors"]
    )
    for threshold in CELL_THRESHOLDS:
        cell_coverage[f"all_donors_ge_{threshold}_cells"] = (
            cell_coverage["minimum_cells_across_donors"] >= threshold
        )
    cell_coverage.to_csv(
        RESULT_DIR / "GSE275938_disease_endothelial_subtype_cell_coverage.tsv",
        sep="\t",
        index=False,
    )

    qc_by_condition = (
        disease_metadata.groupby(["condition", "subtype"], as_index=False)
        .agg(
            donors=("sample", "nunique"),
            total_cells=("cells", "sum"),
            median_cells=("cells", "median"),
            median_library_size=("library_size", "median"),
            median_genes_detected=("genes_detected", "median"),
            median_cell_nCount_RNA=("median_metadata_nCount_RNA", "median"),
            median_cell_nFeature_RNA=(
                "median_metadata_nFeature_RNA",
                "median",
            ),
            median_cell_percent_mt=("median_metadata_percent_mt", "median"),
        )
    )
    qc_wide = qc_by_condition.pivot(
        index="subtype",
        columns="condition",
        values=[
            "median_cells",
            "median_library_size",
            "median_genes_detected",
            "median_cell_nCount_RNA",
            "median_cell_nFeature_RNA",
            "median_cell_percent_mt",
        ],
    )
    qc_wide.columns = [
        f"{metric}_{condition.replace('+', '_plus_')}"
        for metric, condition in qc_wide.columns
    ]
    qc_wide = qc_wide.reset_index()
    qc_wide["BPD_to_BPD_PH_median_cell_nCount_ratio"] = (
        qc_wide["median_cell_nCount_RNA_BPD"]
        / qc_wide["median_cell_nCount_RNA_BPD_plus_PH"]
    )
    qc_wide["BPD_to_BPD_PH_median_cell_nFeature_ratio"] = (
        qc_wide["median_cell_nFeature_RNA_BPD"]
        / qc_wide["median_cell_nFeature_RNA_BPD_plus_PH"]
    )
    qc_wide.to_csv(
        RESULT_DIR / "GSE275938_disease_endothelial_QC_by_condition.tsv",
        sep="\t",
        index=False,
    )

    donor_score_rows: list[dict] = []
    signature_summary_rows: list[dict] = []
    gene_effect_rows: list[dict] = []
    gene_direction_summary_rows: list[dict] = []

    subtype_order = list(dict.fromkeys(expression["subtype"]))
    minimum_cells_map = disease_metadata.groupby("subtype")["cells"].min()

    for subtype in subtype_order:
        subtype_expression = expression.loc[
            expression["subtype"].eq(subtype)
        ].copy()
        eligible_genes = eligible_genes_for_subtype(subtype_expression)

        for normalization_name, expression_column in NORMALIZATIONS.items():
            expression_matrix = subtype_expression.pivot(
                index="sample",
                columns="gene",
                values=expression_column,
            ).reindex(SAMPLE_ORDER)
            raw_matrix = subtype_expression.pivot(
                index="sample",
                columns="gene",
                values="raw_count",
            ).reindex(SAMPLE_ORDER)

            gene_effects = (
                expression_matrix.loc[PH_SAMPLES].mean(axis=0)
                - expression_matrix.loc[BPD_SAMPLES].mean(axis=0)
            )
            for gene in expression_matrix.columns:
                mouse_direction = int(locked.loc[gene, "mouse_direction"])
                pairwise_gene_effects = np.asarray(
                    [
                        expression_matrix.loc[ph_sample, gene]
                        - expression_matrix.loc[bpd_sample, gene]
                        for ph_sample, bpd_sample in product(
                            PH_SAMPLES,
                            BPD_SAMPLES,
                        )
                    ]
                )
                gene_effect_rows.append(
                    {
                        "subtype": subtype,
                        "normalization": normalization_name,
                        "gene": gene,
                        "mouse_bulk_log2FC": locked.loc[
                            gene, "mouse_bulk_log2FC"
                        ],
                        "mouse_concordant_endothelial_subtypes": locked.loc[
                            gene, "concordant_endothelial_subtypes"
                        ],
                        "eligible_expression": gene in eligible_genes,
                        "total_raw_count": int(raw_matrix[gene].sum()),
                        "donors_detected": int((raw_matrix[gene] > 0).sum()),
                        "human_BPD_PH_minus_BPD_log2CPM": float(
                            gene_effects[gene]
                        ),
                        "human_mouse_direction_concordant": bool(
                            np.sign(gene_effects[gene]) == mouse_direction
                        ),
                        "pairwise_mouse_direction_concordant_count": int(
                            np.sum(
                                np.sign(pairwise_gene_effects)
                                == mouse_direction
                            )
                        ),
                        "all_four_pairwise_mouse_direction_concordant": bool(
                            np.all(
                                np.sign(pairwise_gene_effects)
                                == mouse_direction
                            )
                        ),
                    }
                )

            for signature_name, minimum_mouse_subtypes in SIGNATURE_THRESHOLDS.items():
                signature_genes = locked.index[
                    locked["concordant_endothelial_subtypes"].ge(
                        minimum_mouse_subtypes
                    )
                ]
                scored_genes = pd.Index(signature_genes).intersection(
                    eligible_genes
                )
                scored_genes = scored_genes.intersection(expression_matrix.columns)
                directions = locked.loc[scored_genes, "mouse_direction"]
                scores, variable_gene_count = oriented_signature_score(
                    expression_matrix[scored_genes],
                    directions,
                )
                for sample in SAMPLE_ORDER:
                    sample_meta = disease_metadata.loc[
                        disease_metadata["sample"].eq(sample)
                        & disease_metadata["subtype"].eq(subtype)
                    ].iloc[0]
                    donor_score_rows.append(
                        {
                            "sample": sample,
                            "condition": sample_meta["condition"],
                            "subtype": subtype,
                            "normalization": normalization_name,
                            "signature": signature_name,
                            "cells": int(sample_meta["cells"]),
                            "genes_requested": int(len(signature_genes)),
                            "genes_eligible": int(len(scored_genes)),
                            "genes_variable": variable_gene_count,
                            "oriented_signature_score": float(scores[sample]),
                        }
                    )
                summary = donor_contrast_summary(scores)
                signature_summary_rows.append(
                    {
                        "subtype": subtype,
                        "normalization": normalization_name,
                        "signature": signature_name,
                        "minimum_cells_across_donors": int(
                            minimum_cells_map[subtype]
                        ),
                        "genes_requested": int(len(signature_genes)),
                        "genes_eligible": int(len(scored_genes)),
                        "genes_variable": variable_gene_count,
                        **summary,
                    }
                )

            eligible_effects = pd.DataFrame(
                {
                    "human_effect": gene_effects,
                    "mouse_effect": locked.loc[
                        gene_effects.index,
                        "mouse_bulk_log2FC",
                    ],
                    "mouse_breadth": locked.loc[
                        gene_effects.index,
                        "concordant_endothelial_subtypes",
                    ],
                    "eligible": gene_effects.index.isin(eligible_genes),
                }
            )
            for minimum_mouse_subtypes in [1, 2, 3]:
                subset = eligible_effects.loc[
                    eligible_effects["eligible"]
                    & eligible_effects["mouse_breadth"].ge(
                        minimum_mouse_subtypes
                    )
                ]
                direction_concordant = (
                    np.sign(subset["human_effect"])
                    == np.sign(subset["mouse_effect"])
                )
                gene_direction_summary_rows.append(
                    {
                        "subtype": subtype,
                        "normalization": normalization_name,
                        "minimum_mouse_endothelial_subtypes": (
                            minimum_mouse_subtypes
                        ),
                        "minimum_cells_across_donors": int(
                            minimum_cells_map[subtype]
                        ),
                        "genes_evaluated": int(len(subset)),
                        "direction_concordant_genes": int(
                            direction_concordant.sum()
                        ),
                        "direction_concordance_fraction": float(
                            direction_concordant.mean()
                        ),
                        "spearman_mouse_vs_human_effect": spearman(
                            subset["mouse_effect"],
                            subset["human_effect"],
                        ),
                    }
                )

    donor_scores = pd.DataFrame(donor_score_rows)
    signature_summary = pd.DataFrame(signature_summary_rows)
    gene_effects = pd.DataFrame(gene_effect_rows)
    gene_direction_summary = pd.DataFrame(gene_direction_summary_rows)

    donor_scores.to_csv(
        RESULT_DIR / "GSE275938_subtype_locked_signature_donor_scores.tsv",
        sep="\t",
        index=False,
    )
    signature_summary.to_csv(
        RESULT_DIR / "GSE275938_subtype_locked_signature_sensitivity.tsv",
        sep="\t",
        index=False,
    )
    gene_effects.to_csv(
        RESULT_DIR / "GSE275938_subtype_locked_gene_effects.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    gene_direction_summary.to_csv(
        RESULT_DIR / "GSE275938_subtype_gene_direction_summary.tsv",
        sep="\t",
        index=False,
    )

    normalization_comparison = signature_summary.pivot(
        index=["subtype", "signature", "minimum_cells_across_donors"],
        columns="normalization",
        values="BPD_PH_minus_BPD_mean_score",
    ).reset_index()
    normalization_comparison["same_effect_direction"] = (
        np.sign(normalization_comparison["TMM"])
        == np.sign(normalization_comparison["library_CPM"])
    )
    normalization_comparison["TMM_minus_library_CPM_effect"] = (
        normalization_comparison["TMM"]
        - normalization_comparison["library_CPM"]
    )
    normalization_comparison.to_csv(
        RESULT_DIR / "GSE275938_normalization_sensitivity.tsv",
        sep="\t",
        index=False,
    )

    threshold_rows = []
    tmm_summary = signature_summary.loc[
        signature_summary["normalization"].eq("TMM")
    ]
    for signature_name in SIGNATURE_THRESHOLDS:
        signature_subset = tmm_summary.loc[
            tmm_summary["signature"].eq(signature_name)
        ]
        for threshold in CELL_THRESHOLDS:
            passing = signature_subset.loc[
                signature_subset["minimum_cells_across_donors"].ge(threshold)
            ]
            biological = passing.loc[
                passing["subtype"].ne(ALL_ENDOTHELIAL)
            ]
            threshold_rows.append(
                {
                    "signature": signature_name,
                    "minimum_cells_per_donor_threshold": threshold,
                    "passing_biological_subtypes": int(len(biological)),
                    "biological_subtype_names": ";".join(
                        biological["subtype"].tolist()
                    ),
                    "positive_effect_biological_subtypes": int(
                        biological["BPD_PH_minus_BPD_mean_score"].gt(0).sum()
                    ),
                    "all_pairwise_robust_positive_biological_subtypes": int(
                        (
                            biological[
                                "all_four_pairwise_same_direction_as_full"
                            ]
                            & biological[
                                "BPD_PH_minus_BPD_mean_score"
                            ].gt(0)
                        ).sum()
                    ),
                    "all_endothelial_passes": bool(
                        ALL_ENDOTHELIAL in set(passing["subtype"])
                    ),
                    "all_endothelial_effect": (
                        float(
                            passing.loc[
                                passing["subtype"].eq(ALL_ENDOTHELIAL),
                                "BPD_PH_minus_BPD_mean_score",
                            ].iloc[0]
                        )
                        if ALL_ENDOTHELIAL in set(passing["subtype"])
                        else np.nan
                    ),
                }
            )
    threshold_summary = pd.DataFrame(threshold_rows)
    threshold_summary.to_csv(
        RESULT_DIR / "GSE275938_minimum_cell_threshold_sensitivity.tsv",
        sep="\t",
        index=False,
    )

    # Equal-cell downsampling uses cells only to assess sampling sensitivity.
    cell_metadata = pd.read_csv(CELL_METADATA_PATH, sep="\t")
    cell_payload = np.load(CELL_COUNTS_PATH)
    cell_counts = cell_payload["counts"]
    cell_genes = [str(value) for value in cell_payload["genes"]]
    if len(cell_metadata) != cell_counts.shape[0]:
        raise RuntimeError("Cell metadata/count matrix row mismatch.")
    if cell_genes != list(expression["gene"].drop_duplicates()):
        # Order may differ from the long table but the set must be identical.
        if set(cell_genes) != set(expression["gene"]):
            raise RuntimeError("Cell and pseudobulk locked-gene sets differ.")
    cell_gene_to_position = {
        gene: position for position, gene in enumerate(cell_genes)
    }
    random_generator = np.random.default_rng(RNG_SEED)
    downsampling_rows = []

    downsampling_subtypes = [
        subtype
        for subtype in subtype_order
        if int(minimum_cells_map[subtype]) >= 10
    ]
    for subtype in downsampling_subtypes:
        subtype_mask = (
            np.ones(len(cell_metadata), dtype=bool)
            if subtype == ALL_ENDOTHELIAL
            else cell_metadata["subtype"].eq(subtype).to_numpy()
        )
        sample_indices = {
            sample: np.flatnonzero(
                subtype_mask
                & cell_metadata["sample"].eq(sample).to_numpy()
            )
            for sample in SAMPLE_ORDER
        }
        downsample_cells = min(len(indices) for indices in sample_indices.values())
        subtype_expression = expression.loc[
            expression["subtype"].eq(subtype)
        ]
        eligible_genes = eligible_genes_for_subtype(subtype_expression)

        for signature_name, minimum_mouse_subtypes in SIGNATURE_THRESHOLDS.items():
            signature_genes = locked.index[
                locked["concordant_endothelial_subtypes"].ge(
                    minimum_mouse_subtypes
                )
            ]
            scored_genes = (
                pd.Index(signature_genes)
                .intersection(eligible_genes)
                .intersection(cell_genes)
            )
            gene_positions = np.asarray(
                [cell_gene_to_position[gene] for gene in scored_genes],
                dtype=int,
            )
            directions = locked.loc[scored_genes, "mouse_direction"]

            for repeat in range(DOWNSAMPLING_REPEATS):
                donor_expression_rows = []
                for sample in SAMPLE_ORDER:
                    chosen = random_generator.choice(
                        sample_indices[sample],
                        size=downsample_cells,
                        replace=False,
                    )
                    aggregate = cell_counts[
                        chosen[:, None],
                        gene_positions[None, :],
                    ].sum(axis=0)
                    total_umi = float(
                        cell_metadata.loc[chosen, "cell_total_umi"].sum()
                    )
                    donor_expression_rows.append(
                        np.log2(aggregate / total_umi * 1e6 + 0.5)
                    )
                downsampled_expression = pd.DataFrame(
                    donor_expression_rows,
                    index=SAMPLE_ORDER,
                    columns=scored_genes,
                )
                scores, variable_gene_count = oriented_signature_score(
                    downsampled_expression,
                    directions,
                )
                effect = float(
                    scores[PH_SAMPLES].mean() - scores[BPD_SAMPLES].mean()
                )
                downsampling_rows.append(
                    {
                        "subtype": subtype,
                        "signature": signature_name,
                        "repeat": repeat + 1,
                        "cells_per_donor": downsample_cells,
                        "genes_eligible": len(scored_genes),
                        "genes_variable": variable_gene_count,
                        "BPD_PH_minus_BPD_mean_score": effect,
                    }
                )

    downsampling = pd.DataFrame(downsampling_rows)
    downsampling.to_csv(
        RESULT_DIR / "GSE275938_equal_cell_downsampling_all_repeats.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    downsampling_summary = (
        downsampling.groupby(
            ["subtype", "signature", "cells_per_donor"],
            as_index=False,
        )
        .agg(
            repeats=("repeat", "size"),
            median_effect=("BPD_PH_minus_BPD_mean_score", "median"),
            sampling_p2_5=(
                "BPD_PH_minus_BPD_mean_score",
                lambda values: values.quantile(0.025),
            ),
            sampling_p97_5=(
                "BPD_PH_minus_BPD_mean_score",
                lambda values: values.quantile(0.975),
            ),
            positive_effect_fraction=(
                "BPD_PH_minus_BPD_mean_score",
                lambda values: float((values > 0).mean()),
            ),
        )
    )
    downsampling_summary.to_csv(
        RESULT_DIR / "GSE275938_equal_cell_downsampling_summary.tsv",
        sep="\t",
        index=False,
    )

    # Molecule-depth sensitivity: binomially thin each donor pseudobulk to
    # the minimum total-UMI depth within subtype, then recompute the score.
    # This addresses library-depth imbalance without treating cells as
    # independent replicates.
    rarefaction_rows = []
    for subtype in downsampling_subtypes:
        subtype_expression = expression.loc[
            expression["subtype"].eq(subtype)
        ].copy()
        eligible_genes = eligible_genes_for_subtype(subtype_expression)
        raw_matrix = subtype_expression.pivot(
            index="sample",
            columns="gene",
            values="raw_count",
        ).reindex(SAMPLE_ORDER)
        subtype_meta = (
            disease_metadata.loc[disease_metadata["subtype"].eq(subtype)]
            .set_index("sample")
            .reindex(SAMPLE_ORDER)
        )
        library_sizes = subtype_meta["library_size"].astype(float)
        target_depth = int(library_sizes.min())

        for signature_name, minimum_mouse_subtypes in SIGNATURE_THRESHOLDS.items():
            signature_genes = locked.index[
                locked["concordant_endothelial_subtypes"].ge(
                    minimum_mouse_subtypes
                )
            ]
            scored_genes = (
                pd.Index(signature_genes)
                .intersection(eligible_genes)
                .intersection(raw_matrix.columns)
            )
            directions = locked.loc[scored_genes, "mouse_direction"]
            for repeat in range(RAREFACTION_REPEATS):
                donor_rows = []
                for sample in SAMPLE_ORDER:
                    probability = min(
                        1.0,
                        target_depth / float(library_sizes[sample]),
                    )
                    thinned = random_generator.binomial(
                        raw_matrix.loc[sample, scored_genes]
                        .to_numpy(dtype=np.int64),
                        probability,
                    )
                    donor_rows.append(
                        np.log2(thinned / target_depth * 1e6 + 0.5)
                    )
                thinned_expression = pd.DataFrame(
                    donor_rows,
                    index=SAMPLE_ORDER,
                    columns=scored_genes,
                )
                scores, variable_gene_count = oriented_signature_score(
                    thinned_expression,
                    directions,
                )
                rarefaction_rows.append(
                    {
                        "subtype": subtype,
                        "signature": signature_name,
                        "repeat": repeat + 1,
                        "target_total_umi_per_donor": target_depth,
                        "genes_eligible": len(scored_genes),
                        "genes_variable": variable_gene_count,
                        "BPD_PH_minus_BPD_mean_score": float(
                            scores[PH_SAMPLES].mean()
                            - scores[BPD_SAMPLES].mean()
                        ),
                    }
                )
    rarefaction = pd.DataFrame(rarefaction_rows)
    rarefaction.to_csv(
        RESULT_DIR / "GSE275938_library_depth_rarefaction_all_repeats.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    rarefaction_summary = (
        rarefaction.groupby(
            ["subtype", "signature", "target_total_umi_per_donor"],
            as_index=False,
        )
        .agg(
            repeats=("repeat", "size"),
            median_effect=("BPD_PH_minus_BPD_mean_score", "median"),
            sampling_p2_5=(
                "BPD_PH_minus_BPD_mean_score",
                lambda values: values.quantile(0.025),
            ),
            sampling_p97_5=(
                "BPD_PH_minus_BPD_mean_score",
                lambda values: values.quantile(0.975),
            ),
            positive_effect_fraction=(
                "BPD_PH_minus_BPD_mean_score",
                lambda values: float((values > 0).mean()),
            ),
        )
    )
    rarefaction_summary.to_csv(
        RESULT_DIR / "GSE275938_library_depth_rarefaction_summary.tsv",
        sep="\t",
        index=False,
    )

    main_signature = signature_summary.loc[
        signature_summary["normalization"].eq("TMM")
        & signature_summary["signature"].eq("mouse_subtypes_ge3")
    ].copy()
    main_signature = main_signature.sort_values(
        ["minimum_cells_across_donors", "subtype"],
        ascending=[False, True],
    )
    reliable_subtypes = main_signature.loc[
        main_signature["minimum_cells_across_donors"].ge(10)
    ]
    summary = {
        "analysis_version": "1.0",
        "biological_unit": "human donor",
        "disease_donors": {"BPD": 2, "BPD+PH": 2},
        "locked_genes_available": int(expression["gene"].nunique()),
        "signature_sets": SIGNATURE_THRESHOLDS,
        "minimum_cell_threshold_primary": 10,
        "subtypes_with_at_least_10_cells_per_donor": reliable_subtypes[
            "subtype"
        ].tolist(),
        "ge3_signature_positive_effect_subtypes_at_threshold10": reliable_subtypes.loc[
            reliable_subtypes["BPD_PH_minus_BPD_mean_score"].gt(0),
            "subtype",
        ].tolist(),
        "ge3_signature_all_pairwise_robust_subtypes_at_threshold10": (
            reliable_subtypes.loc[
                reliable_subtypes[
                    "all_four_pairwise_same_direction_as_full"
                ],
                "subtype",
            ].tolist()
        ),
        "normalization_same_direction_fraction": float(
            normalization_comparison["same_effect_direction"].mean()
        ),
        "downsampling_repeats": DOWNSAMPLING_REPEATS,
        "library_depth_rarefaction_repeats": RAREFACTION_REPEATS,
        "critical_limitations": [
            "BPD and BPD+PH each contain only two donors.",
            "Endothelial-cell recovery is severely imbalanced across donors.",
            "BPD endothelial cells have substantially greater per-cell UMI and feature counts than BPD+PH cells, so diagnosis is confounded with library/donor quality.",
            "Subtype results with fewer than 10 cells in any donor are excluded from the primary subtype interpretation.",
            "No cell-level p-values or definitive human validation claims are produced.",
        ],
    }
    (RESULT_DIR / "GSE275938_human_subtype_sensitivity_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("\nPrimary high-confidence signature:")
    print(main_signature.to_string(index=False))
    print("\nEqual-cell downsampling:")
    print(
        downsampling_summary.loc[
            downsampling_summary["signature"].eq("mouse_subtypes_ge3")
        ].to_string(index=False)
    )
    print("\nLibrary-depth rarefaction:")
    print(
        rarefaction_summary.loc[
            rarefaction_summary["signature"].eq("mouse_subtypes_ge3")
        ].to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
