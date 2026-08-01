"""Expression/detection-matched random-gene specificity check for human scores.

This is a technical specificity diagnostic, not a replacement for donor-level
replication. Random sets preserve signature size, mouse-direction composition,
mean-expression decile and donor detection count.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from build_ncbi_mouse_human_ortholog_mapping import (
    read_gene_info,
    read_mouse_human_pairs,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROCESSED_DIR = PROJECT_ROOT / "04_processed" / "GSE275938"
RESULT_DIR = PROJECT_ROOT / "07_results" / "human_endothelial_subtypes"
EXTERNAL_DIR = PROJECT_ROOT / "external_data"
COUNTS_PATH = (
    PROCESSED_DIR / "GSE275938_endothelial_subtype_pseudobulk_counts.tsv.gz"
)
METADATA_PATH = (
    PROCESSED_DIR / "GSE275938_endothelial_subtype_pseudobulk_metadata.tsv"
)
NORMALIZATION_QC_PATH = (
    RESULT_DIR / "GSE275938_subtype_TMM_normalization_QC.tsv"
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
SIGNATURE_THRESHOLDS = {
    "mouse_subtypes_ge1": 1,
    "mouse_subtypes_ge2": 2,
    "mouse_subtypes_ge3": 3,
}
PERMUTATIONS = 2000
RNG_SEED = 20260726


def edger_log2_cpm(
    counts: pd.DataFrame,
    effective_library_sizes: pd.Series,
    prior_count: float = 0.5,
) -> pd.DataFrame:
    """Reproduce edgeR::cpm(log=TRUE, prior.count=0.5).

    edgeR scales the prior count by each sample's effective library size and
    adds the scaled prior to both ends of the library.  Keeping this transform
    identical to the primary human analysis prevents the matched-null audit
    from using a subtly different observed effect.
    """
    scaled_prior = (
        prior_count
        * effective_library_sizes
        / effective_library_sizes.mean()
    )
    adjusted_library_sizes = effective_library_sizes + 2 * scaled_prior
    return np.log2(
        counts.add(scaled_prior, axis=0)
        .div(adjusted_library_sizes, axis=0)
        * 1e6
    )


def one_to_one_human_gene_universe() -> set[str]:
    """Return current human symbols with reciprocal one-to-one NCBI orthology."""
    human_info = read_gene_info(
        EXTERNAL_DIR / "Homo_sapiens.gene_info_20260730.gz"
    )
    pairs = read_mouse_human_pairs(
        EXTERNAL_DIR / "gene_orthologs_20260730.gz"
    )
    one_to_one = pairs.loc[
        pairs["mouse_to_human_n"].eq(1)
        & pairs["human_to_mouse_n"].eq(1),
        ["human_GeneID"],
    ].drop_duplicates()
    one_to_one = one_to_one.merge(
        human_info[["GeneID", "Symbol"]],
        left_on="human_GeneID",
        right_on="GeneID",
        how="inner",
    )
    return set(one_to_one["Symbol"].dropna().astype(str))


def bh_adjust(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = numeric.dropna().sort_values()
    if valid.empty:
        return result
    ranks = np.arange(1, len(valid) + 1)
    adjusted = valid.to_numpy() * len(valid) / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[valid.index] = np.minimum(adjusted, 1.0)
    return result


def signature_effect(
    expression: pd.DataFrame,
    directions: pd.Series,
) -> float:
    expression = expression.loc[:, directions.index]
    standard_deviation = expression.std(axis=0, ddof=0)
    variable = standard_deviation.index[standard_deviation > 0]
    if not len(variable):
        return float("nan")
    z_scores = expression[variable].sub(
        expression[variable].mean(axis=0),
        axis=1,
    ).div(standard_deviation[variable], axis=1)
    donor_scores = z_scores.mul(directions[variable], axis=1).mean(axis=1)
    return float(
        donor_scores[PH_SAMPLES].mean()
        - donor_scores[BPD_SAMPLES].mean()
    )


def main() -> int:
    counts = pd.read_csv(COUNTS_PATH, sep="\t", low_memory=False)
    counts = counts.set_index("pseudobulk_id")
    metadata = pd.read_csv(METADATA_PATH, sep="\t")
    metadata = metadata.loc[metadata["sample"].isin(SAMPLE_ORDER)].copy()
    normalization = pd.read_csv(NORMALIZATION_QC_PATH, sep="\t")
    locked = pd.read_csv(LOCKED_EVIDENCE_PATH, sep="\t")
    locked = (
        locked.dropna(subset=["human_gene"])
        .drop_duplicates("human_gene")
        .set_index("human_gene", drop=False)
    )
    locked["mouse_direction"] = np.sign(locked["mouse_bulk_log2FC"]).astype(int)
    all_locked_genes = set(locked.index)
    one_to_one_human_genes = one_to_one_human_gene_universe()
    random_generator = np.random.default_rng(RNG_SEED)

    summary_rows = []
    null_rows = []
    subtype_minimum_cells = metadata.groupby("subtype")["cells"].min()
    tested_subtypes = subtype_minimum_cells.index[
        subtype_minimum_cells.ge(10)
    ]

    for subtype in tested_subtypes:
        subtype_metadata = (
            metadata.loc[metadata["subtype"].eq(subtype)]
            .set_index("sample")
            .loc[SAMPLE_ORDER]
        )
        pseudobulk_ids = subtype_metadata["pseudobulk_id"].tolist()
        subtype_counts = counts.loc[pseudobulk_ids].copy()
        subtype_counts.index = SAMPLE_ORDER

        subtype_normalization = (
            normalization.loc[normalization["subtype"].eq(subtype)]
            .set_index("sample")
            .loc[SAMPLE_ORDER]
        )
        effective_library_sizes = subtype_normalization[
            "effective_library_size"
        ]
        expression = edger_log2_cpm(
            subtype_counts,
            effective_library_sizes,
            prior_count=0.5,
        )

        total_count = subtype_counts.sum(axis=0)
        donors_detected = (subtype_counts > 0).sum(axis=0)
        expression_sd = expression.std(axis=0, ddof=0)
        eligible_universe = total_count.index[
            total_count.ge(10)
            & donors_detected.ge(2)
            & expression_sd.gt(0)
            & total_count.index.isin(one_to_one_human_genes)
        ]
        standardized_gene_effect = (
            expression.loc[PH_SAMPLES, eligible_universe].mean(axis=0)
            - expression.loc[BPD_SAMPLES, eligible_universe].mean(axis=0)
        ).div(expression_sd.loc[eligible_universe])
        eligible_gene_list = eligible_universe.tolist()
        gene_to_position = {
            gene: position for position, gene in enumerate(eligible_gene_list)
        }
        standardized_effect_array = standardized_gene_effect.loc[
            eligible_gene_list
        ].to_numpy()
        mean_expression = expression[eligible_universe].mean(axis=0)
        expression_decile = pd.qcut(
            mean_expression.rank(method="first"),
            q=10,
            labels=False,
        ).astype(int)
        strata = pd.DataFrame(
            {
                "gene": eligible_universe,
                "expression_decile": expression_decile.loc[
                    eligible_universe
                ].to_numpy(),
                "donors_detected": donors_detected.loc[
                    eligible_universe
                ].to_numpy(),
            }
        ).set_index("gene")
        random_universe = strata.loc[
            ~strata.index.isin(all_locked_genes)
        ].copy()
        random_universe["position"] = [
            gene_to_position[gene] for gene in random_universe.index
        ]
        pool_by_stratum = {
            key: group["position"].to_numpy(dtype=int)
            for key, group in random_universe.groupby(
                ["expression_decile", "donors_detected"],
                sort=False,
            )
        }
        pool_by_decile = {
            key: group["position"].to_numpy(dtype=int)
            for key, group in random_universe.groupby(
                "expression_decile",
                sort=False,
            )
        }

        for signature_name, minimum_mouse_subtypes in SIGNATURE_THRESHOLDS.items():
            observed_genes = locked.index[
                locked["concordant_endothelial_subtypes"].ge(
                    minimum_mouse_subtypes
                )
                & locked.index.isin(eligible_universe)
            ]
            observed_directions = locked.loc[
                observed_genes,
                "mouse_direction",
            ]
            observed_positions = np.asarray(
                [gene_to_position[gene] for gene in observed_genes],
                dtype=int,
            )
            observed_effect = float(
                np.mean(
                    standardized_effect_array[observed_positions]
                    * observed_directions.to_numpy(dtype=float)
                )
            )

            observed_strata = strata.loc[observed_genes].copy()
            observed_strata["mouse_direction"] = observed_directions
            stratum_groups = list(
                observed_strata.groupby(
                    ["expression_decile", "donors_detected"],
                    sort=True,
                )
            )
            prepared_groups = []
            for (expression_bin, detection_count), group in stratum_groups:
                primary_pool = pool_by_stratum.get(
                    (expression_bin, detection_count),
                    np.asarray([], dtype=int),
                )
                fallback_pool = pool_by_decile[expression_bin]
                prepared_groups.append(
                    (
                        primary_pool,
                        fallback_pool,
                        group["mouse_direction"].to_numpy(dtype=float),
                    )
                )
            null_effects = np.empty(PERMUTATIONS, dtype=float)

            for permutation in range(PERMUTATIONS):
                selected_positions: list[int] = []
                selected_directions: list[float] = []
                used_positions: set[int] = set()
                for primary_pool, fallback_pool, directions_template in prepared_groups:
                    candidate_pool = primary_pool
                    if len(candidate_pool) < len(directions_template):
                        candidate_pool = fallback_pool
                    if used_positions:
                        candidate_pool = candidate_pool[
                            ~np.isin(
                                candidate_pool,
                                np.fromiter(used_positions, dtype=int),
                            )
                        ]
                    chosen = random_generator.choice(
                        candidate_pool,
                        size=len(directions_template),
                        replace=False,
                    )
                    directions = directions_template.copy()
                    random_generator.shuffle(directions)
                    selected_positions.extend(chosen.tolist())
                    selected_directions.extend(directions.tolist())
                    used_positions.update(chosen.tolist())
                null_effects[permutation] = float(
                    np.mean(
                        standardized_effect_array[
                            np.asarray(selected_positions, dtype=int)
                        ]
                        * np.asarray(selected_directions, dtype=float)
                    )
                )
                null_rows.append(
                    {
                        "subtype": subtype,
                        "signature": signature_name,
                        "permutation": permutation + 1,
                        "random_signature_effect": null_effects[permutation],
                    }
                )

            empirical_two_sided = (
                1 + np.sum(np.abs(null_effects) >= abs(observed_effect))
            ) / (PERMUTATIONS + 1)
            empirical_greater = (
                1 + np.sum(null_effects >= observed_effect)
            ) / (PERMUTATIONS + 1)
            summary_rows.append(
                {
                    "subtype": subtype,
                    "signature": signature_name,
                    "minimum_cells_across_donors": int(
                        subtype_minimum_cells[subtype]
                    ),
                    "observed_genes": len(observed_genes),
                    "eligible_random_gene_universe": int(
                        len(random_universe)
                    ),
                    "mouse_up_genes": int((observed_directions > 0).sum()),
                    "mouse_down_genes": int((observed_directions < 0).sum()),
                    "observed_effect": observed_effect,
                    "null_median": float(np.median(null_effects)),
                    "null_p2_5": float(np.quantile(null_effects, 0.025)),
                    "null_p97_5": float(np.quantile(null_effects, 0.975)),
                    "observed_percentile_in_null": float(
                        np.mean(null_effects <= observed_effect)
                    ),
                    "empirical_p_greater": empirical_greater,
                    "empirical_p_two_sided_abs_effect": empirical_two_sided,
                    "permutations": PERMUTATIONS,
                }
            )

    summary = pd.DataFrame(summary_rows)
    summary["empirical_q_greater_BH_within_signature"] = (
        summary.groupby("signature", group_keys=False)[
            "empirical_p_greater"
        ].apply(bh_adjust)
    )
    summary["empirical_q_two_sided_BH_within_signature"] = (
        summary.groupby("signature", group_keys=False)[
            "empirical_p_two_sided_abs_effect"
        ].apply(bh_adjust)
    )
    null_table = pd.DataFrame(null_rows)
    summary.to_csv(
        RESULT_DIR / "GSE275938_expression_matched_random_set_summary.tsv",
        sep="\t",
        index=False,
    )
    null_table.to_csv(
        RESULT_DIR
        / "GSE275938_expression_matched_random_set_all_permutations.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    result_summary = {
        "permutations_per_subtype_signature": PERMUTATIONS,
        "tested_subtypes": tested_subtypes.tolist(),
        "matching_variables": [
            "mean expression decile",
            "number of disease donors with nonzero counts",
            "signature size",
            "mouse up/down direction composition within matching stratum",
        ],
        "eligible_random_gene_universe": (
            "Genes with total raw count >=10, nonzero counts in >=2 disease "
            "donors, nonzero donor-level expression variance, and a reciprocal "
            "one-to-one mouse-human relationship in NCBI Gene; all formal "
            "mouse-signature genes were excluded from the random pool."
        ),
        "orthology_restriction": (
            "Reciprocal one-to-one mouse-human mappings in NCBI Gene "
            "gene_orthologs downloaded 2026-07-30; human symbols were taken "
            "from Homo_sapiens.gene_info downloaded on the same date."
        ),
        "expression_transform": (
            "TMM-normalized log2 CPM reproduced exactly from "
            "edgeR::cpm(log=TRUE, prior.count=0.5, "
            "normalized.lib.sizes=TRUE), including effective-library-size "
            "scaling of the prior count."
        ),
        "mouse_direction_assignment_for_random_sets": (
            "Random human genes do not receive gene-specific mouse effects. "
            "Within each expression/detection stratum, the observed signature "
            "mouse-direction template is shuffled and assigned to sampled "
            "genes, preserving signature size and up/down composition."
        ),
        "empirical_p_formula": "(1 + number of null effects >= observed effect) / (2000 + 1)",
        "pool_construction": (
            "A separate eligible and matched random-gene pool was constructed "
            "for each tested endothelial population."
        ),
        "sampling_with_replacement": False,
        "multiple_testing": (
            "Benjamini-Hochberg adjustment across all tested endothelial "
            "subtypes separately for each signature definition."
        ),
        "interpretation": (
            "The empirical test evaluates signature specificity relative to "
            "expression-matched random genes; it does not create additional "
            "human donors or establish clinical validation."
        ),
    }
    (
        RESULT_DIR / "GSE275938_expression_matched_random_set_method.json"
    ).write_text(
        json.dumps(result_summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
