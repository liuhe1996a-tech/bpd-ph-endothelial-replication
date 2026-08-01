"""Recover and summarize completed fgsea results after an interrupted tail step."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULT_DIR = PROJECT_ROOT / "07_results" / "pathway_exact"


def joined_unique(values: pd.Series) -> str:
    return ";".join(sorted({str(value) for value in values if pd.notna(value)}))


def main() -> None:
    replication_path = RESULT_DIR / "mouse_pathway_replication_all.tsv.gz"
    replication = pd.read_csv(replication_path, sep="\t", low_memory=False)
    for column in ("both_fdr05", "same_direction", "replicated"):
        if replication[column].dtype != bool:
            replication[column] = (
                replication[column].astype(str).str.upper().eq("TRUE")
            )

    primary = replication.loc[
        replication["model"].eq("age_adjusted_all_ages")
        & replication["replicated"]
    ].copy()
    primary.to_csv(
        RESULT_DIR / "mouse_pathway_replication_primary.tsv",
        sep="\t",
        index=False,
    )

    primary_breadth = (
        primary.groupby(["collection", "pathway"], as_index=False)
        .agg(
            direction=("direction", "first"),
            primary_subtypes_replicated=("cell_type", "nunique"),
            primary_subtypes=("cell_type", joined_unique),
            bulk_NES=("bulk_NES", "first"),
            bulk_padj=("bulk_padj", "first"),
            median_endothelial_NES=("endothelial_NES", "median"),
            minimum_leading_edge_overlap=("leading_edge_overlap_n", "min"),
            median_leading_edge_overlap=("leading_edge_overlap_n", "median"),
        )
    )

    p14 = replication.loc[
        replication["model"].eq("P14_oxygen_only")
        & replication["replicated"]
    ].copy()
    p14_breadth = (
        p14.groupby(["collection", "pathway"], as_index=False)
        .agg(
            p14_subtypes_replicated=("cell_type", "nunique"),
            p14_subtypes=("cell_type", joined_unique),
        )
    )
    breadth = primary_breadth.merge(
        p14_breadth,
        on=["collection", "pathway"],
        how="left",
    )
    breadth["p14_subtypes_replicated"] = (
        breadth["p14_subtypes_replicated"].fillna(0).astype(int)
    )
    breadth["p14_subtypes"] = breadth["p14_subtypes"].fillna("")
    breadth = breadth.sort_values(
        [
            "primary_subtypes_replicated",
            "p14_subtypes_replicated",
            "bulk_padj",
            "median_leading_edge_overlap",
        ],
        ascending=[False, False, True, False],
        kind="stable",
    )
    breadth.to_csv(
        RESULT_DIR / "mouse_pathway_replication_breadth.tsv",
        sep="\t",
        index=False,
    )

    summary = (
        replication.groupby(
            ["model", "cell_type", "collection"], as_index=False
        )
        .agg(
            pathways_compared=("pathway", "size"),
            pathways_both_fdr05=("both_fdr05", "sum"),
            pathways_replicated_same_direction=("replicated", "sum"),
        )
    )
    direction_counts = (
        replication.loc[replication["replicated"]]
        .groupby(
            ["model", "cell_type", "collection", "direction"],
            as_index=False,
        )
        .size()
        .pivot(
            index=["model", "cell_type", "collection"],
            columns="direction",
            values="size",
        )
        .fillna(0)
        .reset_index()
    )
    direction_counts.columns.name = None
    direction_counts = direction_counts.rename(
        columns={
            "up": "pathways_replicated_up",
            "down": "pathways_replicated_down",
        }
    )
    summary = summary.merge(
        direction_counts,
        on=["model", "cell_type", "collection"],
        how="left",
    )
    for column in ("pathways_replicated_up", "pathways_replicated_down"):
        if column not in summary:
            summary[column] = 0
        summary[column] = summary[column].fillna(0).astype(int)
    summary = summary.sort_values(
        ["model", "cell_type", "collection"],
        kind="stable",
    )
    summary.to_csv(
        RESULT_DIR / "mouse_pathway_replication_summary.tsv",
        sep="\t",
        index=False,
    )

    conditions = [
        (
            breadth["primary_subtypes_replicated"].ge(4)
            & breadth["p14_subtypes_replicated"].ge(3)
        ),
        (
            breadth["primary_subtypes_replicated"].ge(3)
            & breadth["p14_subtypes_replicated"].ge(2)
        ),
        (
            breadth["primary_subtypes_replicated"].ge(2)
            & breadth["p14_subtypes_replicated"].ge(1)
        ),
    ]
    breadth["evidence_tier"] = np.select(
        conditions,
        ["Tier_1", "Tier_2", "Tier_3"],
        default="Exploratory",
    )
    priority = breadth.loc[
        breadth["evidence_tier"].ne("Exploratory")
    ].copy()
    tier_order = pd.CategoricalDtype(
        ["Tier_1", "Tier_2", "Tier_3"],
        ordered=True,
    )
    priority["evidence_tier"] = priority["evidence_tier"].astype(tier_order)
    priority = priority.sort_values(
        [
            "evidence_tier",
            "primary_subtypes_replicated",
            "p14_subtypes_replicated",
            "bulk_padj",
        ],
        ascending=[True, False, False, True],
        kind="stable",
    )
    priority.to_csv(
        RESULT_DIR / "mouse_pathway_priority.tsv",
        sep="\t",
        index=False,
    )

    print(f"replication_rows={len(replication):,}")
    print(f"primary_replicated_rows={len(primary):,}")
    print(f"breadth_pathways={len(breadth):,}")
    print(
        priority.groupby(
            ["evidence_tier", "collection"],
            observed=True,
        ).size()
    )


if __name__ == "__main__":
    main()
