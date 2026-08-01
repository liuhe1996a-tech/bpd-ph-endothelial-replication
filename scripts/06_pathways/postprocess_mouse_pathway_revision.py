"""Rebuild pathway replication summaries from completed fgsea results.

The fgsea calculations are written before summary construction in the R
workflow. This deterministic postprocessor avoids repeating those expensive
calculations when only summary-table logic changes.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULT_DIR = PROJECT_ROOT / "07_results" / "pathway_exact"
PRIMARY_MODEL = "age_by_oxygen_P14_contrast"


def joined(values: pd.Series) -> str:
    return ";".join(sorted(set(values.dropna().astype(str))))


def main() -> int:
    replication = pd.read_csv(
        RESULT_DIR / "mouse_pathway_replication_all.tsv.gz",
        sep="\t",
        low_memory=False,
    )
    for column in ["both_fdr05", "same_direction", "replicated"]:
        if replication[column].dtype == object:
            replication[column] = (
                replication[column].astype(str).str.upper().eq("TRUE")
            )

    primary = replication.loc[
        replication["model"].eq(PRIMARY_MODEL)
        & replication["replicated"]
    ].copy()
    primary = primary.sort_values(
        [
            "collection",
            "pathway",
            "leading_edge_overlap_n",
            "cell_type",
        ],
        ascending=[True, True, False, True],
    )
    primary.to_csv(
        RESULT_DIR / "mouse_pathway_replication_primary.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    rows: list[dict] = []
    for (collection, pathway), frame in primary.groupby(
        ["collection", "pathway"],
        sort=False,
    ):
        p14 = replication.loc[
            replication["collection"].eq(collection)
            & replication["pathway"].eq(pathway)
            & replication["model"].eq("P14_oxygen_only")
            & replication["replicated"]
        ]
        age_adjusted = replication.loc[
            replication["collection"].eq(collection)
            & replication["pathway"].eq(pathway)
            & replication["model"].eq("age_adjusted_all_ages")
            & replication["replicated"]
        ]
        rows.append(
            {
                "collection": collection,
                "pathway": pathway,
                "direction": frame["direction"].iloc[0],
                "primary_subtypes_replicated": frame[
                    "cell_type"
                ].nunique(),
                "primary_subtypes": joined(frame["cell_type"]),
                "bulk_NES": frame["bulk_NES"].iloc[0],
                "bulk_padj": frame["bulk_padj"].iloc[0],
                "median_endothelial_NES": frame[
                    "endothelial_NES"
                ].median(),
                "minimum_leading_edge_overlap": frame[
                    "leading_edge_overlap_n"
                ].min(),
                "median_leading_edge_overlap": frame[
                    "leading_edge_overlap_n"
                ].median(),
                "p14_subtypes_replicated": p14["cell_type"].nunique(),
                "p14_subtypes": joined(p14["cell_type"]),
                "age_adjusted_subtypes_replicated": age_adjusted[
                    "cell_type"
                ].nunique(),
                "age_adjusted_subtypes": joined(
                    age_adjusted["cell_type"]
                ),
            }
        )
    breadth = pd.DataFrame(rows).sort_values(
        [
            "primary_subtypes_replicated",
            "p14_subtypes_replicated",
            "age_adjusted_subtypes_replicated",
            "bulk_padj",
        ],
        ascending=[False, False, False, True],
    )
    breadth.to_csv(
        RESULT_DIR / "mouse_pathway_replication_breadth.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary = (
        replication.groupby(
            ["model", "cell_type", "collection"],
            as_index=False,
        )
        .agg(
            pathways_compared=("pathway", "size"),
            pathways_both_fdr05=("both_fdr05", "sum"),
            pathways_replicated_same_direction=("replicated", "sum"),
            pathways_replicated_up=(
                "direction",
                lambda values: (
                    replication.loc[values.index, "replicated"]
                    & values.eq("up")
                ).sum(),
            ),
            pathways_replicated_down=(
                "direction",
                lambda values: (
                    replication.loc[values.index, "replicated"]
                    & values.eq("down")
                ).sum(),
            ),
        )
        .sort_values(["model", "cell_type", "collection"])
    )
    summary.to_csv(
        RESULT_DIR / "mouse_pathway_replication_summary.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    print(
        f"Primary replicated pathways: {breadth.shape[0]}; "
        f"primary replicated rows: {primary.shape[0]}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
