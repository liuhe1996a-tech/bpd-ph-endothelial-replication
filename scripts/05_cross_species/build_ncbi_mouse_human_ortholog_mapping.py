"""Build a version-fixed one-to-one mouse-to-human ortholog projection.

Source: NCBI Gene ``gene_orthologs`` and species ``gene_info`` files
downloaded on 2026-07-30. Mapping cardinality is evaluated across all
mouse-human pairs before the one-to-one restriction is applied.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXTERNAL_DIR = PROJECT_ROOT / "external_data"
MOUSE_DIR = PROJECT_ROOT / "07_results" / "mouse_endothelium_exact"
RESULT_DIR = PROJECT_ROOT / "07_results" / "cross_species"
ACCESS_DATE = "2026-07-30"
MOUSE_TAX_ID = 10090
HUMAN_TAX_ID = 9606


def read_gene_info(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        dtype={"GeneID": "int64", "Symbol": "string"},
        low_memory=False,
    )
    return frame[
        ["GeneID", "Symbol", "Synonyms", "type_of_gene"]
    ].drop_duplicates("GeneID")


def resolve_mouse_symbols(
    signature: pd.DataFrame,
    mouse_info: pd.DataFrame,
) -> pd.DataFrame:
    current = mouse_info.rename(
        columns={
            "GeneID": "mouse_GeneID",
            "Symbol": "mouse_ncbi_symbol",
            "Synonyms": "mouse_ncbi_synonyms",
            "type_of_gene": "mouse_gene_type",
        }
    )
    exact = signature.merge(
        current,
        left_on="gene",
        right_on="mouse_ncbi_symbol",
        how="left",
    )
    exact["mouse_symbol_resolution"] = np.where(
        exact["mouse_GeneID"].notna(),
        "exact current NCBI Symbol",
        "unresolved",
    )

    synonym_rows = []
    for _, row in current.iterrows():
        synonyms = str(row["mouse_ncbi_synonyms"])
        if synonyms in {"-", "nan", "<NA>"}:
            continue
        for synonym in synonyms.split("|"):
            if synonym:
                synonym_rows.append(
                    {
                        "input_symbol": synonym,
                        "mouse_GeneID_synonym": row["mouse_GeneID"],
                        "mouse_ncbi_symbol_synonym": row[
                            "mouse_ncbi_symbol"
                        ],
                        "mouse_ncbi_synonyms_synonym": row[
                            "mouse_ncbi_synonyms"
                        ],
                        "mouse_gene_type_synonym": row["mouse_gene_type"],
                    }
                )
    synonym_map = pd.DataFrame(synonym_rows)
    synonym_counts = synonym_map.groupby("input_symbol")[
        "mouse_GeneID_synonym"
    ].transform("nunique")
    synonym_map = synonym_map.loc[synonym_counts.eq(1)].drop_duplicates(
        "input_symbol"
    )
    exact = exact.merge(
        synonym_map,
        left_on="gene",
        right_on="input_symbol",
        how="left",
    )
    use_synonym = exact["mouse_GeneID"].isna() & exact[
        "mouse_GeneID_synonym"
    ].notna()
    for target, source in {
        "mouse_GeneID": "mouse_GeneID_synonym",
        "mouse_ncbi_symbol": "mouse_ncbi_symbol_synonym",
        "mouse_ncbi_synonyms": "mouse_ncbi_synonyms_synonym",
        "mouse_gene_type": "mouse_gene_type_synonym",
    }.items():
        exact.loc[use_synonym, target] = exact.loc[use_synonym, source]
    exact.loc[
        use_synonym,
        "mouse_symbol_resolution",
    ] = "unique match in NCBI Synonyms"
    return exact.drop(
        columns=[
            "input_symbol",
            "mouse_GeneID_synonym",
            "mouse_ncbi_symbol_synonym",
            "mouse_ncbi_synonyms_synonym",
            "mouse_gene_type_synonym",
        ]
    )


def read_mouse_human_pairs(path: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    columns = [
        "#tax_id",
        "GeneID",
        "relationship",
        "Other_tax_id",
        "Other_GeneID",
    ]
    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        dtype={
            "#tax_id": "int32",
            "GeneID": "int64",
            "relationship": "string",
            "Other_tax_id": "int32",
            "Other_GeneID": "int64",
        },
        chunksize=1_000_000,
    ):
        forward = chunk.loc[
            chunk["#tax_id"].eq(MOUSE_TAX_ID)
            & chunk["Other_tax_id"].eq(HUMAN_TAX_ID),
            columns,
        ].copy()
        if not forward.empty:
            forward = forward.rename(
                columns={
                    "GeneID": "mouse_GeneID",
                    "Other_GeneID": "human_GeneID",
                }
            )
            frames.append(
                forward[
                    ["mouse_GeneID", "human_GeneID", "relationship"]
                ]
            )
        reverse = chunk.loc[
            chunk["#tax_id"].eq(HUMAN_TAX_ID)
            & chunk["Other_tax_id"].eq(MOUSE_TAX_ID),
            columns,
        ].copy()
        if not reverse.empty:
            reverse = reverse.rename(
                columns={
                    "GeneID": "human_GeneID",
                    "Other_GeneID": "mouse_GeneID",
                }
            )
            frames.append(
                reverse[
                    ["mouse_GeneID", "human_GeneID", "relationship"]
                ]
            )
    pairs = pd.concat(frames, ignore_index=True).drop_duplicates(
        ["mouse_GeneID", "human_GeneID"]
    )
    pairs["mouse_to_human_n"] = pairs.groupby("mouse_GeneID")[
        "human_GeneID"
    ].transform("nunique")
    pairs["human_to_mouse_n"] = pairs.groupby("human_GeneID")[
        "mouse_GeneID"
    ].transform("nunique")
    return pairs


def build_mouse_signature() -> pd.DataFrame:
    convergent = pd.read_csv(
        MOUSE_DIR / "mouse_endothelial_convergent_genes.tsv",
        sep="\t",
    )
    convergent = convergent.loc[
        convergent["direction_concordant"]
    ].copy()
    return (
        convergent.groupby("gene", as_index=False)
        .agg(
            mouse_bulk_log2FC=(
                "bulk_log2FC_hyperoxia_vs_air",
                "first",
            ),
            mouse_bulk_fdr=("bulk_fdr_bh", "first"),
            mouse_endothelial_mean_log2FC=(
                "single_cell_logFC_hyperoxia",
                "mean",
            ),
            mouse_endothelial_best_fdr=(
                "single_cell_local_fdr",
                "min",
            ),
            concordant_endothelial_subtypes=(
                "concordant_endothelial_subtypes",
                "max",
            ),
            replicated_subtype_names=(
                "CellType",
                lambda values: ";".join(sorted(set(values))),
            ),
        )
        .sort_values(
            [
                "concordant_endothelial_subtypes",
                "mouse_bulk_fdr",
                "gene",
            ],
            ascending=[False, True, True],
        )
    )


def main() -> int:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    mouse_info = read_gene_info(
        EXTERNAL_DIR / "Mus_musculus.gene_info_20260730.gz"
    ).rename(
        columns={
            "GeneID": "mouse_GeneID",
            "Symbol": "mouse_ncbi_symbol",
            "Synonyms": "mouse_ncbi_synonyms",
            "type_of_gene": "mouse_gene_type",
        }
    )
    human_info = read_gene_info(
        EXTERNAL_DIR / "Homo_sapiens.gene_info_20260730.gz"
    ).rename(
        columns={
            "GeneID": "human_GeneID",
            "Symbol": "human_gene",
            "Synonyms": "human_ncbi_synonyms",
            "type_of_gene": "human_gene_type",
        }
    )
    pairs = read_mouse_human_pairs(
        EXTERNAL_DIR / "gene_orthologs_20260730.gz"
    )
    pairs = (
        pairs.merge(mouse_info, on="mouse_GeneID", how="left")
        .merge(human_info, on="human_GeneID", how="left")
    )

    signature = build_mouse_signature()
    signature = resolve_mouse_symbols(
        signature,
        mouse_info.rename(
            columns={
                "mouse_GeneID": "GeneID",
                "mouse_ncbi_symbol": "Symbol",
                "mouse_ncbi_synonyms": "Synonyms",
                "mouse_gene_type": "type_of_gene",
            }
        ),
    )
    mapped = signature.merge(
        pairs[
            [
                "mouse_GeneID",
                "human_GeneID",
                "relationship",
                "mouse_to_human_n",
                "human_to_mouse_n",
                "human_gene",
                "human_ncbi_synonyms",
                "human_gene_type",
            ]
        ],
        on="mouse_GeneID",
        how="left",
    )
    mapped["orthology_cardinality"] = "unmapped"
    has_pair = mapped["human_GeneID"].notna()
    mapped.loc[
        has_pair
        & mapped["mouse_to_human_n"].eq(1)
        & mapped["human_to_mouse_n"].eq(1),
        "orthology_cardinality",
    ] = "one-to-one"
    mapped.loc[
        has_pair
        & mapped["mouse_to_human_n"].gt(1)
        & mapped["human_to_mouse_n"].eq(1),
        "orthology_cardinality",
    ] = "one-to-many"
    mapped.loc[
        has_pair
        & mapped["mouse_to_human_n"].eq(1)
        & mapped["human_to_mouse_n"].gt(1),
        "orthology_cardinality",
    ] = "many-to-one"
    mapped.loc[
        has_pair
        & mapped["mouse_to_human_n"].gt(1)
        & mapped["human_to_mouse_n"].gt(1),
        "orthology_cardinality",
    ] = "many-to-many"
    mapped["included_one_to_one_projection"] = mapped[
        "orthology_cardinality"
    ].eq("one-to-one")
    mapped["exclusion_reason"] = ""
    mapped.loc[
        mapped["mouse_GeneID"].isna(),
        "exclusion_reason",
    ] = "mouse symbol absent from NCBI Mus musculus gene_info"
    mapped.loc[
        mapped["mouse_GeneID"].notna() & mapped["human_GeneID"].isna(),
        "exclusion_reason",
    ] = "no NCBI mouse-human ortholog"
    mapped.loc[
        mapped["orthology_cardinality"].isin(
            ["one-to-many", "many-to-one", "many-to-many"]
        ),
        "exclusion_reason",
    ] = (
        "non-one-to-one NCBI orthology "
        + mapped["orthology_cardinality"].astype(str)
    )
    mapped["mapping_source"] = "NCBI Gene gene_orthologs"
    mapped["mapping_access_date"] = ACCESS_DATE
    mapped["mouse_gene_info_file"] = (
        "Mus_musculus.gene_info downloaded 2026-07-30"
    )
    mapped["human_gene_info_file"] = (
        "Homo_sapiens.gene_info downloaded 2026-07-30"
    )

    one_to_one = (
        mapped.loc[mapped["included_one_to_one_projection"]]
        .drop_duplicates("gene")
        .copy()
    )
    one_to_one["signature_definition"] = (
        "P14 age-by-oxygen mouse replication; GSE216046 and "
        "GSE151974 FDR<0.05, |observed log2FC|>=1, same direction"
    )
    one_to_one["projection_status"] = (
        "sequentially defined and frozen for current human projection"
    )

    mapped.to_csv(
        RESULT_DIR / "NCBI_mouse_human_ortholog_mapping_audit.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    one_to_one.to_csv(
        RESULT_DIR / "formal_one_to_one_mouse_signature.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    # Downstream scripts use this stable path; it now contains the formal
    # one-to-one projection rather than manual capitalization rules.
    one_to_one.to_csv(
        RESULT_DIR / "locked_mouse_signature_human_BPD_PH_evidence.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    per_gene = mapped.drop_duplicates("gene")
    summary = {
        "mapping_source": "NCBI Gene",
        "ortholog_file": "gene_orthologs downloaded 2026-07-30",
        "access_date": ACCESS_DATE,
        "mouse_signature_genes_before_orthology": int(
            signature["gene"].nunique()
        ),
        "one_to_one_genes": int(one_to_one["gene"].nunique()),
        "genes_without_any_pair": int(
            per_gene["human_GeneID"].isna().sum()
        ),
        "genes_excluded_for_non_one_to_one_mapping": int(
            per_gene["orthology_cardinality"]
            .isin(["one-to-many", "many-to-one", "many-to-many"])
            .sum()
        ),
        "high_confidence_mouse_genes_ge3_subtypes": int(
            signature["concordant_endothelial_subtypes"].ge(3).sum()
        ),
        "high_confidence_one_to_one_genes_ge3_subtypes": int(
            one_to_one["concordant_endothelial_subtypes"].ge(3).sum()
        ),
    }
    (
        RESULT_DIR / "NCBI_mouse_human_ortholog_mapping_summary.json"
    ).write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
