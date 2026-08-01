"""Aggregate GSE151974 raw UMI counts by animal and endothelial subtype.

The 61,839-cell deposited matrix is read in bounded gene chunks. Counts are
summed to animal-level pseudobulks; cells are never treated as replicates.
"""

from __future__ import annotations

import os

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = Path(os.environ.get("GSE151974_RAW_ROOT", PROJECT_ROOT / "02_raw/GSE151974"))
MATRIX_PATH = SOURCE_ROOT / "GSE151974_raw_umi_matrix_postfilter.csv.gz"
METADATA_PATH = SOURCE_ROOT / "GSE151974_cell_metadata_postfilter.csv.gz"
PROCESSED_DIR = PROJECT_ROOT / "04_processed" / "GSE151974"
METADATA_DIR = PROJECT_ROOT / "03_metadata"
RESULT_DIR = PROJECT_ROOT / "07_results" / "input_qc"

ENDOTHELIAL_TYPES = ("Cap", "Cap-a", "Art", "Vein", "Lymph")
ENDOTHELIAL_FRACTION_LABELS = {
    "Cap": "capillary_general_like",
    "Cap-a": "capillary_aerocyte_like",
    "Art": "arterial",
    "Vein": "venous",
    "Lymph": "lymphatic",
}
CHUNK_SIZE = 100


def build_animal_composition(metadata: pd.DataFrame) -> pd.DataFrame:
    """Build the animal-level composition table used by Figure 3 tests."""
    animal_metadata = (
        metadata.groupby("animal_id", as_index=False)
        .agg(
            Age=("Age", "first"),
            Oxygen=("Oxygen", "first"),
            cells_total=("CellType", "size"),
            age_unique=("Age", "nunique"),
            oxygen_unique=("Oxygen", "nunique"),
        )
    )
    if not (
        animal_metadata["age_unique"].eq(1).all()
        and animal_metadata["oxygen_unique"].eq(1).all()
    ):
        raise ValueError("Animal identity maps to inconsistent Age/Oxygen labels.")

    counts = (
        metadata.groupby(["animal_id", "CellType"])
        .size()
        .unstack(fill_value=0)
        .reindex(animal_metadata["animal_id"])
        .fillna(0)
    )
    endothelial_total = counts.reindex(
        columns=ENDOTHELIAL_TYPES,
        fill_value=0,
    ).sum(axis=1)
    animal = animal_metadata.set_index("animal_id")
    animal["endothelial_cells"] = endothelial_total
    animal["endothelial_fraction_all_cells"] = (
        animal["endothelial_cells"] / animal["cells_total"]
    )
    for cell_type, label in ENDOTHELIAL_FRACTION_LABELS.items():
        animal[f"{cell_type}_cells"] = counts.get(cell_type, 0)
        animal[f"{label}_fraction_of_endothelium"] = (
            animal[f"{cell_type}_cells"]
            / animal["endothelial_cells"].replace(0, np.nan)
        )
    return animal.reset_index()


def main() -> int:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    METADATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(METADATA_PATH, index_col=0, low_memory=False)
    if metadata.index.duplicated().any():
        raise ValueError("GSE151974 metadata cell identifiers are not unique.")
    metadata["animal_source"] = metadata["orig.ident"].replace(
        {"P3_P7_1": "P3_P7", "P3_P7_2": "P3_P7"}
    )
    metadata["animal_id"] = (
        metadata["animal_source"].astype(str)
        + "_"
        + metadata["Barcode"].astype(str)
    )
    animal_composition = build_animal_composition(metadata)
    animal_composition.to_csv(
        METADATA_DIR / "GSE151974_animal_level_cell_composition.tsv",
        sep="\t",
        index=False,
    )
    endothelial = metadata.loc[
        metadata["CellType"].isin(ENDOTHELIAL_TYPES)
    ].copy()

    group_metadata = (
        endothelial.groupby(["CellType", "animal_id"], as_index=False)
        .agg(
            Age=("Age", "first"),
            Oxygen=("Oxygen", "first"),
            n_cells=("CellType", "size"),
            age_unique=("Age", "nunique"),
            oxygen_unique=("Oxygen", "nunique"),
        )
        .sort_values(["CellType", "animal_id"])
    )
    if not (
        (group_metadata["age_unique"] == 1).all()
        and (group_metadata["oxygen_unique"] == 1).all()
    ):
        raise ValueError("A pseudobulk maps to inconsistent Age/Oxygen labels.")
    group_metadata["pseudobulk_id"] = (
        group_metadata["CellType"].astype(str)
        + "__"
        + group_metadata["animal_id"].astype(str)
    )
    group_metadata = group_metadata.drop(
        columns=["age_unique", "oxygen_unique"]
    )

    group_lookup = {
        (row.CellType, row.animal_id): row.pseudobulk_id
        for row in group_metadata.itertuples()
    }
    endothelial["pseudobulk_id"] = [
        group_lookup[(cell_type, animal_id)]
        for cell_type, animal_id in zip(
            endothelial["CellType"],
            endothelial["animal_id"],
        )
    ]

    ordered_cells = endothelial.sort_values(
        ["pseudobulk_id"],
        kind="stable",
    ).index.tolist()
    selected_cell_set = set(ordered_cells)
    ordered_group_labels = (
        endothelial.loc[ordered_cells, "pseudobulk_id"].astype(str).to_numpy()
    )
    unique_groups, first_indices = np.unique(
        ordered_group_labels,
        return_index=True,
    )
    group_order = np.argsort(first_indices)
    unique_groups = unique_groups[group_order]
    group_to_index = {
        group: index for index, group in enumerate(unique_groups)
    }
    group_indices = np.array(
        [group_to_index[group] for group in ordered_group_labels],
        dtype=np.int32,
    )
    boundaries = np.flatnonzero(
        np.r_[True, group_indices[1:] != group_indices[:-1]]
    )
    if len(boundaries) != len(unique_groups):
        raise ValueError("Pseudobulk cell ordering is not contiguous.")

    with gzip.open(MATRIX_PATH, "rt", encoding="utf-8-sig") as handle:
        matrix_header = handle.readline().rstrip("\r\n").split(",")[1:]
    if len(matrix_header) != len(metadata):
        raise ValueError(
            f"Matrix has {len(matrix_header)} cells but metadata has "
            f"{len(metadata)} rows."
        )
    if len(set(matrix_header)) != len(matrix_header):
        raise ValueError("GSE151974 matrix cell identifiers are not unique.")
    missing_metadata = set(matrix_header) - set(metadata.index)
    missing_matrix = set(metadata.index) - set(matrix_header)
    if missing_metadata or missing_matrix:
        raise ValueError(
            "Matrix/metadata cell identifiers do not match exactly: "
            f"matrix_only={len(missing_metadata)}, "
            f"metadata_only={len(missing_matrix)}"
        )

    output_counts = (
        PROCESSED_DIR / "GSE151974_endothelial_animal_pseudobulk_counts.tsv.gz"
    )
    output_metadata = (
        PROCESSED_DIR / "GSE151974_endothelial_animal_pseudobulk_metadata.tsv"
    )
    output_qc = RESULT_DIR / "GSE151974_pseudobulk_input_audit.json"

    group_metadata.set_index("pseudobulk_id").loc[
        unique_groups
    ].reset_index().to_csv(
        output_metadata,
        sep="\t",
        index=False,
    )

    genes_processed = 0
    total_umi_by_group = np.zeros(len(unique_groups), dtype=np.int64)
    detected_genes_by_group = np.zeros(len(unique_groups), dtype=np.int64)
    with gzip.open(output_counts, "wt", encoding="utf-8", newline="") as output:
        first_chunk = True
        reader = pd.read_csv(
            MATRIX_PATH,
            index_col=0,
            usecols=lambda column: (
                column == "Unnamed: 0" or column in selected_cell_set
            ),
            chunksize=CHUNK_SIZE,
        )
        for chunk_number, chunk in enumerate(reader, start=1):
            chunk = chunk.loc[:, ordered_cells]
            values = chunk.to_numpy(dtype=np.int64, copy=False)
            if (values < 0).any():
                raise ValueError("Negative UMI counts detected.")
            aggregated = np.add.reduceat(values, boundaries, axis=1)
            aggregated_frame = pd.DataFrame(
                aggregated,
                index=chunk.index.astype(str),
                columns=unique_groups,
            )
            aggregated_frame.index.name = "gene"
            aggregated_frame.to_csv(
                output,
                sep="\t",
                header=first_chunk,
                index=True,
            )
            first_chunk = False
            genes_processed += len(chunk)
            total_umi_by_group += aggregated.sum(axis=0)
            detected_genes_by_group += (aggregated > 0).sum(axis=0)
            if chunk_number % 25 == 0:
                print(
                    f"genes_processed={genes_processed}",
                    flush=True,
                )

    qc_metadata = group_metadata.set_index("pseudobulk_id").loc[
        unique_groups
    ].copy()
    qc_metadata["total_umi"] = total_umi_by_group
    qc_metadata["detected_genes"] = detected_genes_by_group
    qc_metadata.reset_index().to_csv(
        RESULT_DIR / "GSE151974_pseudobulk_sample_QC.tsv",
        sep="\t",
        index=False,
    )

    audit = {
        "matrix_cells": len(matrix_header),
        "metadata_cells": len(metadata),
        "endothelial_cells": len(endothelial),
        "animals": int(metadata["animal_id"].nunique()),
        "pseudobulks": len(unique_groups),
        "genes": genes_processed,
        "matrix_metadata_exact_cell_id_match": True,
        "endothelial_cell_types": group_metadata.groupby("CellType").agg(
            pseudobulks=("pseudobulk_id", "size"),
            cells=("n_cells", "sum"),
            minimum_cells_per_pseudobulk=("n_cells", "min"),
            median_cells_per_pseudobulk=("n_cells", "median"),
            maximum_cells_per_pseudobulk=("n_cells", "max"),
        ).reset_index().to_dict(orient="records"),
        "output_counts": str(output_counts),
        "output_metadata": str(output_metadata),
    }
    output_qc.write_text(
        json.dumps(audit, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
