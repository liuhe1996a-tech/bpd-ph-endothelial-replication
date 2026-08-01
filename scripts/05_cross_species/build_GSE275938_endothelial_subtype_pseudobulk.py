"""Build donor-by-endothelial-subtype pseudobulks from GSE275938 H5 counts.

All inference and sensitivity summaries downstream use donors, not cells, as
the biological unit. Cell-level locked-gene matrices are retained only for
equal-cell downsampling diagnostics.
"""

from __future__ import annotations

import os

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PYTHON_RUNTIME = PROJECT_ROOT / ".python_runtime"
sys.path.insert(0, str(PYTHON_RUNTIME))
import h5py  # noqa: E402


SOURCE_ROOT = Path(os.environ.get("GSE275938_RAW_ROOT", PROJECT_ROOT / "02_raw/GSE275938"))
H5_ROOT = SOURCE_ROOT / "GSE275938_RAW"
METADATA_PATH = SOURCE_ROOT / "GSE275938_cell_metadata.csv.gz"
PROCESSED_DIR = PROJECT_ROOT / "04_processed" / "GSE275938"
RESULT_DIR = PROJECT_ROOT / "07_results" / "human_endothelial_subtypes"
LOCKED_EVIDENCE_PATH = (
    PROJECT_ROOT
    / "07_results"
    / "cross_species"
    / "locked_mouse_signature_human_BPD_PH_evidence.tsv"
)

SAMPLE_MAP = [
    (
        "Acute preterm injury 1",
        "GSM8488395_7466-1-filtered_feature_bc_matrix.h5",
        "Acute preterm injury",
    ),
    ("BPD 1", "GSM8488396_7mo-1-out_cellbender_filtered.h5", "BPD"),
    ("BPD 2", "GSM8488397_7mo-2-out_cellbender_filtered.h5", "BPD"),
    ("BPD+PH 1", "GSM8488398_5796-1-out_cellbender_filtered.h5", "BPD+PH"),
    ("BPD+PH 2", "GSM8488399_5796-2-out_cellbender_filtered.h5", "BPD+PH"),
    ("Term infant 1", "GSM8488400_0d-norm-out_cellbender_filtered.h5", "Term control"),
    ("Term infant 2", "GSM8488401_20d-norm-out_cellbender_filtered.h5", "Term control"),
]

ENDOTHELIAL_SUBTYPES = [
    "gCap",
    "aCap",
    "abCap",
    "Arterial EC",
    "Pulmonary venous EC",
    "Systemic venous EC",
    "Lymphatic",
]
ALL_ENDOTHELIAL = "All Endothelial"
DISEASE_CONDITIONS = {"BPD", "BPD+PH"}


def decode(values: np.ndarray) -> list[str]:
    return [
        value.decode("utf-8")
        if isinstance(value, (bytes, np.bytes_))
        else str(value)
        for value in values
    ]


def sha256_strings(values: list[str]) -> str:
    payload = "\n".join(values).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def collapse_duplicate_symbols(
    feature_names: list[str],
) -> tuple[list[str], np.ndarray]:
    symbol_to_position: dict[str, int] = {}
    unique_symbols: list[str] = []
    feature_to_symbol = np.empty(len(feature_names), dtype=np.int32)
    for feature_index, symbol in enumerate(feature_names):
        position = symbol_to_position.get(symbol)
        if position is None:
            position = len(unique_symbols)
            symbol_to_position[symbol] = position
            unique_symbols.append(symbol)
        feature_to_symbol[feature_index] = position
    return unique_symbols, feature_to_symbol


def main() -> int:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(METADATA_PATH, low_memory=False)
    required_columns = {
        "id",
        "dataset",
        "celltype",
        "celltype_lineage",
        "nCount_RNA",
        "nFeature_RNA",
        "percent.mt",
    }
    missing_columns = required_columns.difference(metadata.columns)
    if missing_columns:
        raise RuntimeError(f"Missing metadata columns: {sorted(missing_columns)}")
    if metadata["id"].duplicated().any():
        raise RuntimeError("Metadata cell IDs are not globally unique.")

    metadata = metadata.copy()
    metadata["base_barcode"] = (
        metadata["id"].astype(str).str.split("_", n=1).str[0]
    )
    duplicate_base_keys = int(
        metadata.duplicated(["dataset", "base_barcode"]).sum()
    )
    if duplicate_base_keys:
        raise RuntimeError(
            f"Duplicate dataset/base-barcode keys: {duplicate_base_keys}"
        )

    locked = pd.read_csv(LOCKED_EVIDENCE_PATH, sep="\t")
    locked = locked.dropna(subset=["human_gene"]).copy()
    locked_genes_requested = list(dict.fromkeys(locked["human_gene"].astype(str)))

    first_h5_path = H5_ROOT / SAMPLE_MAP[0][1]
    with h5py.File(first_h5_path, "r") as handle:
        first_features = decode(handle["matrix/features/name"][:])
        first_feature_ids = decode(handle["matrix/features/id"][:])
    unique_symbols, feature_to_symbol = collapse_duplicate_symbols(first_features)
    symbol_to_position = {
        symbol: position for position, symbol in enumerate(unique_symbols)
    }
    locked_genes_available = [
        gene for gene in locked_genes_requested if gene in symbol_to_position
    ]
    locked_to_position = {
        gene: position for position, gene in enumerate(locked_genes_available)
    }
    feature_to_locked = np.full(len(first_features), -1, dtype=np.int32)
    for feature_index, symbol in enumerate(first_features):
        if symbol in locked_to_position:
            feature_to_locked[feature_index] = locked_to_position[symbol]

    subtype_to_position = {
        subtype: position for position, subtype in enumerate(ENDOTHELIAL_SUBTYPES)
    }
    pseudobulk_counts: list[np.ndarray] = []
    pseudobulk_metadata_rows: list[dict] = []
    recovery_rows: list[dict] = []
    disease_cell_counts: list[np.ndarray] = []
    disease_cell_metadata_rows: list[dict] = []
    feature_hashes: list[dict] = []

    expected_feature_name_hash = sha256_strings(first_features)
    expected_feature_id_hash = sha256_strings(first_feature_ids)

    for sample, filename, condition in SAMPLE_MAP:
        sample_metadata = metadata.loc[metadata["dataset"].eq(sample)].copy()
        metadata_by_barcode = sample_metadata.set_index("base_barcode", drop=False)
        endothelial_metadata = sample_metadata.loc[
            sample_metadata["celltype_lineage"].eq("Endothelial")
        ].copy()
        unexpected_subtypes = sorted(
            set(endothelial_metadata["celltype"]).difference(ENDOTHELIAL_SUBTYPES)
        )
        if unexpected_subtypes:
            raise RuntimeError(
                f"Unexpected endothelial subtypes in {sample}: {unexpected_subtypes}"
            )

        h5_path = H5_ROOT / filename
        with h5py.File(h5_path, "r") as handle:
            matrix = handle["matrix"]
            feature_names = decode(matrix["features/name"][:])
            feature_ids = decode(matrix["features/id"][:])
            feature_name_hash = sha256_strings(feature_names)
            feature_id_hash = sha256_strings(feature_ids)
            if feature_name_hash != expected_feature_name_hash:
                raise RuntimeError(f"Feature-name order mismatch in {filename}")
            if feature_id_hash != expected_feature_id_hash:
                raise RuntimeError(f"Feature-ID order mismatch in {filename}")

            barcodes = decode(matrix["barcodes"][:])
            data = matrix["data"][:]
            indices = matrix["indices"][:]
            indptr = matrix["indptr"][:]

        feature_hashes.append(
            {
                "sample": sample,
                "filename": filename,
                "features": len(feature_names),
                "unique_symbols": len(set(feature_names)),
                "feature_name_sha256": feature_name_hash,
                "feature_id_sha256": feature_id_hash,
            }
        )

        h5_barcode_set = set(barcodes)
        metadata_barcode_set = set(sample_metadata["base_barcode"])
        endothelial_barcode_set = set(endothelial_metadata["base_barcode"])
        matched_metadata_barcodes = metadata_barcode_set.intersection(h5_barcode_set)
        matched_endothelial_barcodes = (
            endothelial_barcode_set.intersection(h5_barcode_set)
        )

        subtype_feature_counts = np.zeros(
            (len(ENDOTHELIAL_SUBTYPES), len(first_features)),
            dtype=np.int64,
        )
        subtype_cells = np.zeros(len(ENDOTHELIAL_SUBTYPES), dtype=np.int64)
        subtype_total_umi = np.zeros(len(ENDOTHELIAL_SUBTYPES), dtype=np.int64)

        for column, barcode in enumerate(barcodes):
            if barcode not in endothelial_barcode_set:
                continue
            row = metadata_by_barcode.loc[barcode]
            subtype = str(row["celltype"])
            subtype_position = subtype_to_position[subtype]
            start, end = int(indptr[column]), int(indptr[column + 1])
            gene_indices = indices[start:end]
            values = data[start:end].astype(np.int64, copy=False)
            np.add.at(
                subtype_feature_counts[subtype_position],
                gene_indices,
                values,
            )
            cell_total_umi = int(values.sum())
            subtype_cells[subtype_position] += 1
            subtype_total_umi[subtype_position] += cell_total_umi

            if condition in DISEASE_CONDITIONS:
                selected_positions = feature_to_locked[gene_indices]
                keep = selected_positions >= 0
                selected_counts = np.zeros(
                    len(locked_genes_available),
                    dtype=np.int64,
                )
                if keep.any():
                    np.add.at(
                        selected_counts,
                        selected_positions[keep],
                        values[keep],
                    )
                disease_cell_counts.append(selected_counts)
                disease_cell_metadata_rows.append(
                    {
                        "sample": sample,
                        "condition": condition,
                        "subtype": subtype,
                        "barcode": barcode,
                        "cell_total_umi": cell_total_umi,
                        "metadata_nCount_RNA": float(row["nCount_RNA"]),
                        "metadata_nFeature_RNA": float(row["nFeature_RNA"]),
                        "metadata_percent_mt": float(row["percent.mt"]),
                    }
                )

        collapsed_subtype_counts = np.zeros(
            (len(ENDOTHELIAL_SUBTYPES), len(unique_symbols)),
            dtype=np.int64,
        )
        for subtype_position in range(len(ENDOTHELIAL_SUBTYPES)):
            np.add.at(
                collapsed_subtype_counts[subtype_position],
                feature_to_symbol,
                subtype_feature_counts[subtype_position],
            )

        for subtype_position, subtype in enumerate(ENDOTHELIAL_SUBTYPES):
            subtype_meta = endothelial_metadata.loc[
                endothelial_metadata["celltype"].eq(subtype)
                & endothelial_metadata["base_barcode"].isin(h5_barcode_set)
            ]
            count_vector = collapsed_subtype_counts[subtype_position]
            pseudobulk_counts.append(count_vector)
            pseudobulk_metadata_rows.append(
                {
                    "pseudobulk_id": f"{sample}::{subtype}",
                    "sample": sample,
                    "condition": condition,
                    "subtype": subtype,
                    "cells": int(subtype_cells[subtype_position]),
                    "library_size": int(count_vector.sum()),
                    "genes_detected": int((count_vector > 0).sum()),
                    "median_metadata_nCount_RNA": float(
                        subtype_meta["nCount_RNA"].median()
                    ),
                    "median_metadata_nFeature_RNA": float(
                        subtype_meta["nFeature_RNA"].median()
                    ),
                    "median_metadata_percent_mt": float(
                        subtype_meta["percent.mt"].median()
                    ),
                    "source_h5": filename,
                }
            )

        all_counts = collapsed_subtype_counts.sum(axis=0)
        all_meta = endothelial_metadata.loc[
            endothelial_metadata["base_barcode"].isin(h5_barcode_set)
        ]
        pseudobulk_counts.append(all_counts)
        pseudobulk_metadata_rows.append(
            {
                "pseudobulk_id": f"{sample}::{ALL_ENDOTHELIAL}",
                "sample": sample,
                "condition": condition,
                "subtype": ALL_ENDOTHELIAL,
                "cells": int(subtype_cells.sum()),
                "library_size": int(all_counts.sum()),
                "genes_detected": int((all_counts > 0).sum()),
                "median_metadata_nCount_RNA": float(
                    all_meta["nCount_RNA"].median()
                ),
                "median_metadata_nFeature_RNA": float(
                    all_meta["nFeature_RNA"].median()
                ),
                "median_metadata_percent_mt": float(
                    all_meta["percent.mt"].median()
                ),
                "source_h5": filename,
            }
        )

        recovery_rows.append(
            {
                "sample": sample,
                "condition": condition,
                "metadata_cells": int(len(sample_metadata)),
                "h5_cells": int(len(barcodes)),
                "metadata_cells_matched_to_h5": int(
                    len(matched_metadata_barcodes)
                ),
                "metadata_to_h5_match_fraction": (
                    len(matched_metadata_barcodes) / len(sample_metadata)
                    if len(sample_metadata)
                    else np.nan
                ),
                "metadata_endothelial_cells": int(len(endothelial_metadata)),
                "metadata_endothelial_cells_matched_to_h5": int(
                    len(matched_endothelial_barcodes)
                ),
                "endothelial_match_fraction": (
                    len(matched_endothelial_barcodes)
                    / len(endothelial_metadata)
                    if len(endothelial_metadata)
                    else np.nan
                ),
                "unmatched_metadata_endothelial_cells": int(
                    len(endothelial_barcode_set.difference(h5_barcode_set))
                ),
            }
        )

        print(
            sample,
            f"metadata={len(sample_metadata):,}",
            f"h5={len(barcodes):,}",
            f"matched_EC={len(matched_endothelial_barcodes):,}",
        )

    counts_frame = pd.DataFrame(
        np.vstack(pseudobulk_counts),
        columns=unique_symbols,
    )
    pseudobulk_metadata = pd.DataFrame(pseudobulk_metadata_rows)
    counts_frame.insert(
        0,
        "pseudobulk_id",
        pseudobulk_metadata["pseudobulk_id"].to_numpy(),
    )
    counts_frame.to_csv(
        PROCESSED_DIR / "GSE275938_endothelial_subtype_pseudobulk_counts.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pseudobulk_metadata.to_csv(
        PROCESSED_DIR / "GSE275938_endothelial_subtype_pseudobulk_metadata.tsv",
        sep="\t",
        index=False,
    )

    disease_cell_metadata = pd.DataFrame(disease_cell_metadata_rows)
    disease_cell_metadata.to_csv(
        PROCESSED_DIR / "GSE275938_disease_endothelial_cell_metadata.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    np.savez_compressed(
        PROCESSED_DIR / "GSE275938_disease_endothelial_locked_gene_counts.npz",
        counts=np.vstack(disease_cell_counts),
        genes=np.asarray(locked_genes_available, dtype=str),
    )

    recovery = pd.DataFrame(recovery_rows)
    recovery.to_csv(
        RESULT_DIR / "GSE275938_cell_recovery_audit.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(feature_hashes).to_csv(
        RESULT_DIR / "GSE275938_H5_feature_identity_audit.tsv",
        sep="\t",
        index=False,
    )

    audit_summary = {
        "metadata_rows": int(len(metadata)),
        "metadata_columns": int(len(metadata.columns)),
        "metadata_id_duplicates": int(metadata["id"].duplicated().sum()),
        "dataset_base_barcode_duplicates": duplicate_base_keys,
        "h5_files": len(SAMPLE_MAP),
        "features": len(first_features),
        "unique_gene_symbols": len(unique_symbols),
        "duplicate_gene_symbols_collapsed": len(first_features) - len(unique_symbols),
        "locked_genes_requested": len(locked_genes_requested),
        "locked_genes_available": len(locked_genes_available),
        "disease_endothelial_cells_saved": int(len(disease_cell_metadata)),
        "feature_identity_consistent_across_h5": True,
    }
    (RESULT_DIR / "GSE275938_pseudobulk_build_summary.json").write_text(
        json.dumps(audit_summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(audit_summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
