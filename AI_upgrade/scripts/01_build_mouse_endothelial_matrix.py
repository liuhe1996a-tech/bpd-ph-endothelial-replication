"""Build a sparse endothelial single-cell matrix for virtual-cell modeling.

The deposited GSE151974 matrix is gene-by-cell CSV.  This script streams it
in bounded chunks, selects the five endothelial annotations used in the
manuscript, and writes a cell-by-gene sparse matrix plus immutable metadata.
Animal identity follows the lane-merging rule of the validated R5 workflow.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ENDOTHELIAL_TYPES = ("Cap", "Cap-a", "Art", "Vein", "Lymph")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=100)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(args.metadata, index_col=0, low_memory=False)
    if metadata.index.duplicated().any():
        raise ValueError("Metadata cell identifiers are not unique.")
    metadata["animal_source"] = metadata["orig.ident"].replace(
        {"P3_P7_1": "P3_P7", "P3_P7_2": "P3_P7"}
    )
    metadata["animal_id"] = (
        metadata["animal_source"].astype(str)
        + "_"
        + metadata["Barcode"].astype(str)
    )
    endothelial = metadata.loc[
        metadata["CellType"].isin(ENDOTHELIAL_TYPES)
    ].copy()
    selected_cells = endothelial.index.astype(str).tolist()
    selected_set = set(selected_cells)

    blocks: list[sparse.csr_matrix] = []
    genes: list[str] = []
    reader = pd.read_csv(
        args.matrix,
        index_col=0,
        usecols=lambda column: column == "Unnamed: 0" or column in selected_set,
        chunksize=args.chunk_size,
    )
    for chunk_number, chunk in enumerate(reader, start=1):
        chunk = chunk.loc[:, selected_cells]
        values = chunk.to_numpy(dtype=np.int32, copy=False)
        if np.any(values < 0):
            raise ValueError("Negative UMI counts detected.")
        blocks.append(sparse.csr_matrix(values))
        genes.extend(chunk.index.astype(str).tolist())
        if chunk_number % 25 == 0:
            print(f"genes_processed={len(genes)}", flush=True)

    if len(genes) != len(set(genes)):
        raise ValueError("Gene identifiers are not unique.")
    gene_by_cell = sparse.vstack(blocks, format="csr")
    counts = gene_by_cell.transpose().tocsr()
    counts.sort_indices()

    matrix_out = args.output_dir / "GSE151974_endothelial_counts_cells_by_genes.npz"
    metadata_out = args.output_dir / "GSE151974_endothelial_cell_metadata.tsv.gz"
    genes_out = args.output_dir / "GSE151974_endothelial_genes.tsv"
    audit_out = args.output_dir / "GSE151974_endothelial_sparse_matrix_audit.json"
    sparse.save_npz(matrix_out, counts, compressed=True)
    endothelial.loc[selected_cells].reset_index(names="cell_id").to_csv(
        metadata_out, sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame({"gene_index": np.arange(len(genes)), "gene": genes}).to_csv(
        genes_out, sep="\t", index=False
    )

    audit = {
        "source_matrix": str(args.matrix),
        "source_metadata": str(args.metadata),
        "source_matrix_sha256": sha256(args.matrix),
        "source_metadata_sha256": sha256(args.metadata),
        "cells": int(counts.shape[0]),
        "genes": int(counts.shape[1]),
        "nonzero_entries": int(counts.nnz),
        "matrix_density": float(counts.nnz / np.prod(counts.shape)),
        "animals": int(endothelial["animal_id"].nunique()),
        "ages": endothelial["Age"].value_counts().sort_index().to_dict(),
        "oxygen": endothelial["Oxygen"].value_counts().sort_index().to_dict(),
        "cell_types": endothelial["CellType"].value_counts().to_dict(),
        "output_matrix_sha256": sha256(matrix_out),
    }
    audit_out.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
