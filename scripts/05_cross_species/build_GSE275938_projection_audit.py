"""Build a row-level mouse-to-human projection audit from raw GSE275938 H5 files.

The audit separates four sequential questions:
1. Was a replicated mouse gene assigned a formal one-to-one human ortholog?
2. Is that human symbol present in the deposited H5 feature matrix?
3. Does the gene belong to the mouse >=3/5 endothelial-subtype signature?
4. Does it pass the prespecified human expression filter in each population?
"""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def decode(values: np.ndarray) -> list[str]:
    return [
        value.decode("utf-8")
        if isinstance(value, (bytes, np.bytes_))
        else str(value)
        for value in values
    ]


def sha256_strings(values: list[str]) -> str:
    return hashlib.sha256("\n".join(values).encode("utf-8")).hexdigest()


def parse_args() -> argparse.Namespace:
    default_project = Path(__file__).resolve().parents[2]
    default_h5 = os.environ.get(
        "GSE275938_H5_ROOT",
        str(default_project / "02_raw" / "GSE275938" / "GSE275938_RAW"),
    )
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=default_project)
    parser.add_argument("--h5-root", type=Path, default=Path(default_h5))
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    h5_root = args.h5_root.resolve()
    output = args.output or (
        project_root
        / "07_results"
        / "human_endothelial_subtypes"
        / "GSE275938_mouse_human_H5_expression_projection_audit.tsv"
    )

    h5_files = sorted(h5_root.glob("*.h5"))
    if len(h5_files) != 7:
        raise RuntimeError(
            f"Expected seven deposited H5 files in {h5_root}; found {len(h5_files)}."
        )

    feature_hashes: dict[str, str] = {}
    reference_features: list[str] | None = None
    for h5_path in h5_files:
        with h5py.File(h5_path, "r") as handle:
            features = decode(handle["matrix/features/name"][:])
        feature_hashes[h5_path.name] = sha256_strings(features)
        if reference_features is None:
            reference_features = features
        elif features != reference_features:
            raise RuntimeError(f"H5 feature order differs in {h5_path.name}.")
    assert reference_features is not None
    h5_symbols = set(reference_features)
    common_hash = next(iter(feature_hashes.values()))

    mapping_path = (
        project_root
        / "07_results"
        / "cross_species"
        / "NCBI_mouse_human_ortholog_mapping_audit.tsv"
    )
    effect_path = (
        project_root
        / "07_results"
        / "human_endothelial_subtypes"
        / "GSE275938_subtype_locked_gene_effects.tsv.gz"
    )
    mapping = pd.read_csv(mapping_path, sep="\t")
    effects = pd.read_csv(effect_path, sep="\t")
    effects = effects.loc[effects["normalization"].eq("TMM")].copy()

    audit = mapping.copy()
    audit.insert(0, "mouse_gene", audit.pop("gene"))
    audit["human_ortholog"] = audit["human_gene"]
    audit["one_to_one_ortholog"] = (
        audit["included_one_to_one_projection"].fillna(False).astype(bool)
    )
    audit["in_33_gene_set"] = audit["concordant_endothelial_subtypes"].ge(3)
    audit["present_in_H5"] = (
        audit["one_to_one_ortholog"]
        & audit["human_ortholog"].fillna("").isin(h5_symbols)
    )

    population_columns = {
        "All Endothelial": "expression_eligible_all_endothelium",
        "gCap": "expression_eligible_gCap",
        "aCap": "expression_eligible_aCap",
    }
    for population, column in population_columns.items():
        population_effects = (
            effects.loc[effects["subtype"].eq(population), ["gene", "eligible_expression"]]
            .drop_duplicates("gene")
            .set_index("gene")["eligible_expression"]
        )
        audit[column] = (
            audit["human_ortholog"].map(population_effects).fillna(False).astype(bool)
        )

    audit["H5_files_checked"] = len(h5_files)
    audit["H5_features_per_file"] = len(reference_features)
    audit["H5_unique_symbols"] = len(h5_symbols)
    audit["H5_feature_name_sha256"] = common_hash
    audit["H5_feature_identity_consistent_across_files"] = (
        len(set(feature_hashes.values())) == 1
    )
    audit["H5_audit_basis"] = (
        "direct read of matrix/features/name from all seven deposited H5 files"
    )

    leading = [
        "mouse_gene",
        "human_ortholog",
        "one_to_one_ortholog",
        "in_33_gene_set",
        "present_in_H5",
        "expression_eligible_all_endothelium",
        "expression_eligible_gCap",
        "expression_eligible_aCap",
        "concordant_endothelial_subtypes",
        "replicated_subtype_names",
        "H5_files_checked",
        "H5_features_per_file",
        "H5_unique_symbols",
        "H5_feature_name_sha256",
        "H5_feature_identity_consistent_across_files",
        "H5_audit_basis",
    ]
    remaining = [column for column in audit.columns if column not in leading]
    audit = audit[leading + remaining]
    output.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output, sep="\t", index=False)

    one_to_one = audit["one_to_one_ortholog"]
    signature = audit["in_33_gene_set"]
    summary = {
        "replicated_mouse_genes": int(len(audit)),
        "one_to_one_orthologs": int(one_to_one.sum()),
        "one_to_one_present_in_H5": int(
            (one_to_one & audit["present_in_H5"]).sum()
        ),
        "signature_genes": int(signature.sum()),
        "signature_present_in_H5": int(
            (signature & audit["present_in_H5"]).sum()
        ),
        **{
            column: int((signature & audit[column]).sum())
            for column in population_columns.values()
        },
    }
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
