"""Animal-level ligand-receptor analysis for GSE151974.

The deposited gene-by-cell UMI matrix is streamed in bounded chunks. Counts
for OmniPath ligand/receptor genes are aggregated to animal-by-compartment
pseudobulks. Hyperoxia effects are then tested with exact label permutations
within each age. This avoids treating cells as biological replicates.

The analysis is deliberately described as expression-supported communication
potential: it does not claim physical ligand-receptor binding or causal
signalling in the absence of perturbation experiments.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns


LR_URL = "https://omnipathdb.org/interactions"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_fdr(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def broad_compartment(cell_type: str) -> str:
    if cell_type == "Cap":
        return "Cap"
    if cell_type == "Cap-a":
        return "Cap-a"
    if cell_type in {"Art", "Vein", "Lymph"}:
        return "Other endothelial"
    if "fibroblast" in cell_type.lower() or cell_type == "Myofibroblast":
        return "Fibroblast"
    if cell_type.startswith("Pericyte") or cell_type == "SMC":
        return "Mural"
    if cell_type.startswith("AT") or cell_type in {"Ciliated", "Club"}:
        return "Epithelial"
    if cell_type in {"Alv Mf", "Int Mf", "Mono", "Neut 1", "Neut 2", "DC1", "DC2"}:
        return "Myeloid"
    if any(token in cell_type for token in ("T cell", "B cell", "NK cell", "ILC2", "gd T")):
        return "Lymphoid"
    if cell_type == "Mast Ba2":
        return "Mast"
    if cell_type == "Mesothelial":
        return "Mesothelial"
    raise ValueError(f"Unmapped cell type: {cell_type}")


def load_or_download_lr(path: Path) -> pd.DataFrame:
    if not path.exists():
        response = requests.get(
            LR_URL,
            params={
                "datasets": "ligrecextra",
                "organisms": 10090,
                "genesymbols": 1,
                "format": "tsv",
            },
            timeout=180,
        )
        response.raise_for_status()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(response.content)
    lr = pd.read_csv(path, sep="\t")
    lr = lr.rename(
        columns={"source_genesymbol": "ligand", "target_genesymbol": "receptor"}
    )
    lr = lr.loc[:, ["ligand", "receptor", "is_directed"]].copy()
    lr = lr[lr["is_directed"].astype(str).str.lower().eq("true")]
    lr = lr.dropna().drop_duplicates(["ligand", "receptor"])
    lr = lr[
        lr["ligand"].str.match(r"^[A-Za-z][A-Za-z0-9._-]+$")
        & lr["receptor"].str.match(r"^[A-Za-z][A-Za-z0-9._-]+$")
    ]
    return lr.reset_index(drop=True)


def exact_two_group_p(scores: np.ndarray, observed_hyper: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return observed mean difference and exact two-sided permutation p values."""
    scores = np.asarray(scores, dtype=float)
    observed_hyper = np.asarray(observed_hyper, dtype=bool)
    n = scores.shape[1]
    n_hyper = int(observed_hyper.sum())
    if n != 12 or n_hyper != 6:
        raise ValueError(f"Expected 12 animals with 6 hyperoxia animals; got {n}, {n_hyper}.")
    observed = scores[:, observed_hyper].mean(1) - scores[:, ~observed_hyper].mean(1)
    permutations = np.array(
        [np.isin(np.arange(n), c) for c in itertools.combinations(range(n), n_hyper)],
        dtype=float,
    )
    hyper_means = scores @ permutations.T / n_hyper
    norm_means = scores @ (1.0 - permutations).T / (n - n_hyper)
    null = hyper_means - norm_means
    p = (np.abs(null) >= np.abs(observed[:, None]) - 1e-12).mean(1)
    return observed, p


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--lr-resource", type=Path, required=True)
    parser.add_argument("--external-cellchat", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=100)
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--min-logcpm", type=float, default=np.log1p(1.0))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    lr = load_or_download_lr(args.lr_resource)
    lr_genes = set(lr["ligand"]) | set(lr["receptor"])

    metadata = pd.read_csv(args.metadata, index_col=0, low_memory=False)
    metadata["animal_source"] = metadata["orig.ident"].replace(
        {"P3_P7_1": "P3_P7", "P3_P7_2": "P3_P7"}
    )
    metadata["animal_id"] = metadata["animal_source"].astype(str) + "_" + metadata["Barcode"].astype(str)
    metadata["compartment"] = metadata["CellType"].map(broad_compartment)
    metadata["group"] = metadata["animal_id"] + "||" + metadata["compartment"]
    cells = metadata.index.astype(str).tolist()

    group_order = sorted(metadata["group"].unique())
    group_to_idx = {g: i for i, g in enumerate(group_order)}
    group_index = metadata.loc[cells, "group"].map(group_to_idx).to_numpy()
    cell_counts = np.bincount(group_index, minlength=len(group_order))
    total_umi = metadata.groupby("group", observed=True)["nCount_RNA"].sum().reindex(group_order).to_numpy(float)

    selected_genes: list[str] = []
    aggregated_rows: list[np.ndarray] = []
    with gzip.open(args.matrix, "rt", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split(",")
        matrix_cells = header[1:]
        if set(matrix_cells) != set(cells):
            missing = len(set(cells) - set(matrix_cells))
            extra = len(set(matrix_cells) - set(cells))
            raise ValueError(f"Matrix/metadata cell mismatch: missing={missing}, extra={extra}")
        header_group_index = metadata.loc[matrix_cells, "group"].map(group_to_idx).to_numpy()
        for gene_number, line in enumerate(handle, start=1):
            gene, separator, payload = line.rstrip("\r\n").partition(",")
            if separator and gene in lr_genes:
                values = np.fromstring(payload, sep=",", dtype=np.int64)
                if values.size != len(matrix_cells):
                    raise ValueError(
                        f"Unexpected cell count for {gene}: {values.size} != {len(matrix_cells)}"
                    )
                aggregated = np.bincount(
                    header_group_index,
                    weights=values,
                    minlength=len(group_order),
                ).astype(np.int64)
                selected_genes.append(gene)
                aggregated_rows.append(aggregated)
            if gene_number % 1000 == 0:
                print(
                    f"genes_scanned={gene_number}; lr_genes_retained={len(selected_genes)}",
                    flush=True,
                )

    counts = pd.DataFrame(np.vstack(aggregated_rows), index=selected_genes, columns=group_order)
    if counts.index.duplicated().any():
        counts = counts.groupby(level=0).sum()
    logcpm = np.log1p(counts.div(total_umi, axis=1) * 1_000_000)

    group_meta = (
        metadata.groupby("group", observed=True)
        .agg(
            animal_id=("animal_id", "first"),
            age=("Age", "first"),
            oxygen=("Oxygen", "first"),
            compartment=("compartment", "first"),
            cells=("CellType", "size"),
        )
        .reindex(group_order)
        .reset_index()
    )
    group_meta["total_umi"] = total_umi
    group_meta.to_csv(args.output_dir / "GSE151974_compartment_pseudobulk_metadata.tsv", sep="\t", index=False)
    logcpm.T.to_csv(args.output_dir / "GSE151974_ligand_receptor_log1pCPM.tsv.gz", sep="\t", compression="gzip")

    present = set(logcpm.index)
    lr = lr[lr["ligand"].isin(present) & lr["receptor"].isin(present)].copy()
    results: list[pd.DataFrame] = []
    for age in ("P3", "P7", "P14"):
        animals = (
            metadata.loc[metadata["Age"].eq(age), ["animal_id", "Oxygen"]]
            .drop_duplicates()
            .sort_values("animal_id")
        )
        animal_ids = animals["animal_id"].tolist()
        hyper = animals["Oxygen"].eq("Hyperoxia").to_numpy()
        for target in ("Cap", "Cap-a"):
            target_groups = [f"{animal}||{target}" for animal in animal_ids]
            if not set(target_groups).issubset(logcpm.columns):
                continue
            receptor_values = logcpm.loc[lr["receptor"], target_groups].to_numpy()
            target_cells = group_meta.set_index("group").loc[target_groups, "cells"].to_numpy()
            for source in sorted(metadata["compartment"].unique()):
                source_groups = [f"{animal}||{source}" for animal in animal_ids]
                if not set(source_groups).issubset(logcpm.columns):
                    continue
                source_cells = group_meta.set_index("group").loc[source_groups, "cells"].to_numpy()
                adequate = (source_cells >= args.min_cells) & (target_cells >= args.min_cells)
                if adequate.sum() != len(animal_ids):
                    continue
                ligand_values = logcpm.loc[lr["ligand"], source_groups].to_numpy()
                expression_ok = (
                    (ligand_values[:, hyper] >= args.min_logcpm).sum(1) >= 4
                ) & (
                    (ligand_values[:, ~hyper] >= args.min_logcpm).sum(1) >= 4
                ) & (
                    (receptor_values[:, hyper] >= args.min_logcpm).sum(1) >= 4
                ) & (
                    (receptor_values[:, ~hyper] >= args.min_logcpm).sum(1) >= 4
                )
                if not expression_ok.any():
                    continue
                local = lr.loc[expression_ok].reset_index(drop=True)
                lv = ligand_values[expression_ok]
                rv = receptor_values[expression_ok]
                scores = np.sqrt(np.clip(lv, 0, None) * np.clip(rv, 0, None))
                delta, p = exact_two_group_p(scores, hyper)
                ligand_delta = lv[:, hyper].mean(1) - lv[:, ~hyper].mean(1)
                receptor_delta = rv[:, hyper].mean(1) - rv[:, ~hyper].mean(1)
                local = local.assign(
                    age=age,
                    source=source,
                    target=target,
                    communication_delta=delta,
                    ligand_delta=ligand_delta,
                    receptor_delta=receptor_delta,
                    exact_p=p,
                    min_source_cells=int(source_cells.min()),
                    min_target_cells=int(target_cells.min()),
                )
                local["fdr_bh"] = bh_fdr(local["exact_p"].to_numpy())
                results.append(local)

    evidence = pd.concat(results, ignore_index=True)
    evidence["direction"] = np.where(evidence["communication_delta"] > 0, "increased", "decreased")
    evidence["pair"] = evidence["ligand"] + "–" + evidence["receptor"]
    evidence.to_csv(args.output_dir / "GSE151974_animal_level_lr_evidence.tsv.gz", sep="\t", index=False, compression="gzip")

    recurrence = (
        evidence.assign(
            nominal=evidence["exact_p"].le(0.05),
            increased=evidence["communication_delta"].gt(0),
            pair=evidence["ligand"] + "–" + evidence["receptor"],
        )
        .groupby(["source", "target", "ligand", "receptor", "pair"], observed=True)
        .agg(
            ages_tested=("age", "nunique"),
            positive_ages=("increased", "sum"),
            nominal_ages=("nominal", "sum"),
            median_delta=("communication_delta", "median"),
            min_p=("exact_p", "min"),
            min_fdr=("fdr_bh", "min"),
        )
        .reset_index()
    )
    recurrence["priority_score"] = (
        recurrence["positive_ages"]
        * recurrence["nominal_ages"]
        * np.maximum(recurrence["median_delta"], 0)
        * -np.log10(np.maximum(recurrence["min_p"], 1 / 924))
    )

    external = pd.read_excel(args.external_cellchat, sheet_name=0, header=0)
    external.columns = [str(c).strip() for c in external.columns]
    external_pairs = set(
        zip(
            external.get("ligand", pd.Series(dtype=str)).dropna().astype(str).str.title(),
            external.get("receptor", pd.Series(dtype=str)).dropna().astype(str).str.replace(r"[_+].*$", "", regex=True).str.title(),
        )
    )
    recurrence["reported_in_external_p53_cellchat"] = [
        (str(ligand).title(), str(receptor).title()) in external_pairs
        for ligand, receptor in zip(recurrence["ligand"], recurrence["receptor"])
    ]
    recurrence = recurrence.sort_values(
        ["reported_in_external_p53_cellchat", "priority_score", "median_delta"],
        ascending=[False, False, False],
    )
    recurrence.to_csv(args.output_dir / "GSE151974_recurrent_lr_summary.tsv", sep="\t", index=False)

    top = recurrence.query("positive_ages >= 2 and nominal_ages >= 1").head(24).copy()
    if top.empty:
        top = recurrence.query("median_delta > 0").head(24).copy()
    plot_data = evidence.merge(top[["source", "target", "ligand", "receptor", "pair"]], on=["source", "target", "ligand", "receptor", "pair"], how="inner")
    plot_data["route"] = plot_data["source"] + " → " + plot_data["target"] + ": " + plot_data["pair"]
    order_routes = top.assign(route=top["source"] + " → " + top["target"] + ": " + top["pair"])["route"].tolist()[::-1]
    fig, ax = plt.subplots(figsize=(10.5, max(5.8, 0.27 * len(order_routes))))
    sns.scatterplot(
        data=plot_data,
        x="communication_delta",
        y="route",
        hue="age",
        size=-np.log10(plot_data["exact_p"].clip(lower=1 / 924)),
        sizes=(30, 130),
        palette={"P3": "#4C78A8", "P7": "#F58518", "P14": "#B279A2"},
        ax=ax,
    )
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_yticks(range(len(order_routes)), order_routes)
    ax.set_xlabel("Hyperoxia minus normoxia communication-potential score")
    ax.set_ylabel("")
    ax.set_title("Animal-level ligand–receptor evidence across postnatal ages")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(args.output_dir / "animal_level_ligand_receptor.png", dpi=300, bbox_inches="tight")
    fig.savefig(args.output_dir / "animal_level_ligand_receptor.pdf", bbox_inches="tight")
    plt.close(fig)

    audit = {
        "source_matrix_sha256": sha256(args.matrix),
        "source_metadata_sha256": sha256(args.metadata),
        "lr_resource_sha256": sha256(args.lr_resource),
        "external_cellchat_sha256": sha256(args.external_cellchat),
        "cells": int(len(metadata)),
        "animals": int(metadata["animal_id"].nunique()),
        "compartments": metadata["compartment"].value_counts().to_dict(),
        "lr_resource_pairs": int(len(lr)),
        "lr_genes_in_matrix": int(len(present)),
        "tests": int(len(evidence)),
        "recurrent_positive_pairs": int(((recurrence["positive_ages"] >= 2) & (recurrence["nominal_ages"] >= 1)).sum()),
        "interpretive_boundary": "Expression-supported communication potential; no causal or physical-interaction claim.",
    }
    (args.output_dir / "ligand_receptor_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
