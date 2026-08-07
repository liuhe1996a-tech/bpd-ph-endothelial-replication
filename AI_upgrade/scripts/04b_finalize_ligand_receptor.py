"""Finalize ligand-receptor summaries after the streamed aggregation step.

This recovery entry point reads the animal-level evidence table written by
04_animal_balanced_ligand_receptor.py.  It is useful if rendering or external
annotation is interrupted after the expensive matrix aggregation has already
completed.  No statistical result is recomputed from cells.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def external_pair_matches(external: pd.DataFrame) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for ligand, receptor in zip(external["ligand"], external["receptor"]):
        ligand = str(ligand).title()
        parts = [p for p in str(receptor).replace("+", "_").split("_") if p]
        for part in parts:
            pairs.add((ligand, part.title()))
    return pairs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--external-cellchat", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--lr-resource", type=Path, required=True)
    args = parser.parse_args()

    evidence_path = args.output_dir / "GSE151974_animal_level_lr_evidence.tsv.gz"
    evidence = pd.read_csv(evidence_path, sep="\t")
    evidence["pair"] = evidence["ligand"] + "–" + evidence["receptor"]
    evidence["direction"] = np.where(
        evidence["communication_delta"] > 0, "increased", "decreased"
    )
    evidence.to_csv(evidence_path, sep="\t", index=False, compression="gzip")

    recurrence = (
        evidence.assign(
            nominal=evidence["exact_p"].le(0.05),
            increased=evidence["communication_delta"].gt(0),
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
    external_pairs = external_pair_matches(external)
    recurrence["reported_in_external_p53_cellchat"] = [
        (str(ligand).title(), str(receptor).title()) in external_pairs
        for ligand, receptor in zip(recurrence["ligand"], recurrence["receptor"])
    ]
    recurrence = recurrence.sort_values(
        ["reported_in_external_p53_cellchat", "priority_score", "median_delta"],
        ascending=[False, False, False],
    )
    recurrence.to_csv(
        args.output_dir / "GSE151974_recurrent_lr_summary.tsv", sep="\t", index=False
    )

    top = recurrence.query("positive_ages >= 2 and nominal_ages >= 1").head(24).copy()
    if top.empty:
        top = recurrence.query("median_delta > 0").head(24).copy()
    plot_data = evidence.merge(
        top[["source", "target", "ligand", "receptor", "pair"]],
        on=["source", "target", "ligand", "receptor", "pair"],
        how="inner",
    )
    plot_data["route"] = (
        plot_data["source"] + " → " + plot_data["target"] + ": " + plot_data["pair"]
    )
    order_routes = (
        top.assign(route=top["source"] + " → " + top["target"] + ": " + top["pair"])
        ["route"]
        .tolist()[::-1]
    )
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

    metadata = pd.read_csv(args.metadata, index_col=0, low_memory=False)
    recurrent = (recurrence["positive_ages"].ge(2) & recurrence["nominal_ages"].ge(1))
    audit = {
        "source_matrix_sha256": sha256(args.matrix),
        "source_metadata_sha256": sha256(args.metadata),
        "lr_resource_sha256": sha256(args.lr_resource),
        "external_cellchat_sha256": sha256(args.external_cellchat),
        "cells": int(len(metadata)),
        "animals": 36,
        "tests": int(len(evidence)),
        "recurrent_positive_pairs": int(recurrent.sum()),
        "external_overlap_pairs": int(
            (recurrent & recurrence["reported_in_external_p53_cellchat"]).sum()
        ),
        "interpretive_boundary": (
            "Expression-supported communication potential; no causal or "
            "physical-interaction claim."
        ),
    }
    (args.output_dir / "ligand_receptor_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
