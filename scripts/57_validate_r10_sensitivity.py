#!/usr/bin/env python3
"""Validate the exact R10 GSE243129 animal-sensitivity release."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    result = root / "results" / "virtual_cell_r10_sensitivity"
    output = (args.output_dir or root / "logs" / "r10_sensitivity_validation").resolve()
    output.mkdir(parents=True, exist_ok=True)

    checks: list[dict[str, object]] = []

    def add(name: str, passed: bool, detail: object) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})

    audit = json.loads(
        (result / "GSE243129_exact_animal_sensitivity_audit.json").read_text(
            encoding="utf-8"
        )
    )
    exact = pd.read_csv(
        result / "GSE243129_all_nonempty_subset_method_ranks.tsv", sep="\t"
    )
    pair = pd.read_csv(result / "GSE243129_one_pair_method_ranks.tsv", sep="\t")
    leave = pd.read_csv(
        result / "GSE243129_leave_one_animal_method_ranks.tsv", sep="\t"
    )
    frequency = pd.read_csv(
        result / "GSE243129_all_nonempty_subset_rank_frequency.tsv", sep="\t"
    )
    methods = {
        "gene_linear", "pca_latent", "vae_latent", "cvae_counterfactual",
        "scgen_adapted", "cpa_adapted", "sinkhorn_ot",
    }
    add("eight deposited animals", audit["animals_total"] == 8, audit["animals_total"])
    add(
        "two animals per age-by-oxygen group",
        audit["animals_per_age_by_oxygen"] == 2,
        audit["animals_per_age_by_oxygen"],
    )
    add("81 exact non-empty configurations", exact["configuration"].nunique() == 81, exact["configuration"].nunique())
    add("16 one-pair configurations", pair["configuration"].nunique() == 16, pair["configuration"].nunique())
    add("8 leave-one-animal configurations", leave["configuration"].nunique() == 8, leave["configuration"].nunique())
    add("seven methods in exact analysis", set(exact["method"]) == methods, sorted(exact["method"].unique()))
    add("exact rank keys unique", not exact.duplicated(["configuration", "method"]).any(), int(exact.duplicated(["configuration", "method"]).sum()))
    add("one-pair rank keys unique", not pair.duplicated(["configuration", "method"]).any(), int(pair.duplicated(["configuration", "method"]).sum()))
    add("leave-one rank keys unique", not leave.duplicated(["configuration", "method"]).any(), int(leave.duplicated(["configuration", "method"]).sum()))
    scgen = frequency.loc[frequency["method"].eq("scgen_adapted")].iloc[0]
    add("scGen-style first in 67 of 81 exact subsets", int(scgen["times_ranked_first"]) == 67, int(scgen["times_ranked_first"]))
    add(
        "500 bootstrap draws reduce to 77 unordered configurations",
        audit["bootstrap_replicates"] == 500
        and audit["unique_unordered_configurations"] == 77,
        f"replicates={audit['bootstrap_replicates']}; unique={audit['unique_unordered_configurations']}",
    )
    add(
        "animal selections jointly applied to Cap and Cap-a",
        bool(audit["joint_application_to_cap_and_cap_a"]),
        audit["joint_application_to_cap_and_cap_a"],
    )
    table = pd.DataFrame(checks)
    table.to_csv(output / "R10_sensitivity_validation_checks.tsv", sep="\t", index=False)
    summary = {
        "checks": int(len(table)),
        "passed": int(table["passed"].sum()),
        "failed": int((~table["passed"]).sum()),
        "all_passed": bool(table["passed"].all()),
    }
    (output / "R10_sensitivity_validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(table.to_string(index=False))
    print(json.dumps(summary, indent=2))
    return 0 if summary["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
