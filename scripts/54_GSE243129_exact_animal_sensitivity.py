#!/usr/bin/env python3
"""Exact animal-subset sensitivity analysis for the small GSE243129 cohort.

GSE243129 contains two wild-type animals in each age-by-oxygen stratum.  This
script does not refit any model.  It reconstructs the held-out animal-level
normoxia and hyperoxia pseudobulks, combines them with the frozen R10
predictions, and evaluates every non-empty animal subset while applying each
age-specific selection jointly to Cap and Cap-a.  It also reports exact
one-pair-per-stratum and leave-one-animal analyses, plus the number of unique
unordered configurations represented by the original 500 bootstrap draws.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


METHODS = (
    "gene_linear",
    "pca_latent",
    "vae_latent",
    "cvae_counterfactual",
    "scgen_adapted",
    "cpa_adapted",
    "sinkhorn_ot",
)
AGES = ("P7", "P14")
CELL_TYPES = ("Cap", "Cap-a")
OXYGENS = ("Normoxia", "Hyperoxia")


def load_core(path: Path):
    spec = importlib.util.spec_from_file_location("virtual_cell_core", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load benchmark helpers from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_counts(path: Path) -> sparse.csr_matrix:
    from scipy.io import mmread

    if path.name.endswith(".npz"):
        return sparse.load_npz(path).tocsr()
    return mmread(path).tocsr()


def nonempty_subsets(n: int) -> list[tuple[int, ...]]:
    return [
        subset
        for size in range(1, n + 1)
        for subset in itertools.combinations(range(n), size)
    ]


def reconstruction(
    core,
    counts: sparse.csr_matrix,
    metadata: pd.DataFrame,
    manifest: pd.DataFrame,
    age: str,
    cell_type: str,
) -> tuple[np.ndarray, list[str], np.ndarray, list[str]]:
    fold = manifest.loc[manifest["heldout_age"].eq(age)].sort_values(
        "model_feature_position"
    )
    selected_idx = fold["gene_index"].to_numpy(dtype=int)
    x = core.normalize_selected_counts(counts, selected_idx)
    norm_mask = (
        metadata["Age"].eq(age)
        & metadata["Oxygen"].eq("Normoxia")
        & metadata["CellType"].eq(cell_type)
    ).to_numpy()
    hyp_mask = (
        metadata["Age"].eq(age)
        & metadata["Oxygen"].eq("Hyperoxia")
        & metadata["CellType"].eq(cell_type)
    ).to_numpy()
    norm, norm_animals = core.animal_means(
        x[norm_mask], metadata.loc[norm_mask].reset_index(drop=True)
    )
    hyp, hyp_animals = core.animal_means(
        x[hyp_mask], metadata.loc[hyp_mask].reset_index(drop=True)
    )
    return norm, norm_animals, hyp, hyp_animals


def prediction_entries(arrays: np.lib.npyio.NpzFile, age: str, cell_type: str):
    prefix = f"{age}_{cell_type.replace('-', 'a')}_"
    entries: list[tuple[str, int, np.ndarray]] = []
    for key in arrays.files:
        if not key.startswith(prefix):
            continue
        method_seed = key[len(prefix):]
        method, seed_text = method_seed.rsplit("_seed", 1)
        if method in METHODS:
            entries.append((method, int(seed_text), arrays[key]))
    return entries


def configuration_rows(
    core,
    arrays: np.lib.npyio.NpzFile,
    task_data: dict[tuple[str, str], dict[str, object]],
    configurations: list[dict[tuple[str, str], tuple[int, ...]]],
    analysis: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    all_indices = np.arange(1800, dtype=int)
    for config_id, config in enumerate(configurations, start=1):
        for age in AGES:
            norm_idx = np.asarray(config[(age, "Normoxia")], dtype=int)
            hyp_idx = np.asarray(config[(age, "Hyperoxia")], dtype=int)
            for cell_type in CELL_TYPES:
                data = task_data[(age, cell_type)]
                norm = data["norm"]
                hyp = data["hyp"]
                for method, seed, prediction in prediction_entries(arrays, age, cell_type):
                    metric = core.metric_row(
                        prediction[norm_idx], hyp[hyp_idx], norm[norm_idx], all_indices
                    )
                    rows.append(
                        {
                            "analysis": analysis,
                            "configuration": config_id,
                            "heldout_age": age,
                            "cell_type": cell_type,
                            "method": method,
                            "seed": seed,
                            "normoxia_subset": ";".join(
                                data["norm_animals"][i] for i in norm_idx
                            ),
                            "hyperoxia_subset": ";".join(
                                data["hyp_animals"][i] for i in hyp_idx
                            ),
                            **metric,
                        }
                    )
    return pd.DataFrame(rows)


def aggregate_and_rank(metrics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    seed_mean = (
        metrics.groupby(
            ["analysis", "configuration", "heldout_age", "cell_type", "method"],
            observed=True,
        )["spearman_effect"]
        .mean()
        .reset_index()
    )
    aggregate = (
        seed_mean.groupby(["analysis", "configuration", "method"], observed=True)[
            "spearman_effect"
        ]
        .mean()
        .reset_index(name="mean_spearman_across_four_tasks")
    )
    aggregate["rank"] = aggregate.groupby(
        ["analysis", "configuration"], observed=True
    )["mean_spearman_across_four_tasks"].rank(method="min", ascending=False)
    n_config = int(aggregate["configuration"].nunique())
    frequency = (
        aggregate.groupby(["analysis", "method"], observed=True)
        .agg(
            configurations=("configuration", "nunique"),
            mean_spearman=("mean_spearman_across_four_tasks", "mean"),
            minimum_spearman=("mean_spearman_across_four_tasks", "min"),
            maximum_spearman=("mean_spearman_across_four_tasks", "max"),
            median_rank=("rank", "median"),
            minimum_rank=("rank", "min"),
            maximum_rank=("rank", "max"),
            times_ranked_first=("rank", lambda x: int((x == 1).sum())),
            times_ranked_top3=("rank", lambda x: int((x <= 3).sum())),
        )
        .reset_index()
    )
    frequency["fraction_ranked_first"] = frequency["times_ranked_first"] / n_config
    frequency["fraction_ranked_top3"] = frequency["times_ranked_top3"] / n_config
    return aggregate, frequency


def bootstrap_configuration_audit(draws: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    signatures: list[dict[str, object]] = []
    for bootstrap, frame in draws.groupby("bootstrap", sort=True):
        parts = []
        for age in AGES:
            for oxygen in OXYGENS:
                animals = sorted(
                    frame.loc[
                        frame["heldout_age"].eq(age)
                        & frame["oxygen_group"].eq(oxygen),
                        "animal_id",
                    ].astype(str)
                )
                parts.append(f"{age}:{oxygen}:{','.join(animals)}")
        signatures.append(
            {"bootstrap": int(bootstrap), "unordered_configuration": "|".join(parts)}
        )
    ledger = pd.DataFrame(signatures)
    counts = (
        ledger.groupby("unordered_configuration", observed=True)
        .size()
        .reset_index(name="bootstrap_replicates")
        .sort_values(["bootstrap_replicates", "unordered_configuration"], ascending=[False, True])
        .reset_index(drop=True)
    )
    audit = {
        "bootstrap_replicates": int(len(ledger)),
        "unique_unordered_configurations": int(len(counts)),
        "maximum_possible_unordered_configurations": 81,
        "interpretation": (
            "The bootstrap quantifies resampling variability from eight animals; "
            "it does not create additional biological replication."
        ),
    }
    return counts, audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = (args.output_dir or root / "results" / "virtual_cell_r10_sensitivity").resolve()
    output.mkdir(parents=True, exist_ok=True)

    core = load_core(root / "scripts" / "27_virtual_cell_leakage_free_benchmark.py")
    raw = root / "raw" / "GSE243129_virtual_cell"
    expanded = root / "results" / "virtual_cell_benchmark_r10"
    counts = load_counts(raw / "GSE243129_WT_capillary_counts_cells_by_genes.mtx.gz")
    metadata = pd.read_csv(
        raw / "GSE243129_WT_capillary_cell_metadata.tsv.gz", sep="\t"
    )
    manifest = pd.read_csv(expanded / "GSE243129_fold_feature_manifest.tsv", sep="\t")
    arrays = np.load(expanded / "GSE243129_capillary_animal_arrays.npz")

    task_data: dict[tuple[str, str], dict[str, object]] = {}
    for age in AGES:
        for cell_type in CELL_TYPES:
            norm, norm_animals, hyp, hyp_animals = reconstruction(
                core, counts, metadata, manifest, age, cell_type
            )
            if len(norm_animals) != 2 or len(hyp_animals) != 2:
                raise ValueError(f"Expected two animals per state for {age} {cell_type}")
            task_data[(age, cell_type)] = {
                "norm": norm,
                "norm_animals": norm_animals,
                "hyp": hyp,
                "hyp_animals": hyp_animals,
            }

    choices = nonempty_subsets(2)
    exhaustive = [
        {
            ("P7", "Normoxia"): p7n,
            ("P7", "Hyperoxia"): p7h,
            ("P14", "Normoxia"): p14n,
            ("P14", "Hyperoxia"): p14h,
        }
        for p7n, p7h, p14n, p14h in itertools.product(choices, repeat=4)
    ]
    exact_metrics = configuration_rows(
        core, arrays, task_data, exhaustive, "all_nonempty_subsets"
    )
    exact_aggregate, exact_frequency = aggregate_and_rank(exact_metrics)

    singleton = [(0,), (1,)]
    one_pair = [
        {
            ("P7", "Normoxia"): p7n,
            ("P7", "Hyperoxia"): p7h,
            ("P14", "Normoxia"): p14n,
            ("P14", "Hyperoxia"): p14h,
        }
        for p7n, p7h, p14n, p14h in itertools.product(singleton, repeat=4)
    ]
    pair_metrics = configuration_rows(
        core, arrays, task_data, one_pair, "one_pair_per_age_oxygen_stratum"
    )
    pair_aggregate, pair_frequency = aggregate_and_rank(pair_metrics)

    full = {(age, oxygen): (0, 1) for age in AGES for oxygen in OXYGENS}
    leave_one: list[dict[tuple[str, str], tuple[int, ...]]] = []
    leave_labels: list[str] = []
    for age in AGES:
        for oxygen in OXYGENS:
            animals = task_data[(age, "Cap")][
                "norm_animals" if oxygen == "Normoxia" else "hyp_animals"
            ]
            for removed_idx, animal in enumerate(animals):
                config = dict(full)
                config[(age, oxygen)] = (1 - removed_idx,)
                leave_one.append(config)
                leave_labels.append(str(animal))
    leave_metrics = configuration_rows(
        core, arrays, task_data, leave_one, "leave_one_animal"
    )
    leave_metrics["removed_animal"] = leave_metrics["configuration"].map(
        {index + 1: animal for index, animal in enumerate(leave_labels)}
    )
    leave_aggregate, leave_frequency = aggregate_and_rank(leave_metrics)
    leave_aggregate["removed_animal"] = leave_aggregate["configuration"].map(
        {index + 1: animal for index, animal in enumerate(leave_labels)}
    )

    exact_metrics.to_csv(
        output / "GSE243129_all_nonempty_subset_metrics.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    exact_aggregate.to_csv(
        output / "GSE243129_all_nonempty_subset_method_ranks.tsv", sep="\t", index=False
    )
    exact_frequency.to_csv(
        output / "GSE243129_all_nonempty_subset_rank_frequency.tsv", sep="\t", index=False
    )
    pair_metrics.to_csv(
        output / "GSE243129_one_pair_metrics.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    pair_aggregate.to_csv(
        output / "GSE243129_one_pair_method_ranks.tsv", sep="\t", index=False
    )
    pair_frequency.to_csv(
        output / "GSE243129_one_pair_rank_frequency.tsv", sep="\t", index=False
    )
    leave_metrics.to_csv(
        output / "GSE243129_leave_one_animal_metrics.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    leave_aggregate.to_csv(
        output / "GSE243129_leave_one_animal_method_ranks.tsv", sep="\t", index=False
    )
    leave_frequency.to_csv(
        output / "GSE243129_leave_one_animal_rank_frequency.tsv", sep="\t", index=False
    )

    draws = pd.read_csv(
        expanded / "GSE243129_joint_capillary_bootstrap_draws.tsv.gz", sep="\t"
    )
    bootstrap_counts, bootstrap_audit = bootstrap_configuration_audit(draws)
    bootstrap_counts.to_csv(
        output / "GSE243129_bootstrap_unordered_configuration_counts.tsv",
        sep="\t", index=False,
    )
    audit = {
        "analysis_version": "R10 exact GSE243129 animal sensitivity",
        "animals_total": int(metadata["animal_id"].nunique()),
        "animals_per_age_by_oxygen": 2,
        "joint_application_to_cap_and_cap_a": True,
        "model_refitting_per_configuration": False,
        "all_nonempty_subset_configurations": int(len(exhaustive)),
        "one_pair_configurations": int(len(one_pair)),
        "leave_one_animal_configurations": int(len(leave_one)),
        **bootstrap_audit,
    }
    (output / "GSE243129_exact_animal_sensitivity_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, indent=2))
    print(exact_frequency.sort_values("median_rank").to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
