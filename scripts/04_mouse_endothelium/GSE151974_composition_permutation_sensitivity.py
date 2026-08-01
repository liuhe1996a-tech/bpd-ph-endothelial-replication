"""Animal-level composition sensitivities using logit, CLR and exact labels."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
INPUT_PATH = (
    PROJECT_ROOT
    / "03_metadata"
    / "GSE151974_animal_level_cell_composition.tsv"
)
OUTPUT_DIR = PROJECT_ROOT / "07_results" / "mouse_endothelium_exact"
AGES = ["P3", "P7", "P14"]
BOOTSTRAP_SEED = 20260731
BOOTSTRAP_REPEATS = 10_000
SUBTYPE_COUNTS = {
    "Cap": "Cap_cells",
    "Cap-a": "Cap-a_cells",
    "Art": "Art_cells",
    "Vein": "Vein_cells",
    "Lymph": "Lymph_cells",
}
METRICS = {
    "endothelial_fraction_all_cells": (
        "endothelial_cells",
        "cells_total",
    ),
    **{
        f"{subtype}_fraction_of_endothelium": (
            column,
            "endothelial_cells",
        )
        for subtype, column in SUBTYPE_COUNTS.items()
    },
}


def bh_adjust(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = numeric.dropna().sort_values()
    ranks = np.arange(1, len(valid) + 1)
    adjusted = valid.to_numpy() * len(valid) / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[valid.index] = np.minimum(adjusted, 1.0)
    return result


def all_group_assignments(n: int, group_size: int) -> list[np.ndarray]:
    assignments = []
    for selected in combinations(range(n), group_size):
        mask = np.zeros(n, dtype=bool)
        mask[list(selected)] = True
        assignments.append(mask)
    return assignments


def exact_mean_difference(
    values: np.ndarray,
    observed_hyperoxia: np.ndarray,
) -> tuple[float, float, int]:
    observed = float(
        values[observed_hyperoxia].mean()
        - values[~observed_hyperoxia].mean()
    )
    assignments = all_group_assignments(
        len(values),
        int(observed_hyperoxia.sum()),
    )
    effects = np.asarray(
        [
            values[mask].mean() - values[~mask].mean()
            for mask in assignments
        ]
    )
    p_value = float(np.mean(np.abs(effects) >= abs(observed) - 1e-15))
    return observed, p_value, len(assignments)


def stratified_bootstrap_mean_difference(
    normoxia: np.ndarray,
    hyperoxia: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float]:
    """Percentile interval for the raw mean difference.

    Animals are resampled with replacement within their observed oxygen group.
    This interval is descriptive; inference and multiplicity control use the
    exact permutation P values below.
    """
    normoxia_indices = rng.integers(
        0,
        len(normoxia),
        size=(BOOTSTRAP_REPEATS, len(normoxia)),
    )
    hyperoxia_indices = rng.integers(
        0,
        len(hyperoxia),
        size=(BOOTSTRAP_REPEATS, len(hyperoxia)),
    )
    effects = (
        hyperoxia[hyperoxia_indices].mean(axis=1)
        - normoxia[normoxia_indices].mean(axis=1)
    )
    low, high = np.quantile(effects, [0.025, 0.975])
    return float(low), float(high)


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    animal = pd.read_csv(INPUT_PATH, sep="\t")
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    logit_rows: list[dict] = []
    clr_rows: list[dict] = []
    clr_donor_rows: list[dict] = []

    for age in AGES:
        frame = animal.loc[animal["Age"].eq(age)].copy().reset_index(drop=True)
        observed_hyperoxia = frame["Oxygen"].eq("Hyperoxia").to_numpy()
        for metric, (numerator, denominator) in METRICS.items():
            raw_fraction = (
                frame[numerator].to_numpy(dtype=float)
                / frame[denominator].to_numpy(dtype=float)
            )
            adjusted_fraction = (
                frame[numerator].to_numpy(dtype=float) + 0.5
            ) / (frame[denominator].to_numpy(dtype=float) + 1.0)
            logit = np.log(adjusted_fraction / (1 - adjusted_fraction))
            effect, p_value, permutations = exact_mean_difference(
                logit,
                observed_hyperoxia,
            )
            ci_low, ci_high = stratified_bootstrap_mean_difference(
                raw_fraction[~observed_hyperoxia],
                raw_fraction[observed_hyperoxia],
                rng,
            )
            logit_rows.append(
                {
                    "Age": age,
                    "metric": metric,
                    "n_normoxia": int((~observed_hyperoxia).sum()),
                    "n_hyperoxia": int(observed_hyperoxia.sum()),
                    "mean_raw_fraction_normoxia": float(
                        raw_fraction[~observed_hyperoxia].mean()
                    ),
                    "mean_raw_fraction_hyperoxia": float(
                        raw_fraction[observed_hyperoxia].mean()
                    ),
                    "raw_fraction_difference_hyperoxia_minus_normoxia": float(
                        raw_fraction[observed_hyperoxia].mean()
                        - raw_fraction[~observed_hyperoxia].mean()
                    ),
                    "raw_difference_ci95_low": ci_low,
                    "raw_difference_ci95_high": ci_high,
                    "raw_difference_interval_method": (
                        "within-oxygen-group nonparametric bootstrap percentile "
                        f"interval; {BOOTSTRAP_REPEATS} repeats"
                    ),
                    "bootstrap_seed": BOOTSTRAP_SEED,
                    "continuity_correction": (
                        "(numerator+0.5)/(denominator+1)"
                    ),
                    "mean_logit_normoxia": float(
                        logit[~observed_hyperoxia].mean()
                    ),
                    "mean_logit_hyperoxia": float(
                        logit[observed_hyperoxia].mean()
                    ),
                    "logit_difference_hyperoxia_minus_normoxia": effect,
                    "exact_two_sided_permutation_p": p_value,
                    "label_assignments": permutations,
                }
            )

        count_matrix = frame[list(SUBTYPE_COUNTS.values())].to_numpy(
            dtype=float
        )
        log_counts = np.log(count_matrix + 0.5)
        clr = log_counts - log_counts.mean(axis=1, keepdims=True)
        observed_vector = (
            clr[observed_hyperoxia].mean(axis=0)
            - clr[~observed_hyperoxia].mean(axis=0)
        )
        observed_distance = float(np.sum(observed_vector**2))
        assignments = all_group_assignments(
            len(frame),
            int(observed_hyperoxia.sum()),
        )
        permuted_distances = np.asarray(
            [
                np.sum(
                    (
                        clr[mask].mean(axis=0)
                        - clr[~mask].mean(axis=0)
                    )
                    ** 2
                )
                for mask in assignments
            ]
        )
        global_p = float(
            np.mean(permuted_distances >= observed_distance - 1e-15)
        )
        for position, subtype in enumerate(SUBTYPE_COUNTS):
            effect, component_p, permutations = exact_mean_difference(
                clr[:, position],
                observed_hyperoxia,
            )
            clr_rows.append(
                {
                    "Age": age,
                    "subtype": subtype,
                    "clr_difference_hyperoxia_minus_normoxia": effect,
                    "component_exact_two_sided_permutation_p": component_p,
                    "global_squared_clr_distance": observed_distance,
                    "global_exact_permutation_p": global_p,
                    "label_assignments": permutations,
                    "pseudocount": 0.5,
                }
            )
        for row_index, animal_id in enumerate(frame["animal_id"]):
            for position, subtype in enumerate(SUBTYPE_COUNTS):
                clr_donor_rows.append(
                    {
                        "animal_id": animal_id,
                        "Age": age,
                        "Oxygen": frame.loc[row_index, "Oxygen"],
                        "subtype": subtype,
                        "clr_value": float(clr[row_index, position]),
                    }
                )

    logit_results = pd.DataFrame(logit_rows)
    logit_results["BH_q_all_18_logit_tests"] = bh_adjust(
        logit_results["exact_two_sided_permutation_p"]
    )
    clr_results = pd.DataFrame(clr_rows)
    clr_results["BH_q_all_15_CLR_component_tests"] = bh_adjust(
        clr_results["component_exact_two_sided_permutation_p"]
    )
    logit_results.to_csv(
        OUTPUT_DIR / "GSE151974_composition_logit_exact_permutation.tsv",
        sep="\t",
        index=False,
    )
    clr_results.to_csv(
        OUTPUT_DIR / "GSE151974_composition_CLR_exact_permutation.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(clr_donor_rows).to_csv(
        OUTPUT_DIR / "GSE151974_composition_CLR_animal_values.tsv",
        sep="\t",
        index=False,
    )
    print("P14 logit sensitivities")
    print(
        logit_results.loc[logit_results["Age"].eq("P14")].to_string(
            index=False
        )
    )
    print("\nCLR global tests")
    print(
        clr_results[
            [
                "Age",
                "global_squared_clr_distance",
                "global_exact_permutation_p",
            ]
        ]
        .drop_duplicates()
        .to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
