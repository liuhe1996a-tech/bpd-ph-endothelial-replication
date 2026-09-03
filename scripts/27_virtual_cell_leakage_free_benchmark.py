"""Leakage-free animal-grouped benchmark for neonatal endothelial virtual cells.

This script replaces the quarantined 2026-08-07 benchmark.  Its design makes
the outer held-out age completely unavailable during feature selection and
model fitting, selects highly variable genes from outer-training cells only,
uses whole animals for early-stopping validation, repeats both deep models over
ten fixed seeds, and applies the same animal bootstrap draws to every model.
The model matrix contains only fold-specific HVGs.  The overlaps of those HVGs
with the internally selected 198-gene and 33-gene sets are reported as
post-selection sensitivity endpoints and never alter the fitted feature space.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from torch import nn
from torch.utils.data import DataLoader, TensorDataset, WeightedRandomSampler


AGES = ("P3", "P7", "P14")
CELL_TYPES = ("Cap", "Cap-a", "Art", "Vein", "Lymph")
CAPILLARY_TYPES = ("Cap", "Cap-a")
SIMPLE_METHODS = ("identity", "gene_linear", "pca_latent")
DEEP_METHODS = ("vae_latent", "cvae_counterfactual")
ALL_METHODS = SIMPLE_METHODS + DEEP_METHODS
METRICS = (
    "pearson_effect",
    "spearman_effect",
    "rmse_effect",
    "direction_accuracy",
)


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(min(8, max(1, os.cpu_count() or 1)))


class VAE(nn.Module):
    def __init__(self, n_input: int, n_latent: int) -> None:
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_input, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(256, 128),
            nn.LayerNorm(128),
            nn.GELU(),
        )
        self.mu = nn.Linear(128, n_latent)
        self.logvar = nn.Linear(128, n_latent)
        self.decoder = nn.Sequential(
            nn.Linear(n_latent, 128),
            nn.LayerNorm(128),
            nn.GELU(),
            nn.Linear(128, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Linear(256, n_input),
        )

    def encode(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        hidden = self.encoder(x)
        return self.mu(hidden), self.logvar(hidden)

    def forward(
        self, x: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mu, logvar = self.encode(x)
        z = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
        return self.decoder(z), mu, logvar


class ConditionalVAE(nn.Module):
    def __init__(self, n_input: int, n_condition: int, n_latent: int) -> None:
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_input + n_condition, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(256, 128),
            nn.LayerNorm(128),
            nn.GELU(),
        )
        self.mu = nn.Linear(128, n_latent)
        self.logvar = nn.Linear(128, n_latent)
        self.decoder = nn.Sequential(
            nn.Linear(n_latent + n_condition, 128),
            nn.LayerNorm(128),
            nn.GELU(),
            nn.Linear(128, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Linear(256, n_input),
        )

    def encode(
        self, x: torch.Tensor, condition: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        hidden = self.encoder(torch.cat([x, condition], dim=1))
        return self.mu(hidden), self.logvar(hidden)

    def decode(self, z: torch.Tensor, condition: torch.Tensor) -> torch.Tensor:
        return self.decoder(torch.cat([z, condition], dim=1))

    def forward(
        self, x: torch.Tensor, condition: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mu, logvar = self.encode(x, condition)
        z = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
        return self.decode(z, condition), mu, logvar


@dataclass
class SelectionFit:
    model: nn.Module
    best_epoch: int
    history: pd.DataFrame


def make_conditions(metadata: pd.DataFrame) -> np.ndarray:
    age_days = metadata["Age"].map({"P3": 3.0, "P7": 7.0, "P14": 14.0})
    age_scaled = ((age_days.to_numpy(dtype=np.float32) - 7.0) / 5.0)[:, None]
    cell_one_hot = np.column_stack(
        [
            metadata["CellType"].eq(cell_type).to_numpy(dtype=np.float32)
            for cell_type in CELL_TYPES
        ]
    )
    oxygen = metadata["Oxygen"].eq("Hyperoxia").to_numpy(dtype=np.float32)[:, None]
    return np.column_stack([age_scaled, cell_one_hot, oxygen]).astype(np.float32)


def build_model(
    kind: str, n_input: int, n_condition: int, n_latent: int
) -> nn.Module:
    if kind == "vae_latent":
        return VAE(n_input, n_latent)
    if kind == "cvae_counterfactual":
        return ConditionalVAE(n_input, n_condition, n_latent)
    raise ValueError(f"Unsupported model kind: {kind}")


def per_sample_loss(
    model: nn.Module,
    kind: str,
    x: torch.Tensor,
    condition: torch.Tensor | None,
    beta: float,
) -> tuple[torch.Tensor, torch.Tensor]:
    if kind == "vae_latent":
        reconstruction, mu, logvar = model(x)
    else:
        if condition is None:
            raise ValueError("Conditional model requires condition matrix")
        reconstruction, mu, logvar = model(x, condition)
    mse = torch.mean((reconstruction - x) ** 2, dim=1)
    kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)
    return mse + beta * kl, reconstruction


def group_balanced_sampler(groups: np.ndarray, seed: int) -> WeightedRandomSampler:
    unique, counts = np.unique(groups, return_counts=True)
    inverse = {group: 1.0 / count for group, count in zip(unique, counts)}
    weights = torch.as_tensor([inverse[group] for group in groups], dtype=torch.double)
    return WeightedRandomSampler(
        weights,
        num_samples=len(groups),
        replacement=True,
        generator=torch.Generator().manual_seed(seed),
    )


def train_with_grouped_validation(
    kind: str,
    x: np.ndarray,
    condition: np.ndarray,
    groups: np.ndarray,
    train_idx: np.ndarray,
    validation_idx: np.ndarray,
    seed: int,
    n_latent: int,
    max_epochs: int,
    batch_size: int,
    beta: float,
) -> SelectionFit:
    set_seed(seed)
    model = build_model(kind, x.shape[1], condition.shape[1], n_latent)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-5)
    sampler = group_balanced_sampler(groups[train_idx], seed)
    if kind == "vae_latent":
        train_dataset = TensorDataset(torch.from_numpy(x[train_idx]))
    else:
        train_dataset = TensorDataset(
            torch.from_numpy(x[train_idx]), torch.from_numpy(condition[train_idx])
        )
    loader = DataLoader(train_dataset, batch_size=batch_size, sampler=sampler)
    validation_groups = groups[validation_idx]
    history: list[dict[str, float | int]] = []
    best_state: dict[str, torch.Tensor] | None = None
    best_epoch = -1
    best_validation = math.inf
    patience = 8

    for epoch in range(1, max_epochs + 1):
        model.train()
        train_sum = 0.0
        train_n = 0
        for batch in loader:
            batch_x = batch[0]
            batch_condition = batch[1] if kind != "vae_latent" else None
            optimizer.zero_grad(set_to_none=True)
            losses, _ = per_sample_loss(model, kind, batch_x, batch_condition, beta)
            loss = losses.mean()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
            train_sum += float(loss.detach()) * len(batch_x)
            train_n += len(batch_x)

        model.eval()
        with torch.no_grad():
            val_x = torch.from_numpy(x[validation_idx])
            val_condition = (
                None
                if kind == "vae_latent"
                else torch.from_numpy(condition[validation_idx])
            )
            val_losses, _ = per_sample_loss(model, kind, val_x, val_condition, beta)
        val_frame = pd.DataFrame(
            {
                "animal_id": validation_groups,
                "loss": val_losses.detach().cpu().numpy(),
            }
        )
        validation_loss = float(val_frame.groupby("animal_id")["loss"].mean().mean())
        history.append(
            {
                "epoch": epoch,
                "train_loss": train_sum / max(1, train_n),
                "animal_grouped_validation_loss": validation_loss,
            }
        )
        if validation_loss < best_validation - 1e-5:
            best_validation = validation_loss
            best_epoch = epoch
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
        elif epoch - best_epoch >= patience:
            break

    if best_state is None or best_epoch < 1:
        raise RuntimeError(f"{kind} did not produce a finite grouped-validation fit")
    model.load_state_dict(best_state)
    model.eval()
    return SelectionFit(model=model, best_epoch=best_epoch, history=pd.DataFrame(history))


def refit_all_outer_training(
    kind: str,
    x: np.ndarray,
    condition: np.ndarray,
    groups: np.ndarray,
    outer_train_idx: np.ndarray,
    seed: int,
    n_latent: int,
    epochs: int,
    batch_size: int,
    beta: float,
) -> nn.Module:
    set_seed(seed)
    model = build_model(kind, x.shape[1], condition.shape[1], n_latent)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-5)
    sampler = group_balanced_sampler(groups[outer_train_idx], seed)
    if kind == "vae_latent":
        dataset = TensorDataset(torch.from_numpy(x[outer_train_idx]))
    else:
        dataset = TensorDataset(
            torch.from_numpy(x[outer_train_idx]),
            torch.from_numpy(condition[outer_train_idx]),
        )
    loader = DataLoader(dataset, batch_size=batch_size, sampler=sampler)
    for _ in range(epochs):
        model.train()
        for batch in loader:
            batch_x = batch[0]
            batch_condition = batch[1] if kind != "vae_latent" else None
            optimizer.zero_grad(set_to_none=True)
            losses, _ = per_sample_loss(model, kind, batch_x, batch_condition, beta)
            losses.mean().backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
    model.eval()
    return model


def encode_vae(model: VAE, x: np.ndarray, batch_size: int = 512) -> np.ndarray:
    output: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(x), batch_size):
            mu, _ = model.encode(torch.from_numpy(x[start : start + batch_size]))
            output.append(mu.detach().cpu().numpy())
    return np.vstack(output)


def decode_vae(model: VAE, z: np.ndarray, batch_size: int = 512) -> np.ndarray:
    output: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(z), batch_size):
            batch = torch.from_numpy(z[start : start + batch_size].astype(np.float32))
            output.append(model.decoder(batch).detach().cpu().numpy())
    return np.vstack(output)


def encode_cvae(
    model: ConditionalVAE,
    x: np.ndarray,
    condition: np.ndarray,
    batch_size: int = 512,
) -> np.ndarray:
    output: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(x), batch_size):
            mu, _ = model.encode(
                torch.from_numpy(x[start : start + batch_size]),
                torch.from_numpy(condition[start : start + batch_size]),
            )
            output.append(mu.detach().cpu().numpy())
    return np.vstack(output)


def decode_cvae(
    model: ConditionalVAE,
    z: np.ndarray,
    condition: np.ndarray,
    batch_size: int = 512,
) -> np.ndarray:
    output: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(z), batch_size):
            output.append(
                model.decode(
                    torch.from_numpy(z[start : start + batch_size].astype(np.float32)),
                    torch.from_numpy(condition[start : start + batch_size]),
                )
                .detach()
                .cpu()
                .numpy()
            )
    return np.vstack(output)


def normalize_selected_counts(
    counts: sparse.csr_matrix, selected_idx: np.ndarray
) -> np.ndarray:
    library = np.asarray(counts.sum(axis=1)).ravel()
    if np.any(library <= 0):
        raise ValueError("Cells with zero library size detected")
    selected = counts[:, selected_idx].multiply((10000.0 / library)[:, None]).tocsr()
    selected = selected.astype(np.float32)
    selected.data = np.log1p(selected.data)
    return selected.toarray().astype(np.float32)


def select_outer_fold_features(
    counts: sparse.csr_matrix,
    genes: np.ndarray,
    outer_train_mask: np.ndarray,
    replicated_genes: set[str],
    signature_genes: set[str],
    n_hvg: int,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    train_counts = counts[outer_train_mask]
    library = np.asarray(train_counts.sum(axis=1)).ravel()
    norm = train_counts.multiply((10000.0 / library)[:, None]).tocsr().astype(np.float32)
    norm.data = np.log1p(norm.data)
    detected = np.asarray((train_counts > 0).sum(axis=0)).ravel()
    mean = np.asarray(norm.mean(axis=0)).ravel()
    second = np.asarray(norm.power(2).mean(axis=0)).ravel()
    variance = np.maximum(second - mean**2, 0)
    eligible = detected >= max(20, int(round(0.005 * train_counts.shape[0])))
    dispersion = np.full(len(genes), -np.inf, dtype=float)
    dispersion[eligible] = variance[eligible] / (mean[eligible] + 1e-6)
    eligible_idx = np.flatnonzero(np.isfinite(dispersion))
    if len(eligible_idx) < n_hvg:
        raise ValueError(f"Only {len(eligible_idx)} genes passed fold detection filter")
    hvg_idx = eligible_idx[np.argsort(dispersion[eligible_idx])[-n_hvg:]]
    # The fitted feature matrix is defined solely from outer-training cells.
    # Disease-informed genes are annotations on the selected HVGs, not forced
    # additions to the model space.
    selected_idx = np.sort(hvg_idx)
    manifest = pd.DataFrame(
        {
            "gene_index": selected_idx,
            "gene": genes[selected_idx],
            "primary_fold_hvg": np.ones(len(selected_idx), dtype=bool),
            "secondary_replicated_198_hvg_intersection": [
                gene in replicated_genes for gene in genes[selected_idx]
            ],
            "secondary_signature_33_hvg_intersection": [
                gene in signature_genes for gene in genes[selected_idx]
            ],
            "outer_training_detected_cells": detected[selected_idx],
            "outer_training_dispersion": dispersion[selected_idx],
        }
    )
    return selected_idx, hvg_idx, manifest


def stratified_validation_animals(
    metadata: pd.DataFrame, outer_train_mask: np.ndarray, seed: int
) -> tuple[set[str], pd.DataFrame]:
    animal_table = (
        metadata.loc[outer_train_mask, ["animal_id", "Age", "Oxygen"]]
        .drop_duplicates()
        .sort_values(["Age", "Oxygen", "animal_id"])
    )
    if animal_table["animal_id"].duplicated().any():
        raise ValueError("Animal appears in more than one age/oxygen stratum")
    rng = np.random.default_rng(seed)
    validation: set[str] = set()
    rows: list[dict[str, object]] = []
    for (age, oxygen), group in animal_table.groupby(["Age", "Oxygen"], sort=True):
        animals = group["animal_id"].astype(str).to_numpy()
        chosen = str(rng.choice(animals, size=1, replace=False)[0])
        validation.add(chosen)
        for animal in animals:
            rows.append(
                {
                    "animal_id": str(animal),
                    "age": age,
                    "oxygen": oxygen,
                    "role": "validation" if str(animal) == chosen else "fit_training",
                }
            )
    return validation, pd.DataFrame(rows)


def animal_means(
    values: np.ndarray, metadata: pd.DataFrame
) -> tuple[np.ndarray, list[str]]:
    animals = sorted(metadata["animal_id"].astype(str).unique())
    animal_vector = metadata["animal_id"].astype(str).to_numpy()
    means = np.vstack([values[animal_vector == animal].mean(axis=0) for animal in animals])
    return means.astype(np.float32), animals


def balanced_effect(
    values: np.ndarray,
    metadata: pd.DataFrame,
    cell_type: str,
    training_ages: tuple[str, ...],
) -> np.ndarray:
    effects: list[np.ndarray] = []
    for age in training_ages:
        subset = metadata["Age"].eq(age) & metadata["CellType"].eq(cell_type)
        age_metadata = metadata.loc[subset].reset_index(drop=True)
        age_values = values[subset.to_numpy()]
        condition_means: dict[str, np.ndarray] = {}
        for oxygen in ("Normoxia", "Hyperoxia"):
            condition = age_metadata["Oxygen"].eq(oxygen)
            pseudobulk, _ = animal_means(
                age_values[condition.to_numpy()], age_metadata.loc[condition]
            )
            condition_means[oxygen] = pseudobulk.mean(axis=0)
        effects.append(condition_means["Hyperoxia"] - condition_means["Normoxia"])
    return np.vstack(effects).mean(axis=0)


def safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    if method == "pearson":
        return float(pearsonr(x, y).statistic)
    return float(spearmanr(x, y).statistic)


def metric_row(
    prediction: np.ndarray,
    observed: np.ndarray,
    normoxia: np.ndarray,
    gene_indices: np.ndarray,
) -> dict[str, float | int]:
    predicted_effect = (
        prediction.mean(axis=0)[gene_indices] - normoxia.mean(axis=0)[gene_indices]
    )
    observed_effect = (
        observed.mean(axis=0)[gene_indices] - normoxia.mean(axis=0)[gene_indices]
    )
    return {
        "n_genes": int(len(gene_indices)),
        "pearson_effect": safe_corr(predicted_effect, observed_effect, "pearson"),
        "spearman_effect": safe_corr(predicted_effect, observed_effect, "spearman"),
        "rmse_effect": float(np.sqrt(np.mean((predicted_effect - observed_effect) ** 2))),
        "direction_accuracy": float(
            np.mean(np.sign(predicted_effect) == np.sign(observed_effect))
        ),
    }


def endpoint_indices(manifest: pd.DataFrame) -> dict[str, np.ndarray]:
    return {
        "primary_fold_hvg": np.flatnonzero(manifest["primary_fold_hvg"].to_numpy()),
        "secondary_replicated_198_hvg_intersection": np.flatnonzero(
            manifest["secondary_replicated_198_hvg_intersection"].to_numpy()
        ),
        "secondary_signature_33_hvg_intersection": np.flatnonzero(
            manifest["secondary_signature_33_hvg_intersection"].to_numpy()
        ),
    }


def summarize_comparisons(bootstrap: pd.DataFrame) -> pd.DataFrame:
    capillary = bootstrap[bootstrap["cell_type"].isin(CAPILLARY_TYPES)].copy()
    rows: list[dict[str, object]] = []
    comparisons = [
        ("vae_latent", "gene_linear"),
        ("vae_latent", "pca_latent"),
        ("cvae_counterfactual", "gene_linear"),
        ("cvae_counterfactual", "pca_latent"),
    ]
    for endpoint, endpoint_data in capillary.groupby("endpoint", sort=False):
        baseline = (
            endpoint_data[endpoint_data["seed"].eq(-1)]
            .groupby(["bootstrap", "method"], observed=True)[list(METRICS)]
            .mean()
            .reset_index()
        )
        # Average across the fixed training seeds within each biological
        # bootstrap replicate.  The percentile interval therefore reflects
        # paired animal resampling, while seed-to-seed variation is reported in
        # the separate training-variance table.
        deep = (
            endpoint_data[endpoint_data["seed"].ge(0)]
            .groupby(["bootstrap", "method"], observed=True)[list(METRICS)]
            .mean()
            .reset_index()
        )
        for deep_method, baseline_method in comparisons:
            if deep_method not in set(deep["method"]):
                continue
            for metric in METRICS:
                deep_values = deep.loc[
                    deep["method"].eq(deep_method),
                    ["bootstrap", metric],
                ]
                base_values = baseline.loc[
                    baseline["method"].eq(baseline_method),
                    ["bootstrap", metric],
                ]
                paired = deep_values.merge(
                    base_values,
                    on="bootstrap",
                    how="inner",
                    suffixes=("_deep", "_baseline"),
                    validate="one_to_one",
                )
                if metric == "rmse_effect":
                    delta = (
                        paired[f"{metric}_baseline"]
                        - paired[f"{metric}_deep"]
                    ).to_numpy()
                else:
                    delta = (
                        paired[f"{metric}_deep"]
                        - paired[f"{metric}_baseline"]
                    ).to_numpy()
                rows.append(
                    {
                        "endpoint": endpoint,
                        "deep_model": deep_method,
                        "baseline": baseline_method,
                        "metric": metric,
                        "median_delta": float(np.nanmedian(delta)),
                        "lower_2_5": float(np.nanquantile(delta, 0.025)),
                        "upper_97_5": float(np.nanquantile(delta, 0.975)),
                        "win_fraction": float(np.nanmean(delta > 0)),
                        "n_paired_bootstrap_values": int(np.isfinite(delta).sum()),
                    }
                )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--heldout-ages", nargs="+", default=list(AGES))
    parser.add_argument("--n-hvg", type=int, default=1800)
    parser.add_argument("--latent-dim", type=int, default=32)
    parser.add_argument("--max-epochs", type=int, default=60)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--beta", type=float, default=1e-4)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--seeds", nargs="+", type=int, default=list(range(20260817, 20260827)))
    parser.add_argument("--base-seed", type=int, default=20260817)
    parser.add_argument(
        "--smoke-test",
        action="store_true",
        help="Permit fewer than ten seeds for a non-reportable pipeline check.",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    counts = sparse.load_npz(args.matrix).tocsr()
    metadata = pd.read_csv(args.metadata, sep="\t")
    genes_frame = pd.read_csv(args.genes, sep="\t")
    genes = genes_frame["gene"].astype(str).to_numpy()
    signature = pd.read_csv(args.signature, sep="\t")
    if counts.shape != (len(metadata), len(genes)):
        raise ValueError("Matrix, metadata and gene dimensions do not match")
    if len(args.seeds) < 10 and not args.smoke_test:
        raise ValueError("At least ten deep-model seeds are required")
    replicated_genes = set(signature["gene"].astype(str).unique())
    signature_genes = set(
        signature.loc[
            signature["concordant_endothelial_subtypes"].ge(3), "gene"
        ].astype(str)
    )
    if len(replicated_genes) != 198 or len(signature_genes) != 33:
        raise ValueError(
            f"Expected locked 198/33 sets, observed {len(replicated_genes)}/{len(signature_genes)}"
        )

    point_rows: list[dict[str, object]] = []
    bootstrap_rows: list[dict[str, object]] = []
    prediction_rows: list[pd.DataFrame] = []
    feature_manifests: list[pd.DataFrame] = []
    split_rows: list[pd.DataFrame] = []
    history_rows: list[pd.DataFrame] = []
    bootstrap_draw_rows: list[dict[str, object]] = []
    array_payload: dict[str, np.ndarray] = {}
    fold_audit: dict[str, object] = {}

    for fold_number, heldout_age in enumerate(args.heldout_ages, start=1):
        if heldout_age not in AGES:
            raise ValueError(f"Unsupported held-out age: {heldout_age}")
        outer_train_mask = metadata["Age"].ne(heldout_age).to_numpy()
        outer_test_mask = ~outer_train_mask
        selected_idx, _, manifest = select_outer_fold_features(
            counts,
            genes,
            outer_train_mask,
            replicated_genes,
            signature_genes,
            args.n_hvg,
        )
        manifest.insert(0, "heldout_age", heldout_age)
        manifest.insert(1, "model_feature_position", np.arange(len(manifest)))
        feature_manifests.append(manifest)
        selected_genes = manifest["gene"].astype(str).to_numpy()
        endpoints = endpoint_indices(manifest)
        x = normalize_selected_counts(counts, selected_idx)
        train_mean = x[outer_train_mask].mean(axis=0)
        train_sd = x[outer_train_mask].std(axis=0)
        train_sd[train_sd < 1e-4] = 1.0
        x_scaled = ((x - train_mean) / train_sd).astype(np.float32)
        conditions = make_conditions(metadata)
        groups = metadata["animal_id"].astype(str).to_numpy()
        outer_train_idx = np.flatnonzero(outer_train_mask)
        training_ages = tuple(age for age in AGES if age != heldout_age)
        pca_components = min(args.latent_dim, x.shape[1] - 1, len(outer_train_idx) - 1)
        pca = PCA(
            n_components=pca_components,
            svd_solver="randomized",
            random_state=args.base_seed + fold_number,
        )
        pca.fit(x_scaled[outer_train_idx])
        pca_latent = pca.transform(x_scaled).astype(np.float32)

        simple_predictions: dict[str, dict[str, np.ndarray]] = {}
        task_observed: dict[str, np.ndarray] = {}
        task_normoxia: dict[str, np.ndarray] = {}
        task_norm_animals: dict[str, list[str]] = {}
        task_obs_animals: dict[str, list[str]] = {}
        for cell_type in CELL_TYPES:
            norm_mask = (
                metadata["Age"].eq(heldout_age)
                & metadata["Oxygen"].eq("Normoxia")
                & metadata["CellType"].eq(cell_type)
            )
            hyp_mask = (
                metadata["Age"].eq(heldout_age)
                & metadata["Oxygen"].eq("Hyperoxia")
                & metadata["CellType"].eq(cell_type)
            )
            if norm_mask.sum() < 5 or hyp_mask.sum() < 5:
                continue
            norm_meta = metadata.loc[norm_mask].reset_index(drop=True)
            hyp_meta = metadata.loc[hyp_mask].reset_index(drop=True)
            norm_pb, norm_animals = animal_means(x[norm_mask.to_numpy()], norm_meta)
            obs_pb, obs_animals = animal_means(x[hyp_mask.to_numpy()], hyp_meta)
            gene_delta = balanced_effect(x, metadata, cell_type, training_ages)
            pca_delta = balanced_effect(pca_latent, metadata, cell_type, training_ages)
            pred_pca_scaled = pca.inverse_transform(
                pca_latent[norm_mask.to_numpy()] + pca_delta
            )
            pred_pca_cells = np.maximum(pred_pca_scaled * train_sd + train_mean, 0)
            pred_pca, pca_animals = animal_means(pred_pca_cells, norm_meta)
            if norm_animals != pca_animals:
                raise ValueError("PCA prediction and normoxia animal order differ")
            task = f"{heldout_age}__{cell_type}"
            task_observed[task] = obs_pb
            task_normoxia[task] = norm_pb
            task_norm_animals[task] = norm_animals
            task_obs_animals[task] = obs_animals
            simple_predictions[task] = {
                "identity": norm_pb.copy(),
                "gene_linear": np.maximum(norm_pb + gene_delta, 0),
                "pca_latent": pred_pca,
            }

        deep_predictions: dict[str, dict[int, dict[str, np.ndarray]]] = {
            task: {seed: {} for seed in args.seeds} for task in simple_predictions
        }
        for seed in args.seeds:
            validation_animals, split = stratified_validation_animals(
                metadata, outer_train_mask, seed + fold_number * 1000
            )
            split.insert(0, "heldout_age", heldout_age)
            split.insert(1, "seed", seed)
            split_rows.append(split)
            fit_train_mask = outer_train_mask & ~metadata["animal_id"].astype(str).isin(
                validation_animals
            ).to_numpy()
            validation_mask = outer_train_mask & metadata["animal_id"].astype(str).isin(
                validation_animals
            ).to_numpy()
            fit_train_idx = np.flatnonzero(fit_train_mask)
            validation_idx = np.flatnonzero(validation_mask)
            if set(groups[fit_train_idx]) & set(groups[validation_idx]):
                raise RuntimeError("Animal leakage between fit and validation cells")
            trained: dict[str, nn.Module] = {}
            for model_offset, kind in enumerate(DEEP_METHODS):
                model_seed = seed + fold_number * 10000 + model_offset * 100000
                selection = train_with_grouped_validation(
                    kind,
                    x_scaled,
                    conditions,
                    groups,
                    fit_train_idx,
                    validation_idx,
                    model_seed,
                    args.latent_dim,
                    args.max_epochs,
                    args.batch_size,
                    args.beta,
                )
                history = selection.history.copy()
                history.insert(0, "heldout_age", heldout_age)
                history.insert(1, "seed", seed)
                history.insert(2, "method", kind)
                history["selected_best_epoch"] = selection.best_epoch
                history_rows.append(history)
                trained[kind] = refit_all_outer_training(
                    kind,
                    x_scaled,
                    conditions,
                    groups,
                    outer_train_idx,
                    model_seed,
                    args.latent_dim,
                    selection.best_epoch,
                    args.batch_size,
                    args.beta,
                )
                print(
                    f"fold={heldout_age} seed={seed} method={kind} "
                    f"best_epoch={selection.best_epoch}",
                    flush=True,
                )

            vae_latent = encode_vae(trained["vae_latent"], x_scaled)
            cvae_latent = encode_cvae(
                trained["cvae_counterfactual"], x_scaled, conditions
            )
            for cell_type in CELL_TYPES:
                task = f"{heldout_age}__{cell_type}"
                if task not in simple_predictions:
                    continue
                norm_mask = (
                    metadata["Age"].eq(heldout_age)
                    & metadata["Oxygen"].eq("Normoxia")
                    & metadata["CellType"].eq(cell_type)
                )
                norm_meta = metadata.loc[norm_mask].reset_index(drop=True)
                vae_delta = balanced_effect(
                    vae_latent, metadata, cell_type, training_ages
                )
                vae_pred_scaled = decode_vae(
                    trained["vae_latent"],
                    vae_latent[norm_mask.to_numpy()] + vae_delta,
                )
                vae_pred_cells = np.maximum(
                    vae_pred_scaled * train_sd + train_mean, 0
                )
                vae_pred, vae_animals = animal_means(vae_pred_cells, norm_meta)
                norm_condition = conditions[norm_mask.to_numpy()].copy()
                hyper_condition = norm_condition.copy()
                hyper_condition[:, -1] = 1.0
                cvae_pred_scaled = decode_cvae(
                    trained["cvae_counterfactual"],
                    cvae_latent[norm_mask.to_numpy()],
                    hyper_condition,
                )
                cvae_pred_cells = np.maximum(
                    cvae_pred_scaled * train_sd + train_mean, 0
                )
                cvae_pred, cvae_animals = animal_means(cvae_pred_cells, norm_meta)
                if not (
                    task_norm_animals[task] == vae_animals == cvae_animals
                ):
                    raise ValueError("Deep prediction and normoxia animal order differ")
                deep_predictions[task][seed] = {
                    "vae_latent": vae_pred,
                    "cvae_counterfactual": cvae_pred,
                }

        capillary_tasks = [
            f"{heldout_age}__{cell_type}"
            for cell_type in CAPILLARY_TYPES
            if f"{heldout_age}__{cell_type}" in simple_predictions
        ]
        common_capillary_normoxia = sorted(
            set.intersection(*(set(task_norm_animals[task]) for task in capillary_tasks))
        )
        common_capillary_hyperoxia = sorted(
            set.intersection(*(set(task_obs_animals[task]) for task in capillary_tasks))
        )
        if not common_capillary_normoxia or not common_capillary_hyperoxia:
            raise RuntimeError(f"No common capillary animals in held-out fold {heldout_age}")
        capillary_rng = np.random.default_rng(args.base_seed + fold_number * 1000)
        capillary_draws = [
            (
                capillary_rng.choice(
                    common_capillary_normoxia,
                    size=len(common_capillary_normoxia),
                    replace=True,
                ),
                capillary_rng.choice(
                    common_capillary_hyperoxia,
                    size=len(common_capillary_hyperoxia),
                    replace=True,
                ),
            )
            for _ in range(args.bootstrap)
        ]
        for replicate, (norm_ids, observed_ids) in enumerate(capillary_draws):
            bootstrap_draw_rows.extend(
                {
                    "heldout_age": heldout_age,
                    "bootstrap": replicate,
                    "oxygen_group": "Normoxia",
                    "draw_position": position,
                    "animal_id": str(animal),
                }
                for position, animal in enumerate(norm_ids)
            )
            bootstrap_draw_rows.extend(
                {
                    "heldout_age": heldout_age,
                    "bootstrap": replicate,
                    "oxygen_group": "Hyperoxia",
                    "draw_position": position,
                    "animal_id": str(animal),
                }
                for position, animal in enumerate(observed_ids)
            )

        for task, method_predictions in simple_predictions.items():
            _, cell_type = task.split("__", 1)
            observed = task_observed[task]
            normoxia = task_normoxia[task]
            for method, prediction in method_predictions.items():
                for endpoint, indices in endpoints.items():
                    row = metric_row(prediction, observed, normoxia, indices)
                    row.update(
                        {
                            "heldout_age": heldout_age,
                            "cell_type": cell_type,
                            "method": method,
                            "seed": -1,
                            "endpoint": endpoint,
                            "n_normoxia_animals": len(normoxia),
                            "n_hyperoxia_animals": len(observed),
                        }
                    )
                    point_rows.append(row)
            for seed in args.seeds:
                for method, prediction in deep_predictions[task][seed].items():
                    for endpoint, indices in endpoints.items():
                        row = metric_row(prediction, observed, normoxia, indices)
                        row.update(
                            {
                                "heldout_age": heldout_age,
                                "cell_type": cell_type,
                                "method": method,
                                "seed": seed,
                                "endpoint": endpoint,
                                "n_normoxia_animals": len(normoxia),
                                "n_hyperoxia_animals": len(observed),
                            }
                        )
                        point_rows.append(row)

            if cell_type in CAPILLARY_TYPES:
                key = task.replace("-", "a").replace("__", "_")
                array_payload[f"{key}_observed"] = observed
                array_payload[f"{key}_normoxia"] = normoxia
                for method, prediction in method_predictions.items():
                    array_payload[f"{key}_{method}"] = prediction
                for seed in args.seeds:
                    for method, prediction in deep_predictions[task][seed].items():
                        array_payload[f"{key}_{method}_seed{seed}"] = prediction

            norm_mean = normoxia.mean(axis=0)
            obs_mean = observed.mean(axis=0)
            for method, prediction in method_predictions.items():
                prediction_rows.append(
                    pd.DataFrame(
                        {
                            "heldout_age": heldout_age,
                            "cell_type": cell_type,
                            "method": method,
                            "seed": -1,
                            "gene": selected_genes,
                            "normoxia_mean": norm_mean,
                            "observed_hyperoxia_mean": obs_mean,
                            "predicted_hyperoxia_mean": prediction.mean(axis=0),
                            "primary_fold_hvg": manifest["primary_fold_hvg"].to_numpy(),
                            "secondary_replicated_198_hvg_intersection": manifest[
                                "secondary_replicated_198_hvg_intersection"
                            ].to_numpy(),
                            "secondary_signature_33_hvg_intersection": manifest[
                                "secondary_signature_33_hvg_intersection"
                            ].to_numpy(),
                        }
                    )
                )
            for seed in args.seeds:
                for method, prediction in deep_predictions[task][seed].items():
                    prediction_rows.append(
                        pd.DataFrame(
                            {
                                "heldout_age": heldout_age,
                                "cell_type": cell_type,
                                "method": method,
                                "seed": seed,
                                "gene": selected_genes,
                                "normoxia_mean": norm_mean,
                                "observed_hyperoxia_mean": obs_mean,
                                "predicted_hyperoxia_mean": prediction.mean(axis=0),
                                "primary_fold_hvg": manifest[
                                    "primary_fold_hvg"
                                ].to_numpy(),
                                "secondary_replicated_198_hvg_intersection": manifest[
                                    "secondary_replicated_198_hvg_intersection"
                                ].to_numpy(),
                                "secondary_signature_33_hvg_intersection": manifest[
                                    "secondary_signature_33_hvg_intersection"
                                ].to_numpy(),
                            }
                        )
                    )

            rng = np.random.default_rng(
                args.base_seed + fold_number * 1000 + CELL_TYPES.index(cell_type) * 100
            )
            for replicate in range(args.bootstrap):
                if cell_type in CAPILLARY_TYPES:
                    norm_ids, observed_ids = capillary_draws[replicate]
                    norm_lookup = {animal: index for index, animal in enumerate(task_norm_animals[task])}
                    observed_lookup = {animal: index for index, animal in enumerate(task_obs_animals[task])}
                    norm_idx = np.asarray([norm_lookup[str(animal)] for animal in norm_ids], dtype=int)
                    observed_idx = np.asarray([observed_lookup[str(animal)] for animal in observed_ids], dtype=int)
                else:
                    norm_idx = rng.integers(0, len(normoxia), len(normoxia))
                    observed_idx = rng.integers(0, len(observed), len(observed))
                shared: list[tuple[str, int, np.ndarray]] = [
                    (method, -1, prediction)
                    for method, prediction in method_predictions.items()
                ]
                shared.extend(
                    (method, seed, prediction)
                    for seed in args.seeds
                    for method, prediction in deep_predictions[task][seed].items()
                )
                for method, seed, prediction in shared:
                    for endpoint, indices in endpoints.items():
                        row = metric_row(
                            prediction[norm_idx],
                            observed[observed_idx],
                            normoxia[norm_idx],
                            indices,
                        )
                        row.update(
                            {
                                "heldout_age": heldout_age,
                                "cell_type": cell_type,
                                "method": method,
                                "seed": seed,
                                "endpoint": endpoint,
                                "bootstrap": replicate,
                            }
                        )
                        bootstrap_rows.append(row)

        fold_audit[heldout_age] = {
            "outer_training_cells": int(outer_train_mask.sum()),
            "outer_test_cells": int(outer_test_mask.sum()),
            "outer_training_ages": list(training_ages),
            "heldout_age_completely_excluded_from_feature_selection_and_fit": True,
            "primary_hvg_count": int(manifest["primary_fold_hvg"].sum()),
            "model_feature_count": int(len(manifest)),
            "replicated_198_hvg_intersection": int(manifest["secondary_replicated_198_hvg_intersection"].sum()),
            "signature_33_hvg_intersection": int(manifest["secondary_signature_33_hvg_intersection"].sum()),
            "common_capillary_normoxia_animals": len(common_capillary_normoxia),
            "common_capillary_hyperoxia_animals": len(common_capillary_hyperoxia),
            "feature_sha256": hashlib.sha256(
                "\n".join(selected_genes).encode("utf-8")
            ).hexdigest(),
        }

    point = pd.DataFrame(point_rows)
    bootstrap = pd.DataFrame(bootstrap_rows)
    predictions = pd.concat(prediction_rows, ignore_index=True)
    features = pd.concat(feature_manifests, ignore_index=True)
    splits = pd.concat(split_rows, ignore_index=True)
    histories = pd.concat(history_rows, ignore_index=True)
    comparisons = summarize_comparisons(bootstrap)
    training_variance = (
        point[
            point["cell_type"].isin(CAPILLARY_TYPES)
            & point["method"].isin(DEEP_METHODS)
        ]
        .groupby(["endpoint", "method", "seed"], observed=True)[list(METRICS)]
        .mean()
        .reset_index()
        .groupby(["endpoint", "method"], observed=True)[list(METRICS)]
        .agg(["mean", "std", "min", "max"])
    )
    training_variance.columns = ["__".join(column) for column in training_variance.columns]
    training_variance = training_variance.reset_index()

    point.to_csv(args.output_dir / "virtual_cell_point_metrics.tsv", sep="\t", index=False)
    bootstrap.to_csv(
        args.output_dir / "virtual_cell_strictly_paired_bootstrap.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(bootstrap_draw_rows).to_csv(
        args.output_dir / "virtual_cell_joint_capillary_bootstrap_draws.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    predictions.to_csv(
        args.output_dir / "virtual_cell_gene_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    features.to_csv(
        args.output_dir / "virtual_cell_fold_feature_manifest.tsv", sep="\t", index=False
    )
    splits.to_csv(
        args.output_dir / "virtual_cell_grouped_validation_splits.tsv", sep="\t", index=False
    )
    histories.to_csv(
        args.output_dir / "virtual_cell_training_histories.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    comparisons.to_csv(
        args.output_dir / "virtual_cell_paired_model_comparisons.tsv",
        sep="\t",
        index=False,
    )
    training_variance.to_csv(
        args.output_dir / "virtual_cell_deep_model_training_variance.tsv",
        sep="\t",
        index=False,
    )
    np.savez_compressed(
        args.output_dir / "virtual_cell_capillary_animal_level_arrays.npz",
        **array_payload,
    )
    audit = {
        "analysis_version": "2026-08-18 R8 HVG-only benchmark",
        "input_cells": int(len(metadata)),
        "input_animals": int(metadata["animal_id"].nunique()),
        "outer_folds": args.heldout_ages,
        "deep_model_seeds": args.seeds,
        "n_deep_model_seeds": len(args.seeds),
        "reportable_run": bool(not args.smoke_test and len(args.seeds) >= 10),
        "bootstrap_replicates": args.bootstrap,
        "bootstrap_pairing": (
            "one animal-ID draw per held-out age and oxygen group is shared by Cap and Cap-a "
            "and by every method; predicted profiles remain paired to starting normoxia animals"
        ),
        "comparison_interval_policy": (
            "deep-model metrics are averaged across training seeds within each biological "
            "bootstrap replicate; training-seed variance is reported separately"
        ),
        "model_feature_space": "exactly 1,800 fold-specific HVGs selected without outer-test cells; no disease-informed genes are forced into the matrix",
        "primary_endpoint": "the complete fold-specific HVG model space",
        "secondary_endpoints": [
            "intersection of the fold-specific HVGs with the internally selected 198-gene set",
            "intersection of the fold-specific HVGs with the internally selected 33-gene set"
        ],
        "validation_split": "whole animals, one animal per outer-training age-by-oxygen stratum",
        "model_selection": (
            "fixed architecture and hyperparameters; grouped validation selected epoch only; "
            "final models refit on all outer-training animals for the selected epoch"
        ),
        "folds": fold_audit,
        "output_rows": {
            "point_metrics": int(len(point)),
            "bootstrap_metrics": int(len(bootstrap)),
            "gene_predictions": int(len(predictions)),
        },
    }
    (args.output_dir / "virtual_cell_leakage_free_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(
        point[
            point["cell_type"].isin(CAPILLARY_TYPES)
            & point["endpoint"].eq("primary_fold_hvg")
        ]
        .groupby(["method", "seed"], observed=True)["spearman_effect"]
        .mean()
        .to_string(),
        flush=True,
    )
    print(comparisons.to_string(index=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
