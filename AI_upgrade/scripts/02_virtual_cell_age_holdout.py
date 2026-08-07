"""Benchmark a scGen-like virtual-cell model with age-held-out validation.

The target age-by-hyperoxia cells are excluded from model fitting.  Hyperoxia
vectors are learned from the other ages using animal-balanced latent means and
applied to normoxic cells from the held-out age.  Evaluation is performed on
animal-level pseudobulks and includes identity, gene-space linear-shift and PCA
latent-shift baselines.  This is intentionally a benchmark, not a claim that a
deep model is intrinsically superior.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from torch import nn
from torch.utils.data import DataLoader, TensorDataset, WeightedRandomSampler


AGES = ("P3", "P7", "P14")
CELL_TYPES = ("Cap", "Cap-a", "Art", "Vein", "Lymph")
METHODS = (
    "identity",
    "gene_linear",
    "pca_latent",
    "vae_latent",
    "cvae_counterfactual",
)


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(min(12, max(1, torch.get_num_threads())))


class VAE(nn.Module):
    def __init__(self, n_input: int, n_latent: int = 32) -> None:
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
        h = self.encoder(x)
        return self.mu(h), self.logvar(h)

    def forward(
        self, x: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mu, logvar = self.encode(x)
        std = torch.exp(0.5 * logvar)
        z = mu + torch.randn_like(std) * std
        return self.decoder(z), mu, logvar


class ConditionalVAE(nn.Module):
    """Conditional VAE for direct normoxia-to-hyperoxia counterfactuals."""

    def __init__(self, n_input: int, n_condition: int, n_latent: int = 32) -> None:
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
        h = self.encoder(torch.cat([x, condition], dim=1))
        return self.mu(h), self.logvar(h)

    def decode(self, z: torch.Tensor, condition: torch.Tensor) -> torch.Tensor:
        return self.decoder(torch.cat([z, condition], dim=1))

    def forward(
        self, x: torch.Tensor, condition: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mu, logvar = self.encode(x, condition)
        std = torch.exp(0.5 * logvar)
        z = mu + torch.randn_like(std) * std
        return self.decode(z, condition), mu, logvar


@dataclass
class FitResult:
    model: VAE
    history: pd.DataFrame
    best_epoch: int


@dataclass
class ConditionalFitResult:
    model: ConditionalVAE
    history: pd.DataFrame
    best_epoch: int


def train_vae(
    x: np.ndarray,
    groups: np.ndarray,
    seed: int,
    max_epochs: int,
    batch_size: int,
    beta: float,
) -> FitResult:
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(x))
    n_val = max(batch_size, int(round(0.1 * len(x))))
    val_idx = order[:n_val]
    train_idx = order[n_val:]
    train_groups = groups[train_idx]
    unique, counts = np.unique(train_groups, return_counts=True)
    inverse = {group: 1.0 / count for group, count in zip(unique, counts)}
    weights = torch.as_tensor(
        [inverse[group] for group in train_groups], dtype=torch.double
    )
    sampler = WeightedRandomSampler(
        weights, num_samples=len(train_idx), replacement=True,
        generator=torch.Generator().manual_seed(seed),
    )
    train_tensor = torch.from_numpy(x[train_idx])
    val_tensor = torch.from_numpy(x[val_idx])
    loader = DataLoader(
        TensorDataset(train_tensor), batch_size=batch_size, sampler=sampler
    )
    val_loader = DataLoader(
        TensorDataset(val_tensor), batch_size=batch_size, shuffle=False
    )
    model = VAE(x.shape[1])
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-5)
    history: list[dict[str, float | int]] = []
    best_state: dict[str, torch.Tensor] | None = None
    best_epoch = -1
    best_val = math.inf
    patience = 10

    for epoch in range(1, max_epochs + 1):
        model.train()
        train_loss = 0.0
        train_n = 0
        for (batch,) in loader:
            optimizer.zero_grad(set_to_none=True)
            recon, mu, logvar = model(batch)
            mse = torch.mean((recon - batch) ** 2)
            kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
            loss = mse + beta * kl
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
            train_loss += float(loss.detach()) * len(batch)
            train_n += len(batch)
        model.eval()
        val_loss = 0.0
        val_n = 0
        with torch.no_grad():
            for (batch,) in val_loader:
                recon, mu, logvar = model(batch)
                mse = torch.mean((recon - batch) ** 2)
                kl = -0.5 * torch.mean(
                    1 + logvar - mu.pow(2) - logvar.exp()
                )
                loss = mse + beta * kl
                val_loss += float(loss) * len(batch)
                val_n += len(batch)
        row = {
            "epoch": epoch,
            "train_loss": train_loss / train_n,
            "validation_loss": val_loss / val_n,
        }
        history.append(row)
        if row["validation_loss"] < best_val - 1e-5:
            best_val = float(row["validation_loss"])
            best_epoch = epoch
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
        elif epoch - best_epoch >= patience:
            break
        if epoch == 1 or epoch % 10 == 0:
            print(
                f"epoch={epoch} train={row['train_loss']:.5f} "
                f"val={row['validation_loss']:.5f}",
                flush=True,
            )
    if best_state is None:
        raise RuntimeError("VAE training did not produce a finite model.")
    model.load_state_dict(best_state)
    model.eval()
    return FitResult(model, pd.DataFrame(history), best_epoch)


def train_cvae(
    x: np.ndarray,
    condition: np.ndarray,
    groups: np.ndarray,
    seed: int,
    max_epochs: int,
    batch_size: int,
    beta: float,
) -> ConditionalFitResult:
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(x))
    n_val = max(batch_size, int(round(0.1 * len(x))))
    val_idx = order[:n_val]
    train_idx = order[n_val:]
    train_groups = groups[train_idx]
    unique, counts = np.unique(train_groups, return_counts=True)
    inverse = {group: 1.0 / count for group, count in zip(unique, counts)}
    weights = torch.as_tensor(
        [inverse[group] for group in train_groups], dtype=torch.double
    )
    sampler = WeightedRandomSampler(
        weights, num_samples=len(train_idx), replacement=True,
        generator=torch.Generator().manual_seed(seed),
    )
    train_x = torch.from_numpy(x[train_idx])
    train_c = torch.from_numpy(condition[train_idx])
    val_x = torch.from_numpy(x[val_idx])
    val_c = torch.from_numpy(condition[val_idx])
    loader = DataLoader(
        TensorDataset(train_x, train_c), batch_size=batch_size, sampler=sampler
    )
    val_loader = DataLoader(
        TensorDataset(val_x, val_c), batch_size=batch_size, shuffle=False
    )
    model = ConditionalVAE(x.shape[1], condition.shape[1])
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-5)
    history: list[dict[str, float | int]] = []
    best_state: dict[str, torch.Tensor] | None = None
    best_epoch = -1
    best_val = math.inf
    patience = 10
    for epoch in range(1, max_epochs + 1):
        model.train()
        train_loss = 0.0
        train_n = 0
        for batch, cond in loader:
            optimizer.zero_grad(set_to_none=True)
            recon, mu, logvar = model(batch, cond)
            mse = torch.mean((recon - batch) ** 2)
            kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
            loss = mse + beta * kl
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
            train_loss += float(loss.detach()) * len(batch)
            train_n += len(batch)
        model.eval()
        val_loss = 0.0
        val_n = 0
        with torch.no_grad():
            for batch, cond in val_loader:
                recon, mu, logvar = model(batch, cond)
                mse = torch.mean((recon - batch) ** 2)
                kl = -0.5 * torch.mean(
                    1 + logvar - mu.pow(2) - logvar.exp()
                )
                loss = mse + beta * kl
                val_loss += float(loss) * len(batch)
                val_n += len(batch)
        row = {
            "epoch": epoch,
            "train_loss": train_loss / train_n,
            "validation_loss": val_loss / val_n,
        }
        history.append(row)
        if row["validation_loss"] < best_val - 1e-5:
            best_val = float(row["validation_loss"])
            best_epoch = epoch
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
        elif epoch - best_epoch >= patience:
            break
        if epoch == 1 or epoch % 10 == 0:
            print(
                f"cvae_epoch={epoch} train={row['train_loss']:.5f} "
                f"val={row['validation_loss']:.5f}",
                flush=True,
            )
    if best_state is None:
        raise RuntimeError("Conditional VAE training did not produce a finite model.")
    model.load_state_dict(best_state)
    model.eval()
    return ConditionalFitResult(model, pd.DataFrame(history), best_epoch)


def encode(model: VAE, x: np.ndarray, batch_size: int = 512) -> np.ndarray:
    output: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(x), batch_size):
            batch = torch.from_numpy(x[start : start + batch_size])
            mu, _ = model.encode(batch)
            output.append(mu.numpy())
    return np.vstack(output)


def decode(model: VAE, z: np.ndarray, batch_size: int = 512) -> np.ndarray:
    output: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(z), batch_size):
            batch = torch.from_numpy(z[start : start + batch_size].astype(np.float32))
            output.append(model.decoder(batch).numpy())
    return np.vstack(output)


def cvae_encode(
    model: ConditionalVAE,
    x: np.ndarray,
    condition: np.ndarray,
    batch_size: int = 512,
) -> np.ndarray:
    output: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(x), batch_size):
            batch = torch.from_numpy(x[start : start + batch_size])
            cond = torch.from_numpy(condition[start : start + batch_size])
            mu, _ = model.encode(batch, cond)
            output.append(mu.numpy())
    return np.vstack(output)


def cvae_decode(
    model: ConditionalVAE,
    z: np.ndarray,
    condition: np.ndarray,
    batch_size: int = 512,
) -> np.ndarray:
    output: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(z), batch_size):
            latent = torch.from_numpy(z[start : start + batch_size].astype(np.float32))
            cond = torch.from_numpy(condition[start : start + batch_size])
            output.append(model.decode(latent, cond).numpy())
    return np.vstack(output)


def make_conditions(metadata: pd.DataFrame) -> np.ndarray:
    age_days = metadata["Age"].map({"P3": 3.0, "P7": 7.0, "P14": 14.0})
    age_scaled = ((age_days.to_numpy(dtype=np.float32) - 7.0) / 5.0)[:, None]
    cell_one_hot = np.column_stack(
        [metadata["CellType"].eq(cell_type).to_numpy(dtype=np.float32) for cell_type in CELL_TYPES]
    )
    oxygen = metadata["Oxygen"].eq("Hyperoxia").to_numpy(dtype=np.float32)[:, None]
    return np.column_stack([age_scaled, cell_one_hot, oxygen]).astype(np.float32)


def animal_means(values: np.ndarray, meta: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    animals = sorted(meta["animal_id"].unique())
    means = np.vstack(
        [values[meta["animal_id"].to_numpy() == animal].mean(axis=0) for animal in animals]
    )
    return means, animals


def balanced_effect(
    values: np.ndarray,
    meta: pd.DataFrame,
    cell_type: str,
    exclude_age: str,
) -> np.ndarray:
    effects: list[np.ndarray] = []
    for age in AGES:
        if age == exclude_age:
            continue
        subset = meta["Age"].eq(age) & meta["CellType"].eq(cell_type)
        age_meta = meta.loc[subset].reset_index(drop=True)
        age_values = values[subset.to_numpy()]
        condition_means: dict[str, np.ndarray] = {}
        for oxygen in ("Normoxia", "Hyperoxia"):
            cond = age_meta["Oxygen"].eq(oxygen)
            pb, _ = animal_means(age_values[cond.to_numpy()], age_meta.loc[cond])
            condition_means[oxygen] = pb.mean(axis=0)
        effects.append(condition_means["Hyperoxia"] - condition_means["Normoxia"])
    return np.vstack(effects).mean(axis=0)


def safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    if method == "pearson":
        return float(pearsonr(x, y).statistic)
    return float(spearmanr(x, y).statistic)


def metric_row(
    pred: np.ndarray,
    observed: np.ndarray,
    normoxia: np.ndarray,
    gene_indices: np.ndarray,
) -> dict[str, float | int]:
    pred_mean = pred.mean(axis=0)[gene_indices]
    obs_mean = observed.mean(axis=0)[gene_indices]
    norm_mean = normoxia.mean(axis=0)[gene_indices]
    pred_effect = pred_mean - norm_mean
    obs_effect = obs_mean - norm_mean
    return {
        "n_genes": int(len(gene_indices)),
        "pearson_expression": safe_corr(pred_mean, obs_mean, "pearson"),
        "spearman_expression": safe_corr(pred_mean, obs_mean, "spearman"),
        "rmse_expression": float(np.sqrt(np.mean((pred_mean - obs_mean) ** 2))),
        "pearson_effect": safe_corr(pred_effect, obs_effect, "pearson"),
        "spearman_effect": safe_corr(pred_effect, obs_effect, "spearman"),
        "rmse_effect": float(np.sqrt(np.mean((pred_effect - obs_effect) ** 2))),
        "direction_accuracy": float(np.mean(np.sign(pred_effect) == np.sign(obs_effect))),
    }


def bootstrap_metrics(
    pred: np.ndarray,
    observed: np.ndarray,
    normoxia: np.ndarray,
    gene_indices: np.ndarray,
    n_bootstrap: int,
    seed: int,
) -> list[dict[str, float | int]]:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, float | int]] = []
    for replicate in range(n_bootstrap):
        pred_i = rng.integers(0, len(pred), len(pred))
        obs_i = rng.integers(0, len(observed), len(observed))
        norm_i = rng.integers(0, len(normoxia), len(normoxia))
        row = metric_row(
            pred[pred_i], observed[obs_i], normoxia[norm_i], gene_indices
        )
        row["bootstrap"] = replicate
        rows.append(row)
    return rows


def choose_genes(
    counts: sparse.csr_matrix,
    genes: np.ndarray,
    signature: pd.DataFrame,
    n_hvg: int,
) -> tuple[sparse.csr_matrix, np.ndarray, pd.DataFrame]:
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    library = np.asarray(counts.sum(axis=1)).ravel()
    if np.any(library <= 0):
        raise ValueError("Cells with zero library size detected.")
    norm = counts.multiply((10000.0 / library)[:, None]).tocsr().astype(np.float32)
    norm.data = np.log1p(norm.data)
    mean = np.asarray(norm.mean(axis=0)).ravel()
    second = np.asarray(norm.power(2).mean(axis=0)).ravel()
    variance = np.maximum(second - mean**2, 0)
    eligible = detected >= max(20, int(round(0.005 * counts.shape[0])))
    dispersion = np.full(len(genes), -np.inf)
    dispersion[eligible] = variance[eligible] / (mean[eligible] + 1e-6)
    hvg_idx = np.argsort(dispersion)[-n_hvg:]
    required = set(signature["gene"].astype(str))
    required.update({
        "Trp53", "Cdkn1a", "Zmat3", "Phlda3", "Gdf15", "Pgf", "Angpt2",
        "Osgin1", "Mdm2", "Bax", "Ccnd1", "Inhba", "Emp2", "Itgb5",
    })
    required_idx = np.flatnonzero(np.isin(genes, sorted(required)))
    selected_idx = np.unique(np.r_[hvg_idx, required_idx])
    selected_idx.sort()
    selected = pd.DataFrame({
        "gene_index_original": selected_idx,
        "gene": genes[selected_idx],
        "detected_cells": detected[selected_idx],
        "mean_log1p_10k": mean[selected_idx],
        "variance_log1p_10k": variance[selected_idx],
        "dispersion": dispersion[selected_idx],
        "selected_as_hvg": np.isin(selected_idx, hvg_idx),
        "selected_as_required": np.isin(selected_idx, required_idx),
    })
    return norm[:, selected_idx].tocsr(), genes[selected_idx], selected


def plot_fold(
    metrics: pd.DataFrame,
    predictions: pd.DataFrame,
    output: Path,
    heldout_age: str,
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    primary = metrics.loc[
        metrics["gene_set"].eq("replicated_198")
        & metrics["cell_type"].isin(["Cap", "Cap-a"])
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    sns.barplot(
        data=primary,
        x="method",
        y="spearman_effect",
        hue="cell_type",
        ax=axes[0],
        palette="Set2",
    )
    axes[0].set_title(f"{heldout_age} held-out effect prediction")
    axes[0].set_xlabel("")
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].set_ylabel("Spearman correlation\n(predicted vs observed effect)")

    scatter = predictions.loc[
        predictions["cell_type"].eq("Cap")
        & predictions["gene_set_33"].eq(True)
        & predictions["method"].eq("cvae_counterfactual")
    ]
    axes[1].axhline(0, color="#999999", lw=0.8)
    axes[1].axvline(0, color="#999999", lw=0.8)
    axes[1].scatter(
        scatter["observed_effect"], scatter["predicted_effect"],
        s=38, color="#355C7D", alpha=0.85,
    )
    for row in scatter.nlargest(6, "observed_effect").itertuples():
        axes[1].annotate(row.gene, (row.observed_effect, row.predicted_effect), fontsize=8)
    axes[1].set_xlabel("Observed hyperoxia effect")
    axes[1].set_ylabel("Conditional-VAE predicted effect")
    axes[1].set_title("33-gene signature: Cap")
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--heldout-ages", nargs="+", default=["P14"])
    parser.add_argument("--n-hvg", type=int, default=1800)
    parser.add_argument("--latent-dim", type=int, default=32)
    parser.add_argument("--max-epochs", type=int, default=80)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--beta", type=float, default=1e-4)
    parser.add_argument("--bootstrap", type=int, default=250)
    parser.add_argument("--seed", type=int, default=20260807)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    set_seed(args.seed)

    counts = sparse.load_npz(args.matrix).tocsr()
    metadata = pd.read_csv(args.metadata, sep="\t")
    genes_frame = pd.read_csv(args.genes, sep="\t")
    genes = genes_frame["gene"].astype(str).to_numpy()
    signature = pd.read_csv(args.signature, sep="\t")
    if counts.shape != (len(metadata), len(genes)):
        raise ValueError("Matrix, metadata, and gene dimensions do not match.")
    norm_sparse, selected_genes, selected = choose_genes(
        counts, genes, signature, args.n_hvg
    )
    selected.to_csv(args.output_dir / "virtual_cell_selected_genes.tsv", sep="\t", index=False)
    x = norm_sparse.toarray().astype(np.float32)
    conditions = make_conditions(metadata)
    gene_to_idx = {gene: i for i, gene in enumerate(selected_genes)}
    sets = {
        "all_model_genes": np.arange(len(selected_genes), dtype=int),
        "replicated_198": np.array(
            [gene_to_idx[g] for g in signature["gene"] if g in gene_to_idx], dtype=int
        ),
        "signature_33": np.array(
            [
                gene_to_idx[row.gene]
                for row in signature.itertuples()
                if row.included_high_confidence_ge3 and row.gene in gene_to_idx
            ],
            dtype=int,
        ),
    }
    all_metrics: list[dict[str, object]] = []
    all_bootstrap: list[dict[str, object]] = []
    all_predictions: list[pd.DataFrame] = []
    audit: dict[str, object] = {
        "seed": args.seed,
        "model": "animal-balanced scGen-like VAE latent shift",
        "baselines": list(METHODS[:-1]),
        "input_cells": int(len(metadata)),
        "modeled_genes": int(len(selected_genes)),
        "heldout_ages": args.heldout_ages,
        "evaluation_unit": "animal pseudobulk",
    }

    for fold_number, heldout_age in enumerate(args.heldout_ages, start=1):
        if heldout_age not in AGES:
            raise ValueError(f"Unsupported age: {heldout_age}")
        fold_dir = args.output_dir / f"heldout_{heldout_age}"
        fold_dir.mkdir(parents=True, exist_ok=True)
        heldout_hyper = metadata["Age"].eq(heldout_age) & metadata["Oxygen"].eq("Hyperoxia")
        train_mask = ~heldout_hyper
        train_mean = x[train_mask].mean(axis=0)
        train_sd = x[train_mask].std(axis=0)
        train_sd[train_sd < 1e-4] = 1.0
        x_scaled = ((x - train_mean) / train_sd).astype(np.float32)
        groups = (
            metadata.loc[train_mask, "animal_id"].astype(str)
            + "__" + metadata.loc[train_mask, "CellType"].astype(str)
        ).to_numpy()
        print(f"training heldout={heldout_age} cells={int(train_mask.sum())}", flush=True)
        fit = train_vae(
            x_scaled[train_mask.to_numpy()], groups, args.seed + fold_number,
            args.max_epochs, args.batch_size, args.beta,
        )
        fit.history.to_csv(fold_dir / "vae_training_history.tsv", sep="\t", index=False)
        torch.save(
            {
                "state_dict": fit.model.state_dict(),
                "selected_genes": selected_genes.tolist(),
                "train_mean": train_mean,
                "train_sd": train_sd,
                "best_epoch": fit.best_epoch,
                "heldout_age": heldout_age,
                "seed": args.seed + fold_number,
            },
            fold_dir / "vae_model.pt",
        )
        cfit = train_cvae(
            x_scaled[train_mask.to_numpy()],
            conditions[train_mask.to_numpy()],
            groups,
            args.seed + 100 + fold_number,
            args.max_epochs,
            args.batch_size,
            args.beta,
        )
        cfit.history.to_csv(
            fold_dir / "conditional_vae_training_history.tsv", sep="\t", index=False
        )
        torch.save(
            {
                "state_dict": cfit.model.state_dict(),
                "selected_genes": selected_genes.tolist(),
                "condition_columns": [
                    "age_scaled", *[f"cell_type_{x}" for x in CELL_TYPES], "hyperoxia"
                ],
                "train_mean": train_mean,
                "train_sd": train_sd,
                "best_epoch": cfit.best_epoch,
                "heldout_age": heldout_age,
                "seed": args.seed + 100 + fold_number,
            },
            fold_dir / "conditional_vae_model.pt",
        )

        pca = PCA(n_components=args.latent_dim, svd_solver="randomized", random_state=args.seed)
        pca_latent = np.full((len(x), args.latent_dim), np.nan, dtype=np.float32)
        pca.fit(x_scaled[train_mask.to_numpy()])
        pca_latent[:] = pca.transform(x_scaled).astype(np.float32)
        vae_latent = encode(fit.model, x_scaled)
        cvae_latent = cvae_encode(cfit.model, x_scaled, conditions)
        fold_predictions: list[pd.DataFrame] = []

        for cell_type in CELL_TYPES:
            target_norm_mask = (
                metadata["Age"].eq(heldout_age)
                & metadata["Oxygen"].eq("Normoxia")
                & metadata["CellType"].eq(cell_type)
            )
            target_hyp_mask = (
                metadata["Age"].eq(heldout_age)
                & metadata["Oxygen"].eq("Hyperoxia")
                & metadata["CellType"].eq(cell_type)
            )
            if target_norm_mask.sum() < 5 or target_hyp_mask.sum() < 5:
                continue
            norm_meta = metadata.loc[target_norm_mask].reset_index(drop=True)
            hyp_meta = metadata.loc[target_hyp_mask].reset_index(drop=True)
            norm_pb, norm_animals = animal_means(x[target_norm_mask.to_numpy()], norm_meta)
            obs_pb, obs_animals = animal_means(x[target_hyp_mask.to_numpy()], hyp_meta)

            gene_delta = balanced_effect(x, metadata, cell_type, heldout_age)
            pca_delta = balanced_effect(pca_latent, metadata, cell_type, heldout_age)
            vae_delta = balanced_effect(vae_latent, metadata, cell_type, heldout_age)

            pred_identity = norm_pb.copy()
            pred_gene = np.maximum(norm_pb + gene_delta, 0)
            norm_pca_cells = pca_latent[target_norm_mask.to_numpy()]
            pred_pca_cells_scaled = pca.inverse_transform(norm_pca_cells + pca_delta)
            pred_pca_cells = np.maximum(
                pred_pca_cells_scaled * train_sd + train_mean, 0
            )
            pred_pca, pca_animals = animal_means(pred_pca_cells, norm_meta)
            norm_vae_cells = vae_latent[target_norm_mask.to_numpy()]
            pred_vae_scaled = decode(fit.model, norm_vae_cells + vae_delta)
            pred_vae_cells = np.maximum(pred_vae_scaled * train_sd + train_mean, 0)
            pred_vae, vae_animals = animal_means(pred_vae_cells, norm_meta)
            norm_condition = conditions[target_norm_mask.to_numpy()].copy()
            hyper_condition = norm_condition.copy()
            hyper_condition[:, -1] = 1.0
            pred_cvae_scaled = cvae_decode(
                cfit.model,
                cvae_latent[target_norm_mask.to_numpy()],
                hyper_condition,
            )
            pred_cvae_cells = np.maximum(
                pred_cvae_scaled * train_sd + train_mean, 0
            )
            pred_cvae, cvae_animals = animal_means(pred_cvae_cells, norm_meta)
            if not (norm_animals == pca_animals == vae_animals == cvae_animals):
                raise ValueError("Predicted animal ordering mismatch.")
            prediction_by_method = {
                "identity": pred_identity,
                "gene_linear": pred_gene,
                "pca_latent": pred_pca,
                "vae_latent": pred_vae,
                "cvae_counterfactual": pred_cvae,
            }
            observed_mean = obs_pb.mean(axis=0)
            norm_mean = norm_pb.mean(axis=0)
            for method, prediction in prediction_by_method.items():
                pred_mean = prediction.mean(axis=0)
                pred_frame = pd.DataFrame({
                    "heldout_age": heldout_age,
                    "cell_type": cell_type,
                    "method": method,
                    "gene": selected_genes,
                    "normoxia_mean": norm_mean,
                    "observed_hyperoxia_mean": observed_mean,
                    "predicted_hyperoxia_mean": pred_mean,
                    "observed_effect": observed_mean - norm_mean,
                    "predicted_effect": pred_mean - norm_mean,
                    "gene_set_198": np.isin(np.arange(len(selected_genes)), sets["replicated_198"]),
                    "gene_set_33": np.isin(np.arange(len(selected_genes)), sets["signature_33"]),
                })
                fold_predictions.append(pred_frame)
                for set_name, indices in sets.items():
                    row = metric_row(prediction, obs_pb, norm_pb, indices)
                    row.update({
                        "heldout_age": heldout_age,
                        "cell_type": cell_type,
                        "method": method,
                        "gene_set": set_name,
                        "n_normoxia_animals": len(norm_pb),
                        "n_hyperoxia_animals": len(obs_pb),
                    })
                    all_metrics.append(row)
                    boot = bootstrap_metrics(
                        prediction, obs_pb, norm_pb, indices,
                        args.bootstrap,
                        args.seed + 1000 * fold_number + 100 * CELL_TYPES.index(cell_type) + METHODS.index(method),
                    )
                    for boot_row in boot:
                        boot_row.update({
                            "heldout_age": heldout_age,
                            "cell_type": cell_type,
                            "method": method,
                            "gene_set": set_name,
                        })
                    all_bootstrap.extend(boot)

        fold_prediction_frame = pd.concat(fold_predictions, ignore_index=True)
        fold_prediction_frame.to_csv(
            fold_dir / "virtual_cell_gene_predictions.tsv.gz",
            sep="\t", index=False, compression="gzip",
        )
        all_predictions.append(fold_prediction_frame)

    metrics = pd.DataFrame(all_metrics)
    bootstrap = pd.DataFrame(all_bootstrap)
    predictions = pd.concat(all_predictions, ignore_index=True)
    metrics.to_csv(args.output_dir / "virtual_cell_benchmark_metrics.tsv", sep="\t", index=False)
    bootstrap.to_csv(
        args.output_dir / "virtual_cell_benchmark_bootstrap.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    predictions.to_csv(
        args.output_dir / "virtual_cell_all_gene_predictions.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    for heldout_age in args.heldout_ages:
        plot_fold(
            metrics.loc[metrics["heldout_age"].eq(heldout_age)],
            predictions.loc[predictions["heldout_age"].eq(heldout_age)],
            args.output_dir / f"virtual_cell_benchmark_{heldout_age}",
            heldout_age,
        )
    audit["results_rows"] = int(len(metrics))
    audit["bootstrap_rows"] = int(len(bootstrap))
    audit["best_epochs"] = {
        age: int(pd.read_csv(
            args.output_dir / f"heldout_{age}" / "vae_training_history.tsv", sep="\t"
        )["validation_loss"].idxmin() + 1)
        for age in args.heldout_ages
    }
    (args.output_dir / "virtual_cell_benchmark_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    summary = metrics.loc[
        metrics["cell_type"].isin(["Cap", "Cap-a"])
        & metrics["gene_set"].eq("replicated_198"),
        ["heldout_age", "cell_type", "method", "spearman_effect", "rmse_effect", "direction_accuracy"],
    ]
    print(summary.to_string(index=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
