#!/usr/bin/env python3
"""Add protocol-adapted scGen, CPA and Sinkhorn OT to the R8 benchmark.

The outer biological split, fold-only HVG selection, animal-grouped validation,
fixed seeds and animal-level bootstrap are identical across methods.  For the
primary GSE151974 cohort, the script reuses and verifies the frozen R8 feature,
validation and bootstrap ledgers.  It can also run the complete seven-method
benchmark in an independent cohort such as GSE243129.

The scGen and CPA implementations preserve the defining algorithms described
in their papers and official source code, but are adapted to Python 3.12 and to
the animal-grouped protocol.  They are deliberately labelled protocol-adapted
implementations rather than executions of the version-incompatible legacy
packages.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import random
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy import sparse
from scipy.io import mmread
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

try:
    import ot
except ImportError as exc:  # pragma: no cover - explicit environment error
    raise ImportError("POT is required for the Sinkhorn OT baseline") from exc


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "scripts" / "27_virtual_cell_leakage_free_benchmark.py"
SPEC = importlib.util.spec_from_file_location("r8_virtual_cell_core", CORE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Unable to import R8 benchmark core from {CORE_PATH}")
core = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = core
SPEC.loader.exec_module(core)

METRICS = core.METRICS
ADDED_METHODS = ("scgen_adapted", "cpa_adapted", "sinkhorn_ot")
FULL_METHODS = (
    "identity",
    "gene_linear",
    "pca_latent",
    "vae_latent",
    "cvae_counterfactual",
) + ADDED_METHODS


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(min(8, max(1, __import__("os").cpu_count() or 1)))


def load_counts(path: Path) -> sparse.csr_matrix:
    if path.name.endswith(".npz"):
        return sparse.load_npz(path).tocsr()
    if path.name.endswith(".mtx") or path.name.endswith(".mtx.gz"):
        return sparse.csr_matrix(mmread(path))
    raise ValueError(f"Unsupported matrix file: {path}")


def global_balanced_effect(
    values: np.ndarray,
    metadata: pd.DataFrame,
    training_ages: tuple[str, ...],
    cell_types: tuple[str, ...],
) -> np.ndarray:
    """Animal-, age- and subtype-balanced perturbation vector for scGen."""
    effects: list[np.ndarray] = []
    for age in training_ages:
        for cell_type in cell_types:
            subset = metadata["Age"].eq(age) & metadata["CellType"].eq(cell_type)
            frame = metadata.loc[subset].reset_index(drop=True)
            matrix = values[subset.to_numpy()]
            condition_means: dict[str, np.ndarray] = {}
            valid = True
            for oxygen in ("Normoxia", "Hyperoxia"):
                keep = frame["Oxygen"].eq(oxygen)
                if not keep.any():
                    valid = False
                    break
                animal_values, _ = core.animal_means(
                    matrix[keep.to_numpy()], frame.loc[keep]
                )
                condition_means[oxygen] = animal_values.mean(axis=0)
            if valid:
                effects.append(condition_means["Hyperoxia"] - condition_means["Normoxia"])
    if not effects:
        raise RuntimeError("No complete training strata were available for scGen arithmetic")
    return np.vstack(effects).mean(axis=0).astype(np.float32)


class GradientReversal(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x: torch.Tensor, strength: float) -> torch.Tensor:
        ctx.strength = strength
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output: torch.Tensor) -> tuple[torch.Tensor, None]:
        return -ctx.strength * grad_output, None


class ProtocolCPA(nn.Module):
    """Protocol-adapted CPA-style model for a binary oxygen perturbation."""

    def __init__(self, n_input: int, n_cell_types: int, n_latent: int) -> None:
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_input, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(256, 128),
            nn.LayerNorm(128),
            nn.GELU(),
            nn.Linear(128, n_latent),
        )
        self.perturbation = nn.Embedding(2, n_latent)
        self.cell_context = nn.Embedding(n_cell_types, n_latent)
        self.decoder = nn.Sequential(
            nn.Linear(n_latent, 128),
            nn.LayerNorm(128),
            nn.GELU(),
            nn.Linear(128, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Linear(256, n_input),
        )
        self.oxygen_adversary = nn.Linear(n_latent, 2)
        self.cell_adversary = nn.Linear(n_latent, n_cell_types)

    def forward(
        self,
        x: torch.Tensor,
        oxygen: torch.Tensor,
        cell_type: torch.Tensor,
        adversarial_strength: float = 0.1,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        basal = self.encoder(x)
        after = basal + self.perturbation(oxygen) + self.cell_context(cell_type)
        reconstruction = self.decoder(after)
        reversed_basal = GradientReversal.apply(basal, adversarial_strength)
        return (
            reconstruction,
            basal,
            self.oxygen_adversary(reversed_basal),
            self.cell_adversary(reversed_basal),
        )

    def predict(
        self, x: torch.Tensor, target_oxygen: torch.Tensor, cell_type: torch.Tensor
    ) -> torch.Tensor:
        basal = self.encoder(x)
        after = basal + self.perturbation(target_oxygen) + self.cell_context(cell_type)
        return self.decoder(after)


@dataclass
class CPASelection:
    model: ProtocolCPA
    best_epoch: int
    history: pd.DataFrame


def cpa_loss(
    model: ProtocolCPA,
    x: torch.Tensor,
    oxygen: torch.Tensor,
    cell_type: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    reconstruction, _, oxygen_logits, cell_logits = model(x, oxygen, cell_type)
    mse = torch.mean((reconstruction - x) ** 2, dim=1)
    oxygen_ce = nn.functional.cross_entropy(oxygen_logits, oxygen, reduction="none")
    cell_ce = nn.functional.cross_entropy(cell_logits, cell_type, reduction="none")
    return mse + 0.05 * oxygen_ce + 0.05 * cell_ce, mse


def train_cpa_with_grouped_validation(
    x: np.ndarray,
    oxygen: np.ndarray,
    cell_type: np.ndarray,
    groups: np.ndarray,
    train_idx: np.ndarray,
    validation_idx: np.ndarray,
    seed: int,
    n_latent: int,
    max_epochs: int,
    batch_size: int,
) -> CPASelection:
    set_seed(seed)
    model = ProtocolCPA(x.shape[1], int(cell_type.max()) + 1, n_latent)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-5)
    dataset = TensorDataset(
        torch.from_numpy(x[train_idx]),
        torch.from_numpy(oxygen[train_idx]),
        torch.from_numpy(cell_type[train_idx]),
    )
    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        sampler=core.group_balanced_sampler(groups[train_idx], seed),
    )
    history: list[dict[str, float | int]] = []
    best_state: dict[str, torch.Tensor] | None = None
    best_epoch = -1
    best_validation = math.inf
    patience = 8
    validation_groups = groups[validation_idx]
    for epoch in range(1, max_epochs + 1):
        model.train()
        total = 0.0
        n_total = 0
        for batch_x, batch_oxygen, batch_cell_type in loader:
            optimizer.zero_grad(set_to_none=True)
            losses, _ = cpa_loss(model, batch_x, batch_oxygen, batch_cell_type)
            loss = losses.mean()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            total += float(loss.detach()) * len(batch_x)
            n_total += len(batch_x)
        model.eval()
        with torch.no_grad():
            _, validation_mse = cpa_loss(
                model,
                torch.from_numpy(x[validation_idx]),
                torch.from_numpy(oxygen[validation_idx]),
                torch.from_numpy(cell_type[validation_idx]),
            )
        frame = pd.DataFrame(
            {"animal_id": validation_groups, "loss": validation_mse.cpu().numpy()}
        )
        validation_loss = float(frame.groupby("animal_id")["loss"].mean().mean())
        history.append(
            {
                "epoch": epoch,
                "train_loss": total / max(n_total, 1),
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
        raise RuntimeError("CPA did not produce a finite grouped-validation fit")
    model.load_state_dict(best_state)
    model.eval()
    return CPASelection(model, best_epoch, pd.DataFrame(history))


def refit_cpa(
    x: np.ndarray,
    oxygen: np.ndarray,
    cell_type: np.ndarray,
    groups: np.ndarray,
    outer_train_idx: np.ndarray,
    seed: int,
    n_latent: int,
    epochs: int,
    batch_size: int,
) -> ProtocolCPA:
    set_seed(seed)
    model = ProtocolCPA(x.shape[1], int(cell_type.max()) + 1, n_latent)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-5)
    dataset = TensorDataset(
        torch.from_numpy(x[outer_train_idx]),
        torch.from_numpy(oxygen[outer_train_idx]),
        torch.from_numpy(cell_type[outer_train_idx]),
    )
    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        sampler=core.group_balanced_sampler(groups[outer_train_idx], seed),
    )
    for _ in range(epochs):
        model.train()
        for batch_x, batch_oxygen, batch_cell_type in loader:
            optimizer.zero_grad(set_to_none=True)
            losses, _ = cpa_loss(model, batch_x, batch_oxygen, batch_cell_type)
            losses.mean().backward()
            nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
    model.eval()
    return model


def predict_cpa(
    model: ProtocolCPA,
    x: np.ndarray,
    cell_type: np.ndarray,
    batch_size: int = 512,
) -> np.ndarray:
    output: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(x), batch_size):
            batch = torch.from_numpy(x[start : start + batch_size])
            target = torch.ones(len(batch), dtype=torch.long)
            context = torch.from_numpy(cell_type[start : start + batch_size])
            output.append(model.predict(batch, target, context).cpu().numpy())
    return np.vstack(output)


def animal_balanced_cells(
    values: np.ndarray,
    metadata: pd.DataFrame,
    seed: int,
    cap_per_animal: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    animals = sorted(metadata["animal_id"].astype(str).unique())
    animal_vector = metadata["animal_id"].astype(str).to_numpy()
    for animal in animals:
        candidates = np.flatnonzero(animal_vector == animal)
        n = min(len(candidates), cap_per_animal)
        chosen = rng.choice(candidates, size=n, replace=False)
        selected.append(values[chosen])
        weights.append(np.full(n, 1.0 / (len(animals) * n), dtype=float))
    return np.vstack(selected), np.concatenate(weights)


def sinkhorn_map_one_age(
    source: np.ndarray,
    source_weights: np.ndarray,
    target: np.ndarray,
    target_weights: np.ndarray,
    query: np.ndarray,
    regularization: float,
    neighbors: int,
) -> np.ndarray:
    cost = ot.dist(source, target, metric="sqeuclidean")
    positive = cost[cost > 0]
    scale = float(np.median(positive)) if len(positive) else 1.0
    cost = cost / max(scale, 1e-8)
    coupling = ot.sinkhorn(
        source_weights,
        target_weights,
        cost,
        reg=regularization,
        method="sinkhorn_log",
        numItermax=2000,
        stopThr=1e-7,
        warn=False,
    )
    row_mass = np.maximum(coupling.sum(axis=1), 1e-12)
    mapped_source = (coupling @ target) / row_mass[:, None]
    displacement = mapped_source - source
    n_neighbors = min(neighbors, len(source))
    nearest = NearestNeighbors(n_neighbors=n_neighbors).fit(source)
    distance, index = nearest.kneighbors(query)
    weight = 1.0 / np.maximum(distance, 1e-6)
    weight /= weight.sum(axis=1, keepdims=True)
    query_displacement = np.sum(displacement[index] * weight[:, :, None], axis=1)
    return query + query_displacement


def sinkhorn_prediction(
    pca_latent: np.ndarray,
    metadata: pd.DataFrame,
    query_mask: np.ndarray,
    cell_type: str,
    training_ages: tuple[str, ...],
    seed: int,
    cap_per_animal: int,
    regularization: float,
    neighbors: int,
) -> np.ndarray:
    predictions: list[np.ndarray] = []
    query = pca_latent[query_mask]
    for age_index, age in enumerate(training_ages):
        source_mask = (
            metadata["Age"].eq(age)
            & metadata["Oxygen"].eq("Normoxia")
            & metadata["CellType"].eq(cell_type)
        )
        target_mask = (
            metadata["Age"].eq(age)
            & metadata["Oxygen"].eq("Hyperoxia")
            & metadata["CellType"].eq(cell_type)
        )
        source_meta = metadata.loc[source_mask].reset_index(drop=True)
        target_meta = metadata.loc[target_mask].reset_index(drop=True)
        source, source_weights = animal_balanced_cells(
            pca_latent[source_mask.to_numpy()],
            source_meta,
            seed + age_index * 1009,
            cap_per_animal,
        )
        target, target_weights = animal_balanced_cells(
            pca_latent[target_mask.to_numpy()],
            target_meta,
            seed + age_index * 1009 + 503,
            cap_per_animal,
        )
        predictions.append(
            sinkhorn_map_one_age(
                source,
                source_weights,
                target,
                target_weights,
                query,
                regularization,
                neighbors,
            )
        )
    return np.stack(predictions).mean(axis=0).astype(np.float32)


def endpoint_indices(manifest: pd.DataFrame) -> dict[str, np.ndarray]:
    return core.endpoint_indices(manifest)


def feature_manifest(
    counts: sparse.csr_matrix,
    genes: np.ndarray,
    metadata: pd.DataFrame,
    heldout_age: str,
    replicated_genes: set[str],
    signature_genes: set[str],
    n_hvg: int,
    reuse_dir: Path | None,
) -> tuple[np.ndarray, pd.DataFrame]:
    outer_train_mask = metadata["Age"].ne(heldout_age).to_numpy()
    selected_idx, _, recomputed = core.select_outer_fold_features(
        counts,
        genes,
        outer_train_mask,
        replicated_genes,
        signature_genes,
        n_hvg,
    )
    if reuse_dir is None:
        manifest = recomputed.copy()
        manifest.insert(0, "heldout_age", heldout_age)
        manifest.insert(1, "model_feature_position", np.arange(len(manifest)))
        return selected_idx, manifest
    frozen = pd.read_csv(
        reuse_dir / "virtual_cell_fold_feature_manifest.tsv", sep="\t"
    )
    manifest = frozen[frozen["heldout_age"].eq(heldout_age)].copy()
    manifest = manifest.sort_values("model_feature_position").reset_index(drop=True)
    if manifest["gene"].astype(str).tolist() != recomputed["gene"].astype(str).tolist():
        raise RuntimeError(f"Recomputed features differ from frozen R8 ledger for {heldout_age}")
    return selected_idx, manifest


def validation_animals_for_seed(
    metadata: pd.DataFrame,
    outer_train_mask: np.ndarray,
    heldout_age: str,
    seed: int,
    fold_number: int,
    reuse_dir: Path | None,
) -> tuple[set[str], pd.DataFrame]:
    if reuse_dir is None:
        return core.stratified_validation_animals(
            metadata, outer_train_mask, seed + fold_number * 1000
        )
    ledger = pd.read_csv(
        reuse_dir / "virtual_cell_grouped_validation_splits.tsv", sep="\t"
    )
    split = ledger[
        ledger["heldout_age"].eq(heldout_age) & ledger["seed"].eq(seed)
    ].copy()
    if split.empty:
        raise RuntimeError(f"Missing frozen R8 validation split for {heldout_age}, {seed}")
    validation = set(split.loc[split["role"].eq("validation"), "animal_id"].astype(str))
    return validation, split.drop(columns=["heldout_age", "seed"])


def build_bootstrap_draws(
    heldout_age: str,
    fold_number: int,
    task_norm_animals: dict[str, list[str]],
    task_obs_animals: dict[str, list[str]],
    capillary_types: tuple[str, ...],
    n_bootstrap: int,
    base_seed: int,
    reuse_dir: Path | None,
) -> tuple[list[tuple[np.ndarray, np.ndarray]], list[dict[str, object]]]:
    capillary_tasks = [f"{heldout_age}__{cell_type}" for cell_type in capillary_types]
    common_norm = sorted(
        set.intersection(*(set(task_norm_animals[task]) for task in capillary_tasks))
    )
    common_obs = sorted(
        set.intersection(*(set(task_obs_animals[task]) for task in capillary_tasks))
    )
    rows: list[dict[str, object]] = []
    draws: list[tuple[np.ndarray, np.ndarray]] = []
    if reuse_dir is not None:
        ledger = pd.read_csv(
            reuse_dir / "virtual_cell_joint_capillary_bootstrap_draws.tsv.gz", sep="\t"
        )
        ledger = ledger[ledger["heldout_age"].eq(heldout_age)]
        for replicate in range(n_bootstrap):
            subset = ledger[ledger["bootstrap"].eq(replicate)]
            norm = (
                subset[subset["oxygen_group"].eq("Normoxia")]
                .sort_values("draw_position")["animal_id"]
                .astype(str)
                .to_numpy()
            )
            obs = (
                subset[subset["oxygen_group"].eq("Hyperoxia")]
                .sort_values("draw_position")["animal_id"]
                .astype(str)
                .to_numpy()
            )
            draws.append((norm, obs))
        return draws, ledger.to_dict("records")
    rng = np.random.default_rng(base_seed + fold_number * 1000)
    for replicate in range(n_bootstrap):
        norm = rng.choice(common_norm, size=len(common_norm), replace=True)
        obs = rng.choice(common_obs, size=len(common_obs), replace=True)
        draws.append((norm, obs))
        for group, identifiers in (("Normoxia", norm), ("Hyperoxia", obs)):
            rows.extend(
                {
                    "heldout_age": heldout_age,
                    "bootstrap": replicate,
                    "oxygen_group": group,
                    "draw_position": position,
                    "animal_id": str(animal),
                }
                for position, animal in enumerate(identifiers)
            )
    return draws, rows


def run(args: argparse.Namespace) -> int:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    counts = load_counts(args.matrix)
    metadata = pd.read_csv(args.metadata, sep="\t")
    genes = pd.read_csv(args.genes, sep="\t")["gene"].astype(str).to_numpy()
    signature = pd.read_csv(args.signature, sep="\t")
    if counts.shape != (len(metadata), len(genes)):
        raise ValueError(f"Dimension mismatch: {counts.shape}, {len(metadata)}, {len(genes)}")
    replicated_genes = set(signature["gene"].astype(str).unique())
    signature_genes = set(
        signature.loc[
            signature["concordant_endothelial_subtypes"].ge(3), "gene"
        ].astype(str)
    )
    if len(replicated_genes) != 198 or len(signature_genes) != 33:
        raise ValueError("The locked 198/33 sets were not recovered")
    cell_types = tuple(args.cell_types)
    model_cell_types = tuple(
        cell_type
        for cell_type in core.CELL_TYPES
        if cell_type in set(metadata["CellType"].astype(str))
    )
    model_cell_types += tuple(
        sorted(set(metadata["CellType"].astype(str)) - set(model_cell_types))
    )
    capillary_types = tuple(ct for ct in ("Cap", "Cap-a") if ct in cell_types)
    missing = set(cell_types) - set(metadata["CellType"].astype(str))
    if missing:
        raise ValueError(f"Requested cell types absent from metadata: {sorted(missing)}")
    logical_seeds = list(args.seeds)
    if len(logical_seeds) < 10 and not args.smoke_test:
        raise ValueError("Ten fixed seeds are required for a reportable run")

    point_rows: list[dict[str, object]] = []
    bootstrap_rows: list[dict[str, object]] = []
    prediction_rows: list[pd.DataFrame] = []
    feature_rows: list[pd.DataFrame] = []
    split_rows: list[pd.DataFrame] = []
    history_rows: list[pd.DataFrame] = []
    bootstrap_draw_rows: list[dict[str, object]] = []
    array_payload: dict[str, np.ndarray] = {}
    fold_audit: dict[str, object] = {}

    for fold_number, heldout_age in enumerate(args.heldout_ages, start=1):
        outer_train_mask = metadata["Age"].ne(heldout_age).to_numpy()
        outer_train_idx = np.flatnonzero(outer_train_mask)
        training_ages = tuple(
            age for age in sorted(metadata["Age"].astype(str).unique()) if age != heldout_age
        )
        selected_idx, manifest = feature_manifest(
            counts,
            genes,
            metadata,
            heldout_age,
            replicated_genes,
            signature_genes,
            args.n_hvg,
            args.reuse_r8_dir,
        )
        feature_rows.append(manifest)
        selected_genes = manifest["gene"].astype(str).to_numpy()
        endpoints = endpoint_indices(manifest)
        x = core.normalize_selected_counts(counts, selected_idx)
        train_mean = x[outer_train_mask].mean(axis=0)
        train_sd = x[outer_train_mask].std(axis=0)
        train_sd[train_sd < 1e-4] = 1.0
        x_scaled = ((x - train_mean) / train_sd).astype(np.float32)
        groups = metadata["animal_id"].astype(str).to_numpy()
        oxygen = metadata["Oxygen"].eq("Hyperoxia").to_numpy(dtype=np.int64)
        cell_map = {cell_type: index for index, cell_type in enumerate(model_cell_types)}
        cell_index = metadata["CellType"].map(cell_map).to_numpy(dtype=np.int64)
        conditions = core.make_conditions(metadata)
        pca_components = min(args.latent_dim, len(selected_idx) - 1, len(outer_train_idx) - 1)
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
        task_norm_masks: dict[str, np.ndarray] = {}
        for cell_type in cell_types:
            norm_mask = (
                metadata["Age"].eq(heldout_age)
                & metadata["Oxygen"].eq("Normoxia")
                & metadata["CellType"].eq(cell_type)
            ).to_numpy()
            hyp_mask = (
                metadata["Age"].eq(heldout_age)
                & metadata["Oxygen"].eq("Hyperoxia")
                & metadata["CellType"].eq(cell_type)
            ).to_numpy()
            if norm_mask.sum() < 4 or hyp_mask.sum() < 4:
                continue
            norm_meta = metadata.loc[norm_mask].reset_index(drop=True)
            hyp_meta = metadata.loc[hyp_mask].reset_index(drop=True)
            norm_pb, norm_animals = core.animal_means(x[norm_mask], norm_meta)
            obs_pb, obs_animals = core.animal_means(x[hyp_mask], hyp_meta)
            gene_delta = core.balanced_effect(x, metadata, cell_type, training_ages)
            pca_delta = core.balanced_effect(
                pca_latent, metadata, cell_type, training_ages
            )
            pca_pred_scaled = pca.inverse_transform(pca_latent[norm_mask] + pca_delta)
            pca_pred_cells = np.maximum(pca_pred_scaled * train_sd + train_mean, 0)
            pca_pred, _ = core.animal_means(pca_pred_cells, norm_meta)
            task = f"{heldout_age}__{cell_type}"
            task_observed[task] = obs_pb
            task_normoxia[task] = norm_pb
            task_norm_animals[task] = norm_animals
            task_obs_animals[task] = obs_animals
            task_norm_masks[task] = norm_mask
            simple_predictions[task] = {
                "identity": norm_pb.copy(),
                "gene_linear": np.maximum(norm_pb + gene_delta, 0),
                "pca_latent": pca_pred,
            }

        methods_to_train = ["scgen_adapted", "cpa_adapted"]
        if args.include_existing_models:
            methods_to_train = ["vae_latent", "cvae_counterfactual"] + methods_to_train
        deep_predictions: dict[str, dict[int, dict[str, np.ndarray]]] = {
            task: {seed: {} for seed in logical_seeds} for task in simple_predictions
        }
        for logical_seed in logical_seeds:
            validation_animals, split = validation_animals_for_seed(
                metadata,
                outer_train_mask,
                heldout_age,
                logical_seed,
                fold_number,
                args.reuse_r8_dir,
            )
            split = split.copy()
            split.insert(0, "heldout_age", heldout_age)
            split.insert(1, "seed", logical_seed)
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
                raise RuntimeError("Animal leakage between fit and validation")
            fitted: dict[str, nn.Module] = {}
            best_epochs: dict[str, int] = {}
            for method_index, method in enumerate(methods_to_train):
                model_seed = logical_seed + fold_number * 10000 + (method_index + 1) * 200000
                if method in ("vae_latent", "cvae_counterfactual", "scgen_adapted"):
                    kind = "vae_latent" if method == "scgen_adapted" else method
                    selection = core.train_with_grouped_validation(
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
                    fitted[method] = core.refit_all_outer_training(
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
                    history = selection.history.copy()
                    best_epochs[method] = selection.best_epoch
                else:
                    selection = train_cpa_with_grouped_validation(
                        x_scaled,
                        oxygen,
                        cell_index,
                        groups,
                        fit_train_idx,
                        validation_idx,
                        model_seed,
                        args.latent_dim,
                        args.max_epochs,
                        args.batch_size,
                    )
                    fitted[method] = refit_cpa(
                        x_scaled,
                        oxygen,
                        cell_index,
                        groups,
                        outer_train_idx,
                        model_seed,
                        args.latent_dim,
                        selection.best_epoch,
                        args.batch_size,
                    )
                    history = selection.history.copy()
                    best_epochs[method] = selection.best_epoch
                history.insert(0, "heldout_age", heldout_age)
                history.insert(1, "seed", logical_seed)
                history.insert(2, "method", method)
                history["selected_best_epoch"] = best_epochs[method]
                history_rows.append(history)
                print(
                    f"cohort={args.cohort} fold={heldout_age} seed={logical_seed} "
                    f"method={method} best_epoch={best_epochs[method]}",
                    flush=True,
                )

            vae_latent = (
                core.encode_vae(fitted["vae_latent"], x_scaled)
                if "vae_latent" in fitted
                else None
            )
            cvae_latent = (
                core.encode_cvae(fitted["cvae_counterfactual"], x_scaled, conditions)
                if "cvae_counterfactual" in fitted
                else None
            )
            scgen_latent = core.encode_vae(fitted["scgen_adapted"], x_scaled)
            scgen_delta = global_balanced_effect(
                scgen_latent, metadata, training_ages, model_cell_types
            )
            for cell_type in cell_types:
                task = f"{heldout_age}__{cell_type}"
                if task not in simple_predictions:
                    continue
                norm_mask = task_norm_masks[task]
                norm_meta = metadata.loc[norm_mask].reset_index(drop=True)
                predictions: dict[str, np.ndarray] = {}
                if vae_latent is not None:
                    vae_delta = core.balanced_effect(
                        vae_latent, metadata, cell_type, training_ages
                    )
                    decoded = core.decode_vae(
                        fitted["vae_latent"], vae_latent[norm_mask] + vae_delta
                    )
                    pred_cells = np.maximum(decoded * train_sd + train_mean, 0)
                    predictions["vae_latent"], _ = core.animal_means(pred_cells, norm_meta)
                if cvae_latent is not None:
                    hyper_condition = conditions[norm_mask].copy()
                    hyper_condition[:, -1] = 1.0
                    decoded = core.decode_cvae(
                        fitted["cvae_counterfactual"],
                        cvae_latent[norm_mask],
                        hyper_condition,
                    )
                    pred_cells = np.maximum(decoded * train_sd + train_mean, 0)
                    predictions["cvae_counterfactual"], _ = core.animal_means(
                        pred_cells, norm_meta
                    )
                decoded = core.decode_vae(
                    fitted["scgen_adapted"], scgen_latent[norm_mask] + scgen_delta
                )
                pred_cells = np.maximum(decoded * train_sd + train_mean, 0)
                predictions["scgen_adapted"], _ = core.animal_means(pred_cells, norm_meta)
                decoded = predict_cpa(
                    fitted["cpa_adapted"],
                    x_scaled[norm_mask],
                    cell_index[norm_mask],
                )
                pred_cells = np.maximum(decoded * train_sd + train_mean, 0)
                predictions["cpa_adapted"], _ = core.animal_means(pred_cells, norm_meta)
                ot_latent = sinkhorn_prediction(
                    pca_latent,
                    metadata,
                    norm_mask,
                    cell_type,
                    training_ages,
                    logical_seed + fold_number * 10000 + 900000,
                    args.ot_cap_per_animal,
                    args.ot_regularization,
                    args.ot_neighbors,
                )
                ot_scaled = pca.inverse_transform(ot_latent)
                ot_cells = np.maximum(ot_scaled * train_sd + train_mean, 0)
                predictions["sinkhorn_ot"], _ = core.animal_means(ot_cells, norm_meta)
                deep_predictions[task][logical_seed] = predictions

        draws, draw_rows = build_bootstrap_draws(
            heldout_age,
            fold_number,
            task_norm_animals,
            task_obs_animals,
            capillary_types,
            args.bootstrap,
            args.base_seed,
            args.reuse_r8_dir,
        )
        bootstrap_draw_rows.extend(draw_rows)

        for task, simple in simple_predictions.items():
            _, cell_type = task.split("__", 1)
            observed = task_observed[task]
            normoxia = task_normoxia[task]
            methods: list[tuple[str, int, np.ndarray]] = []
            if args.include_existing_models:
                methods.extend((method, -1, pred) for method, pred in simple.items())
            methods.extend(
                (method, seed, prediction)
                for seed in logical_seeds
                for method, prediction in deep_predictions[task][seed].items()
            )
            for method, seed, prediction in methods:
                for endpoint, indices in endpoints.items():
                    row = core.metric_row(prediction, observed, normoxia, indices)
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
                prediction_rows.append(
                    pd.DataFrame(
                        {
                            "heldout_age": heldout_age,
                            "cell_type": cell_type,
                            "method": method,
                            "seed": seed,
                            "gene": selected_genes,
                            "normoxia_mean": normoxia.mean(axis=0),
                            "observed_hyperoxia_mean": observed.mean(axis=0),
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
                if cell_type in capillary_types:
                    key = f"{heldout_age}_{cell_type.replace('-', 'a')}"
                    array_payload[f"{key}_{method}_seed{seed}"] = prediction

            if cell_type in capillary_types:
                norm_lookup = {
                    animal: index for index, animal in enumerate(task_norm_animals[task])
                }
                obs_lookup = {
                    animal: index for index, animal in enumerate(task_obs_animals[task])
                }
                for replicate, (norm_ids, obs_ids) in enumerate(draws):
                    norm_idx = np.asarray([norm_lookup[str(x)] for x in norm_ids], dtype=int)
                    obs_idx = np.asarray([obs_lookup[str(x)] for x in obs_ids], dtype=int)
                    for method, seed, prediction in methods:
                        for endpoint, indices in endpoints.items():
                            row = core.metric_row(
                                prediction[norm_idx],
                                observed[obs_idx],
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
            "training_ages": list(training_ages),
            "heldout_age_excluded_from_feature_selection_fit_and_tuning": True,
            "n_model_features": int(len(manifest)),
            "n_training_animals": int(metadata.loc[outer_train_mask, "animal_id"].nunique()),
            "n_test_animals": int(metadata.loc[~outer_train_mask, "animal_id"].nunique()),
            "feature_sha256": hashlib.sha256(
                "\n".join(selected_genes).encode("utf-8")
            ).hexdigest(),
        }

    point = pd.DataFrame(point_rows)
    bootstrap = pd.DataFrame(bootstrap_rows)
    predictions = pd.concat(prediction_rows, ignore_index=True)
    features = pd.concat(feature_rows, ignore_index=True)
    splits = pd.concat(split_rows, ignore_index=True)
    histories = pd.concat(history_rows, ignore_index=True)
    point.to_csv(args.output_dir / "expanded_point_metrics.tsv", sep="\t", index=False)
    bootstrap.to_csv(
        args.output_dir / "expanded_strictly_paired_bootstrap.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    predictions.to_csv(
        args.output_dir / "expanded_gene_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    features.to_csv(args.output_dir / "expanded_fold_feature_manifest.tsv", sep="\t", index=False)
    splits.to_csv(args.output_dir / "expanded_grouped_validation_splits.tsv", sep="\t", index=False)
    histories.to_csv(
        args.output_dir / "expanded_training_histories.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(bootstrap_draw_rows).to_csv(
        args.output_dir / "expanded_joint_capillary_bootstrap_draws.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    np.savez_compressed(args.output_dir / "expanded_capillary_animal_arrays.npz", **array_payload)
    audit = {
        "analysis_version": "R9 expanded perturbation benchmark",
        "cohort": args.cohort,
        "protocol_adapted_methods": {
            "scgen_adapted": "VAE plus animal-, age- and subtype-balanced latent vector arithmetic",
            "cpa_adapted": "basal encoder plus additive perturbation/context embeddings and adversarial disentanglement",
            "sinkhorn_ot": "animal-balanced entropic OT in fold-trained PCA space with k-nearest out-of-sample displacement",
        },
        "legacy_package_compatibility_note": (
            "The official scGen 2.1 and CPA 0.8 packages target older scvi/Python/Torch APIs. "
            "Protocol-adapted implementations preserving their defining perturbation operators "
            "were used so every model shares Python 3.12, the same features, animal splits and "
            "training budget; these runs are not presented as executions of the official packages."
        ),
        "outer_folds": args.heldout_ages,
        "logical_training_seeds": logical_seeds,
        "n_training_seeds": len(logical_seeds),
        "bootstrap_replicates": args.bootstrap,
        "reused_frozen_r8_ledgers": bool(args.reuse_r8_dir is not None),
        "model_feature_space": f"exactly {args.n_hvg} outer-training HVGs per fold",
        "ot_parameters": {
            "regularization_relative_to_median_cost": args.ot_regularization,
            "neighbors": args.ot_neighbors,
            "maximum_cells_per_animal": args.ot_cap_per_animal,
        },
        "folds": fold_audit,
        "output_rows": {
            "point_metrics": int(len(point)),
            "bootstrap_metrics": int(len(bootstrap)),
            "gene_predictions": int(len(predictions)),
        },
    }
    (args.output_dir / "expanded_benchmark_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(
        point[
            point["cell_type"].isin(capillary_types)
            & point["endpoint"].eq("primary_fold_hvg")
        ]
        .groupby("method", observed=True)["spearman_effect"]
        .mean()
        .sort_values(ascending=False)
        .to_string(),
        flush=True,
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--heldout-ages", nargs="+", required=True)
    parser.add_argument("--cell-types", nargs="+", required=True)
    parser.add_argument("--reuse-r8-dir", type=Path)
    parser.add_argument("--include-existing-models", action="store_true")
    parser.add_argument("--n-hvg", type=int, default=1800)
    parser.add_argument("--latent-dim", type=int, default=32)
    parser.add_argument("--max-epochs", type=int, default=60)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--beta", type=float, default=1e-4)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--seeds", nargs="+", type=int, default=list(range(20260817, 20260827)))
    parser.add_argument("--base-seed", type=int, default=20260817)
    parser.add_argument("--ot-regularization", type=float, default=0.05)
    parser.add_argument("--ot-neighbors", type=int, default=15)
    parser.add_argument("--ot-cap-per-animal", type=int, default=64)
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
