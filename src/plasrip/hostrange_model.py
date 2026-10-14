"""Model loading and host scoring for PlasRiP host-range prediction."""

from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

from .exceptions import PlasRipModelError
from .utils.validate import validate_model_dir


class PlasmidTower(nn.Module):
    """k-mer branch for the plasmid tower."""

    def __init__(self, in_dim: int, hidden_dims: list[int], out_dim: int, dropout: float):
        super().__init__()
        layers: list[nn.Module] = []
        last = in_dim
        for h in hidden_dims:
            layers.extend([nn.Linear(last, h), nn.ReLU(), nn.Dropout(dropout)])
            last = h
        layers.append(nn.Linear(last, out_dim))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        z = self.net(x)
        return F.normalize(z, dim=1)


class TwoTowerModel(nn.Module):
    """Symmetric two-tower model matching the training script architecture."""

    def __init__(
        self,
        kmer_dim: int,
        plasmid_feat_dim: int,
        num_hosts: int,
        embed_dim: int,
        hidden_dims: list[int],
        dropout: float,
        host_features_tensor: torch.Tensor,
        plasmid_feat_hidden: int,
        host_feat_hidden: int,
    ):
        super().__init__()
        self.plasmid_tower = PlasmidTower(kmer_dim, hidden_dims, embed_dim, dropout)
        self.plasmid_feat_proj = nn.Sequential(
            nn.Linear(plasmid_feat_dim, plasmid_feat_hidden),
            nn.ReLU(),
            nn.Linear(plasmid_feat_hidden, embed_dim),
        )
        self.plasmid_fuse = nn.Sequential(
            nn.Linear(embed_dim * 2, embed_dim),
            nn.ReLU(),
            nn.Linear(embed_dim, embed_dim),
        )

        self.host_emb = nn.Embedding(num_hosts, embed_dim)
        self.register_buffer("host_features", host_features_tensor)
        self.host_feat_proj = nn.Sequential(
            nn.Linear(host_features_tensor.shape[1], host_feat_hidden),
            nn.ReLU(),
            nn.Linear(host_feat_hidden, embed_dim),
        )
        self.host_fuse = nn.Sequential(
            nn.Linear(embed_dim * 2, embed_dim),
            nn.ReLU(),
            nn.Linear(embed_dim, embed_dim),
        )

    def encode_plasmid(self, x_kmer: torch.Tensor, x_codon_feat: torch.Tensor) -> torch.Tensor:
        e_kmer = self.plasmid_tower(x_kmer)
        e_feat = self.plasmid_feat_proj(x_codon_feat)
        e_fused = self.plasmid_fuse(torch.cat([e_kmer, e_feat], dim=-1))
        return F.normalize(e_fused, dim=-1)

    def encode_host(self, host_ids: torch.Tensor) -> torch.Tensor:
        e_id = self.host_emb(host_ids)
        h_feat = self.host_features[host_ids]
        e_feat = self.host_feat_proj(h_feat)
        e_fused = self.host_fuse(torch.cat([e_id, e_feat], dim=-1))
        return F.normalize(e_fused, dim=-1)

    def score_all_hosts(self, z_plasmid: torch.Tensor) -> torch.Tensor:
        all_ids = torch.arange(self.host_emb.num_embeddings, device=z_plasmid.device)
        all_host_emb = self.encode_host(all_ids)
        return torch.matmul(z_plasmid, all_host_emb.t())


def _pick_validate_key(files: dict[str, str | None], *keys: str) -> str | None:
    for key in keys:
        if key in files and files[key]:
            return str(files[key])
    return None


def _infer_architecture_from_state_dict(state_dict: dict[str, torch.Tensor]) -> dict[str, int | list[int] | float]:
    try:
        kmer_dim = int(state_dict["plasmid_tower.net.0.weight"].shape[1])
        h1 = int(state_dict["plasmid_tower.net.0.weight"].shape[0])
        h2 = int(state_dict["plasmid_tower.net.3.weight"].shape[0])
        embed_dim = int(state_dict["plasmid_tower.net.6.weight"].shape[0])
        plasmid_feat_dim = int(state_dict["plasmid_feat_proj.0.weight"].shape[1])
        plasmid_feat_hidden = int(state_dict["plasmid_feat_proj.0.weight"].shape[0])
        host_feat_hidden = int(state_dict["host_feat_proj.0.weight"].shape[0])
        num_hosts = int(state_dict["host_emb.weight"].shape[0])
        host_feat_dim = int(state_dict["host_features"].shape[1])
    except KeyError as e:
        raise PlasRipModelError(f"Model state dict is missing expected key: {e}") from e

    return {
        "kmer_dim": kmer_dim,
        "hidden_dims": [h1, h2],
        "embed_dim": embed_dim,
        "plasmid_feat_dim": plasmid_feat_dim,
        "plasmid_feat_hidden": plasmid_feat_hidden,
        "host_feat_hidden": host_feat_hidden,
        "num_hosts": num_hosts,
        "host_feat_dim": host_feat_dim,
        "dropout": 0.2,
    }


def load_hostrange_model(model_dir: str | None):
    """Load the trained two-tower model and related resources."""
    files = validate_model_dir(model_dir)

    model_path = _pick_validate_key(files, "model_path", "model/two_tower_best_3.2.pth")
    host_vocab_path = _pick_validate_key(files, "host_vocab_path", "cache/host_vocab.joblib")
    host_profiles_path = _pick_validate_key(files, "host_profiles_path", "cache/host_profiles.joblib")

    if model_path is None:
        raise PlasRipModelError("No model checkpoint path was resolved from model_dir.")
    if host_vocab_path is None:
        raise PlasRipModelError("No host_vocab.joblib path was resolved from model_dir.")

    state_dict = torch.load(model_path, map_location="cpu")
    if not isinstance(state_dict, dict):
        raise PlasRipModelError(f"Unsupported checkpoint type at {model_path}: {type(state_dict)}")

    arch = _infer_architecture_from_state_dict(state_dict)

    if "host_features" not in state_dict:
        raise PlasRipModelError("Checkpoint is missing the host_features buffer.")
    host_features_tensor = state_dict["host_features"].detach().clone().float()

    model = TwoTowerModel(
        kmer_dim=int(arch["kmer_dim"]),
        plasmid_feat_dim=int(arch["plasmid_feat_dim"]),
        num_hosts=int(arch["num_hosts"]),
        embed_dim=int(arch["embed_dim"]),
        hidden_dims=list(arch["hidden_dims"]),
        dropout=float(arch["dropout"]),
        host_features_tensor=host_features_tensor,
        plasmid_feat_hidden=int(arch["plasmid_feat_hidden"]),
        host_feat_hidden=int(arch["host_feat_hidden"]),
    )
    model.load_state_dict(state_dict, strict=True)
    model.eval()

    host_vocab = joblib.load(host_vocab_path)
    if not isinstance(host_vocab, dict) or "id2host_key" not in host_vocab:
        raise PlasRipModelError("host_vocab.joblib has an unexpected structure.")

    host_profiles = None
    if host_profiles_path is not None and Path(host_profiles_path).is_file():
        host_profiles = joblib.load(host_profiles_path)

    root = _pick_validate_key(files, "root") or str(Path(model_path).resolve().parent)

    return {
        "model_dir": str(Path(root).resolve()),
        "model_path": str(Path(model_path).resolve()),
        "model_name": Path(model_path).name,
        "files": files,
        "state_dict": state_dict,
        "arch": arch,
        "model": model,
        "host_vocab": host_vocab,
        "host_profiles": host_profiles,
    }


def load_host_reference(model_dir: str | None) -> dict:
    """Load host reference resources used at inference time."""
    handle = load_hostrange_model(model_dir)
    return {
        "model_dir": handle["model_dir"],
        "host_vocab": handle["host_vocab"],
        "host_profiles": handle.get("host_profiles"),
        "arch": handle["arch"],
    }


def batch_score_hosts(
    plasmid_features,
    host_reference,
    model,
    batch_size: int = 256,
):
    """Score all candidate hosts for all plasmids in batches."""
    if plasmid_features is None:
        raise PlasRipModelError("plasmid_features is None")

    plasmid_ids = list(plasmid_features.get("plasmid_ids", []))
    if not plasmid_ids:
        raise PlasRipModelError("No plasmid IDs were provided for scoring.")

    kmer_features = np.asarray(plasmid_features.get("kmer_features"), dtype=np.float32)
    codon_features = np.asarray(plasmid_features.get("codon_features"), dtype=np.float32)

    if kmer_features.ndim != 2 or codon_features.ndim != 2:
        raise PlasRipModelError("Feature arrays must be 2D.")
    if kmer_features.shape[0] != len(plasmid_ids) or codon_features.shape[0] != len(plasmid_ids):
        raise PlasRipModelError("Feature row counts do not match the number of plasmids.")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    model.eval()

    host_vocab = host_reference["host_vocab"]
    id2host_key = host_vocab["id2host_key"]
    if isinstance(id2host_key, list):
        id2host_key = {i: x for i, x in enumerate(id2host_key)}
    else:
        id2host_key = {int(k): v for k, v in id2host_key.items()}

    rows: list[dict[str, object]] = []
    n = len(plasmid_ids)
    with torch.no_grad():
        for start in range(0, n, batch_size):
            end = min(n, start + batch_size)
            x_kmer = torch.from_numpy(kmer_features[start:end]).to(device).float()
            x_codon = torch.from_numpy(codon_features[start:end]).to(device).float()
            z = model.encode_plasmid(x_kmer, x_codon)
            scores = model.score_all_hosts(z).detach().cpu().numpy()

            for i, plasmid_id in enumerate(plasmid_ids[start:end]):
                for host_id in range(scores.shape[1]):
                    rows.append(
                        {
                            "PLASMID": plasmid_id,
                            "Host_id": int(host_id),
                            "Host_key": id2host_key.get(host_id, str(host_id)),
                            "Score": float(scores[i, host_id]),
                        }
                    )

    return pd.DataFrame(rows)


def standardize_score_output(score_df: pd.DataFrame) -> pd.DataFrame:
    """Ensure score output columns exist and are ordered."""
    expected = ["PLASMID", "Host_id", "Host_key", "Score"]
    out = score_df.copy()
    for col in expected:
        if col not in out.columns:
            out[col] = None
    ordered = expected + [c for c in out.columns if c not in expected]
    return out.loc[:, ordered]


def score_hosts(
    features: dict,
    model_dir: str | None,
    threads: int = 1,
) -> pd.DataFrame:
    """Score candidate hosts for each plasmid using the trained two-tower model."""
    if not features.get("plasmid_ids"):
        raise PlasRipModelError("No plasmid features were provided for scoring.")

    handle = load_hostrange_model(model_dir)
    model = handle["model"]
    host_reference = {
        "host_vocab": handle["host_vocab"],
        "host_profiles": handle.get("host_profiles"),
        "model_dir": handle["model_dir"],
    }

    batch_size = 256
    try:
        if threads and int(threads) > 0:
            batch_size = max(64, min(1024, int(threads) * 8))
    except Exception:
        batch_size = 256

    score_df = batch_score_hosts(features, host_reference, model, batch_size=batch_size)
    return standardize_score_output(score_df)
