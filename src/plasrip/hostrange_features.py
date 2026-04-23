"""Feature extraction for PlasRiP v2 host-range prediction.

This module ports the preprocessing logic used by the two-tower host-range
model. It computes, for each input plasmid FASTA:

1. k-mer frequency features (default: 6-mer)
2. Prodigal-derived CDS features with fixed ``-p meta`` mode
3. Codon usage features: RSCU + GC3 + ENC
4. Application of the training-time plasmid codon scaler

Notes
-----
- Prodigal is invoked with ``-p meta`` by design, matching the user's request.
- If Prodigal or downstream codon parsing fails for an input plasmid, a default
  codon feature vector is used so the whole workflow can continue.
- The model directory may use either a structured or flat layout; both are
  handled here.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
from Bio import SeqIO

from .exceptions import PlasRipDependencyError
from .utils.fasta import get_plasmid_id_from_fasta, validate_fastas
from .utils.validate import validate_model_dir


# ---------------------------------------------------------------------------
# Genetic code and codon features
# ---------------------------------------------------------------------------
GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

CODONS = sorted([c for c, aa in GENETIC_CODE.items() if aa != "*"])
AA2CODONS: dict[str, list[str]] = defaultdict(list)
for _codon, _aa in GENETIC_CODE.items():
    if _aa != "*":
        AA2CODONS[_aa].append(_codon)

AA_DEGENERACY = {
    "F": 2, "Y": 2, "C": 2, "H": 2, "Q": 2, "N": 2, "K": 2, "D": 2, "E": 2,
    "I": 3,
    "V": 4, "P": 4, "T": 4, "A": 4, "G": 4,
    "L": 6, "S": 6, "R": 6,
    "M": 1, "W": 1,
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _get_validate_key(files: dict[str, str | None], *candidates: str) -> str | None:
    """Return the first existing key from a validate_model_dir result."""
    for key in candidates:
        if key in files and files[key]:
            return str(files[key])
    return None


def load_feature_preprocessors(model_dir: str) -> dict[str, Any]:
    """Load feature preprocessors/resources from a model directory.

    Returns
    -------
    dict
        Keys:
            plasmid_codon_scaler
            host_vocab
            host_profiles
            files
    """
    files = validate_model_dir(model_dir)

    scaler_path = _get_validate_key(
        files,
        "plasmid_codon_scaler_path",
        "cache/plasmid_codon_scaler.joblib",
    )
    host_vocab_path = _get_validate_key(
        files,
        "host_vocab_path",
        "cache/host_vocab.joblib",
    )
    host_profiles_path = _get_validate_key(
        files,
        "host_profiles_path",
        "cache/host_profiles.joblib",
    )

    out: dict[str, Any] = {"files": files}
    if scaler_path:
        out["plasmid_codon_scaler"] = joblib.load(scaler_path)
    else:
        raise KeyError("plasmid_codon_scaler_path")

    if host_vocab_path:
        out["host_vocab"] = joblib.load(host_vocab_path)
    if host_profiles_path:
        out["host_profiles"] = joblib.load(host_profiles_path)
    return out


def default_codon_feature_vector() -> np.ndarray:
    """Default codon feature vector when CDS extraction/parsing fails."""
    rscu = [1.0] * len(CODONS)
    gc3 = 0.5
    enc = 40.0
    return np.array(rscu + [gc3, enc], dtype=np.float32)


def default_host_feature_vector() -> np.ndarray:
    """Default host-side feature vector with the same dimensionality."""
    return default_codon_feature_vector()


def normalize_name(name: str) -> str:
    """Normalize host/species strings for downstream lookup."""
    x = os.path.splitext(os.path.basename(str(name)))[0]
    x = re.sub(r"_\d+$", "", x)
    x = x.strip().lower()
    x = re.sub(r"[^a-z0-9]+", "_", x)
    x = re.sub(r"_+", "_", x).strip("_")
    return x


# ---------------------------------------------------------------------------
# k-mer features
# ---------------------------------------------------------------------------
def _product(bases: list[str], repeat: int):
    pools = [bases] * repeat
    res = [[]]
    for pool in pools:
        res = [x + [y] for x in res for y in pool]
    for prod in res:
        yield tuple(prod)


def calculate_kmer_freqs(sequence: str, k: int) -> np.ndarray:
    """Compute relative k-mer frequencies for a nucleotide sequence."""
    bases = ["A", "C", "G", "T"]
    all_kmers = ["".join(p) for p in _product(bases, repeat=k)]
    kmer2idx = {km: i for i, km in enumerate(all_kmers)}

    seq = (sequence or "").upper().replace("U", "T")
    if len(seq) < k:
        return np.zeros((len(all_kmers),), dtype=np.float32)

    counts = np.zeros((len(all_kmers),), dtype=np.float32)
    total = 0
    for i in range(len(seq) - k + 1):
        km = seq[i : i + k]
        idx = kmer2idx.get(km)
        if idx is not None:
            counts[idx] += 1.0
            total += 1

    if total <= 0:
        return counts
    return counts / float(total)


def prepare_kmer_features(
    input_fastas: list[str],
    cache_dir: str,
    model_dir: str,
    kmer_k: int = 6,
) -> pd.DataFrame:
    """Compute plasmid k-mer features and optionally cache them to disk."""
    fastas = validate_fastas(input_fastas)
    rows = []
    for fasta in fastas:
        plasmid_id = get_plasmid_id_from_fasta(fasta)
        try:
            rec = next(SeqIO.parse(fasta, "fasta"))
            seq = str(rec.seq)
        except Exception:
            seq = ""
        feat = calculate_kmer_freqs(seq, kmer_k)
        rows.append({
            "PLASMID": plasmid_id,
            "PLASMID_FASTA": fasta,
            "kmer_features": feat,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Prodigal + codon features (fixed -p meta)
# ---------------------------------------------------------------------------
def check_prodigal_available() -> str:
    """Check whether prodigal is available on PATH."""
    exe = shutil.which("prodigal")
    if exe is None:
        raise PlasRipDependencyError("Required dependency not found on PATH: prodigal")
    return exe


def get_plasmid_cds_path(fasta_path: str, cds_dir: str) -> str:
    base = os.path.splitext(os.path.basename(fasta_path))[0]
    return os.path.join(cds_dir, base + ".cds.fna")


def run_prodigal(in_fna: str, out_cds_fna: str, quiet: bool = True) -> str:
    """Run prodigal for CDS extraction using fixed ``-p meta`` mode."""
    check_prodigal_available()
    os.makedirs(os.path.dirname(out_cds_fna), exist_ok=True)
    cmd = ["prodigal", "-i", in_fna, "-d", out_cds_fna, "-p", "meta"]
    if quiet:
        cmd.append("-q")
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return "meta"


def get_codons_from_seq(seq: str) -> list[str]:
    seq = (seq or "").upper().replace("U", "T")
    if "N" in seq:
        return []
    out: list[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if len(codon) == 3 and codon in CODONS:
            out.append(codon)
    return out


def compute_counts_from_fasta(path: str):
    codon_counts: Counter[str] = Counter()
    total = 0
    for rec in SeqIO.parse(path, "fasta"):
        seq = str(rec.seq).upper()
        cods = get_codons_from_seq(seq)
        if cods:
            codon_counts.update(cods)
            total += len(cods)
    return codon_counts, total


def compute_rscu_from_counts(codon_counts: Counter[str], min_count_threshold: int = 1) -> dict[str, float]:
    rscu: dict[str, float] = {}
    for aa, cods in AA2CODONS.items():
        total = sum(codon_counts.get(c, 0) for c in cods)
        if total < min_count_threshold:
            for c in cods:
                rscu[c] = 0.0
        else:
            mean = total / len(cods)
            for c in cods:
                rscu[c] = (codon_counts.get(c, 0) + 1) / (mean if mean > 0 else 1.0)
    for c in CODONS:
        rscu.setdefault(c, 0.0)
    return rscu


def compute_gc3_from_counts(codon_counts: Counter[str], total_codons: int) -> float:
    if total_codons == 0:
        return np.nan
    third_gc, third_total = 0, 0
    for codon, cnt in codon_counts.items():
        if cnt <= 0:
            continue
        third_total += cnt
        if codon[2] in ("G", "C"):
            third_gc += cnt
    return third_gc / third_total if third_total else np.nan


def _f2(n1, n2):
    n = n1 + n2
    return (n1 * n1 + n2 * n2 - n) / (n * (n - 1)) if n > 1 else np.nan


def _f3(n1, n2, n3):
    n = n1 + n2 + n3
    return (n1 * n1 + n2 * n2 + n3 * n3 - n) / (n * (n - 1)) if n > 1 else np.nan


def _f4(n1, n2, n3, n4):
    n = n1 + n2 + n3 + n4
    return (n1 * n1 + n2 * n2 + n3 * n3 + n4 * n4 - n) / (n * (n - 1)) if n > 1 else np.nan


def compute_enc_wright1990(codon_counts: Counter[str]) -> float:
    aa_codon_counts = {aa: [codon_counts.get(c, 0) for c in AA2CODONS[aa]] for aa in AA2CODONS}
    f2_list, f3_list, f4_list = [], [], []

    for aa, cods in AA2CODONS.items():
        deg = AA_DEGENERACY.get(aa)
        if deg is None or deg == 1:
            continue

        counts = aa_codon_counts[aa]
        if sum(counts) <= 1:
            continue

        if aa in ("L", "S", "R"):
            codon_list = AA2CODONS[aa]
            if aa == "L":
                idx2 = [codon_list.index("TTA"), codon_list.index("TTG")]
                idx4 = [codon_list.index("CTT"), codon_list.index("CTC"), codon_list.index("CTA"), codon_list.index("CTG")]
            elif aa == "S":
                idx2 = [codon_list.index("AGT"), codon_list.index("AGC")]
                idx4 = [codon_list.index("TCT"), codon_list.index("TCC"), codon_list.index("TCA"), codon_list.index("TCG")]
            else:
                idx2 = [codon_list.index("AGA"), codon_list.index("AGG")]
                idx4 = [codon_list.index("CGT"), codon_list.index("CGC"), codon_list.index("CGA"), codon_list.index("CGG")]

            c2 = [counts[i] for i in idx2]
            c4 = [counts[i] for i in idx4]
            f2 = _f2(c2[0], c2[1])
            f4 = _f4(c4[0], c4[1], c4[2], c4[3])
            if not np.isnan(f2):
                f2_list.append(f2)
            if not np.isnan(f4):
                f4_list.append(f4)
            continue

        if deg == 2:
            f = _f2(counts[0], counts[1])
            if not np.isnan(f):
                f2_list.append(f)
        elif deg == 3:
            f = _f3(counts[0], counts[1], counts[2])
            if not np.isnan(f):
                f3_list.append(f)
        elif deg == 4:
            f = _f4(counts[0], counts[1], counts[2], counts[3])
            if not np.isnan(f):
                f4_list.append(f)

    F2 = float(np.mean(f2_list)) if f2_list else np.nan
    F3 = float(np.mean(f3_list)) if f3_list else np.nan
    F4 = float(np.mean(f4_list)) if f4_list else np.nan

    eps = 1e-9
    if np.isnan(F2):
        F2 = eps
    if np.isnan(F3):
        F3 = eps
    if np.isnan(F4):
        F4 = eps

    enc = 2.0 + 9.0 / (F2 + eps) + 1.0 / (F3 + eps) + 5.0 / (F4 + eps)
    return float(max(20.0, min(61.0, enc)))


def compute_codon_feature_vector_from_counts(codon_counts: Counter[str], total_codons: int) -> np.ndarray:
    if codon_counts is None or total_codons <= 0:
        return default_codon_feature_vector()

    rscu_dict = compute_rscu_from_counts(codon_counts)
    gc3 = compute_gc3_from_counts(codon_counts, total_codons)
    enc = compute_enc_wright1990(codon_counts)

    rscu = [rscu_dict.get(c, 0.0) for c in CODONS]
    feat = np.array(rscu + [gc3, enc], dtype=np.float32)
    feat = np.nan_to_num(feat, nan=0.0)
    return feat


def compute_plasmid_codon_feature_from_fasta(fasta_path: str, cds_dir: str) -> np.ndarray:
    """Extract plasmid codon features, using fixed ``-p meta`` for prodigal."""
    os.makedirs(cds_dir, exist_ok=True)
    cds_path = get_plasmid_cds_path(fasta_path, cds_dir)

    if not (os.path.exists(cds_path) and os.path.getsize(cds_path) > 0):
        try:
            run_prodigal(fasta_path, cds_path, quiet=True)
            logging.info("[Prodigal] success mode=meta input=%s", os.path.basename(fasta_path))
        except Exception as exc:
            logging.warning(
                "[Prodigal] failed mode=meta input=%s; using default codon feature. error=%s",
                fasta_path,
                exc,
            )
            return default_codon_feature_vector()
    else:
        logging.info("[Prodigal] reuse cached CDS: %s", cds_path)

    try:
        counts, total = compute_counts_from_fasta(cds_path)
        if counts is None or total <= 0:
            logging.warning("[CodonFeat] empty/invalid CDS counts for %s; using default feature.", fasta_path)
            return default_codon_feature_vector()
        return compute_codon_feature_vector_from_counts(counts, total)
    except Exception as exc:
        logging.warning(
            "[CodonFeat] failed to compute codon feature for %s; using default feature. error=%s",
            fasta_path,
            exc,
        )
        return default_codon_feature_vector()


def prepare_orf_features(
    input_fastas: list[str],
    cache_dir: str,
    threads: int = 1,
):
    """Prepare CDS files for plasmids using prodigal ``-p meta``.

    Returns a list of per-plasmid records documenting the cache paths and
    whether CDS extraction succeeded.
    """
    fastas = validate_fastas(input_fastas)
    cds_dir = Path(cache_dir).resolve() / "plasmid_prodigal_orfs"
    cds_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, Any]] = []
    for fasta in fastas:
        cds_path = get_plasmid_cds_path(fasta, str(cds_dir))
        success = False
        used_mode = "meta"

        if os.path.exists(cds_path) and os.path.getsize(cds_path) > 0:
            success = True
            used_mode = "cached_meta"
        else:
            try:
                run_prodigal(fasta, cds_path, quiet=True)
                success = True
            except Exception as exc:
                logging.warning("[Prodigal] ORF preparation failed for %s: %s", fasta, exc)
                success = False

        records.append(
            {
                "PLASMID": get_plasmid_id_from_fasta(fasta),
                "PLASMID_FASTA": fasta,
                "CDS_FASTA": cds_path,
                "PRODIGAL_MODE": used_mode,
                "PRODIGAL_SUCCESS": success,
            }
        )

    return records


# ---------------------------------------------------------------------------
# Feature preprocessing / scaling
# ---------------------------------------------------------------------------
def apply_plasmid_codon_scaler(feature: np.ndarray, scaler: dict[str, Any]) -> np.ndarray:
    """Apply the saved training-time GC3/ENC scaler to a plasmid feature vector."""
    feat = np.asarray(feature, dtype=np.float32).copy()
    gc3_idx = len(CODONS)
    enc_idx = len(CODONS) + 1

    feat[gc3_idx] = (feat[gc3_idx] - scaler["gc3_mean"]) / (scaler["gc3_std"] + 1e-9)
    feat[enc_idx] = (feat[enc_idx] - scaler["enc_mean"]) / (scaler["enc_std"] + 1e-9)
    feat = np.nan_to_num(feat, nan=0.0)
    return feat.astype(np.float32)


def prepare_codon_features(
    input_fastas: list[str],
    cache_dir: str,
    model_dir: str,
) -> pd.DataFrame:
    """Compute plasmid codon features and apply the saved scaler."""
    fastas = validate_fastas(input_fastas)
    preprocessors = load_feature_preprocessors(model_dir)
    scaler = preprocessors["plasmid_codon_scaler"]

    cds_dir = Path(cache_dir).resolve() / "plasmid_prodigal_orfs"
    cds_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for fasta in fastas:
        plasmid_id = get_plasmid_id_from_fasta(fasta)
        feat = compute_plasmid_codon_feature_from_fasta(fasta, str(cds_dir))
        feat = apply_plasmid_codon_scaler(feat, scaler)
        rows.append(
            {
                "PLASMID": plasmid_id,
                "PLASMID_FASTA": fasta,
                "codon_features": feat,
            }
        )
    return pd.DataFrame(rows)


def build_feature_index(plasmid_ids: list[str]) -> pd.DataFrame:
    """Build a simple feature index dataframe."""
    return pd.DataFrame({"PLASMID": plasmid_ids})


# ---------------------------------------------------------------------------
# Top-level extraction
# ---------------------------------------------------------------------------
def extract_plasmid_features(
    input_fastas: list[str],
    model_dir: str,
    cache_dir: str | None = None,
    threads: int = 1,
    keep_intermediate: bool = False,
) -> dict[str, Any]:
    """Extract all model-ready plasmid features for host-range prediction."""
    fastas = validate_fastas(input_fastas)
    cache_root = Path(cache_dir).resolve() if cache_dir else None
    if cache_root is not None:
        cache_root.mkdir(parents=True, exist_ok=True)

    kmer = prepare_kmer_features(
        input_fastas=fastas,
        cache_dir=str(cache_root) if cache_root else str(Path.cwd()),
        model_dir=model_dir,
        kmer_k=6,
    )
    codon = prepare_codon_features(fastas, str(cache_root) if cache_root else str(Path.cwd()), model_dir)

    merged = kmer.merge(codon[["PLASMID", "codon_features"]], on="PLASMID", how="inner")
    plasmid_ids = merged["PLASMID"].tolist()

    feature_index = build_feature_index(plasmid_ids)
    return {
        "input_fastas": fastas,
        "plasmid_ids": plasmid_ids,
        "model_dir": str(Path(model_dir).resolve()),
        "cache_dir": str(cache_root) if cache_root else None,
        "keep_intermediate": keep_intermediate,
        "threads": threads,
        "feature_index": feature_index,
        "kmer": kmer,
        "codon": codon,
        "feature_table": merged,
        "kmer_features": np.stack(merged["kmer_features"].tolist(), axis=0).astype(np.float32) if len(merged) else np.empty((0, 4 ** 6), dtype=np.float32),
        "codon_features": np.stack(merged["codon_features"].tolist(), axis=0).astype(np.float32) if len(merged) else np.empty((0, len(CODONS) + 2), dtype=np.float32),
    }
