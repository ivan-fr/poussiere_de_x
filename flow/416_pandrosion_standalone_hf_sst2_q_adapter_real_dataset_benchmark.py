#!/usr/bin/env python3
"""SST-2 real-dataset benchmark for Pandrosion q-adapters.

416 implements the next research steps after 415:

  1. Real open-source dataset: nyu-mll/glue, config sst2.
  2. Equal-parameter comparisons against LoRA and dense adapters.
  3. Multi-seed runs.
  4. Two model sizes: BERT hidden 128 and hidden 256.
  5. JSON, CSV, and publication-style benchmark curves.

The benchmark trains only adapter/classifier parameters on frozen BERT features
for speed, then measures end-to-end BERT+adapter latency separately.  This is a
frozen-backbone adapter benchmark, not a full fine-tuning benchmark.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
import random
import sys
import time
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from datasets import load_dataset
from transformers import AutoModel, AutoTokenizer


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "flow" / "412_pandrosion_standalone_local_jet_geometry_composite_cached_irp_engine.py"
OUT_JSON = Path("/private/tmp/416_hf_sst2_q_adapter_real_dataset_benchmark.json")
OUT_CSV = Path("/private/tmp/416_hf_sst2_q_adapter_real_dataset_rows.csv")
FIG_DIR = Path("/private/tmp/416_hf_sst2_q_adapter_figures")
DATASET_REPO = "nyu-mll/glue"
DATASET_CONFIG = "sst2"
MODEL_IDS = [
    "google/bert_uncased_L-2_H-128_A-2",
    "google/bert_uncased_L-4_H-256_A-4",
]
FP32_ADAPTER_ACCEPT = 1e-6


def sync() -> None:
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    if torch.backends.mps.is_available():
        torch.mps.synchronize()


def bench(fn: Any, repeats: int = 12, warmup: int = 4) -> float:
    for _ in range(max(0, warmup)):
        fn()
    sync()
    vals: list[float] = []
    for _ in range(max(1, repeats)):
        sync()
        t0 = time.perf_counter()
        fn()
        sync()
        vals.append(time.perf_counter() - t0)
    return float(min(vals))


def selected_device() -> Any:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def set_seed(seed: int) -> None:
    random.seed(int(seed))
    np.random.seed(int(seed) & 0xFFFFFFFF)
    torch.manual_seed(int(seed) & 0x7FFFFFFF)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(int(seed) & 0x7FFFFFFF)


def load_engine() -> Any:
    spec = importlib.util.spec_from_file_location("pandrosion412_sst2", ENGINE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {ENGINE}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def configure(mod: Any) -> None:
    mod.configure_matrix_backend(
        argparse.Namespace(
            matrix_backend="torch",
            matrix_algorithm="adaptive-directional-sketch",
            batch_kernel="auto",
            torch_device="auto",
            torch_complex_dtype="auto",
            torch_real_dtype="auto",
            ns_iters=8,
            ns_damping_floor=1e-5,
            sketch_dim=0,
            sketch_factor=2.75,
            sketch_min_n=1,
            sketch_seed=412,
            sketch_mode="cached-rademacher",
            sketch_solver="svd",
            sketch_basis_cache=True,
            sketch_basis_cache_max=8192,
            directional_jet=True,
            directional_jet_min_n=1,
            directional_jet_factor=2.75,
            directional_diff="auto",
            directional_auto_central_cap=0.35,
            directional_coded_probe=True,
            directional_coded_factor=2.0,
            directional_coded_min=8,
            directional_coded_max=0,
            directional_fast_projector=True,
            directional_fast_projector_cap=0.05,
            directional_fast_projector_diagnostics=False,
            directional_basis_reuse=True,
        )
    )


class DenseAdapter(torch.nn.Module):
    def __init__(self, hidden_size: int) -> None:
        super().__init__()
        self.linear = torch.nn.Linear(int(hidden_size), int(hidden_size), bias=False)

    def forward(self, x: Any) -> Any:
        return self.linear(x)


class RandomLoRAAdapter(torch.nn.Module):
    def __init__(self, hidden_size: int, rank: int) -> None:
        super().__init__()
        self.down = torch.nn.Linear(int(hidden_size), int(rank), bias=False)
        self.up = torch.nn.Linear(int(rank), int(hidden_size), bias=False)

    def forward(self, x: Any) -> Any:
        return self.up(self.down(x))


class PandrosionQAdapter(torch.nn.Module):
    def __init__(self, hidden_size: int, q_basis_np: Any, rank: int) -> None:
        super().__init__()
        C = np.asarray(q_basis_np, dtype=np.float32)[:, : int(rank)]
        if C.ndim != 2 or int(C.shape[0]) != int(hidden_size):
            raise ValueError("bad q basis shape")
        self.down = torch.nn.Linear(int(hidden_size), int(rank), bias=False)
        self.register_buffer("fixed_q_basis_t", torch.as_tensor(C.T.copy(), dtype=torch.float32))
        self.rank = int(rank)

    def forward(self, x: Any) -> Any:
        return self.down(x) @ self.fixed_q_basis_t


class FeatureClassifier(torch.nn.Module):
    def __init__(self, hidden_size: int, adapter: torch.nn.Module | None, scope: str, adapter_scale: float = 0.35) -> None:
        super().__init__()
        self.adapter = adapter
        self.scope = str(scope)
        self.adapter_scale = float(adapter_scale)
        self.classifier = torch.nn.Linear(int(hidden_size), 2)

    def pooled(self, hidden: Any, attention_mask: Any) -> tuple[Any, Any]:
        adapter_delta: Any | None = None
        if self.adapter is None:
            pooled = hidden[:, 0, :]
            adapter_delta = pooled
        elif self.scope == "all_tokens":
            adapter_delta = self.adapter(hidden)
            hidden = hidden + self.adapter_scale * adapter_delta
            mask = attention_mask.to(dtype=hidden.dtype)[:, :, None]
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp_min(1.0)
        else:
            pooled = hidden[:, 0, :]
            adapter_delta = self.adapter(pooled)
            pooled = pooled + self.adapter_scale * adapter_delta
        return pooled, adapter_delta

    def forward(self, hidden: Any, attention_mask: Any, return_adapter: bool = False) -> Any:
        pooled, adapter_delta = self.pooled(hidden, attention_mask)
        logits = self.classifier(pooled)
        if return_adapter:
            return logits, adapter_delta
        return logits


class EndToEndClassifier(torch.nn.Module):
    def __init__(self, bert: torch.nn.Module, feature_model: FeatureClassifier) -> None:
        super().__init__()
        self.bert = bert
        self.feature_model = feature_model
        for p in self.bert.parameters():
            p.requires_grad_(False)
        self.bert.eval()

    def forward(self, input_ids: Any, attention_mask: Any, token_type_ids: Any | None = None) -> Any:
        with torch.no_grad():
            kwargs = {"input_ids": input_ids, "attention_mask": attention_mask}
            if token_type_ids is not None:
                kwargs["token_type_ids"] = token_type_ids
            hidden = self.bert(**kwargs).last_hidden_state.detach()
        return self.feature_model(hidden, attention_mask)


class VectorTarget:
    def __init__(self, root: Any) -> None:
        self.root = np.asarray(root, dtype=np.complex128)
        self.evals = 0

    def eval(self, z: Any) -> Any:
        self.evals += 1
        return np.asarray(z, dtype=np.complex128) - self.root

    def eval_batch(self, Z: Any) -> Any:
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        self.evals += int(ZZ.shape[0])
        return ZZ - self.root[None, :]

    def residual(self, z: Any) -> float:
        return float(np.linalg.norm(self.eval(z)))

    def residuals_batch(self, Z: Any) -> Any:
        return np.linalg.norm(self.eval_batch(Z), axis=1)


def normalize_vector(x: Any, scale: float = 0.25) -> Any:
    v = np.asarray(x, dtype=np.float64).reshape(-1)
    nrm = float(np.linalg.norm(v))
    if not math.isfinite(nrm) or nrm <= 0.0:
        raise ValueError("empty_or_zero_vector")
    return (float(scale) * v / nrm).astype(np.complex128)


def q_dimensions(mod: Any, n: int) -> tuple[int, int]:
    k = max(1, min(int(n), int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))))
    q = int(mod.directional_coded_dim_for(int(n), int(k)))
    return int(k), int(q)


def q_salt(seed: int) -> int:
    return int(seed) + 0x412D1A


def q_basis(mod: Any, n: int, seed: int) -> tuple[Any, int, int]:
    k, q = q_dimensions(mod, int(n))
    C = np.asarray(mod._coded_composite_basis_np(int(n), k, q, salt=q_salt(seed)), dtype=np.complex128)
    return C.real.astype(np.float64, copy=False), k, q


def subspace_energy(mod: Any, root: Any, seed: int) -> dict[str, Any]:
    r = np.asarray(root, dtype=np.complex128).reshape(-1)
    n = int(r.size)
    C, k, q = q_basis(mod, n, seed)
    P = np.asarray(mod._sketch_basis_np(n, k, salt=q_salt(seed)), dtype=np.complex128)
    Cc = C.astype(np.complex128, copy=False)
    den = max(1e-300, float(np.vdot(r, r).real))
    q_coeff = Cc.conj().T @ r
    p_coeff = P.conj().T @ r
    return {
        "n": n,
        "k": k,
        "q": q,
        "q_energy_fraction": float(np.vdot(q_coeff, q_coeff).real / den),
        "full_sketch_energy_fraction": float(np.vdot(p_coeff, p_coeff).real / den),
        "expected_random_q_fraction": float(q / max(1, n)),
        "expected_random_full_sketch_fraction": float(k / max(1, n)),
    }


def run_corrector(mod: Any, root: Any, seed: int) -> dict[str, Any]:
    root_np = np.asarray(root, dtype=np.complex128)
    target = VectorTarget(root_np)
    y0 = np.zeros_like(root_np)
    t0 = time.perf_counter()
    rec = mod.hypercube_matrixjet_corrector(
        target,
        y0,
        max_epochs=2,
        tol=1e-8,
        accept=FP32_ADAPTER_ACCEPT,
        trial_timeout=0.0,
        line_search=6,
        line_grid=(1.0, 0.5, 0.25, 0.125),
        direction_seed=int(seed),
        cloud_nodes=0,
        lm_damping=0.0,
        trust_radius=0.0,
        matrix_rcond=1e-12,
        matrix_condition_cap=1e12,
    )
    ms = 1e3 * (time.perf_counter() - t0)
    y = np.asarray(rec.get("y", y0), dtype=np.complex128)
    return {
        "corrector_ms": float(ms),
        "status": rec.get("status"),
        "accepted": bool(rec.get("accepted")),
        "residual": float(rec.get("residual", float("inf"))),
        "root_relative_error": float(np.linalg.norm(y - root_np) / max(1e-300, float(np.linalg.norm(root_np)))),
        "probe_kind": rec.get("directional_probe_kind"),
        "nodes": int(rec.get("hypercube_nodes", 0) or 0),
        "coded_dim": rec.get("directional_coded_dim"),
        "parent_k": rec.get("directional_parent_sketch_dim"),
        "linear_relative_residual": rec.get("directional_linear_relative_residual"),
        "q_active": bool(
            rec.get("directional_probe_kind") == "coded-probe"
            and rec.get("hypercube_nodes") == rec.get("directional_coded_dim")
        ),
        **subspace_energy(mod, root_np, seed),
    }


def trainable_params(module: torch.nn.Module) -> int:
    return int(sum(p.numel() for p in module.parameters() if p.requires_grad))


def buffer_params(module: torch.nn.Module) -> int:
    return int(sum(b.numel() for b in module.buffers()))


def dataset_subset(limit_train: int, limit_val: int, seed: int) -> tuple[list[str], list[int], list[str], list[int]]:
    ds = load_dataset(DATASET_REPO, DATASET_CONFIG)
    train = ds["train"].shuffle(seed=int(seed)).select(range(min(int(limit_train), len(ds["train"]))))
    val = ds["validation"].shuffle(seed=int(seed) + 17).select(range(min(int(limit_val), len(ds["validation"]))))
    train_text = [str(x) for x in train["sentence"]]
    train_y = [int(x) for x in train["label"]]
    val_text = [str(x) for x in val["sentence"]]
    val_y = [int(x) for x in val["label"]]
    return train_text, train_y, val_text, val_y


def tokenize(tokenizer: Any, texts: list[str], labels: list[int], device: Any) -> dict[str, Any]:
    enc = tokenizer(texts, return_tensors="pt", padding=True, truncation=True, max_length=64)
    out = {k: v.to(device) for k, v in enc.items()}
    out["labels"] = torch.as_tensor(labels, dtype=torch.long, device=device)
    return out


def precompute_hidden(bert: torch.nn.Module, data: dict[str, Any], batch_size: int) -> dict[str, Any]:
    hidden_parts: list[Any] = []
    n = int(data["labels"].shape[0])
    for start in range(0, n, int(batch_size)):
        end = min(n, start + int(batch_size))
        kwargs = {
            "input_ids": data["input_ids"][start:end],
            "attention_mask": data["attention_mask"][start:end],
        }
        if "token_type_ids" in data:
            kwargs["token_type_ids"] = data["token_type_ids"][start:end]
        with torch.no_grad():
            hidden_parts.append(bert(**kwargs).last_hidden_state.detach().cpu())
    return {
        "hidden": torch.cat(hidden_parts, dim=0),
        "attention_mask": data["attention_mask"].detach().cpu(),
        "labels": data["labels"].detach().cpu(),
        "input_ids": data["input_ids"].detach().cpu(),
        **({"token_type_ids": data["token_type_ids"].detach().cpu()} if "token_type_ids" in data else {}),
    }


def to_device_batch(batch: dict[str, Any], device: Any) -> dict[str, Any]:
    return {k: v.to(device) for k, v in batch.items()}


def iter_feature_batches(data: dict[str, Any], batch_size: int, seed: int, shuffle: bool, device: Any) -> list[dict[str, Any]]:
    n = int(data["labels"].shape[0])
    idx = np.arange(n)
    if shuffle:
        rng = np.random.default_rng(int(seed))
        rng.shuffle(idx)
    batches: list[dict[str, Any]] = []
    for start in range(0, n, int(batch_size)):
        sel_np = idx[start : start + int(batch_size)]
        batch = {
            "hidden": data["hidden"][sel_np],
            "attention_mask": data["attention_mask"][sel_np],
            "labels": data["labels"][sel_np],
        }
        batches.append(to_device_batch(batch, device))
    return batches


def evaluate(model: FeatureClassifier, data: dict[str, Any], batch_size: int, device: Any) -> dict[str, float]:
    model.eval()
    loss_fn = torch.nn.CrossEntropyLoss(reduction="sum")
    total_loss = 0.0
    correct = 0
    total = 0
    for batch in iter_feature_batches(data, batch_size, 0, False, device):
        with torch.no_grad():
            logits = model(batch["hidden"], batch["attention_mask"])
            loss = loss_fn(logits, batch["labels"])
            pred = torch.argmax(logits, dim=-1)
        total_loss += float(loss.detach().cpu().item())
        correct += int((pred == batch["labels"]).sum().detach().cpu().item())
        total += int(batch["labels"].numel())
    return {"loss": total_loss / max(1, total), "accuracy": correct / max(1, total)}


def train_model(model: FeatureClassifier, train_data: dict[str, Any], val_data: dict[str, Any], seed: int, epochs: int, batch_size: int, device: Any) -> dict[str, Any]:
    set_seed(seed)
    model.train()
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad], lr=3e-3, weight_decay=1e-3)
    loss_fn = torch.nn.CrossEntropyLoss()
    history: list[dict[str, float]] = []
    t0 = time.perf_counter()
    initial = evaluate(model, val_data, batch_size, device)
    for ep in range(int(epochs)):
        total = 0.0
        count = 0
        for batch in iter_feature_batches(train_data, batch_size, seed + 1009 * ep, True, device):
            model.train()
            opt.zero_grad(set_to_none=True)
            logits = model(batch["hidden"], batch["attention_mask"])
            loss = loss_fn(logits, batch["labels"])
            loss.backward()
            torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad], 5.0)
            opt.step()
            total += float(loss.detach().cpu().item()) * int(batch["labels"].numel())
            count += int(batch["labels"].numel())
        val = evaluate(model, val_data, batch_size, device)
        history.append({"epoch": float(ep + 1), "train_loss": total / max(1, count), "val_loss": val["loss"], "val_accuracy": val["accuracy"]})
    final = evaluate(model, val_data, batch_size, device)
    return {
        "epochs": int(epochs),
        "batch_size": int(batch_size),
        "initial_val_loss": initial["loss"],
        "initial_val_accuracy": initial["accuracy"],
        "final_val_loss": final["loss"],
        "final_val_accuracy": final["accuracy"],
        "train_seconds": float(time.perf_counter() - t0),
        "history": history,
    }


def make_adapter(kind: str, hidden: int, C: Any, q_full: int, rank: int) -> tuple[torch.nn.Module | None, int]:
    if kind == "none":
        return None, 0
    if kind == "dense":
        return DenseAdapter(hidden), hidden
    if kind == "lora":
        return RandomLoRAAdapter(hidden, int(rank)), int(rank)
    if kind == "pandrosion_q":
        return PandrosionQAdapter(hidden, C, int(rank)), int(rank)
    raise ValueError(f"unknown adapter kind {kind!r}")


def validation_roots(model: FeatureClassifier, val_data: dict[str, Any], device: Any, max_roots: int = 3) -> list[Any]:
    batch = {
        "hidden": val_data["hidden"][:max_roots].to(device),
        "attention_mask": val_data["attention_mask"][:max_roots].to(device),
    }
    model.eval()
    with torch.no_grad():
        _, delta = model(batch["hidden"], batch["attention_mask"], return_adapter=True)
    arr = delta.detach().cpu().numpy()
    roots: list[Any] = []
    if arr.ndim == 2:
        for i in range(arr.shape[0]):
            roots.append(normalize_vector(arr[i]))
    elif arr.ndim == 3:
        mask = batch["attention_mask"].detach().cpu().numpy()
        for i in range(arr.shape[0]):
            token_ids = np.where(mask[i] > 0)[0]
            for j in token_ids[:1]:
                roots.append(normalize_vector(arr[i, j]))
                if len(roots) >= max_roots:
                    return roots
    return roots[:max_roots]


def validate_q_path(mod: Any, model: FeatureClassifier, val_data: dict[str, Any], seed: int, device: Any) -> dict[str, Any]:
    roots = validation_roots(model, val_data, device, max_roots=3)
    recs = [run_corrector(mod, root, seed) for root in roots]
    if not recs:
        return {"q_validation_roots": 0}
    return {
        "q_validation_roots": len(recs),
        "q_active_rate": float(np.mean([bool(r["q_active"]) for r in recs])),
        "q_accepted_rate": float(np.mean([bool(r["accepted"]) for r in recs])),
        "q_energy_mean": float(np.mean([float(r["q_energy_fraction"]) for r in recs])),
        "q_energy_min": float(np.min([float(r["q_energy_fraction"]) for r in recs])),
        "q_root_error_max": float(np.max([float(r["root_relative_error"]) for r in recs])),
        "q_residual_max": float(np.max([float(r["residual"]) for r in recs])),
        "q_corrector_ms_mean": float(np.mean([float(r["corrector_ms"]) for r in recs])),
        "q_probe_k": int(recs[0]["k"]),
        "q_probe_q": int(recs[0]["q"]),
    }


def end_to_end_latency(bert: torch.nn.Module, feature_model: FeatureClassifier, tokenized_val: dict[str, Any], device: Any, batch_size: int) -> float:
    e2e = EndToEndClassifier(bert, feature_model).to(device=device, dtype=torch.float32).eval()
    n = min(int(batch_size), int(tokenized_val["labels"].shape[0]))
    batch = {k: v[:n].to(device) for k, v in tokenized_val.items() if k != "labels"}

    def once() -> Any:
        with torch.no_grad():
            return e2e(**batch)

    return 1e3 * bench(once, repeats=16, warmup=4)


def stats(vals: list[float]) -> dict[str, float]:
    a = np.asarray(vals, dtype=float)
    return {
        "mean": float(np.mean(a)),
        "std": float(np.std(a, ddof=1)) if a.size > 1 else 0.0,
        "min": float(np.min(a)),
        "median": float(np.median(a)),
        "max": float(np.max(a)),
    }


def summarize(rows: list[dict[str, Any]], keys: list[str]) -> dict[str, Any]:
    groups = sorted({tuple(str(r[k]) for k in keys) for r in rows})
    out: dict[str, Any] = {}
    for group in groups:
        label = " / ".join(group)
        rs = [r for r in rows if tuple(str(r[k]) for k in keys) == group]
        out[label] = {
            "rows": len(rs),
            "final_val_accuracy": stats([float(r["final_val_accuracy"]) for r in rs]),
            "final_val_loss": stats([float(r["final_val_loss"]) for r in rs]),
            "forward_ms_end_to_end": stats([float(r["forward_ms_end_to_end"]) for r in rs]),
            "trainable_params": stats([float(r["trainable_params"]) for r in rs]),
            "adapter_trainable_params": stats([float(r["adapter_trainable_params"]) for r in rs]),
            "q_active_rate": stats([float(r.get("q_active_rate", 0.0)) for r in rs]),
            "q_accepted_rate": stats([float(r.get("q_accepted_rate", 0.0)) for r in rs]),
            "q_energy_mean": stats([float(r.get("q_energy_mean", 0.0)) for r in rs]),
            "q_root_error_max": stats([float(r.get("q_root_error_max", 0.0)) for r in rs]),
        }
    return out


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "model_id",
        "hidden_size",
        "seed",
        "variant",
        "kind",
        "scope",
        "ablation",
        "adapter_rank",
        "engine_q_dim",
        "equal_lora_rank",
        "trainable_params",
        "adapter_trainable_params",
        "final_val_accuracy",
        "final_val_loss",
        "forward_ms_end_to_end",
        "q_active_rate",
        "q_accepted_rate",
        "q_energy_mean",
        "q_root_error_max",
        "q_corrector_ms_mean",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def plot_figures(rows: list[dict[str, Any]], fig_dir: Path) -> dict[str, str]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    figs: dict[str, str] = {}
    colors = {
        "none": "#777777",
        "dense": "#1f77b4",
        "lora": "#ff7f0e",
        "pandrosion_q": "#2ca02c",
    }

    plt.figure(figsize=(8, 5))
    for r in rows:
        marker = "o" if int(r["hidden_size"]) == 128 else "s"
        plt.scatter(
            float(r["trainable_params"]),
            float(r["final_val_accuracy"]),
            c=colors.get(str(r["kind"]), "#333333"),
            marker=marker,
            alpha=0.78,
            s=54,
        )
    plt.xscale("log")
    plt.xlabel("Trainable parameters (log scale)")
    plt.ylabel("SST-2 validation accuracy")
    plt.title("416 SST-2 accuracy vs parameters")
    plt.grid(True, alpha=0.25)
    path = fig_dir / "416_accuracy_vs_params.png"
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.savefig(path.with_suffix(".pdf"))
    plt.close()
    figs["accuracy_vs_params_png"] = str(path)
    figs["accuracy_vs_params_pdf"] = str(path.with_suffix(".pdf"))

    labels = sorted({str(r["variant"]) for r in rows})
    means = [float(np.mean([float(r.get("q_active_rate", 0.0)) for r in rows if str(r["variant"]) == lab])) for lab in labels]
    plt.figure(figsize=(10, 5))
    plt.bar(range(len(labels)), means, color="#2ca02c")
    plt.xticks(range(len(labels)), labels, rotation=35, ha="right")
    plt.ylim(0.0, 1.05)
    plt.ylabel("Mean q-active rate")
    plt.title("416 q-path activation by adapter")
    plt.grid(True, axis="y", alpha=0.25)
    path = fig_dir / "416_q_active_by_variant.png"
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.savefig(path.with_suffix(".pdf"))
    plt.close()
    figs["q_active_by_variant_png"] = str(path)
    figs["q_active_by_variant_pdf"] = str(path.with_suffix(".pdf"))

    plt.figure(figsize=(8, 5))
    for kind in sorted({str(r["kind"]) for r in rows}):
        rs = [r for r in rows if str(r["kind"]) == kind]
        plt.scatter(
            [float(r["forward_ms_end_to_end"]) for r in rs],
            [float(r["final_val_accuracy"]) for r in rs],
            label=kind,
            c=colors.get(kind, "#333333"),
            alpha=0.78,
            s=54,
        )
    plt.xlabel("End-to-end forward latency, ms")
    plt.ylabel("SST-2 validation accuracy")
    plt.title("416 latency vs accuracy")
    plt.grid(True, alpha=0.25)
    plt.legend()
    path = fig_dir / "416_latency_vs_accuracy.png"
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.savefig(path.with_suffix(".pdf"))
    plt.close()
    figs["latency_vs_accuracy_png"] = str(path)
    figs["latency_vs_accuracy_pdf"] = str(path.with_suffix(".pdf"))

    plt.figure(figsize=(9, 5))
    for variant in sorted({str(r["variant"]) for r in rows}):
        histories = [r["history"] for r in rows if str(r["variant"]) == variant]
        if not histories:
            continue
        max_ep = max(len(h) for h in histories)
        xs = list(range(1, max_ep + 1))
        ys = []
        for ep in range(max_ep):
            vals = [float(h[ep]["val_loss"]) for h in histories if len(h) > ep]
            ys.append(float(np.mean(vals)))
        plt.plot(xs, ys, label=variant, linewidth=1.4)
    plt.xlabel("Epoch")
    plt.ylabel("Mean validation loss")
    plt.title("416 validation loss curves")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=7)
    path = fig_dir / "416_val_loss_curves.png"
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.savefig(path.with_suffix(".pdf"))
    plt.close()
    figs["val_loss_curves_png"] = str(path)
    figs["val_loss_curves_pdf"] = str(path.with_suffix(".pdf"))
    return figs


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--train-limit", type=int, default=512)
    p.add_argument("--val-limit", type=int, default=256)
    p.add_argument("--epochs", type=int, default=4)
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--seeds", type=str, default="41601,41602,41603,41604,41605")
    p.add_argument("--models", type=str, default=",".join(MODEL_IDS))
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    mod = load_engine()
    configure(mod)
    device = selected_device()
    seeds = [int(x.strip()) for x in str(args.seeds).split(",") if x.strip()]
    model_ids = [x.strip() for x in str(args.models).split(",") if x.strip()]
    train_text, train_y, val_text, val_y = dataset_subset(int(args.train_limit), int(args.val_limit), 416)
    rows: list[dict[str, Any]] = []
    model_meta: list[dict[str, Any]] = []

    for model_id in model_ids:
        tokenizer = AutoTokenizer.from_pretrained(model_id, local_files_only=True)
        bert = AutoModel.from_pretrained(model_id, local_files_only=True).to(device=device, dtype=torch.float32).eval()
        for p in bert.parameters():
            p.requires_grad_(False)
        hidden = int(bert.config.hidden_size)
        train_tok = tokenize(tokenizer, train_text, train_y, device)
        val_tok = tokenize(tokenizer, val_text, val_y, device)
        train_features = precompute_hidden(bert, train_tok, int(args.batch_size))
        val_features = precompute_hidden(bert, val_tok, int(args.batch_size))
        _, parent_k0, q_full0 = q_basis(mod, hidden, 4160000 + hidden)
        equal_lora_rank = max(1, int(math.ceil(q_full0 / 2.0)))
        variants = [
            {"variant": "bert_head_only", "kind": "none", "scope": "cls", "rank": 0, "ablation": "baseline"},
            {"variant": "dense_adapter_cls", "kind": "dense", "scope": "cls", "rank": hidden, "ablation": "dense"},
            {"variant": "lora_equal_params_cls", "kind": "lora", "scope": "cls", "rank": equal_lora_rank, "ablation": "lora_equal_trainable"},
            {"variant": "pandrosion_q_cls_rank_half", "kind": "pandrosion_q", "scope": "cls", "rank": equal_lora_rank, "ablation": "rank_half"},
            {"variant": "pandrosion_q_cls_rank_full", "kind": "pandrosion_q", "scope": "cls", "rank": q_full0, "ablation": "rank_full"},
            {"variant": "pandrosion_q_all_tokens_rank_full", "kind": "pandrosion_q", "scope": "all_tokens", "rank": q_full0, "ablation": "all_tokens"},
        ]
        model_meta.append({"model_id": model_id, "hidden_size": hidden, "parent_k": parent_k0, "q_full": q_full0, "equal_lora_rank": equal_lora_rank})
        for seed in seeds:
            C, parent_k, q_dim = q_basis(mod, hidden, seed)
            for spec in variants:
                set_seed(seed)
                adapter, adapter_rank = make_adapter(spec["kind"], hidden, C, q_dim, int(spec["rank"]))
                feature_model = FeatureClassifier(hidden, adapter, spec["scope"]).to(device=device, dtype=torch.float32)
                train = train_model(feature_model, train_features, val_features, seed, int(args.epochs), int(args.batch_size), device)
                qval = validate_q_path(mod, feature_model, val_features, seed, device)
                fwd_ms = end_to_end_latency(bert, feature_model, val_tok, device, int(args.batch_size))
                row = {
                    "model_id": model_id,
                    "hidden_size": hidden,
                    "seed": seed,
                    "variant": spec["variant"],
                    "kind": spec["kind"],
                    "scope": spec["scope"],
                    "ablation": spec["ablation"],
                    "adapter_rank": adapter_rank,
                    "parent_k": parent_k,
                    "engine_q_dim": q_dim,
                    "equal_lora_rank": equal_lora_rank,
                    "trainable_params": trainable_params(feature_model),
                    "adapter_trainable_params": trainable_params(adapter) if adapter is not None else 0,
                    "buffer_params": buffer_params(feature_model),
                    "forward_ms_end_to_end": fwd_ms,
                    **train,
                    **qval,
                }
                rows.append(row)
                print(
                    f"{Path(model_id).name} {row['variant']} seed={seed} "
                    f"acc={row['final_val_accuracy']:.3f} loss={row['final_val_loss']:.3f} "
                    f"fwd={row['forward_ms_end_to_end']:.3f}ms params={row['trainable_params']} "
                    f"q_rate={row.get('q_active_rate', 0.0):.2f}"
                )

    write_csv(rows, OUT_CSV)
    figures = plot_figures(rows, FIG_DIR)
    result = {
        "engine": str(ENGINE),
        "dataset_repo": DATASET_REPO,
        "dataset_config": DATASET_CONFIG,
        "output_json": str(OUT_JSON),
        "output_csv": str(OUT_CSV),
        "figures": figures,
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "torch": torch.__version__,
        "transformers": __import__("transformers").__version__,
        "datasets": __import__("datasets").__version__,
        "device": str(device),
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "train_examples": len(train_text),
        "val_examples": len(val_text),
        "epochs": int(args.epochs),
        "batch_size": int(args.batch_size),
        "seeds": seeds,
        "models": model_meta,
        "rows": rows,
        "summary_by_model_variant": summarize(rows, ["model_id", "variant"]),
        "summary_by_variant": summarize(rows, ["variant"]),
        "summary_by_kind": summarize(rows, ["kind"]),
        "interpretation": {
            "training": "Frozen BERT hidden states are precomputed for adapter/classifier training; BERT weights are not fine-tuned.",
            "latency": "forward_ms_end_to_end is measured through the actual BERT model plus the trained adapter/classifier.",
            "equal_params": "LoRA rank is ceil(q/2), giving approximately the same adapter trainable parameters as full Pandrosion q.",
            "q_validation": "q_active_rate validates held-out adapter outputs with the standalone 412 Pandrosion IRP corrector.",
            "scope": "CLS scope adapts only pooled CLS; all_tokens adapts token states and uses masked mean pooling.",
        },
    }
    OUT_JSON.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"json={OUT_JSON}")
    print(f"csv={OUT_CSV}")
    for key, value in figures.items():
        print(f"{key}={value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
