#!/usr/bin/env python3
"""End-to-end Hugging Face adapter benchmark for Pandrosion IRP.

415 covers the next five research steps:

  1. Put adapters inside a BERT forward path.
  2. Measure a small text-classification task: loss, accuracy, latency.
  3. Compare against a LoRA adapter with approximately equal trainable params.
  4. Repeat over multiple seeds.
  5. Run ablations: q-rank and CLS-only vs all-token adapter placement.

The model is the locally cached google/bert_uncased_L-2_H-128_A-2.  The dataset
is a tiny built-in sentiment task used only as a micro-benchmark; it is not a
GLUE/SST-2 replacement.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import random
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "flow" / "412_pandrosion_standalone_local_jet_geometry_composite_cached_irp_engine.py"
OUT = Path("/private/tmp/415_hf_q_adapter_end_to_end_benchmark.json")
HF_MODEL_ID = "google/bert_uncased_L-2_H-128_A-2"
FP32_ADAPTER_ACCEPT = 1e-6


def sync() -> None:
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    if torch.backends.mps.is_available():
        torch.mps.synchronize()


def bench(fn: Any, repeats: int = 20, warmup: int = 5) -> float:
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
    spec = importlib.util.spec_from_file_location("pandrosion412_e2e", ENGINE)
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
            sketch_basis_cache_max=4096,
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
        C = np.asarray(q_basis_np, dtype=np.float32)
        C = C[:, : int(rank)]
        if C.ndim != 2 or int(C.shape[0]) != int(hidden_size):
            raise ValueError("bad q basis shape")
        self.down = torch.nn.Linear(int(hidden_size), int(rank), bias=False)
        self.register_buffer("fixed_q_basis_t", torch.as_tensor(C.T.copy(), dtype=torch.float32))
        self.rank = int(rank)

    def forward(self, x: Any) -> Any:
        return self.down(x) @ self.fixed_q_basis_t


class AdapterClassifier(torch.nn.Module):
    def __init__(
        self,
        bert: torch.nn.Module,
        hidden_size: int,
        adapter: torch.nn.Module | None,
        scope: str,
        adapter_scale: float = 0.35,
    ) -> None:
        super().__init__()
        self.bert = bert
        self.adapter = adapter
        self.scope = str(scope)
        self.adapter_scale = float(adapter_scale)
        self.classifier = torch.nn.Linear(int(hidden_size), 2)
        for p in self.bert.parameters():
            p.requires_grad_(False)
        self.bert.eval()

    def forward(self, input_ids: Any, attention_mask: Any, token_type_ids: Any | None = None, return_adapter: bool = False) -> Any:
        with torch.no_grad():
            kwargs = {"input_ids": input_ids, "attention_mask": attention_mask}
            if token_type_ids is not None:
                kwargs["token_type_ids"] = token_type_ids
            hidden = self.bert(**kwargs).last_hidden_state.detach()
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
        logits = self.classifier(pooled)
        if return_adapter:
            return logits, adapter_delta
        return logits


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


def dataset() -> tuple[list[tuple[str, int]], list[tuple[str, int]]]:
    positive = [
        "the product feels reliable and the result is excellent",
        "this film is warm clever and genuinely enjoyable",
        "the service was fast helpful and pleasant",
        "i trust this tool because the output is accurate",
        "the update made the app smoother and easier to use",
        "the explanation is clear and the example works well",
        "a strong result with stable performance and good value",
        "the interface is clean responsive and satisfying",
        "the model answer was useful precise and readable",
        "the experiment succeeded and the metrics improved",
        "the team delivered a robust and elegant solution",
        "the benchmark is promising and the residual is tiny",
        "the adapter converged quickly with excellent accuracy",
        "the documentation is concise practical and correct",
        "the training run is stable and the loss decreases",
        "the feature is simple fast and dependable",
    ]
    negative = [
        "the product feels unreliable and the result is poor",
        "this film is dull confusing and hard to enjoy",
        "the service was slow unhelpful and unpleasant",
        "i distrust this tool because the output is inaccurate",
        "the update made the app slower and harder to use",
        "the explanation is vague and the example fails",
        "a weak result with unstable performance and bad value",
        "the interface is cluttered laggy and frustrating",
        "the model answer was useless imprecise and unreadable",
        "the experiment failed and the metrics degraded",
        "the team delivered a fragile and messy solution",
        "the benchmark is disappointing and the residual is large",
        "the adapter diverged slowly with terrible accuracy",
        "the documentation is verbose impractical and wrong",
        "the training run is unstable and the loss increases",
        "the feature is complex slow and unreliable",
    ]
    train = [(x, 1) for x in positive[:12]] + [(x, 0) for x in negative[:12]]
    val = [(x, 1) for x in positive[12:]] + [(x, 0) for x in negative[12:]]
    return train, val


def encode(tokenizer: Any, rows: list[tuple[str, int]], device: Any) -> dict[str, Any]:
    texts = [x for x, _ in rows]
    labels = torch.as_tensor([int(y) for _, y in rows], dtype=torch.long, device=device)
    enc = tokenizer(texts, return_tensors="pt", padding=True, truncation=True, max_length=32)
    out = {k: v.to(device) for k, v in enc.items()}
    out["labels"] = labels
    return out


def iter_batches(data: dict[str, Any], batch_size: int, seed: int, shuffle: bool) -> list[dict[str, Any]]:
    n = int(data["labels"].shape[0])
    idx = np.arange(n)
    if shuffle:
        rng = np.random.default_rng(int(seed))
        rng.shuffle(idx)
    batches: list[dict[str, Any]] = []
    for start in range(0, n, int(batch_size)):
        sel = torch.as_tensor(idx[start : start + int(batch_size)], dtype=torch.long, device=data["labels"].device)
        batch = {k: v.index_select(0, sel) for k, v in data.items()}
        batches.append(batch)
    return batches


def trainable_params(module: torch.nn.Module) -> int:
    return int(sum(p.numel() for p in module.parameters() if p.requires_grad))


def buffer_params(module: torch.nn.Module) -> int:
    return int(sum(b.numel() for b in module.buffers()))


def evaluate(model: AdapterClassifier, data: dict[str, Any], batch_size: int) -> dict[str, float]:
    model.eval()
    loss_fn = torch.nn.CrossEntropyLoss(reduction="sum")
    total_loss = 0.0
    correct = 0
    total = 0
    for batch in iter_batches(data, batch_size, 0, False):
        with torch.no_grad():
            logits = model(
                input_ids=batch["input_ids"],
                attention_mask=batch["attention_mask"],
                token_type_ids=batch.get("token_type_ids"),
            )
            loss = loss_fn(logits, batch["labels"])
            pred = torch.argmax(logits, dim=-1)
        total_loss += float(loss.detach().cpu().item())
        correct += int((pred == batch["labels"]).sum().detach().cpu().item())
        total += int(batch["labels"].numel())
    return {"loss": total_loss / max(1, total), "accuracy": correct / max(1, total)}


def train_model(model: AdapterClassifier, train_data: dict[str, Any], val_data: dict[str, Any], seed: int) -> dict[str, Any]:
    set_seed(seed)
    model.train()
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad], lr=4e-3, weight_decay=1e-3)
    loss_fn = torch.nn.CrossEntropyLoss()
    history: list[dict[str, float]] = []
    t0 = time.perf_counter()
    epochs = 12
    batch_size = 8
    initial = evaluate(model, val_data, batch_size)
    for ep in range(epochs):
        total = 0.0
        count = 0
        for batch in iter_batches(train_data, batch_size, seed + 1009 * ep, True):
            model.train()
            opt.zero_grad(set_to_none=True)
            logits = model(
                input_ids=batch["input_ids"],
                attention_mask=batch["attention_mask"],
                token_type_ids=batch.get("token_type_ids"),
            )
            loss = loss_fn(logits, batch["labels"])
            loss.backward()
            torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad], 5.0)
            opt.step()
            total += float(loss.detach().cpu().item()) * int(batch["labels"].numel())
            count += int(batch["labels"].numel())
        val = evaluate(model, val_data, batch_size)
        history.append({"epoch": float(ep + 1), "train_loss": total / max(1, count), "val_loss": val["loss"], "val_accuracy": val["accuracy"]})
    final = evaluate(model, val_data, batch_size)
    return {
        "epochs": epochs,
        "batch_size": batch_size,
        "initial_val_loss": initial["loss"],
        "initial_val_accuracy": initial["accuracy"],
        "final_val_loss": final["loss"],
        "final_val_accuracy": final["accuracy"],
        "train_seconds": float(time.perf_counter() - t0),
        "history": history,
    }


def make_adapter(kind: str, hidden: int, q_full_basis: Any, q_full: int, rank: int) -> tuple[torch.nn.Module | None, int]:
    if kind == "none":
        return None, 0
    if kind == "dense":
        return DenseAdapter(hidden), hidden
    if kind == "lora":
        return RandomLoRAAdapter(hidden, rank), rank
    if kind == "pandrosion_q":
        return PandrosionQAdapter(hidden, q_full_basis, rank), rank
    raise ValueError(f"unknown adapter kind {kind!r}")


def validation_roots(model: AdapterClassifier, data: dict[str, Any], max_roots: int = 3) -> list[Any]:
    batch = {k: v[: min(max_roots, int(data["labels"].shape[0]))] for k, v in data.items()}
    model.eval()
    with torch.no_grad():
        _, delta = model(
            input_ids=batch["input_ids"],
            attention_mask=batch["attention_mask"],
            token_type_ids=batch.get("token_type_ids"),
            return_adapter=True,
        )
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


def validate_q_path(mod: Any, model: AdapterClassifier, val_data: dict[str, Any], seed: int) -> dict[str, Any]:
    roots = validation_roots(model, val_data, max_roots=3)
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


def forward_latency_ms(model: AdapterClassifier, data: dict[str, Any]) -> float:
    batch = {k: v for k, v in data.items() if k != "labels"}

    def once() -> Any:
        with torch.no_grad():
            return model(**batch)

    return 1e3 * bench(once, repeats=30, warmup=5)


def stats(vals: list[float]) -> dict[str, float]:
    a = np.asarray(vals, dtype=float)
    return {
        "mean": float(np.mean(a)),
        "std": float(np.std(a, ddof=1)) if a.size > 1 else 0.0,
        "min": float(np.min(a)),
        "median": float(np.median(a)),
        "max": float(np.max(a)),
    }


def summarize(rows: list[dict[str, Any]], key: str) -> dict[str, Any]:
    groups = sorted({str(r[key]) for r in rows})
    out: dict[str, Any] = {}
    for g in groups:
        rs = [r for r in rows if str(r[key]) == g]
        out[g] = {
            "rows": len(rs),
            "final_val_accuracy": stats([float(r["final_val_accuracy"]) for r in rs]),
            "final_val_loss": stats([float(r["final_val_loss"]) for r in rs]),
            "forward_ms": stats([float(r["forward_ms"]) for r in rs]),
            "trainable_params": stats([float(r["trainable_params"]) for r in rs]),
            "q_active_rate": stats([float(r.get("q_active_rate", 0.0)) for r in rs]),
            "q_accepted_rate": stats([float(r.get("q_accepted_rate", 0.0)) for r in rs]),
            "q_energy_mean": stats([float(r.get("q_energy_mean", 0.0)) for r in rs]),
            "q_root_error_max": stats([float(r.get("q_root_error_max", 0.0)) for r in rs]),
        }
    return out


def main() -> int:
    mod = load_engine()
    configure(mod)
    device = selected_device()
    from transformers import AutoModel, AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained(HF_MODEL_ID, local_files_only=True)
    bert = AutoModel.from_pretrained(HF_MODEL_ID, local_files_only=True).to(device=device, dtype=torch.float32).eval()
    for p in bert.parameters():
        p.requires_grad_(False)
    hidden = int(bert.config.hidden_size)
    train_rows, val_rows = dataset()
    train_data = encode(tokenizer, train_rows, device)
    val_data = encode(tokenizer, val_rows, device)
    _, _, q_full = q_basis(mod, hidden, 4150000)
    equal_lora_rank = max(1, int(math.ceil(q_full / 2.0)))
    q_half_rank = max(1, equal_lora_rank)
    seeds = [41501, 41502, 41503]
    variants = [
        {"variant": "bert_head_only", "kind": "none", "scope": "cls", "rank": 0, "ablation": "baseline"},
        {"variant": "dense_adapter_cls", "kind": "dense", "scope": "cls", "rank": hidden, "ablation": "dense"},
        {"variant": "lora_equal_params_cls", "kind": "lora", "scope": "cls", "rank": equal_lora_rank, "ablation": "lora_equal_trainable"},
        {"variant": "pandrosion_q_cls_rank_half", "kind": "pandrosion_q", "scope": "cls", "rank": q_half_rank, "ablation": "rank_half"},
        {"variant": "pandrosion_q_cls_rank_full", "kind": "pandrosion_q", "scope": "cls", "rank": q_full, "ablation": "rank_full"},
        {"variant": "pandrosion_q_all_tokens_rank_full", "kind": "pandrosion_q", "scope": "all_tokens", "rank": q_full, "ablation": "all_tokens"},
    ]

    rows: list[dict[str, Any]] = []
    for seed in seeds:
        C, parent_k, q_dim = q_basis(mod, hidden, seed)
        for spec in variants:
            set_seed(seed)
            adapter, adapter_rank = make_adapter(spec["kind"], hidden, C, q_dim, int(spec["rank"]))
            model = AdapterClassifier(bert, hidden, adapter, spec["scope"]).to(device=device, dtype=torch.float32)
            train = train_model(model, train_data, val_data, seed)
            qval = validate_q_path(mod, model, val_data, seed)
            fwd_ms = forward_latency_ms(model, val_data)
            row = {
                "seed": seed,
                "variant": spec["variant"],
                "kind": spec["kind"],
                "scope": spec["scope"],
                "ablation": spec["ablation"],
                "hidden_size": hidden,
                "adapter_rank": adapter_rank,
                "parent_k": parent_k,
                "engine_q_dim": q_dim,
                "equal_lora_rank": equal_lora_rank,
                "trainable_params": trainable_params(model),
                "adapter_trainable_params": trainable_params(adapter) if adapter is not None else 0,
                "buffer_params": buffer_params(model),
                "forward_ms": fwd_ms,
                **train,
                **qval,
            }
            rows.append(row)
            print(
                f"{row['variant']} seed={seed} acc={row['final_val_accuracy']:.3f} "
                f"loss={row['final_val_loss']:.3f} forward={row['forward_ms']:.3f}ms "
                f"params={row['trainable_params']} q_rate={row.get('q_active_rate', 0.0):.2f} "
                f"q_err={row.get('q_root_error_max', 0.0):.2e}"
            )

    result = {
        "engine": str(ENGINE),
        "model_id": HF_MODEL_ID,
        "output": str(OUT),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "torch": torch.__version__,
        "transformers": __import__("transformers").__version__,
        "device": str(device),
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "dataset": {
            "type": "built_in_tiny_sentiment_micro_benchmark",
            "train_examples": len(train_rows),
            "val_examples": len(val_rows),
            "classes": 2,
            "note": "Micro-benchmark only; not a replacement for GLUE/SST-2.",
        },
        "seeds": seeds,
        "variants": variants,
        "rows": rows,
        "summary_by_variant": summarize(rows, "variant"),
        "summary_by_kind": summarize(rows, "kind"),
        "interpretation": {
            "end_to_end": "Every row runs the cached BERT model plus a classifier forward; BERT is frozen and adapter/classifier parameters are trained.",
            "lora_equal_params": "LoRA rank is ceil(q/2), giving approximately the same trainable adapter parameter count as a full-rank Pandrosion q adapter.",
            "q_validation": "q_active_rate validates adapter outputs with the standalone 412 Pandrosion IRP corrector on held-out examples.",
            "ablation": "Rank-half, rank-full, and all-token Pandrosion adapter placements are included.",
        },
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
