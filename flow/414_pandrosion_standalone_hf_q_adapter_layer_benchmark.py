#!/usr/bin/env python3
"""Real PyTorch q-adapter layer benchmark for Pandrosion IRP.

414 turns the 413 construction into actual torch.nn.Module adapters:

  - DenseAdapter: a standard dense linear correction.
  - RandomLoRAAdapter: a low-rank adapter with a random output basis.
  - PandrosionQAdapter: a low-rank adapter whose fixed output basis is the
    deterministic 412 Pandrosion q-basis.

The benchmark uses the cached Hugging Face model
google/bert_uncased_L-2_H-128_A-2.  It does not download during the run.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "flow" / "412_pandrosion_standalone_local_jet_geometry_composite_cached_irp_engine.py"
OUT = Path("/private/tmp/414_hf_pandrosion_q_adapter_layer_benchmark.json")
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


def load_engine() -> Any:
    spec = importlib.util.spec_from_file_location("pandrosion412_adapter_layer", ENGINE)
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
            sketch_basis_cache_max=2048,
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
    def __init__(self, feature_dim: int, out_dim: int) -> None:
        super().__init__()
        self.linear = torch.nn.Linear(int(feature_dim), int(out_dim), bias=False)

    def forward(self, x: Any) -> Any:
        return self.linear(x)


class RandomLoRAAdapter(torch.nn.Module):
    def __init__(self, feature_dim: int, out_dim: int, rank: int) -> None:
        super().__init__()
        self.down = torch.nn.Linear(int(feature_dim), int(rank), bias=False)
        self.up = torch.nn.Linear(int(rank), int(out_dim), bias=False)

    def forward(self, x: Any) -> Any:
        return self.up(self.down(x))


class PandrosionQAdapter(torch.nn.Module):
    def __init__(self, feature_dim: int, q_basis_np: Any) -> None:
        super().__init__()
        C = np.asarray(q_basis_np, dtype=np.float32)
        if C.ndim != 2:
            raise ValueError("q_basis must be a 2D array")
        out_dim, rank = int(C.shape[0]), int(C.shape[1])
        self.down = torch.nn.Linear(int(feature_dim), rank, bias=False)
        self.register_buffer("fixed_q_basis_t", torch.as_tensor(C.T.copy(), dtype=torch.float32))
        self.out_dim = out_dim
        self.rank = rank

    def forward(self, x: Any) -> Any:
        return self.down(x) @ self.fixed_q_basis_t.to(dtype=x.dtype, device=x.device)


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


def q_basis(mod: Any, n: int, seed: int) -> tuple[Any, int, int, float]:
    k, q = q_dimensions(mod, int(n))
    C = np.asarray(mod._coded_composite_basis_np(int(n), k, q, salt=q_salt(seed)), dtype=np.complex128)
    imag_ratio = float(np.linalg.norm(C.imag) / max(1e-300, float(np.linalg.norm(C))))
    return C.real.astype(np.float64, copy=False), k, q, imag_ratio


def subspace_energy(mod: Any, root: Any, seed: int) -> dict[str, Any]:
    r = np.asarray(root, dtype=np.complex128).reshape(-1)
    n = int(r.size)
    C, k, q, imag_ratio = q_basis(mod, n, seed)
    P = np.asarray(mod._sketch_basis_np(n, k, salt=q_salt(seed)), dtype=np.complex128)
    Cc = C.astype(np.complex128, copy=False)
    den = max(1e-300, float(np.vdot(r, r).real))
    q_coeff = Cc.conj().T @ r
    p_coeff = P.conj().T @ r
    return {
        "n": n,
        "k": k,
        "q": q,
        "q_basis_imag_ratio": imag_ratio,
        "q_energy_fraction": float(np.vdot(q_coeff, q_coeff).real / den),
        "full_sketch_energy_fraction": float(np.vdot(p_coeff, p_coeff).real / den),
        "expected_random_q_fraction": float(q / max(1, n)),
        "expected_random_full_sketch_fraction": float(k / max(1, n)),
    }


def run_corrector(mod: Any, root: Any, seed: int, repeats: int = 3) -> dict[str, Any]:
    root_np = np.asarray(root, dtype=np.complex128)

    def once() -> dict[str, Any]:
        target = VectorTarget(root_np)
        y0 = np.zeros_like(root_np)
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
        y = np.asarray(rec.get("y", y0), dtype=np.complex128)
        return {
            "status": rec.get("status"),
            "accepted": bool(rec.get("accepted")),
            "residual": float(rec.get("residual", float("inf"))),
            "root_relative_error": float(np.linalg.norm(y - root_np) / max(1e-300, float(np.linalg.norm(root_np)))),
            "probe_kind": rec.get("directional_probe_kind"),
            "nodes": int(rec.get("hypercube_nodes", 0) or 0),
            "coded_dim": rec.get("directional_coded_dim"),
            "parent_k": rec.get("directional_parent_sketch_dim"),
            "solver": rec.get("directional_reduced_solver"),
            "linear_relative_residual": rec.get("directional_linear_relative_residual"),
            "coded_fallback": rec.get("directional_coded_fallback"),
            "forward_fallback": rec.get("directional_forward_fallback"),
            "target_evals": int(target.evals),
            "q_active": bool(
                rec.get("directional_probe_kind") == "coded-probe"
                and rec.get("hypercube_nodes") == rec.get("directional_coded_dim")
            ),
        }

    rec0 = once()
    ms = 1e3 * bench(once, repeats=max(1, int(repeats)), warmup=1)
    return {"corrector_ms": ms, **rec0, **subspace_energy(mod, root_np, seed)}


def load_hf_features(model_id: str, prompts: list[str], device: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    from transformers import AutoModel, AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained(model_id, local_files_only=True)
    model = AutoModel.from_pretrained(model_id, local_files_only=True).to(device=device, dtype=torch.float32).eval()
    runs: list[dict[str, Any]] = []
    features: list[dict[str, Any]] = []
    for i, prompt in enumerate(prompts):
        x = tokenizer(prompt, return_tensors="pt", padding="max_length", truncation=True, max_length=16)
        x = {k: v.to(device) for k, v in x.items()}

        def forward_once() -> Any:
            with torch.no_grad():
                return model(**x, output_hidden_states=True)

        forward_ms = 1e3 * bench(forward_once, repeats=8, warmup=2)
        y = forward_once()
        sync()
        last = y.last_hidden_state.detach()
        mask = x.get("attention_mask", torch.ones(last.shape[:2], dtype=torch.float32, device=device)).to(last.dtype)
        mean_last = (last * mask[:, :, None]).sum(dim=1) / mask.sum(dim=1).clamp_min(1.0)[:, None]
        run = {
            "prompt_index": i,
            "prompt": prompt,
            "model_id": model_id,
            "bert_forward_ms": float(forward_ms),
            "hidden_size": int(last.shape[-1]),
            "seq_len": int(last.shape[1]),
            "device": str(device),
        }
        runs.append(run)
        samples = {
            "hf_cls_last": last[0, 0, :].detach().cpu().numpy().astype(np.float64),
            "hf_mean_last": mean_last[0].detach().cpu().numpy().astype(np.float64),
            "hf_flat_last_16x128": last[0].reshape(-1).detach().cpu().numpy().astype(np.float64),
        }
        for name, feature in samples.items():
            features.append({"activation": name, "feature": feature, **run})
    return runs, features


def init_modules(feature_dim: int, out_dim: int, rank: int, C: Any, seed: int, device: Any) -> dict[str, Any]:
    torch.manual_seed(int(seed) & 0x7FFFFFFF)
    dense = DenseAdapter(feature_dim, out_dim)
    lora = RandomLoRAAdapter(feature_dim, out_dim, rank)
    q_adapter = PandrosionQAdapter(feature_dim, C)
    with torch.no_grad():
        shared_down = torch.randn((int(rank), int(feature_dim)), dtype=torch.float32) / math.sqrt(float(max(1, feature_dim)))
        lora.down.weight.copy_(shared_down)
        q_adapter.down.weight.copy_(shared_down)
    return {
        "dense_adapter": dense.to(device=device, dtype=torch.float32).eval(),
        "random_lora_adapter": lora.to(device=device, dtype=torch.float32).eval(),
        "pandrosion_q_adapter": q_adapter.to(device=device, dtype=torch.float32).eval(),
    }


def trainable_params(module: torch.nn.Module) -> int:
    return int(sum(p.numel() for p in module.parameters() if p.requires_grad))


def buffer_params(module: torch.nn.Module) -> int:
    return int(sum(b.numel() for b in module.buffers()))


def module_output_and_time(module: torch.nn.Module, feature: Any, device: Any) -> tuple[Any, float]:
    x = torch.as_tensor(np.asarray(feature, dtype=np.float32)[None, :], dtype=torch.float32, device=device)

    def once() -> Any:
        with torch.no_grad():
            return module(x)

    ms = 1e3 * bench(once, repeats=80, warmup=10)
    y = once().detach().cpu().numpy().reshape(-1)
    return normalize_vector(y), ms


def matmul_times(n: int) -> dict[str, Any]:
    device = selected_device()
    out: dict[str, Any] = {"matmul_device": str(device)}
    dtypes = [torch.float16, torch.float32] if device.type != "cpu" else [torch.float32]
    for dtype in dtypes:
        a = torch.randn((int(n), int(n)), dtype=dtype, device=device)
        b = torch.randn((int(n), int(n)), dtype=dtype, device=device)
        repeats = 12 if n <= 128 else (8 if n <= 1024 else 5)
        try:
            t = bench(lambda: torch.matmul(a, b), repeats=repeats, warmup=3)
            name = str(dtype).replace("torch.", "")
            out[f"square_matmul_{name}_ms"] = 1e3 * t
            out[f"square_matmul_{name}_gflops"] = (2.0 * int(n) ** 3) / t / 1e9
        except Exception as exc:
            name = str(dtype).replace("torch.", "")
            out[f"square_matmul_{name}_error"] = f"{type(exc).__name__}: {exc}"
    return out


def stats(vals: list[float]) -> dict[str, float]:
    a = np.asarray(vals, dtype=float)
    return {
        "mean": float(np.mean(a)),
        "std": float(np.std(a, ddof=1)) if a.size > 1 else 0.0,
        "min": float(np.min(a)),
        "median": float(np.median(a)),
        "max": float(np.max(a)),
    }


def summarize_group(name: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "name": name,
        "rows": len(rows),
        "q_active_rate": float(np.mean([bool(r["q_active"]) for r in rows])) if rows else 0.0,
        "accepted_rate": float(np.mean([bool(r["accepted"]) for r in rows])) if rows else 0.0,
        "forward_ms": stats([float(r["adapter_forward_ms"]) for r in rows]) if rows else {},
        "corrector_ms": stats([float(r["corrector_ms"]) for r in rows]) if rows else {},
        "q_energy_fraction": stats([float(r["q_energy_fraction"]) for r in rows]) if rows else {},
        "root_relative_error": stats([float(r["root_relative_error"]) for r in rows]) if rows else {},
        "trainable_params": stats([float(r["trainable_params"]) for r in rows]) if rows else {},
    }


def main() -> int:
    mod = load_engine()
    configure(mod)
    device = selected_device()
    prompts = [
        "Pandrosion IRP benchmark on a tiny open source BERT model.",
        "Matrix multiplication kernels dominate modern AI inference workloads.",
        "Local jet geometry probes a small correction subspace.",
    ]
    hf_runs, features = load_hf_features(HF_MODEL_ID, prompts, device)
    rows: list[dict[str, Any]] = []
    matmul_cache: dict[int, dict[str, Any]] = {}

    for idx, item in enumerate(features):
        feature = np.asarray(item["feature"], dtype=np.float64).reshape(-1)
        out_dim = int(feature.size)
        feature_dim = int(feature.size)
        seed = 4140000 + 1000003 * out_dim + 7919 * int(item["prompt_index"]) + 101 * idx
        C, k, q, imag_ratio = q_basis(mod, out_dim, seed)
        modules = init_modules(feature_dim, out_dim, q, C, seed, device)
        if out_dim not in matmul_cache:
            matmul_cache[out_dim] = matmul_times(out_dim)
        dense_ms = None
        for kind, module in modules.items():
            root, forward_ms = module_output_and_time(module, feature, device)
            rec = run_corrector(mod, root, seed, repeats=4 if kind == "pandrosion_q_adapter" else 2)
            if kind == "dense_adapter":
                dense_ms = forward_ms
            row = {
                "kind": kind,
                "activation": item["activation"],
                "prompt_index": item["prompt_index"],
                "prompt": item["prompt"],
                "model_id": item["model_id"],
                "hidden_size": item["hidden_size"],
                "seq_len": item["seq_len"],
                "bert_forward_ms": item["bert_forward_ms"],
                "feature_dim": feature_dim,
                "out_dim": out_dim,
                "rank": q,
                "parent_k": k,
                "probe_seed": seed,
                "adapter_forward_ms": forward_ms,
                "trainable_params": trainable_params(module),
                "buffer_params": buffer_params(module),
                "q_basis_imag_ratio_input": imag_ratio,
                **rec,
                **matmul_cache[out_dim],
            }
            if dense_ms is not None:
                row["adapter_speedup_vs_dense_forward"] = float(dense_ms / max(1e-300, forward_ms))
            if "square_matmul_float16_ms" in row:
                row["corrector_speedup_vs_square_matmul_fp16"] = float(row["square_matmul_float16_ms"] / row["corrector_ms"])
            if "square_matmul_float32_ms" in row:
                row["corrector_speedup_vs_square_matmul_fp32"] = float(row["square_matmul_float32_ms"] / row["corrector_ms"])
            rows.append(row)
            print(
                f"{kind} {item['activation']} n={out_dim} q/k={q}/{k} "
                f"forward={forward_ms:.4f}ms q_energy={row['q_energy_fraction']:.3f} "
                f"q_active={row['q_active']} accepted={row['accepted']} "
                f"corrector={row['corrector_ms']:.3f}ms err={row['root_relative_error']:.2e}"
            )

    groups = sorted({str(r["kind"]) for r in rows})
    result = {
        "engine": str(ENGINE),
        "model_id": HF_MODEL_ID,
        "output": str(OUT),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "torch": torch.__version__,
        "transformers": __import__("transformers").__version__,
        "cpu_count": os.cpu_count(),
        "device": str(device),
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "hf_runs": hf_runs,
        "rows": rows,
        "summary": {g: summarize_group(g, [r for r in rows if r["kind"] == g]) for g in groups},
        "interpretation": {
            "dense_adapter": "Standard dense adapter; no q alignment is enforced.",
            "random_lora_adapter": "Low-rank adapter with random output basis; low rank alone should not imply 412 q-path activation.",
            "pandrosion_q_adapter": "Real torch.nn.Module with trainable down projection and fixed Pandrosion q output basis; outputs are designed to be 412 q-validable.",
        },
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
