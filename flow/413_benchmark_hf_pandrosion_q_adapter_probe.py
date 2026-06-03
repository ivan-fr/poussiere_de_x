#!/usr/bin/env python3
"""Hugging Face q-adapter probe for Pandrosion IRP 412.

This is the constructive step after the open-source AI negative control:

  - natural BERT activations are not expected to activate the 412 q-path;
  - a random low-rank/LoRA-style output is low-rank, but not aligned with the
    deterministic Pandrosion q-basis, so it should usually be rejected too;
  - a Pandrosion q-adapter emits corrections directly in the validated q-basis,
    so 412 should accept the small q-path and recover the target accurately.

The script uses the locally cached Hugging Face model
google/bert_uncased_L-2_H-128_A-2 and does not download during the benchmark.
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
OUT = Path("/private/tmp/413_hf_pandrosion_q_adapter_probe_benchmark.json")
HF_MODEL_ID = "google/bert_uncased_L-2_H-128_A-2"


def sync() -> None:
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    if torch.backends.mps.is_available():
        torch.mps.synchronize()


def bench(fn: Any, repeats: int = 5, warmup: int = 2) -> float:
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
    spec = importlib.util.spec_from_file_location("pandrosion412_q_adapter", ENGINE)
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
    C = mod._coded_composite_basis_np(int(n), k, q, salt=q_salt(seed))
    return np.asarray(C, dtype=np.complex128), k, q


def subspace_energy(mod: Any, root: Any, seed: int) -> dict[str, Any]:
    r = np.asarray(root, dtype=np.complex128).reshape(-1)
    n = int(r.size)
    C, k, q = q_basis(mod, n, seed)
    P = np.asarray(mod._sketch_basis_np(n, k, salt=q_salt(seed)), dtype=np.complex128)
    den = max(1e-300, float(np.vdot(r, r).real))
    q_coeff = C.conj().T @ r
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


def run_corrector(mod: Any, root: Any, seed: int, repeats: int = 3) -> dict[str, Any]:
    root_np = np.asarray(root, dtype=np.complex128)

    def once() -> dict[str, Any]:
        target = VectorTarget(root_np)
        y0 = np.zeros_like(root_np)
        rec = mod.hypercube_matrixjet_corrector(
            target,
            y0,
            max_epochs=2,
            tol=1e-12,
            accept=1e-9,
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
    return {"ms": ms, **rec0, **subspace_energy(mod, root_np, seed)}


def load_hf_activations(model_id: str, prompts: list[str], device: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    from transformers import AutoModel, AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained(model_id, local_files_only=True)
    model = AutoModel.from_pretrained(model_id, local_files_only=True).to(device=device, dtype=torch.float32).eval()
    runs: list[dict[str, Any]] = []
    activations: list[dict[str, Any]] = []
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
            "forward_ms": float(forward_ms),
            "hidden_size": int(last.shape[-1]),
            "seq_len": int(last.shape[1]),
            "device": str(device),
        }
        runs.append(run)
        samples = {
            "hf_cls_last": normalize_vector(last[0, 0, :].cpu().numpy()),
            "hf_mean_last": normalize_vector(mean_last[0].cpu().numpy()),
            "hf_flat_last_16x128": normalize_vector(last[0].reshape(-1).cpu().numpy()),
        }
        for name, root in samples.items():
            activations.append({"activation": name, "root": root, **run})
    return runs, activations


def random_lora_output(n: int, feature: Any, q: int, seed: int, scale: float = 0.25) -> Any:
    rng = np.random.default_rng(int(seed) ^ 0x41310A)
    basis_raw = rng.standard_normal((int(n), int(q)))
    basis, _ = np.linalg.qr(basis_raw, mode="reduced")
    feat = np.asarray(feature, dtype=np.float64).reshape(-1)
    proj = rng.standard_normal((int(q), int(feat.size))) / math.sqrt(float(max(1, feat.size)))
    coeff = proj @ feat
    out = basis @ coeff
    return normalize_vector(out, scale=scale)


def pandrosion_q_adapter_output(mod: Any, n: int, feature: Any, seed: int, scale: float = 0.25) -> Any:
    C, _, q = q_basis(mod, int(n), int(seed))
    rng = np.random.default_rng(int(seed) ^ 0x413ADA)
    feat = np.asarray(feature, dtype=np.float64).reshape(-1)
    proj = rng.standard_normal((int(q), int(feat.size))) / math.sqrt(float(max(1, feat.size)))
    coeff = proj @ feat
    nrm = float(np.linalg.norm(coeff))
    if not math.isfinite(nrm) or nrm <= 0.0:
        coeff = rng.standard_normal(int(q))
        nrm = float(np.linalg.norm(coeff))
    coeff = float(scale) * coeff / max(1e-300, nrm)
    return C @ coeff.astype(np.complex128, copy=False)


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
            out[f"matmul_{name}_ms"] = 1e3 * t
            out[f"matmul_{name}_gflops"] = (2.0 * int(n) ** 3) / t / 1e9
        except Exception as exc:
            name = str(dtype).replace("torch.", "")
            out[f"matmul_{name}_error"] = f"{type(exc).__name__}: {exc}"
    return out


def adapter_forward_times(n: int, feature_dim: int, q: int) -> dict[str, Any]:
    device = selected_device()
    out: dict[str, Any] = {"adapter_device": str(device), "adapter_feature_dim": int(feature_dim)}
    dtypes = [torch.float16, torch.float32] if device.type != "cpu" else [torch.float32]
    for dtype in dtypes:
        x = torch.randn((int(feature_dim),), dtype=dtype, device=device)
        dense = torch.randn((int(n), int(feature_dim)), dtype=dtype, device=device)
        down = torch.randn((int(q), int(feature_dim)), dtype=dtype, device=device)
        up = torch.randn((int(n), int(q)), dtype=dtype, device=device)
        repeats = 100 if n <= 128 else 30
        dense_t = bench(lambda: dense @ x, repeats=repeats, warmup=5)
        lowrank_t = bench(lambda: up @ (down @ x), repeats=repeats, warmup=5)
        name = str(dtype).replace("torch.", "")
        out[f"dense_linear_{name}_ms"] = 1e3 * dense_t
        out[f"lowrank_adapter_{name}_ms"] = 1e3 * lowrank_t
        out[f"lowrank_speedup_vs_dense_{name}"] = float(dense_t / max(1e-300, lowrank_t))
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
        "ms": stats([float(r["ms"]) for r in rows]) if rows else {},
        "q_energy_fraction": stats([float(r["q_energy_fraction"]) for r in rows]) if rows else {},
        "root_relative_error": stats([float(r["root_relative_error"]) for r in rows]) if rows else {},
        "residual": stats([float(r["residual"]) for r in rows]) if rows else {},
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
    hf_runs, activations = load_hf_activations(HF_MODEL_ID, prompts, device)
    rows: list[dict[str, Any]] = []
    matmul_cache: dict[int, dict[str, Any]] = {}
    adapter_cache: dict[tuple[int, int, int], dict[str, Any]] = {}

    for idx, item in enumerate(activations):
        natural = np.asarray(item["root"], dtype=np.complex128)
        n = int(natural.size)
        probe_seed = 4130000 + 1000003 * n + 7919 * int(item["prompt_index"]) + 101 * idx
        _, _, q = q_basis(mod, n, probe_seed)
        feature = natural.real
        candidates = [
            ("natural_hf_activation", natural),
            ("random_lora_output", random_lora_output(n, feature, q, probe_seed)),
            ("pandrosion_q_adapter_output", pandrosion_q_adapter_output(mod, n, feature, probe_seed)),
        ]
        if n not in matmul_cache:
            matmul_cache[n] = matmul_times(n)
        cache_key = (n, int(feature.size), q)
        if cache_key not in adapter_cache:
            adapter_cache[cache_key] = adapter_forward_times(n, int(feature.size), q)
        for kind, root in candidates:
            rec = run_corrector(mod, root, probe_seed, repeats=4 if kind == "pandrosion_q_adapter_output" else 2)
            row = {
                "kind": kind,
                "activation": item["activation"],
                "prompt_index": item["prompt_index"],
                "prompt": item["prompt"],
                "model_id": item["model_id"],
                "hidden_size": item["hidden_size"],
                "seq_len": item["seq_len"],
                "hf_forward_ms": item["forward_ms"],
                "probe_seed": probe_seed,
                **rec,
                **matmul_cache[n],
                **adapter_cache[cache_key],
            }
            if "matmul_float16_ms" in row:
                row["speedup_vs_square_matmul_fp16"] = float(row["matmul_float16_ms"] / row["ms"])
            if "matmul_float32_ms" in row:
                row["speedup_vs_square_matmul_fp32"] = float(row["matmul_float32_ms"] / row["ms"])
            rows.append(row)
            print(
                f"{kind} {item['activation']} n={n} q/k={row['q']}/{row['k']} "
                f"q_energy={row['q_energy_fraction']:.3f} q_active={row['q_active']} "
                f"accepted={row['accepted']} 412={row['ms']:.3f}ms "
                f"err={row['root_relative_error']:.2e}"
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
            "natural_hf_activation": "Real cached Hugging Face BERT activations; q should only activate if the model naturally aligns with the deterministic q-basis.",
            "random_lora_output": "A classical low-rank output using a random output basis; low rank alone is not enough for the 412 q-path.",
            "pandrosion_q_adapter_output": "An adapter whose output basis is the deterministic Pandrosion q-basis for the layer/seed; this is the constructive q-compatible AI path.",
        },
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
