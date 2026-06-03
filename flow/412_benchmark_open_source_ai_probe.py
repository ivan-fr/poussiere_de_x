#!/usr/bin/env python3
"""Probe 412 Pandrosion IRP on local open-source-style PyTorch AI workloads.

This benchmark intentionally separates two questions:

1. Do natural activations from a PyTorch Transformer block activate the 412 q-path?
2. When an AI-like vector is explicitly aligned with the 412 q-basis, does the
   same engine still validate and accelerate the correction?

No model download is required.  If Hugging Face transformers is installed in the
environment, the script records that fact, but it only depends on stdlib, NumPy,
PyTorch, and the standalone 412 engine.
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
OUT = Path("/private/tmp/412_open_source_ai_probe_benchmark.json")
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
    return min(vals)


def load_engine() -> Any:
    spec = importlib.util.spec_from_file_location("pandrosion412_ai_probe", ENGINE)
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
            sketch_basis_cache_max=1024,
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


def try_import_version(name: str) -> dict[str, Any]:
    try:
        mod = __import__(name)
        return {"available": True, "version": str(getattr(mod, "__version__", "unknown"))}
    except Exception as exc:
        return {"available": False, "error": f"{type(exc).__name__}: {exc}"}


class TinyTransformerBlock(torch.nn.Module):
    """Small decoder-style block using standard torch.nn building blocks."""

    def __init__(self, d_model: int, nhead: int, mlp_mult: int = 4) -> None:
        super().__init__()
        hidden = int(d_model) * int(mlp_mult)
        self.ln1 = torch.nn.LayerNorm(int(d_model))
        self.attn = torch.nn.MultiheadAttention(
            int(d_model),
            int(nhead),
            dropout=0.0,
            bias=True,
            batch_first=True,
        )
        self.ln2 = torch.nn.LayerNorm(int(d_model))
        self.fc1 = torch.nn.Linear(int(d_model), hidden)
        self.fc2 = torch.nn.Linear(hidden, int(d_model))

    def forward(self, x: Any) -> tuple[Any, dict[str, Any]]:
        z = self.ln1(x)
        attn, _ = self.attn(z, z, z, need_weights=False)
        x = x + attn
        hidden = torch.nn.functional.gelu(self.fc1(self.ln2(x)))
        out = x + self.fc2(hidden)
        return out, {
            "attn_token": attn[:, 0, :],
            "block_token": out[:, 0, :],
            "mlp_hidden": hidden[:, 0, :],
        }


class VectorTarget:
    """Linear target F(z)=z-root for q-path activation probing."""

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
        raise ValueError("empty_or_zero_activation")
    return (float(scale) * v / nrm).astype(np.complex128)


def q_dimensions(mod: Any, n: int) -> tuple[int, int]:
    k = max(1, min(int(n), int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))))
    q = int(mod.directional_coded_dim_for(int(n), int(k)))
    return int(k), int(q)


def q_salt(seed: int) -> int:
    return int(seed) + 0x412D1A


def subspace_energy(mod: Any, root: Any, seed: int) -> dict[str, Any]:
    r = np.asarray(root, dtype=np.complex128).reshape(-1)
    n = int(r.size)
    k, q = q_dimensions(mod, n)
    salt = q_salt(seed)
    c_basis = mod._coded_composite_basis_np(n, k, q, salt=salt)
    p_basis = mod._sketch_basis_np(n, k, salt=salt)
    den = max(1e-300, float(np.vdot(r, r).real))
    q_coeff = c_basis.conj().T @ r
    p_coeff = p_basis.conj().T @ r
    return {
        "n": n,
        "k": k,
        "q": q,
        "q_energy_fraction": float(np.vdot(q_coeff, q_coeff).real / den),
        "full_sketch_energy_fraction": float(np.vdot(p_coeff, p_coeff).real / den),
        "expected_random_q_fraction": float(q / max(1, n)),
        "expected_random_full_sketch_fraction": float(k / max(1, n)),
    }


def make_q_aligned_root(mod: Any, n: int, seed: int, scale: float = 0.25) -> Any:
    k, q = q_dimensions(mod, int(n))
    c_basis = mod._coded_composite_basis_np(int(n), k, q, salt=q_salt(seed))
    rng = np.random.default_rng(int(seed) ^ 0x412A17)
    coeff = rng.standard_normal(q) + 1j * rng.standard_normal(q)
    coeff = float(scale) * coeff / max(1e-300, float(np.linalg.norm(coeff)))
    return c_basis @ coeff


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
    return {"ms": float(ms), **rec0, **subspace_energy(mod, root_np, seed)}


def selected_device() -> Any:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def make_transformer_samples(
    d_model: int,
    seq_len: int,
    batch: int,
    nhead: int,
    seed: int,
    device: Any,
) -> tuple[float, dict[str, Any]]:
    torch.manual_seed(int(seed))
    block = TinyTransformerBlock(int(d_model), int(nhead)).to(device=device, dtype=torch.float32).eval()
    x = torch.randn((int(batch), int(seq_len), int(d_model)), dtype=torch.float32, device=device)

    def forward_once() -> Any:
        with torch.no_grad():
            return block(x)

    forward_ms = 1e3 * bench(forward_once, repeats=8 if d_model <= 512 else 5, warmup=2)
    out, acts = forward_once()
    sync()
    del out
    samples = {
        "attn_token": normalize_vector(acts["attn_token"][0].detach().cpu().numpy()),
        "block_token": normalize_vector(acts["block_token"][0].detach().cpu().numpy()),
    }
    if d_model <= 512:
        samples["mlp_hidden"] = normalize_vector(acts["mlp_hidden"][0].detach().cpu().numpy())
    return float(forward_ms), samples


def make_huggingface_samples(model_id: str, prompt: str, device: Any) -> tuple[float, dict[str, Any], dict[str, Any]]:
    from transformers import AutoModel, AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained(model_id, local_files_only=True)
    model = AutoModel.from_pretrained(model_id, local_files_only=True).to(device=device, dtype=torch.float32).eval()
    x = tokenizer(
        prompt,
        return_tensors="pt",
        padding="max_length",
        truncation=True,
        max_length=16,
    )
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
    samples = {
        "hf_cls_last": normalize_vector(last[0, 0, :].cpu().numpy()),
        "hf_mean_last": normalize_vector(mean_last[0].cpu().numpy()),
        "hf_flat_last_16x128": normalize_vector(last[0].reshape(-1).cpu().numpy()),
    }
    meta = {
        "model_id": model_id,
        "hidden_size": int(last.shape[-1]),
        "seq_len": int(last.shape[1]),
        "forward_device": str(device),
    }
    return float(forward_ms), samples, meta


def matmul_times(n: int) -> dict[str, Any]:
    device = selected_device()
    out: dict[str, Any] = {"device": str(device)}
    dtypes = [torch.float16, torch.float32] if device.type != "cpu" else [torch.float32]
    for dtype in dtypes:
        a = torch.randn((int(n), int(n)), dtype=dtype, device=device)
        b = torch.randn((int(n), int(n)), dtype=dtype, device=device)
        repeats = 8 if n <= 1024 else (5 if n <= 2048 else 3)
        try:
            t = bench(lambda: torch.matmul(a, b), repeats=repeats, warmup=2)
            name = str(dtype).replace("torch.", "")
            out[f"matmul_{name}_ms"] = 1e3 * t
            out[f"matmul_{name}_gflops"] = (2.0 * int(n) ** 3) / t / 1e9
        except Exception as exc:
            name = str(dtype).replace("torch.", "")
            out[f"matmul_{name}_error"] = f"{type(exc).__name__}: {exc}"
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
        "root_relative_error": stats([float(r["root_relative_error"]) for r in rows]) if rows else {},
        "q_energy_fraction": stats([float(r["q_energy_fraction"]) for r in rows]) if rows else {},
        "full_sketch_energy_fraction": stats([float(r["full_sketch_energy_fraction"]) for r in rows]) if rows else {},
    }


def main() -> int:
    mod = load_engine()
    configure(mod)
    device = selected_device()
    seeds = [412101, 412102, 412103]
    configs = [
        {"d_model": 512, "seq_len": 64, "batch": 2, "nhead": 8},
        {"d_model": 1024, "seq_len": 32, "batch": 1, "nhead": 8},
    ]

    natural_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    transformer_runs: list[dict[str, Any]] = []
    huggingface_runs: list[dict[str, Any]] = []
    matmul_cache: dict[int, dict[str, Any]] = {}

    for cfg in configs:
        for seed in seeds:
            try:
                forward_ms, samples = make_transformer_samples(
                    cfg["d_model"],
                    cfg["seq_len"],
                    cfg["batch"],
                    cfg["nhead"],
                    seed,
                    device,
                )
                transformer_runs.append({"seed": seed, "device": str(device), "forward_ms": forward_ms, **cfg})
            except Exception as exc:
                transformer_runs.append({"seed": seed, "device": str(device), "error": f"{type(exc).__name__}: {exc}", **cfg})
                continue
            for name, root in samples.items():
                n = int(np.asarray(root).size)
                if n not in matmul_cache:
                    matmul_cache[n] = matmul_times(n)
                probe_seed = int(seed) + 1000003 * n + 7919 * (len(natural_rows) + 1)
                row = {
                    "kind": "natural_transformer_activation",
                    "activation": name,
                    "seed": seed,
                    "probe_seed": probe_seed,
                    "transformer_forward_ms": forward_ms,
                    **cfg,
                    **run_corrector(mod, root, probe_seed, repeats=2),
                    **matmul_cache[n],
                }
                natural_rows.append(row)
                print(
                    f"natural {name} n={n} seed={seed} q/k={row['q']}/{row['k']} "
                    f"q_energy={row['q_energy_fraction']:.3f} q_active={row['q_active']} "
                    f"accepted={row['accepted']} err={row['root_relative_error']:.3f}"
                )

    hf_prompts = [
        "Pandrosion IRP benchmark on a tiny open source BERT model.",
        "Matrix multiplication kernels dominate modern AI inference workloads.",
        "Local jet geometry probes a small correction subspace.",
    ]
    for i, prompt in enumerate(hf_prompts):
        try:
            forward_ms, samples, hf_meta = make_huggingface_samples(HF_MODEL_ID, prompt, device)
            huggingface_runs.append({"prompt_index": i, "prompt": prompt, "forward_ms": forward_ms, **hf_meta})
        except Exception as exc:
            if device.type != "cpu":
                try:
                    cpu_device = torch.device("cpu")
                    forward_ms, samples, hf_meta = make_huggingface_samples(HF_MODEL_ID, prompt, cpu_device)
                    huggingface_runs.append({"prompt_index": i, "prompt": prompt, "forward_ms": forward_ms, "device_fallback": "cpu", **hf_meta})
                except Exception as exc2:
                    huggingface_runs.append({"prompt_index": i, "prompt": prompt, "error": f"{type(exc).__name__}: {exc}; cpu: {type(exc2).__name__}: {exc2}"})
                    continue
            else:
                huggingface_runs.append({"prompt_index": i, "prompt": prompt, "error": f"{type(exc).__name__}: {exc}"})
                continue
        for name, root in samples.items():
            n = int(np.asarray(root).size)
            if n not in matmul_cache:
                matmul_cache[n] = matmul_times(n)
            probe_seed = 91357 + 1000003 * n + 7919 * (i + 1) + 101 * len(name)
            row = {
                "kind": "huggingface_tiny_bert_activation",
                "activation": name,
                "prompt_index": i,
                "prompt": prompt,
                "probe_seed": probe_seed,
                "transformer_forward_ms": forward_ms,
                **hf_meta,
                **run_corrector(mod, root, probe_seed, repeats=2),
                **matmul_cache[n],
            }
            natural_rows.append(row)
            print(
                f"hf {name} n={n} prompt={i} q/k={row['q']}/{row['k']} "
                f"q_energy={row['q_energy_fraction']:.3f} q_active={row['q_active']} "
                f"accepted={row['accepted']} err={row['root_relative_error']:.3f}"
            )

    for n in sorted(matmul_cache):
        for seed in seeds:
            probe_seed = int(seed) + 2000003 * n
            root = make_q_aligned_root(mod, n, probe_seed)
            row = {
                "kind": "q_aligned_ai_vector_control",
                "seed": seed,
                "probe_seed": probe_seed,
                **run_corrector(mod, root, probe_seed, repeats=4 if n <= 1024 else 3),
                **matmul_cache[n],
            }
            if "matmul_float16_ms" in row:
                row["speedup_vs_matmul_fp16"] = float(row["matmul_float16_ms"] / row["ms"])
            if "matmul_float32_ms" in row:
                row["speedup_vs_matmul_fp32"] = float(row["matmul_float32_ms"] / row["ms"])
            control_rows.append(row)
            print(
                f"q-control n={n} seed={seed} q/k={row['q']}/{row['k']} "
                f"412={row['ms']:.3f}ms q_active={row['q_active']} "
                f"fp16_speed={row.get('speedup_vs_matmul_fp16', 0.0):.2f} "
                f"err={row['root_relative_error']:.2e}"
            )

    result = {
        "engine": str(ENGINE),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cpu_count": os.cpu_count(),
        "device": str(device),
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "transformers": try_import_version("transformers"),
        "torchvision": try_import_version("torchvision"),
        "configs": configs,
        "seeds": seeds,
        "transformer_runs": transformer_runs,
        "huggingface_runs": huggingface_runs,
        "natural_rows": natural_rows,
        "control_rows": control_rows,
        "summary": {
            "natural_all": summarize_group("natural_all", natural_rows),
            "local_tiny_transformer_activation": summarize_group("local_tiny_transformer_activation", [r for r in natural_rows if r.get("kind") == "natural_transformer_activation"]),
            "huggingface_tiny_bert_activation": summarize_group("huggingface_tiny_bert_activation", [r for r in natural_rows if r.get("kind") == "huggingface_tiny_bert_activation"]),
            "q_aligned_ai_vector_control": summarize_group("q_aligned_ai_vector_control", control_rows),
        },
        "interpretation": {
            "natural": "A natural activation is a real AI workload probe. q_active=false means 412 correctly rejects the small q-path for that activation.",
            "control": "The q-aligned control checks the same 412 path under the condition where Pandrosion IRP has a validated small q basis.",
        },
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
