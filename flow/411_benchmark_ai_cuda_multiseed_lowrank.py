#!/usr/bin/env python3
"""Extra 411 benchmarks: CUDA status, AI workloads, multi-seed stats, sketch baselines."""
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
import torch.nn.functional as F


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "flow" / "411_pandrosion_standalone_local_jet_geometry_cached_general_irp_engine.py"
OUT = Path("/private/tmp/411_ai_cuda_multiseed_lowrank_benchmark.json")


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
    spec = importlib.util.spec_from_file_location("pandrosion411_extra_bench", ENGINE)
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
            sketch_seed=411,
            sketch_mode="cached-rademacher",
            sketch_solver="svd",
            sketch_basis_cache=True,
            sketch_basis_cache_max=512,
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
            directional_basis_reuse=True,
        )
    )


class SparseQuadraticTarget:
    def __init__(self, root: Any, beta: float = 1e-8, gamma: float = 0.0) -> None:
        self.root = np.asarray(root, dtype=np.complex128)
        self.beta = float(beta)
        self.gamma = float(gamma)
        self.c = self._poly(self.root)
        self.evals = 0

    @property
    def n(self) -> int:
        return int(self.root.size)

    def _poly(self, Z: Any) -> Any:
        ZZ = np.asarray(Z, dtype=np.complex128)
        return ZZ + self.beta * ZZ * ZZ + self.gamma * ZZ * np.roll(ZZ, -1, axis=-1)

    def eval(self, z: Any) -> Any:
        self.evals += 1
        return self._poly(z) - self.c

    def eval_batch(self, Z: Any) -> Any:
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        self.evals += int(ZZ.shape[0])
        return self._poly(ZZ) - self.c[None, :]

    def residual(self, z: Any) -> float:
        return float(np.linalg.norm(self.eval(z)))

    def residuals_batch(self, Z: Any) -> Any:
        return np.linalg.norm(self.eval_batch(Z), axis=1)


def make_coded_root(mod: Any, n: int, seed: int, scale: float = 0.25) -> tuple[Any, int, int]:
    k = max(1, min(n, int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))))
    q = mod.directional_coded_dim_for(n, k)
    salt = int(seed) + 0x411D1A
    P = mod._sketch_basis_np(n, k, salt=salt)
    R = mod._fast_projected_coded_basis_np(k, q, salt=salt + 0xC0DEC0DE)
    rng = np.random.default_rng(seed ^ 0x411A1)
    coeff = scale * (rng.standard_normal(q) + 1j * rng.standard_normal(q)) / math.sqrt(2.0 * q)
    return (P @ R) @ coeff, k, q


def run_411_corrector(mod: Any, root: Any, seed: int, beta: float = 1e-8, gamma: float = 0.0) -> dict[str, Any]:
    target = SparseQuadraticTarget(root, beta, gamma)
    y0 = np.zeros(target.n, dtype=np.complex128)
    rec = mod.hypercube_matrixjet_corrector(
        target,
        y0,
        max_epochs=4,
        tol=1e-12,
        accept=1e-9,
        trial_timeout=0.0,
        line_search=8,
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
        "root_relative_error": float(np.linalg.norm(y - root) / max(1e-300, float(np.linalg.norm(root)))),
        "probe_kind": rec.get("directional_probe_kind"),
        "nodes": int(rec.get("hypercube_nodes", 0) or 0),
        "coded_dim": rec.get("directional_coded_dim"),
        "parent_k": rec.get("directional_parent_sketch_dim"),
        "solver": rec.get("directional_reduced_solver"),
        "linear_rel": rec.get("directional_linear_relative_residual"),
        "target_evals": int(target.evals),
        "q_active": bool(
            rec.get("directional_probe_kind") == "coded-probe"
            and rec.get("hypercube_nodes") == rec.get("directional_coded_dim")
        ),
    }


def matmul_baseline(n: int) -> dict[str, Any]:
    device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))
    out: dict[str, Any] = {"device": str(device)}
    for dtype in ([torch.float16, torch.float32] if device.type != "cpu" else [torch.float32]):
        a = torch.randn((n, n), dtype=dtype, device=device)
        b = torch.randn((n, n), dtype=dtype, device=device)
        repeats = 8 if n <= 1024 else (5 if n <= 2048 else 3)
        ms = bench(lambda: torch.matmul(a, b), repeats=repeats, warmup=3) * 1e3
        out[f"matmul_{str(dtype).replace('torch.', '')}_ms"] = ms
    return out


def multiseed_q_stats(mod: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    dims = [512, 1024, 2048, 4096]
    seeds = [41101, 41102, 41103, 41104, 41105]
    for n in dims:
        mm = matmul_baseline(n)
        samples: list[dict[str, Any]] = []
        for seed0 in seeds:
            seed = seed0 + n
            root, k, q = make_coded_root(mod, n, seed)

            def run() -> dict[str, Any]:
                return run_411_corrector(mod, root, seed, beta=1e-8, gamma=0.0)

            rec = run()
            repeats = 4 if n <= 1024 else (3 if n <= 2048 else 2)
            ms = bench(run, repeats=repeats, warmup=1) * 1e3
            samples.append({"seed": seed, "ms": ms, "k": k, "q": q, **rec})
        times = np.asarray([s["ms"] for s in samples], dtype=float)
        residuals = np.asarray([s["residual"] for s in samples], dtype=float)
        root_errors = np.asarray([s["root_relative_error"] for s in samples], dtype=float)
        row: dict[str, Any] = {
            "n": n,
            "q": int(samples[0]["q"]),
            "k": int(samples[0]["k"]),
            "seeds": len(samples),
            "411_ms_mean": float(np.mean(times)),
            "411_ms_std": float(np.std(times, ddof=1)) if len(times) > 1 else 0.0,
            "411_ms_min": float(np.min(times)),
            "411_ms_max": float(np.max(times)),
            "q_active_rate": float(np.mean([bool(s["q_active"]) for s in samples])),
            "accepted_rate": float(np.mean([bool(s["accepted"]) for s in samples])),
            "residual_median": float(np.median(residuals)),
            "residual_max": float(np.max(residuals)),
            "root_error_median": float(np.median(root_errors)),
            "root_error_max": float(np.max(root_errors)),
            "samples": samples,
            **mm,
        }
        if "matmul_float16_ms" in mm:
            row["mean_speedup_vs_matmul_fp16"] = float(mm["matmul_float16_ms"] / row["411_ms_mean"])
        if "matmul_float32_ms" in mm:
            row["mean_speedup_vs_matmul_fp32"] = float(mm["matmul_float32_ms"] / row["411_ms_mean"])
        rows.append(row)
    return rows


def ai_workload_baselines() -> list[dict[str, Any]]:
    device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))
    rows: list[dict[str, Any]] = []
    dtypes = [torch.float16, torch.float32] if device.type != "cpu" else [torch.float32]
    for n in [512, 1024, 2048]:
        for dtype in dtypes:
            b = 32
            x = torch.randn((b, n), dtype=dtype, device=device)
            w = torch.randn((n, n), dtype=dtype, device=device)
            ms = bench(lambda: F.linear(x, w), repeats=10 if n <= 1024 else 6, warmup=3) * 1e3
            rows.append({"workload": "linear_batch32", "n": n, "dtype": str(dtype).replace("torch.", ""), "device": str(device), "ms": ms})
            h = min(4 * n, 4096)
            w1 = torch.randn((h, n), dtype=dtype, device=device)
            w2 = torch.randn((n, h), dtype=dtype, device=device)
            ms = bench(lambda: F.linear(F.gelu(F.linear(x, w1)), w2), repeats=8 if n <= 1024 else 4, warmup=2) * 1e3
            rows.append({"workload": "mlp_batch32_hidden_capped", "n": n, "hidden": h, "dtype": str(dtype).replace("torch.", ""), "device": str(device), "ms": ms})
    for n in [512, 1024]:
        for dtype in dtypes:
            seq, batch = 128, 1
            x = torch.randn((batch, seq, n), dtype=dtype, device=device)
            wq = torch.randn((n, n), dtype=dtype, device=device)
            wk = torch.randn((n, n), dtype=dtype, device=device)
            wv = torch.randn((n, n), dtype=dtype, device=device)

            def attention_projection() -> Any:
                q = F.linear(x, wq)
                k = F.linear(x, wk)
                v = F.linear(x, wv)
                scores = torch.matmul(q, k.transpose(-2, -1)) / math.sqrt(float(n))
                return torch.matmul(torch.softmax(scores, dim=-1), v)

            ms = bench(attention_projection, repeats=6, warmup=2) * 1e3
            rows.append({"workload": "single_head_attention_seq128", "n": n, "seq": seq, "dtype": str(dtype).replace("torch.", ""), "device": str(device), "ms": ms})
    return rows


def classical_sketch_baselines(mod: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for n in [512, 1024, 2048]:
        seed = 512000 + n
        root, k, q = make_coded_root(mod, n, seed)
        rec = run_411_corrector(mod, root, seed, beta=1e-8, gamma=0.0)
        t411 = bench(lambda: run_411_corrector(mod, root, seed, beta=1e-8, gamma=0.0), repeats=4 if n <= 1024 else 3, warmup=1) * 1e3
        rng = np.random.default_rng(seed ^ 0xC1A551C)

        def random_projection(dim: int) -> tuple[float, float]:
            t0 = time.perf_counter()
            raw = rng.choice([-1.0, 1.0], size=(n, dim)).astype(np.float64) / math.sqrt(float(n))
            Q, _ = np.linalg.qr(raw, mode="reduced")
            delta = Q @ (Q.conj().T @ root)
            return 1e3 * (time.perf_counter() - t0), float(np.linalg.norm(delta - root) / max(1e-300, float(np.linalg.norm(root))))

        def coordinate_top(dim: int) -> tuple[float, float]:
            t0 = time.perf_counter()
            idx = np.argpartition(np.abs(root), -dim)[-dim:]
            delta = np.zeros_like(root)
            delta[idx] = root[idx]
            return 1e3 * (time.perf_counter() - t0), float(np.linalg.norm(delta - root) / max(1e-300, float(np.linalg.norm(root))))

        rq_ms, rq_err = random_projection(q)
        rk_ms, rk_err = random_projection(k)
        cq_ms, cq_err = coordinate_top(q)
        rows.append(
            {
                "n": n,
                "q": q,
                "k": k,
                "411_ms": t411,
                "411_residual": rec["residual"],
                "411_root_error": rec["root_relative_error"],
                "411_q_active": rec["q_active"],
                "random_projection_q_ms": rq_ms,
                "random_projection_q_error": rq_err,
                "random_projection_k_ms": rk_ms,
                "random_projection_k_error": rk_err,
                "coordinate_top_q_ms": cq_ms,
                "coordinate_top_q_error": cq_err,
            }
        )
    return rows


def cuda_tensor_core_status() -> dict[str, Any]:
    out: dict[str, Any] = {
        "cuda_available": bool(torch.cuda.is_available()),
        "cuda_device_count": int(torch.cuda.device_count()) if torch.cuda.is_available() else 0,
        "note": "",
        "rows": [],
    }
    if not torch.cuda.is_available():
        out["note"] = "CUDA/Tensor Cores are not available on this machine; benchmark script records this as unavailable."
        return out
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.set_float32_matmul_precision("high")
    device = torch.device("cuda")
    out["device_name"] = torch.cuda.get_device_name(0)
    rows: list[dict[str, Any]] = []
    for n in [1024, 2048, 4096]:
        for dtype in [torch.float16, torch.float32]:
            a = torch.randn((n, n), dtype=dtype, device=device)
            b = torch.randn((n, n), dtype=dtype, device=device)
            ms = bench(lambda: torch.matmul(a, b), repeats=8 if n <= 2048 else 5, warmup=4) * 1e3
            rows.append({"n": n, "dtype": str(dtype).replace("torch.", ""), "cuda_matmul_ms": ms})
    out["rows"] = rows
    return out


def main() -> int:
    mod = load_engine()
    configure(mod)
    result = {
        "engine": str(ENGINE),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cpu_count": os.cpu_count(),
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "cuda_tensor_cores": cuda_tensor_core_status(),
        "multiseed_q_active": multiseed_q_stats(mod),
        "ai_workload_baselines": ai_workload_baselines(),
        "classical_sketch_baselines": classical_sketch_baselines(mod),
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print("CUDA:", result["cuda_tensor_cores"].get("note") or result["cuda_tensor_cores"].get("device_name"))
    print("multi-seed q-active:")
    for r in result["multiseed_q_active"]:
        print(
            f"  n={r['n']} q/k={r['q']}/{r['k']} seeds={r['seeds']} "
            f"411={r['411_ms_mean']:.3f}+/-{r['411_ms_std']:.3f}ms "
            f"q_rate={r['q_active_rate']:.2f} fp16_speed={r.get('mean_speedup_vs_matmul_fp16', 0):.2f} "
            f"res_max={r['residual_max']:.2e}"
        )
    print("AI workload baselines:")
    for r in result["ai_workload_baselines"]:
        print(f"  {r['workload']} n={r['n']} dtype={r['dtype']} ms={r['ms']:.3f}")
    print("classic sketch baselines:")
    for r in result["classical_sketch_baselines"]:
        print(
            f"  n={r['n']} q/k={r['q']}/{r['k']} 411_err={r['411_root_error']:.2e} "
            f"rand_q_err={r['random_projection_q_error']:.2e} "
            f"rand_k_err={r['random_projection_k_error']:.2e} "
            f"coord_q_err={r['coordinate_top_q_error']:.2e}"
        )
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
