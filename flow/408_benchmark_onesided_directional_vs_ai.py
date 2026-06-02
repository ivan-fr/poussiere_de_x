#!/usr/bin/env python3
"""Benchmark 408 one-sided directional Pandrosion IRP against full jets and AI matmul."""
from __future__ import annotations

import importlib.util
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import torch


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "flow" / "408_pandrosion_standalone_local_jet_geometry_onesided_directional_irp_engine.py"
OUT = Path("/private/tmp/408_pandrosion_onesided_directional_benchmark.json")


def sync() -> None:
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    if torch.backends.mps.is_available():
        torch.mps.synchronize()


def bench(fn: Any, repeats: int = 5, warmup: int = 2) -> float:
    for _ in range(max(0, warmup)):
        fn()
    sync()
    t0 = time.perf_counter()
    for _ in range(max(1, repeats)):
        fn()
    sync()
    return (time.perf_counter() - t0) / max(1, repeats)


def load_engine() -> Any:
    spec = importlib.util.spec_from_file_location("pandrosion408", ENGINE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {ENGINE}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pandrosion408"] = mod
    spec.loader.exec_module(mod)
    return mod


def configure(
    mod: Any,
    algorithm: str,
    *,
    directional: bool = True,
    directional_diff: str = "auto",
    sketch_solver: str = "svd",
    sketch_factor: float = 2.75,
    directional_factor: float = 2.75,
    sketch_min_n: int = 96,
    ns_iters: int = 8,
) -> None:
    args = SimpleNamespace(
        matrix_backend="torch",
        matrix_algorithm=algorithm,
        batch_kernel="auto",
        torch_device="auto",
        torch_complex_dtype="auto",
        torch_real_dtype="auto",
        ns_iters=ns_iters,
        ns_damping_floor=1e-5,
        sketch_dim=0,
        sketch_factor=sketch_factor,
        sketch_min_n=sketch_min_n,
        sketch_seed=408,
        sketch_mode="rademacher",
        sketch_solver=sketch_solver,
        directional_jet=directional,
        directional_jet_min_n=sketch_min_n,
        directional_jet_factor=directional_factor,
        directional_diff=directional_diff,
        directional_auto_central_cap=0.35,
    )
    mod.configure_matrix_backend(args)


class LinearTarget:
    def __init__(self, A: Any, b: Any) -> None:
        self.A = np.asarray(A, dtype=np.complex128)
        self.b = np.asarray(b, dtype=np.complex128)
        self.evals = 0

    def eval(self, y: Any) -> Any:
        yy = np.asarray(y, dtype=np.complex128)
        self.evals += 1
        return yy @ self.A.T - self.b

    def eval_batch(self, Y: Any) -> Any:
        YY = np.asarray(Y, dtype=np.complex128)
        if YY.ndim == 1:
            return self.eval(YY)[None, :]
        self.evals += int(YY.shape[0])
        return YY @ self.A.T - self.b[None, :]

    def residuals_batch(self, Y: Any) -> Any:
        return np.linalg.norm(self.eval_batch(Y), axis=1)


def directional_variant_bench(mod: Any) -> list[dict[str, Any]]:
    rng = np.random.default_rng(408)
    rows: list[dict[str, Any]] = []
    for n in [128, 256, 512, 1024]:
        configure(mod, "adaptive-directional-sketch", directional=True, directional_diff="auto", sketch_solver="svd", sketch_min_n=1)
        k = int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))
        seed = 9100 + n
        ep = 0
        salt = seed + 65537 * ep + 0x408D1A
        P = mod._sketch_basis_np(n, k, salt=salt)
        A = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / math.sqrt(n)
        u = rng.normal(size=k) + 1j * rng.normal(size=k)
        x_true = P @ u
        b = A @ x_true

        variants: list[dict[str, Any]] = []
        for diff in ["forward", "auto", "central"]:
            configure(mod, "adaptive-directional-sketch", directional=True, directional_diff=diff, sketch_solver="svd", sketch_min_n=1)

            def run_directional() -> tuple[Any, dict[str, Any], int]:
                target = LinearTarget(A, b)
                y0 = np.zeros(n, dtype=np.complex128)
                f0 = target.eval(y0)
                delta, meta = mod.hypercube_delta(target, y0, f0, ep, seed, 0, 0.0, 0.0, 1e-12, 1e12)
                return delta, meta, target.evals

            delta, meta, evals = run_directional()
            repeats = 5 if n <= 512 else 4
            t_dir = bench(run_directional, repeats=repeats, warmup=2)
            variants.append(
                {
                    "requested_diff": diff,
                    "used_diff": str(meta.get("directional_diff")),
                    "fallback": meta.get("directional_forward_fallback"),
                    "ms": 1e3 * t_dir,
                    "nodes": int(meta.get("hypercube_nodes", 0)),
                    "evals": int(evals),
                    "nodes_avoided": int(meta.get("directional_full_hypercube_samples_avoided", 0)),
                    "rel_step_error": float(np.linalg.norm(delta - x_true) / max(1e-300, np.linalg.norm(x_true))),
                    "linear_residual": float(np.linalg.norm(A @ delta - b) / max(1e-300, np.linalg.norm(b))),
                    "solver": str(meta.get("directional_reduced_solver", meta.get("directional_sketch_solver"))),
                }
            )

        full_result: dict[str, Any] | None = None
        if n <= 512:
            configure(mod, "gemm-ns", directional=False, directional_diff="central", sketch_solver="svd", sketch_min_n=1)

            def run_full() -> tuple[Any, dict[str, Any], int]:
                target = LinearTarget(A, b)
                y0 = np.zeros(n, dtype=np.complex128)
                f0 = target.eval(y0)
                delta_full, meta_full = mod.hypercube_delta(target, y0, f0, ep, seed, 0, 0.0, 0.0, 1e-12, 1e12)
                return delta_full, meta_full, target.evals

            d_full, m_full, e_full = run_full()
            t_full = bench(run_full, repeats=3 if n <= 256 else 2, warmup=1)
            full_result = {
                "full_ms": 1e3 * t_full,
                "full_nodes": int(m_full.get("hypercube_nodes", 0)),
                "full_evals": int(e_full),
                "full_rel_step_error": float(np.linalg.norm(d_full - x_true) / max(1e-300, np.linalg.norm(x_true))),
                "full_linear_residual": float(np.linalg.norm(A @ d_full - b) / max(1e-300, np.linalg.norm(b))),
            }

        device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))
        aa = torch.randn((n, n), dtype=torch.float32, device=device)
        bb = torch.randn((n, n), dtype=torch.float32, device=device)
        t_matmul = bench(lambda: torch.matmul(aa, bb), repeats=8 if n <= 512 else 5, warmup=3)

        row: dict[str, Any] = {
            "n": n,
            "k": k,
            "full_hypercube_nodes_nominal": int(max(2 * n + 4, 16)),
            "variants": variants,
            "best_directional_ms": min(v["ms"] for v in variants),
            "forward_vs_central_node_reduction": (
                next(v["nodes"] for v in variants if v["requested_diff"] == "central")
                / max(1, next(v["nodes"] for v in variants if v["requested_diff"] == "forward"))
            ),
            "matmul_float32_ms": 1e3 * t_matmul,
            "matmul_float32_gflops": (2.0 * n**3) / t_matmul / 1e9,
        }
        if full_result is not None:
            row.update(full_result)
            forward = next(v for v in variants if v["requested_diff"] == "forward")
            row["forward_vs_full_speedup"] = full_result["full_ms"] / forward["ms"] if forward["ms"] > 0 else None
            row["forward_oracle_reduction"] = full_result["full_evals"] / max(1, forward["evals"])
        rows.append(row)
    return rows


def ai_matmul_bench() -> dict[str, Any]:
    device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))
    out: dict[str, Any] = {"device": str(device), "matmul": [], "bmm": []}
    for dtype in [torch.float32, torch.float16]:
        if device.type == "cpu" and dtype == torch.float16:
            continue
        for n in [512, 1024, 2048, 4096]:
            a = torch.randn((n, n), dtype=dtype, device=device)
            b = torch.randn((n, n), dtype=dtype, device=device)
            try:
                t = bench(lambda: torch.matmul(a, b), repeats=5 if n <= 2048 else 3, warmup=3)
                out["matmul"].append({"dtype": str(dtype).replace("torch.", ""), "n": n, "ms": 1e3 * t, "gflops": (2.0 * n**3) / t / 1e9})
            except Exception as exc:
                out["matmul"].append({"dtype": str(dtype).replace("torch.", ""), "n": n, "error": f"{type(exc).__name__}:{exc}"})
    for dtype in [torch.float32, torch.float16]:
        if device.type == "cpu" and dtype == torch.float16:
            continue
        for batch, n in [(64, 128), (128, 128), (256, 128), (64, 256), (128, 256)]:
            a = torch.randn((batch, n, n), dtype=dtype, device=device)
            b = torch.randn((batch, n, n), dtype=dtype, device=device)
            try:
                t = bench(lambda: torch.bmm(a, b), repeats=8, warmup=3)
                out["bmm"].append({"dtype": str(dtype).replace("torch.", ""), "batch": batch, "n": n, "ms": 1e3 * t, "gflops": (2.0 * batch * n**3) / t / 1e9})
            except Exception as exc:
                out["bmm"].append({"dtype": str(dtype).replace("torch.", ""), "batch": batch, "n": n, "error": f"{type(exc).__name__}:{exc}"})
    return out


def run_engine_case(args: list[str]) -> dict[str, Any]:
    cmd = [sys.executable, str(ENGINE), *args]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True, check=False)
    out_path = args[args.index("--out") + 1] if "--out" in args else None
    summary = None
    params = None
    if out_path and Path(out_path).exists():
        data = json.loads(Path(out_path).read_text(encoding="utf-8"))
        summary = data.get("summary")
        params = data.get("parameters")
    return {
        "cmd": cmd,
        "returncode": proc.returncode,
        "wall_seconds": time.perf_counter() - t0,
        "stdout_tail": proc.stdout.splitlines()[-8:],
        "stderr_tail": proc.stderr.splitlines()[-8:],
        "summary": summary,
        "parameters": params,
    }


def end_to_end_bench() -> list[dict[str, Any]]:
    cases = [
        ["--self-test", "--out", "/private/tmp/408_bench_self.json"],
        ["--system-source", "poly", "--cases", "2,2", "--polys", "x^2 - 1; y^2 - 4", "--variables", "x,y", "--starts", "-2,-3; -2,3; 2,-3; 2,3", "--count", "4", "--pool", "4", "--accept", "1e-8", "--out", "/private/tmp/408_bench_poly2d.json"],
        ["--cases", "2,4", "--count", "2", "--pool", "32", "--accept", "1e-8", "--out", "/private/tmp/408_bench_ks2x4.json"],
        ["--cases", "2,4", "--count", "2", "--pool", "32", "--accept", "1e-8", "--matrix-algorithm", "gemm-ns", "--batch-kernel", "complex", "--no-directional-jet", "--out", "/private/tmp/408_bench_ks2x4_full.json"],
    ]
    return [{"args": args, "run": run_engine_case(args)} for args in cases]


def main() -> int:
    mod = load_engine()
    result = {
        "engine": str(ENGINE),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "torch": torch.__version__,
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "torch_threads": int(torch.get_num_threads()),
        "cpu_count": os.cpu_count(),
        "directional_variants": directional_variant_bench(mod),
        "ai_matmul": ai_matmul_bench(),
        "end_to_end": end_to_end_bench(),
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
