#!/usr/bin/env python3
"""Benchmark 405 Pandrosion matrix kernels against AI-like Torch matmul kernels."""
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
ENGINE = ROOT / "flow" / "405_pandrosion_standalone_local_jet_geometry_realpack_batched_irp_gemm_engine.py"
OUT = Path("/private/tmp/405_pandrosion_benchmark_matrix_vs_ai.json")


def sync() -> None:
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    if torch.backends.mps.is_available():
        torch.mps.synchronize()


def bench(fn: Any, repeats: int = 10, warmup: int = 3) -> float:
    for _ in range(max(0, warmup)):
        fn()
    sync()
    t0 = time.perf_counter()
    for _ in range(max(1, repeats)):
        fn()
    sync()
    return (time.perf_counter() - t0) / max(1, repeats)


def load_engine() -> Any:
    spec = importlib.util.spec_from_file_location("pandrosion405", ENGINE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {ENGINE}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pandrosion405"] = mod
    spec.loader.exec_module(mod)
    return mod


def configure(mod: Any, algorithm: str, batch_kernel: str = "auto", real_dtype: str = "auto", ns_iters: int = 8) -> None:
    args = SimpleNamespace(
        matrix_backend="torch",
        matrix_algorithm=algorithm,
        batch_kernel=batch_kernel,
        torch_device="auto",
        torch_complex_dtype="auto",
        torch_real_dtype=real_dtype,
        ns_iters=ns_iters,
        ns_damping_floor=1e-5,
    )
    mod.configure_matrix_backend(args)


def kernel_quality(A: Any, b: Any, x: Any) -> float:
    return float(np.linalg.norm(np.asarray(A) @ np.asarray(x) - np.asarray(b)) / max(1e-300, np.linalg.norm(b)))


def matrix_kernel_bench(mod: Any) -> list[dict[str, Any]]:
    rng = np.random.default_rng(405)
    rows: list[dict[str, Any]] = []
    for n in [16, 32, 64, 128, 256, 512]:
        A = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / math.sqrt(n)
        b = rng.normal(size=n) + 1j * rng.normal(size=n)
        repeats = max(3, 512 // n)

        configure(mod, "adaptive-ns")
        x_adapt, meta_adapt = mod.matrix_svd_step(A, b)
        t_adapt = bench(lambda: mod.matrix_svd_step(A, b), repeats=repeats)

        configure(mod, "gemm-ns")
        x_complex, _meta_complex = mod.matrix_svd_step(A, b)
        t_complex = bench(lambda: mod.matrix_svd_step(A, b), repeats=repeats)

        configure(mod, "realpack-ns", batch_kernel="realpack")
        x_real, _meta_real = mod.matrix_svd_step(A, b)
        t_real = bench(lambda: mod.matrix_svd_step(A, b), repeats=repeats)

        configure(mod, "svd")
        x_svd, _meta_svd = mod.matrix_svd_step(A, b)
        t_svd = bench(lambda: mod.matrix_svd_step(A, b), repeats=max(2, repeats // 2))

        rows.append({
            "n": n,
            "repeats": repeats,
            "adaptive_resolved": meta_adapt.get("matrix_algorithm", meta_adapt.get("matrix_method")),
            "adaptive_ms": 1e3 * t_adapt,
            "complex_404style_ms": 1e3 * t_complex,
            "realpack_405_ms": 1e3 * t_real,
            "torch_svd_ms": 1e3 * t_svd,
            "adaptive_rel_residual": kernel_quality(A, b, x_adapt),
            "complex_rel_residual": kernel_quality(A, b, x_complex),
            "realpack_rel_residual": kernel_quality(A, b, x_real),
            "svd_rel_residual": kernel_quality(A, b, x_svd),
        })
    return rows


def prepass_bench(mod: Any) -> list[dict[str, Any]]:
    rng = np.random.default_rng(406)
    lambdas = np.asarray([1.0, 0.75, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03125, 0.02], dtype=float)
    rows: list[dict[str, Any]] = []
    for batch, n in [(8, 16), (32, 16), (64, 32), (128, 32), (64, 64), (64, 128)]:
        m = max(2 * n + 4, 16)
        y = (rng.normal(size=(batch, n)) + 1j * rng.normal(size=(batch, n))) / math.sqrt(n)
        dy = (rng.normal(size=(batch, m, n)) + 1j * rng.normal(size=(batch, m, n))) * 1e-5
        # Make dF locally consistent enough to avoid benchmarking only bad conditioning.
        jac = (rng.normal(size=(batch, n, n)) + 1j * rng.normal(size=(batch, n, n))) / math.sqrt(n)
        df = np.einsum("bmn,bkn->bmk", dy, jac)
        f = (rng.normal(size=(batch, n)) + 1j * rng.normal(size=(batch, n))) / math.sqrt(n)
        yn = np.maximum(1.0, np.linalg.norm(y, axis=1))
        repeats = max(3, 256 // max(1, batch))

        configure(mod, "gemm-ns", batch_kernel="complex")
        t_complex = bench(lambda: mod.torch_batched_gemm_candidates_405(dy, df, f, y, lambdas, 0.0, 0.0, yn), repeats=repeats)

        configure(mod, "realpack-ns", batch_kernel="realpack")
        t_real = bench(lambda: mod.torch_batched_realpack_candidates_405(dy, df, f, y, lambdas, 0.0, 0.0, yn), repeats=repeats)

        configure(mod, "adaptive-ns", batch_kernel="auto")
        kernel = mod.selected_batch_kernel_405()
        t_adapt = bench(
            lambda: (
                mod.torch_batched_realpack_candidates_405(dy, df, f, y, lambdas, 0.0, 0.0, yn)
                if kernel == "realpack"
                else mod.torch_batched_gemm_candidates_405(dy, df, f, y, lambdas, 0.0, 0.0, yn)
            ),
            repeats=repeats,
        )

        real_ops = float(batch) * (8.0 * m * n * n + 2.0 * 8 * n**3 + 8.0 * n**3)
        rows.append({
            "batch": batch,
            "n": n,
            "jet_samples": m,
            "repeats": repeats,
            "adaptive_kernel": kernel,
            "adaptive_ms": 1e3 * t_adapt,
            "complex_404style_ms": 1e3 * t_complex,
            "realpack_405_ms": 1e3 * t_real,
            "realpack_vs_complex_speedup": t_complex / t_real if t_real > 0 else None,
            "approx_realpack_gflops": real_ops / t_real / 1e9 if t_real > 0 else None,
        })
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
            repeats = 5 if n <= 2048 else 3
            try:
                t = bench(lambda: torch.matmul(a, b), repeats=repeats, warmup=3)
                out["matmul"].append({
                    "dtype": str(dtype).replace("torch.", ""),
                    "n": n,
                    "ms": 1e3 * t,
                    "gflops": (2.0 * n**3) / t / 1e9,
                })
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
                out["bmm"].append({
                    "dtype": str(dtype).replace("torch.", ""),
                    "batch": batch,
                    "n": n,
                    "ms": 1e3 * t,
                    "gflops": (2.0 * batch * n**3) / t / 1e9,
                })
            except Exception as exc:
                out["bmm"].append({"dtype": str(dtype).replace("torch.", ""), "batch": batch, "n": n, "error": f"{type(exc).__name__}:{exc}"})
    return out


def run_engine_case(args: list[str]) -> dict[str, Any]:
    cmd = [sys.executable, str(ENGINE), *args]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True, check=False)
    return {
        "cmd": cmd,
        "returncode": proc.returncode,
        "wall_seconds": time.perf_counter() - t0,
        "stdout_tail": proc.stdout.splitlines()[-8:],
        "stderr_tail": proc.stderr.splitlines()[-8:],
    }


def end_to_end_bench() -> list[dict[str, Any]]:
    cases = [
        ["--self-test", "--out", "/private/tmp/405_bench_self_adaptive.json"],
        ["--system-source", "poly", "--cases", "2,2", "--polys", "x^2 - 1; y^2 - 4", "--variables", "x,y", "--starts", "-2,-3; -2,3; 2,-3; 2,3", "--count", "4", "--pool", "4", "--accept", "1e-8", "--out", "/private/tmp/405_bench_poly2d_adaptive.json"],
        ["--cases", "2,4", "--count", "2", "--pool", "32", "--accept", "1e-8", "--out", "/private/tmp/405_bench_ks2x4_adaptive.json"],
        ["--cases", "2,4", "--count", "2", "--pool", "32", "--accept", "1e-8", "--matrix-algorithm", "realpack-ns", "--batch-kernel", "realpack", "--out", "/private/tmp/405_bench_ks2x4_realpack.json"],
    ]
    rows: list[dict[str, Any]] = []
    for args in cases:
        run = run_engine_case(args)
        out_path = None
        if "--out" in args:
            out_path = args[args.index("--out") + 1]
        summary = None
        params = None
        if out_path and Path(out_path).exists():
            data = json.loads(Path(out_path).read_text(encoding="utf-8"))
            summary = data.get("summary")
            params = data.get("parameters")
        rows.append({"args": args, "run": run, "summary": summary, "parameters": params})
    return rows


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
        "matrix_kernels": matrix_kernel_bench(mod),
        "prepass_kernels": prepass_bench(mod),
        "ai_matmul": ai_matmul_bench(),
        "end_to_end": end_to_end_bench(),
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
