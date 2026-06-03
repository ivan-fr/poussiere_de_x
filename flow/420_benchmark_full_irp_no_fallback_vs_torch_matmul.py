#!/usr/bin/env python3
"""Benchmark Pandrosion 420 full IRP against torch.matmul baselines.

The comparison is intentionally split:

  - q_active identity roots: 420 should stop in the small q stage.
  - random identity roots: 420 has no fallback and must go to full-n IRP.
  - square GEMM A @ B: dense AI matmul throughput proxy.
  - matrix-vector A @ x: closer to the vector correction shape of 420.

This does not claim that a root-correction and a dense GEMM are the same
operation.  It quantifies when 420's IRP-only cascade is faster or slower than
torch.matmul-shaped dense kernels.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "flow" / "420_pandrosion_standalone_full_irp_no_fallback_engine.py"
OUT = Path("/private/tmp/420_full_irp_no_fallback_vs_torch_matmul_benchmark.json")


def sync(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize()
    elif device.type == "mps":
        torch.mps.synchronize()


def bench(fn: Any, device: torch.device | None = None, repeats: int = 7, warmup: int = 3) -> tuple[float, Any]:
    last: Any = None
    for _ in range(max(0, warmup)):
        last = fn()
    if device is not None:
        sync(device)
    vals: list[float] = []
    for _ in range(max(1, repeats)):
        if device is not None:
            sync(device)
        t0 = time.perf_counter()
        last = fn()
        if device is not None:
            sync(device)
        vals.append(time.perf_counter() - t0)
    return 1e3 * min(vals), last


def stats(vals: list[float]) -> dict[str, float]:
    a = np.asarray(vals, dtype=float)
    return {
        "mean": float(np.mean(a)),
        "std": float(np.std(a, ddof=1)) if a.size > 1 else 0.0,
        "min": float(np.min(a)),
        "median": float(np.median(a)),
        "max": float(np.max(a)),
    }


def device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def load_engine() -> Any:
    spec = importlib.util.spec_from_file_location("pandrosion420_bench", ENGINE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {ENGINE}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def make_cfg(mod: Any, seed: int) -> Any:
    return mod.Config(
        seed=int(seed),
        sketch_factor=2.75,
        coded_factor=2.0,
        coded_min=8,
        growth=2.0,
        epochs=2,
        tol=1e-12,
        accept=1e-9,
        linear_accept=0.15,
        line_grid=(1.0,),
    )


def run_420_sample(mod: Any, n: int, seed: int, family: str) -> dict[str, Any]:
    cfg = make_cfg(mod, seed)
    q_aligned = family == "q_active"
    target, y0 = mod.make_identity(int(n), int(seed), q_aligned=q_aligned, cfg=cfg)
    ms, rec = bench(lambda: mod.corrector(target, y0, cfg), device=None, repeats=7 if n <= 1024 else 5, warmup=2)
    y = np.asarray(rec["y"], dtype=np.complex128)
    root = np.asarray(target.root, dtype=np.complex128)
    err = float(mod.finite_norm(y - root) / max(1e-300, mod.finite_norm(root)))
    return {
        "ms_420": float(ms),
        "accepted": bool(rec["accepted"]),
        "residual": float(rec["residual"]),
        "relative_error": float(err),
        "accepted_dim": int(rec.get("last_dim", n) or n),
        "tried_dims": rec.get("last_tried_dims"),
        "stage": int(rec.get("last_stage", 0) or 0),
        "no_external_fallback": bool(rec.get("no_external_fallback", False)),
    }


def matmul_baselines(n: int, dev: torch.device) -> dict[str, Any]:
    out: dict[str, Any] = {}
    dtypes = [torch.float16, torch.float32] if dev.type != "cpu" else [torch.float32]
    for dtype in dtypes:
        name = str(dtype).replace("torch.", "")
        repeats = 16 if n <= 512 else (10 if n <= 2048 else 5)
        warmup = 5 if n <= 2048 else 3
        gen = torch.Generator(device="cpu")
        gen.manual_seed(420000 + int(n))
        a = torch.randn((int(n), int(n)), generator=gen, dtype=torch.float32).to(device=dev, dtype=dtype)
        b = torch.randn((int(n), int(n)), generator=gen, dtype=torch.float32).to(device=dev, dtype=dtype)
        x = torch.randn((int(n),), generator=gen, dtype=torch.float32).to(device=dev, dtype=dtype)
        gemm_ms, _ = bench(lambda: torch.matmul(a, b), device=dev, repeats=repeats, warmup=warmup)
        matvec_ms, _ = bench(lambda: torch.matmul(a, x), device=dev, repeats=max(20, repeats), warmup=warmup)
        out[f"square_matmul_{name}_ms"] = float(gemm_ms)
        out[f"square_matmul_{name}_gflops"] = float((2.0 * int(n) ** 3) / (gemm_ms / 1e3) / 1e9)
        out[f"matrix_vector_{name}_ms"] = float(matvec_ms)
        out[f"matrix_vector_{name}_gflops"] = float((2.0 * int(n) ** 2) / (matvec_ms / 1e3) / 1e9)
    return out


def summarize(family: str, n: int, samples: list[dict[str, Any]], baselines: dict[str, Any]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "family": family,
        "n": int(n),
        "seeds": len(samples),
        "q": int(samples[0]["q"]),
        "k": int(samples[0]["k"]),
        "ms_420": stats([float(s["ms_420"]) for s in samples]),
        "accepted_rate": float(np.mean([bool(s["accepted"]) for s in samples])),
        "relative_error": stats([float(s["relative_error"]) for s in samples]),
        "residual": stats([float(s["residual"]) for s in samples]),
        "accepted_dim": stats([float(s["accepted_dim"]) for s in samples]),
        "stage": stats([float(s["stage"]) for s in samples]),
        "no_external_fallback_rate": float(np.mean([bool(s["no_external_fallback"]) for s in samples])),
        "samples": samples,
    }
    for key, val in baselines.items():
        row[key] = {"mean": float(val), "std": 0.0, "min": float(val), "median": float(val), "max": float(val)}
    for dtype in ("float16", "float32"):
        gkey = f"square_matmul_{dtype}_ms"
        vkey = f"matrix_vector_{dtype}_ms"
        if gkey in baselines:
            row[f"speedup_vs_square_matmul_{dtype}"] = stats([float(baselines[gkey]) / max(1e-300, float(s["ms_420"])) for s in samples])
        if vkey in baselines:
            row[f"speedup_vs_matrix_vector_{dtype}"] = stats([float(baselines[vkey]) / max(1e-300, float(s["ms_420"])) for s in samples])
    return row


def main() -> int:
    mod = load_engine()
    dev = device()
    dims = [128, 512, 1024, 2048, 4096]
    seeds = [42001, 42002, 42003]
    baselines_by_n: dict[int, dict[str, Any]] = {}
    rows: list[dict[str, Any]] = []
    for n in dims:
        baselines_by_n[n] = matmul_baselines(n, dev)
    for family in ("q_active", "random"):
        for n in dims:
            samples: list[dict[str, Any]] = []
            for seed in seeds:
                cfg = make_cfg(mod, seed)
                rec = run_420_sample(mod, n, seed + 1009 * n, family)
                rec.update({"seed": int(seed), "q": int(mod.coded_dim(n, mod.sketch_dim(n, cfg), cfg)), "k": int(mod.sketch_dim(n, cfg))})
                samples.append(rec)
            row = summarize(family, n, samples, baselines_by_n[n])
            rows.append(row)
            print(
                f"{family} n={n} q/k={row['q']}/{row['k']} 420={row['ms_420']['mean']:.4f}ms "
                f"accepted_dim={row['accepted_dim']['mean']:.1f} "
                f"fp16_gemm_speed={row.get('speedup_vs_square_matmul_float16', {}).get('mean', 0.0):.2f}x "
                f"fp16_mv_speed={row.get('speedup_vs_matrix_vector_float16', {}).get('mean', 0.0):.2f}x "
                f"err_max={row['relative_error']['max']:.2e}",
                flush=True,
            )
    result = {
        "engine": str(ENGINE),
        "output": str(OUT),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cpu_count": os.cpu_count(),
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "device": str(dev),
        "dims": dims,
        "seeds": seeds,
        "interpretation": {
            "q_active": "420 solves F(z)=z-root and stops in the small q IRP stage.",
            "random": "420 has no fallback and therefore grows to full-n IRP.",
            "square_matmul": "Dense n-by-n torch.matmul(A, B), the common AI GEMM proxy.",
            "matrix_vector": "Dense torch.matmul(A, x), closer to the vector shape corrected by 420.",
            "claim_scope": "420 is not a dense GEMM replacement. It wins only when the IRP cascade stops early; the no-fallback random case is slower.",
        },
        "rows": rows,
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"out={OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
