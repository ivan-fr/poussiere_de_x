#!/usr/bin/env python3
"""Benchmark 421 random full-n IRP against torch.matmul.

This is the negative-control / stress case for the q cascade:

  - roots are random dense vectors, so the cascade reaches full n;
  - 421 remains no-fallback and uses its fast full-n identity IRP stage;
  - baselines are torch.matmul(A, B) square GEMM and torch.matmul(A, x) matvec.

The operations are not mathematically identical.  The benchmark measures latency
against dense matmul-shaped kernels to quantify the full-n 421 regime.
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
ENGINE = ROOT / "flow" / "421_pandrosion_standalone_full_irp_fast_fulln_engine.py"
OUT = Path("/private/tmp/421_random_fulln_vs_torch_matmul_benchmark.json")


def sync(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize()
    elif device.type == "mps":
        torch.mps.synchronize()


def bench(fn: Any, *, device: torch.device | None = None, repeats: int = 9, warmup: int = 3) -> tuple[float, Any]:
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
    spec = importlib.util.spec_from_file_location("pandrosion421_random_fulln", ENGINE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {ENGINE}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def cfg(mod: Any, seed: int) -> Any:
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


def random_root(mod: Any, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(mod.splitmix64(int(seed) + 1009 * int(n) + 0x421BAD) & 0xFFFFFFFF)
    root = rng.standard_normal(int(n)) + 1j * rng.standard_normal(int(n))
    return (root / max(1e-300, float(np.linalg.norm(root)))).astype(np.complex128)


def run_421(mod: Any, n: int, seed: int) -> dict[str, Any]:
    c = cfg(mod, seed)
    root = random_root(mod, n, seed)
    target = mod.IdentityTarget(root)
    y0 = np.zeros(int(n), dtype=np.complex128)
    repeats = 13 if n <= 1024 else 9
    ms, rec = bench(lambda: mod.corrector(target, y0, c), repeats=repeats, warmup=4)
    y = np.asarray(rec["y"], dtype=np.complex128)
    err = float(mod.finite_norm(y - root) / max(1e-300, mod.finite_norm(root)))
    return {
        "seed": int(seed),
        "ms_421": float(ms),
        "accepted": bool(rec["accepted"]),
        "residual": float(rec["residual"]),
        "relative_error": float(err),
        "accepted_dim": int(rec.get("last_dim", n) or n),
        "tried_dims": rec.get("last_tried_dims"),
        "stage": int(rec.get("last_stage", 0) or 0),
        "solver": rec.get("last_solver"),
        "fast_full_n_identity": bool(rec.get("last_fast_full_n_identity", False)),
        "no_external_fallback": bool(rec.get("no_external_fallback", False)),
    }


def matmul_baselines(n: int, dev: torch.device) -> dict[str, float]:
    out: dict[str, float] = {}
    dtypes = [torch.float16, torch.float32] if dev.type != "cpu" else [torch.float32]
    for dtype in dtypes:
        name = str(dtype).replace("torch.", "")
        repeats = 18 if n <= 512 else (12 if n <= 2048 else 6)
        warmup = 6 if n <= 2048 else 4
        gen = torch.Generator(device="cpu")
        gen.manual_seed(421000 + int(n))
        a = torch.randn((int(n), int(n)), generator=gen, dtype=torch.float32).to(device=dev, dtype=dtype)
        b = torch.randn((int(n), int(n)), generator=gen, dtype=torch.float32).to(device=dev, dtype=dtype)
        x = torch.randn((int(n),), generator=gen, dtype=torch.float32).to(device=dev, dtype=dtype)
        gemm_ms, _ = bench(lambda: torch.matmul(a, b), device=dev, repeats=repeats, warmup=warmup)
        matvec_ms, _ = bench(lambda: torch.matmul(a, x), device=dev, repeats=max(24, repeats), warmup=warmup)
        out[f"square_matmul_{name}_ms"] = float(gemm_ms)
        out[f"square_matmul_{name}_gflops"] = float((2.0 * int(n) ** 3) / (gemm_ms / 1e3) / 1e9)
        out[f"matrix_vector_{name}_ms"] = float(matvec_ms)
        out[f"matrix_vector_{name}_gflops"] = float((2.0 * int(n) ** 2) / (matvec_ms / 1e3) / 1e9)
    return out


def summarize(mod: Any, n: int, samples: list[dict[str, Any]], baselines: dict[str, float]) -> dict[str, Any]:
    c = cfg(mod, 42101)
    row: dict[str, Any] = {
        "family": "random_full_n",
        "n": int(n),
        "q": int(mod.coded_dim(n, mod.sketch_dim(n, c), c)),
        "k": int(mod.sketch_dim(n, c)),
        "seeds": len(samples),
        "ms_421": stats([float(s["ms_421"]) for s in samples]),
        "accepted_dim": stats([float(s["accepted_dim"]) for s in samples]),
        "stage": stats([float(s["stage"]) for s in samples]),
        "accepted_rate": float(np.mean([bool(s["accepted"]) for s in samples])),
        "fast_full_n_identity_rate": float(np.mean([bool(s["fast_full_n_identity"]) for s in samples])),
        "no_external_fallback_rate": float(np.mean([bool(s["no_external_fallback"]) for s in samples])),
        "relative_error": stats([float(s["relative_error"]) for s in samples]),
        "residual": stats([float(s["residual"]) for s in samples]),
        "samples": samples,
    }
    for key, val in baselines.items():
        row[key] = {"mean": float(val), "std": 0.0, "min": float(val), "median": float(val), "max": float(val)}
    for dtype in ("float16", "float32"):
        for kind in ("square_matmul", "matrix_vector"):
            key = f"{kind}_{dtype}_ms"
            if key in baselines:
                row[f"speedup_vs_{kind}_{dtype}"] = stats([float(baselines[key]) / max(1e-300, float(s["ms_421"])) for s in samples])
    return row


def main() -> int:
    mod = load_engine()
    dev = device()
    dims = [128, 512, 1024, 2048, 4096]
    seeds = [42101, 42102, 42103, 42104, 42105]
    rows: list[dict[str, Any]] = []
    for n in dims:
        baselines = matmul_baselines(n, dev)
        samples = [run_421(mod, n, seed + 1009 * n) for seed in seeds]
        row = summarize(mod, n, samples, baselines)
        rows.append(row)
        print(
            f"random_full_n n={n} q/k={row['q']}/{row['k']} "
            f"421={row['ms_421']['mean']:.4f}ms accepted_dim={row['accepted_dim']['mean']:.1f} "
            f"fp16_gemm={row.get('speedup_vs_square_matmul_float16', {}).get('mean', 0.0):.2f}x "
            f"fp16_mv={row.get('speedup_vs_matrix_vector_float16', {}).get('mean', 0.0):.2f}x "
            f"fp32_gemm={row.get('speedup_vs_square_matmul_float32', {}).get('mean', 0.0):.2f}x "
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
            "random_full_n": "Dense random identity roots force the cascade to the final full-n IRP stage.",
            "421_fast_full_n": "For validated F(z)=z-root targets, 421 uses the exact full-n IRP step delta=-F(y), not an external matmul fallback.",
            "square_matmul": "Dense n-by-n torch.matmul(A, B), the common AI GEMM proxy.",
            "matrix_vector": "Dense torch.matmul(A, x), closer to vector correction shape.",
            "claim_scope": "This is a latency comparison, not an equivalence between vector IRP correction and dense GEMM.",
        },
        "rows": rows,
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"out={OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
