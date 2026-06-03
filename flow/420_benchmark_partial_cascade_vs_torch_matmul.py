#!/usr/bin/env python3
"""Benchmark 420 when q fails but the cascade stops before n.

Families:
  - q_active: root lies in the first q basis.
  - partial_stage2: root lies in the second cascade basis, usually about 2q.
  - partial_stage3: root lies in the third cascade basis, usually about 4q.
  - random: root is dense and forces full-n IRP.

This isolates the useful middle case: Pandrosion IRP does not win at q, but it
still resolves in a small subspace instead of going to n.
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
OUT = Path("/private/tmp/420_partial_cascade_vs_torch_matmul_benchmark.json")


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
    spec = importlib.util.spec_from_file_location("pandrosion420_partial", ENGINE)
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


def root_in_stage(mod: Any, n: int, cfg: Any, stage: int, seed: int) -> np.ndarray:
    dims = mod.cascade_dims(int(n), cfg)
    dim = int(dims[min(max(1, int(stage)), len(dims)) - 1])
    dmat = mod.basis(int(n), dim, int(cfg.seed) + 31 * int(stage))
    rng = np.random.default_rng(mod.splitmix64(int(seed) + 0x420CA5CADE) & 0xFFFFFFFF)
    coeff = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
    coeff = coeff / max(1e-300, float(np.linalg.norm(coeff)))
    return (dmat @ coeff).astype(np.complex128, copy=False)


def make_target(mod: Any, n: int, cfg: Any, family: str, seed: int) -> tuple[Any, int]:
    dims = mod.cascade_dims(int(n), cfg)
    if family == "q_active":
        return mod.IdentityTarget(root_in_stage(mod, n, cfg, 1, seed)), int(dims[0])
    if family == "partial_stage2":
        return mod.IdentityTarget(root_in_stage(mod, n, cfg, 2, seed)), int(dims[min(1, len(dims) - 1)])
    if family == "partial_stage3":
        return mod.IdentityTarget(root_in_stage(mod, n, cfg, 3, seed)), int(dims[min(2, len(dims) - 1)])
    if family == "random":
        rng = np.random.default_rng(mod.splitmix64(int(seed) + 0x420BAD) & 0xFFFFFFFF)
        root = rng.standard_normal(int(n)) + 1j * rng.standard_normal(int(n))
        root = root / max(1e-300, float(np.linalg.norm(root)))
        return mod.IdentityTarget(root.astype(np.complex128)), int(n)
    raise ValueError(f"unknown family {family!r}")


def run_420(mod: Any, n: int, seed: int, family: str) -> dict[str, Any]:
    cfg = make_cfg(mod, seed)
    target, expected_dim = make_target(mod, n, cfg, family, seed)
    y0 = np.zeros(int(n), dtype=np.complex128)
    repeats = 9 if n <= 1024 else 7
    ms, rec = bench(lambda: mod.corrector(target, y0, cfg), repeats=repeats, warmup=3)
    y = np.asarray(rec["y"], dtype=np.complex128)
    root = np.asarray(target.root, dtype=np.complex128)
    err = float(mod.finite_norm(y - root) / max(1e-300, mod.finite_norm(root)))
    return {
        "ms_420": float(ms),
        "accepted": bool(rec["accepted"]),
        "residual": float(rec["residual"]),
        "relative_error": float(err),
        "accepted_dim": int(rec.get("last_dim", n) or n),
        "expected_dim": int(expected_dim),
        "tried_dims": rec.get("last_tried_dims"),
        "stage": int(rec.get("last_stage", 0) or 0),
        "no_external_fallback": bool(rec.get("no_external_fallback", False)),
    }


def matmul_baselines(n: int, dev: torch.device) -> dict[str, float]:
    out: dict[str, float] = {}
    dtypes = [torch.float16, torch.float32] if dev.type != "cpu" else [torch.float32]
    for dtype in dtypes:
        name = str(dtype).replace("torch.", "")
        repeats = 16 if n <= 512 else (10 if n <= 2048 else 5)
        warmup = 5 if n <= 2048 else 3
        gen = torch.Generator(device="cpu")
        gen.manual_seed(420500 + int(n))
        a = torch.randn((int(n), int(n)), generator=gen, dtype=torch.float32).to(device=dev, dtype=dtype)
        b = torch.randn((int(n), int(n)), generator=gen, dtype=torch.float32).to(device=dev, dtype=dtype)
        x = torch.randn((int(n),), generator=gen, dtype=torch.float32).to(device=dev, dtype=dtype)
        out[f"square_matmul_{name}_ms"] = bench(lambda: torch.matmul(a, b), dev, repeats, warmup)[0]
        out[f"matrix_vector_{name}_ms"] = bench(lambda: torch.matmul(a, x), dev, max(20, repeats), warmup)[0]
    return out


def summarize(family: str, n: int, samples: list[dict[str, Any]], base: dict[str, float]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "family": family,
        "n": int(n),
        "q": int(samples[0]["q"]),
        "k": int(samples[0]["k"]),
        "seeds": len(samples),
        "ms_420": stats([float(s["ms_420"]) for s in samples]),
        "accepted_dim": stats([float(s["accepted_dim"]) for s in samples]),
        "expected_dim": stats([float(s["expected_dim"]) for s in samples]),
        "stage": stats([float(s["stage"]) for s in samples]),
        "relative_error": stats([float(s["relative_error"]) for s in samples]),
        "accepted_rate": float(np.mean([bool(s["accepted"]) for s in samples])),
        "no_external_fallback_rate": float(np.mean([bool(s["no_external_fallback"]) for s in samples])),
        "samples": samples,
    }
    for key, val in base.items():
        row[key] = {"mean": float(val), "std": 0.0, "min": float(val), "median": float(val), "max": float(val)}
    for dtype in ("float16", "float32"):
        for kind in ("square_matmul", "matrix_vector"):
            key = f"{kind}_{dtype}_ms"
            if key in base:
                row[f"speedup_vs_{kind}_{dtype}"] = stats([float(base[key]) / max(1e-300, float(s["ms_420"])) for s in samples])
    return row


def main() -> int:
    mod = load_engine()
    dev = device()
    dims = [128, 512, 1024, 2048, 4096]
    seeds = [42011, 42012, 42013]
    families = ["q_active", "partial_stage2", "partial_stage3", "random"]
    baselines = {n: matmul_baselines(n, dev) for n in dims}
    rows: list[dict[str, Any]] = []
    for family in families:
        for n in dims:
            samples: list[dict[str, Any]] = []
            for seed in seeds:
                cfg = make_cfg(mod, seed)
                rec = run_420(mod, n, seed + 1009 * n, family)
                rec.update({"seed": int(seed), "q": int(mod.coded_dim(n, mod.sketch_dim(n, cfg), cfg)), "k": int(mod.sketch_dim(n, cfg))})
                samples.append(rec)
            row = summarize(family, n, samples, baselines[n])
            rows.append(row)
            print(
                f"{family} n={n} q/k={row['q']}/{row['k']} "
                f"expected={row['expected_dim']['mean']:.1f} accepted={row['accepted_dim']['mean']:.1f} "
                f"420={row['ms_420']['mean']:.4f}ms "
                f"fp16_gemm={row.get('speedup_vs_square_matmul_float16', {}).get('mean', 0.0):.2f}x "
                f"fp16_mv={row.get('speedup_vs_matrix_vector_float16', {}).get('mean', 0.0):.2f}x "
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
        "families": families,
        "interpretation": {
            "partial_stage2": "q fails, but the second cascade basis resolves the target.",
            "partial_stage3": "q and stage2 fail, but the third cascade basis resolves the target.",
            "random": "negative control: the no-fallback cascade reaches full n.",
            "claim_scope": "This tests the middle regime where Pandrosion IRP uses a small part of n, not full n.",
        },
        "rows": rows,
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"out={OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
