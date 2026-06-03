#!/usr/bin/env python3
"""Benchmark 422 adaptive dense jump against 421 fast full-n."""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
E421 = ROOT / "flow" / "421_pandrosion_standalone_full_irp_fast_fulln_engine.py"
E422 = ROOT / "flow" / "422_pandrosion_standalone_full_irp_adaptive_dense_jump_engine.py"
OUT = Path("/private/tmp/422_adaptive_dense_jump_vs_421_benchmark.json")


def load(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def bench(fn: Any, repeats: int = 11, warmup: int = 3) -> tuple[float, Any]:
    last: Any = None
    for _ in range(warmup):
        last = fn()
    vals: list[float] = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        last = fn()
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


def cfg(mod: Any, seed: int) -> Any:
    return mod.Config(seed=int(seed), epochs=2, line_grid=(1.0,), tol=1e-12, accept=1e-9)


def root_for(mod: Any, n: int, seed: int, family: str) -> np.ndarray:
    c = cfg(mod, seed)
    rng = np.random.default_rng(mod.splitmix64(seed + 1009 * n) & 0xFFFFFFFF)
    if family == "q_active":
        q = mod.coded_dim(n, mod.sketch_dim(n, c), c)
        dmat = mod.basis(n, q, int(c.seed) + 31)
        coeff = rng.standard_normal(q) + 1j * rng.standard_normal(q)
        return (dmat @ (coeff / max(1e-300, np.linalg.norm(coeff)))).astype(np.complex128)
    root = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    return (root / max(1e-300, np.linalg.norm(root))).astype(np.complex128)


def run_one(mod: Any, root: np.ndarray, seed: int, repeats: int) -> tuple[float, dict[str, Any]]:
    c = cfg(mod, seed)
    target = mod.IdentityTarget(root)
    y0 = np.zeros(root.size, dtype=np.complex128)
    ms, rec = bench(lambda: mod.corrector(target, y0, c), repeats=repeats, warmup=3)
    y = np.asarray(rec["y"], dtype=np.complex128)
    err = float(mod.finite_norm(y - root) / max(1e-300, mod.finite_norm(root)))
    return ms, {
        "accepted": bool(rec["accepted"]),
        "relative_error": err,
        "accepted_dim": int(rec.get("last_dim", root.size) or root.size),
        "tried_dims": rec.get("last_tried_dims"),
        "solver": rec.get("last_solver"),
        "fast_full_n_identity": bool(rec.get("last_fast_full_n_identity", False)),
        "adaptive_dense_jump": bool(rec.get("last_adaptive_dense_jump", False)),
        "adaptive_dense_jump_from_dim": rec.get("last_adaptive_dense_jump_from_dim"),
    }


def main() -> int:
    m421 = load(E421, "p421")
    m422 = load(E422, "p422")
    dims = [128, 512, 1024, 2048, 4096]
    seeds = [42201, 42202, 42203, 42204, 42205]
    rows: list[dict[str, Any]] = []
    for family in ("q_active", "random"):
        for n in dims:
            samples: list[dict[str, Any]] = []
            repeats = 13 if n <= 1024 else 9
            for seed in seeds:
                root = root_for(m421, n, seed, family)
                ms421, r421 = run_one(m421, root, seed, repeats)
                ms422, r422 = run_one(m422, root, seed, repeats)
                samples.append({
                    "seed": int(seed),
                    "ms421": float(ms421),
                    "ms422": float(ms422),
                    "speedup_422_vs_421": float(ms421 / max(1e-300, ms422)),
                    "421": r421,
                    "422": r422,
                })
            c = cfg(m421, seeds[0])
            row = {
                "family": family,
                "n": int(n),
                "q": int(m421.coded_dim(n, m421.sketch_dim(n, c), c)),
                "k": int(m421.sketch_dim(n, c)),
                "ms421": stats([s["ms421"] for s in samples]),
                "ms422": stats([s["ms422"] for s in samples]),
                "speedup_422_vs_421": stats([s["speedup_422_vs_421"] for s in samples]),
                "accepted_dim_422": stats([float(s["422"]["accepted_dim"]) for s in samples]),
                "adaptive_dense_jump_rate": float(np.mean([bool(s["422"]["adaptive_dense_jump"]) for s in samples])),
                "relative_error_422": stats([float(s["422"]["relative_error"]) for s in samples]),
                "samples": samples,
            }
            rows.append(row)
            print(
                f"{family} n={n} q/k={row['q']}/{row['k']} "
                f"421={row['ms421']['mean']:.4f}ms 422={row['ms422']['mean']:.4f}ms "
                f"speed={row['speedup_422_vs_421']['mean']:.2f}x jump={row['adaptive_dense_jump_rate']:.2f} "
                f"err={row['relative_error_422']['max']:.2e}",
                flush=True,
            )
    OUT.write_text(json.dumps({"engine_421": str(E421), "engine_422": str(E422), "rows": rows}, indent=2), encoding="utf-8")
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
