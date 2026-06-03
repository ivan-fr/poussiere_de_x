#!/usr/bin/env python3
"""Direct benchmark of 421 fast full-n IRP against 420.

Both engines are standalone; this benchmark imports them only to run the same
identity targets with the same seeds.  It measures the q-active case and the
random full-n case where 421 should avoid 420's dense identity projection.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
E420 = ROOT / "flow" / "420_pandrosion_standalone_full_irp_no_fallback_engine.py"
E421 = ROOT / "flow" / "421_pandrosion_standalone_full_irp_fast_fulln_engine.py"
OUT = Path("/private/tmp/421_fast_fulln_vs_420_benchmark.json")


def load(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def bench(fn: Any, repeats: int = 9, warmup: int = 3) -> tuple[float, Any]:
    last: Any = None
    for _ in range(max(0, warmup)):
        last = fn()
    vals: list[float] = []
    for _ in range(max(1, repeats)):
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
    return float(ms), {
        "accepted": bool(rec["accepted"]),
        "residual": float(rec["residual"]),
        "relative_error": err,
        "accepted_dim": int(rec.get("last_dim", root.size) or root.size),
        "tried_dims": rec.get("last_tried_dims"),
        "solver": rec.get("last_solver"),
        "fast_full_n_identity": bool(rec.get("last_fast_full_n_identity", False)),
    }


def main() -> int:
    m420 = load(E420, "p420")
    m421 = load(E421, "p421")
    dims = [128, 512, 1024, 2048, 4096]
    seeds = [42101, 42102, 42103]
    rows: list[dict[str, Any]] = []
    for family in ("q_active", "random"):
        for n in dims:
            samples: list[dict[str, Any]] = []
            repeats = 11 if n <= 1024 else 7
            for seed in seeds:
                root = root_for(m420, n, seed, family)
                ms420, r420 = run_one(m420, root, seed, repeats)
                ms421, r421 = run_one(m421, root, seed, repeats)
                samples.append({
                    "seed": int(seed),
                    "ms420": ms420,
                    "ms421": ms421,
                    "speedup_421_vs_420": ms420 / max(1e-300, ms421),
                    "420": r420,
                    "421": r421,
                })
            c = cfg(m420, seeds[0])
            row = {
                "family": family,
                "n": int(n),
                "q": int(m420.coded_dim(n, m420.sketch_dim(n, c), c)),
                "k": int(m420.sketch_dim(n, c)),
                "ms420": stats([s["ms420"] for s in samples]),
                "ms421": stats([s["ms421"] for s in samples]),
                "speedup_421_vs_420": stats([s["speedup_421_vs_420"] for s in samples]),
                "accepted_dim_421": stats([float(s["421"]["accepted_dim"]) for s in samples]),
                "fast_full_n_identity_rate": float(np.mean([bool(s["421"]["fast_full_n_identity"]) for s in samples])),
                "relative_error_421": stats([float(s["421"]["relative_error"]) for s in samples]),
                "samples": samples,
            }
            rows.append(row)
            print(
                f"{family} n={n} q/k={row['q']}/{row['k']} "
                f"420={row['ms420']['mean']:.4f}ms 421={row['ms421']['mean']:.4f}ms "
                f"speed={row['speedup_421_vs_420']['mean']:.2f}x "
                f"fast_fulln={row['fast_full_n_identity_rate']:.2f} err={row['relative_error_421']['max']:.2e}",
                flush=True,
            )
    result = {
        "engine_420": str(E420),
        "engine_421": str(E421),
        "output": str(OUT),
        "dims": dims,
        "seeds": seeds,
        "rows": rows,
        "interpretation": {
            "q_active": "421 should match 420 because both stop in q.",
            "random": "421 should improve full-n identity cases by avoiding dense identity projection.",
        },
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
