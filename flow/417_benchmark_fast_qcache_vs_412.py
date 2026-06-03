#!/usr/bin/env python3
"""Benchmark 417 fast-qcache against 412 on validated q-adapter targets."""
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
ENGINE_412 = ROOT / "flow" / "412_pandrosion_standalone_local_jet_geometry_composite_cached_irp_engine.py"
ENGINE_417 = ROOT / "flow" / "417_pandrosion_standalone_local_jet_geometry_fast_qcache_irp_engine.py"
OUT = Path("/private/tmp/417_fast_qcache_vs_412_benchmark.json")


def sync() -> None:
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    if torch.backends.mps.is_available():
        torch.mps.synchronize()


def bench(fn: Any, repeats: int = 7, warmup: int = 3) -> tuple[float, dict[str, Any]]:
    last: dict[str, Any] = {}
    for _ in range(max(0, warmup)):
        last = fn()
    sync()
    vals: list[float] = []
    for _ in range(max(1, repeats)):
        sync()
        t0 = time.perf_counter()
        last = fn()
        sync()
        vals.append(time.perf_counter() - t0)
    return 1e3 * min(vals), last


def load_engine(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def configure(mod: Any, *, fastpath: bool, sketch_seed: int = 412) -> None:
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
            sketch_seed=int(sketch_seed),
            sketch_mode="cached-rademacher",
            sketch_solver="svd",
            sketch_basis_cache=True,
            sketch_basis_cache_max=4096,
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
            directional_identity_target_fastpath=bool(fastpath),
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


def q_dims(mod: Any, n: int) -> tuple[int, int]:
    k = max(1, min(int(n), int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))))
    q = int(mod.directional_coded_dim_for(int(n), int(k)))
    return k, q


def make_q_root(mod: Any, n: int, seed: int, scale: float = 0.25) -> tuple[Any, int, int]:
    k, q = q_dims(mod, n)
    salt = int(seed) + 0x412D1A
    C = mod._coded_composite_basis_np(int(n), k, q, salt=salt)
    rng = np.random.default_rng(int(seed) ^ 0x417A11)
    coeff = rng.standard_normal(q) + 1j * rng.standard_normal(q)
    coeff = float(scale) * coeff / max(1e-300, float(np.linalg.norm(coeff)))
    return C @ coeff.astype(np.complex128, copy=False), k, q


def run_corrector(mod: Any, root: Any, seed: int) -> dict[str, Any]:
    target = VectorTarget(root)
    y0 = np.zeros_like(np.asarray(root, dtype=np.complex128))
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
    identity_fast = bool(rec.get("directional_reduced_identity_target_fastpath", False))
    q_active = bool(
        rec.get("directional_probe_kind") == "coded-probe"
        and (rec.get("hypercube_nodes") == rec.get("directional_coded_dim") or identity_fast)
    )
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
        "identity_fastpath": identity_fast,
        "line_search_evals": int(rec.get("line_search_evals", 0) or 0),
        "early_unit_step": bool(rec.get("line_search_early_unit_step", False)),
        "target_evals": int(target.evals),
        "q_active": q_active,
    }


def run_case(mod412: Any, mod417: Any, n: int, seeds: list[int]) -> dict[str, Any]:
    samples: list[dict[str, Any]] = []
    for base_seed in seeds:
        seed = int(base_seed) + 1009 * int(n)
        configure(mod412, fastpath=False, sketch_seed=412)
        configure(mod417, fastpath=True, sketch_seed=412)
        root, k, q = make_q_root(mod412, n, seed)

        configure(mod412, fastpath=False, sketch_seed=412)
        ms412, rec412 = bench(lambda: run_corrector(mod412, root, seed), repeats=7, warmup=3)

        configure(mod417, fastpath=False, sketch_seed=412)
        ms417_no, rec417_no = bench(lambda: run_corrector(mod417, root, seed), repeats=7, warmup=3)

        configure(mod417, fastpath=True, sketch_seed=412)
        ms417, rec417 = bench(lambda: run_corrector(mod417, root, seed), repeats=9, warmup=4)

        samples.append(
            {
                "seed": seed,
                "n": int(n),
                "q": int(q),
                "k": int(k),
                "ms_412": ms412,
                "ms_417_no_identity_fastpath": ms417_no,
                "ms_417": ms417,
                "speedup_417_vs_412": float(ms412 / ms417),
                "speedup_417_no_identity_vs_412": float(ms412 / ms417_no),
                "rec_412": rec412,
                "rec_417_no_identity_fastpath": rec417_no,
                "rec_417": rec417,
            }
        )
    return summarize_case(n, samples)


def stats(vals: list[float]) -> dict[str, float]:
    arr = np.asarray(vals, dtype=float)
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0,
        "min": float(np.min(arr)),
        "median": float(np.median(arr)),
        "max": float(np.max(arr)),
    }


def summarize_case(n: int, samples: list[dict[str, Any]]) -> dict[str, Any]:
    row = {
        "n": int(n),
        "seeds": len(samples),
        "q": int(samples[0]["q"]),
        "k": int(samples[0]["k"]),
        "ms_412": stats([float(s["ms_412"]) for s in samples]),
        "ms_417_no_identity_fastpath": stats([float(s["ms_417_no_identity_fastpath"]) for s in samples]),
        "ms_417": stats([float(s["ms_417"]) for s in samples]),
        "speedup_417_vs_412": stats([float(s["speedup_417_vs_412"]) for s in samples]),
        "speedup_417_no_identity_vs_412": stats([float(s["speedup_417_no_identity_vs_412"]) for s in samples]),
        "q_active_412": float(np.mean([bool(s["rec_412"]["q_active"]) for s in samples])),
        "q_active_417": float(np.mean([bool(s["rec_417"]["q_active"]) for s in samples])),
        "accepted_412": float(np.mean([bool(s["rec_412"]["accepted"]) for s in samples])),
        "accepted_417": float(np.mean([bool(s["rec_417"]["accepted"]) for s in samples])),
        "identity_fastpath_rate_417": float(np.mean([bool(s["rec_417"]["identity_fastpath"]) for s in samples])),
        "line_evals_412": stats([float(s["rec_412"]["line_search_evals"]) for s in samples]),
        "line_evals_417": stats([float(s["rec_417"]["line_search_evals"]) for s in samples]),
        "target_evals_412": stats([float(s["rec_412"]["target_evals"]) for s in samples]),
        "target_evals_417": stats([float(s["rec_417"]["target_evals"]) for s in samples]),
        "residual_417_max": float(np.max([float(s["rec_417"]["residual"]) for s in samples])),
        "root_error_417_max": float(np.max([float(s["rec_417"]["root_relative_error"]) for s in samples])),
        "samples": samples,
    }
    print(
        f"n={n} q/k={row['q']}/{row['k']} 412={row['ms_412']['mean']:.4f}ms "
        f"417={row['ms_417']['mean']:.4f}ms speed={row['speedup_417_vs_412']['mean']:.2f}x "
        f"q_rate={row['q_active_417']:.2f} err={row['root_error_417_max']:.2e}"
    )
    return row


def main() -> int:
    mod412 = load_engine(ENGINE_412, "pandrosion412_for_417_bench")
    mod417 = load_engine(ENGINE_417, "pandrosion417_for_bench")
    seeds = [41701, 41702, 41703, 41704, 41705]
    dims = [128, 512, 1024, 2048, 4096]
    rows = [run_case(mod412, mod417, n, seeds) for n in dims]
    result = {
        "engine_412": str(ENGINE_412),
        "engine_417": str(ENGINE_417),
        "output": str(OUT),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cpu_count": os.cpu_count(),
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "seeds": seeds,
        "dims": dims,
        "rows": rows,
        "interpretation": {
            "417_no_identity_fastpath": "Same 417 engine with the affine identity target fast path disabled; mostly measures early unit-step line-search exit.",
            "417": "Fast q-cache path: validates F(z)=z-root consistency, projects onto coded q basis directly, then still checks exact target residual.",
            "scope": "This benchmark targets AI adapter validation/vector correction workloads. Polynomial/nonlinear fallback behavior remains inherited from 412.",
        },
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
