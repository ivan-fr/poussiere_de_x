#!/usr/bin/env python3
"""Complete 417 fast-qcache benchmark against torch.matmul baselines.

The benchmark reports both:

  - square GEMM: A @ B with A,B in R^{n x n}, the usual dense AI matmul proxy;
  - matrix-vector: A @ x with A in R^{n x n}, x in R^n, closer to the vector
    correction shape solved by Pandrosion 417.

Rows are split into q-active roots and random q-rejected roots.  This keeps the
claim bounded: 417 is fast when Pandrosion IRP validates a small q path.
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
ENGINE = ROOT / "flow" / "417_pandrosion_standalone_local_jet_geometry_fast_qcache_irp_engine.py"
OUT = Path("/private/tmp/417_fast_qcache_vs_matmul_benchmark.json")


def sync() -> None:
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    if torch.backends.mps.is_available():
        torch.mps.synchronize()


def bench(fn: Any, repeats: int = 7, warmup: int = 3) -> tuple[float, Any]:
    last: Any = None
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


def load_engine() -> Any:
    spec = importlib.util.spec_from_file_location("pandrosion417_vs_matmul", ENGINE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {ENGINE}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def configure(mod: Any, *, fastpath: bool = True) -> None:
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
            sketch_seed=417,
            sketch_mode="cached-rademacher",
            sketch_solver="svd",
            sketch_basis_cache=True,
            sketch_basis_cache_max=8192,
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


def device() -> Any:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def stats(vals: list[float]) -> dict[str, float]:
    a = np.asarray(vals, dtype=float)
    return {
        "mean": float(np.mean(a)),
        "std": float(np.std(a, ddof=1)) if a.size > 1 else 0.0,
        "min": float(np.min(a)),
        "median": float(np.median(a)),
        "max": float(np.max(a)),
    }


def q_dims(mod: Any, n: int) -> tuple[int, int]:
    k = max(1, min(int(n), int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))))
    q = int(mod.directional_coded_dim_for(int(n), int(k)))
    return k, q


def subspace_energy(mod: Any, root: Any, seed: int) -> dict[str, Any]:
    r = np.asarray(root, dtype=np.complex128).reshape(-1)
    n = int(r.size)
    k, q = q_dims(mod, n)
    salt = int(seed) + 0x412D1A
    C = np.asarray(mod._coded_composite_basis_np(n, k, q, salt=salt), dtype=np.complex128)
    P = np.asarray(mod._sketch_basis_np(n, k, salt=salt), dtype=np.complex128)
    den = max(1e-300, float(np.vdot(r, r).real))
    q_coeff = C.conj().T @ r
    p_coeff = P.conj().T @ r
    return {
        "q_energy_fraction": float(np.vdot(q_coeff, q_coeff).real / den),
        "full_sketch_energy_fraction": float(np.vdot(p_coeff, p_coeff).real / den),
        "expected_random_q_fraction": float(q / max(1, n)),
        "expected_random_full_sketch_fraction": float(k / max(1, n)),
    }


def make_q_root(mod: Any, n: int, seed: int, scale: float = 0.25) -> tuple[Any, int, int]:
    k, q = q_dims(mod, n)
    C = mod._coded_composite_basis_np(int(n), k, q, salt=int(seed) + 0x412D1A)
    rng = np.random.default_rng(int(seed) ^ 0x417A11)
    coeff = rng.standard_normal(q) + 1j * rng.standard_normal(q)
    coeff = float(scale) * coeff / max(1e-300, float(np.linalg.norm(coeff)))
    return C @ coeff.astype(np.complex128, copy=False), k, q


def make_random_root(n: int, seed: int, scale: float = 0.25) -> Any:
    rng = np.random.default_rng(int(seed) ^ 0x417BAD)
    v = rng.standard_normal(int(n)) + 1j * rng.standard_normal(int(n))
    return float(scale) * v.astype(np.complex128, copy=False) / max(1e-300, float(np.linalg.norm(v)))


def run_417(mod: Any, root: Any, seed: int, *, repeats: int) -> tuple[float, dict[str, Any]]:
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
        identity_fast = bool(rec.get("directional_reduced_identity_target_fastpath", False))
        q_active = bool(
            rec.get("directional_probe_kind") == "coded-probe"
            and (rec.get("hypercube_nodes") == rec.get("directional_coded_dim") or identity_fast)
        )
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
            "identity_fastpath": identity_fast,
            "line_search_evals": int(rec.get("line_search_evals", 0) or 0),
            "early_unit_step": bool(rec.get("line_search_early_unit_step", False)),
            "target_evals": int(target.evals),
            "q_active": q_active,
            **subspace_energy(mod, root_np, seed),
        }

    ms, rec = bench(once, repeats=max(1, repeats), warmup=3)
    return ms, rec


def matmul_baselines(n: int) -> dict[str, Any]:
    dev = device()
    out: dict[str, Any] = {"device": str(dev)}
    dtypes = [torch.float16, torch.float32] if dev.type != "cpu" else [torch.float32]
    for dtype in dtypes:
        name = str(dtype).replace("torch.", "")
        repeats = 14 if n <= 512 else (10 if n <= 2048 else 5)
        warmup = 5 if n <= 2048 else 3
        a = torch.randn((int(n), int(n)), dtype=dtype, device=dev)
        b = torch.randn((int(n), int(n)), dtype=dtype, device=dev)
        x = torch.randn((int(n),), dtype=dtype, device=dev)
        ms_gemm, _ = bench(lambda: torch.matmul(a, b), repeats=repeats, warmup=warmup)
        ms_mv, _ = bench(lambda: torch.matmul(a, x), repeats=max(repeats, 20), warmup=warmup)
        out[f"square_matmul_{name}_ms"] = float(ms_gemm)
        out[f"square_matmul_{name}_gflops"] = float((2.0 * int(n) ** 3) / (ms_gemm / 1e3) / 1e9)
        out[f"matrix_vector_{name}_ms"] = float(ms_mv)
        out[f"matrix_vector_{name}_gflops"] = float((2.0 * int(n) ** 2) / (ms_mv / 1e3) / 1e9)
    return out


def run_family(mod: Any, family: str, dims: list[int], seeds: list[int], matmul_cache: dict[int, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for n in dims:
        if n not in matmul_cache:
            matmul_cache[n] = matmul_baselines(n)
        samples: list[dict[str, Any]] = []
        for base_seed in seeds:
            seed = int(base_seed) + 1009 * int(n)
            if family == "q_active":
                root, k, q = make_q_root(mod, n, seed)
                repeats = 11 if n <= 1024 else 7
            elif family == "q_rejected_random":
                k, q = q_dims(mod, n)
                root = make_random_root(n, seed)
                repeats = 5 if n <= 1024 else 3
            else:
                raise ValueError(f"unknown family {family!r}")
            configure(mod, fastpath=True)
            ms, rec = run_417(mod, root, seed, repeats=repeats)
            sample = {"seed": seed, "n": n, "q": q, "k": k, "ms_417": float(ms), **rec}
            for key, val in matmul_cache[n].items():
                sample[key] = val
            if "square_matmul_float16_ms" in sample:
                sample["speedup_vs_square_matmul_fp16"] = float(sample["square_matmul_float16_ms"] / max(1e-300, ms))
            if "square_matmul_float32_ms" in sample:
                sample["speedup_vs_square_matmul_fp32"] = float(sample["square_matmul_float32_ms"] / max(1e-300, ms))
            if "matrix_vector_float16_ms" in sample:
                sample["speedup_vs_matrix_vector_fp16"] = float(sample["matrix_vector_float16_ms"] / max(1e-300, ms))
            if "matrix_vector_float32_ms" in sample:
                sample["speedup_vs_matrix_vector_fp32"] = float(sample["matrix_vector_float32_ms"] / max(1e-300, ms))
            samples.append(sample)

        row = summarize_family(family, n, samples)
        rows.append(row)
        print(
            f"{family} n={n} q/k={row['q']}/{row['k']} 417={row['ms_417']['mean']:.4f}ms "
            f"q_rate={row['q_active_rate']:.2f} fp16_gemm_speed={row.get('speedup_vs_square_matmul_fp16', {}).get('mean', 0.0):.2f} "
            f"fp16_mv_speed={row.get('speedup_vs_matrix_vector_fp16', {}).get('mean', 0.0):.3f} "
            f"err_max={row['root_relative_error']['max']:.2e}"
        )
    return rows


def summarize_family(family: str, n: int, samples: list[dict[str, Any]]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "family": family,
        "n": int(n),
        "seeds": len(samples),
        "q": int(samples[0]["q"]),
        "k": int(samples[0]["k"]),
        "ms_417": stats([float(s["ms_417"]) for s in samples]),
        "q_active_rate": float(np.mean([bool(s["q_active"]) for s in samples])),
        "accepted_rate": float(np.mean([bool(s["accepted"]) for s in samples])),
        "identity_fastpath_rate": float(np.mean([bool(s["identity_fastpath"]) for s in samples])),
        "residual": stats([float(s["residual"]) for s in samples]),
        "root_relative_error": stats([float(s["root_relative_error"]) for s in samples]),
        "q_energy_fraction": stats([float(s["q_energy_fraction"]) for s in samples]),
        "target_evals": stats([float(s["target_evals"]) for s in samples]),
        "line_search_evals": stats([float(s["line_search_evals"]) for s in samples]),
        "samples": samples,
    }
    baseline_keys = [
        "square_matmul_float16_ms",
        "square_matmul_float32_ms",
        "matrix_vector_float16_ms",
        "matrix_vector_float32_ms",
        "square_matmul_float16_gflops",
        "square_matmul_float32_gflops",
        "matrix_vector_float16_gflops",
        "matrix_vector_float32_gflops",
        "speedup_vs_square_matmul_fp16",
        "speedup_vs_square_matmul_fp32",
        "speedup_vs_matrix_vector_fp16",
        "speedup_vs_matrix_vector_fp32",
    ]
    for key in baseline_keys:
        vals = [float(s[key]) for s in samples if key in s and math.isfinite(float(s[key]))]
        if vals:
            row[key] = stats(vals)
    return row


def main() -> int:
    mod = load_engine()
    configure(mod, fastpath=True)
    seeds = [41701, 41702, 41703, 41704, 41705]
    dims = [128, 512, 1024, 2048, 4096]
    matmul_cache: dict[int, dict[str, Any]] = {}
    rows: list[dict[str, Any]] = []
    rows.extend(run_family(mod, "q_active", dims, seeds, matmul_cache))
    rows.extend(run_family(mod, "q_rejected_random", [128, 512, 1024, 2048], seeds, matmul_cache))
    result = {
        "engine": str(ENGINE),
        "output": str(OUT),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cpu_count": os.cpu_count(),
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "device": str(device()),
        "seeds": seeds,
        "dims_q_active": dims,
        "dims_q_rejected_random": [128, 512, 1024, 2048],
        "rows": rows,
        "interpretation": {
            "square_matmul": "Dense n-by-n GEMM proxy often used to characterize AI matrix throughput.",
            "matrix_vector": "Dense n-by-n by vector baseline, closer to a single vector correction shape.",
            "q_active": "Constructed roots live in the Pandrosion coded q basis and validate the fast path.",
            "q_rejected_random": "Random roots are not q-aligned; this is the negative control where 417 should not claim a q-path win.",
            "claim_scope": "417 can beat dense square matmul when q is active because it avoids dense GEMM. It is not a universal replacement for dense matrix multiplication.",
        },
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
