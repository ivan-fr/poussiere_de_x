#!/usr/bin/env python3
"""Multi-seed q-path benchmark for 412 composite-cached Pandrosion IRP."""
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
ENGINE = ROOT / "flow" / "412_pandrosion_standalone_local_jet_geometry_composite_cached_irp_engine.py"
OUT = Path("/private/tmp/412_multiseed_qpath_benchmark.json")


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
    spec = importlib.util.spec_from_file_location("pandrosion412_multiseed", ENGINE)
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
            sketch_seed=412,
            sketch_mode="cached-rademacher",
            sketch_solver="svd",
            sketch_basis_cache=True,
            sketch_basis_cache_max=1024,
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

    def _poly(self, z: Any) -> Any:
        zz = np.asarray(z, dtype=np.complex128)
        return zz + self.beta * zz * zz + self.gamma * zz * np.roll(zz, -1, axis=-1)

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
    salt = int(seed) + 0x412D1A
    C = mod._coded_composite_basis_np(n, k, q, salt=salt)
    rng = np.random.default_rng(seed ^ 0x412A1)
    coeff = scale * (rng.standard_normal(q) + 1j * rng.standard_normal(q)) / math.sqrt(2.0 * q)
    return C @ coeff, k, q


def run_corrector(mod: Any, root: Any, seed: int, *, beta: float, gamma: float) -> dict[str, Any]:
    target = SparseQuadraticTarget(root, beta, gamma)
    y0 = np.zeros_like(root)
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


def matmul_times(n: int) -> dict[str, Any]:
    device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))
    out: dict[str, Any] = {"device": str(device)}
    for dtype in ([torch.float16, torch.float32] if device.type != "cpu" else [torch.float32]):
        a = torch.randn((n, n), dtype=dtype, device=device)
        b = torch.randn((n, n), dtype=dtype, device=device)
        repeats = 8 if n <= 1024 else (5 if n <= 2048 else 3)
        t = bench(lambda: torch.matmul(a, b), repeats=repeats, warmup=3)
        name = str(dtype).replace("torch.", "")
        out[f"matmul_{name}_ms"] = 1e3 * t
        out[f"matmul_{name}_gflops"] = (2.0 * n**3) / t / 1e9
    return out


def stats(vals: list[float]) -> dict[str, float]:
    a = np.asarray(vals, dtype=float)
    return {
        "mean": float(np.mean(a)),
        "std": float(np.std(a, ddof=1)) if a.size > 1 else 0.0,
        "min": float(np.min(a)),
        "median": float(np.median(a)),
        "max": float(np.max(a)),
    }


def run_family(mod: Any, family: str, beta: float, gamma: float, dims: list[int], seeds: list[int]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for n in dims:
        mm = matmul_times(n)
        samples: list[dict[str, Any]] = []
        for base_seed in seeds:
            seed = int(base_seed) + 1009 * n + (17 if gamma else 0)
            root, k, q = make_coded_root(mod, n, seed)

            def run_once() -> dict[str, Any]:
                return run_corrector(mod, root, seed, beta=beta, gamma=gamma)

            rec = run_once()
            repeats = 5 if n <= 1024 else (4 if n <= 2048 else 3)
            ms = 1e3 * bench(run_once, repeats=repeats, warmup=2)
            samples.append({"seed": seed, "ms": ms, "q": q, "k": k, **rec})
        times = [float(s["ms"]) for s in samples]
        residuals = [float(s["residual"]) for s in samples]
        root_errors = [float(s["root_relative_error"]) for s in samples]
        row: dict[str, Any] = {
            "family": family,
            "n": n,
            "seeds": len(samples),
            "q": int(samples[0]["q"]),
            "k": int(samples[0]["k"]),
            "ms": stats(times),
            "q_active_rate": float(np.mean([bool(s["q_active"]) for s in samples])),
            "accepted_rate": float(np.mean([bool(s["accepted"]) for s in samples])),
            "residual": stats(residuals),
            "root_relative_error": stats(root_errors),
            "samples": samples,
            **mm,
        }
        if "matmul_float16_ms" in mm:
            row["speedup_mean_vs_matmul_fp16"] = float(mm["matmul_float16_ms"] / row["ms"]["mean"])
            row["speedup_min_vs_matmul_fp16"] = float(mm["matmul_float16_ms"] / row["ms"]["max"])
        if "matmul_float32_ms" in mm:
            row["speedup_mean_vs_matmul_fp32"] = float(mm["matmul_float32_ms"] / row["ms"]["mean"])
        rows.append(row)
        print(
            f"{family} n={n} q/k={row['q']}/{row['k']} seeds={len(samples)} "
            f"412={row['ms']['mean']:.3f}+/-{row['ms']['std']:.3f}ms "
            f"q_rate={row['q_active_rate']:.2f} fp16_speed={row.get('speedup_mean_vs_matmul_fp16', 0.0):.2f} "
            f"res_max={row['residual']['max']:.2e}"
        )
    return rows


def main() -> int:
    mod = load_engine()
    configure(mod)
    seeds = [41201, 41202, 41203, 41204, 41205, 41206, 41207]
    rows: list[dict[str, Any]] = []
    rows.extend(run_family(mod, "near_linear_quadratic", 1e-8, 0.0, [512, 1024, 2048, 4096], seeds))
    rows.extend(run_family(mod, "near_linear_ring_quadratic", 1e-8, 1e-8, [1024, 2048], seeds))
    result = {
        "engine": str(ENGINE),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cpu_count": os.cpu_count(),
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "seeds": seeds,
        "rows": rows,
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
