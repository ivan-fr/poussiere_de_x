#!/usr/bin/env python3
"""Benchmark 410 on large polynomial systems where the q-coded path is active."""
from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import torch


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "flow" / "410_pandrosion_standalone_local_jet_geometry_fast_projected_coded_irp_engine.py"
OUT = Path("/private/tmp/410_large_polynomial_qpath_benchmark.json")


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
    spec = importlib.util.spec_from_file_location("pandrosion410_poly_qpath", ENGINE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {ENGINE}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pandrosion410_poly_qpath"] = mod
    spec.loader.exec_module(mod)
    return mod


def configure(mod: Any, *, coded: bool = True, fast_projector: bool = True) -> None:
    mod.configure_matrix_backend(
        SimpleNamespace(
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
            sketch_seed=410,
            sketch_mode="coordinate",
            sketch_solver="svd",
            directional_jet=True,
            directional_jet_min_n=1,
            directional_jet_factor=2.75,
            directional_diff="auto",
            directional_auto_central_cap=0.35,
            directional_coded_probe=coded,
            directional_coded_factor=2.0,
            directional_coded_min=8,
            directional_coded_max=0,
            directional_fast_projector=fast_projector,
            directional_fast_projector_cap=0.05,
        )
    )


class SparseQuadraticPolynomialTarget:
    """Vector polynomial: z + beta*z^2 + gamma*z*roll(z,-1) - F(root)."""

    def __init__(self, root: Any, beta: float, gamma: float) -> None:
        self.root = np.asarray(root, dtype=np.complex128)
        self.beta = float(beta)
        self.gamma = float(gamma)
        self.c = self._poly(self.root)
        self.evals = 0

    @property
    def n(self) -> int:
        return int(self.root.size)

    def _poly(self, Z: Any) -> Any:
        ZZ = np.asarray(Z, dtype=np.complex128)
        return ZZ + self.beta * ZZ * ZZ + self.gamma * ZZ * np.roll(ZZ, -1, axis=-1)

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
    k = int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))
    q = mod.directional_coded_dim_for(n, k)
    salt = int(seed) + 0x410D1A
    P = mod._sketch_basis_np(n, k, salt=salt)
    R = mod._fast_projected_coded_basis_np(k, q, salt=salt + 0xC0DEC0DE)
    rng = np.random.default_rng(seed ^ 0x410410)
    coeff = scale * (rng.standard_normal(q) + 1j * rng.standard_normal(q)) / math.sqrt(2.0 * q)
    return (P @ R) @ coeff, k, q


def run_corrector(mod: Any, target: Any, seed: int, *, coded: bool, fast_projector: bool) -> dict[str, Any]:
    configure(mod, coded=coded, fast_projector=fast_projector)
    y0 = np.zeros(target.n, dtype=np.complex128)
    t0 = time.perf_counter()
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
    elapsed = time.perf_counter() - t0
    y = np.asarray(rec.get("y", y0), dtype=np.complex128)
    root_rel = float(np.linalg.norm(y - target.root) / max(1e-300, np.linalg.norm(target.root)))
    return {
        "wall_ms": 1e3 * elapsed,
        "accepted": bool(rec.get("accepted")),
        "status": str(rec.get("status")),
        "epochs": int(rec.get("epochs", 0)),
        "residual": float(rec.get("residual", float("inf"))),
        "root_relative_error": root_rel,
        "target_evals": int(target.evals),
        "hypercube_nodes": int(rec.get("hypercube_nodes", 0) or 0),
        "hypercube_total_evals": int(rec.get("hypercube_total_evals", 0) or 0),
        "hypercube_used_count": int(rec.get("hypercube_used_count", 0) or 0),
        "directional_probe_kind": rec.get("directional_probe_kind"),
        "directional_sketch_dim": rec.get("directional_sketch_dim"),
        "directional_coded_dim": rec.get("directional_coded_dim"),
        "directional_parent_sketch_dim": rec.get("directional_parent_sketch_dim"),
        "directional_oracle_samples": rec.get("directional_oracle_samples"),
        "directional_solver": rec.get("directional_reduced_solver"),
        "directional_fast_projector_accepted": rec.get("directional_reduced_fast_projector_accepted"),
        "directional_fast_projector_fallback": rec.get("directional_reduced_fast_projector_fallback"),
        "directional_linear_relative_residual": rec.get("directional_linear_relative_residual"),
        "q_path_active": bool(rec.get("directional_probe_kind") == "coded-probe" and rec.get("hypercube_nodes") == rec.get("directional_coded_dim")),
    }


def matmul_times(n: int) -> dict[str, Any]:
    device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))
    out: dict[str, Any] = {"device": str(device)}
    for dtype in [torch.float16, torch.float32]:
        if dtype == torch.float16 and device.type == "cpu":
            continue
        a = torch.randn((n, n), dtype=dtype, device=device)
        b = torch.randn((n, n), dtype=dtype, device=device)
        repeats = 8 if n <= 1024 else (5 if n <= 2048 else 3)
        t = bench(lambda: torch.matmul(a, b), repeats=repeats, warmup=3)
        name = str(dtype).replace("torch.", "")
        out[f"matmul_{name}_ms"] = 1e3 * t
        out[f"matmul_{name}_gflops"] = (2.0 * n**3) / t / 1e9
    return out


def polynomial_qpath_bench(mod: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    cases = [
        ("affine_linear_polynomial", 0.0, 0.0, [128, 256, 512, 1024, 2048, 4096]),
        ("near_linear_quadratic", 1e-8, 0.0, [128, 256, 512, 1024, 2048, 4096]),
        ("near_linear_ring_quadratic", 1e-8, 1e-8, [512, 1024, 2048]),
    ]
    for family, beta, gamma, dims in cases:
        for n in dims:
            seed = 410910 + 17 * n + (1 if gamma else 0)
            configure(mod, coded=True, fast_projector=True)
            root, k, q = make_coded_root(mod, n, seed)

            def make_target() -> SparseQuadraticPolynomialTarget:
                return SparseQuadraticPolynomialTarget(root, beta, gamma)

            fast = run_corrector(mod, make_target(), seed, coded=True, fast_projector=True)
            no_fast = run_corrector(mod, make_target(), seed, coded=True, fast_projector=False)
            no_coded = run_corrector(mod, make_target(), seed, coded=False, fast_projector=False)

            repeats = 5 if n <= 1024 else (3 if n <= 2048 else 2)
            fast_ms = 1e3 * bench(lambda: run_corrector(mod, make_target(), seed, coded=True, fast_projector=True), repeats=repeats, warmup=1)
            no_fast_ms = 1e3 * bench(lambda: run_corrector(mod, make_target(), seed, coded=True, fast_projector=False), repeats=max(2, repeats // 2), warmup=1)
            no_coded_ms = 1e3 * bench(lambda: run_corrector(mod, make_target(), seed, coded=False, fast_projector=False), repeats=max(1, repeats // 2), warmup=1)
            mm = matmul_times(n)
            row: dict[str, Any] = {
                "family": family,
                "degree": 2,
                "n": n,
                "k": k,
                "q": q,
                "beta": beta,
                "gamma": gamma,
                "410_fast_ms": fast_ms,
                "410_no_fast_ms": no_fast_ms,
                "410_no_coded_ms": no_coded_ms,
                "speedup_fast_vs_no_fast": no_fast_ms / fast_ms if fast_ms > 0 else None,
                "speedup_fast_vs_no_coded": no_coded_ms / fast_ms if fast_ms > 0 else None,
                "fast": fast,
                "no_fast": no_fast,
                "no_coded": no_coded,
                **mm,
            }
            if "matmul_float16_ms" in mm:
                row["speedup_fast_vs_matmul_fp16"] = mm["matmul_float16_ms"] / fast_ms
                row["beats_matmul_fp16"] = bool(row["speedup_fast_vs_matmul_fp16"] > 1.0)
            if "matmul_float32_ms" in mm:
                row["speedup_fast_vs_matmul_fp32"] = mm["matmul_float32_ms"] / fast_ms
            rows.append(row)
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
        "polynomial_qpath": polynomial_qpath_bench(mod),
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
