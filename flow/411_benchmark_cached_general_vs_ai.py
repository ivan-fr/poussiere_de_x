#!/usr/bin/env python3
"""Complete benchmark for 411 cached Pandrosion IRP against dense AI matmul."""
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
ENGINE_410 = ROOT / "flow" / "410_pandrosion_standalone_local_jet_geometry_fast_projected_coded_irp_engine.py"
ENGINE_411 = ROOT / "flow" / "411_pandrosion_standalone_local_jet_geometry_cached_general_irp_engine.py"
OUT = Path("/private/tmp/411_pandrosion_cached_general_complete_benchmark.json")


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


def load_engine(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def configure_410(mod: Any, *, coded: bool = True, fast_projector: bool = True, sketch_mode: str = "rademacher") -> None:
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
            sketch_mode=sketch_mode,
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


def configure_411(
    mod: Any,
    *,
    coded: bool = True,
    fast_projector: bool = True,
    sketch_mode: str = "cached-rademacher",
    cache: bool = True,
    basis_reuse: bool = True,
) -> None:
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
            sketch_seed=411,
            sketch_mode=sketch_mode,
            sketch_solver="svd",
            sketch_basis_cache=cache,
            sketch_basis_cache_max=256,
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
            directional_basis_reuse=basis_reuse,
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


def matmul_times(n: int) -> dict[str, Any]:
    device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))
    out: dict[str, Any] = {"device": str(device), "n": int(n)}
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


def salt_for(mod: Any, engine: str, seed: int, ep: int = 0) -> int:
    if engine == "411":
        salt_ep = 0 if bool(getattr(mod, "DIRECTIONAL_BASIS_REUSE", True)) else int(ep)
        return int(seed) + 65537 * int(salt_ep) + 0x411D1A
    return int(seed) + 65537 * int(ep) + 0x410D1A


def coded_r_salt(engine: str, salt: int) -> int:
    return int(salt) + 0xC0DEC0DE


def make_coded_root(mod: Any, engine: str, n: int, seed: int, scale: float = 0.25) -> tuple[Any, int, int]:
    k = int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))
    q = mod.directional_coded_dim_for(n, k)
    salt = salt_for(mod, engine, seed, 0)
    P = mod._sketch_basis_np(n, k, salt=salt)
    R = mod._fast_projected_coded_basis_np(k, q, salt=coded_r_salt(engine, salt))
    rng = np.random.default_rng(seed ^ (0x411411 if engine == "411" else 0x410410))
    coeff = scale * (rng.standard_normal(q) + 1j * rng.standard_normal(q)) / math.sqrt(2.0 * q)
    return (P @ R) @ coeff, k, q


def run_corrector(mod: Any, target: Any, seed: int) -> dict[str, Any]:
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
        "basis_cache_hits": rec.get("directional_basis_cache_hits"),
        "basis_cache_misses": rec.get("directional_basis_cache_misses"),
        "basis_cache_entries": rec.get("directional_basis_cache_entries"),
    }


def basis_setup_bench(m410: Any, m411: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for n in [128, 256, 512, 1024, 2048, 4096, 8192]:
        configure_410(m410, sketch_mode="rademacher")
        configure_411(m411, sketch_mode="cached-rademacher", cache=True, basis_reuse=True)
        k = int(math.ceil(2.75 * math.sqrt(float(n))))
        q = m411.directional_coded_dim_for(n, k)
        salt410 = 410000 + n + 0x410D1A
        salt411 = 411000 + n + 0x411D1A
        t410_p = bench(lambda: m410._sketch_basis_np(n, k, salt=salt410), repeats=5, warmup=1)
        t410_r = bench(lambda: m410._fast_projected_coded_basis_np(k, q, salt=salt410 + 0xC0DEC0DE), repeats=5, warmup=1)
        m411._sketch_basis_np(n, k, salt=salt411)
        m411._fast_projected_coded_basis_np(k, q, salt=salt411 + 0xC0DEC0DE)
        t411_p_hot = bench(lambda: m411._sketch_basis_np(n, k, salt=salt411), repeats=8, warmup=1)
        t411_r_hot = bench(lambda: m411._fast_projected_coded_basis_np(k, q, salt=salt411 + 0xC0DEC0DE), repeats=8, warmup=1)
        configure_411(m411, sketch_mode="sparse-sign", cache=False, basis_reuse=True)
        t411_sparse_p = bench(lambda: m411._sketch_basis_np(n, k, salt=salt411 + 17), repeats=5, warmup=1)
        rows.append(
            {
                "n": int(n),
                "k": int(k),
                "q": int(q),
                "410_P_rademacher_qr_ms": 1e3 * t410_p,
                "410_R_qr_ms": 1e3 * t410_r,
                "411_P_cached_hot_ms": 1e3 * t411_p_hot,
                "411_R_cached_hot_ms": 1e3 * t411_r_hot,
                "411_P_sparse_sign_uncached_ms": 1e3 * t411_sparse_p,
                "speedup_411_cached_P_vs_410_P": t410_p / t411_p_hot if t411_p_hot > 0 else None,
                "speedup_411_cached_R_vs_410_R": t410_r / t411_r_hot if t411_r_hot > 0 else None,
                "speedup_411_sparse_P_vs_410_P": t410_p / t411_sparse_p if t411_sparse_p > 0 else None,
            }
        )
    return rows


def polynomial_qpath_bench(m410: Any, m411: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    cases = [
        ("affine_linear_polynomial", 0.0, 0.0, [128, 256, 512, 1024, 2048, 4096]),
        ("near_linear_quadratic", 1e-8, 0.0, [128, 256, 512, 1024, 2048, 4096]),
        ("near_linear_ring_quadratic", 1e-8, 1e-8, [512, 1024, 2048]),
    ]
    matmul_cache: dict[int, dict[str, Any]] = {}
    for family, beta, gamma, dims in cases:
        for n in dims:
            seed = 411910 + 17 * n + (1 if gamma else 0)
            configure_411(m411, sketch_mode="cached-rademacher", cache=True, basis_reuse=True)
            root411, k411, q411 = make_coded_root(m411, "411", n, seed)

            def target411() -> SparseQuadraticPolynomialTarget:
                return SparseQuadraticPolynomialTarget(root411, beta, gamma)

            def run411_cached() -> dict[str, Any]:
                configure_411(m411, sketch_mode="cached-rademacher", cache=True, basis_reuse=True)
                return run_corrector(m411, target411(), seed)

            def run411_nocache() -> dict[str, Any]:
                configure_411(m411, sketch_mode="cached-rademacher", cache=False, basis_reuse=False)
                return run_corrector(m411, target411(), seed)

            def run411_sparse() -> dict[str, Any]:
                configure_411(m411, sketch_mode="sparse-sign", cache=True, basis_reuse=True)
                root_sparse, _, _ = make_coded_root(m411, "411", n, seed + 29)
                return run_corrector(m411, SparseQuadraticPolynomialTarget(root_sparse, beta, gamma), seed + 29)

            configure_410(m410, sketch_mode="rademacher")
            root410, k410, q410 = make_coded_root(m410, "410", n, seed)

            def run410() -> dict[str, Any]:
                configure_410(m410, sketch_mode="rademacher")
                return run_corrector(m410, SparseQuadraticPolynomialTarget(root410, beta, gamma), seed)

            r411 = run411_cached()
            r411_nc = run411_nocache()
            r411_sp = run411_sparse()
            r410 = run410()
            repeats = 6 if n <= 512 else (5 if n <= 1024 else (3 if n <= 2048 else 2))
            t411 = bench(run411_cached, repeats=repeats, warmup=2)
            t411_nc = bench(run411_nocache, repeats=max(2, repeats // 2), warmup=1)
            t411_sp = bench(run411_sparse, repeats=max(2, repeats // 2), warmup=1)
            t410 = bench(run410, repeats=max(2, repeats // 2), warmup=1)
            if n not in matmul_cache:
                matmul_cache[n] = matmul_times(n)
            mm = matmul_cache[n]
            row: dict[str, Any] = {
                "family": family,
                "degree": 2,
                "n": int(n),
                "beta": beta,
                "gamma": gamma,
                "k": int(k411),
                "q": int(q411),
                "410_k": int(k410),
                "410_q": int(q410),
                "411_cached_ms": 1e3 * t411,
                "411_no_cache_no_reuse_ms": 1e3 * t411_nc,
                "411_sparse_sign_ms": 1e3 * t411_sp,
                "410_rademacher_ms": 1e3 * t410,
                "speedup_411_cached_vs_410": t410 / t411 if t411 > 0 else None,
                "speedup_411_cached_vs_411_no_cache": t411_nc / t411 if t411 > 0 else None,
                "speedup_411_sparse_vs_cached": t411 / t411_sp if t411_sp > 0 else None,
                "411_cached": r411,
                "411_no_cache_no_reuse": r411_nc,
                "411_sparse_sign": r411_sp,
                "410_rademacher": r410,
                **mm,
            }
            if "matmul_float16_ms" in mm:
                row["speedup_411_cached_vs_matmul_fp16"] = mm["matmul_float16_ms"] / row["411_cached_ms"]
                row["beats_matmul_fp16"] = bool(row["speedup_411_cached_vs_matmul_fp16"] > 1.0)
            if "matmul_float32_ms" in mm:
                row["speedup_411_cached_vs_matmul_fp32"] = mm["matmul_float32_ms"] / row["411_cached_ms"]
            rows.append(row)
            print(
                f"{family} n={n} q/k={q411}/{k411} "
                f"411={row['411_cached_ms']:.3f}ms 410={row['410_rademacher_ms']:.3f}ms "
                f"fp16={row.get('matmul_float16_ms', float('nan')):.3f}ms "
                f"q_active={r411['q_path_active']}",
                flush=True,
            )
    return rows


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    wins = [r for r in rows if r.get("beats_matmul_fp16")]
    q_active = [r for r in rows if r.get("411_cached", {}).get("q_path_active")]
    return {
        "rows": len(rows),
        "q_path_active_rows": len(q_active),
        "beats_matmul_fp16_rows": len(wins),
        "best_speedup_vs_matmul_fp16": max([float(r.get("speedup_411_cached_vs_matmul_fp16", 0.0)) for r in rows] or [0.0]),
        "best_speedup_vs_410": max([float(r.get("speedup_411_cached_vs_410", 0.0)) for r in rows] or [0.0]),
    }


def main() -> int:
    m410 = load_engine(ENGINE_410, "pandrosion410_for_411_bench")
    m411 = load_engine(ENGINE_411, "pandrosion411_for_bench")
    rows_basis = basis_setup_bench(m410, m411)
    rows_poly = polynomial_qpath_bench(m410, m411)
    result = {
        "engine_411": str(ENGINE_411),
        "engine_410_reference": str(ENGINE_410),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "torch": torch.__version__,
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "torch_threads": int(torch.get_num_threads()),
        "cpu_count": os.cpu_count(),
        "basis_setup": rows_basis,
        "polynomial_qpath": rows_poly,
        "summary": summarize(rows_poly),
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result["summary"], indent=2), flush=True)
    print(f"out={OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
