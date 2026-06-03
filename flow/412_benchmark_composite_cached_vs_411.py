#!/usr/bin/env python3
"""Targeted benchmark for 412 composite-cached q-path against 411 and matmul."""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch


ROOT = Path(__file__).resolve().parents[1]
ENGINE_411 = ROOT / "flow" / "411_pandrosion_standalone_local_jet_geometry_cached_general_irp_engine.py"
ENGINE_412 = ROOT / "flow" / "412_pandrosion_standalone_local_jet_geometry_composite_cached_irp_engine.py"
OUT = Path("/private/tmp/412_composite_cached_vs_411_benchmark.json")


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


def configure(mod: Any, seed: int) -> None:
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
            sketch_seed=seed,
            sketch_mode="cached-rademacher",
            sketch_solver="svd",
            sketch_basis_cache=True,
            sketch_basis_cache_max=512,
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


class NearLinearQuadraticTarget:
    def __init__(self, root: Any) -> None:
        self.root = np.asarray(root, dtype=np.complex128)
        self.c = self._poly(self.root)
        self.evals = 0

    @property
    def n(self) -> int:
        return int(self.root.size)

    def _poly(self, z: Any) -> Any:
        zz = np.asarray(z, dtype=np.complex128)
        return zz + 1e-8 * zz * zz

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


def make_root(mod: Any, engine: str, n: int, seed: int) -> tuple[Any, int, int]:
    k = max(1, min(n, int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))))
    q = mod.directional_coded_dim_for(n, k)
    salt = int(seed) + (0x412D1A if engine == "412" else 0x411D1A)
    if engine == "412":
        C = mod._coded_composite_basis_np(n, k, q, salt=salt)
    else:
        P = mod._sketch_basis_np(n, k, salt=salt)
        R = mod._fast_projected_coded_basis_np(k, q, salt=salt + 0xC0DEC0DE)
        C = P @ R
    rng = np.random.default_rng(seed ^ 0x412412)
    coeff = 0.25 * (rng.standard_normal(q) + 1j * rng.standard_normal(q)) / math.sqrt(2.0 * q)
    return C @ coeff, k, q


def run_corrector(mod: Any, root: Any, seed: int) -> dict[str, Any]:
    target = NearLinearQuadraticTarget(root)
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
        "residual": float(rec.get("residual", float("inf"))),
        "root_relative_error": float(np.linalg.norm(y - root) / max(1e-300, float(np.linalg.norm(root)))),
        "probe_kind": rec.get("directional_probe_kind"),
        "nodes": int(rec.get("hypercube_nodes", 0) or 0),
        "solver": rec.get("directional_reduced_solver"),
        "q_active": bool(rec.get("directional_probe_kind") == "coded-probe" and rec.get("hypercube_nodes") == rec.get("directional_coded_dim")),
    }


def matmul_fp16_ms(n: int) -> float | None:
    device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))
    if device.type == "cpu":
        return None
    a = torch.randn((n, n), dtype=torch.float16, device=device)
    b = torch.randn((n, n), dtype=torch.float16, device=device)
    repeats = 8 if n <= 1024 else (5 if n <= 2048 else 3)
    return 1e3 * bench(lambda: torch.matmul(a, b), repeats=repeats, warmup=3)


def main() -> int:
    m411 = load_engine(ENGINE_411, "pandrosion411_for_412_bench_file")
    m412 = load_engine(ENGINE_412, "pandrosion412_for_bench_file")
    rows: list[dict[str, Any]] = []
    for n in [1024, 2048, 4096]:
        seed = 812000 + n
        configure(m411, 411)
        root411, k411, q411 = make_root(m411, "411", n, seed)
        configure(m412, 412)
        root412, k412, q412 = make_root(m412, "412", n, seed)

        def f411() -> dict[str, Any]:
            configure(m411, 411)
            return run_corrector(m411, root411, seed)

        def f412() -> dict[str, Any]:
            configure(m412, 412)
            return run_corrector(m412, root412, seed)

        r411 = f411()
        r412 = f412()
        repeats = 5 if n <= 2048 else 3
        t411 = 1e3 * bench(f411, repeats=repeats, warmup=2)
        t412 = 1e3 * bench(f412, repeats=repeats, warmup=2)
        mm = matmul_fp16_ms(n)
        row: dict[str, Any] = {
            "n": n,
            "q": q412,
            "k": k412,
            "411_ms": t411,
            "412_ms": t412,
            "speedup_412_vs_411": t411 / t412 if t412 > 0 else None,
            "matmul_fp16_ms": mm,
            "speedup_412_vs_matmul_fp16": (mm / t412 if mm is not None and t412 > 0 else None),
            "411": r411,
            "412": r412,
            "411_q": q411,
            "411_k": k411,
        }
        rows.append(row)
        print(
            f"n={n} q/k={q412}/{k412} 411={t411:.3f}ms 412={t412:.3f}ms "
            f"speedup={row['speedup_412_vs_411']:.2f}x fp16={mm if mm is not None else float('nan'):.3f}ms "
            f"q_active={r412['q_active']} residual={r412['residual']:.2e}"
        )
    OUT.write_text(json.dumps({"rows": rows}, indent=2), encoding="utf-8")
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
