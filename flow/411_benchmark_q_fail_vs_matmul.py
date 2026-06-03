#!/usr/bin/env python3
"""Small 411 benchmark where the q-path is forced to fail, compared to matmul."""
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
ENGINE = ROOT / "flow" / "411_pandrosion_standalone_local_jet_geometry_cached_general_irp_engine.py"
OUT = Path("/private/tmp/411_q_fail_vs_matmul_small_benchmark.json")


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
    spec = importlib.util.spec_from_file_location("pandrosion411_qfail", ENGINE)
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
            sketch_seed=411,
            sketch_mode="cached-rademacher",
            sketch_solver="svd",
            sketch_basis_cache=True,
            sketch_basis_cache_max=256,
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
            directional_basis_reuse=True,
        )
    )


class IdentityTarget:
    def __init__(self, root: Any) -> None:
        self.root = np.asarray(root, dtype=np.complex128)
        self.evals = 0

    @property
    def n(self) -> int:
        return int(self.root.size)

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


def make_root_outside_q(mod: Any, n: int, seed: int) -> tuple[Any, int, int]:
    k = max(1, min(n, int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))))
    q = mod.directional_coded_dim_for(n, k)
    salt = int(seed) + 0x411D1A
    P = mod._sketch_basis_np(n, k, salt=salt)
    R = mod._fast_projected_coded_basis_np(k, q, salt=salt + 0xC0DEC0DE)
    rng = np.random.default_rng(seed ^ 0xBAD411)
    u = rng.standard_normal(k) + 1j * rng.standard_normal(k)
    u = u - R @ (R.conj().T @ u)
    u = 0.25 * u / max(1e-300, float(np.linalg.norm(u)))
    return P @ u, k, q


def run_corrector(mod: Any, root: Any, seed: int) -> dict[str, Any]:
    target = IdentityTarget(root)
    y0 = np.zeros_like(root)
    rec = mod.hypercube_matrixjet_corrector(
        target,
        y0,
        max_epochs=3,
        tol=1e-12,
        accept=1e-9,
        trial_timeout=0.0,
        line_search=6,
        line_grid=(1.0, 0.5, 0.25),
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
        "coded_fallback": rec.get("directional_coded_fallback"),
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
        repeats = 8 if n <= 1024 else 5
        ms = bench(lambda: torch.matmul(a, b), repeats=repeats, warmup=3) * 1e3
        out[f"matmul_{str(dtype).replace('torch.', '')}_ms"] = ms
    return out


def main() -> int:
    mod = load_engine()
    configure(mod)
    rows: list[dict[str, Any]] = []
    for n in [512, 1024, 2048]:
        seed = 941100 + n
        root, k, q = make_root_outside_q(mod, n, seed)
        rec = run_corrector(mod, root, seed)
        repeats = 5 if n <= 1024 else 3
        ms = bench(lambda: run_corrector(mod, root, seed), repeats=repeats, warmup=2) * 1e3
        mm = matmul_times(n)
        row = {"n": n, "k": k, "q": q, "411_q_fail_ms": ms, **rec, **mm}
        if "matmul_float16_ms" in mm:
            row["speedup_vs_matmul_fp16"] = mm["matmul_float16_ms"] / ms
            row["beats_matmul_fp16"] = bool(row["speedup_vs_matmul_fp16"] > 1.0)
        if "matmul_float32_ms" in mm:
            row["speedup_vs_matmul_fp32"] = mm["matmul_float32_ms"] / ms
        rows.append(row)
    OUT.write_text(json.dumps({"rows": rows}, indent=2), encoding="utf-8")
    for r in rows:
        print(
            f"n={r['n']} q/k={r['q']}/{r['k']} q_active={r['q_active']} "
            f"probe={r['probe_kind']} nodes={r['nodes']} "
            f"411={r['411_q_fail_ms']:.3f}ms fp16={r.get('matmul_float16_ms', float('nan')):.3f}ms "
            f"speed_fp16={r.get('speedup_vs_matmul_fp16', 0.0):.2f} "
            f"residual={r['residual']:.2e} err={r['root_relative_error']:.2e}"
        )
        print(f"  fallback={r['coded_fallback']} solver={r['solver']}")
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
