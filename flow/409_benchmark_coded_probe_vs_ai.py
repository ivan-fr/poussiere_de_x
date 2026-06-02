#!/usr/bin/env python3
"""Complete benchmark for 409 coded Pandrosion IRP against dense AI matmul."""
from __future__ import annotations

import importlib.util
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import torch


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "flow" / "409_pandrosion_standalone_local_jet_geometry_coded_probe_irp_engine.py"
OUT = Path("/private/tmp/409_pandrosion_coded_probe_complete_benchmark.json")


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
    spec = importlib.util.spec_from_file_location("pandrosion409", ENGINE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {ENGINE}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pandrosion409"] = mod
    spec.loader.exec_module(mod)
    return mod


def configure(
    mod: Any,
    *,
    coded: bool = True,
    sketch_mode: str = "coordinate",
    sketch_solver: str = "svd",
    coded_factor: float = 2.0,
    auto_cap: float = 0.35,
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
            sketch_seed=409,
            sketch_mode=sketch_mode,
            sketch_solver=sketch_solver,
            directional_jet=True,
            directional_jet_min_n=1,
            directional_jet_factor=2.75,
            directional_diff="auto",
            directional_auto_central_cap=auto_cap,
            directional_coded_probe=coded,
            directional_coded_factor=coded_factor,
            directional_coded_min=8,
            directional_coded_max=0,
        )
    )


class IdentityTarget:
    def __init__(self, xstar: Any) -> None:
        self.xstar = np.asarray(xstar, dtype=np.complex128)
        self.evals = 0

    def eval(self, y: Any) -> Any:
        self.evals += 1
        return np.asarray(y, dtype=np.complex128) - self.xstar

    def eval_batch(self, Y: Any) -> Any:
        YY = np.asarray(Y, dtype=np.complex128)
        if YY.ndim == 1:
            return self.eval(YY)[None, :]
        self.evals += int(YY.shape[0])
        return YY - self.xstar[None, :]


class DenseLinearTarget:
    def __init__(self, A: Any, b: Any) -> None:
        self.A = np.asarray(A, dtype=np.complex128)
        self.b = np.asarray(b, dtype=np.complex128)
        self.evals = 0

    def eval(self, y: Any) -> Any:
        self.evals += 1
        return np.asarray(y, dtype=np.complex128) @ self.A.T - self.b

    def eval_batch(self, Y: Any) -> Any:
        YY = np.asarray(Y, dtype=np.complex128)
        if YY.ndim == 1:
            return self.eval(YY)[None, :]
        self.evals += int(YY.shape[0])
        return YY @ self.A.T - self.b[None, :]


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


def make_coded_case(mod: Any, n: int, seed: int, ep: int, rng: Any) -> tuple[Any, Any, int, int]:
    k = int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))
    q = mod.directional_coded_dim_for(n, k)
    salt = seed + 65537 * ep + 0x409D1A
    P = mod._sketch_basis_np(n, k, salt=salt)
    R = mod._coded_probe_basis_np(k, q, salt=salt + 0xC0DEC0DE)
    C = P @ R
    v = (rng.standard_normal(q) + 1j * rng.standard_normal(q)) / math.sqrt(2.0 * q)
    return C @ v, P, k, q


def run_delta(mod: Any, xstar: Any, seed: int, ep: int, target_kind: str = "identity", A: Any | None = None) -> tuple[Any, dict[str, Any], int]:
    if target_kind == "dense":
        if A is None:
            raise ValueError("A is required for dense target")
        target = DenseLinearTarget(A, A @ np.asarray(xstar, dtype=np.complex128))
    else:
        target = IdentityTarget(xstar)
    y0 = np.zeros(len(xstar), dtype=np.complex128)
    f0 = target.eval(y0)
    delta, meta = mod.hypercube_delta(target, y0, f0, ep, seed, 0, 0.0, 0.0, 1e-12, 1e12)
    return delta, meta, int(target.evals)


def coded_identity_bench(mod: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(40901)
    matmul_cache: dict[int, dict[str, Any]] = {}
    for n in [128, 256, 512, 1024, 2048, 4096]:
        seed, ep = 409000 + n, 0
        configure(mod, coded=True, sketch_mode="coordinate", sketch_solver="svd")
        xstar, _P, k, q = make_coded_case(mod, n, seed, ep, rng)

        def run_coded() -> tuple[Any, dict[str, Any], int]:
            configure(mod, coded=True, sketch_mode="coordinate", sketch_solver="svd")
            return run_delta(mod, xstar, seed, ep)

        def run_no_coded() -> tuple[Any, dict[str, Any], int]:
            configure(mod, coded=False, sketch_mode="coordinate", sketch_solver="svd")
            return run_delta(mod, xstar, seed, ep)

        d_coded, m_coded, e_coded = run_coded()
        d_full, m_full, e_full = run_no_coded()
        repeats = 8 if n <= 1024 else (5 if n <= 2048 else 3)
        t_coded = bench(run_coded, repeats=repeats, warmup=2)
        t_full = bench(run_no_coded, repeats=max(3, repeats // 2), warmup=1)
        matmul_cache[n] = matmul_times(n)
        mm = matmul_cache[n]
        row = {
            "n": n,
            "k": k,
            "q": q,
            "case": "coded_identity_result",
            "409_coded_ms": 1e3 * t_coded,
            "409_no_coded_ms": 1e3 * t_full,
            "409_coded_probe_kind": m_coded.get("directional_probe_kind"),
            "409_coded_nodes": int(m_coded.get("hypercube_nodes", 0)),
            "409_coded_evals": int(e_coded),
            "409_coded_rel_error": float(np.linalg.norm(d_coded - xstar) / max(1e-300, np.linalg.norm(xstar))),
            "409_coded_linear_rel": float(m_coded.get("directional_linear_relative_residual", float("nan"))),
            "409_no_coded_nodes": int(m_full.get("hypercube_nodes", 0)),
            "409_no_coded_evals": int(e_full),
            "409_no_coded_rel_error": float(np.linalg.norm(d_full - xstar) / max(1e-300, np.linalg.norm(xstar))),
            "speedup_coded_vs_no_coded": t_full / t_coded if t_coded > 0 else None,
            **mm,
        }
        if "matmul_float16_ms" in mm:
            row["speedup_coded_vs_matmul_fp16"] = mm["matmul_float16_ms"] / row["409_coded_ms"]
            row["beats_matmul_fp16"] = bool(row["speedup_coded_vs_matmul_fp16"] > 1.0)
        if "matmul_float32_ms" in mm:
            row["speedup_coded_vs_matmul_fp32"] = mm["matmul_float32_ms"] / row["409_coded_ms"]
        rows.append(row)
    return rows


def fallback_correctness_bench(mod: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(40902)
    for n in [128, 256, 512, 1024]:
        seed, ep = 409500 + n, 0
        configure(mod, coded=True, sketch_mode="coordinate", sketch_solver="svd")
        k = int(math.ceil(float(mod.DIRECTIONAL_JET_FACTOR) * math.sqrt(float(n))))
        q = mod.directional_coded_dim_for(n, k)
        salt = seed + 65537 * ep + 0x409D1A
        P = mod._sketch_basis_np(n, k, salt=salt)
        u = (rng.standard_normal(k) + 1j * rng.standard_normal(k)) / math.sqrt(2.0 * k)
        xstar = P @ u

        def run() -> tuple[Any, dict[str, Any], int]:
            configure(mod, coded=True, sketch_mode="coordinate", sketch_solver="svd")
            return run_delta(mod, xstar, seed, ep)

        delta, meta, evals = run()
        t = bench(run, repeats=6 if n <= 512 else 4, warmup=2)
        rows.append(
            {
                "n": n,
                "k": k,
                "q": q,
                "case": "full_sketch_result_auto_fallback",
                "ms": 1e3 * t,
                "probe_kind": str(meta.get("directional_probe_kind")),
                "diff": str(meta.get("directional_diff")),
                "nodes": int(meta.get("hypercube_nodes", 0)),
                "evals": int(evals),
                "coded_fallback": meta.get("directional_coded_fallback"),
                "rel_error": float(np.linalg.norm(delta - xstar) / max(1e-300, np.linalg.norm(xstar))),
                "linear_rel": float(meta.get("directional_linear_relative_residual", float("nan"))),
            }
        )
    return rows


def dense_oracle_limit_bench(mod: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(40903)
    for n in [128, 256, 512]:
        seed, ep = 409900 + n, 0
        configure(mod, coded=True, sketch_mode="coordinate", sketch_solver="svd")
        xstar, _P, k, q = make_coded_case(mod, n, seed, ep, rng)
        A = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / math.sqrt(float(n))

        def run() -> tuple[Any, dict[str, Any], int]:
            configure(mod, coded=True, sketch_mode="coordinate", sketch_solver="svd")
            return run_delta(mod, xstar, seed, ep, target_kind="dense", A=A)

        delta, meta, evals = run()
        t = bench(run, repeats=4, warmup=1)
        mm = matmul_times(n)
        row = {
            "n": n,
            "k": k,
            "q": q,
            "case": "coded_dense_linear_oracle",
            "409_ms": 1e3 * t,
            "probe_kind": str(meta.get("directional_probe_kind")),
            "nodes": int(meta.get("hypercube_nodes", 0)),
            "evals": int(evals),
            "rel_error": float(np.linalg.norm(delta - xstar) / max(1e-300, np.linalg.norm(xstar))),
            "linear_rel": float(np.linalg.norm(A @ delta - A @ xstar) / max(1e-300, np.linalg.norm(A @ xstar))),
            **mm,
        }
        if "matmul_float16_ms" in mm:
            row["speedup_vs_matmul_fp16"] = mm["matmul_float16_ms"] / row["409_ms"]
        rows.append(row)
    return rows


def run_engine_case(args: list[str]) -> dict[str, Any]:
    cmd = [sys.executable, str(ENGINE), *args]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True, check=False)
    out_path = args[args.index("--out") + 1] if "--out" in args else None
    summary = None
    params = None
    if out_path and Path(out_path).exists():
        data = json.loads(Path(out_path).read_text(encoding="utf-8"))
        summary = data.get("summary")
        params = data.get("parameters")
    return {
        "cmd": cmd,
        "returncode": proc.returncode,
        "wall_seconds": time.perf_counter() - t0,
        "stdout_tail": proc.stdout.splitlines()[-8:],
        "stderr_tail": proc.stderr.splitlines()[-8:],
        "summary": summary,
        "parameters": params,
    }


def end_to_end_bench() -> list[dict[str, Any]]:
    cases = [
        ["--self-test", "--out", "/private/tmp/409_complete_self.json"],
        ["--system-source", "poly", "--cases", "2,2", "--polys", "x^2 - 1; y^2 - 4", "--variables", "x,y", "--starts", "-2,-3; -2,3; 2,-3; 2,3", "--count", "4", "--pool", "4", "--accept", "1e-8", "--out", "/private/tmp/409_complete_poly2d.json"],
        ["--cases", "2,4", "--count", "2", "--pool", "32", "--accept", "1e-8", "--out", "/private/tmp/409_complete_ks2x4.json"],
    ]
    return [{"args": args, "run": run_engine_case(args)} for args in cases]


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
        "coded_identity": coded_identity_bench(mod),
        "fallback_correctness": fallback_correctness_bench(mod),
        "dense_oracle_limit": dense_oracle_limit_bench(mod),
        "end_to_end": end_to_end_bench(),
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"out={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
