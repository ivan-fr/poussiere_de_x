#!/usr/bin/env python3
"""
411_pandrosion_standalone_local_jet_geometry_cached_general_irp_engine.py

Standalone PyTorch/NumPy engine that reuses the 317 local-jet geometry and keeps
Pandrosion IRP, with a batched GEMM prepass before each IRP wave.

Mathematical idea:
  1. Keep 317's faithful geometric oracle: around the current point z, sample the
     true system F on a small hypercube and fit a local jet

         F(z + delta) ~= F(z) + J delta + higher local defects.

  2. Do not call a multivariate algebraic solver, Groebner/resultant machinery, or a
     fragile square solve as the main numerical primitive.  Each correction is the
     solution of a matrix least-squares problem on the sampled jet:

         min_delta ||J delta + F(z)||_2^2 + mu ||delta||_2^2.

     The 410 default is adaptive-directional-sketch with a coded one-sided
     directional jet: for small systems it keeps the fastest stable
     GEMM/Newton-Schulz path; for larger local jets it first samples q coded
     directions, where q ~= O(log n), before falling back to the safer
     k ~= O(sqrt(n)) directional basis.  The goal is not to beat a dense matmul
     kernel on its own terms, but to avoid constructing and multiplying a dense
     n-by-n local jet when Pandrosion IRP exposes a smaller direction space.

  3. 410 prepares waves of starts with a fused batched direct prepass.  The prepass
     can use the complex 404-style kernel, the 405 real-packed kernel, or the new
     sketch kernel, then hands every candidate to the normal Pandrosion lazy IRP
     chart corrector.

This is still a numerical root extractor for F=0, but the engine surface is local
matrix calculus on faithful sampled jets rather than direct multivariate system
resolution.  The 410 difference from 406 is that the local jet itself can be
directional:

    delta = P u,       min_u ||J P u + F||^2 + mu ||u||^2,

where J P is estimated directly from black-box oracle samples.  407 used central
differences

    F(y +/- h P[:,j])

while 408 first tried the one-sided SPSA-like/secant form

    F(y + h P[:,j]) - F(y)

410 adds an even smaller coded probe basis D = P R, q << k, and solves

    delta = D v,       min_v ||J D v + F||^2 + mu ||v||^2.

410 adds a fast projected reduced solver before SVD/Newton-Schulz:

    v ~= diag((JD)*JD + mu I)^-1 (JD)*(-F).

This is exact enough when the coded probe is close to an isometry, and it removes
the CPU SVD overhead seen on MPS.  If the projected step has a poor reduced
residual, 410 falls back to the 409/408 reduced solver.  This can beat dense
matmul only by changing the problem size and avoiding the dense jet.  The full
hypercube, complex GEMM, real-packed GEMM, and SVD paths remain available as
diagnostics/fallbacks.  The full matrix kernel remains a GEMM-first pseudo-inverse:

411 keeps the same Pandrosion IRP surface and attacks the fixed Python/NumPy
overhead that dominated the smaller benchmarks:

    - deterministic cache for the sketch basis P and coded basis R;
    - optional reuse of the same directional basis across corrector epochs;
    - sparse-sign sketch mode, which avoids dense QR construction when speed is
      more important than an orthonormal randomized basis.

This still does not claim a universal dense-GEMM replacement.  General matmul is
already close to hardware peak through BLAS/CUTLASS-style tiling, packing, and
micro-kernels.  411 wins when Pandrosion IRP can reduce the effective dimension
or when repeated local jets amortize the cached basis construction.

    (J*J + mu I) delta = -J*F

is applied through a Newton-Schulz inverse approximation using mostly matrix
multiplications (`torch.matmul` / `torch.bmm`).  This is deliberately closer to the
matrix kernels used by AI workloads than SVD.  The SVD path remains available as a
robust fallback/diagnostic.

Systems supported as the geometric oracle:
  - exact user polynomial systems (--system-source poly --polys "...");
  - KS/Kostlan systems (dense / lazy-feature / geometry-kernel backends).

Dependencies: Python stdlib + NumPy + PyTorch only.  No local flow imports.
"""
from __future__ import annotations

import argparse
import ast
import dataclasses
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
import torch


MATRIX_BACKEND = "torch"
MATRIX_ALGORITHM = "adaptive-directional-sketch"
BATCH_KERNEL = "auto"
TORCH_DEVICE = "auto"
TORCH_COMPLEX_DTYPE = "auto"
TORCH_REAL_DTYPE = "auto"
TORCH_RESOLVED_DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"
TORCH_RESOLVED_DTYPE = "complex64"
TORCH_RESOLVED_REAL_DTYPE = "float32"
NS_ITERS = 8
NS_DAMPING_FLOOR = 1e-5
SKETCH_DIM = 0
SKETCH_FACTOR = 2.75
SKETCH_MIN_N = 96
SKETCH_SEED = 411
SKETCH_MODE = "cached-rademacher"
SKETCH_SOLVER = "svd"
SKETCH_BASIS_CACHE = True
SKETCH_BASIS_CACHE_MAX = 128
DIRECTIONAL_JET = True
DIRECTIONAL_JET_MIN_N = 96
DIRECTIONAL_JET_FACTOR = 2.75
DIRECTIONAL_DIFF_MODE = "auto"
DIRECTIONAL_AUTO_CENTRAL_CAP = 0.35
DIRECTIONAL_CODED_PROBE = True
DIRECTIONAL_CODED_FACTOR = 2.0
DIRECTIONAL_CODED_MIN = 8
DIRECTIONAL_CODED_MAX = 0
DIRECTIONAL_FAST_PROJECTOR = True
DIRECTIONAL_FAST_PROJECTOR_CAP = 0.05
DIRECTIONAL_BASIS_REUSE = True
_BASIS_CACHE: dict[tuple[Any, ...], Any] = {}
_BASIS_CACHE_HITS = 0
_BASIS_CACHE_MISSES = 0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def now() -> float:
    return time.time()


def cjson(z: complex) -> list[float]:
    w = complex(z)
    return [float(w.real), float(w.imag)]


def root_to_json(z: Sequence[complex]) -> list[list[float]]:
    return [cjson(v) for v in z]


def parse_case(raw: str) -> tuple[int, int]:
    s = str(raw).strip().lower().replace("x", ",").replace(":", ",")
    p = [q.strip() for q in s.split(",") if q.strip()]
    if len(p) != 2:
        raise ValueError(f"case must be n,d, got {raw!r}")
    n, d = int(p[0]), int(p[1])
    if n <= 0 or d <= 0:
        raise ValueError("n,d must be positive")
    return n, d


def parse_float_list(raw: Optional[str], default: Sequence[float], positive: bool = False) -> list[float]:
    if raw is None or str(raw).strip() == "":
        return list(default)
    out: list[float] = []
    for part in str(raw).replace(";", ",").split(","):
        try:
            x = float(part.strip())
        except Exception:
            continue
        if math.isfinite(x) and ((not positive) or x > 0):
            out.append(x)
    return out or list(default)


def splitmix64(x: int) -> int:
    x = (int(x) + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return (x ^ (x >> 31)) & 0xFFFFFFFFFFFFFFFF


def u01(x: int) -> float:
    return ((splitmix64(x) >> 11) & ((1 << 53) - 1)) / float(1 << 53)


def phase(theta: float) -> complex:
    return complex(math.cos(theta), math.sin(theta))


def stable_seed(n: int, d: int, seed_index: int = 0, salt: int = 0) -> int:
    return int(splitmix64(0x50414E44524F5349 + 1000003 * n + 9176 * d + 97 * seed_index + salt) & 0x7FFFFFFF)


def exact_kostlan_terms(n: int, d: int) -> int:
    return int(math.comb(int(n) + int(d), int(d)))


def auto_lazy_feature_count(n: int, d: int, cap: int) -> int:
    base = int(32 * math.sqrt(max(1, (int(n) + 1) * (int(d) + 1))))
    return int(max(int(n) + 8, min(max(1, int(cap)), max(512, base))))


def auto_geometry_anchor_count(n: int, d: int, cap: int) -> int:
    base = int(24 * math.sqrt(max(1, (int(n) + 1) * (int(d) + 1))))
    return int(max(int(n) + 16, min(max(1, int(cap)), max(512, base))))


def finite_norm(v: Any) -> float:
    try:
        r = float(np.linalg.norm(v))
        return r if math.isfinite(r) else float("inf")
    except Exception:
        return float("inf")


def log_energy(y: Sequence[complex]) -> float:
    yy = np.asarray(y, dtype=np.complex128)
    if yy.size == 0:
        return 0.0
    return float(np.mean(np.abs(np.log(np.maximum(np.abs(yy), 1e-300)))))


def realness(z: Sequence[complex]) -> float:
    zz = np.asarray(z, dtype=np.complex128)
    den = max(1e-300, float(np.linalg.norm(zz)))
    return float(np.linalg.norm(zz.imag) / den)


def score_root(residual: float, realness_value: float, cond: Optional[float]) -> float:
    c = float(cond) if cond is not None and math.isfinite(float(cond)) else 1e300
    return float(math.log1p(max(0.0, residual)) + 0.01 * math.log1p(max(0.0, realness_value)) + 0.001 * math.log1p(max(0.0, c)))


def cluster_index(roots: list[dict[str, Any]], z: Any, sep: float) -> Optional[int]:
    zz = np.asarray(z, dtype=np.complex128)
    for i, r in enumerate(roots):
        if float(np.linalg.norm(zz - np.asarray(r["z_complex"], dtype=np.complex128))) <= float(sep):
            return i
    return None


def parse_complex_token(raw: str) -> complex:
    s = str(raw).strip().replace("i", "j")
    if not s:
        raise ValueError("empty complex token")
    return complex(s)


def parse_start_points(raw: Optional[str], n: int) -> list[Any]:
    if raw is None or str(raw).strip() == "":
        return []
    text = str(raw).strip()
    if int(n) == 1 and ";" not in text and "|" not in text:
        return [np.asarray([parse_complex_token(p)], dtype=np.complex128) for p in text.split(",") if p.strip()]
    pts = []
    for part in text.replace("|", ";").split(";"):
        if not part.strip():
            continue
        coords = [p.strip() for p in part.split(",") if p.strip()]
        if len(coords) != int(n):
            raise ValueError(f"start point {part!r} has {len(coords)} coord(s), expected {n}")
        pts.append(np.asarray([parse_complex_token(p) for p in coords], dtype=np.complex128))
    return pts


def safe_norms(F: Any) -> Any:
    R = np.linalg.norm(np.asarray(F, dtype=np.complex128), axis=1)
    R = np.asarray(R, dtype=float)
    R[~np.isfinite(R)] = np.inf
    return R


def resolve_torch_device(raw: str) -> str:
    mode = str(raw or "auto").strip().lower()
    if mode == "auto":
        if torch.cuda.is_available():
            return "cuda"
        if torch.backends.mps.is_available():
            return "mps"
        return "cpu"
    if mode == "cuda" and not torch.cuda.is_available():
        return "cpu"
    if mode == "mps" and not torch.backends.mps.is_available():
        return "cpu"
    return mode if mode in {"cpu", "cuda", "mps"} else "cpu"


def resolve_torch_complex_dtype(raw: str, device: str) -> str:
    mode = str(raw or "auto").strip().lower()
    if mode == "auto":
        return "complex64" if str(device) == "mps" else "complex128"
    if mode in {"complex64", "c64", "float32"}:
        return "complex64"
    if mode in {"complex128", "c128", "float64"}:
        return "complex64" if str(device) == "mps" else "complex128"
    return "complex64" if str(device) == "mps" else "complex128"


def resolve_torch_real_dtype(raw: str, device: str) -> str:
    mode = str(raw or "auto").strip().lower()
    if mode == "auto":
        return "float32"
    if mode in {"float16", "fp16", "half"}:
        return "float16"
    if mode in {"bfloat16", "bf16"}:
        return "bfloat16"
    if mode in {"float64", "fp64", "double"}:
        return "float32" if str(device) == "mps" else "float64"
    return "float32"


def configure_matrix_backend(args: argparse.Namespace) -> None:
    global MATRIX_BACKEND, MATRIX_ALGORITHM, BATCH_KERNEL, TORCH_DEVICE, TORCH_COMPLEX_DTYPE, TORCH_REAL_DTYPE, TORCH_RESOLVED_DEVICE, TORCH_RESOLVED_DTYPE, TORCH_RESOLVED_REAL_DTYPE, NS_ITERS, NS_DAMPING_FLOOR, SKETCH_DIM, SKETCH_FACTOR, SKETCH_MIN_N, SKETCH_SEED, SKETCH_MODE, SKETCH_SOLVER, SKETCH_BASIS_CACHE, SKETCH_BASIS_CACHE_MAX, DIRECTIONAL_JET, DIRECTIONAL_JET_MIN_N, DIRECTIONAL_JET_FACTOR, DIRECTIONAL_DIFF_MODE, DIRECTIONAL_AUTO_CENTRAL_CAP, DIRECTIONAL_CODED_PROBE, DIRECTIONAL_CODED_FACTOR, DIRECTIONAL_CODED_MIN, DIRECTIONAL_CODED_MAX, DIRECTIONAL_FAST_PROJECTOR, DIRECTIONAL_FAST_PROJECTOR_CAP, DIRECTIONAL_BASIS_REUSE
    backend = str(getattr(args, "matrix_backend", "torch") or "torch").strip().lower()
    MATRIX_BACKEND = backend if backend in {"torch", "numpy", "auto"} else "torch"
    algorithm = str(getattr(args, "matrix_algorithm", "adaptive-directional-sketch") or "adaptive-directional-sketch").strip().lower().replace("_", "-")
    MATRIX_ALGORITHM = algorithm if algorithm in {"adaptive-directional-sketch", "directional-sketch", "adaptive-sketch", "sketch-ns", "adaptive-ns", "realpack-ns", "gemm-ns", "svd", "auto"} else "adaptive-directional-sketch"
    batch_kernel = str(getattr(args, "batch_kernel", "auto") or "auto").strip().lower().replace("_", "-")
    BATCH_KERNEL = batch_kernel if batch_kernel in {"auto", "sketch", "realpack", "complex"} else "auto"
    TORCH_DEVICE = str(getattr(args, "torch_device", "auto") or "auto").strip().lower()
    TORCH_COMPLEX_DTYPE = str(getattr(args, "torch_complex_dtype", "auto") or "auto").strip().lower()
    TORCH_REAL_DTYPE = str(getattr(args, "torch_real_dtype", "auto") or "auto").strip().lower()
    TORCH_RESOLVED_DEVICE = resolve_torch_device(TORCH_DEVICE)
    TORCH_RESOLVED_DTYPE = resolve_torch_complex_dtype(TORCH_COMPLEX_DTYPE, TORCH_RESOLVED_DEVICE)
    TORCH_RESOLVED_REAL_DTYPE = resolve_torch_real_dtype(TORCH_REAL_DTYPE, TORCH_RESOLVED_DEVICE)
    NS_ITERS = max(1, int(getattr(args, "ns_iters", 8)))
    NS_DAMPING_FLOOR = max(0.0, float(getattr(args, "ns_damping_floor", 1e-5)))
    SKETCH_DIM = max(0, int(getattr(args, "sketch_dim", 0)))
    SKETCH_FACTOR = max(1.0, float(getattr(args, "sketch_factor", 2.75)))
    SKETCH_MIN_N = max(1, int(getattr(args, "sketch_min_n", 96)))
    SKETCH_SEED = int(getattr(args, "sketch_seed", 411))
    mode = str(getattr(args, "sketch_mode", "cached-rademacher") or "cached-rademacher").strip().lower().replace("_", "-")
    SKETCH_MODE = mode if mode in {"rademacher", "cached-rademacher", "coordinate", "sparse-sign", "auto"} else "rademacher"
    solver = str(getattr(args, "sketch_solver", "svd") or "svd").strip().lower().replace("_", "-")
    SKETCH_SOLVER = solver if solver in {"svd", "ns", "auto"} else "svd"
    SKETCH_BASIS_CACHE = bool(getattr(args, "sketch_basis_cache", True))
    SKETCH_BASIS_CACHE_MAX = max(0, int(getattr(args, "sketch_basis_cache_max", 128)))
    DIRECTIONAL_JET = bool(getattr(args, "directional_jet", True))
    DIRECTIONAL_JET_MIN_N = max(1, int(getattr(args, "directional_jet_min_n", SKETCH_MIN_N)))
    DIRECTIONAL_JET_FACTOR = max(1.0, float(getattr(args, "directional_jet_factor", SKETCH_FACTOR)))
    diff_mode = str(getattr(args, "directional_diff", "auto") or "auto").strip().lower().replace("_", "-")
    DIRECTIONAL_DIFF_MODE = diff_mode if diff_mode in {"auto", "forward", "central"} else "auto"
    DIRECTIONAL_AUTO_CENTRAL_CAP = max(0.0, float(getattr(args, "directional_auto_central_cap", 0.35)))
    DIRECTIONAL_CODED_PROBE = bool(getattr(args, "directional_coded_probe", True))
    DIRECTIONAL_CODED_FACTOR = max(0.25, float(getattr(args, "directional_coded_factor", 2.0)))
    DIRECTIONAL_CODED_MIN = max(1, int(getattr(args, "directional_coded_min", 8)))
    DIRECTIONAL_CODED_MAX = max(0, int(getattr(args, "directional_coded_max", 0)))
    DIRECTIONAL_FAST_PROJECTOR = bool(getattr(args, "directional_fast_projector", True))
    DIRECTIONAL_FAST_PROJECTOR_CAP = max(0.0, float(getattr(args, "directional_fast_projector_cap", 0.05)))
    DIRECTIONAL_BASIS_REUSE = bool(getattr(args, "directional_basis_reuse", True))


def torch_complex_dtype() -> Any:
    return torch.complex64 if TORCH_RESOLVED_DTYPE == "complex64" else torch.complex128


def torch_real_dtype() -> Any:
    if TORCH_RESOLVED_REAL_DTYPE == "float16":
        return torch.float16
    if TORCH_RESOLVED_REAL_DTYPE == "bfloat16":
        return torch.bfloat16
    if TORCH_RESOLVED_REAL_DTYPE == "float64":
        return torch.float64
    return torch.float32


def torch_device() -> Any:
    return torch.device(TORCH_RESOLVED_DEVICE)


def selected_matrix_algorithm_410(dim: int = 0) -> str:
    if MATRIX_ALGORITHM == "directional-sketch":
        return "sketch-ns"
    if MATRIX_ALGORITHM not in {"adaptive-directional-sketch", "adaptive-sketch", "adaptive-ns"}:
        return MATRIX_ALGORITHM
    if MATRIX_ALGORITHM in {"adaptive-directional-sketch", "adaptive-sketch"} and int(dim) >= int(SKETCH_MIN_N):
        return "sketch-ns"
    # On current MPS, complex64 GEMM beats the doubled real block for the small
    # local jets used by Pandrosion. CUDA can profit from the real GEMM path on
    # larger blocks, while retaining the same IRP/oracle surface.
    if TORCH_RESOLVED_DEVICE == "cuda" and int(dim) >= 128:
        return "realpack-ns"
    return "gemm-ns"


def _svd_threshold(s: Any, rcond: float, condition_cap: float) -> tuple[float, float, float]:
    ss = np.asarray(s, dtype=float)
    if ss.size == 0 or not np.any(np.isfinite(ss)):
        return 0.0, 0.0, 0.0
    smax = float(np.nanmax(ss))
    if not math.isfinite(smax) or smax <= 0.0:
        return 0.0, 0.0, 0.0
    rel = max(0.0, float(rcond)) * smax
    if math.isfinite(float(condition_cap)) and float(condition_cap) > 0.0:
        rel = max(rel, smax / float(condition_cap))
    return float(rel), float(smax), float(np.nanmin(ss[np.isfinite(ss)]))


def _torch_ns_inverse(B: Any, iters: int) -> Any:
    n = int(B.shape[-1])
    I = torch.eye(n, dtype=B.dtype, device=B.device)
    norm1 = torch.max(torch.sum(torch.abs(B), dim=-2))
    norminf = torch.max(torch.sum(torch.abs(B), dim=-1))
    alpha = torch.clamp(norm1 * norminf, min=torch.finfo(torch.float32).tiny)
    X = B.conj().transpose(-2, -1) / alpha
    for _ in range(max(1, int(iters))):
        X = X @ (2.0 * I - B @ X)
    return X


def sketch_dim_for(n: int) -> int:
    nn = max(1, int(n))
    if int(SKETCH_DIM) > 0:
        return max(1, min(nn, int(SKETCH_DIM)))
    if nn < int(SKETCH_MIN_N):
        return nn
    return max(1, min(nn, int(math.ceil(float(SKETCH_FACTOR) * math.sqrt(float(nn))))))


def directional_coded_dim_for(n: int, k: int) -> int:
    kk = max(1, int(k))
    q = int(math.ceil(float(DIRECTIONAL_CODED_FACTOR) * math.log2(float(max(3, int(n) + 1)))))
    q = max(int(DIRECTIONAL_CODED_MIN), q)
    if int(DIRECTIONAL_CODED_MAX) > 0:
        q = min(q, int(DIRECTIONAL_CODED_MAX))
    return max(1, min(kk, q))


def _basis_cache_lookup(key: tuple[Any, ...]) -> Optional[Any]:
    global _BASIS_CACHE_HITS, _BASIS_CACHE_MISSES
    if not bool(SKETCH_BASIS_CACHE) or int(SKETCH_BASIS_CACHE_MAX) <= 0:
        _BASIS_CACHE_MISSES += 1
        return None
    val = _BASIS_CACHE.get(key)
    if val is None:
        _BASIS_CACHE_MISSES += 1
        return None
    _BASIS_CACHE_HITS += 1
    return val


def _basis_cache_store(key: tuple[Any, ...], value: Any) -> Any:
    if not bool(SKETCH_BASIS_CACHE) or int(SKETCH_BASIS_CACHE_MAX) <= 0:
        return value
    if len(_BASIS_CACHE) >= int(SKETCH_BASIS_CACHE_MAX):
        _BASIS_CACHE.clear()
    arr = np.asarray(value, dtype=np.complex128)
    arr.setflags(write=False)
    _BASIS_CACHE[key] = arr
    return arr


def basis_cache_stats() -> dict[str, Any]:
    return {
        "basis_cache_enabled": bool(SKETCH_BASIS_CACHE),
        "basis_cache_entries": int(len(_BASIS_CACHE)),
        "basis_cache_hits": int(_BASIS_CACHE_HITS),
        "basis_cache_misses": int(_BASIS_CACHE_MISSES),
        "basis_cache_max": int(SKETCH_BASIS_CACHE_MAX),
    }


def _sketch_basis_np(n: int, k: int, salt: int = 0) -> Any:
    nn, kk = int(n), int(k)
    if kk >= nn:
        return np.eye(nn, dtype=np.complex128)
    mode = "cached-rademacher" if SKETCH_MODE == "auto" else str(SKETCH_MODE)
    key = ("P", mode, int(SKETCH_SEED), nn, kk, int(salt))
    cached = _basis_cache_lookup(key)
    if cached is not None:
        return cached
    seed = int(splitmix64(int(SKETCH_SEED) + 1000003 * nn + 9176 * kk + 37 * int(salt)) & 0xFFFFFFFF)
    rng = np.random.default_rng(seed)
    if mode == "coordinate":
        idx = np.sort(rng.choice(nn, size=kk, replace=False))
        P = np.zeros((nn, kk), dtype=np.complex128)
        P[idx, np.arange(kk)] = 1.0
        return _basis_cache_store(key, P)
    if mode == "sparse-sign":
        # CountSketch-like columns: cheap to build and normalized, but not as
        # orthogonal as the dense QR Rademacher basis. Fallbacks remain enabled.
        width = max(4, min(nn, int(math.ceil(math.sqrt(float(nn))))))
        P = np.zeros((nn, kk), dtype=np.float64)
        scale = 1.0 / math.sqrt(float(width))
        for j in range(kk):
            idx = rng.choice(nn, size=width, replace=False)
            P[idx, j] = rng.choice([-scale, scale], size=width)
        return _basis_cache_store(key, P.astype(np.complex128, copy=False))
    raw = rng.choice([-1.0, 1.0], size=(nn, kk)).astype(np.float64) / math.sqrt(float(nn))
    # Real orthonormal columns keep P.T @ P ~= I, which matches the complex-linear
    # local-jet convention used by the polynomial oracle while still allowing
    # complex coefficients in delta=P u.
    qmat, _ = np.linalg.qr(raw, mode="reduced")
    return _basis_cache_store(key, qmat.astype(np.complex128, copy=False))


def _fast_projected_coded_basis_np(k: int, q: int, salt: int = 0) -> Any:
    kk, qq = int(k), int(q)
    if qq >= kk:
        return np.eye(kk, dtype=np.complex128)
    key = ("R", int(SKETCH_SEED), kk, qq, int(salt))
    cached = _basis_cache_lookup(key)
    if cached is not None:
        return cached
    seed = int(splitmix64(int(SKETCH_SEED) + 1000003 * kk + 9176 * qq + 53 * int(salt) + 0x411C0DE) & 0xFFFFFFFF)
    rng = np.random.default_rng(seed)
    raw = rng.choice([-1.0, 1.0], size=(kk, qq)).astype(np.float64) / math.sqrt(float(kk))
    qmat, _ = np.linalg.qr(raw, mode="reduced")
    return _basis_cache_store(key, qmat.astype(np.complex128, copy=False))


def _torch_sketch_basis(n: int, k: int, dtype: Any, device: Any, salt: int = 0) -> Any:
    return torch.as_tensor(_sketch_basis_np(int(n), int(k), int(salt)), dtype=dtype, device=device)


def _numpy_fast_projected_apply(A: Any, rhs: Any, damping: float = 0.0, condition_cap: float = 1e12) -> tuple[Any, dict[str, Any]]:
    AA = np.asarray(A, dtype=np.complex128)
    bb = np.asarray(rhs, dtype=np.complex128)
    if AA.ndim != 2 or bb.ndim != 1 or int(AA.shape[0]) != int(bb.shape[0]):
        raise ValueError("bad_fast_projected_shapes")
    col_norm2 = np.sum(np.abs(AA) ** 2, axis=0).real
    if col_norm2.size == 0 or not np.all(np.isfinite(col_norm2)):
        raise ValueError("nonfinite_fast_projected_norms")
    diag_mean = float(np.mean(col_norm2[col_norm2 > 0])) if np.any(col_norm2 > 0) else 1.0
    floor = max(1e-12, (1.0 / max(float(condition_cap), 1.0)) ** 2)
    mu = max(float(damping), floor) * max(diag_mean, 1e-300)
    denom = np.maximum(col_norm2 + mu, 1e-300)
    u = (AA.conj().T @ bb) / denom
    if not np.all(np.isfinite(u)):
        raise ValueError("nonfinite_fast_projected_solution")
    if int(AA.shape[1]) <= 192:
        G = AA.conj().T @ AA
        off = G - np.diag(np.diag(G))
        off_ratio = finite_norm(off) / max(1e-300, finite_norm(np.diag(np.diag(G))))
        diag_spread = float(np.max(col_norm2) / max(1e-300, np.min(col_norm2[col_norm2 > 0]))) if np.any(col_norm2 > 0) else float("inf")
    else:
        off_ratio = None
        diag_spread = None
    return u.astype(np.complex128, copy=False), {
        "reduced_solver": "fast-projected-diagonal-normal",
        "reduced_rank": int(AA.shape[1]),
        "reduced_damping_mu": float(mu),
        "reduced_fast_projector": True,
        "reduced_fast_projector_offdiag_ratio": off_ratio,
        "reduced_fast_projector_diag_spread": diag_spread,
    }


def _torch_reduced_svd_apply(A: Any, rhs: Any, rcond: float, condition_cap: float, damping: float = 0.0) -> tuple[Any, dict[str, Any]]:
    U, s, Vh = torch.linalg.svd(A, full_matrices=False)
    s_np = s.detach().cpu().numpy()
    cutoff, smax, smin = _svd_threshold(s_np, rcond, condition_cap)
    keep = torch.isfinite(s) & (s > float(cutoff))
    rank = int(keep.sum().detach().cpu().item())
    filt = torch.zeros_like(s)
    mu = 0.0
    if rank:
        if float(damping) > 0.0:
            mu_t = float(damping) * (torch.mean(s[keep] * s[keep]) + torch.as_tensor(1e-300, dtype=s.dtype, device=s.device))
            filt[keep] = s[keep] / (s[keep] * s[keep] + mu_t)
            mu = float(mu_t.detach().cpu().item())
        else:
            filt[keep] = 1.0 / s[keep]
        x = Vh.conj().transpose(-2, -1) @ (filt * (U.conj().transpose(-2, -1) @ rhs))
    else:
        x = torch.zeros(A.shape[1], dtype=A.dtype, device=A.device)
    return x, {
        "reduced_solver": "svd",
        "reduced_rank": int(rank),
        "reduced_smax": float(smax),
        "reduced_smin": float(smin),
        "reduced_cutoff": float(cutoff),
        "reduced_damping_mu": float(mu),
    }


def _torch_reduced_ns_apply(A: Any, rhs: Any, rcond: float, condition_cap: float, damping: float = 0.0) -> tuple[Any, dict[str, Any]]:
    Ah = A.conj().transpose(-2, -1)
    G = Ah @ A
    g = Ah @ rhs
    k = int(G.shape[-1])
    diag_mean = torch.mean(torch.real(torch.diagonal(G))).clamp_min(torch.finfo(torch.float32).tiny)
    floor = max(float(NS_DAMPING_FLOOR), (1.0 / max(float(condition_cap), 1.0)) ** 2)
    mu = (max(float(damping), floor) * diag_mean).to(A.dtype)
    x = _torch_ns_inverse(G + mu * torch.eye(k, dtype=A.dtype, device=A.device), NS_ITERS) @ g
    return x, {
        "reduced_solver": "newton-schulz",
        "reduced_rank": int(k),
        "reduced_damping_mu": float(torch.real(mu).detach().cpu().item()),
    }


def _torch_reduced_apply(A: Any, rhs: Any, rcond: float, condition_cap: float, damping: float = 0.0) -> tuple[Any, dict[str, Any]]:
    if SKETCH_SOLVER == "ns":
        return _torch_reduced_ns_apply(A, rhs, rcond, condition_cap, damping)
    try:
        return _torch_reduced_svd_apply(A, rhs, rcond, condition_cap, damping)
    except Exception as exc:
        if SKETCH_SOLVER != "auto":
            raise
        x, meta = _torch_reduced_ns_apply(A, rhs, rcond, condition_cap, damping)
        meta["reduced_solver_fallback"] = f"svd-error:{type(exc).__name__}"
        return x, meta


def torch_sketch_matrix_step(A: Any, rhs: Any, rcond: float = 1e-12, condition_cap: float = 1e12, damping: float = 0.0) -> tuple[Any, dict[str, Any]]:
    AA_np = np.asarray(A, dtype=np.complex128)
    bb_np = np.asarray(rhs, dtype=np.complex128)
    n = int(AA_np.shape[1])
    k = sketch_dim_for(n)
    if k >= n:
        x, meta = torch_gemm_matrix_step(AA_np, bb_np, rcond, condition_cap, damping)
        meta["matrix_algorithm"] = "sketch-bypass-full-gemm-newton-schulz"
        meta["matrix_sketch_dim"] = int(k)
        meta["matrix_sketch_mode"] = "full"
        return x, meta
    device, dtype = torch_device(), torch_complex_dtype()
    AA = torch.as_tensor(AA_np, dtype=dtype, device=device)
    bb = torch.as_tensor(bb_np, dtype=dtype, device=device)
    P = _torch_sketch_basis(n, k, dtype, device, salt=int(AA_np.shape[0]))
    AP = AA @ P
    z, rmeta = _torch_reduced_apply(AP, bb, float(rcond), float(condition_cap), float(damping))
    x = P @ z
    x_np = x.detach().cpu().numpy().astype(np.complex128, copy=False)
    lin = AA_np @ x_np - bb_np
    normal_res = AA_np.conj().T @ lin
    meta = {
        "matrix_backend": "torch",
        "matrix_algorithm": "sketch-newton-schulz",
        "torch_device": TORCH_RESOLVED_DEVICE,
        "torch_complex_dtype": TORCH_RESOLVED_DTYPE,
        "matrix_method": "subspace-sketch-normal-equations-newton-schulz",
        "matrix_dim": int(n),
        "matrix_sketch_dim": int(k),
        "matrix_sketch_ratio": float(k / max(1, n)),
        "matrix_sketch_mode": str(SKETCH_MODE),
        "matrix_sketch_seed": int(SKETCH_SEED),
        "matrix_sketch_solver": str(SKETCH_SOLVER),
        "matrix_ns_iters": int(NS_ITERS),
        "matrix_damping_mu": float(rmeta.get("reduced_damping_mu", 0.0)),
        "matrix_damping_floor": float(NS_DAMPING_FLOOR),
        "matrix_rcond": float(rcond),
        "matrix_condition_cap": float(condition_cap),
        "matrix_rhs_norm": finite_norm(bb_np),
        "matrix_step_norm": finite_norm(x_np),
        "matrix_linear_residual": finite_norm(lin),
        "matrix_normal_residual": finite_norm(normal_res),
        **rmeta,
    }
    return x_np, meta


def torch_sketch_least_squares_fit(X: Any, Y: Any, rcond: float = 1e-12, condition_cap: float = 1e12) -> tuple[Any, dict[str, Any]]:
    XX_np = np.asarray(X, dtype=np.complex128)
    YY_np = np.asarray(Y, dtype=np.complex128)
    n = int(XX_np.shape[1])
    k = sketch_dim_for(n)
    if k >= n:
        B, meta = torch_gemm_least_squares_fit(XX_np, YY_np, rcond, condition_cap)
        meta["matrix_algorithm"] = "sketch-bypass-full-gemm-newton-schulz"
        meta["matrix_fit_sketch_dim"] = int(k)
        meta["matrix_fit_sketch_mode"] = "full"
        return B, meta
    device, dtype = torch_device(), torch_complex_dtype()
    XX = torch.as_tensor(XX_np, dtype=dtype, device=device)
    YY = torch.as_tensor(YY_np, dtype=dtype, device=device)
    P = _torch_sketch_basis(n, k, dtype, device, salt=int(XX_np.shape[0]) + 11)
    XP = XX @ P
    if SKETCH_SOLVER == "ns":
        Xh = XP.conj().transpose(-2, -1)
        G = Xh @ XP
        RHS = Xh @ YY
        diag_mean = torch.mean(torch.real(torch.diagonal(G))).clamp_min(torch.finfo(torch.float32).tiny)
        floor = max(float(NS_DAMPING_FLOOR), (1.0 / max(float(condition_cap), 1.0)) ** 2)
        mu = (floor * diag_mean).to(dtype)
        C = _torch_ns_inverse(G + mu * torch.eye(k, dtype=dtype, device=device), NS_ITERS) @ RHS
        rmeta = {"reduced_solver": "newton-schulz", "reduced_damping_mu": float(torch.real(mu).detach().cpu().item())}
    else:
        U, s, Vh = torch.linalg.svd(XP, full_matrices=False)
        s_np = s.detach().cpu().numpy()
        cutoff, smax, smin = _svd_threshold(s_np, rcond, condition_cap)
        keep = torch.isfinite(s) & (s > float(cutoff))
        inv = torch.zeros_like(s)
        inv[keep] = 1.0 / s[keep]
        C = Vh.conj().transpose(-2, -1) @ (inv[:, None] * (U.conj().transpose(-2, -1) @ YY))
        rmeta = {"reduced_solver": "svd", "reduced_rank": int(keep.sum().detach().cpu().item()), "reduced_smax": float(smax), "reduced_smin": float(smin), "reduced_cutoff": float(cutoff), "reduced_damping_mu": 0.0}
    B = P @ C
    B_np = B.detach().cpu().numpy().astype(np.complex128, copy=False)
    fit_res = XX_np @ B_np - YY_np
    meta = {
        "matrix_backend": "torch",
        "matrix_algorithm": "sketch-newton-schulz",
        "torch_device": TORCH_RESOLVED_DEVICE,
        "torch_complex_dtype": TORCH_RESOLVED_DTYPE,
        "matrix_fit_method": "subspace-sketch-normal-equations-newton-schulz",
        "matrix_fit_rows": int(XX_np.shape[0]),
        "matrix_fit_cols": int(n),
        "matrix_fit_sketch_dim": int(k),
        "matrix_fit_sketch_ratio": float(k / max(1, n)),
        "matrix_fit_sketch_mode": str(SKETCH_MODE),
        "matrix_fit_sketch_solver": str(SKETCH_SOLVER),
        "matrix_fit_ns_iters": int(NS_ITERS),
        "matrix_fit_damping_mu": float(rmeta.get("reduced_damping_mu", 0.0)),
        "matrix_fit_residual": finite_norm(fit_res),
        **rmeta,
    }
    return B_np, meta


def _torch_real_tiny(dtype: Any) -> float:
    if dtype in {torch.float16, torch.bfloat16}:
        return 1e-6
    return float(torch.finfo(dtype).tiny)


def _real_tensor_to_numpy(x: Any) -> Any:
    if x.dtype in {torch.float16, torch.bfloat16}:
        x = x.to(dtype=torch.float32)
    return x.detach().cpu().numpy().astype(np.float64, copy=False)


def _torch_pack_complex_operator(A: Any, dtype: Any, device: Any) -> Any:
    """Column-vector complex operator A x -> real block [[Re,-Im],[Im,Re]]."""
    AA = np.asarray(A, dtype=np.complex128)
    R = torch.as_tensor(AA.real, dtype=dtype, device=device)
    I = torch.as_tensor(AA.imag, dtype=dtype, device=device)
    return torch.cat((torch.cat((R, -I), dim=-1), torch.cat((I, R), dim=-1)), dim=-2)


def _torch_pack_complex_design(X: Any, dtype: Any, device: Any) -> Any:
    """Row-vector design X B -> Y, stacked as [ReY; ImY]."""
    XX = np.asarray(X, dtype=np.complex128)
    R = torch.as_tensor(XX.real, dtype=dtype, device=device)
    I = torch.as_tensor(XX.imag, dtype=dtype, device=device)
    return torch.cat((torch.cat((R, -I), dim=-1), torch.cat((I, R), dim=-1)), dim=-2)


def _torch_pack_complex_operator_batched(A: Any, dtype: Any, device: Any) -> Any:
    AA = np.asarray(A, dtype=np.complex128)
    R = torch.as_tensor(AA.real, dtype=dtype, device=device)
    I = torch.as_tensor(AA.imag, dtype=dtype, device=device)
    return torch.cat((torch.cat((R, -I), dim=-1), torch.cat((I, R), dim=-1)), dim=-2)


def _torch_pack_complex_design_batched(X: Any, dtype: Any, device: Any) -> Any:
    XX = np.asarray(X, dtype=np.complex128)
    R = torch.as_tensor(XX.real, dtype=dtype, device=device)
    I = torch.as_tensor(XX.imag, dtype=dtype, device=device)
    return torch.cat((torch.cat((R, -I), dim=-1), torch.cat((I, R), dim=-1)), dim=-2)


def _torch_realpack_normal_step(A: Any, rhs: Any, damping: float = 0.0, condition_cap: float = 1e12) -> tuple[Any, Any, Any]:
    Ah = A.transpose(-2, -1)
    G = Ah @ A
    g = Ah @ rhs
    n2 = int(G.shape[-1])
    diag_mean = torch.mean(torch.diagonal(G)).clamp_min(_torch_real_tiny(A.dtype))
    floor = max(float(NS_DAMPING_FLOOR), (1.0 / max(float(condition_cap), 1.0)) ** 2)
    mu = max(float(damping), floor) * diag_mean
    B = G + mu * torch.eye(n2, dtype=A.dtype, device=A.device)
    return _torch_ns_inverse(B, NS_ITERS) @ g, mu, G


def _torch_realpack_batched_normal_step(A: Any, rhs: Any, damping: float = 0.0) -> tuple[Any, Any]:
    Ah = A.transpose(-2, -1)
    G = torch.bmm(Ah, A)
    g = torch.bmm(Ah, rhs[:, :, None])
    n2 = int(G.shape[-1])
    diag = torch.diagonal(G, dim1=-2, dim2=-1).mean(dim=-1).clamp_min(_torch_real_tiny(A.dtype))
    mu = (max(float(damping), float(NS_DAMPING_FLOOR)) * diag)[:, None, None]
    I = torch.eye(n2, dtype=A.dtype, device=A.device)[None, :, :]
    return torch.bmm(_torch_batched_ns_inverse(G + mu * I, NS_ITERS), g)[:, :, 0], mu


def torch_realpack_matrix_step(A: Any, rhs: Any, rcond: float = 1e-12, condition_cap: float = 1e12, damping: float = 0.0) -> tuple[Any, dict[str, Any]]:
    AA_np = np.asarray(A, dtype=np.complex128)
    bb_np = np.asarray(rhs, dtype=np.complex128)
    device, dtype = torch_device(), torch_real_dtype()
    Ar = _torch_pack_complex_operator(AA_np, dtype, device)
    br = torch.as_tensor(np.concatenate((bb_np.real, bb_np.imag)), dtype=dtype, device=device)
    xr, mu, _G = _torch_realpack_normal_step(Ar, br, float(damping), float(condition_cap))
    x_pack = _real_tensor_to_numpy(xr)
    n = int(AA_np.shape[1])
    x_np = x_pack[:n] + 1j * x_pack[n:]
    lin = AA_np @ x_np - bb_np
    normal_res = (AA_np.conj().T @ lin) + complex(float(_real_tensor_to_numpy(mu.reshape(1))[0])) * x_np
    meta = {
        "matrix_backend": "torch",
        "matrix_algorithm": "realpack-newton-schulz",
        "torch_device": TORCH_RESOLVED_DEVICE,
        "torch_real_dtype": TORCH_RESOLVED_REAL_DTYPE,
        "matrix_method": "real-packed-normal-equations-newton-schulz",
        "matrix_rank": int(2 * n),
        "matrix_dim": int(n),
        "matrix_real_dim": int(2 * n),
        "matrix_ns_iters": int(NS_ITERS),
        "matrix_damping_mu": float(_real_tensor_to_numpy(mu.reshape(1))[0]),
        "matrix_damping_floor": float(NS_DAMPING_FLOOR),
        "matrix_rcond": float(rcond),
        "matrix_condition_cap": float(condition_cap),
        "matrix_rhs_norm": finite_norm(bb_np),
        "matrix_step_norm": finite_norm(x_np),
        "matrix_linear_residual": finite_norm(lin),
        "matrix_normal_residual": finite_norm(normal_res),
    }
    return x_np, meta


def torch_realpack_least_squares_fit(X: Any, Y: Any, rcond: float = 1e-12, condition_cap: float = 1e12) -> tuple[Any, dict[str, Any]]:
    XX_np = np.asarray(X, dtype=np.complex128)
    YY_np = np.asarray(Y, dtype=np.complex128)
    device, dtype = torch_device(), torch_real_dtype()
    Xr = _torch_pack_complex_design(XX_np, dtype, device)
    Yr = torch.as_tensor(np.concatenate((YY_np.real, YY_np.imag), axis=0), dtype=dtype, device=device)
    Xh = Xr.transpose(-2, -1)
    G = Xh @ Xr
    RHS = Xh @ Yr
    n2 = int(G.shape[-1])
    n = n2 // 2
    diag_mean = torch.mean(torch.diagonal(G)).clamp_min(_torch_real_tiny(dtype))
    floor = max(float(NS_DAMPING_FLOOR), (1.0 / max(float(condition_cap), 1.0)) ** 2)
    mu = floor * diag_mean
    coef = _torch_ns_inverse(G + mu * torch.eye(n2, dtype=dtype, device=device), NS_ITERS) @ RHS
    C = _real_tensor_to_numpy(coef)
    B_np = C[:n, :] + 1j * C[n:, :]
    fit_res = XX_np @ B_np - YY_np
    meta = {
        "matrix_backend": "torch",
        "matrix_algorithm": "realpack-newton-schulz",
        "torch_device": TORCH_RESOLVED_DEVICE,
        "torch_real_dtype": TORCH_RESOLVED_REAL_DTYPE,
        "matrix_fit_method": "real-packed-normal-equations-newton-schulz",
        "matrix_fit_rank": int(n2),
        "matrix_fit_rows": int(XX_np.shape[0]),
        "matrix_fit_cols": int(XX_np.shape[1]),
        "matrix_fit_real_rows": int(2 * XX_np.shape[0]),
        "matrix_fit_real_cols": int(n2),
        "matrix_fit_ns_iters": int(NS_ITERS),
        "matrix_fit_damping_mu": float(_real_tensor_to_numpy(mu.reshape(1))[0]),
        "matrix_fit_residual": finite_norm(fit_res),
    }
    return B_np, meta


def torch_gemm_matrix_step(A: Any, rhs: Any, rcond: float = 1e-12, condition_cap: float = 1e12, damping: float = 0.0) -> tuple[Any, dict[str, Any]]:
    AA_np = np.asarray(A, dtype=np.complex128)
    bb_np = np.asarray(rhs, dtype=np.complex128)
    device = torch_device()
    dtype = torch_complex_dtype()
    AA = torch.as_tensor(AA_np, dtype=dtype, device=device)
    bb = torch.as_tensor(bb_np, dtype=dtype, device=device)
    Ah = AA.conj().transpose(-2, -1)
    G = Ah @ AA
    g = Ah @ bb
    n = int(G.shape[-1])
    diag_mean = torch.mean(torch.real(torch.diagonal(G))).clamp_min(torch.finfo(torch.float32).tiny)
    floor = max(float(NS_DAMPING_FLOOR), (1.0 / max(float(condition_cap), 1.0)) ** 2)
    mu = (max(float(damping), floor) * diag_mean).to(dtype)
    I = torch.eye(n, dtype=dtype, device=device)
    B = G + mu * I
    Binv = _torch_ns_inverse(B, NS_ITERS)
    x = Binv @ g
    x_np = x.detach().cpu().numpy().astype(np.complex128, copy=False)
    lin = AA_np @ x_np - bb_np
    normal_res = (AA_np.conj().T @ lin) + complex(float(torch.real(mu).detach().cpu().item())) * x_np
    meta = {
        "matrix_backend": "torch",
        "matrix_algorithm": "gemm-newton-schulz",
        "torch_device": TORCH_RESOLVED_DEVICE,
        "torch_complex_dtype": TORCH_RESOLVED_DTYPE,
        "matrix_method": "normal-equations-newton-schulz",
        "matrix_rank": int(n),
        "matrix_dim": int(n),
        "matrix_ns_iters": int(NS_ITERS),
        "matrix_damping_mu": float(torch.real(mu).detach().cpu().item()),
        "matrix_damping_floor": float(NS_DAMPING_FLOOR),
        "matrix_rcond": float(rcond),
        "matrix_condition_cap": float(condition_cap),
        "matrix_rhs_norm": finite_norm(bb_np),
        "matrix_step_norm": finite_norm(x_np),
        "matrix_linear_residual": finite_norm(lin),
        "matrix_normal_residual": finite_norm(normal_res),
    }
    return x_np, meta


def torch_gemm_least_squares_fit(X: Any, Y: Any, rcond: float = 1e-12, condition_cap: float = 1e12) -> tuple[Any, dict[str, Any]]:
    XX_np = np.asarray(X, dtype=np.complex128)
    YY_np = np.asarray(Y, dtype=np.complex128)
    device = torch_device()
    dtype = torch_complex_dtype()
    XX = torch.as_tensor(XX_np, dtype=dtype, device=device)
    YY = torch.as_tensor(YY_np, dtype=dtype, device=device)
    Xh = XX.conj().transpose(-2, -1)
    G = Xh @ XX
    RHS = Xh @ YY
    n = int(G.shape[-1])
    diag_mean = torch.mean(torch.real(torch.diagonal(G))).clamp_min(torch.finfo(torch.float32).tiny)
    floor = max(float(NS_DAMPING_FLOOR), (1.0 / max(float(condition_cap), 1.0)) ** 2)
    mu = (floor * diag_mean).to(dtype)
    I = torch.eye(n, dtype=dtype, device=device)
    Binv = _torch_ns_inverse(G + mu * I, NS_ITERS)
    B = Binv @ RHS
    B_np = B.detach().cpu().numpy().astype(np.complex128, copy=False)
    fit_res = XX_np @ B_np - YY_np
    meta = {
        "matrix_backend": "torch",
        "matrix_algorithm": "gemm-newton-schulz",
        "torch_device": TORCH_RESOLVED_DEVICE,
        "torch_complex_dtype": TORCH_RESOLVED_DTYPE,
        "matrix_fit_method": "normal-equations-newton-schulz",
        "matrix_fit_rank": int(n),
        "matrix_fit_rows": int(XX_np.shape[0]),
        "matrix_fit_cols": int(XX_np.shape[1]),
        "matrix_fit_ns_iters": int(NS_ITERS),
        "matrix_fit_damping_mu": float(torch.real(mu).detach().cpu().item()),
        "matrix_fit_residual": finite_norm(fit_res),
    }
    return B_np, meta


def torch_matrix_svd_step(A: Any, rhs: Any, rcond: float = 1e-12, condition_cap: float = 1e12, damping: float = 0.0) -> tuple[Any, dict[str, Any]]:
    AA_np = np.asarray(A, dtype=np.complex128)
    bb_np = np.asarray(rhs, dtype=np.complex128)
    device = torch_device()
    dtype = torch_complex_dtype()
    AA = torch.as_tensor(AA_np, dtype=dtype, device=device)
    bb = torch.as_tensor(bb_np, dtype=dtype, device=device)
    U, s, Vh = torch.linalg.svd(AA, full_matrices=False)
    s_np = s.detach().cpu().numpy()
    cutoff, smax, smin = _svd_threshold(s_np, rcond, condition_cap)
    keep = torch.isfinite(s) & (s > float(cutoff))
    rank = int(keep.sum().detach().cpu().item())
    mu = 0.0
    filt = torch.zeros_like(s)
    if rank:
        if float(damping) > 0.0:
            mu_t = float(damping) * (torch.mean(s[keep] * s[keep]) + torch.as_tensor(1e-300, dtype=s.dtype, device=s.device))
            filt[keep] = s[keep] / (s[keep] * s[keep] + mu_t)
            mu = float(mu_t.detach().cpu().item())
        else:
            filt[keep] = 1.0 / s[keep]
        x = Vh.conj().transpose(-2, -1) @ (filt * (U.conj().transpose(-2, -1) @ bb))
    else:
        x = torch.zeros(AA.shape[1], dtype=dtype, device=device)
    x_np = x.detach().cpu().numpy().astype(np.complex128, copy=False)
    cond = float(smax / max(smin, 1e-300)) if smax > 0.0 else float("inf")
    kept = s_np[np.isfinite(s_np) & (s_np > cutoff)]
    kept_cond = float(smax / max(float(np.min(kept)), 1e-300)) if kept.size else float("inf")
    meta = {
        "matrix_backend": "torch",
        "torch_device": TORCH_RESOLVED_DEVICE,
        "torch_complex_dtype": TORCH_RESOLVED_DTYPE,
        "matrix_method": "torch-svd-filter",
        "matrix_rank": int(rank),
        "matrix_dim": int(AA_np.shape[1]),
        "matrix_smax": float(smax),
        "matrix_smin": float(smin),
        "matrix_condition": cond,
        "matrix_kept_condition": kept_cond,
        "matrix_cutoff": float(cutoff),
        "matrix_rcond": float(rcond),
        "matrix_condition_cap": float(condition_cap),
        "matrix_damping_mu": float(mu),
        "matrix_rhs_norm": finite_norm(bb_np),
        "matrix_step_norm": finite_norm(x_np),
        "matrix_linear_residual": finite_norm(AA_np @ x_np - bb_np),
    }
    return x_np, meta


def torch_matrix_least_squares_fit(X: Any, Y: Any, rcond: float = 1e-12, condition_cap: float = 1e12) -> tuple[Any, dict[str, Any]]:
    XX_np = np.asarray(X, dtype=np.complex128)
    YY_np = np.asarray(Y, dtype=np.complex128)
    device = torch_device()
    dtype = torch_complex_dtype()
    XX = torch.as_tensor(XX_np, dtype=dtype, device=device)
    YY = torch.as_tensor(YY_np, dtype=dtype, device=device)
    U, s, Vh = torch.linalg.svd(XX, full_matrices=False)
    s_np = s.detach().cpu().numpy()
    cutoff, smax, smin = _svd_threshold(s_np, rcond, condition_cap)
    keep = torch.isfinite(s) & (s > float(cutoff))
    inv = torch.zeros_like(s)
    inv[keep] = 1.0 / s[keep]
    B = Vh.conj().transpose(-2, -1) @ (inv[:, None] * (U.conj().transpose(-2, -1) @ YY))
    B_np = B.detach().cpu().numpy().astype(np.complex128, copy=False)
    rank = int(keep.sum().detach().cpu().item())
    cond = float(smax / max(smin, 1e-300)) if smax > 0.0 else float("inf")
    meta = {
        "matrix_backend": "torch",
        "torch_device": TORCH_RESOLVED_DEVICE,
        "torch_complex_dtype": TORCH_RESOLVED_DTYPE,
        "matrix_fit_method": "torch-svd-pseudoinverse",
        "matrix_fit_rank": int(rank),
        "matrix_fit_rows": int(XX_np.shape[0]),
        "matrix_fit_cols": int(XX_np.shape[1]),
        "matrix_fit_smax": float(smax),
        "matrix_fit_smin": float(smin),
        "matrix_fit_condition": cond,
        "matrix_fit_cutoff": float(cutoff),
    }
    return B_np, meta


def matrix_svd_step(A: Any, rhs: Any, rcond: float = 1e-12, condition_cap: float = 1e12, damping: float = 0.0) -> tuple[Any, dict[str, Any]]:
    """Apply the configured matrix pseudo-inverse to rhs.

    In 410 the default path is real-packed GEMM/Newton-Schulz on regularised
    normal equations. ``--matrix-algorithm gemm-ns`` preserves the 404 complex
    GEMM path and ``--matrix-algorithm svd`` preserves the 402 SVD filter.
    """
    if MATRIX_BACKEND in {"torch", "auto"}:
        try:
            selected_algorithm = selected_matrix_algorithm_410(int(np.asarray(A).shape[-1]))
            if selected_algorithm == "sketch-ns":
                return torch_sketch_matrix_step(A, rhs, rcond, condition_cap, damping)
            if selected_algorithm in {"realpack-ns", "auto"}:
                return torch_realpack_matrix_step(A, rhs, rcond, condition_cap, damping)
            if selected_algorithm == "gemm-ns":
                return torch_gemm_matrix_step(A, rhs, rcond, condition_cap, damping)
            return torch_matrix_svd_step(A, rhs, rcond, condition_cap, damping)
        except Exception as exc:
            if MATRIX_ALGORITHM in {"auto", "adaptive-sketch", "adaptive-ns"}:
                if selected_algorithm != "gemm-ns":
                    try:
                        x, meta = torch_gemm_matrix_step(A, rhs, rcond, condition_cap, damping)
                        meta["matrix_algorithm_fallback"] = f"{selected_algorithm}-error:{type(exc).__name__}"
                        return x, meta
                    except Exception:
                        pass
                try:
                    x, meta = torch_matrix_svd_step(A, rhs, rcond, condition_cap, damping)
                    meta["matrix_algorithm_fallback"] = f"{selected_algorithm}-error:{type(exc).__name__}"
                    return x, meta
                except Exception:
                    pass
            if MATRIX_BACKEND == "torch":
                raise
            x, meta = numpy_matrix_svd_step(A, rhs, rcond, condition_cap, damping)
            meta["matrix_backend_fallback"] = f"torch-error:{type(exc).__name__}"
            return x, meta
    return numpy_matrix_svd_step(A, rhs, rcond, condition_cap, damping)


def numpy_matrix_svd_step(A: Any, rhs: Any, rcond: float = 1e-12, condition_cap: float = 1e12, damping: float = 0.0) -> tuple[Any, dict[str, Any]]:
    AA = np.asarray(A, dtype=np.complex128)
    bb = np.asarray(rhs, dtype=np.complex128)
    U, s, Vh = np.linalg.svd(AA, full_matrices=False)
    cutoff, smax, smin = _svd_threshold(s, rcond, condition_cap)
    keep = np.asarray(np.isfinite(s) & (s > cutoff), dtype=bool)
    rank = int(np.sum(keep))
    mu = 0.0
    filt = np.zeros_like(s, dtype=np.float64)
    if rank:
        if float(damping) > 0.0:
            mu = float(damping) * (float(np.mean(s[keep] * s[keep])) + 1e-300)
            filt[keep] = s[keep] / (s[keep] * s[keep] + mu)
        else:
            filt[keep] = 1.0 / s[keep]
        x = Vh.conj().T @ (filt * (U.conj().T @ bb))
    else:
        x = np.zeros(AA.shape[1], dtype=np.complex128)
    cond = float(smax / max(smin, 1e-300)) if smax > 0.0 else float("inf")
    kept_cond = float(smax / max(float(np.min(s[keep])), 1e-300)) if rank else float("inf")
    meta = {
        "matrix_backend": "numpy",
        "matrix_method": "svd-filter",
        "matrix_rank": int(rank),
        "matrix_dim": int(AA.shape[1]),
        "matrix_smax": float(smax),
        "matrix_smin": float(smin),
        "matrix_condition": cond,
        "matrix_kept_condition": kept_cond,
        "matrix_cutoff": float(cutoff),
        "matrix_rcond": float(rcond),
        "matrix_condition_cap": float(condition_cap),
        "matrix_damping_mu": float(mu),
        "matrix_rhs_norm": finite_norm(bb),
        "matrix_step_norm": finite_norm(x),
        "matrix_linear_residual": finite_norm(AA @ x - bb) if AA.ndim == 2 else float("inf"),
    }
    return x, meta


def matrix_least_squares_fit(X: Any, Y: Any, rcond: float = 1e-12, condition_cap: float = 1e12) -> tuple[Any, dict[str, Any]]:
    """Fit B in X @ B ~= Y by the configured matrix pseudo-inverse."""
    if MATRIX_BACKEND in {"torch", "auto"}:
        try:
            selected_algorithm = selected_matrix_algorithm_410(int(np.asarray(X).shape[-1]))
            if selected_algorithm == "sketch-ns":
                return torch_sketch_least_squares_fit(X, Y, rcond, condition_cap)
            if selected_algorithm in {"realpack-ns", "auto"}:
                return torch_realpack_least_squares_fit(X, Y, rcond, condition_cap)
            if selected_algorithm == "gemm-ns":
                return torch_gemm_least_squares_fit(X, Y, rcond, condition_cap)
            return torch_matrix_least_squares_fit(X, Y, rcond, condition_cap)
        except Exception as exc:
            if MATRIX_ALGORITHM in {"auto", "adaptive-sketch", "adaptive-ns"}:
                if selected_algorithm != "gemm-ns":
                    try:
                        B, meta = torch_gemm_least_squares_fit(X, Y, rcond, condition_cap)
                        meta["matrix_algorithm_fallback"] = f"{selected_algorithm}-error:{type(exc).__name__}"
                        return B, meta
                    except Exception:
                        pass
                try:
                    B, meta = torch_matrix_least_squares_fit(X, Y, rcond, condition_cap)
                    meta["matrix_algorithm_fallback"] = f"{selected_algorithm}-error:{type(exc).__name__}"
                    return B, meta
                except Exception:
                    pass
            if MATRIX_BACKEND == "torch":
                raise
            B, meta = numpy_matrix_least_squares_fit(X, Y, rcond, condition_cap)
            meta["matrix_backend_fallback"] = f"torch-error:{type(exc).__name__}"
            return B, meta
    return numpy_matrix_least_squares_fit(X, Y, rcond, condition_cap)


def numpy_matrix_least_squares_fit(X: Any, Y: Any, rcond: float = 1e-12, condition_cap: float = 1e12) -> tuple[Any, dict[str, Any]]:
    XX = np.asarray(X, dtype=np.complex128)
    YY = np.asarray(Y, dtype=np.complex128)
    U, s, Vh = np.linalg.svd(XX, full_matrices=False)
    cutoff, smax, smin = _svd_threshold(s, rcond, condition_cap)
    keep = np.asarray(np.isfinite(s) & (s > cutoff), dtype=bool)
    inv = np.zeros_like(s, dtype=np.float64)
    inv[keep] = 1.0 / s[keep]
    B = Vh.conj().T @ (inv[:, None] * (U.conj().T @ YY))
    rank = int(np.sum(keep))
    cond = float(smax / max(smin, 1e-300)) if smax > 0.0 else float("inf")
    meta = {
        "matrix_backend": "numpy",
        "matrix_fit_method": "svd-pseudoinverse",
        "matrix_fit_rank": int(rank),
        "matrix_fit_rows": int(XX.shape[0]),
        "matrix_fit_cols": int(XX.shape[1]),
        "matrix_fit_smax": float(smax),
        "matrix_fit_smin": float(smin),
        "matrix_fit_condition": cond,
        "matrix_fit_cutoff": float(cutoff),
    }
    return B, meta


def finite_slope_matrix(system: Any, a: Sequence[complex], b: Sequence[complex]) -> Any:
    t0 = now()
    aa = np.asarray(a, dtype=np.complex128)
    bb = np.asarray(b, dtype=np.complex128)
    cur = aa.copy()
    f_prev = system.eval(cur)
    Q = np.zeros((int(system.n), int(system.n)), dtype=np.complex128)
    for j in range(int(system.n)):
        old = cur[j]
        cur[j] = bb[j]
        f_next = system.eval(cur)
        dz = bb[j] - old
        if abs(dz) > 1e-300:
            Q[:, j] = (f_next - f_prev) / dz
        else:
            h = 1e-6 * max(1.0, abs(old))
            plus, minus = cur.copy(), cur.copy()
            plus[j] = old + h
            minus[j] = old - h
            Q[:, j] = (system.eval(plus) - system.eval(minus)) / (2.0 * h)
        f_prev = f_next
    system.slope_count = int(getattr(system, "slope_count", 0)) + 1
    system.seconds_slope = float(getattr(system, "seconds_slope", 0.0)) + (now() - t0)
    return Q


# ---------------------------------------------------------------------------
# KS systems (used only as a black-box geometric oracle)
# ---------------------------------------------------------------------------

def compositions_leq(d: int, n: int) -> Any:
    out: list[tuple[int, ...]] = []

    def rec(pos: int, rem: int, cur: list[int]) -> None:
        if pos == n - 1:
            for k in range(rem + 1):
                out.append(tuple(cur + [k]))
            return
        for k in range(rem + 1):
            cur.append(k)
            rec(pos + 1, rem - k, cur)
            cur.pop()

    rec(0, int(d), [])
    return np.asarray(out, dtype=np.int16 if d < 32767 else np.int32)


def multinomial_kostlan_weights(exps: Any, d: int) -> Any:
    totals = np.sum(exps, axis=1).astype(np.int64)
    logfac = np.zeros(int(d) + 1, dtype=np.float64)
    for k in range(1, int(d) + 1):
        logfac[k] = logfac[k - 1] + math.log(k)
    logs = logfac[d] - logfac[d - totals]
    for j in range(exps.shape[1]):
        logs -= logfac[exps[:, j].astype(np.int64)]
    return np.exp(0.5 * logs)


@dataclasses.dataclass
class DenseKostlanSystem:
    n: int
    d: int
    seed: int
    exps: Any
    coeff: Any
    weights: Any
    equation_normalize: bool = True
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, seed_index: int = 0, equation_normalize: bool = True) -> "DenseKostlanSystem":
        t0 = now()
        exps = compositions_leq(d, n)
        weights = multinomial_kostlan_weights(exps, d)
        seed = stable_seed(n, d, seed_index)
        rng = np.random.default_rng(seed)
        coeff = (rng.standard_normal((n, exps.shape[0])) + 1j * rng.standard_normal((n, exps.shape[0]))) / math.sqrt(2.0)
        coeff = coeff * weights[None, :]
        if equation_normalize:
            row = np.linalg.norm(coeff, axis=1)
            coeff = coeff / np.where(row > 0, row, 1.0)[:, None]
        obj = cls(int(n), int(d), seed, exps, coeff.astype(np.complex128), weights, bool(equation_normalize))
        obj._generation_seconds = now() - t0
        return obj

    @property
    def terms_per_poly(self) -> int:
        return int(self.exps.shape[0])

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return float(getattr(self, "_generation_seconds", 0.0))

    def monomials_batch(self, Z: Any) -> Any:
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        B, M = int(ZZ.shape[0]), self.terms_per_poly
        mon = np.ones((B, M), dtype=np.complex128)
        for j in range(self.n):
            p = np.empty((B, self.d + 1), dtype=np.complex128)
            p[:, 0] = 1.0
            if self.d > 0:
                p[:, 1] = ZZ[:, j]
                for k in range(2, self.d + 1):
                    p[:, k] = p[:, k - 1] * ZZ[:, j]
            mon *= p[:, self.exps[:, j]]
        return mon

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = (self.monomials_batch(np.asarray(z, dtype=np.complex128))[0] @ self.coeff.T)
        self.eval_count += 1
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        out = self.monomials_batch(ZZ) @ self.coeff.T
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms}


@dataclasses.dataclass
class LazyFeatureKostlanSystem:
    n: int
    d: int
    seed: int
    feature_exps: Any
    feature_exps_t: Any
    feature_log_scales: Any
    coeff: Any
    lazy_features: int
    projective_normalize: bool = True
    dynamic_normalize: bool = True
    equation_normalize: bool = True
    eval_block: int = 128
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, seed_index: int = 0, equation_normalize: bool = True, lazy_features: int = 0, lazy_feature_cap: int = 8192, projective_normalize: bool = True, dynamic_normalize: bool = True, eval_block: int = 128) -> "LazyFeatureKostlanSystem":
        t0 = now()
        n, d = int(n), int(d)
        m = max(n + 1, int(lazy_features) if int(lazy_features) > 0 else auto_lazy_feature_count(n, d, int(lazy_feature_cap)))
        seed = stable_seed(n, d, seed_index, salt=0x314314)
        rng = np.random.default_rng(seed)
        exps = np.zeros((m, n), dtype=np.int16 if d < 32767 else np.int32)
        degrees = np.zeros(m, dtype=np.int64)
        idx = 1
        for j in range(n):
            if idx >= m:
                break
            exps[idx, j] = 1 if d >= 1 else 0
            degrees[idx] = 1 if d >= 1 else 0
            idx += 1
        probs = np.full(n, 1.0 / max(1, n))
        while idx < m:
            q = (idx - (n + 1) + 0.5) / max(1, m - (n + 1))
            k = int(min(d, max(0, math.floor(q * (d + 1)))))
            if d > 0 and idx % 7 == 0:
                k = int(rng.integers(0, d + 1))
            if k:
                exps[idx, :] = rng.multinomial(k, probs).astype(exps.dtype)
            degrees[idx] = k
            idx += 1
        log_scales = np.empty(m, dtype=np.float64)
        log_m, log_deg, log_n = math.log(max(1, m)), math.log(max(1, d + 1)), math.log(max(1, n))
        for i, k in enumerate(degrees):
            lc = math.lgamma(d + 1) - math.lgamma(int(k) + 1) - math.lgamma(d - int(k) + 1)
            log_scales[i] = 0.5 * (lc + log_deg + int(k) * log_n - log_m)
        coeff = (rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m))) / math.sqrt(2.0)
        if equation_normalize:
            row = np.linalg.norm(coeff, axis=1)
            coeff = coeff / np.where(row > 0, row, 1.0)[:, None] * math.sqrt(float(m))
        obj = cls(n, d, seed, exps, exps.astype(np.float64).T, log_scales, coeff.astype(np.complex128), m, bool(projective_normalize), bool(dynamic_normalize), bool(equation_normalize), max(1, int(eval_block)))
        obj._generation_seconds = now() - t0
        return obj

    @property
    def terms_per_poly(self) -> int:
        return exact_kostlan_terms(self.n, self.d)

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return float(getattr(self, "_generation_seconds", 0.0))

    def _eval_block(self, ZZ: Any) -> Any:
        ZZ = np.asarray(ZZ, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        az = np.abs(ZZ)
        log_amp = np.log(np.maximum(az, 1e-300)) @ self.feature_exps_t + self.feature_log_scales[None, :]
        phase_arg = np.angle(ZZ) @ self.feature_exps_t
        if self.projective_normalize:
            log_amp -= (0.5 * float(self.d) * np.log1p(np.sum(az * az, axis=1)))[:, None]
        if self.dynamic_normalize:
            shift = np.max(np.where(np.isfinite(log_amp), log_amp, -np.inf), axis=1)
            log_amp = log_amp - np.where(np.isfinite(shift), shift, 0.0)[:, None]
            log_amp = np.clip(log_amp, -745.0, 0.0)
        else:
            log_amp = np.clip(log_amp, -745.0, 700.0)
        Phi = np.exp(log_amp) * np.exp(1j * phase_arg)
        out = Phi @ self.coeff.T
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        return out

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = self._eval_block(z)[0]
        self.eval_count += 1
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        chunks = [self._eval_block(ZZ[i:i + self.eval_block]) for i in range(0, int(ZZ.shape[0]), self.eval_block)]
        out = np.vstack(chunks) if chunks else np.empty((0, self.n), dtype=np.complex128)
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms, "lazy_features": int(self.lazy_features)}


@dataclasses.dataclass
class GeometryKernelKostlanSystem:
    n: int
    d: int
    seed: int
    anchors: Any
    anchor_conj_t: Any
    anchor_den: Any
    coeff: Any
    geometry_anchors: int
    anchor_scales: list[float]
    dynamic_normalize: bool = True
    self_normalize: bool = True
    equation_normalize: bool = True
    eval_block: int = 128
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, seed_index: int = 0, equation_normalize: bool = True, geometry_anchors: int = 0, geometry_anchor_cap: int = 4106, geometry_anchor_scales: Optional[Sequence[float]] = None, dynamic_normalize: bool = True, self_normalize: bool = True, eval_block: int = 128) -> "GeometryKernelKostlanSystem":
        t0 = now()
        n, d = int(n), int(d)
        m = max(n + 2, int(geometry_anchors) if int(geometry_anchors) > 0 else auto_geometry_anchor_count(n, d, int(geometry_anchor_cap)))
        seed = stable_seed(n, d, seed_index, salt=0x314C0DE)
        rng = np.random.default_rng(seed)
        scales = [float(x) for x in (geometry_anchor_scales or []) if math.isfinite(float(x)) and float(x) > 0] or [0.25, 0.5, 1.0, 2.0, 4.0]
        anchors = np.zeros((m, n), dtype=np.complex128)
        idx = 1
        axis_scale = scales[min(2, len(scales) - 1)]
        for j in range(n):
            if idx >= m:
                break
            anchors[idx, j] = complex(axis_scale)
            idx += 1
        sqrt_n = math.sqrt(max(1, n))
        while idx < m:
            shell = scales[(idx - 1 - n) % len(scales)]
            v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
            v = v / max(1e-300, float(np.linalg.norm(v))) * (shell * sqrt_n)
            anchors[idx, :] = v
            idx += 1
        coeff = (rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m))) / math.sqrt(2.0 * max(1, m))
        if equation_normalize:
            row = np.linalg.norm(coeff, axis=1)
            coeff = coeff / np.where(row > 0, row, 1.0)[:, None]
        anchor_norm2 = np.sum(np.abs(anchors) ** 2, axis=1)
        obj = cls(n, d, seed, anchors, np.conjugate(anchors).T, np.sqrt(1.0 + anchor_norm2), coeff.astype(np.complex128), m, [float(x) for x in scales], bool(dynamic_normalize), bool(self_normalize), bool(equation_normalize), max(1, int(eval_block)))
        obj._generation_seconds = now() - t0
        return obj

    @property
    def terms_per_poly(self) -> int:
        return exact_kostlan_terms(self.n, self.d)

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return float(getattr(self, "_generation_seconds", 0.0))

    def _kernel_block(self, ZZ: Any) -> Any:
        ZZ = np.asarray(ZZ, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        dot = ZZ @ self.anchor_conj_t
        zn = np.sqrt(1.0 + np.sum(np.abs(ZZ) ** 2, axis=1))
        base = (1.0 + dot) / np.maximum(1e-300, zn[:, None] * self.anchor_den[None, :])
        mag = np.abs(base)
        base = np.where(mag > 1.0, base / np.maximum(mag, 1e-300), base)
        mag = np.minimum(mag, 1.0)
        log_amp = float(self.d) * np.log(np.maximum(mag, 1e-300))
        phase_arg = float(self.d) * np.angle(base)
        if self.dynamic_normalize:
            shift = np.max(np.where(np.isfinite(log_amp), log_amp, -np.inf), axis=1)
            log_amp = log_amp - np.where(np.isfinite(shift), shift, 0.0)[:, None]
        K = np.exp(np.clip(log_amp, -745.0, 0.0)) * np.exp(1j * phase_arg)
        if self.self_normalize:
            row = np.sqrt(np.mean(np.abs(K) ** 2, axis=1))
            K = K / np.maximum(row[:, None], 1e-300)
        return K

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = (self._kernel_block(z) @ self.coeff.T)[0]
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        self.eval_count += 1
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        chunks = []
        for i in range(0, int(ZZ.shape[0]), self.eval_block):
            F = self._kernel_block(ZZ[i:i + self.eval_block]) @ self.coeff.T
            F[~np.isfinite(F)] = complex(1e300, 0.0)
            chunks.append(F)
        out = np.vstack(chunks) if chunks else np.empty((0, self.n), dtype=np.complex128)
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms, "geometry_anchors": int(self.geometry_anchors), "geometry_anchor_scales": list(self.anchor_scales)}


def make_kostlan_base(args: argparse.Namespace, n: int, d: int) -> Any:
    mode = str(getattr(args, "system_mode", "auto")).strip().lower().replace("_", "-")
    dense_terms = exact_kostlan_terms(n, d)
    if mode == "dense" or (mode == "auto" and dense_terms <= int(args.dense_max_terms)):
        return DenseKostlanSystem.make(n, d, int(args.seed_index), bool(args.equation_normalize))
    if mode in {"auto", "geometry", "geometry-kernel", "kernel", "projective-kernel"}:
        return GeometryKernelKostlanSystem.make(n, d, int(args.seed_index), bool(args.equation_normalize), int(args.geometry_anchors), int(args.geometry_anchor_cap), parse_float_list(args.geometry_anchor_scales, [0.25, 0.5, 1.0, 2.0, 4.0], positive=True), bool(args.geometry_dynamic_normalize), bool(args.geometry_self_normalize), int(args.geometry_eval_block))
    if mode not in {"lazy", "lazy-feature", "feature", "stream"}:
        raise ValueError(f"unknown system mode {mode!r}")
    return LazyFeatureKostlanSystem.make(n, d, int(args.seed_index), bool(args.equation_normalize), int(args.lazy_features), int(args.lazy_feature_cap), bool(args.lazy_projective_normalize), bool(args.lazy_dynamic_normalize), int(args.lazy_eval_block))


# ---------------------------------------------------------------------------
# Exact expression polynomial systems (used only as a black-box geometric oracle)
# ---------------------------------------------------------------------------

class SafeExpression:
    _allowed = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub, ast.UAdd, ast.Load)

    def __init__(self, raw: str) -> None:
        self.raw = str(raw).strip()
        tree = ast.parse(self.raw.replace("^", "**"), mode="eval")
        self._validate(tree)
        self.tree = tree

    def _validate(self, node: ast.AST) -> None:
        if not isinstance(node, self._allowed):
            raise ValueError(f"unsupported expression node: {type(node).__name__}")
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
            if not isinstance(node.right, ast.Constant):
                raise ValueError("polynomial powers must be numeric constants")
            if isinstance(node.right.value, complex) or float(node.right.value) != int(node.right.value):
                raise ValueError("polynomial powers must be integer constants")
        for child in ast.iter_child_nodes(node):
            self._validate(child)

    def eval(self, env: dict[str, Any]) -> Any:
        return self._eval(self.tree.body, env)

    def _eval(self, node: ast.AST, env: dict[str, Any]) -> Any:
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            if node.id not in env:
                raise ValueError(f"unknown variable {node.id!r}")
            return env[node.id]
        if isinstance(node, ast.UnaryOp):
            v = self._eval(node.operand, env)
            return -v if isinstance(node.op, ast.USub) else v
        if isinstance(node, ast.BinOp):
            a, b = self._eval(node.left, env), self._eval(node.right, env)
            if isinstance(node.op, ast.Add):
                return a + b
            if isinstance(node.op, ast.Sub):
                return a - b
            if isinstance(node.op, ast.Mult):
                return a * b
            if isinstance(node.op, ast.Div):
                return a / b
            if isinstance(node.op, ast.Pow):
                return a ** int(b)
        raise ValueError(f"unsupported expression node: {type(node).__name__}")


@dataclasses.dataclass
class ExpressionPolynomialSystem:
    n: int
    d: int
    seed: int
    expressions_raw: list[str]
    expressions: list[SafeExpression]
    variable_names: list[str]
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, raw: str, variable_names: Optional[Sequence[str]] = None, seed_index: int = 0) -> "ExpressionPolynomialSystem":
        parts = [p.strip() for p in str(raw).replace("|", ";").split(";") if p.strip()]
        if not parts:
            raise ValueError("--polys/--poly is required for polynomial systems")
        names = [str(p).strip() for p in (variable_names or []) if str(p).strip()] or (["x"] if int(n) == 1 else [f"x{i + 1}" for i in range(int(n))])
        if len(names) != int(n) or len(parts) != int(n):
            raise ValueError(f"expected {n} variables and {n} polynomials")
        return cls(int(n), int(d), stable_seed(n, d, seed_index, salt=0x315EAC7), parts, [SafeExpression(p) for p in parts], names)

    @property
    def terms_per_poly(self) -> int:
        return exact_kostlan_terms(self.n, self.d)

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return 0.0

    def _env(self, ZZ: Any) -> dict[str, Any]:
        env: dict[str, Any] = {"I": 1j, "j": 1j}
        for i, name in enumerate(self.variable_names):
            env[name] = ZZ[:, i]
        if self.n == 1:
            env.setdefault("x", ZZ[:, 0])
            env.setdefault("z", ZZ[:, 0])
        for i in range(self.n):
            env.setdefault(f"x{i + 1}", ZZ[:, i])
            env.setdefault(f"z{i + 1}", ZZ[:, i])
        return env

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = self.eval_batch(np.asarray(z, dtype=np.complex128)[None, :])[0]
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        env = self._env(ZZ)
        cols = []
        for expr in self.expressions:
            val = expr.eval(env)
            arr = np.asarray(val, dtype=np.complex128)
            if arr.ndim == 0:
                arr = np.full(int(ZZ.shape[0]), complex(arr), dtype=np.complex128)
            cols.append(arr.reshape(int(ZZ.shape[0])))
        out = np.stack(cols, axis=1)
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms, "expressions": list(self.expressions_raw), "variables": list(self.variable_names)}


# ---------------------------------------------------------------------------
# The local-jet geometry -- inherited from 317
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class LocalJetGeometry:
    """The geometry that faithfully represents a system: its field of local jets.

    There is NO global precompute (generation is O(1)).  The base system is treated
    as a black-box geometric oracle; around any point the geometry is described by a
    locally sampled jet (value + affine Jacobian + optional curvature) built from a
    small hypercube of ~2n+4 oracle queries.  Per correction step the oracle is hit
    O(n) times -- O(log) of the global comb(n+d, d) monomial budget.

    Because jets are sampled from the TRUE system, eval()/residual() are faithful;
    geometric acceptance therefore coincides with exact acceptance.  Optionally, a
    small jet cache reuses jets across nearby queries.
    """

    base: Any
    n: int
    d: int
    seed: int
    use_quadratic: bool = True
    cache_enabled: bool = True
    cache_decimals: int = 9
    jet_radius: float = 1e-5
    eval_count: int = 0
    slope_count: int = 0
    oracle_samples: int = 0
    jets_built: int = 0
    cache_hits: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, base: Any, n: int, d: int, args: argparse.Namespace) -> "LocalJetGeometry":
        obj = cls(
            base, int(n), int(d), int(getattr(base, "seed", stable_seed(n, d, int(args.seed_index)))),
            bool(args.jet_quadratic), bool(args.jet_cache), max(1, int(args.jet_cache_decimals)),
            max(1e-12, float(args.jet_radius)),
        )
        obj._cache = {}
        obj._build_seconds = 0.0
        return obj

    @property
    def terms_per_poly(self) -> int:
        return int(self.base.terms_per_poly)

    @property
    def total_terms(self) -> int:
        return int(self.base.total_terms)

    @property
    def bezout(self) -> int:
        return int(self.base.bezout)

    @property
    def generation_seconds(self) -> float:
        return float(getattr(self.base, "generation_seconds", 0.0))

    @property
    def samples_per_jet(self) -> int:
        return int(max(2 * self.n + 4, 16))

    # The oracle is queried only through eval/eval_batch -- never via the algebra.
    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = self.base.eval(np.asarray(z, dtype=np.complex128))
        self.eval_count += 1
        self.oracle_samples += 1
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        out = self.base.eval_batch(ZZ)
        self.eval_count += int(ZZ.shape[0])
        self.oracle_samples += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def local_jet(self, z: Sequence[complex], radius: Optional[float] = None) -> dict[str, Any]:
        """Sample the oracle on a small hypercube around z and fit a local jet.

        Returns {f0, J, (Q,) z}, the explicit local geometry at z.  This is the same
        object the hypercube matrix-jet corrector builds; exposing it makes the
        jet-field nature of the geometry concrete and supports a jet-Newton polish.
        """
        zz = np.asarray(z, dtype=np.complex128)
        key = None
        if self.cache_enabled:
            key = tuple(np.round(np.concatenate([zz.real, zz.imag]), self.cache_decimals))
            hit = self._cache.get(key)
            if hit is not None:
                self.cache_hits += 1
                return hit
        h = float(radius) if radius is not None else self.jet_radius * max(1.0, float(np.linalg.norm(zz)))
        f0 = self.eval(zz)
        eye = np.eye(self.n, dtype=np.complex128)
        Pp = zz[None, :] + h * eye
        Pm = zz[None, :] - h * eye
        Fp = self.eval_batch(Pp)
        Fm = self.eval_batch(Pm)
        J = ((Fp - Fm) / (2.0 * h)).T
        jet: dict[str, Any] = {"z": zz, "f0": f0, "J": J, "h": float(h)}
        if self.use_quadratic:
            jet["Qdiag"] = ((Fp - 2.0 * f0[None, :] + Fm) / (h * h)).T
        self.jets_built += 1
        if key is not None:
            self._cache[key] = jet
        return jet

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        jet = self.local_jet(np.asarray(a, dtype=np.complex128))
        self.slope_count += 1
        return jet["J"]

    def stats(self) -> dict[str, Any]:
        return {
            "eval_count": int(self.eval_count), "slope_count": int(self.slope_count),
            "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope),
            "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms,
            "geometry": "local-jet-field", "oracle_samples": int(self.oracle_samples),
            "jets_built": int(self.jets_built), "jet_cache_hits": int(self.cache_hits),
            "samples_per_jet": int(self.samples_per_jet), "use_quadratic": bool(self.use_quadratic),
            "base_stats": self.base.stats(),
        }


def describe_base_backend(base: Any) -> str:
    if isinstance(base, ExpressionPolynomialSystem):
        return "polynomial-exact"
    if isinstance(base, DenseKostlanSystem):
        return "ks-dense"
    if isinstance(base, LazyFeatureKostlanSystem):
        return "ks-lazy-feature"
    if isinstance(base, GeometryKernelKostlanSystem):
        return "ks-geometry-kernel"
    return "unknown"


def make_system_410(args: argparse.Namespace, n: int, d: int) -> tuple[LocalJetGeometry, Any, str, str]:
    source = str(args.system_source).strip().lower().replace("_", "-")
    if source in {"poly", "polynomial", "expr", "expression"}:
        vars_raw = str(args.variables or "").strip()
        variables = [p.strip() for p in vars_raw.replace(";", ",").split(",") if p.strip()] if vars_raw else None
        base: Any = ExpressionPolynomialSystem.make(n, d, str(args.polys or ""), variables, int(args.seed_index))
        source = "polynomial"
    elif source in {"ks", "kostlan"}:
        base = make_kostlan_base(args, n, d)
        source = "ks"
    else:
        raise ValueError(f"unknown --system-source {source!r}")
    geometry = LocalJetGeometry.make(base, n, d, args)
    return geometry, base, source, describe_base_backend(base)


# ---------------------------------------------------------------------------
# Charts, starts, and correctors (operate on the geometry through a generic target)
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class LinearChart:
    A: Any
    Ainv: Any

    @classmethod
    def identity(cls, n: int, scale: float = 1.0) -> "LinearChart":
        A = np.eye(int(n), dtype=np.complex128) * complex(scale)
        Ainv = np.eye(int(n), dtype=np.complex128) / complex(scale)
        return cls(A, Ainv)

    def z_from_y(self, y: Sequence[complex]) -> Any:
        return self.A @ np.asarray(y, dtype=np.complex128)

    def y_from_z(self, z: Sequence[complex]) -> Any:
        return self.Ainv @ np.asarray(z, dtype=np.complex128)


@dataclasses.dataclass
class TargetTrack:
    system: Any
    chart: LinearChart

    def eval(self, y: Sequence[complex]) -> Any:
        return self.system.eval(self.chart.z_from_y(y))

    def eval_batch(self, Y: Any) -> Any:
        YY = np.asarray(Y, dtype=np.complex128)
        if YY.ndim == 1:
            return self.eval(YY)[None, :]
        return self.system.eval_batch(YY @ np.asarray(self.chart.A, dtype=np.complex128).T)

    def residual(self, y: Sequence[complex]) -> float:
        return finite_norm(self.eval(y))

    def residuals_batch(self, Y: Any) -> Any:
        try:
            return safe_norms(self.eval_batch(Y))
        except Exception:
            YY = np.asarray(Y, dtype=np.complex128)
            return np.full(1 if YY.ndim == 1 else int(YY.shape[0]), np.inf, dtype=float)


DEFAULT_POWERS = [2.0 ** k for k in range(-20, 21)] + [3.0 * (2.0 ** k) for k in range(-12, 18)] + [10.0 ** k for k in range(-6, 7)]
DEFAULT_ANGLES_DEG = [0, 6, 12, 18, 24, 32, 40, 48, 56, 64, 72, 80, 86, 89, 90, 91, 94, 100, 108, 116, 128, 140, 152, 164, 172]


def monomial_scale_ladder(subdivisions: int, octaves: int, base: float = 2.0) -> list[float]:
    """The 'complete logarithmic ladder': equally spaced logarithmic scale offsets.

    Between consecutive powers ``base**m`` the ladder inserts the ``subdivisions - 1`` geometric
    (proportional) means ``base**(m + k/subdivisions)``, k = 1..subdivisions-1.  For base 2 and
    subdivisions 3 those means are the cube-root proportional means 2**(1/3), 2**(2/3) -- exactly
    the x**(1/p), x**(2/p), ... ladder of the Pandrosion construction (the diagram that yields
    cbrt(2) yields cbrt(4) = cbrt(2)**2 in the same figure).

    The monomial scale-palette theorem proves that equally spaced logarithmic offsets minimise the
    worst-case raw residual multiplier 1 - p/S_p(q**(1/(2K))).  Using the ladder as the multi-start
    homothety palette (and as the IRP chart palette) therefore makes |log y| uniformly small before
    iteration, so a start is more likely to land inside a convergence basin.
    """
    p = max(1, int(subdivisions))
    m = max(1, int(octaves))
    b = float(base) if float(base) > 1.0 else 2.0
    return [b ** (j / p) for j in range(-m * p, m * p + 1)]
DEFAULT_RADII = [0.025, 0.04, 0.06, 0.08, 0.12, 0.18, 0.27, 0.4, 0.6, 0.85, 1.15, 1.55, 2.05, 2.75, 3.6, 4.8, 6.4]
DEFAULT_GAINS = [0.035, 0.055, 0.085, 0.12, 0.18, 0.27, 0.4, 0.58, 0.78, 1.0, 1.28, 1.65, 2.2, 3.0, 4.2, 6.0, 8.5, 12.0]


def raw_direction(n: int, trial: int, seed: int) -> Any:
    vals = []
    for j in range(int(n)):
        h1 = splitmix64(seed + 0xD1A5E + 1000003 * trial + 4109 * (j + 1))
        h2 = splitmix64(seed + 0xBADC0DE + 1000033 * trial + 9176 * (j + 1))
        vals.append(math.exp(0.45 * (2.0 * u01(h2) - 1.0)) * phase(2.0 * math.pi * u01(h1)))
    v = np.asarray(vals, dtype=np.complex128)
    return v / max(1e-300, float(np.linalg.norm(v))) * math.sqrt(max(1, int(n)))


def universal_atlas_start(target: TargetTrack, n: int, trial: int, seed: int, powers: Sequence[float], angles: Sequence[float], radii: Sequence[float], cap: float, roots_found: int, duplicates: int, failures: int, target_count: int, universal_cells: int = 16, universal_shells: int = 5, atlas_selection: str = "diverse-shell", **_: Any) -> tuple[Any, dict[str, Any]]:
    cells = max(1, int(universal_cells))
    shells = max(1, int(universal_shells))
    pows = [min(max(float(x), 1e-300), float(cap)) for x in powers if float(x) > 0] or [1.0]
    mode = str(atlas_selection or "diverse-shell").strip().lower().replace("_", "-")
    if mode in {"diverse", "diverse-shell", "shell", "raw-shell"}:
        rr = [float(x) for x in radii if math.isfinite(float(x)) and float(x) > 0] or DEFAULT_RADII
        idx = int(trial)
        rad = rr[(idx + 3 * roots_found + failures) % len(rr)]
        y = float(rad) * raw_direction(n, idx, seed + 65537 * idx + 104729 * duplicates)
        r0 = target.residual(y)
        meta = {
            "chart": "410-standalone-universal-diverse-shell-atlas",
            "atlas_mode": "410-diverse-shell",
            "atlas_selection": mode,
            "atlas_startopt_bypass_recommended": True,
            "atlas_cells_tested": 1,
            "atlas_admissible_cells": int(math.isfinite(r0)),
            "atlas_cell_residual": float(r0),
            "atlas_selected_index": int(idx),
            "atlas_selected_layer": int(idx % shells),
            "homothety": float(rad),
            "base_homothety": 1.0,
            "thales_thrust": 1.0,
            "theta_deg": None,
            "base_radius": float(rad),
            "dup_pressure": float((duplicates + 1.0) / (roots_found + 1.0)),
            "fail_pressure": float((failures + 1.0) / (trial + 1.0)),
            "progress": float(min(1.0, roots_found / max(1.0, float(target_count)))),
        }
        return np.asarray(y, dtype=np.complex128).copy(), meta

    candidates = []
    metas = []
    for k in range(cells):
        idx = trial * cells + k
        layer = (idx // cells) % shells
        pwr = pows[(idx * 37 + 11 * layer + roots_found) % len(pows)]
        rad = radii[(idx * 13 + failures) % len(radii)] if radii else 1.0
        theta = angles[(idx * 19 + duplicates) % len(angles)] if angles else 0.0
        thrust = [1.0, 1.6, 2.5, 4.0, 6.5, 10.0, 16.0][idx % 7]
        lam = min(float(cap), pwr * (thrust ** (0.2 + 0.8 * min(1.0, roots_found / max(1.0, float(target_count))))))
        d = raw_direction(n, idx, seed + 65537 * k)
        y = rad * d
        c, s = math.cos(theta), math.sin(theta)
        pole = phase(0.37 * (idx + 1))
        w = y / pole
        den = -s * w + c
        den = np.where(np.abs(den) < 1e-12, den + 1e-12, den)
        y = lam * pole * ((c * w + s) / den)
        candidates.append(np.asarray(y, dtype=np.complex128))
        metas.append({"homothety": float(lam), "base_homothety": float(pwr), "thales_thrust": float(thrust), "theta_deg": float(math.degrees(theta)), "base_radius": float(rad), "atlas_selected_index": int(idx), "atlas_selected_layer": int(layer)})
    Y = np.asarray(candidates, dtype=np.complex128)
    R = target.residuals_batch(Y)
    idx_best = int(np.nanargmin(R)) if np.any(np.isfinite(R)) else 0
    meta = dict(metas[idx_best])
    meta.update({"chart": "410-standalone-universal-mobius-kernel-atlas", "atlas_mode": "410-compact-batched-score", "atlas_selection": mode, "atlas_startopt_bypass_recommended": False, "atlas_cells_tested": int(cells), "atlas_admissible_cells": int(np.sum(np.isfinite(R))), "atlas_cell_residual": float(R[idx_best]) if idx_best < len(R) else float("inf"), "dup_pressure": float((duplicates + 1.0) / (roots_found + 1.0)), "fail_pressure": float((failures + 1.0) / (trial + 1.0)), "progress": float(min(1.0, roots_found / max(1.0, float(target_count))))})
    return Y[idx_best].copy(), meta


def origin_affine_start(target: TargetTrack, n: int, h: float, max_norm: float, matrix_rcond: float = 1e-12, matrix_condition_cap: float = 1e12) -> tuple[Optional[Any], dict[str, Any]]:
    h0 = max(1e-12, float(h))
    y0 = np.zeros(int(n), dtype=np.complex128)
    try:
        f0 = target.eval(y0)
        P = np.eye(int(n), dtype=np.complex128) * h0
        J = ((target.eval_batch(P) - f0[None, :]) / h0).T
        y, mmeta = matrix_svd_step(J, -f0, float(matrix_rcond), float(matrix_condition_cap), 0.0)
        if not np.all(np.isfinite(y)):
            return None, {"origin_seed_enabled": True, "origin_seed_status": "nonfinite"}
        cap = float(max_norm)
        yn = float(np.linalg.norm(y))
        if math.isfinite(cap) and cap > 0 and yn > cap:
            y = y * (cap / max(yn, 1e-300))
        return y, {
            "origin_seed_enabled": True,
            "origin_seed_status": "ok",
            "origin_seed_method": "matrix-svd-filter",
            "origin_seed_h": float(h0),
            "origin_seed_r0": finite_norm(f0),
            "origin_seed_r1": target.residual(y),
            "origin_seed_norm": float(np.linalg.norm(y)),
            **{f"origin_{k}": v for k, v in mmeta.items()},
        }
    except Exception as exc:
        return None, {"origin_seed_enabled": True, "origin_seed_status": f"error:{type(exc).__name__}", "origin_seed_h": float(h0)}


def startopt(target: TargetTrack, y0: Any, trial: int, seed: int, steps: int, candidates: int, gains: Sequence[float], micro_epochs: int) -> tuple[Any, dict[str, Any]]:
    best = np.asarray(y0, dtype=np.complex128).copy()
    best_r = target.residual(best)
    initial = best_r
    evals = 1
    chosen_gain = 1.0
    for step in range(max(0, int(steps))):
        pool = [(1.0, best)]
        for c in range(max(0, int(candidates) - 1)):
            gain = float(gains[(trial + 3 * step + c) % len(gains)])
            cand = gain * best
            if c % 3 == 1:
                noise = raw_direction(len(best), trial + 31 * step + c, seed)
                cand = cand + 0.08 * max(1.0, float(np.linalg.norm(best))) * noise
            elif c % 3 == 2:
                cand = cand * phase(0.11 * (c + 1))
            pool.append((gain, np.asarray(cand, dtype=np.complex128)))
        Y = np.asarray([p[1] for p in pool], dtype=np.complex128)
        R = target.residuals_batch(Y)
        evals += len(pool)
        idx = int(np.nanargmin(R)) if np.any(np.isfinite(R)) else 0
        if float(R[idx]) < best_r:
            best, best_r, chosen_gain = Y[idx].copy(), float(R[idx]), float(pool[idx][0])
    return best, {"startopt_enabled": bool(steps > 0), "startopt_r0": float(initial), "startopt_r1": float(best_r), "startopt_ratio": (float(best_r / initial) if math.isfinite(best_r) and math.isfinite(initial) and initial > 0 else None), "startopt_steps": int(max(0, steps)), "startopt_evals": int(evals), "startopt_micro_epochs": int(max(0, micro_epochs)), "startopt_gain": float(chosen_gain), "startopt_batch_numpy": True}


def _line_lambdas(line_search: int, line_grid: Sequence[float]) -> list[float]:
    vals = [float(x) for x in line_grid if math.isfinite(float(x)) and float(x) > 0]
    if vals:
        return vals
    return [1.0, 0.75, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03125, 0.02][:max(1, int(line_search))]


def directional_sketch_delta(target: TargetTrack, y: Any, f: Any, ep: int, seed: int, lm_damping: float = 0.0, trust_radius: float = 0.0, matrix_rcond: float = 1e-12, matrix_condition_cap: float = 1e12) -> Optional[tuple[Any, dict[str, Any]]]:
    """Estimate only J D from coded directional oracle samples and return delta=D v."""
    if not bool(DIRECTIONAL_JET):
        return None
    n = int(len(y))
    if n < int(DIRECTIONAL_JET_MIN_N) or selected_matrix_algorithm_410(n) != "sketch-ns":
        return None
    base_k = int(SKETCH_DIM) if int(SKETCH_DIM) > 0 else int(math.ceil(float(DIRECTIONAL_JET_FACTOR) * math.sqrt(float(max(1, n)))))
    k = max(1, min(n, base_k))
    if k >= n:
        return None
    yy = np.asarray(y, dtype=np.complex128)
    ff = np.asarray(f, dtype=np.complex128)
    yn = max(1.0, float(np.linalg.norm(yy)))
    h = 1e-5 * yn
    salt_ep = 0 if bool(DIRECTIONAL_BASIS_REUSE) else int(ep)
    salt = int(seed) + 65537 * int(salt_ep) + 0x411D1A
    P = _sketch_basis_np(n, k, salt=salt)
    device, dtype = torch_device(), torch_complex_dtype()
    rhs = torch.as_tensor(-ff, dtype=dtype, device=device)

    def estimate_and_solve(basis: Any, mode: str) -> tuple[Any, Any, dict[str, Any], int, Any]:
        dirs = np.asarray(basis, dtype=np.complex128).T
        qloc = int(dirs.shape[0])
        if mode == "central":
            Fp = target.eval_batch(yy[None, :] + h * dirs)
            Fm = target.eval_batch(yy[None, :] - h * dirs)
            JD = ((Fp - Fm) / (2.0 * h)).T
            samples = int(2 * qloc)
        else:
            Fp = target.eval_batch(yy[None, :] + h * dirs)
            JD = ((Fp - ff[None, :]) / h).T
            samples = int(qloc)
        if not np.all(np.isfinite(JD)):
            raise ValueError("nonfinite_directional_jp")
        fast_fallback: Optional[str] = None
        if bool(DIRECTIONAL_FAST_PROJECTOR):
            try:
                u_fast, fast_meta = _numpy_fast_projected_apply(JD, -ff, float(lm_damping), float(matrix_condition_cap))
                lin_fast = JD @ u_fast + ff
                fast_rel = finite_norm(lin_fast) / max(1e-300, finite_norm(ff))
                fast_meta["reduced_fast_projector_relative_residual"] = float(fast_rel)
                if fast_rel <= float(DIRECTIONAL_FAST_PROJECTOR_CAP):
                    fast_meta["reduced_fast_projector_accepted"] = True
                    return u_fast, lin_fast, fast_meta, samples, JD
                fast_fallback = f"projected-linear-residual:{fast_rel:.3e}"
            except Exception as exc:
                fast_fallback = f"projected-error:{type(exc).__name__}"
        AP = torch.as_tensor(JD, dtype=dtype, device=device)
        u, rmeta = _torch_reduced_apply(AP, rhs, float(matrix_rcond), float(matrix_condition_cap), float(lm_damping))
        if fast_fallback is not None:
            rmeta["reduced_fast_projector"] = False
            rmeta["reduced_fast_projector_attempted"] = True
            rmeta["reduced_fast_projector_fallback"] = fast_fallback
            rmeta["reduced_fast_projector_cap"] = float(DIRECTIONAL_FAST_PROJECTOR_CAP)
        u_np = u.detach().cpu().numpy().astype(np.complex128, copy=False)
        lin = JD @ u_np + ff
        return u_np, lin, rmeta, samples, JD

    def solve_or_raise(basis: Any, mode: str) -> tuple[Any, Any, dict[str, Any], int, Any]:
        u_loc, lin_loc, meta_loc, samples_loc, jd_loc = estimate_and_solve(basis, mode)
        if not np.all(np.isfinite(u_loc)):
            raise ValueError("nonfinite_directional_coefficients")
        return u_loc, lin_loc, meta_loc, samples_loc, jd_loc

    full_diff = "central" if DIRECTIONAL_DIFF_MODE == "central" else "forward"
    diff_used = full_diff
    probe_kind = "full-sketch"
    active_basis = P
    active_dim = k
    coded_q = directional_coded_dim_for(n, k)
    fallback_reason: Optional[str] = None
    coded_fallback_reason: Optional[str] = None
    u_np: Any
    lin: Any
    rmeta: dict[str, Any]
    samples: int
    JD: Any

    accepted = False
    if bool(DIRECTIONAL_CODED_PROBE) and coded_q < k:
        R = _fast_projected_coded_basis_np(k, coded_q, salt=salt + 0xC0DEC0DE)
        C = P @ R
        try:
            u_np, lin, rmeta, samples, JD = solve_or_raise(C, full_diff)
            lin_rel = finite_norm(lin) / max(1e-300, finite_norm(ff))
            accepted = DIRECTIONAL_DIFF_MODE != "auto" or lin_rel <= float(DIRECTIONAL_AUTO_CENTRAL_CAP)
            if accepted:
                active_basis, active_dim, probe_kind, diff_used = C, coded_q, "coded-probe", full_diff
            elif DIRECTIONAL_DIFF_MODE == "auto":
                coded_fallback_reason = f"coded-forward-linear-residual:{lin_rel:.3e}"
                try:
                    u_np, lin, rmeta, samples, JD = solve_or_raise(C, "central")
                    lin_rel = finite_norm(lin) / max(1e-300, finite_norm(ff))
                    accepted = lin_rel <= float(DIRECTIONAL_AUTO_CENTRAL_CAP)
                    if accepted:
                        active_basis, active_dim, probe_kind, diff_used = C, coded_q, "coded-probe", "central"
                        fallback_reason = coded_fallback_reason
                    else:
                        coded_fallback_reason = f"coded-central-linear-residual:{lin_rel:.3e}"
                except Exception as exc:
                    coded_fallback_reason = f"coded-central-error:{type(exc).__name__}"
        except Exception as exc:
            if DIRECTIONAL_DIFF_MODE != "auto":
                return np.zeros_like(yy), {"hypercube_error": f"directional_coded_failed:{type(exc).__name__}", "hypercube_order": 0, "hypercube_nodes": int(coded_q if full_diff == "forward" else 2 * coded_q), "directional_jet": True, "directional_probe_kind": "coded-probe", "directional_parent_sketch_dim": int(k), "directional_sketch_dim": int(coded_q), "directional_diff": full_diff}
            coded_fallback_reason = f"coded-forward-error:{type(exc).__name__}"

    if not accepted:
        try:
            u_np, lin, rmeta, samples, JD = solve_or_raise(P, full_diff)
            lin_rel = finite_norm(lin) / max(1e-300, finite_norm(ff))
            if DIRECTIONAL_DIFF_MODE == "auto" and lin_rel > float(DIRECTIONAL_AUTO_CENTRAL_CAP):
                fallback_reason = coded_fallback_reason or f"forward-linear-residual:{lin_rel:.3e}"
                diff_used = "central"
                u_np, lin, rmeta, samples, JD = solve_or_raise(P, "central")
            else:
                fallback_reason = coded_fallback_reason
                diff_used = full_diff
            active_basis, active_dim, probe_kind = P, k, "full-sketch"
        except Exception as exc:
            if DIRECTIONAL_DIFF_MODE == "auto" and full_diff == "forward":
                fallback_reason = coded_fallback_reason or f"forward-error:{type(exc).__name__}"
                try:
                    diff_used = "central"
                    u_np, lin, rmeta, samples, JD = solve_or_raise(P, "central")
                    active_basis, active_dim, probe_kind = P, k, "full-sketch"
                except Exception as exc2:
                    return np.zeros_like(yy), {"hypercube_error": f"directional_reduced_failed:{type(exc2).__name__}", "hypercube_order": 0, "hypercube_nodes": int(2 * k), "directional_jet": True, "directional_probe_kind": "full-sketch", "directional_parent_sketch_dim": int(k), "directional_sketch_dim": int(k), "directional_diff": "central", "directional_forward_fallback": fallback_reason, "directional_coded_fallback": coded_fallback_reason}
            else:
                return np.zeros_like(yy), {"hypercube_error": f"directional_reduced_failed:{type(exc).__name__}", "hypercube_order": 0, "hypercube_nodes": int(k if full_diff == "forward" else 2 * k), "directional_jet": True, "directional_probe_kind": "full-sketch", "directional_parent_sketch_dim": int(k), "directional_sketch_dim": int(k), "directional_diff": full_diff, "directional_coded_fallback": coded_fallback_reason}

    delta = active_basis @ u_np
    if float(trust_radius) > 0.0:
        lim = float(trust_radius) * yn
        dn = float(np.linalg.norm(delta))
        if math.isfinite(dn) and dn > lim > 0:
            delta = delta * (lim / max(dn, 1e-300))
    lin = JD @ u_np + ff
    return delta, {
        "hypercube_order": 1,
        "hypercube_nodes": int(samples),
        "hypercube_delta1_norm": finite_norm(delta),
        "directional_jet": True,
        "directional_jet_algorithm": "411-cached-coded-probe-or-full-directional-oracle-sampled-jp-subspace",
        "directional_diff": str(diff_used),
        "directional_auto_central_cap": float(DIRECTIONAL_AUTO_CENTRAL_CAP),
        "directional_forward_fallback": fallback_reason,
        "directional_coded_fallback": coded_fallback_reason,
        "directional_probe_kind": str(probe_kind),
        "directional_fast_projected_coded": bool(probe_kind == "coded-probe"),
        "directional_coded_requested": bool(DIRECTIONAL_CODED_PROBE),
        "directional_parent_sketch_dim": int(k),
        "directional_coded_dim": int(coded_q),
        "directional_fast_projector": bool(DIRECTIONAL_FAST_PROJECTOR),
        "directional_fast_projector_cap": float(DIRECTIONAL_FAST_PROJECTOR_CAP),
        "directional_sketch_dim": int(active_dim),
        "directional_sketch_ratio": float(active_dim / max(1, n)),
        "directional_oracle_samples": int(samples),
        "directional_full_hypercube_samples_avoided": int(max(2 * n + 4, 16) - int(samples)),
        "directional_linear_residual": finite_norm(lin),
        "directional_linear_relative_residual": finite_norm(lin) / max(1e-300, finite_norm(ff)),
        "directional_h": float(h),
        "directional_sketch_mode": str(SKETCH_MODE),
        "directional_sketch_solver": str(SKETCH_SOLVER),
        "directional_basis_reuse": bool(DIRECTIONAL_BASIS_REUSE),
        **{f"directional_{kk}": vv for kk, vv in basis_cache_stats().items()},
        "directional_lm_damping": float(lm_damping),
        "directional_trust_radius": float(trust_radius),
        **{f"directional_{kk}": vv for kk, vv in rmeta.items()},
    }


def hypercube_delta(target: TargetTrack, y: Any, f: Any, ep: int, seed: int, cloud_nodes: int = 0, lm_damping: float = 0.0, trust_radius: float = 0.0, matrix_rcond: float = 1e-12, matrix_condition_cap: float = 1e12) -> tuple[Any, dict[str, Any]]:
    direct = directional_sketch_delta(target, y, f, ep, seed, lm_damping, trust_radius, matrix_rcond, matrix_condition_cap)
    if direct is not None:
        return direct
    n = len(y)
    M = max(int(cloud_nodes), 2 * n + 4, 16)
    yn = max(1.0, float(np.linalg.norm(y)))
    h = 1e-5 * yn
    rng = np.random.default_rng(int(seed) + int(ep) * 1337)
    signs = rng.choice([-1.0, 1.0], size=(M, n))
    dY = h * signs
    dF = target.eval_batch(y[None, :] + dY) - f[None, :]
    B, fit_meta = matrix_least_squares_fit(dY, dF, float(matrix_rcond), float(matrix_condition_cap))
    A = B.T
    if not np.all(np.isfinite(A)):
        return np.zeros_like(y), {"hypercube_error": "nonfinite_jacobian", "hypercube_order": 0, "hypercube_nodes": M}
    # Matrix-jet inversion of A.  A is not solved as an exact square system; every
    # right-hand side is filtered through the same regularised SVD pseudo-inverse.
    tr = float(trust_radius) > 0.0

    lim = float(trust_radius) * yn if tr else float("inf")

    def cap(v: Any) -> Any:
        if tr:
            nv = float(np.linalg.norm(v))
            if math.isfinite(nv) and nv > lim > 0.0:
                return v * (lim / nv)
        return v

    extra = {
        "hypercube_lm_damping": float(lm_damping),
        "hypercube_trust_radius": float(trust_radius),
        "hypercube_matrix_rcond": float(matrix_rcond),
        "hypercube_matrix_condition_cap": float(matrix_condition_cap),
        **{f"hypercube_fit_{k}": v for k, v in fit_meta.items()},
    }
    try:
        d1_raw, step_meta = matrix_svd_step(A, -f, float(matrix_rcond), float(matrix_condition_cap), float(lm_damping))
        d1 = cap(d1_raw)
    except Exception:
        return np.zeros_like(y), {"hypercube_error": "matrix_filter_failed", "hypercube_order": 0, "hypercube_nodes": M, **extra}
    extra.update({f"hypercube_{k}": v for k, v in step_meta.items()})
    if int(step_meta.get("matrix_rank", 0)) <= 0:
        return np.zeros_like(y), {"hypercube_error": "zero_rank_matrix", "hypercube_order": 0, "hypercube_nodes": M, **extra}
    nd1 = float(np.linalg.norm(d1))
    if not math.isfinite(nd1) or nd1 < 1e-14 or (not tr and nd1 > 1e4 * yn):
        return d1, {"hypercube_order": 1, "hypercube_delta1_norm": nd1, "hypercube_nodes": M, **extra}
    hh = max(1e-5, min(1e-2, 0.05 * yn / nd1))
    p1, m1 = target.eval(y + hh * d1), target.eval(y - hh * d1)
    p2, m2 = target.eval(y + 2.0 * hh * d1), target.eval(y - 2.0 * hh * d1)
    D2 = (-p2 + 16.0 * p1 - 30.0 * f + 16.0 * m1 - m2) / (12.0 * hh**2)
    D3 = (p2 - 2.0 * p1 + 2.0 * m1 - m2) / (2.0 * hh**3)
    try:
        d2 = cap(matrix_svd_step(A, -0.5 * D2, float(matrix_rcond), float(matrix_condition_cap), float(lm_damping))[0])
    except Exception:
        return d1, {"hypercube_order": 1, "hypercube_error": "matrix_delta2_failed", "hypercube_nodes": M, "hypercube_delta1_norm": nd1, **extra}
    nd2 = float(np.linalg.norm(d2))
    if not math.isfinite(nd2) or (not tr and nd2 > 2.0 * max(nd1, yn)):
        return d1, {"hypercube_order": 1, "hypercube_error": "delta2_rejected", "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2, **extra}
    try:
        f2 = target.eval(y + hh * d2)
        f12 = target.eval(y + hh * d1 + hh * d2)
        cross = (f12 - p1 - f2 + f) / (hh**2)
        d3 = cap(matrix_svd_step(A, -(cross + (1.0 / 6.0) * D3), float(matrix_rcond), float(matrix_condition_cap), float(lm_damping))[0])
    except Exception:
        return d1 + d2, {"hypercube_order": 2, "hypercube_error": "matrix_delta3_failed", "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2, **extra}
    nd3 = float(np.linalg.norm(d3))
    if not math.isfinite(nd3) or (not tr and nd3 > 2.0 * max(nd1, yn)):
        return d1 + d2, {"hypercube_order": 2, "hypercube_error": "delta3_rejected", "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2, "hypercube_delta3_norm": nd3, **extra}
    return cap(d1 + d2 + d3), {"hypercube_order": 4, "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2, "hypercube_delta3_norm": nd3, **extra}


def hypercube_matrixjet_corrector(target: TargetTrack, y0: Sequence[complex], max_epochs: int, tol: float, accept: float, trial_timeout: float, line_search: int = 12, line_grid: Sequence[float] = (), direction_seed: int = 0, cloud_nodes: int = 0, lm_damping: float = 0.0, trust_radius: float = 0.0, matrix_rcond: float = 1e-12, matrix_condition_cap: float = 1e12) -> dict[str, Any]:
    t0 = now()
    deadline = t0 + float(trial_timeout) if trial_timeout and trial_timeout > 0 else None
    y = np.asarray(y0, dtype=np.complex128).copy()
    best_y, best_r = y.copy(), target.residual(y)
    status, ok, epochs = "started", False, 0
    total_line, total_hyper, used, last = 0, 0, 0, {}
    lambdas = _line_lambdas(line_search, line_grid)
    for ep in range(max(1, int(max_epochs))):
        if deadline is not None and now() > deadline:
            status = "timeout"
            break
        f = target.eval(y)
        r = finite_norm(f)
        if r < best_r:
            best_y, best_r = y.copy(), r
        if r <= max(float(tol), float(accept)) and (accept <= 0 or r < accept):
            ok, status = True, "converged"
            break
        delta, meta = hypercube_delta(target, y, f, ep, int(direction_seed), int(cloud_nodes), float(lm_damping), float(trust_radius), float(matrix_rcond), float(matrix_condition_cap))
        last = dict(meta)
        total_hyper += int(meta.get("hypercube_nodes", max(2 * len(y) + 4, 16))) + 7
        used += 1
        if not np.all(np.isfinite(delta)):
            status = "nonfinite-matrix-step"
            break
        L = np.asarray(lambdas, dtype=float)
        Y = y[None, :] + L[:, None] * delta[None, :]
        R = target.residuals_batch(Y)
        total_line += int(len(Y))
        better = np.isfinite(R) & ((R < r) | (R < best_r))
        if not np.any(better):
            status = "no-matrixjet-decrease"
            break
        idx = int(np.nanargmin(np.where(better, R, np.inf)))
        y = Y[idx].copy()
        if float(R[idx]) < best_r:
            best_y, best_r = y.copy(), float(R[idx])
        epochs = ep + 1
        last["hypercube_chosen_lambda"] = float(L[idx])
    else:
        status = "max-epochs"
    final_r = target.residual(best_y)
    if final_r <= max(float(tol), float(accept)) and (accept <= 0 or final_r < accept):
        ok, status, best_r = True, "converged", final_r
    return {"accepted": bool(ok if accept <= 0 else (math.isfinite(best_r) and best_r < accept)), "ok": bool(ok), "status": status, "epochs": int(epochs), "residual": float(best_r), "y": best_y, "seconds": float(now() - t0), "slope_cond": None, "corrector": "411-cached-hypercube-pandrosion-irp-gemm-matrixjet", "line_search_evals": int(total_line), "line_lambdas": [float(x) for x in lambdas], "hypercube_total_evals": int(total_hyper), "hypercube_used_count": int(used), "matrix_rcond": float(matrix_rcond), "matrix_condition_cap": float(matrix_condition_cap), "halley_enabled": False, "halley_total_evals": 0, "halley_used_count": 0, **last}


class NonlinearChartTarget:
    def __init__(self, base: TargetTrack, scale: complex, reciprocal: bool) -> None:
        self.base, self.scale, self.reciprocal = base, complex(scale), bool(reciprocal)

    def to_base(self, u: Any) -> Any:
        uu = np.asarray(u, dtype=np.complex128)
        if not self.reciprocal:
            return self.scale * uu
        den = self.scale * uu
        den = np.where(np.abs(den) < 1e-14, den + 1e-14, den)
        return 1.0 / den

    def from_base(self, y: Any) -> Any:
        yy = np.asarray(y, dtype=np.complex128)
        if not self.reciprocal:
            return yy / self.scale
        den = self.scale * yy
        den = np.where(np.abs(den) < 1e-14, den + 1e-14, den)
        return 1.0 / den

    def eval(self, u: Sequence[complex]) -> Any:
        return self.base.eval(self.to_base(u))

    def eval_batch(self, U: Any) -> Any:
        UU = np.asarray(U, dtype=np.complex128)
        return self.base.eval_batch(self.to_base(UU))

    def residual(self, u: Sequence[complex]) -> float:
        return finite_norm(self.eval(u))

    def residuals_batch(self, U: Any) -> Any:
        return safe_norms(self.eval_batch(U))


def complex_scale_palette(gains: Sequence[float], phases: Sequence[float], top: int) -> list[complex]:
    out: list[complex] = []
    for g in gains:
        for ph in phases:
            out.append(float(g) * phase(float(ph)))
            if len(out) >= int(top):
                return out
    return out or [1.0 + 0.0j]


def chart_candidates(base: TargetTrack, y: Any, scales: Sequence[complex], inversion: bool, top: int) -> list[tuple[float, NonlinearChartTarget, Any, dict[str, Any]]]:
    recs = []
    for s in list(scales)[:max(1, int(top))]:
        for recip in ([False, True] if inversion else [False]):
            ct = NonlinearChartTarget(base, complex(s), recip)
            u0 = ct.from_base(y)
            r = ct.residual(u0)
            e = log_energy(u0)
            if math.isfinite(r) and math.isfinite(e):
                recs.append((r * (1.0 + 0.002 * e), ct, u0, {"kind": "reciprocal" if recip else "direct", "scale": s, "score": r, "log_energy": e}))
    recs.sort(key=lambda x: x[0])
    return recs


def lazy_irp_hypercube_torch_gemm_matrixjet_corrector_410(base: TargetTrack, y0: Sequence[complex], max_epochs: int, tol: float, accept: float, trial_timeout: float, line_search: int = 12, line_grid: Sequence[float] = (), direction_seed: int = 0, cloud_nodes: int = 0, irp_layers: int = 2, irp_inner_epochs: int = 2, irp_scales: Optional[Sequence[complex]] = None, irp_chart_top: int = 2, irp_inversion: bool = True, collapse: bool = True, collapse_residual: float = 1e-4, collapse_drop: float = 0.42, collapse_rel_step: float = 0.35, collapse_after: int = 2, local_inner_epochs: int = 3, lazy_direct_epochs: int = 1, lazy_trigger_drop: float = 0.82, lazy_trigger_after: int = 1, lazy_bad_cond: float = 1e10, lazy_log_energy: float = 8.0, eager_irp: bool = False, rescue_collapsed: bool = False, lm_damping: float = 0.0, trust_radius: float = 0.0, matrix_rcond: float = 1e-12, matrix_condition_cap: float = 1e12) -> dict[str, Any]:
    del lazy_bad_cond
    t0 = now()
    deadline = t0 + float(trial_timeout) if trial_timeout and trial_timeout > 0 else None
    y = np.asarray(y0, dtype=np.complex128).copy()
    best_y, best_r = y.copy(), base.residual(y)
    status, ok, epochs_done = "started", False, 0
    scales = list(irp_scales or [1.0 + 0.0j])
    total_line = total_hyper = hyper_used = triggers = rescues = direct_steps = skipped = chart_switches = recip_uses = direct_uses = 0
    collapsed, locality_hits, stagnation_hits = False, 0, 0
    collapse_epoch = collapse_reason = None
    last: dict[str, Any] = {}
    last_chart: dict[str, Any] = {"kind": "direct", "scale": 1.0 + 0j, "score": best_r, "log_energy": log_energy(y)}

    def absorb(loc: dict[str, Any]) -> None:
        nonlocal total_line, total_hyper, hyper_used
        total_line += int(loc.get("line_search_evals", 0) or 0)
        total_hyper += int(loc.get("hypercube_total_evals", 0) or 0)
        hyper_used += int(loc.get("hypercube_used_count", 0) or 0)

    for ep in range(max(1, int(max_epochs))):
        if deadline is not None and now() > deadline:
            status = "timeout"
            break
        r0 = base.residual(y)
        if r0 < best_r:
            best_y, best_r = y.copy(), r0
        if r0 <= max(float(tol), float(accept)) and (accept <= 0 or r0 < accept):
            ok, status = True, "converged"
            break
        loc = hypercube_matrixjet_corrector(base, y, max(1, int(local_inner_epochs if collapsed else lazy_direct_epochs)), tol, accept, 0.0 if deadline is None else max(0.0, deadline - now()), line_search, line_grid, int(direction_seed) + 1000003 * ep + 17, cloud_nodes, lm_damping, trust_radius, matrix_rcond, matrix_condition_cap)
        absorb(loc)
        last = dict(loc)
        yd = np.asarray(loc.get("y", y), dtype=np.complex128)
        rd = base.residual(yd)
        improved = math.isfinite(rd) and rd <= r0 * (1.0 - 1e-14)
        if improved:
            direct_steps += 1
            stagnation_hits = 0
            if rd < best_r:
                best_y, best_r = yd.copy(), rd
        else:
            stagnation_hits += 1
        strong = bool(improved and rd / max(r0, 1e-300) <= float(lazy_trigger_drop) and log_energy(yd) <= float(lazy_log_energy))
        trigger = bool(eager_irp and not collapsed) or (not improved and stagnation_hits >= max(1, int(lazy_trigger_after))) or (improved and not strong)
        if collapsed and not rescue_collapsed:
            trigger = False
        if not trigger:
            if improved:
                step = float(np.linalg.norm(yd - y)) / max(1.0, float(np.linalg.norm(y)))
                drop = rd / max(r0, 1e-300)
                y = yd.copy()
                direct_uses += 1
                epochs_done = ep + 1
                if collapse and not collapsed:
                    if rd <= float(collapse_residual):
                        collapsed, collapse_epoch, collapse_reason = True, ep, "residual-threshold"
                    elif drop <= float(collapse_drop) and step <= float(collapse_rel_step):
                        locality_hits += 1
                        if locality_hits >= max(1, int(collapse_after)):
                            collapsed, collapse_epoch, collapse_reason = True, ep, "observed-local-contraction"
                status = "local-collapsed" if collapsed else "direct-410-pandrosion-irp-gemm-running"
                continue
            status = "no-lazy-direct-decrease"
            break
        triggers += 1
        rescues += 1 if not improved else 0
        base_y, base_r = (yd.copy(), rd) if improved else (y.copy(), r0)
        irp_y, irp_r = base_y.copy(), base_r
        any_irp = False
        for layer in range(max(1, int(irp_layers))):
            cands = chart_candidates(base, irp_y, scales, bool(irp_inversion), int(irp_chart_top))
            if not cands:
                status = "no-admissible-lazy-irp-chart"
                break
            best_cy = None
            best_cr = irp_r
            best_cloc: dict[str, Any] = {}
            best_chart: dict[str, Any] = {}
            for _, ct, u0, cm in cands:
                loc2 = hypercube_matrixjet_corrector(ct, u0, max(1, int(irp_inner_epochs)), tol, accept, 0.0, line_search, line_grid, int(direction_seed) + 1000003 * ep + 65537 * layer, cloud_nodes, lm_damping, trust_radius, matrix_rcond, matrix_condition_cap)
                absorb(loc2)
                yb = ct.to_base(loc2.get("y", u0))
                rb = base.residual(yb)
                if math.isfinite(rb) and rb < best_cr:
                    best_cy, best_cr, best_cloc, best_chart = yb.copy(), rb, dict(loc2), dict(cm)
            if best_cy is None:
                status = "lazy-irp-no-extra-decrease"
                break
            irp_y, irp_r, last, last_chart = best_cy, best_cr, best_cloc, best_chart
            any_irp = True
            chart_switches += 1
            recip_uses += 1 if best_chart.get("kind") == "reciprocal" else 0
            direct_uses += 0 if best_chart.get("kind") == "reciprocal" else 1
            if irp_r < best_r:
                best_y, best_r = irp_y.copy(), irp_r
            if irp_r <= max(float(tol), float(accept)) and (accept <= 0 or irp_r < accept):
                ok, status = True, "converged"
                break
        if ok:
            break
        if any_irp and irp_r < r0:
            y = irp_y.copy()
            epochs_done = ep + 1
            status = "lazy-irp-running"
            continue
        if improved:
            y = yd.copy()
            epochs_done = ep + 1
            skipped += 1
            status = "direct-after-weak-irp"
            continue
        status = status if status != "started" else "no-lazy-irp-or-direct-decrease"
        break
    final_r = base.residual(best_y)
    if final_r <= max(float(tol), float(accept)) and (accept <= 0 or final_r < accept):
        ok, status, best_r = True, "converged", final_r
    return {"accepted": bool(ok if accept <= 0 else (math.isfinite(best_r) and best_r < accept)), "ok": bool(ok), "status": status, "epochs": int(epochs_done), "residual": float(best_r), "y": best_y, "seconds": float(now() - t0), "slope_cond": None, "corrector": "411-pandrosion-lazy-irp-plus-cached-hypercube-gemm-matrixjet", "line_search_evals": int(total_line), "hypercube_total_evals": int(total_hyper), "hypercube_used_count": int(hyper_used), "matrix_rcond": float(matrix_rcond), "matrix_condition_cap": float(matrix_condition_cap), "irp_layers_completed": int(chart_switches), "irp_chart_switches": int(chart_switches), "irp_collapsed": bool(collapsed), "irp_collapse_epoch": collapse_epoch, "irp_collapse_reason": collapse_reason, "irp_reciprocal_uses": int(recip_uses), "irp_direct_uses": int(direct_uses), "irp_lazy_triggers": int(triggers), "irp_lazy_rescues": int(rescues), "irp_lazy_direct_steps": int(direct_steps), "irp_lazy_direct_good": int(direct_steps - triggers if direct_steps >= triggers else 0), "irp_lazy_direct_weak": int(triggers), "irp_lazy_skipped": int(skipped), "irp_last_chart_kind": last_chart.get("kind"), "irp_last_chart_scale": cjson(last_chart.get("scale", 1.0 + 0j)), "irp_last_chart_log_energy": last_chart.get("log_energy"), "halley_enabled": False, "halley_total_evals": 0, "halley_used_count": 0, **{k: v for k, v in last.items() if k != "y"}}


def polish_geometric_root(geometry: Any, chart: LinearChart, y: Any, accept: float, max_steps: int, matrix_rcond: float = 1e-12, matrix_condition_cap: float = 1e12, lm_damping: float = 0.0) -> tuple[Any, Any, float, dict[str, Any]]:
    """Matrix polish on the GEOMETRY (finite-difference slope of F_geo)."""
    z = chart.z_from_y(y)
    r = finite_norm(geometry.eval(z))
    meta = {"geometric_polish_enabled": bool(max_steps > 0), "geometric_polish_steps": 0, "geometric_polish_r0": float(r), "geometric_polish_r1": float(r)}
    if max_steps <= 0 or not math.isfinite(r):
        return z, np.asarray(y, dtype=np.complex128), r, meta
    best_z, best_f, best_r = np.asarray(z, dtype=np.complex128).copy(), geometry.eval(z), r
    for step in range(int(max_steps)):
        if best_r <= float(accept):
            break
        try:
            J = geometry.slope_matrix(best_z, best_z)
            delta, mmeta = matrix_svd_step(np.asarray(J, dtype=np.complex128), -best_f, float(matrix_rcond), float(matrix_condition_cap), float(lm_damping))
        except Exception:
            meta["geometric_polish_status"] = "matrix-filter-failed"
            break
        meta.update({f"geometric_polish_{k}": v for k, v in mmeta.items()})
        improved = False
        for gain in (1.0, 0.5, 0.25, 0.125, 0.0625):
            cand_z = best_z + gain * delta
            cand_f = geometry.eval(cand_z)
            cand_r = finite_norm(cand_f)
            if cand_r < best_r:
                best_z, best_f, best_r = cand_z, cand_f, cand_r
                meta["geometric_polish_steps"] = step + 1
                improved = True
                break
        if not improved:
            meta["geometric_polish_status"] = "no-decrease"
            break
    meta["geometric_polish_r1"] = float(best_r)
    return best_z, chart.y_from_z(best_z), best_r, meta


def _torch_batched_ns_inverse(B: Any, iters: int) -> Any:
    n = int(B.shape[-1])
    I = torch.eye(n, dtype=B.dtype, device=B.device)[None, :, :]
    norm1 = torch.max(torch.sum(torch.abs(B), dim=-2), dim=-1).values[:, None, None]
    norminf = torch.max(torch.sum(torch.abs(B), dim=-1), dim=-1).values[:, None, None]
    alpha = torch.clamp(norm1 * norminf, min=torch.finfo(torch.float32).tiny)
    X = B.conj().transpose(-2, -1) / alpha
    for _ in range(max(1, int(iters))):
        X = torch.bmm(X, 2.0 * I - torch.bmm(B, X))
    return X


def torch_batched_gemm_fit_410(X: Any, Y: Any) -> Any:
    XX_np = np.asarray(X, dtype=np.complex128)
    YY_np = np.asarray(Y, dtype=np.complex128)
    device, dtype = torch_device(), torch_complex_dtype()
    XX = torch.as_tensor(XX_np, dtype=dtype, device=device)
    YY = torch.as_tensor(YY_np, dtype=dtype, device=device)
    Xh = XX.conj().transpose(-2, -1)
    G = torch.bmm(Xh, XX)
    RHS = torch.bmm(Xh, YY)
    n = int(G.shape[-1])
    diag = torch.real(torch.diagonal(G, dim1=-2, dim2=-1)).mean(dim=-1).clamp_min(torch.finfo(torch.float32).tiny)
    mu = (float(NS_DAMPING_FLOOR) * diag).to(dtype)[:, None, None]
    I = torch.eye(n, dtype=dtype, device=device)[None, :, :]
    Binv = _torch_batched_ns_inverse(G + mu * I, NS_ITERS)
    return torch.bmm(Binv, RHS).detach().cpu().numpy().astype(np.complex128, copy=False)


def torch_batched_gemm_step_410(A: Any, rhs: Any, damping: float = 0.0) -> Any:
    AA_np = np.asarray(A, dtype=np.complex128)
    bb_np = np.asarray(rhs, dtype=np.complex128)
    device, dtype = torch_device(), torch_complex_dtype()
    AA = torch.as_tensor(AA_np, dtype=dtype, device=device)
    bb = torch.as_tensor(bb_np, dtype=dtype, device=device)
    Ah = AA.conj().transpose(-2, -1)
    G = torch.bmm(Ah, AA)
    g = torch.bmm(Ah, bb[:, :, None])
    n = int(G.shape[-1])
    diag = torch.real(torch.diagonal(G, dim1=-2, dim2=-1)).mean(dim=-1).clamp_min(torch.finfo(torch.float32).tiny)
    mu = (max(float(damping), float(NS_DAMPING_FLOOR)) * diag).to(dtype)[:, None, None]
    I = torch.eye(n, dtype=dtype, device=device)[None, :, :]
    Binv = _torch_batched_ns_inverse(G + mu * I, NS_ITERS)
    return torch.bmm(Binv, g)[:, :, 0].detach().cpu().numpy().astype(np.complex128, copy=False)


def torch_batched_gemm_candidates_410(dY: Any, dF: Any, F: Any, Y: Any, lambdas: Any, damping: float = 0.0, trust_radius: float = 0.0, yn: Optional[Any] = None) -> Any:
    """Fused 410 batch prepass kernel: fit local jets, step, and candidates in Torch."""
    dY_np = np.asarray(dY, dtype=np.complex128)
    dF_np = np.asarray(dF, dtype=np.complex128)
    F_np = np.asarray(F, dtype=np.complex128)
    Y_np = np.asarray(Y, dtype=np.complex128)
    L_np = np.asarray(lambdas, dtype=np.float32)
    device, dtype = torch_device(), torch_complex_dtype()
    X = torch.as_tensor(dY_np, dtype=dtype, device=device)
    Z = torch.as_tensor(dF_np, dtype=dtype, device=device)
    Fv = torch.as_tensor(F_np, dtype=dtype, device=device)
    Yt = torch.as_tensor(Y_np, dtype=dtype, device=device)
    L = torch.as_tensor(L_np, dtype=torch.float32, device=device).to(dtype)

    Xh = X.conj().transpose(-2, -1)
    Gfit = torch.bmm(Xh, X)
    RHSfit = torch.bmm(Xh, Z)
    n = int(Gfit.shape[-1])
    diag_fit = torch.real(torch.diagonal(Gfit, dim1=-2, dim2=-1)).mean(dim=-1).clamp_min(torch.finfo(torch.float32).tiny)
    mu_fit = (float(NS_DAMPING_FLOOR) * diag_fit).to(dtype)[:, None, None]
    I = torch.eye(n, dtype=dtype, device=device)[None, :, :]
    coef = torch.bmm(_torch_batched_ns_inverse(Gfit + mu_fit * I, NS_ITERS), RHSfit)
    A = coef.transpose(-2, -1)

    Ah = A.conj().transpose(-2, -1)
    G = torch.bmm(Ah, A)
    rhs = -Fv
    g = torch.bmm(Ah, rhs[:, :, None])
    diag = torch.real(torch.diagonal(G, dim1=-2, dim2=-1)).mean(dim=-1).clamp_min(torch.finfo(torch.float32).tiny)
    mu = (max(float(damping), float(NS_DAMPING_FLOOR)) * diag).to(dtype)[:, None, None]
    delta = torch.bmm(_torch_batched_ns_inverse(G + mu * I, NS_ITERS), g)[:, :, 0]
    if float(trust_radius) > 0.0:
        yy = torch.as_tensor(np.asarray(yn, dtype=np.float32), dtype=torch.float32, device=device) if yn is not None else torch.linalg.vector_norm(Yt, dim=1).real.clamp_min(1.0)
        lim = float(trust_radius) * yy
        dn = torch.linalg.vector_norm(delta, dim=1).real
        scale = torch.where((dn > lim) & (dn > 0), lim / torch.clamp(dn, min=torch.finfo(torch.float32).tiny), torch.ones_like(dn))
        delta = delta * scale.to(dtype)[:, None]
    C = Yt[:, None, :] + L[None, :, None] * delta[:, None, :]
    return C.detach().cpu().numpy().astype(np.complex128, copy=False)


def torch_batched_realpack_candidates_410(dY: Any, dF: Any, F: Any, Y: Any, lambdas: Any, damping: float = 0.0, trust_radius: float = 0.0, yn: Optional[Any] = None) -> Any:
    """Fused 410 real-packed prepass kernel using real GEMM/Newton-Schulz."""
    dY_np = np.asarray(dY, dtype=np.complex128)
    dF_np = np.asarray(dF, dtype=np.complex128)
    F_np = np.asarray(F, dtype=np.complex128)
    Y_np = np.asarray(Y, dtype=np.complex128)
    L_np = np.asarray(lambdas, dtype=np.float32)
    B, _M, n = int(dY_np.shape[0]), int(dY_np.shape[1]), int(dY_np.shape[2])
    device, dtype = torch_device(), torch_real_dtype()
    X = _torch_pack_complex_design_batched(dY_np, dtype, device)
    Zr = torch.as_tensor(dF_np.real, dtype=dtype, device=device)
    Zi = torch.as_tensor(dF_np.imag, dtype=dtype, device=device)
    Z = torch.cat((Zr, Zi), dim=1)
    Fr = torch.as_tensor(F_np.real, dtype=dtype, device=device)
    Fi = torch.as_tensor(F_np.imag, dtype=dtype, device=device)
    Yr = torch.as_tensor(Y_np.real, dtype=dtype, device=device)
    Yi = torch.as_tensor(Y_np.imag, dtype=dtype, device=device)
    L = torch.as_tensor(L_np, dtype=dtype, device=device)

    Xh = X.transpose(-2, -1)
    Gfit = torch.bmm(Xh, X)
    RHSfit = torch.bmm(Xh, Z)
    n2 = int(Gfit.shape[-1])
    diag_fit = torch.diagonal(Gfit, dim1=-2, dim2=-1).mean(dim=-1).clamp_min(_torch_real_tiny(dtype))
    mu_fit = (float(NS_DAMPING_FLOOR) * diag_fit)[:, None, None]
    I2 = torch.eye(n2, dtype=dtype, device=device)[None, :, :]
    coef_stack = torch.bmm(_torch_batched_ns_inverse(Gfit + mu_fit * I2, NS_ITERS), RHSfit)

    P = coef_stack[:, :n, :]
    Q = coef_stack[:, n:, :]
    Are = P.transpose(-2, -1)
    Aim = Q.transpose(-2, -1)
    A = torch.cat((torch.cat((Are, -Aim), dim=-1), torch.cat((Aim, Are), dim=-1)), dim=-2)
    rhs = -torch.cat((Fr, Fi), dim=-1)
    delta_pack, _mu = _torch_realpack_batched_normal_step(A, rhs, float(damping))
    dr = delta_pack[:, :n]
    di = delta_pack[:, n:]
    if float(trust_radius) > 0.0:
        if yn is not None:
            yy = torch.as_tensor(np.asarray(yn, dtype=np.float32), dtype=dtype, device=device)
        else:
            yy = torch.sqrt(torch.sum(Yr * Yr + Yi * Yi, dim=1)).clamp_min(1.0)
        lim = float(trust_radius) * yy
        dn = torch.sqrt(torch.sum(dr * dr + di * di, dim=1))
        scale = torch.where((dn > lim) & (dn > 0), lim / torch.clamp(dn, min=_torch_real_tiny(dtype)), torch.ones_like(dn))
        dr = dr * scale[:, None]
        di = di * scale[:, None]
    Cr = Yr[:, None, :] + L[None, :, None] * dr[:, None, :]
    Ci = Yi[:, None, :] + L[None, :, None] * di[:, None, :]
    return _real_tensor_to_numpy(Cr) + 1j * _real_tensor_to_numpy(Ci)


def torch_batched_sketch_candidates_410(dY: Any, dF: Any, F: Any, Y: Any, lambdas: Any, damping: float = 0.0, trust_radius: float = 0.0, yn: Optional[Any] = None) -> Any:
    """Fused 410 subspace-sketch prepass kernel.

    Fit B=P C from dF ~= dY B, then solve the local Newton step as delta=P u.
    The heavy normal equations are k-by-k with k=sketch_dim_for(n).
    """
    dY_np = np.asarray(dY, dtype=np.complex128)
    dF_np = np.asarray(dF, dtype=np.complex128)
    F_np = np.asarray(F, dtype=np.complex128)
    Y_np = np.asarray(Y, dtype=np.complex128)
    L_np = np.asarray(lambdas, dtype=np.float32)
    B, _M, n = int(dY_np.shape[0]), int(dY_np.shape[1]), int(dY_np.shape[2])
    k = sketch_dim_for(n)
    if k >= n:
        return torch_batched_gemm_candidates_410(dY_np, dF_np, F_np, Y_np, L_np, damping, trust_radius, yn)
    device, dtype = torch_device(), torch_complex_dtype()
    X = torch.as_tensor(dY_np, dtype=dtype, device=device)
    Z = torch.as_tensor(dF_np, dtype=dtype, device=device)
    Fv = torch.as_tensor(F_np, dtype=dtype, device=device)
    Yt = torch.as_tensor(Y_np, dtype=dtype, device=device)
    L = torch.as_tensor(L_np, dtype=torch.float32, device=device).to(dtype)
    P = _torch_sketch_basis(n, k, dtype, device, salt=int(B) + int(dY_np.shape[1]) + 23)
    XP = torch.matmul(X, P)

    Xh = XP.conj().transpose(-2, -1)
    Gfit = torch.bmm(Xh, XP)
    RHSfit = torch.bmm(Xh, Z)
    diag_fit = torch.real(torch.diagonal(Gfit, dim1=-2, dim2=-1)).mean(dim=-1).clamp_min(torch.finfo(torch.float32).tiny)
    mu_fit = (float(NS_DAMPING_FLOOR) * diag_fit).to(dtype)[:, None, None]
    Ik = torch.eye(k, dtype=dtype, device=device)[None, :, :]
    Cfit = torch.bmm(_torch_batched_ns_inverse(Gfit + mu_fit * Ik, NS_ITERS), RHSfit)

    # AP is the full local operator restricted to the sketch subspace.
    PtP = P.transpose(-2, -1) @ P
    AP = torch.matmul(Cfit.transpose(-2, -1), PtP)
    Ah = AP.conj().transpose(-2, -1)
    G = torch.bmm(Ah, AP)
    rhs = -Fv
    g = torch.bmm(Ah, rhs[:, :, None])
    diag = torch.real(torch.diagonal(G, dim1=-2, dim2=-1)).mean(dim=-1).clamp_min(torch.finfo(torch.float32).tiny)
    mu = (max(float(damping), float(NS_DAMPING_FLOOR)) * diag).to(dtype)[:, None, None]
    u = torch.bmm(_torch_batched_ns_inverse(G + mu * Ik, NS_ITERS), g)[:, :, 0]
    delta = torch.matmul(P[None, :, :], u[:, :, None])[:, :, 0]
    if float(trust_radius) > 0.0:
        yy = torch.as_tensor(np.asarray(yn, dtype=np.float32), dtype=torch.float32, device=device) if yn is not None else torch.linalg.vector_norm(Yt, dim=1).real.clamp_min(1.0)
        lim = float(trust_radius) * yy
        dn = torch.linalg.vector_norm(delta, dim=1).real
        scale = torch.where((dn > lim) & (dn > 0), lim / torch.clamp(dn, min=torch.finfo(torch.float32).tiny), torch.ones_like(dn))
        delta = delta * scale.to(dtype)[:, None]
    C = Yt[:, None, :] + L[None, :, None] * delta[:, None, :]
    return C.detach().cpu().numpy().astype(np.complex128, copy=False)


def selected_batch_kernel_410() -> str:
    if BATCH_KERNEL in {"sketch", "realpack", "complex"}:
        return BATCH_KERNEL
    selected_algorithm = selected_matrix_algorithm_410(128)
    if selected_algorithm == "sketch-ns":
        return "sketch"
    return "realpack" if selected_algorithm in {"realpack-ns", "auto"} else "complex"


def batch_direct_gemm_prepass_410(target: TargetTrack, Y0: Any, max_epochs: int, tol: float, accept: float, line_search: int, line_grid: Sequence[float], direction_seed: int, cloud_nodes: int, lm_damping: float, trust_radius: float) -> list[dict[str, Any]]:
    """Batched direct matrix prepass before the normal Pandrosion IRP corrector.

    This does not replace IRP. It only moves a wave of starts with the same
    hypercube/GEMM local jet operation so the subsequent IRP calls start closer.
    """
    Y = np.asarray(Y0, dtype=np.complex128).copy()
    if Y.ndim != 2 or Y.shape[0] == 0:
        return []
    B, n = int(Y.shape[0]), int(Y.shape[1])
    best_y = Y.copy()
    best_r = target.residuals_batch(Y)
    lambdas = np.asarray(_line_lambdas(line_search, line_grid), dtype=float)
    total_line = 0
    total_hyper = 0
    moved = np.zeros(B, dtype=np.int64)
    last_meta: list[dict[str, Any]] = [{} for _ in range(B)]
    configured_kernel = selected_batch_kernel_410()
    fallback_error: Optional[str] = None
    for ep in range(max(0, int(max_epochs))):
        F = target.eval_batch(Y)
        R = safe_norms(F)
        improve_best = np.isfinite(R) & (R < best_r)
        best_y[improve_best] = Y[improve_best]
        best_r[improve_best] = R[improve_best]
        active = np.where(np.isfinite(R) & (R > max(float(tol), float(accept))))[0]
        if active.size == 0:
            break
        Ya = Y[active]
        Fa = F[active]
        M = max(int(cloud_nodes), 2 * n + 4, 16)
        yn = np.maximum(1.0, np.linalg.norm(Ya, axis=1).astype(float))
        h = 1e-5 * yn
        rng = np.random.default_rng(int(direction_seed) + 104729 * ep + 17 * B)
        signs = rng.choice([-1.0, 1.0], size=(active.size, M, n))
        dY = h[:, None, None] * signs
        P = Ya[:, None, :] + dY
        dF = target.eval_batch(P.reshape(-1, n)).reshape(active.size, M, n) - Fa[:, None, :]
        kernel_used = configured_kernel
        try:
            if configured_kernel == "sketch":
                C = torch_batched_sketch_candidates_410(dY, dF, Fa, Ya, lambdas, float(lm_damping), float(trust_radius), yn)
            elif configured_kernel == "realpack":
                C = torch_batched_realpack_candidates_410(dY, dF, Fa, Ya, lambdas, float(lm_damping), float(trust_radius), yn)
            else:
                C = torch_batched_gemm_candidates_410(dY, dF, Fa, Ya, lambdas, float(lm_damping), float(trust_radius), yn)
        except Exception as exc:
            if BATCH_KERNEL != "auto" or configured_kernel not in {"sketch", "realpack"}:
                raise
            fallback_error = f"{type(exc).__name__}:{exc}"
            kernel_used = "complex-fallback"
            C = torch_batched_gemm_candidates_410(dY, dF, Fa, Ya, lambdas, float(lm_damping), float(trust_radius), yn)
        Rc = target.residuals_batch(C.reshape(-1, n)).reshape(active.size, len(lambdas))
        total_line += int(C.shape[0] * C.shape[1])
        total_hyper += int(active.size * M)
        for row, idx in enumerate(active):
            base_r = float(R[idx])
            scores = np.where(np.isfinite(Rc[row]) & (Rc[row] < base_r), Rc[row], np.inf)
            if np.any(np.isfinite(scores)):
                j = int(np.nanargmin(scores))
                Y[idx] = C[row, j]
                moved[idx] += 1
                if float(scores[j]) < best_r[idx]:
                    best_y[idx] = Y[idx]
                    best_r[idx] = float(scores[j])
                last_meta[idx] = {"batch_prepass_kernel": str(kernel_used), "batch_prepass_last_lambda": float(lambdas[j]), "batch_prepass_last_residual": float(scores[j])}
    out = []
    for i in range(B):
        out.append({
            "y": best_y[i].copy(),
            "batch_prepass_enabled": True,
            "batch_prepass_epochs": int(max_epochs),
            "batch_prepass_wave": int(B),
            "batch_prepass_moved": int(moved[i]),
            "batch_prepass_residual": float(best_r[i]),
            "batch_prepass_line_evals_total": int(total_line),
            "batch_prepass_hypercube_evals_total": int(total_hyper),
            "batch_prepass_algorithm": f"410-{configured_kernel}-fused-gemm-batched-direct-before-irp",
            "batch_prepass_kernel": str(configured_kernel),
            "batch_prepass_fallback_error": fallback_error,
            **last_meta[i],
        })
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    t_case = now()
    n, d = parse_case(case_raw)
    geometry, base, source, backend = make_system_410(args, n, d)
    chart = LinearChart.identity(n, float(args.linear_scale))
    target = TargetTrack(geometry, chart)  # the WHOLE pandrosion stack runs on the geometry
    starts = parse_start_points(args.starts, n)
    # The 'complete logarithmic ladder' (proportional-means / x^{k/p} ladder) is the default
    # scale geometry inherited from 317. It gives even log coverage so a start lands in an easy basin.
    ladder_on = bool(args.scale_ladder)
    default_powers = monomial_scale_ladder(int(args.ladder_subdiv), int(args.ladder_octaves), float(args.ladder_base)) if ladder_on else DEFAULT_POWERS
    powers = sorted(set(round(float(x), 16) for x in parse_float_list(args.powers, default_powers, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    angles = [math.radians(x) for x in parse_float_list(args.angles, DEFAULT_ANGLES_DEG)]
    radii = parse_float_list(args.rays, DEFAULT_RADII, positive=True)
    default_startopt_gains = sorted(monomial_scale_ladder(int(args.ladder_subdiv), min(3, max(1, int(args.ladder_octaves))), float(args.ladder_base)), key=lambda g: abs(math.log(g))) if ladder_on else DEFAULT_GAINS
    gains = parse_float_list(args.startopt_gains, default_startopt_gains, positive=True)
    startopt_gains_source = "user" if args.startopt_gains not in (None, "") else ("pandrosion-ladder" if ladder_on else "legacy")
    # IRP chart palette: a tighter proportional-means ladder around 1, ordered by |log s| so the
    # near-basin scales are tried first within --irp-top, else the legacy dyadic palette.
    default_irp_gains = sorted(monomial_scale_ladder(int(args.ladder_subdiv), min(3, max(1, int(args.ladder_octaves))), float(args.ladder_base)), key=lambda g: abs(math.log(g))) if ladder_on else [1.0, 0.5, 2.0, 0.25, 4.0, 0.125, 8.0]
    irp_gain_values = parse_float_list(args.irp_gains, default_irp_gains, positive=True)
    irp_gains_source = "user" if args.irp_gains not in (None, "") else ("pandrosion-ladder" if ladder_on else "legacy")
    irp_scales = complex_scale_palette(irp_gain_values, parse_float_list(args.irp_phases, [0.0, 0.08, -0.08, 0.19, -0.19]), int(args.irp_top))
    roots: list[dict[str, Any]] = []
    trials: list[dict[str, Any]] = []
    failures = duplicates = 0
    atlas_calls = 0
    t_extract = now()

    def build_trial_start(trial: int) -> dict[str, Any]:
        nonlocal atlas_calls
        explicit = trial < len(starts)
        if explicit:
            y_raw = np.asarray(starts[trial], dtype=np.complex128)
            geom = {"chart": "410-explicit-start", "atlas_mode": "explicit-start", "homothety": 1.0, "theta_deg": None}
        else:
            y_raw = None
            geom: dict[str, Any] = {}
            if bool(args.origin_seed) and (trial == 0 if int(args.origin_seed_period) <= 0 else trial % max(1, int(args.origin_seed_period)) == 0):
                cap = float(args.origin_seed_max_norm) if float(args.origin_seed_max_norm) > 0 else 2.0 * math.sqrt(max(1, n))
                y_origin, ometa = origin_affine_start(target, n, float(args.origin_seed_h), cap, float(args.matrix_rcond), float(args.matrix_condition_cap))
                if y_origin is not None:
                    y_raw = y_origin
                    geom = {"chart": "410-origin-affine-gemm-matrix-seed", "atlas_mode": "origin-affine-gemm-matrix-seed", "homothety": 1.0, "theta_deg": None, **ometa}
            if y_raw is None:
                y_raw, geom = universal_atlas_start(target, n, atlas_calls, geometry.seed + 0x113000, powers, angles, radii, float(args.power_cap), len(roots), duplicates, failures, int(args.count), int(args.universal_cells), int(args.universal_shells), atlas_selection=str(args.atlas_selection), cell_probe_radius=float(args.cell_probe_radius), cell_descent_min=float(args.cell_descent_min), cell_equal_gap_min=float(args.cell_equal_gap_min), cell_log_max=float(args.cell_log_max), universal_cycle=bool(args.universal_cycle))
                atlas_calls += 1
        bypass_startopt = bool(args.atlas_bypass_startopt) and bool(geom.get("atlas_startopt_bypass_recommended")) and args.startopt_gains in (None, "")
        if bypass_startopt:
            y0 = np.asarray(y_raw, dtype=np.complex128).copy()
            r0 = target.residual(y0)
            smeta = {"startopt_enabled": False, "startopt_skipped": "atlas-diversity", "startopt_r0": float(r0), "startopt_r1": float(r0), "startopt_ratio": 1.0 if math.isfinite(r0) and r0 > 0 else None, "startopt_steps": 0, "startopt_evals": 1, "startopt_micro_epochs": 0, "startopt_gain": 1.0, "startopt_batch_numpy": False}
        else:
            y0, smeta = startopt(target, y_raw, trial, geometry.seed + 0x112555, int(args.startopt_steps), int(args.startopt_candidates), gains, int(args.startopt_micro_epochs))
        return {"trial": int(trial), "y0": y0, "geom": geom, "smeta": smeta}

    def finish_trial(item: dict[str, Any]) -> None:
        nonlocal failures, duplicates
        trial = int(item["trial"])
        y0 = np.asarray(item["y0"], dtype=np.complex128)
        geom = dict(item["geom"])
        smeta = dict(item["smeta"])
        loc = lazy_irp_hypercube_torch_gemm_matrixjet_corrector_410(target, y0, int(args.epochs), float(args.tol), float(args.accept), float(args.trial_timeout), int(args.line_search), parse_float_list(args.line_grid, []), geometry.seed + 7919 * trial, int(args.hypercube_nodes), int(args.irp_layers), int(args.irp_inner_epochs), irp_scales, int(args.irp_chart_top), bool(args.irp_inversion), bool(args.collapse), float(args.collapse_residual), float(args.collapse_drop), float(args.collapse_rel_step), int(args.collapse_after), int(args.local_inner_epochs), int(args.lazy_direct_epochs), float(args.lazy_trigger_drop), int(args.lazy_trigger_after), float(args.lazy_bad_cond), float(args.lazy_log_energy), bool(args.eager_irp), bool(args.rescue_collapsed), float(args.lm_damping), float(args.trust_radius), float(args.matrix_rcond), float(args.matrix_condition_cap))
        z, y_final, geo_residual, polish_meta = polish_geometric_root(geometry, chart, loc["y"], float(args.accept), int(args.geometric_polish_steps), float(args.matrix_rcond), float(args.matrix_condition_cap), float(args.lm_damping))
        accepted = bool(math.isfinite(geo_residual) and geo_residual < float(args.accept))
        rec = {"trial": int(trial), "accepted": accepted, "status": loc.get("status"), "r1": float(geo_residual), "geometric_residual": float(geo_residual), "epochs": int(loc.get("epochs", 0)), "seconds": float(loc.get("seconds", 0.0)), **{k: v for k, v in loc.items() if k != "y"}, **geom, **smeta, **polish_meta}
        if bool(args.verbose_trials):
            rec["z"] = root_to_json(z)
            rec["y0"] = root_to_json(y0)
        if not accepted:
            failures += 1
            trials.append(rec)
            return
        dup = cluster_index(roots, z, float(args.cluster_sep))
        if dup is not None:
            duplicates += 1
            rec["status"] = "duplicate"
            rec["cluster"] = int(dup)
            trials.append(rec)
            return
        root = {"id": len(roots), "source": "410-standalone-pandrosion-irp-batched-local-jet-gemm-matrix", "trial": int(trial), "z_complex": np.asarray(z, dtype=np.complex128).copy(), "y_complex": np.asarray(y_final, dtype=np.complex128).copy(), "residual": float(geo_residual), "geometric_residual": float(geo_residual), "realness": realness(z), "cond": None, "epochs": int(loc.get("epochs", 0)), "seconds": float(loc.get("seconds", 0.0)), **{k: v for k, v in loc.items() if k != "y"}, **geom, **smeta, **polish_meta}
        root["score"] = score_root(float(root["residual"]), float(root["realness"]), root["cond"])
        roots.append(root)
        rec["status"] = "new-root"
        rec["root_id"] = int(root["id"])
        trials.append(rec)

    trial = 0
    batch_wave = max(1, int(getattr(args, "batch_wave", 1)))
    while trial < int(args.pool) and len(roots) < int(args.count):
        wave: list[dict[str, Any]] = []
        for _ in range(batch_wave):
            if trial >= int(args.pool) or len(roots) + len(wave) >= int(args.count):
                break
            wave.append(build_trial_start(trial))
            trial += 1
        if not wave:
            break
        use_prepass = bool(args.batch_prepass) and len(wave) > 1 and int(args.batch_prepass_epochs) > 0 and MATRIX_BACKEND in {"torch", "auto"} and MATRIX_ALGORITHM in {"adaptive-directional-sketch", "directional-sketch", "adaptive-sketch", "sketch-ns", "adaptive-ns", "realpack-ns", "gemm-ns", "auto"}
        if use_prepass:
            pre = batch_direct_gemm_prepass_410(target, np.asarray([w["y0"] for w in wave], dtype=np.complex128), int(args.batch_prepass_epochs), float(args.tol), float(args.accept), int(args.line_search), parse_float_list(args.line_grid, []), geometry.seed + 0x410000 + 8191 * trial, int(args.hypercube_nodes), float(args.lm_damping), float(args.trust_radius))
            for w, pmeta in zip(wave, pre):
                w["y0"] = pmeta.pop("y")
                w["smeta"] = {**w["smeta"], **pmeta}
        else:
            for w in wave:
                w["smeta"] = {**w["smeta"], "batch_prepass_enabled": False, "batch_prepass_wave": int(len(wave))}
        for item in wave:
            if len(roots) >= int(args.count):
                break
            finish_trial(item)
    encoded_roots = []
    for root in sorted(roots, key=lambda q: (float(q.get("score", float("inf"))), int(q.get("id", 0)))):
        rr = dict(root)
        rr["z"] = root_to_json(rr.pop("z_complex"))
        rr["y"] = root_to_json(rr.pop("y_complex"))
        encoded_roots.append(rr)
    gstats = geometry.stats()
    extract_seconds = float(now() - t_extract)
    samples_per_trial = (float(geometry.oracle_samples) / max(1, len(trials))) if trials else 0.0
    result = {
        "script": Path(__file__).name, "autonomous": True, "dependencies": {"python_scripts": [], "numpy": True, "torch": True},
        "mode": "411-standalone-pandrosion-irp-cached-directional-sketch-local-jet-matrix",
        "flow_formula": "any system as black-box geometric oracle -> lazily-sampled local-jet field -> universal Mobius atlas starts -> StartOpt -> sketch/realpack/full GEMM direct prepass -> Pandrosion lazy IRP charts -> directional JP oracle samples when k<<n, otherwise full hypercube matrix-jet -> matrix polish -> faithful residual acceptance",
        "case": f"{n},{d}", "family": source, "system_source": source, "base_backend": backend, "fully_geometric": True,
        "geometry_kind": "local-jet-field", "residual_is_faithful": True,
        "seed_index": int(args.seed_index), "seed": int(geometry.seed), "n": int(n), "degree": int(d),
        "terms_per_poly": int(geometry.terms_per_poly), "terms": int(geometry.total_terms), "bezout": int(geometry.bezout),
        "equation_normalize": bool(args.equation_normalize),
        "linear_A": [[cjson(chart.A[i, j]) for j in range(n)] for i in range(n)],
        "geometry": {"kind": "local-jet-field", "use_quadratic": bool(args.jet_quadratic), "jet_radius": float(args.jet_radius), "jet_cache": bool(args.jet_cache), "samples_per_jet": int(geometry.samples_per_jet), "jets_built": int(gstats.get("jets_built", 0)), "jet_cache_hits": int(gstats.get("jet_cache_hits", 0)), "oracle_samples_total": int(geometry.oracle_samples), "construction_complexity": "O(1) global; full path O(n) oracle samples/step; 411 cached coded path O(q) samples/step with q ~= directional_coded_factor*log2(n), fallback O(k) with k ~= sketch_factor*sqrt(n)"},
        "basis_cache": basis_cache_stats(),
        "parameters": {"system_source": str(args.system_source), "polys": str(args.polys or ""), "variables": str(args.variables or ""), "jet_quadratic": bool(args.jet_quadratic), "jet_radius": float(args.jet_radius), "jet_cache": bool(args.jet_cache), "geometric_polish_steps": int(args.geometric_polish_steps), "batch_wave": int(args.batch_wave), "batch_prepass": bool(args.batch_prepass), "batch_prepass_epochs": int(args.batch_prepass_epochs), "batch_kernel": BATCH_KERNEL, "resolved_batch_kernel": selected_batch_kernel_410(), "matrix_backend": MATRIX_BACKEND, "matrix_algorithm": MATRIX_ALGORITHM, "resolved_matrix_algorithm": selected_matrix_algorithm_410(n), "ns_iters": int(NS_ITERS), "ns_damping_floor": float(NS_DAMPING_FLOOR), "sketch_dim": int(SKETCH_DIM), "resolved_sketch_dim": int(sketch_dim_for(n)), "sketch_factor": float(SKETCH_FACTOR), "sketch_min_n": int(SKETCH_MIN_N), "sketch_mode": str(SKETCH_MODE), "sketch_solver": str(SKETCH_SOLVER), "sketch_seed": int(SKETCH_SEED), "directional_jet": bool(DIRECTIONAL_JET), "directional_jet_min_n": int(DIRECTIONAL_JET_MIN_N), "directional_jet_factor": float(DIRECTIONAL_JET_FACTOR), "directional_diff": str(DIRECTIONAL_DIFF_MODE), "directional_auto_central_cap": float(DIRECTIONAL_AUTO_CENTRAL_CAP), "directional_coded_probe": bool(DIRECTIONAL_CODED_PROBE), "directional_fast_projector": bool(DIRECTIONAL_FAST_PROJECTOR), "directional_fast_projector_cap": float(DIRECTIONAL_FAST_PROJECTOR_CAP), "directional_coded_factor": float(DIRECTIONAL_CODED_FACTOR), "directional_coded_min": int(DIRECTIONAL_CODED_MIN), "directional_coded_max": int(DIRECTIONAL_CODED_MAX), "directional_resolved_dim": int(max(1, min(n, int(math.ceil(float(DIRECTIONAL_JET_FACTOR) * math.sqrt(float(max(1, n)))))))), "directional_resolved_coded_dim": int(directional_coded_dim_for(n, max(1, min(n, int(math.ceil(float(DIRECTIONAL_JET_FACTOR) * math.sqrt(float(max(1, n))))))))), "torch_device": TORCH_DEVICE, "torch_resolved_device": TORCH_RESOLVED_DEVICE, "torch_complex_dtype": TORCH_COMPLEX_DTYPE, "torch_resolved_complex_dtype": TORCH_RESOLVED_DTYPE, "torch_real_dtype": TORCH_REAL_DTYPE, "torch_resolved_real_dtype": TORCH_RESOLVED_REAL_DTYPE, "torch_version": str(torch.__version__), "matrix_rcond": float(args.matrix_rcond), "matrix_condition_cap": float(args.matrix_condition_cap), "starts": str(args.starts or ""), "system_mode": str(args.system_mode), "count": int(args.count), "pool": int(args.pool), "accept": float(args.accept), "tol": float(args.tol), "epochs": int(args.epochs), "cluster_sep": float(args.cluster_sep), "line_search": int(args.line_search), "hypercube_nodes": int(args.hypercube_nodes), "atlas_selection": str(args.atlas_selection), "atlas_bypass_startopt": bool(args.atlas_bypass_startopt), "atlas_calls": int(atlas_calls), "startopt_steps": int(args.startopt_steps), "startopt_candidates": int(args.startopt_candidates), "startopt_gains_source": startopt_gains_source, "startopt_gains_count": len(gains), "scale_ladder": bool(args.scale_ladder), "ladder_subdiv": int(args.ladder_subdiv), "ladder_octaves": int(args.ladder_octaves), "ladder_base": float(args.ladder_base), "powers_count": len(powers), "irp_gains_source": irp_gains_source, "irp_gains_count": len(irp_gain_values), "irp_scales_count": len(irp_scales)},
        "roots": encoded_roots,
        "trials": trials if bool(args.verbose_trials) else trials[: min(len(trials), int(args.keep_trials))],
        "summary": {"requested_roots": int(args.count), "unique_roots": len(roots), "success": bool(len(roots) >= int(args.count)), "trials_used": len(trials), "duplicates": int(duplicates), "failures": int(failures), "generation_seconds": float(geometry.generation_seconds), "extract_seconds": extract_seconds, "total_seconds": float(now() - t_case), "oracle_samples_total": int(geometry.oracle_samples), "oracle_samples_per_trial": samples_per_trial, "eval_stats": gstats},
    }
    return result


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="411 standalone Pandrosion IRP cached directional-sketch local-jet matrix engine.")
    p.add_argument("--cases", default="2,4")
    p.add_argument("--seed-index", type=int, default=0)
    p.add_argument("--equation-normalize", action="store_true", default=False)
    p.add_argument("--no-equation-normalize", dest="equation_normalize", action="store_false")
    p.add_argument("--system-mode", choices=["auto", "dense", "geometry-kernel", "geometry", "kernel", "projective-kernel", "lazy-feature", "lazy", "feature", "stream"], default="auto")
    p.add_argument("--dense-max-terms", type=int, default=250000)
    # KS-base geometry-kernel backend knobs (only used when the base itself is a large KS system).
    p.add_argument("--geometry-anchors", type=int, default=0)
    p.add_argument("--geometry-anchor-cap", type=int, default=4106)
    p.add_argument("--geometry-anchor-scales", default="0.25,0.5,1,2,4")
    p.add_argument("--geometry-dynamic-normalize", action="store_true", default=True)
    p.add_argument("--no-geometry-dynamic-normalize", dest="geometry_dynamic_normalize", action="store_false")
    p.add_argument("--geometry-self-normalize", action="store_true", default=True)
    p.add_argument("--no-geometry-self-normalize", dest="geometry_self_normalize", action="store_false")
    p.add_argument("--geometry-eval-block", type=int, default=128)
    p.add_argument("--lazy-features", type=int, default=0)
    p.add_argument("--lazy-feature-cap", type=int, default=8192)
    p.add_argument("--lazy-projective-normalize", action="store_true", default=True)
    p.add_argument("--no-lazy-projective-normalize", dest="lazy_projective_normalize", action="store_false")
    p.add_argument("--lazy-dynamic-normalize", action="store_true", default=True)
    p.add_argument("--no-lazy-dynamic-normalize", dest="lazy_dynamic_normalize", action="store_false")
    p.add_argument("--lazy-eval-block", type=int, default=128)
    # The local-jet geometry inherited from 317.
    p.add_argument("--jet-quadratic", action="store_true", default=True, help="Also sample local curvature (second-order jet).")
    p.add_argument("--no-jet-quadratic", dest="jet_quadratic", action="store_false")
    p.add_argument("--jet-radius", type=float, default=1e-5, help="Relative hypercube radius for local-jet sampling.")
    p.add_argument("--jet-cache", action="store_true", default=True, help="Reuse local jets across nearby queries.")
    p.add_argument("--no-jet-cache", dest="jet_cache", action="store_false")
    p.add_argument("--jet-cache-decimals", type=int, default=9)
    p.add_argument("--geometric-polish-steps", type=int, default=6)
    p.add_argument("--batch-wave", type=int, default=8, help="Number of starts prepared in a wave before individual Pandrosion IRP correction.")
    p.add_argument("--batch-prepass", action="store_true", default=True, help="Run a batched GEMM direct prepass before the normal Pandrosion IRP corrector.")
    p.add_argument("--no-batch-prepass", dest="batch_prepass", action="store_false")
    p.add_argument("--batch-prepass-epochs", type=int, default=1, help="Batched direct GEMM epochs before handing each start to Pandrosion IRP.")
    p.add_argument("--batch-kernel", choices=["auto", "sketch", "realpack", "complex"], default="auto", help="Batched prepass kernel. auto follows the selected 411 matrix algorithm.")
    p.add_argument("--origin-seed", action="store_true", default=True)
    p.add_argument("--no-origin-seed", dest="origin_seed", action="store_false")
    p.add_argument("--origin-seed-h", type=float, default=1e-5)
    p.add_argument("--origin-seed-max-norm", type=float, default=0.0)
    p.add_argument("--origin-seed-period", type=int, default=0)
    p.add_argument("--linear-scale", type=float, default=1.0)
    p.add_argument("--count", type=int, default=8)
    p.add_argument("--pool", type=int, default=4106)
    p.add_argument("--epochs", type=int, default=24)
    p.add_argument("--tol", type=float, default=1e-12)
    p.add_argument("--accept", "--residual-accept", type=float, default=1e-8)
    p.add_argument("--cluster-sep", type=float, default=1e-8)
    p.add_argument("--trial-timeout", type=float, default=0.0)
    p.add_argument("--line-search", type=int, default=12)
    p.add_argument("--hypercube-nodes", type=int, default=0)
    p.add_argument("--irp-layers", type=int, default=2)
    p.add_argument("--irp-inner-epochs", type=int, default=2)
    p.add_argument("--local-inner-epochs", type=int, default=3)
    p.add_argument("--lazy-direct-epochs", type=int, default=1)
    p.add_argument("--lazy-trigger-drop", type=float, default=0.82)
    p.add_argument("--lazy-trigger-after", type=int, default=1)
    p.add_argument("--lazy-bad-cond", type=float, default=1e10)
    p.add_argument("--lazy-log-energy", type=float, default=8.0)
    p.add_argument("--eager-irp", action="store_true", default=False)
    p.add_argument("--rescue-collapsed", action="store_true", default=False)
    # Regularised matrix-jet inversion (tames the giant near-singular step in high n).
    p.add_argument("--matrix-backend", choices=["torch", "numpy", "auto"], default="torch", help="Matrix kernel backend. torch+adaptive-directional-sketch is the 411 default; numpy preserves the SVD fallback behavior.")
    p.add_argument("--matrix-algorithm", choices=["adaptive-directional-sketch", "directional-sketch", "adaptive-sketch", "sketch-ns", "adaptive-ns", "realpack-ns", "gemm-ns", "svd", "auto"], default="adaptive-directional-sketch", help="Matrix correction algorithm. adaptive-directional-sketch samples JP directly for large local jets and uses full GEMM for small systems.")
    p.add_argument("--ns-iters", type=int, default=8, help="Newton-Schulz iterations for GEMM/Newton-Schulz algorithms.")
    p.add_argument("--ns-damping-floor", type=float, default=1e-5, help="Minimum dimensionless normal-equation damping for GEMM/Newton-Schulz.")
    p.add_argument("--torch-device", choices=["auto", "cpu", "mps", "cuda"], default="auto", help="Torch device for matrix kernels. auto uses CUDA if available, then MPS, then CPU.")
    p.add_argument("--torch-complex-dtype", choices=["auto", "complex64", "complex128"], default="auto", help="Torch complex dtype. auto uses complex64 on MPS and complex128 elsewhere.")
    p.add_argument("--torch-real-dtype", choices=["auto", "float16", "bfloat16", "float32", "float64"], default="auto", help="Torch real dtype for 411 real-packed kernels. auto uses float32.")
    p.add_argument("--sketch-dim", type=int, default=0, help="Subspace dimension for --matrix-algorithm sketch-ns (0 = sqrt(n)-scaled auto).")
    p.add_argument("--sketch-factor", type=float, default=2.75, help="Auto sketch dimension multiplier: k ~= sketch_factor*sqrt(n).")
    p.add_argument("--sketch-min-n", type=int, default=96, help="Minimum dimension before adaptive-sketch switches away from full GEMM.")
    p.add_argument("--sketch-seed", type=int, default=411, help="Deterministic seed for the randomized sketch basis.")
    p.add_argument("--sketch-mode", choices=["auto", "rademacher", "cached-rademacher", "coordinate", "sparse-sign"], default="cached-rademacher", help="Sketch basis type. cached-rademacher is safer; sparse-sign is faster and more aggressive.")
    p.add_argument("--sketch-basis-cache", action="store_true", default=True, help="Cache deterministic P/R sketch bases to reduce Python/NumPy setup overhead.")
    p.add_argument("--no-sketch-basis-cache", dest="sketch_basis_cache", action="store_false")
    p.add_argument("--sketch-basis-cache-max", type=int, default=128, help="Maximum cached sketch/coded bases before the cache is cleared.")
    p.add_argument("--sketch-solver", choices=["svd", "ns", "auto"], default="svd", help="Reduced subspace solver. svd is more accurate; ns stays GEMM/Newton-Schulz inside the sketch.")
    p.add_argument("--directional-jet", action="store_true", default=True, help="For sketch-sized systems, sample J P directly from the oracle instead of building a full local jet.")
    p.add_argument("--no-directional-jet", dest="directional_jet", action="store_false")
    p.add_argument("--directional-jet-min-n", type=int, default=96, help="Minimum dimension before directional JP sampling activates.")
    p.add_argument("--directional-jet-factor", type=float, default=2.75, help="Directional JP dimension multiplier: k ~= factor*sqrt(n).")
    p.add_argument("--directional-diff", choices=["auto", "forward", "central"], default="auto", help="Directional JP finite difference. forward uses k oracle samples; central uses 2k; auto tries forward and falls back to central on poor reduced residual.")
    p.add_argument("--directional-auto-central-cap", type=float, default=0.35, help="In --directional-diff auto, switch from forward to central if reduced relative residual is above this cap.")
    p.add_argument("--directional-coded-probe", action="store_true", default=True, help="Try the 411 coded probe D=P R with q<<k before the full directional sketch.")
    p.add_argument("--no-directional-coded-probe", dest="directional_coded_probe", action="store_false")
    p.add_argument("--directional-coded-factor", type=float, default=2.0, help="Coded probe multiplier: q ~= factor*log2(n).")
    p.add_argument("--directional-coded-min", type=int, default=8, help="Minimum coded probe dimension q.")
    p.add_argument("--directional-coded-max", type=int, default=0, help="Maximum coded probe dimension q (0 = no cap).")
    p.add_argument("--directional-fast-projector", action="store_true", default=True, help="Try the 411 diagonal projected reduced solver before SVD/Newton-Schulz.")
    p.add_argument("--no-directional-fast-projector", dest="directional_fast_projector", action="store_false")
    p.add_argument("--directional-fast-projector-cap", type=float, default=0.05, help="Accept the fast projected reduced solver only below this relative residual.")
    p.add_argument("--directional-basis-reuse", action="store_true", default=True, help="Reuse deterministic directional basis across corrector epochs to amortize setup overhead.")
    p.add_argument("--no-directional-basis-reuse", dest="directional_basis_reuse", action="store_false")
    p.add_argument("--matrix-rcond", type=float, default=1e-12, help="Relative SVD cutoff for matrix pseudo-inverse corrections.")
    p.add_argument("--matrix-condition-cap", type=float, default=1e12, help="Drop singular modes beyond this effective condition number.")
    p.add_argument("--lm-damping", type=float, default=0.0, help="Dimensionless LM/Tikhonov damping for matrix-jet inversion (0 = off, try 1e-2).")
    p.add_argument("--trust-radius", type=float, default=0.0, help="Cap the matrix-jet step to trust_radius*||y|| (0 = off, try 1.0).")
    p.add_argument("--trust-region", action="store_true", default=False, help="Convenience: enable LM damping + trust radius with sensible defaults.")
    p.add_argument("--irp-chart-top", type=int, default=2)
    p.add_argument("--irp-gains", default=None, help="Homothety gains for the IRP chart palette. Defaults to the Pandrosion ladder when --scale-ladder is enabled, otherwise legacy dyadic gains.")
    p.add_argument("--irp-phases", default="0,0.08,-0.08,0.19,-0.19")
    p.add_argument("--irp-top", type=int, default=14)
    p.add_argument("--irp-inversion", action="store_true", default=True)
    p.add_argument("--no-irp-inversion", dest="irp_inversion", action="store_false")
    p.add_argument("--collapse", action="store_true", default=True)
    p.add_argument("--no-collapse", dest="collapse", action="store_false")
    p.add_argument("--collapse-residual", type=float, default=1e-4)
    p.add_argument("--collapse-drop", type=float, default=0.42)
    p.add_argument("--collapse-rel-step", type=float, default=0.35)
    p.add_argument("--collapse-after", type=int, default=2)
    p.add_argument("--line-grid", default="1,0.75,0.5,0.35,0.25,0.18,0.125,0.09,0.0625,0.045,0.03125,0.02")
    p.add_argument("--powers", default=None)
    p.add_argument("--power-cap", type=float, default=1048576.0)
    # The 'complete logarithmic ladder' (Pandrosion proportional-means x^{k/p} scaling).
    p.add_argument("--scale-ladder", dest="scale_ladder", action="store_true", default=True, help="Use an equally-spaced logarithmic (proportional-means) ladder for the start homothety palette and IRP charts.")
    p.add_argument("--no-scale-ladder", dest="scale_ladder", action="store_false", help="Use the legacy dyadic/decadic start palette and dyadic IRP chart gains.")
    p.add_argument("--ladder-subdiv", type=int, default=3, help="Geometric means per octave (p): 3 inserts the cube-root means base^{k/3}; the x^{k/p} ladder.")
    p.add_argument("--ladder-octaves", type=int, default=12, help="Half-span of the start ladder in octaves of --ladder-base.")
    p.add_argument("--ladder-base", type=float, default=2.0, help="Octave base of the logarithmic ladder (>1).")
    p.add_argument("--angles", default=None)
    p.add_argument("--rays", default=None)
    p.add_argument("--startopt-steps", type=int, default=1)
    p.add_argument("--startopt-candidates", type=int, default=12)
    p.add_argument("--startopt-gains", default=None)
    p.add_argument("--startopt-micro-epochs", type=int, default=0)
    p.add_argument("--universal-cells", type=int, default=16)
    p.add_argument("--universal-shells", type=int, default=5)
    p.add_argument("--atlas-selection", choices=["diverse-shell", "compact-score"], default="diverse-shell", help="Automatic start selection. diverse-shell cycles deterministic logarithmic shells; compact-score keeps the older residual-min atlas cell selection.")
    p.add_argument("--atlas-bypass-startopt", action="store_true", default=True, help="Do not run StartOpt on diverse-shell atlas starts, preserving basin diversity for all-root extraction.")
    p.add_argument("--no-atlas-bypass-startopt", dest="atlas_bypass_startopt", action="store_false")
    p.add_argument("--cell-probe-radius", type=float, default=0.14)
    p.add_argument("--cell-descent-min", type=float, default=1.02)
    p.add_argument("--cell-equal-gap-min", type=float, default=1e-10)
    p.add_argument("--cell-log-max", type=float, default=80.0)
    p.add_argument("--universal-cycle", action="store_true", default=True)
    p.add_argument("--no-universal-cycle", dest="universal_cycle", action="store_false")
    p.add_argument("--out", default=None)
    p.add_argument("--outdir", default="/private/tmp/411_standalone_cached_pandrosion_irp_gemm_out")
    p.add_argument("--keep-trials", type=int, default=160)
    p.add_argument("--verbose-trials", action="store_true")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--system-source", choices=["ks", "kostlan", "polynomial", "poly", "expr", "expression"], default="ks")
    p.add_argument("--polys", "--poly", default=None)
    p.add_argument("--variables", default=None)
    p.add_argument("--starts", default=None)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if bool(args.trust_region):
        if float(args.lm_damping) <= 0.0:
            args.lm_damping = 1e-2
        if float(args.trust_radius) <= 0.0:
            args.trust_radius = 1.0
    if bool(args.self_test):
        args.system_source = "polynomial"
        args.system_mode = "auto"
        args.cases = "1,2"
        args.polys = "x^2 - 3*x - 10"
        args.starts = "-8,4"
        args.count = 2
        args.pool = min(int(args.pool), 64)
        args.epochs = min(int(args.epochs), 24)
        args.accept = min(float(args.accept), 1e-8)
        args.out = args.out or "/private/tmp/411_standalone_cached_pandrosion_irp_gemm_out/self_test_411.json"
    configure_matrix_backend(args)
    cases = [c.strip() for c in str(args.cases).replace("|", ";").split(";") if c.strip()]
    outputs = [run_case(args, c) for c in cases]
    final = outputs[0] if len(outputs) == 1 else {"script": Path(__file__).name, "autonomous": True, "cases": outputs}
    out = Path(args.out) if args.out else Path(args.outdir) / f"411_cached_pandrosion_irp_gemm_{cases[0].replace(',', 'x') if cases else 'case'}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(final, indent=2), encoding="utf-8")
    print("=" * 120, flush=True)
    print("411 STANDALONE PANDROSION IRP -- CACHED GENERAL CODED MATRIX-JET CALCULUS", flush=True)
    print("PyTorch + NumPy + stdlib only; no local flow imports. The geometry is 317's lazily-sampled local-jet field;", flush=True)
    first_n = parse_case(cases[0])[0] if cases else 0
    directional_k = int(max(1, min(max(1, first_n), int(math.ceil(float(DIRECTIONAL_JET_FACTOR) * math.sqrt(float(max(1, first_n)))))))) if first_n else 0
    directional_q = directional_coded_dim_for(first_n, directional_k) if first_n else 0
    print(f"matrix_backend={MATRIX_BACKEND}, algorithm={MATRIX_ALGORITHM}->{selected_matrix_algorithm_410(first_n)}, torch={torch.__version__}, device={TORCH_RESOLVED_DEVICE}, complex_dtype={TORCH_RESOLVED_DTYPE}, real_dtype={TORCH_RESOLVED_REAL_DTYPE}, batch_kernel={selected_batch_kernel_410()}, sketch_dim={sketch_dim_for(first_n) if first_n else 0}, sketch_mode={SKETCH_MODE}, basis_cache={bool(SKETCH_BASIS_CACHE)}, basis_reuse={bool(DIRECTIONAL_BASIS_REUSE)}, directional_jet={bool(DIRECTIONAL_JET)}, directional_diff={DIRECTIONAL_DIFF_MODE}, coded_probe={bool(DIRECTIONAL_CODED_PROBE)}, fast_projector={bool(DIRECTIONAL_FAST_PROJECTOR)}, directional_q={directional_q}, directional_k={directional_k}, batch_wave={int(args.batch_wave)}, irp=on", flush=True)
    print("=" * 120, flush=True)
    for r in outputs:
        s = r["summary"]
        g = r.get("geometry", {})
        print(f"case={r.get('family')}({r['n']},{r['degree']}), base_backend={r.get('base_backend')}, geometry=local-jet-field (quad={g.get('use_quadratic')}, samples/jet={g.get('samples_per_jet')}), seed={r['seed']}", flush=True)
        print(f"roots={s['unique_roots']}/{s['requested_roots']} success={s['success']} trials={s['trials_used']} duplicates={s['duplicates']} failures={s['failures']}", flush=True)
        print(f"seconds: extract={s['extract_seconds']:.2f}, total={s['total_seconds']:.2f}; oracle_samples={s['oracle_samples_total']} ({s['oracle_samples_per_trial']:.0f}/trial)", flush=True)
        if r.get("roots"):
            best = r["roots"][0]
            print(f"best_root: residual={float(best.get('geometric_residual', float('inf'))):.3e} (faithful), trial={best.get('trial')}", flush=True)
    print(f"out={out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
