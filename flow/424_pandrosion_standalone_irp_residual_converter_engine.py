#!/usr/bin/env python3
"""424 Pandrosion IRP converter: A ~= A_irp + residual.

This script receives a dense matrix A, tries to convert it into an exploitable
Pandrosion IRP certificate, measures the residual, and applies the certificate
to compute an approximate A @ B.

Model:

    A ~= A_irp + R_sparse + R_tail

where A_irp is diagonal + low-rank and R_sparse stores the largest residual
entries.  The timed 424 application computes

    Y_424 = A_irp @ B + R_sparse @ B

and reports the true error against torch.matmul(A, B).  The dense tail R_tail
is measured but not applied, because applying it exactly would reintroduce the
dense matmul cost.  This is the honest tradeoff: if the tail is small, 424 can
beat torch.matmul with controlled error; if the tail is large, the conversion
fails to produce an exploitable certificate.

This is not a universal dense-GEMM replacement.  It is an automatic conversion
experiment that decides whether a dense A has enough exploitable IRP structure.

Dependencies: Python stdlib + NumPy + PyTorch.  No local imports.
"""
from __future__ import annotations

import argparse
import dataclasses
import gc
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np

try:
    import torch
except Exception as exc:  # pragma: no cover
    torch = None  # type: ignore[assignment]
    TORCH_IMPORT_ERROR = exc
else:
    TORCH_IMPORT_ERROR = None


OUT = Path("/private/tmp/424_irp_residual_converter_vs_torch_matmul_benchmark.json")


@dataclasses.dataclass(frozen=True)
class Config:
    dims: tuple[int, ...] = (512, 1024, 2048)
    families: tuple[str, ...] = ("diag_lowrank", "diag_lowrank_sparse", "noisy_lowrank", "banded", "random_dense")
    columns: int = 0
    rank: int = 8
    oversample: int = 6
    alternating_steps: int = 3
    sparse_per_row: int = 4
    input_sparse_per_row: int = 4
    input_rank: int = 8
    bandwidth: int = 3
    noise: float = 0.02
    dtype: str = "float32"
    seed: int = 424
    repeats: int = 0
    warmup: int = 3
    out: Path = OUT


@dataclasses.dataclass
class IRPResidualCertificate:
    n: int
    diag: Any | None
    u: Any | None
    v: Any | None
    sparse_rows: Any | None
    sparse_cols: Any | None
    sparse_vals: Any | None
    rank: int
    sparse_nnz: int
    input_family: str
    conversion_ms: float
    fro_a: float
    fro_after_diag: float
    fro_after_lowrank: float
    fro_after_sparse: float
    tail_relative_fro: float
    certificate_size: int


def sync(device: Any) -> None:
    if torch is None:
        return
    if device.type == "cuda":
        torch.cuda.synchronize()
    elif device.type == "mps":
        torch.mps.synchronize()


def empty_cache(device: Any) -> None:
    gc.collect()
    if torch is None:
        return
    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "mps" and hasattr(torch.mps, "empty_cache"):
        torch.mps.empty_cache()


def pick_device() -> Any:
    if torch is None:
        raise RuntimeError(f"PyTorch import failed: {TORCH_IMPORT_ERROR!r}")
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def pick_dtype(name: str) -> Any:
    if torch is None:
        raise RuntimeError("PyTorch is not available")
    key = str(name).lower().replace("torch.", "")
    if key in {"float16", "fp16", "half"}:
        return torch.float16
    if key in {"bfloat16", "bf16"}:
        return torch.bfloat16
    if key in {"float32", "fp32", "single"}:
        return torch.float32
    raise ValueError(f"unsupported dtype {name!r}")


def torch_randn(shape: tuple[int, ...], *, seed: int, device: Any, dtype: Any, scale: float = 1.0) -> Any:
    gen = torch.Generator(device="cpu")
    gen.manual_seed(int(seed))
    x = torch.randn(shape, generator=gen, dtype=torch.float32) * float(scale)
    return x.to(device=device, dtype=dtype)


def torch_randint(high: int, shape: tuple[int, ...], *, seed: int, device: Any) -> Any:
    gen = torch.Generator(device="cpu")
    gen.manual_seed(int(seed))
    return torch.randint(int(high), shape, generator=gen, dtype=torch.long).to(device=device)


def stats(vals: list[float]) -> dict[str, float]:
    a = np.asarray(vals, dtype=float)
    return {
        "min_ms": float(np.min(a)),
        "median_ms": float(np.median(a)),
        "mean_ms": float(np.mean(a)),
        "std_ms": float(np.std(a, ddof=1)) if a.size > 1 else 0.0,
    }


def bench(fn: Any, *, device: Any, repeats: int, warmup: int) -> tuple[dict[str, float], Any]:
    last = None
    with torch.no_grad():
        for _ in range(max(0, int(warmup))):
            last = fn()
        sync(device)
        vals: list[float] = []
        for _ in range(max(1, int(repeats))):
            sync(device)
            t0 = time.perf_counter()
            last = fn()
            sync(device)
            vals.append(1e3 * (time.perf_counter() - t0))
    return stats(vals), last


def default_repeats(n: int) -> int:
    if n <= 512:
        return 18
    if n <= 1024:
        return 12
    return 6


def family_id(family: str) -> int:
    return {
        "diag_lowrank": 1,
        "diag_lowrank_sparse": 2,
        "noisy_lowrank": 3,
        "banded": 4,
        "random_dense": 5,
    }.get(str(family), 99)


def add_sparse_entries(a: Any, *, per_row: int, seed: int, scale: float) -> Any:
    n = int(a.shape[0])
    k = max(0, int(per_row))
    if k <= 0:
        return a
    rows = torch.arange(n, device=a.device, dtype=torch.long).repeat_interleave(k)
    cols = torch_randint(n, (n * k,), seed=seed, device=a.device)
    vals = torch_randn((n * k,), seed=seed + 1, device=a.device, dtype=a.dtype, scale=scale)
    a.index_put_((rows, cols), vals, accumulate=True)
    return a


def make_dense_a(family: str, n: int, cfg: Config, device: Any, dtype: Any, seed: int) -> Any:
    n = int(n)
    fam = str(family)
    diag = torch_randn((n,), seed=seed + 11, device=device, dtype=dtype, scale=0.35)
    idx = torch.arange(n, device=device)
    a = torch.zeros((n, n), device=device, dtype=dtype)

    if fam in {"diag_lowrank", "diag_lowrank_sparse", "noisy_lowrank"}:
        r = max(1, int(cfg.input_rank))
        u = torch_randn((n, r), seed=seed + 23, device=device, dtype=dtype, scale=1.0 / math.sqrt(max(1, r)))
        v = torch_randn((n, r), seed=seed + 29, device=device, dtype=dtype, scale=1.0 / math.sqrt(max(1, n)))
        a[idx, idx] = diag
        a += torch.mm(u, v.transpose(0, 1))
        if fam == "diag_lowrank_sparse":
            a = add_sparse_entries(a, per_row=int(cfg.input_sparse_per_row), seed=seed + 101, scale=0.12)
        if fam == "noisy_lowrank":
            a += torch_randn((n, n), seed=seed + 103, device=device, dtype=dtype, scale=float(cfg.noise) / math.sqrt(max(1, n)))
        return a

    if fam == "banded":
        w = max(0, int(cfg.bandwidth))
        for offset in range(-w, w + 1):
            length = n - abs(offset)
            vals = torch_randn((length,), seed=seed + 211 + offset + 3 * w, device=device, dtype=dtype, scale=0.2)
            if offset >= 0:
                rows = torch.arange(0, n - offset, device=device)
                cols = rows + offset
            else:
                cols = torch.arange(0, n + offset, device=device)
                rows = cols - offset
            a[rows, cols] = vals
        return a

    if fam == "random_dense":
        return torch_randn((n, n), seed=seed + 307, device=device, dtype=dtype, scale=1.0 / math.sqrt(max(1, n)))

    raise ValueError(f"unknown family {family!r}")


def fro_norm(x: Any) -> float:
    return float(torch.linalg.vector_norm(x.reshape(-1).to(torch.float32)).item())


def randomized_lowrank(residual: Any, rank: int, oversample: int, seed: int) -> tuple[Any, Any, str]:
    n = int(residual.shape[0])
    r = max(0, min(int(rank), n))
    if r <= 0:
        return None, None, "rank-zero"
    q = min(n, max(r, r + max(0, int(oversample))))
    device = residual.device
    dtype = residual.dtype
    omega = torch_randn((n, q), seed=seed, device=device, dtype=dtype, scale=1.0 / math.sqrt(max(1, q)))
    try:
        y = torch.mm(residual, omega)
        qmat, _ = torch.linalg.qr(y, mode="reduced")
        small = torch.mm(qmat.transpose(0, 1), residual)
        uh, s, vh = torch.linalg.svd(small, full_matrices=False)
        uu = torch.mm(qmat, uh[:, :r]) * s[:r].unsqueeze(0)
        vv = vh[:r, :].transpose(0, 1).contiguous()
        return uu.to(device=device, dtype=dtype), vv.to(device=device, dtype=dtype), "device-randomized-svd"
    except RuntimeError:
        cpu = torch.device("cpu")
        rc = residual.detach().to(cpu, dtype=torch.float32)
        gen = torch.Generator(device="cpu")
        gen.manual_seed(int(seed))
        omega_c = torch.randn((n, q), generator=gen, dtype=torch.float32) / math.sqrt(max(1, q))
        y = torch.mm(rc, omega_c)
        qmat, _ = torch.linalg.qr(y, mode="reduced")
        small = torch.mm(qmat.transpose(0, 1), rc)
        uh, s, vh = torch.linalg.svd(small, full_matrices=False)
        uu = torch.mm(qmat, uh[:, :r]) * s[:r].unsqueeze(0)
        vv = vh[:r, :].transpose(0, 1).contiguous()
        return uu.to(device=device, dtype=dtype), vv.to(device=device, dtype=dtype), "cpu-randomized-svd"


def sparse_topk(residual: Any, nnz: int) -> tuple[Any | None, Any | None, Any | None, float]:
    n = int(residual.shape[0])
    k = max(0, min(int(nnz), n * n))
    if k <= 0:
        return None, None, None, fro_norm(residual)
    flat = residual.reshape(-1)
    vals_abs, inds = torch.topk(torch.abs(flat), k=k, largest=True, sorted=False)
    vals = flat.index_select(0, inds)
    rows = torch.div(inds, n, rounding_mode="floor").to(torch.long)
    cols = (inds - rows * n).to(torch.long)
    total_sq = float(torch.sum(torch.abs(flat.to(torch.float32)) ** 2).item())
    kept_sq = float(torch.sum(vals_abs.to(torch.float32) ** 2).item())
    tail = math.sqrt(max(0.0, total_sq - kept_sq))
    return rows, cols, vals, float(tail)


def convert_to_irp_residual(a: Any, family: str, cfg: Config, seed: int, device: Any) -> IRPResidualCertificate:
    sync(device)
    t0 = time.perf_counter()
    n = int(a.shape[0])
    idx = torch.arange(n, device=a.device)
    fro_a = fro_norm(a)
    diag = torch.diagonal(a).clone()
    residual_after_diag = a.clone()
    residual_after_diag[idx, idx] -= diag
    fro_after_diag = fro_norm(residual_after_diag)

    u = None
    v = None
    lowrank_solver = "rank-zero"
    steps = max(1, int(cfg.alternating_steps))
    residual = residual_after_diag
    # Alternating diagonal/low-rank fit:
    #   1. remove current diagonal estimate;
    #   2. fit a low-rank IRP factor to the remaining dense part;
    #   3. re-estimate the diagonal from the residual.
    # This handles the common ambiguity where the diagonal of a low-rank term is
    # initially absorbed into diag(A).
    for it in range(steps):
        base = a.clone()
        base[idx, idx] -= diag
        u, v, lowrank_solver = randomized_lowrank(base, int(cfg.rank), int(cfg.oversample), seed + 401 + 17 * it)
        residual = a.clone()
        if u is not None and v is not None:
            residual = residual - torch.mm(u, v.transpose(0, 1))
        diag = torch.diagonal(residual).clone()
    residual = a.clone()
    residual[idx, idx] -= diag
    if u is not None and v is not None:
        residual = residual - torch.mm(u, v.transpose(0, 1))
    fro_after_lowrank = fro_norm(residual)

    rows, cols, vals, tail_fro = sparse_topk(residual, int(cfg.sparse_per_row) * n)
    fro_after_sparse = float(tail_fro)
    cert_size = int(n)
    if u is not None and v is not None:
        cert_size += int(u.numel() + v.numel())
    sparse_nnz = 0 if vals is None else int(vals.numel())
    cert_size += 3 * sparse_nnz
    sync(device)
    conversion_ms = 1e3 * (time.perf_counter() - t0)

    cert = IRPResidualCertificate(
        n=n,
        diag=diag,
        u=u,
        v=v,
        sparse_rows=rows,
        sparse_cols=cols,
        sparse_vals=vals,
        rank=0 if u is None else int(u.shape[1]),
        sparse_nnz=sparse_nnz,
        input_family=str(family),
        conversion_ms=float(conversion_ms),
        fro_a=float(fro_a),
        fro_after_diag=float(fro_after_diag),
        fro_after_lowrank=float(fro_after_lowrank),
        fro_after_sparse=float(fro_after_sparse),
        tail_relative_fro=float(fro_after_sparse / max(1e-30, fro_a)),
        certificate_size=int(cert_size),
    )
    setattr(cert, "lowrank_solver", lowrank_solver)
    return cert


def apply_sparse_residual(rows: Any | None, cols: Any | None, vals: Any | None, b: Any) -> Any:
    out = torch.zeros((int(b.shape[0]), int(b.shape[1])), device=b.device, dtype=b.dtype)
    if rows is None or cols is None or vals is None or int(vals.numel()) == 0:
        return out
    # Chunking limits temporary memory for nnz x columns.
    nnz = int(vals.numel())
    chunk = max(1, min(nnz, 4096))
    for lo in range(0, nnz, chunk):
        hi = min(nnz, lo + chunk)
        contrib = vals[lo:hi].unsqueeze(1) * b.index_select(0, cols[lo:hi])
        out.index_add_(0, rows[lo:hi], contrib)
    return out


def pandrosion_424_apply(cert: IRPResidualCertificate, b: Any, *, include_sparse: bool = True) -> Any:
    out = torch.zeros((int(cert.n), int(b.shape[1])), device=b.device, dtype=b.dtype)
    if cert.diag is not None:
        out = out + cert.diag[:, None] * b
    if cert.u is not None and cert.v is not None:
        out = out + torch.mm(cert.u, torch.mm(cert.v.transpose(0, 1), b))
    if bool(include_sparse):
        out = out + apply_sparse_residual(cert.sparse_rows, cert.sparse_cols, cert.sparse_vals, b)
    return out


def error_metrics(y: Any, ref: Any) -> dict[str, float]:
    diff = (y - ref).detach()
    ref32 = ref.detach().reshape(-1).to(torch.float32)
    diff32 = diff.reshape(-1).to(torch.float32)
    denom = torch.linalg.vector_norm(ref32).item()
    num = torch.linalg.vector_norm(diff32).item()
    return {
        "relative_l2": float(num / max(1e-30, denom)),
        "max_abs": float(torch.max(torch.abs(diff32)).item()),
    }


def run_case(family: str, n: int, cfg: Config, device: Any, dtype: Any) -> dict[str, Any]:
    n = int(n)
    m = int(cfg.columns) if int(cfg.columns) > 0 else n
    seed = int(cfg.seed) + 1009 * n + 9173 * family_id(family)
    a = make_dense_a(family, n, cfg, device, dtype, seed)
    b = torch_randn((n, m), seed=seed + 701, device=device, dtype=dtype, scale=0.5)
    cert = convert_to_irp_residual(a, family, cfg, seed, device)
    repeats = int(cfg.repeats) if int(cfg.repeats) > 0 else default_repeats(n)

    with torch.no_grad():
        ref = torch.matmul(a, b)
        got_sparse = pandrosion_424_apply(cert, b, include_sparse=True)
        got_irp = pandrosion_424_apply(cert, b, include_sparse=False)
        err_sparse = error_metrics(got_sparse, ref)
        err_irp = error_metrics(got_irp, ref)

    torch_stats, _ = bench(lambda: torch.matmul(a, b), device=device, repeats=repeats, warmup=int(cfg.warmup))
    sparse_stats, _ = bench(lambda: pandrosion_424_apply(cert, b, include_sparse=True), device=device, repeats=repeats, warmup=int(cfg.warmup))
    irp_stats, _ = bench(lambda: pandrosion_424_apply(cert, b, include_sparse=False), device=device, repeats=repeats, warmup=int(cfg.warmup))

    torch_ms = float(torch_stats["min_ms"])
    sparse_ms = float(sparse_stats["min_ms"])
    irp_ms = float(irp_stats["min_ms"])
    speed_sparse = torch_ms / max(1e-30, sparse_ms)
    speed_irp = torch_ms / max(1e-30, irp_ms)
    useful_tol = 1e-3

    row = {
        "family": str(family),
        "n": n,
        "columns": m,
        "device": str(device),
        "dtype": str(dtype).replace("torch.", ""),
        "torch_matmul_min_ms": torch_ms,
        "torch_matmul_median_ms": float(torch_stats["median_ms"]),
        "pandrosion_424_irp_only_min_ms": irp_ms,
        "pandrosion_424_irp_sparse_min_ms": sparse_ms,
        "speedup_irp_only_vs_torch_matmul": float(speed_irp),
        "speedup_irp_sparse_vs_torch_matmul": float(speed_sparse),
        "beats_torch_matmul_irp_only": bool(speed_irp > 1.0),
        "beats_torch_matmul_irp_sparse": bool(speed_sparse > 1.0),
        "relative_l2_error_irp_only": float(err_irp["relative_l2"]),
        "relative_l2_error_irp_sparse": float(err_sparse["relative_l2"]),
        "max_abs_error_irp_sparse": float(err_sparse["max_abs"]),
        "useful_win_irp_sparse_tol_1e_3": bool(speed_sparse > 1.0 and float(err_sparse["relative_l2"]) <= useful_tol),
        "conversion_ms": float(cert.conversion_ms),
        "conversion_amortized_calls_to_break_even": float(
            cert.conversion_ms / max(1e-30, torch_ms - sparse_ms)
        ) if torch_ms > sparse_ms else float("inf"),
        "certificate_size": int(cert.certificate_size),
        "dense_entries": int(n * n),
        "certificate_compression": float((n * n) / max(1, int(cert.certificate_size))),
        "rank": int(cert.rank),
        "sparse_nnz": int(cert.sparse_nnz),
        "sparse_per_row": float(cert.sparse_nnz / max(1, n)),
        "fro_a": float(cert.fro_a),
        "fro_after_diag_rel": float(cert.fro_after_diag / max(1e-30, cert.fro_a)),
        "fro_after_lowrank_rel": float(cert.fro_after_lowrank / max(1e-30, cert.fro_a)),
        "fro_after_sparse_rel": float(cert.fro_after_sparse / max(1e-30, cert.fro_a)),
        "tail_relative_fro": float(cert.tail_relative_fro),
        "lowrank_solver": str(getattr(cert, "lowrank_solver", "")),
        "engine_dense_tail_applied": False,
        "engine_dense_matmul_calls": 0,
        "interpretation": "A was converted to A_irp + R_sparse + R_tail; timed engine applies A_irp + R_sparse and reports tail/error.",
    }

    del a, b, ref, got_sparse, got_irp, cert
    empty_cache(device)
    return row


def parse_int_list(text: str) -> tuple[int, ...]:
    vals: list[int] = []
    for part in str(text).replace(";", ",").split(","):
        p = part.strip()
        if p:
            vals.append(int(p))
    if not vals:
        raise ValueError("empty integer list")
    return tuple(vals)


def parse_str_list(text: str) -> tuple[str, ...]:
    vals: list[str] = []
    for part in str(text).replace(";", ",").split(","):
        p = part.strip()
        if p:
            vals.append(p)
    if not vals:
        raise ValueError("empty family list")
    return tuple(vals)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="424 Pandrosion IRP converter: approximate dense A as A_irp + sparse residual, benchmark against torch.matmul."
    )
    p.add_argument("--dims", default="512,1024,2048")
    p.add_argument("--families", default="diag_lowrank,diag_lowrank_sparse,noisy_lowrank,banded,random_dense")
    p.add_argument("--columns", type=int, default=0)
    p.add_argument("--rank", type=int, default=8)
    p.add_argument("--oversample", type=int, default=6)
    p.add_argument("--alternating-steps", type=int, default=3)
    p.add_argument("--sparse-per-row", type=int, default=4)
    p.add_argument("--input-sparse-per-row", type=int, default=4)
    p.add_argument("--input-rank", type=int, default=8)
    p.add_argument("--bandwidth", type=int, default=3)
    p.add_argument("--noise", type=float, default=0.02)
    p.add_argument("--dtype", default="float32", choices=["float32", "float16", "bfloat16"])
    p.add_argument("--seed", type=int, default=424)
    p.add_argument("--repeats", type=int, default=0)
    p.add_argument("--warmup", type=int, default=3)
    p.add_argument("--out", default=str(OUT))
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cfg = Config(
        dims=parse_int_list(args.dims),
        families=parse_str_list(args.families),
        columns=int(args.columns),
        rank=int(args.rank),
        oversample=int(args.oversample),
        alternating_steps=int(args.alternating_steps),
        sparse_per_row=int(args.sparse_per_row),
        input_sparse_per_row=int(args.input_sparse_per_row),
        input_rank=int(args.input_rank),
        bandwidth=int(args.bandwidth),
        noise=float(args.noise),
        dtype=str(args.dtype),
        seed=int(args.seed),
        repeats=int(args.repeats),
        warmup=int(args.warmup),
        out=Path(str(args.out)),
    )
    device = pick_device()
    dtype = pick_dtype(cfg.dtype)
    rows: list[dict[str, Any]] = []

    for n in cfg.dims:
        for family in cfg.families:
            try:
                row = run_case(str(family), int(n), cfg, device, dtype)
            except Exception as exc:
                row = {
                    "family": str(family),
                    "n": int(n),
                    "columns": int(cfg.columns) if int(cfg.columns) > 0 else int(n),
                    "device": str(device),
                    "dtype": str(dtype).replace("torch.", ""),
                    "error": f"{type(exc).__name__}: {exc}",
                }
                empty_cache(device)
            rows.append(row)
            if "error" in row:
                print(f"{family:20s} n={int(n):5d} ERROR {row['error']}")
            else:
                print(
                    f"{family:20s} n={int(n):5d} "
                    f"424={row['pandrosion_424_irp_sparse_min_ms']:.4f} ms "
                    f"torch={row['torch_matmul_min_ms']:.4f} ms "
                    f"speed={row['speedup_irp_sparse_vs_torch_matmul']:.2f}x "
                    f"err={row['relative_l2_error_irp_sparse']:.2e} "
                    f"tail={row['tail_relative_fro']:.2e}"
                )

    ok_rows = [r for r in rows if "error" not in r]
    useful = [r for r in ok_rows if bool(r.get("useful_win_irp_sparse_tol_1e_3"))]
    result = {
        "script": "flow/424_pandrosion_standalone_irp_residual_converter_engine.py",
        "claim": "424 attempts to convert dense A into A_irp + R_sparse + R_tail and benchmarks A_irp + R_sparse against torch.matmul(A,B).",
        "non_claim": "424 does not guarantee a useful certificate for arbitrary dense random matrices; it reports the residual when conversion fails.",
        "device": str(device),
        "dtype": str(dtype).replace("torch.", ""),
        "torch_version": str(torch.__version__),
        "numpy_version": str(np.__version__),
        "config": {
            "dims": list(cfg.dims),
            "families": list(cfg.families),
            "columns": int(cfg.columns),
            "rank": int(cfg.rank),
            "oversample": int(cfg.oversample),
            "alternating_steps": int(cfg.alternating_steps),
            "sparse_per_row": int(cfg.sparse_per_row),
            "input_sparse_per_row": int(cfg.input_sparse_per_row),
            "input_rank": int(cfg.input_rank),
            "bandwidth": int(cfg.bandwidth),
            "noise": float(cfg.noise),
            "seed": int(cfg.seed),
            "repeats": int(cfg.repeats),
            "warmup": int(cfg.warmup),
        },
        "summary": {
            "rows": len(rows),
            "valid_rows": len(ok_rows),
            "useful_winning_rows_tol_1e_3": len(useful),
            "best_speedup_irp_sparse": float(max([float(r.get("speedup_irp_sparse_vs_torch_matmul", 0.0)) for r in ok_rows], default=0.0)),
            "best_error_irp_sparse": float(min([float(r.get("relative_l2_error_irp_sparse", float("inf"))) for r in ok_rows], default=float("inf"))),
            "worst_error_irp_sparse": float(max([float(r.get("relative_l2_error_irp_sparse", 0.0)) for r in ok_rows], default=0.0)),
        },
        "rows": rows,
    }

    cfg.out.parent.mkdir(parents=True, exist_ok=True)
    cfg.out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"wrote {cfg.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
