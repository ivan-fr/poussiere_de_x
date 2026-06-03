#!/usr/bin/env python3
"""423 Pandrosion full-IRP structured matmul engine.

This script is the honest continuation of 422 for "beat torch.matmul":

    compute the same dense result Y = A @ B,
    but represent A by a Pandrosion IRP certificate and apply that certificate
    without calling torch.matmul inside the 423 engine.

The 423 solve is a full-n identity IRP step on the matrix target

    F(Y) = Y - Apply_A(B).

Starting from Y0=0, the full local-jet correction is exact:

    delta = -F(0) = Apply_A(B).

This does not read a hidden root.  The right-hand side is produced by the
declared action of A.  When A has exploitable structure, this can beat dense
torch.matmul(A, B) while returning the same Y.

Supported A certificates:
  - diagonal
  - scaled permutation
  - diagonal plus low rank
  - circulant
  - banded

Non-claim:
  This does not beat torch.matmul for arbitrary dense unstructured A.  See 440
  for that negative control.

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


OUT = Path("/private/tmp/423_full_irp_structured_matmul_vs_torch_benchmark.json")


@dataclasses.dataclass(frozen=True)
class Config:
    dims: tuple[int, ...] = (512, 1024, 2048, 4096)
    families: tuple[str, ...] = ("diagonal", "scaled_permutation", "diag_lowrank", "circulant")
    columns: int = 0
    rank: int = 8
    bandwidth: int = 3
    reduced_contractions: bool = True
    dtype: str = "float32"
    seed: int = 423
    repeats: int = 0
    warmup: int = 4
    out: Path = OUT


@dataclasses.dataclass
class IRPCertificate:
    family: str
    n: int
    data: dict[str, Any]
    certificate_size: int
    dense_entries: int
    full_irp_dimension: int


@dataclasses.dataclass
class MatrixIdentityTarget:
    cert: IRPCertificate
    b: Any
    evals: int = 0
    applications: int = 0

    def apply(self) -> Any:
        self.applications += 1
        return apply_irp_certificate(self.cert, self.b)

    def eval(self, y: Any) -> Any:
        self.evals += 1
        return y - self.apply()

    def jvp_identity(self, d: Any) -> Any:
        return d


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


def torch_randperm(n: int, *, seed: int, device: Any) -> Any:
    gen = torch.Generator(device="cpu")
    gen.manual_seed(int(seed))
    return torch.randperm(int(n), generator=gen).to(device=device)


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


def default_repeats(n: int, family: str) -> int:
    if family == "circulant":
        return 12 if n <= 1024 else (7 if n <= 2048 else 4)
    if n <= 1024:
        return 18
    if n <= 2048:
        return 10
    return 5


def family_id(family: str) -> int:
    return {
        "diagonal": 1,
        "scaled_permutation": 2,
        "diag_lowrank": 3,
        "circulant": 4,
        "banded": 5,
    }.get(str(family), 99)


def make_certificate(family: str, n: int, cfg: Config, device: Any, dtype: Any, seed: int) -> IRPCertificate:
    n = int(n)
    fam = str(family)
    if fam == "diagonal":
        diag = torch_randn((n,), seed=seed + 11, device=device, dtype=dtype, scale=0.5)
        return IRPCertificate(fam, n, {"diag": diag}, n, n * n, n)

    if fam == "scaled_permutation":
        perm = torch_randperm(n, seed=seed + 17, device=device)
        scale = torch_randn((n,), seed=seed + 19, device=device, dtype=dtype, scale=0.5)
        return IRPCertificate(fam, n, {"perm": perm, "scale": scale}, 2 * n, n * n, n)

    if fam == "diag_lowrank":
        r = max(1, int(cfg.rank))
        diag = torch_randn((n,), seed=seed + 23, device=device, dtype=dtype, scale=0.25)
        u = torch_randn((n, r), seed=seed + 29, device=device, dtype=dtype, scale=1.0 / math.sqrt(max(1, r)))
        v = torch_randn((n, r), seed=seed + 31, device=device, dtype=dtype, scale=1.0 / math.sqrt(max(1, n)))
        return IRPCertificate(
            fam,
            n,
            {
                "diag": diag,
                "u": u,
                "v": v,
                "rank": r,
                "reduced_contractions": bool(cfg.reduced_contractions),
            },
            n + 2 * n * r,
            n * n,
            n,
        )

    if fam == "circulant":
        first_col = torch_randn((n,), seed=seed + 37, device=device, dtype=dtype, scale=1.0 / math.sqrt(max(1, n)))
        return IRPCertificate(fam, n, {"first_col": first_col}, n, n * n, n)

    if fam == "banded":
        w = max(0, int(cfg.bandwidth))
        diags: list[tuple[int, Any]] = []
        cert_size = 0
        for offset in range(-w, w + 1):
            length = n - abs(offset)
            vals = torch_randn((length,), seed=seed + 101 + offset + 3 * w, device=device, dtype=dtype, scale=0.2)
            diags.append((offset, vals))
            cert_size += length
        return IRPCertificate(fam, n, {"bandwidth": w, "diags": diags}, cert_size, n * n, n)

    raise ValueError(f"unknown family {family!r}")


def apply_diag_lowrank_irp(diag: Any, u: Any, v: Any, b: Any, *, reduced_contractions: bool) -> Any:
    out = diag[:, None] * b
    if bool(reduced_contractions):
        # IRP reduced contractions: rank r is small, so this avoids dense n-by-n
        # GEMM while using optimized kernels for the two skinny reductions.
        return out + torch.mm(u, torch.mm(v.transpose(0, 1), b))
    r = int(u.shape[1])
    # Pure elementwise IRP fallback: sum_i u_i outer (v_i dot B).
    for i in range(r):
        coeff = torch.sum(v[:, i].unsqueeze(1) * b, dim=0)
        out = out + u[:, i].unsqueeze(1) * coeff.unsqueeze(0)
    return out


def apply_circulant_fft(first_col: Any, b: Any) -> Any:
    # Circulant with first column c: y_i = sum_j c[(i-j) mod n] b_j.
    c64 = first_col.to(torch.float32)
    b64 = b.to(torch.float32)
    fc = torch.fft.fft(c64.to(torch.complex64), dim=0).unsqueeze(1)
    fb = torch.fft.fft(b64.to(torch.complex64), dim=0)
    y = torch.fft.ifft(fc * fb, dim=0).real
    return y.to(dtype=b.dtype)


def apply_irp_certificate(cert: IRPCertificate, b: Any) -> Any:
    fam = cert.family
    if fam == "diagonal":
        return cert.data["diag"][:, None] * b

    if fam == "scaled_permutation":
        return cert.data["scale"][:, None] * b.index_select(0, cert.data["perm"])

    if fam == "diag_lowrank":
        return apply_diag_lowrank_irp(
            cert.data["diag"],
            cert.data["u"],
            cert.data["v"],
            b,
            reduced_contractions=bool(cert.data.get("reduced_contractions", True)),
        )

    if fam == "circulant":
        return apply_circulant_fft(cert.data["first_col"], b)

    if fam == "banded":
        n = int(cert.n)
        out = torch.zeros_like(b)
        for offset, vals in cert.data["diags"]:
            k = int(offset)
            if k >= 0:
                out[: n - k, :] += vals[:, None] * b[k:n, :]
            else:
                out[-k:n, :] += vals[:, None] * b[: n + k, :]
        return out

    raise ValueError(f"unknown family {fam!r}")


def materialize_dense_a(cert: IRPCertificate, *, device: Any, dtype: Any) -> Any:
    n = int(cert.n)
    fam = cert.family
    a = torch.zeros((n, n), device=device, dtype=dtype)

    if fam == "diagonal":
        idx = torch.arange(n, device=device)
        a[idx, idx] = cert.data["diag"]
        return a

    if fam == "scaled_permutation":
        rows = torch.arange(n, device=device)
        a[rows, cert.data["perm"]] = cert.data["scale"]
        return a

    if fam == "diag_lowrank":
        idx = torch.arange(n, device=device)
        a[idx, idx] = cert.data["diag"]
        for i in range(int(cert.data["rank"])):
            a += cert.data["u"][:, i].unsqueeze(1) * cert.data["v"][:, i].unsqueeze(0)
        return a

    if fam == "circulant":
        rows = torch.arange(n, device=device).unsqueeze(1)
        cols = torch.arange(n, device=device).unsqueeze(0)
        return cert.data["first_col"][(rows - cols) % n].to(dtype=dtype)

    if fam == "banded":
        for offset, vals in cert.data["diags"]:
            k = int(offset)
            if k >= 0:
                rows = torch.arange(0, n - k, device=device)
                cols = rows + k
            else:
                cols = torch.arange(0, n + k, device=device)
                rows = cols - k
            a[rows, cols] = vals
        return a

    raise ValueError(f"unknown family {fam!r}")


def pandrosion_423_full_irp(target: MatrixIdentityTarget, *, verify_residual: bool = False) -> tuple[Any, dict[str, Any]]:
    # For F(Y)=Y-Apply_A(B), the full-n IRP identity correction from Y0=0 is
    # exactly delta=-F(0)=Apply_A(B).  We do not need to allocate Y0 or recompute
    # the residual on the timed path.
    y = target.apply()
    rel = 0.0
    if bool(verify_residual):
        residual = target.eval(y)
        rel_num = torch.linalg.vector_norm(residual.reshape(-1).to(torch.float32)).item()
        denom = torch.linalg.vector_norm(y.reshape(-1).to(torch.float32)).item()
        rel = float(rel_num / max(1e-30, denom))
    return y, {
        "solver": "423-full-n-identity-irp-structured-operator",
        "full_irp": True,
        "no_external_fallback": True,
        "engine_torch_matmul_calls": 0,
        "engine_reduced_contractions": int(2 if target.cert.family == "diag_lowrank" and bool(target.cert.data.get("reduced_contractions", True)) else 0),
        "target_evals": int(target.evals),
        "target_applications": int(target.applications),
        "relative_identity_residual": float(rel),
    }


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
    cert = make_certificate(family, n, cfg, device, dtype, seed)
    b = torch_randn((n, m), seed=seed + 313, device=device, dtype=dtype, scale=0.5)
    dense_a = materialize_dense_a(cert, device=device, dtype=dtype)
    target = MatrixIdentityTarget(cert, b)

    repeats = int(cfg.repeats) if int(cfg.repeats) > 0 else default_repeats(n, family)
    with torch.no_grad():
        ref = torch.matmul(dense_a, b)
        got, meta = pandrosion_423_full_irp(target, verify_residual=True)
        err = error_metrics(got, ref)

    irp_stats, _ = bench(lambda: pandrosion_423_full_irp(MatrixIdentityTarget(cert, b), verify_residual=False)[0], device=device, repeats=repeats, warmup=int(cfg.warmup))
    torch_stats, _ = bench(lambda: torch.matmul(dense_a, b), device=device, repeats=repeats, warmup=int(cfg.warmup))

    irp_ms = float(irp_stats["min_ms"])
    torch_ms = float(torch_stats["min_ms"])
    speed = torch_ms / max(1e-30, irp_ms)
    same_tol = 5e-3 if dtype in {torch.float16, torch.bfloat16} else 5e-5
    dense_work = 2.0 * n * n * m
    if family == "diagonal":
        irp_work = float(n * m)
    elif family == "scaled_permutation":
        irp_work = float(n * m)
    elif family == "diag_lowrank":
        r = int(cert.data["rank"])
        irp_work = float(n * m + 2 * n * r * m)
    elif family == "circulant":
        irp_work = float(8.0 * n * m * math.log2(max(2, n)))
    elif family == "banded":
        irp_work = float(cert.certificate_size * m)
    else:
        irp_work = dense_work

    row = {
        "family": str(family),
        "n": n,
        "columns": m,
        "device": str(device),
        "dtype": str(dtype).replace("torch.", ""),
        "same_result_as_torch_matmul": bool(float(err["relative_l2"]) < same_tol),
        "relative_l2_error": float(err["relative_l2"]),
        "max_abs_error": float(err["max_abs"]),
        "pandrosion_423_min_ms": irp_ms,
        "pandrosion_423_median_ms": float(irp_stats["median_ms"]),
        "torch_matmul_min_ms": torch_ms,
        "torch_matmul_median_ms": float(torch_stats["median_ms"]),
        "speedup_vs_torch_matmul": float(speed),
        "beats_torch_matmul": bool(speed > 1.0),
        "certificate_size": int(cert.certificate_size),
        "dense_entries": int(cert.dense_entries),
        "certificate_compression": float(cert.dense_entries / max(1, cert.certificate_size)),
        "full_irp_dimension": int(cert.full_irp_dimension),
        "engine_torch_matmul_calls": int(meta["engine_torch_matmul_calls"]),
        "engine_reduced_contractions": int(meta["engine_reduced_contractions"]),
        "no_external_fallback": bool(meta["no_external_fallback"]),
        "target_evals": int(meta["target_evals"]),
        "target_applications": int(meta["target_applications"]),
        "identity_residual": float(meta["relative_identity_residual"]),
        "nominal_dense_work": float(dense_work),
        "nominal_irp_work": float(irp_work),
        "nominal_work_reduction": float(dense_work / max(1.0, irp_work)),
        "interpretation": "Same Y=A@B as dense torch.matmul for a certified structured A; not arbitrary dense.",
    }

    del dense_a, b, ref, got, cert, target
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
        description="423 full Pandrosion IRP structured matmul: same A@B as torch.matmul, faster when A has an IRP certificate."
    )
    p.add_argument("--dims", default="512,1024,2048,4096")
    p.add_argument("--families", default="diagonal,scaled_permutation,diag_lowrank,circulant")
    p.add_argument("--columns", type=int, default=0, help="B columns; default 0 means square B.")
    p.add_argument("--rank", type=int, default=8)
    p.add_argument("--bandwidth", type=int, default=3)
    p.add_argument("--reduced-contractions", action="store_true", default=True)
    p.add_argument("--no-reduced-contractions", dest="reduced_contractions", action="store_false")
    p.add_argument("--dtype", default="float32", choices=["float32", "float16", "bfloat16"])
    p.add_argument("--seed", type=int, default=423)
    p.add_argument("--repeats", type=int, default=0)
    p.add_argument("--warmup", type=int, default=4)
    p.add_argument("--out", default=str(OUT))
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cfg = Config(
        dims=parse_int_list(args.dims),
        families=parse_str_list(args.families),
        columns=int(args.columns),
        rank=int(args.rank),
        bandwidth=int(args.bandwidth),
        reduced_contractions=bool(args.reduced_contractions),
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
                print(f"{family:18s} n={int(n):5d} ERROR {row['error']}")
            else:
                print(
                    f"{family:18s} n={int(n):5d} "
                    f"423={row['pandrosion_423_min_ms']:.4f} ms "
                    f"torch={row['torch_matmul_min_ms']:.4f} ms "
                    f"speed={row['speedup_vs_torch_matmul']:.2f}x "
                    f"err={row['relative_l2_error']:.2e}"
                )

    ok_rows = [r for r in rows if "error" not in r]
    wins = [r for r in ok_rows if bool(r.get("beats_torch_matmul")) and bool(r.get("same_result_as_torch_matmul"))]
    result = {
        "script": "flow/423_pandrosion_standalone_full_irp_structured_matmul_engine.py",
        "claim": "423 computes the same Y=A@B as torch.matmul through a full Pandrosion IRP identity correction when A has a certified exploitable structure.",
        "non_claim": "423 does not beat torch.matmul on arbitrary dense unstructured matrices.",
        "device": str(device),
        "dtype": str(dtype).replace("torch.", ""),
        "torch_version": str(torch.__version__),
        "numpy_version": str(np.__version__),
        "config": {
            "dims": list(cfg.dims),
            "families": list(cfg.families),
            "columns": int(cfg.columns),
            "rank": int(cfg.rank),
            "bandwidth": int(cfg.bandwidth),
            "reduced_contractions": bool(cfg.reduced_contractions),
            "seed": int(cfg.seed),
            "repeats": int(cfg.repeats),
            "warmup": int(cfg.warmup),
        },
        "summary": {
            "rows": len(rows),
            "valid_rows": len(ok_rows),
            "winning_same_result_rows": len(wins),
            "best_speedup": float(max([float(r.get("speedup_vs_torch_matmul", 0.0)) for r in ok_rows], default=0.0)),
            "worst_relative_l2_error": float(max([float(r.get("relative_l2_error", 0.0)) for r in ok_rows], default=0.0)),
        },
        "rows": rows,
    }

    cfg.out.parent.mkdir(parents=True, exist_ok=True)
    cfg.out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"wrote {cfg.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
