#!/usr/bin/env python3
"""
121_pandrosion_mlx_vs_118_benchmark.py

Benchmark harness for comparing:
  - 118: NumPy probe-aware PURE Thales Pandrosion engine
  - 121: the same 118 flow with MLX kernels for dense F(z) and Q(a,b)

The comparison keeps the same case, seed, root target, trial pool, epoch budget,
tolerances, and start/corrector formulas.  Only the polynomial evaluation and
telescopic finite-slope matrix kernels are replaced by MLX.

MLX currently stores complex arrays as complex64, so 121 is an accelerator
experiment rather than a bitwise-identical complex128 port of 118.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence


def _load_118() -> Any:
    path = Path(__file__).resolve().with_name("118_pandrosion_probe_aware_pure_thales_engine.py")
    spec = importlib.util.spec_from_file_location("pandrosion_118_for_121", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load 118 engine from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


ENGINE118 = _load_118()
np = ENGINE118.np
ORIGINAL_DENSE = ENGINE118.DenseKostlanSystem

try:
    import mlx.core as mx
except Exception as exc:  # pragma: no cover
    mx = None
    _MLX_IMPORT_ERROR = exc
else:
    _MLX_IMPORT_ERROR = None


def now() -> float:
    return time.perf_counter()


def parse_cases(raw: str) -> list[str]:
    return [c.strip() for c in str(raw).replace("|", ";").split(";") if c.strip()]


def mlx_device(name: str) -> Any:
    if mx is None:
        raise RuntimeError(f"MLX is not installed or cannot be imported: {_MLX_IMPORT_ERROR!r}")
    key = str(name).strip().lower()
    if key == "cpu":
        return mx.Device(mx.DeviceType.cpu, 0)
    if key == "gpu":
        return mx.Device(mx.DeviceType.gpu, 0)
    raise ValueError(f"unsupported MLX device {name!r}; use 'gpu' or 'cpu'")


def ensure_mlx(device_name: str) -> None:
    ENGINE118.ensure_numpy()
    if mx is None:
        raise RuntimeError(f"MLX is required for 121. Import error: {_MLX_IMPORT_ERROR!r}")
    mx.set_default_device(mlx_device(device_name))
    try:
        probe = mx.array(np.asarray([1.0 + 0.0j], dtype=np.complex64))
        mx.eval(probe)
    except Exception as exc:
        raise RuntimeError(
            "MLX could not create an array on the selected device. "
            "On macOS this usually means the process cannot access Metal/GPU; "
            "run from the normal venv session or allow the command outside the sandbox."
        ) from exc


def warmup_mlx(device_name: str) -> None:
    ensure_mlx(device_name)
    a = mx.ones((8, 8), dtype=mx.complex64)
    b = mx.ones((8, 8), dtype=mx.complex64)
    out = mx.matmul(a, b)
    mx.eval(out)


class MlxDenseKostlanSystem(ORIGINAL_DENSE):
    """118 DenseKostlanSystem with MLX-backed eval and slope kernels."""

    mlx_dtype = "complex64"

    @classmethod
    def make(
        cls,
        n: int,
        d: int,
        seed_index: int = 0,
        equation_normalize: bool = True,
    ) -> "MlxDenseKostlanSystem":
        base = ORIGINAL_DENSE.make(
            n,
            d,
            seed_index=seed_index,
            equation_normalize=equation_normalize,
        )
        obj = cls(
            n=base.n,
            d=base.d,
            seed=base.seed,
            exps=base.exps,
            coeff=base.coeff,
            weights=base.weights,
            equation_normalize=base.equation_normalize,
        )
        obj._generation_seconds = float(getattr(base, "_generation_seconds", 0.0))
        obj._prepare_mlx()
        return obj

    def _prepare_mlx(self) -> None:
        if mx is None:
            raise RuntimeError(f"MLX is required for 121. Import error: {_MLX_IMPORT_ERROR!r}")
        self._mx_coeff = mx.array(np.asarray(self.coeff, dtype=np.complex64))
        self._mx_exp_cols = [
            mx.array(np.asarray(self.exps[:, j], dtype=np.int32))
            for j in range(int(self.n))
        ]
        self._mx_terms = int(self.exps.shape[0])
        mx.eval(self._mx_coeff)

    def _mx_powers(self, z: Sequence[complex]) -> list[Any]:
        zz = mx.array(np.asarray(z, dtype=np.complex64))
        tables: list[Any] = []
        for j in range(int(self.n)):
            zj = zz[j]
            vals = [mx.array(np.complex64(1.0 + 0.0j))]
            cur = vals[0]
            for _ in range(1, int(self.d) + 1):
                cur = cur * zj
                vals.append(cur)
            tables.append(mx.stack(vals))
        return tables

    def monomials(self, z: Sequence[complex]) -> Any:
        pows = self._mx_powers(z)
        mon = mx.ones((int(self.terms_per_poly),), dtype=mx.complex64)
        for j in range(int(self.n)):
            mon = mon * pows[j][self._mx_exp_cols[j]]
        return mon

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = ENGINE118.now()
        mon = self.monomials(z)
        f = mx.matmul(self._mx_coeff, mon)
        mx.eval(f)
        out = np.asarray(f, dtype=np.complex64).astype(np.complex128)
        self.eval_count += 1
        self.seconds_eval += ENGINE118.now() - t0
        return out

    def _slope_power_table(self, aj: complex, bj: complex, pows_b_j: Any) -> Any:
        vals = [mx.array(np.complex64(0.0 + 0.0j))]
        acc = vals[0]
        ajx = mx.array(np.complex64(aj))
        for m in range(1, int(self.d) + 1):
            acc = pows_b_j[m - 1] + ajx * acc
            vals.append(acc)
        return mx.stack(vals)

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        t0 = ENGINE118.now()
        pows_a = self._mx_powers(a)
        pows_b = self._mx_powers(b)
        n = int(self.n)
        m_terms = int(self.terms_per_poly)

        pa_cols = [pows_a[j][self._mx_exp_cols[j]] for j in range(n)]
        pb_cols = [pows_b[j][self._mx_exp_cols[j]] for j in range(n)]

        prefix_b: list[Any] = [None] * (n + 1)
        suffix_a: list[Any] = [None] * (n + 1)
        prefix_b[0] = mx.ones((m_terms,), dtype=mx.complex64)
        for j in range(n):
            prefix_b[j + 1] = prefix_b[j] * pb_cols[j]
        suffix_a[n] = mx.ones((m_terms,), dtype=mx.complex64)
        for j in range(n - 1, -1, -1):
            suffix_a[j] = suffix_a[j + 1] * pa_cols[j]

        aa = np.asarray(a, dtype=np.complex64)
        bb = np.asarray(b, dtype=np.complex64)
        cols: list[Any] = []
        for j in range(n):
            slope_table = self._slope_power_table(aa[j], bb[j], pows_b[j])
            cols.append(prefix_b[j] * suffix_a[j + 1] * slope_table[self._mx_exp_cols[j]])
        terms = mx.stack(cols, axis=1)
        q = mx.matmul(self._mx_coeff, terms)
        mx.eval(q)

        out = np.asarray(q, dtype=np.complex64).astype(np.complex128)
        self.slope_count = int(getattr(self, "slope_count", 0)) + 1
        self.seconds_slope = float(getattr(self, "seconds_slope", 0.0)) + (ENGINE118.now() - t0)
        return out

    def stats(self) -> dict[str, Any]:
        stats = super().stats()
        stats.update(
            {
                "accelerator": "mlx",
                "mlx_device": str(mx.default_device()) if mx is not None else None,
                "mlx_dtype": self.mlx_dtype,
                "mlx_note": "complex128 coefficients are evaluated through MLX complex64 kernels",
            }
        )
        return stats


def run_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    ensure_mlx(str(args.mlx_device))
    old_dense = ENGINE118.DenseKostlanSystem
    ENGINE118.DenseKostlanSystem = MlxDenseKostlanSystem
    try:
        result = ENGINE118.run_case(args, case_raw)
    finally:
        ENGINE118.DenseKostlanSystem = old_dense
    result["script"] = "121_pandrosion_mlx_vs_118_benchmark.py"
    result["source_engine"] = "flow/118_pandrosion_probe_aware_pure_thales_engine.py"
    result["mode"] = str(result.get("mode", "")) + "/mlx-eval-slope-kernels"
    result["dependencies"] = {
        "python_scripts": ["flow/118_pandrosion_probe_aware_pure_thales_engine.py"],
        "numpy": bool(np is not None),
        "mlx": bool(mx is not None),
    }
    result["mlx"] = {
        "device": str(mx.default_device()) if mx is not None else None,
        "dtype": MlxDenseKostlanSystem.mlx_dtype,
        "complex_precision": "complex64",
    }
    return result


def residual_stats(values: Sequence[float]) -> dict[str, Optional[float]]:
    finite = [float(v) for v in values if math.isfinite(float(v))]
    if not finite:
        return {"residual_min": None, "residual_mean": None, "residual_max": None}
    return {
        "residual_min": float(min(finite)),
        "residual_mean": float(sum(finite) / len(finite)),
        "residual_max": float(max(finite)),
    }


def probe_eval_sum(result: dict[str, Any]) -> Optional[int]:
    total = 0
    seen = False
    for trial in result.get("trials", []):
        value = trial.get("probe_total_evals")
        if value is None:
            continue
        try:
            total += int(value)
            seen = True
        except Exception:
            pass
    return total if seen else None


def summarize_result(engine: str, case_raw: str, result: dict[str, Any], seconds: float) -> dict[str, Any]:
    summary = result.get("summary", {})
    eval_stats = summary.get("eval_stats", {})
    residuals = [float(r.get("residual", float("nan"))) for r in result.get("roots", [])]
    out = {
        "case": case_raw,
        "engine": engine,
        "status": "ok",
        "success": bool(summary.get("success")),
        "requested_roots": int(summary.get("requested_roots", 0)),
        "unique_roots": int(summary.get("unique_roots", 0)),
        "trials_used": int(summary.get("trials_used", 0)),
        "duplicates": int(summary.get("duplicates", 0)),
        "failures": int(summary.get("failures", 0)),
        "seconds": float(seconds),
        "reported_seconds": float(summary.get("total_seconds", seconds)),
        "extract_seconds": float(summary.get("extract_seconds", 0.0)),
        "generation_seconds": float(summary.get("generation_seconds", 0.0)),
        "eval_used": int(eval_stats.get("eval_count", 0)) if eval_stats.get("eval_count") is not None else None,
        "slope_calls": int(eval_stats.get("slope_count", 0)) if eval_stats.get("slope_count") is not None else None,
        "seconds_eval": float(eval_stats.get("seconds_eval", 0.0)),
        "seconds_slope": float(eval_stats.get("seconds_slope", 0.0)),
        "probe_eval_sum": probe_eval_sum(result),
        "seed": int(result.get("seed", 0)),
        "terms": int(result.get("terms", 0)),
        "bezout": int(result.get("bezout", 0)),
        "error": None,
    }
    out.update(residual_stats(residuals))
    return out


def error_row(engine: str, case_raw: str, exc: BaseException, seconds: float) -> dict[str, Any]:
    return {
        "case": case_raw,
        "engine": engine,
        "status": "error",
        "success": False,
        "requested_roots": None,
        "unique_roots": 0,
        "trials_used": None,
        "duplicates": None,
        "failures": None,
        "seconds": float(seconds),
        "reported_seconds": None,
        "extract_seconds": None,
        "generation_seconds": None,
        "eval_used": None,
        "slope_calls": None,
        "seconds_eval": None,
        "seconds_slope": None,
        "probe_eval_sum": None,
        "seed": None,
        "terms": None,
        "bezout": None,
        "residual_min": None,
        "residual_mean": None,
        "residual_max": None,
        "error": f"{type(exc).__name__}: {exc}",
    }


def run_engine(engine: str, fn: Any, args: argparse.Namespace, case_raw: str) -> tuple[dict[str, Any], dict[str, Any]]:
    t0 = now()
    try:
        result = fn(args, case_raw)
    except Exception as exc:
        return error_row(engine, case_raw, exc, now() - t0), {}
    return summarize_result(engine, case_raw, result, now() - t0), result


def compare_pair(row118: dict[str, Any], row121: dict[str, Any]) -> dict[str, Any]:
    roots118 = int(row118.get("unique_roots") or 0)
    roots121 = int(row121.get("unique_roots") or 0)
    seconds118 = float(row118.get("seconds") or 0.0)
    seconds121 = float(row121.get("seconds") or 0.0)
    if roots121 > roots118:
        winner = "121(118+MLX)"
    elif roots118 > roots121:
        winner = "118"
    elif bool(row121.get("success")) and not bool(row118.get("success")):
        winner = "121(118+MLX)"
    elif bool(row118.get("success")) and not bool(row121.get("success")):
        winner = "118"
    elif seconds121 > 0 and seconds118 > 0:
        winner = "121(118+MLX)" if seconds121 < seconds118 else "118"
    else:
        winner = "tie"
    return {
        "case": row118.get("case"),
        "roots_118": roots118,
        "roots_121": roots121,
        "root_delta_121_minus_118": roots121 - roots118,
        "success_118": bool(row118.get("success")),
        "success_121": bool(row121.get("success")),
        "seconds_118": seconds118,
        "seconds_121": seconds121,
        "speedup_121_vs_118": (seconds118 / seconds121 if seconds121 > 0 else None),
        "evals_118": row118.get("eval_used"),
        "evals_121": row121.get("eval_used"),
        "slopes_118": row118.get("slope_calls"),
        "slopes_121": row121.get("slope_calls"),
        "residual_min_118": row118.get("residual_min"),
        "residual_min_121": row121.get("residual_min"),
        "winner": winner,
    }


def fmt_float(value: Any, digits: int = 4) -> str:
    if value is None:
        return ""
    try:
        x = float(value)
    except Exception:
        return ""
    if not math.isfinite(x):
        return ""
    return f"{x:.{digits}f}"


def fmt_sci(value: Any) -> str:
    if value is None:
        return ""
    try:
        x = float(value)
    except Exception:
        return ""
    if not math.isfinite(x):
        return ""
    return f"{x:.3e}"


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    fields = [
        "case",
        "engine",
        "status",
        "success",
        "requested_roots",
        "unique_roots",
        "trials_used",
        "duplicates",
        "failures",
        "seconds",
        "reported_seconds",
        "extract_seconds",
        "generation_seconds",
        "eval_used",
        "slope_calls",
        "seconds_eval",
        "seconds_slope",
        "probe_eval_sum",
        "seed",
        "terms",
        "bezout",
        "residual_min",
        "residual_mean",
        "residual_max",
        "error",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_pair_csv(path: Path, pairs: Sequence[dict[str, Any]]) -> None:
    fields = [
        "case",
        "roots_118",
        "roots_121",
        "root_delta_121_minus_118",
        "success_118",
        "success_121",
        "seconds_118",
        "seconds_121",
        "speedup_121_vs_118",
        "evals_118",
        "evals_121",
        "slopes_118",
        "slopes_121",
        "residual_min_118",
        "residual_min_121",
        "winner",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in pairs:
            writer.writerow(row)


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    p = payload["parameters"]
    s = payload["summary"]
    lines = [
        "# 121 MLX vs 118 Benchmark",
        "",
        f"- cases: `{', '.join(payload['cases'])}`",
        f"- count/pool/epochs: `{p['count']}/{p['pool']}/{p['epochs']}`",
        f"- MLX device: `{p['mlx_device']}`",
        f"- MLX complex dtype: `complex64`",
        "",
        "## Summary",
        "",
        f"- 118 roots: `{s['roots_118']}/{s['requested_118']}`; complete cases: `{s['successes_118']}/{s['pairs']}`",
        f"- 121 roots: `{s['roots_121']}/{s['requested_121']}`; complete cases: `{s['successes_121']}/{s['pairs']}`",
        f"- wall seconds: 118=`{s['seconds_118']:.4f}`, 121=`{s['seconds_121']:.4f}`, total=`{s['seconds']:.4f}`",
        f"- root delta 121-118: `{s['root_delta_121_minus_118']}`",
        "",
        "## Pair Results",
        "",
        "| case | 118 roots | 121 roots | delta | 118 sec | 121 sec | 121 speedup | 118 evals | 121 evals | winner |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in payload["pairs"]:
        lines.append(
            f"| `{row.get('case')}` | {row.get('roots_118')} | {row.get('roots_121')} | "
            f"{row.get('root_delta_121_minus_118')} | {fmt_float(row.get('seconds_118'), 4)} | "
            f"{fmt_float(row.get('seconds_121'), 4)} | {fmt_float(row.get('speedup_121_vs_118'), 3)} | "
            f"{row.get('evals_118') or ''} | {row.get('evals_121') or ''} | `{row.get('winner')}` |"
        )
    lines.extend(
        [
            "",
            "## Residual Minima",
            "",
            "| case | 118 min residual | 121 min residual |",
            "|---|---:|---:|",
        ]
    )
    for row in payload["pairs"]:
        lines.append(
            f"| `{row.get('case')}` | {fmt_sci(row.get('residual_min_118'))} | "
            f"{fmt_sci(row.get('residual_min_121'))} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    p = ENGINE118.build_parser()
    p.description = "121 benchmark: compare 118 against the same 118 flow accelerated with MLX kernels."
    p.set_defaults(
        cases="2,34",
        count=8,
        pool=512,
        epochs=24,
        outdir="verification/121_mlx_vs_118_benchmark",
    )
    p.add_argument("--mlx-device", choices=["gpu", "cpu"], default="gpu")
    p.add_argument("--warmup-mlx", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--run-order", choices=["118-first", "121-first"], default="118-first")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ENGINE118.ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)
    cases = parse_cases(args.cases)
    if bool(args.warmup_mlx):
        warmup_mlx(str(args.mlx_device))

    rows: list[dict[str, Any]] = []
    pairs: list[dict[str, Any]] = []
    full_results: dict[str, Any] = {}
    t_all = now()

    print("=" * 120, flush=True)
    print("121 benchmark: 121(118 + MLX eval/slope kernels) vs 118 NumPy", flush=True)
    print(
        f"cases={cases} count={args.count} pool={args.pool} epochs={args.epochs} "
        f"mlx_device={args.mlx_device} run_order={args.run_order}",
        flush=True,
    )
    print("=" * 120, flush=True)

    for case_raw in cases:
        print(f"case={case_raw}", flush=True)
        if args.run_order == "121-first":
            row121, result121 = run_engine("121(118+MLX)", run_case, args, case_raw)
            row118, result118 = run_engine("118", ENGINE118.run_case, args, case_raw)
        else:
            row118, result118 = run_engine("118", ENGINE118.run_case, args, case_raw)
            row121, result121 = run_engine("121(118+MLX)", run_case, args, case_raw)

        rows.extend([row118, row121])
        print(
            f"  118: status={row118.get('status')} roots={row118.get('unique_roots')}/{row118.get('requested_roots')} "
            f"evals={row118.get('eval_used')} slopes={row118.get('slope_calls')} seconds={float(row118.get('seconds') or 0.0):.4f}",
            flush=True,
        )
        print(
            f"  121(118+MLX): status={row121.get('status')} roots={row121.get('unique_roots')}/{row121.get('requested_roots')} "
            f"evals={row121.get('eval_used')} slopes={row121.get('slope_calls')} seconds={float(row121.get('seconds') or 0.0):.4f}",
            flush=True,
        )
        pair = compare_pair(row118, row121)
        pairs.append(pair)
        print(
            f"  winner={pair['winner']} delta121-118={pair['root_delta_121_minus_118']} "
            f"speedup121={fmt_float(pair.get('speedup_121_vs_118'), 3)}",
            flush=True,
        )
        full_results[case_raw.replace(",", "x")] = {"118": result118, "121(118+MLX)": result121}

    seconds_all = now() - t_all
    summary = {
        "pairs": len(pairs),
        "runs": len(rows),
        "successes_118": int(sum(1 for row in rows if row.get("engine") == "118" and row.get("success") is True)),
        "successes_121": int(sum(1 for row in rows if row.get("engine") == "121(118+MLX)" and row.get("success") is True)),
        "roots_118": int(sum(int(row.get("unique_roots") or 0) for row in rows if row.get("engine") == "118")),
        "roots_121": int(sum(int(row.get("unique_roots") or 0) for row in rows if row.get("engine") == "121(118+MLX)")),
        "requested_118": int(sum(int(row.get("requested_roots") or 0) for row in rows if row.get("engine") == "118")),
        "requested_121": int(sum(int(row.get("requested_roots") or 0) for row in rows if row.get("engine") == "121(118+MLX)")),
        "seconds_118": float(sum(float(row.get("seconds") or 0.0) for row in rows if row.get("engine") == "118")),
        "seconds_121": float(sum(float(row.get("seconds") or 0.0) for row in rows if row.get("engine") == "121(118+MLX)")),
        "root_delta_121_minus_118": int(sum(int(pair.get("root_delta_121_minus_118") or 0) for pair in pairs)),
        "wins_118": int(sum(1 for pair in pairs if pair.get("winner") == "118")),
        "wins_121": int(sum(1 for pair in pairs if pair.get("winner") == "121(118+MLX)")),
        "ties": int(sum(1 for pair in pairs if pair.get("winner") == "tie")),
        "seconds": float(seconds_all),
    }

    payload = {
        "script": "121_pandrosion_mlx_vs_118_benchmark.py",
        "engines": {
            "118": "flow/118_pandrosion_probe_aware_pure_thales_engine.py",
            "121(118+MLX)": "flow/118_pandrosion_probe_aware_pure_thales_engine.py with MLX eval/slope kernels",
        },
        "cases": cases,
        "parameters": {
            "count": int(args.count),
            "pool": int(args.pool),
            "epochs": int(args.epochs),
            "accept": float(args.accept),
            "tol": float(args.tol),
            "cluster_sep": float(args.cluster_sep),
            "line_search": int(args.line_search),
            "probe_scale": float(args.probe_scale),
            "probe_candidates": int(args.probe_candidates),
            "probe_radii": str(args.probe_radii),
            "probe_self": bool(args.probe_self),
            "seed_index": int(args.seed_index),
            "mlx_device": str(args.mlx_device),
            "run_order": str(args.run_order),
        },
        "summary": summary,
        "pairs": pairs,
        "runs": rows,
        "full_results": full_results if bool(args.verbose_trials) else {},
    }

    out_json = Path(args.out) if args.out else Path(args.outdir) / "121_mlx_vs_118_benchmark.json"
    out_csv = out_json.with_suffix(".csv")
    out_pairs_csv = out_json.with_name(out_json.stem + "_pairs.csv")
    out_md = out_json.with_suffix(".md")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_csv(out_csv, rows)
    write_pair_csv(out_pairs_csv, pairs)
    write_markdown(out_md, payload)

    print("=" * 120, flush=True)
    print(
        f"summary: 118 roots={summary['roots_118']}/{summary['requested_118']} "
        f"successes={summary['successes_118']}/{summary['pairs']} seconds={summary['seconds_118']:.4f}",
        flush=True,
    )
    print(
        f"summary: 121 roots={summary['roots_121']}/{summary['requested_121']} "
        f"successes={summary['successes_121']}/{summary['pairs']} seconds={summary['seconds_121']:.4f}",
        flush=True,
    )
    print(f"root_delta_121_minus_118={summary['root_delta_121_minus_118']}", flush=True)
    print(f"json={out_json}", flush=True)
    print(f"csv={out_csv}", flush=True)
    print(f"pairs_csv={out_pairs_csv}", flush=True)
    print(f"report={out_md}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    finally:
        try:
            sys.stdout.flush()
            sys.stderr.flush()
        except Exception:
            pass
