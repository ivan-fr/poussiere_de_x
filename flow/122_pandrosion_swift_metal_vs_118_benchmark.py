#!/usr/bin/env python3
"""
122_pandrosion_swift_metal_vs_118_benchmark.py

Swift + Metal feasibility benchmark for the 118 dense Kostlan evaluator.

This does not claim to be a complete Swift port of the 118 corrector.  It tests
the GPU-friendly part that can actually benefit from Metal: evaluating the same
118 polynomial system on a large batch of complex points without per-point GPU
round trips.  The report includes a normal 118 solver run for context, plus a
direct evaluator comparison against 118's current NumPy eval loop.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
FLOW = Path(__file__).resolve().parent


def _load_118() -> Any:
    path = FLOW / "118_pandrosion_probe_aware_pure_thales_engine.py"
    spec = importlib.util.spec_from_file_location("pandrosion_118_for_122", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load 118 engine from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


ENGINE118 = _load_118()
np = ENGINE118.np


def now() -> float:
    return time.perf_counter()


def parse_cases(raw: str) -> list[str]:
    return [c.strip() for c in str(raw).replace("|", ";").split(";") if c.strip()]


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


def build_swift_helper(source: Path, binary: Path, force: bool = False) -> dict[str, Any]:
    binary.parent.mkdir(parents=True, exist_ok=True)
    needs_build = force or not binary.exists() or source.stat().st_mtime > binary.stat().st_mtime
    out: dict[str, Any] = {"source": str(source), "binary": str(binary), "rebuilt": bool(needs_build)}
    if not needs_build:
        return out
    cmd = [
        "swiftc",
        "-O",
        str(source),
        "-o",
        str(binary),
        "-framework",
        "Foundation",
        "-framework",
        "Metal",
    ]
    t0 = now()
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True)
    out["build_seconds"] = now() - t0
    out["build_stdout"] = proc.stdout
    out["build_stderr"] = proc.stderr
    if proc.returncode != 0:
        raise RuntimeError(f"swiftc failed with code {proc.returncode}: {proc.stderr}")
    return out


def make_points(n: int, count: int, seed: int) -> Any:
    rng = np.random.default_rng(int(seed) ^ 0x122B17)
    raw = rng.standard_normal((int(count), int(n))) + 1j * rng.standard_normal((int(count), int(n)))
    norms = np.linalg.norm(raw, axis=1)
    norms = np.where(norms > 0, norms, 1.0)
    dirs = raw / norms[:, None] * math.sqrt(max(1, int(n)))
    radii = rng.uniform(0.05, 0.92, size=(int(count), 1))
    phases = np.exp(1j * rng.uniform(0.0, 2.0 * math.pi, size=(int(count), 1)))
    points = dirs * radii * phases
    return points.astype(np.complex128)


def write_metal_input(path: Path, system: Any, points: Any) -> None:
    coeff = np.asarray(system.coeff, dtype=np.complex64)
    pts = np.asarray(points, dtype=np.complex64)
    payload = {
        "n": int(system.n),
        "equations": int(system.n),
        "terms": int(system.terms_per_poly),
        "points": int(pts.shape[0]),
        "exps": [int(x) for x in np.asarray(system.exps, dtype=np.int32).reshape(-1).tolist()],
        "coeffRe": [float(x) for x in coeff.real.reshape(-1).tolist()],
        "coeffIm": [float(x) for x in coeff.imag.reshape(-1).tolist()],
        "pointsRe": [float(x) for x in pts.real.reshape(-1).tolist()],
        "pointsIm": [float(x) for x in pts.imag.reshape(-1).tolist()],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")


def benchmark_118_eval(system: Any, points: Any, loops: int) -> dict[str, Any]:
    checksum_re = 0.0
    checksum_im = 0.0
    checksum_norm = 0.0
    t0 = now()
    for _ in range(max(1, int(loops))):
        for z in points:
            f = system.eval(z)
            if checksum_norm == 0.0:
                checksum_re += float(np.sum(f.real))
                checksum_im += float(np.sum(f.imag))
                checksum_norm += float(np.sum(np.abs(f) ** 2))
    seconds = now() - t0
    point_evals = int(points.shape[0]) * int(system.n) * max(1, int(loops))
    return {
        "seconds": float(seconds),
        "point_evals": int(point_evals),
        "point_evals_per_second": float(point_evals / max(seconds, 1e-300)),
        "checksum_re": float(checksum_re),
        "checksum_im": float(checksum_im),
        "checksum_norm": float(checksum_norm),
        "eval_count": int(system.eval_count),
        "seconds_eval_reported": float(system.seconds_eval),
    }


def run_swift_metal(binary: Path, input_path: Path, loops: int) -> dict[str, Any]:
    t0 = now()
    proc = subprocess.run(
        [str(binary), str(input_path), str(max(1, int(loops)))],
        cwd=str(ROOT),
        text=True,
        capture_output=True,
    )
    wall = now() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"Swift Metal helper failed with code {proc.returncode}: {proc.stderr}")
    payload = json.loads(proc.stdout)
    payload["process_wall_seconds"] = float(wall)
    payload["stderr"] = proc.stderr
    return payload


def run_118_solver(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    t0 = now()
    result = ENGINE118.run_case(args, case_raw)
    result["_wall_seconds"] = now() - t0
    return result


def summarize_solver(result: dict[str, Any]) -> dict[str, Any]:
    summary = result.get("summary", {})
    eval_stats = summary.get("eval_stats", {})
    return {
        "requested_roots": int(summary.get("requested_roots", 0)),
        "unique_roots": int(summary.get("unique_roots", 0)),
        "success": bool(summary.get("success")),
        "trials_used": int(summary.get("trials_used", 0)),
        "seconds": float(result.get("_wall_seconds", summary.get("total_seconds", 0.0))),
        "eval_used": int(eval_stats.get("eval_count", 0)),
        "slope_calls": int(eval_stats.get("slope_count", 0)),
    }


def run_case(args: argparse.Namespace, case_raw: str, helper_build: dict[str, Any], binary: Path) -> dict[str, Any]:
    n, d = ENGINE118.parse_case(case_raw)
    if d <= 33:
        raise ValueError(f"case {case_raw!r} is not degree > 33")

    system_for_points = ENGINE118.DenseKostlanSystem.make(
        n,
        d,
        seed_index=int(args.seed_index),
        equation_normalize=bool(args.equation_normalize),
    )
    points = make_points(n, int(args.metal_points), int(system_for_points.seed))
    input_path = Path(args.outdir) / "inputs" / f"122_{n}x{d}_points{args.metal_points}.json"
    write_metal_input(input_path, system_for_points, points)

    system_for_numpy = ENGINE118.DenseKostlanSystem.make(
        n,
        d,
        seed_index=int(args.seed_index),
        equation_normalize=bool(args.equation_normalize),
    )
    numpy_eval = benchmark_118_eval(system_for_numpy, points, int(args.metal_loops))
    metal_eval = run_swift_metal(binary, input_path, int(args.metal_loops))

    solver_result = run_118_solver(args, case_raw) if bool(args.run_118_solver) else {}
    solver_summary = summarize_solver(solver_result) if solver_result else {}

    return {
        "case": case_raw,
        "n": int(n),
        "degree": int(d),
        "seed": int(system_for_points.seed),
        "terms_per_poly": int(system_for_points.terms_per_poly),
        "metal_points": int(args.metal_points),
        "metal_loops": int(args.metal_loops),
        "input": str(input_path),
        "helper_build": helper_build,
        "numpy_118_eval_loop": numpy_eval,
        "swift_metal_eval_batch": metal_eval,
        "eval_speedup_metal_vs_118_loop": float(
            metal_eval["point_evals_per_second"] / max(numpy_eval["point_evals_per_second"], 1e-300)
        ),
        "metal_kernel_speedup_vs_118_seconds": float(
            numpy_eval["seconds"] / max(float(metal_eval["kernel_seconds"]), 1e-300)
        ),
        "metal_process_speedup_vs_118_seconds": float(
            numpy_eval["seconds"] / max(float(metal_eval["process_wall_seconds"]), 1e-300)
        ),
        "full_118_solver": solver_summary,
    }


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    fields = [
        "case",
        "degree",
        "terms_per_poly",
        "points",
        "loops",
        "numpy_seconds",
        "metal_kernel_seconds",
        "metal_process_seconds",
        "numpy_point_evals_per_second",
        "metal_point_evals_per_second",
        "eval_speedup_metal_vs_118_loop",
        "metal_kernel_speedup_vs_118_seconds",
        "metal_process_speedup_vs_118_seconds",
        "solver_118_roots",
        "solver_118_requested",
        "solver_118_seconds",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            solver = row.get("full_118_solver", {})
            numpy_eval = row["numpy_118_eval_loop"]
            metal_eval = row["swift_metal_eval_batch"]
            writer.writerow(
                {
                    "case": row["case"],
                    "degree": row["degree"],
                    "terms_per_poly": row["terms_per_poly"],
                    "points": row["metal_points"],
                    "loops": row["metal_loops"],
                    "numpy_seconds": numpy_eval["seconds"],
                    "metal_kernel_seconds": metal_eval["kernel_seconds"],
                    "metal_process_seconds": metal_eval["process_wall_seconds"],
                    "numpy_point_evals_per_second": numpy_eval["point_evals_per_second"],
                    "metal_point_evals_per_second": metal_eval["point_evals_per_second"],
                    "eval_speedup_metal_vs_118_loop": row["eval_speedup_metal_vs_118_loop"],
                    "metal_kernel_speedup_vs_118_seconds": row["metal_kernel_speedup_vs_118_seconds"],
                    "metal_process_speedup_vs_118_seconds": row["metal_process_speedup_vs_118_seconds"],
                    "solver_118_roots": solver.get("unique_roots"),
                    "solver_118_requested": solver.get("requested_roots"),
                    "solver_118_seconds": solver.get("seconds"),
                }
            )


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    p = payload["parameters"]
    lines = [
        "# 122 Swift Metal vs 118 Benchmark",
        "",
        f"- cases: `{', '.join(payload['cases'])}`",
        f"- metal points/loops: `{p['metal_points']}/{p['metal_loops']}`",
        "- scope: Swift + Metal batch `F(z)` evaluator for the same 118 dense Kostlan system",
        "",
        "## Results",
        "",
        "| case | terms | 118 eval sec | Metal kernel sec | Metal process sec | kernel speedup | process speedup | full 118 roots | full 118 sec |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["runs"]:
        solver = row.get("full_118_solver", {})
        full_roots = ""
        if solver:
            full_roots = f"{solver.get('unique_roots')}/{solver.get('requested_roots')}"
        lines.append(
            f"| `{row['case']}` | {row['terms_per_poly']} | "
            f"{fmt_float(row['numpy_118_eval_loop']['seconds'], 4)} | "
            f"{fmt_float(row['swift_metal_eval_batch']['kernel_seconds'], 4)} | "
            f"{fmt_float(row['swift_metal_eval_batch']['process_wall_seconds'], 4)} | "
            f"{fmt_float(row['metal_kernel_speedup_vs_118_seconds'], 3)} | "
            f"{fmt_float(row['metal_process_speedup_vs_118_seconds'], 3)} | "
            f"{full_roots} | {fmt_float(solver.get('seconds'), 4)} |"
        )
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- 122 uses `float2` complex arithmetic in Metal, so it is a batch evaluator prototype, not the final `complex128` corrector.",
            "- The useful number is the process speedup when input/output and command dispatch are included.",
            "- A complete Swift/Metal solver would need to keep probe scoring and candidate selection on the GPU, then refine accepted candidates in CPU `complex128`.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    p = ENGINE118.build_parser()
    p.description = "122 benchmark: Swift + Metal batch evaluator against 118 NumPy eval loop."
    p.set_defaults(
        cases="2,34",
        count=8,
        pool=512,
        epochs=24,
        outdir="verification/122_swift_metal_vs_118_benchmark",
    )
    p.add_argument("--metal-points", type=int, default=8192)
    p.add_argument("--metal-loops", type=int, default=4)
    p.add_argument("--swift-source", default=str(FLOW / "122_swift_metal_eval.swift"))
    p.add_argument("--swift-bin", default=None)
    p.add_argument("--rebuild-swift", action="store_true")
    p.add_argument("--run-118-solver", action=argparse.BooleanOptionalAction, default=True)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ENGINE118.ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)
    cases = parse_cases(args.cases)
    out_json = Path(args.out) if args.out else Path(args.outdir) / "122_swift_metal_vs_118_benchmark.json"
    binary = Path(args.swift_bin) if args.swift_bin else Path(args.outdir) / "bin" / "122_swift_metal_eval"
    helper_build = build_swift_helper(Path(args.swift_source), binary, force=bool(args.rebuild_swift))

    t0 = now()
    runs = []
    print("=" * 120, flush=True)
    print("122 benchmark: Swift + Metal batch evaluator vs 118 NumPy evaluator", flush=True)
    print(f"cases={cases} metal_points={args.metal_points} metal_loops={args.metal_loops}", flush=True)
    print("=" * 120, flush=True)
    for case_raw in cases:
        print(f"case={case_raw}", flush=True)
        row = run_case(args, case_raw, helper_build, binary)
        runs.append(row)
        solver = row.get("full_118_solver", {})
        print(
            f"  118 eval loop: seconds={row['numpy_118_eval_loop']['seconds']:.4f} "
            f"eval/s={row['numpy_118_eval_loop']['point_evals_per_second']:.1f}",
            flush=True,
        )
        print(
            f"  122 Swift+Metal: kernel={row['swift_metal_eval_batch']['kernel_seconds']:.4f}s "
            f"process={row['swift_metal_eval_batch']['process_wall_seconds']:.4f}s "
            f"eval/s={row['swift_metal_eval_batch']['point_evals_per_second']:.1f}",
            flush=True,
        )
        print(
            f"  speedup kernel={row['metal_kernel_speedup_vs_118_seconds']:.3f} "
            f"process={row['metal_process_speedup_vs_118_seconds']:.3f}",
            flush=True,
        )
        if solver:
            print(
                f"  full 118 context: roots={solver.get('unique_roots')}/{solver.get('requested_roots')} "
                f"seconds={float(solver.get('seconds') or 0.0):.4f}",
                flush=True,
            )

    payload = {
        "script": "122_pandrosion_swift_metal_vs_118_benchmark.py",
        "cases": cases,
        "parameters": {
            "metal_points": int(args.metal_points),
            "metal_loops": int(args.metal_loops),
            "seed_index": int(args.seed_index),
            "equation_normalize": bool(args.equation_normalize),
            "run_118_solver": bool(args.run_118_solver),
            "count": int(args.count),
            "pool": int(args.pool),
            "epochs": int(args.epochs),
            "accept": float(args.accept),
        },
        "seconds": float(now() - t0),
        "runs": runs,
    }
    out_csv = out_json.with_suffix(".csv")
    out_md = out_json.with_suffix(".md")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_csv(out_csv, runs)
    write_markdown(out_md, payload)

    print("=" * 120, flush=True)
    print(f"json={out_json}", flush=True)
    print(f"csv={out_csv}", flush=True)
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
