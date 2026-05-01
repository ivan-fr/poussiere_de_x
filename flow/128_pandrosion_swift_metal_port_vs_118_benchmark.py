#!/usr/bin/env python3
"""
128_pandrosion_swift_metal_port_vs_118_benchmark.py

Benchmark the first end-to-end Swift-side port of the 118 solve loop:

  - Python 118 still generates the exact dense Kostlan system so the comparison
    is against the same coefficients.
  - Swift + Metal does batch start generation/selection.
  - Swift + Metal does the n=2 pre-polish.
  - Swift Double CPU does the final pure Pandrosion finite-slope correction and
    clustering.

This is intentionally a port milestone, not yet a fully independent Swift
Kostlan generator.
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
    spec = importlib.util.spec_from_file_location("pandrosion_118_for_128", path)
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


def load_system(proc: Any, system: Any) -> None:
    if proc.stdin is None or proc.stdout is None:
        raise RuntimeError("Swift server has no stdin/stdout pipe")
    coeff64 = np.asarray(system.coeff, dtype=np.complex128)
    coeff32 = np.asarray(system.coeff, dtype=np.complex64)
    payload = {
        "op": "load",
        "n": int(system.n),
        "equations": int(system.n),
        "terms": int(system.terms_per_poly),
        "exps": [int(x) for x in np.asarray(system.exps, dtype=np.int32).reshape(-1).tolist()],
        "coeffRe": [float(x) for x in coeff32.real.reshape(-1).tolist()],
        "coeffIm": [float(x) for x in coeff32.imag.reshape(-1).tolist()],
        "coeffRe64": [float(x) for x in coeff64.real.reshape(-1).tolist()],
        "coeffIm64": [float(x) for x in coeff64.imag.reshape(-1).tolist()],
    }
    proc.stdin.write(json.dumps(payload, separators=(",", ":")) + "\n")
    proc.stdin.flush()
    response = proc.stdout.readline()
    if not response:
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        raise RuntimeError(f"Swift server load failed: {stderr}")
    result = json.loads(response)
    if result.get("error") or not result.get("ok"):
        raise RuntimeError(f"Swift server load error: {result}")


def solve_payload(args: argparse.Namespace, system: Any) -> dict[str, Any]:
    powers = sorted(set(round(float(x), 16) for x in ENGINE118.parse_float_list(args.powers, ENGINE118.DEFAULT_POWERS, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    angles = [math.radians(x) for x in ENGINE118.parse_float_list(args.angles, ENGINE118.DEFAULT_ANGLES_DEG)]
    radii = ENGINE118.parse_float_list(args.rays, ENGINE118.DEFAULT_RADII, positive=True)
    probe_radii = ENGINE118.parse_float_list(args.probe_radii, [0.0, 0.35, 0.7, 1.0, 1.6, 2.6, 4.2], positive=False)
    probe_radii = [r for r in probe_radii if r >= 0] or [0.0, 1.0]
    gains = ENGINE118.parse_float_list(args.startopt_gains, ENGINE118.DEFAULT_GAINS, positive=True)
    selector_top = int(args.selector_top) if int(args.selector_top) > 0 else int(args.refine_top) * 8
    return {
        "op": "solve2",
        "candidates": int(args.metal_candidates),
        "targetCount": int(args.count),
        "topK": int(selector_top),
        "refineTop": int(args.refine_top),
        "count": int(args.count),
        "seed": int(system.seed) + 0x113000,
        "powerCap": float(args.power_cap),
        "powers": [float(x) for x in powers],
        "angles": [float(x) for x in angles],
        "radii": [float(x) for x in radii],
        "diversitySep": float(args.diversity_sep),
        "epochs": int(args.epochs),
        "tol": float(args.tol),
        "accept": float(args.accept),
        "clusterSep": float(args.cluster_sep),
        "lineSearch": int(args.line_search),
        "probeScale": float(args.probe_scale),
        "probeCandidates": int(args.probe_candidates),
        "probeRadii": [float(x) for x in probe_radii],
        "probeSelf": bool(args.probe_self),
        "startoptSteps": int(args.startopt_steps),
        "startoptCandidates": int(args.startopt_candidates),
        "startoptGains": [float(x) for x in gains],
        "startoptMicroEpochs": int(args.startopt_micro_epochs),
        "metalPolish2": bool(args.metal_polish2),
        "metalPolishEpochs": int(args.metal_polish_epochs),
        "metalPolishProbes": int(args.metal_polish_probes),
        "metalPolishLineSearch": int(args.metal_polish_line_search),
        "metalPolishProbeScale": float(args.metal_polish_probe_scale),
    }


def run_swift_solve(proc: Any, payload: dict[str, Any]) -> dict[str, Any]:
    if proc.stdin is None or proc.stdout is None:
        raise RuntimeError("Swift server has no stdin/stdout pipe")
    t0 = now()
    proc.stdin.write(json.dumps(payload, separators=(",", ":")) + "\n")
    proc.stdin.flush()
    response = proc.stdout.readline()
    wall = now() - t0
    if not response:
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        raise RuntimeError(f"Swift server solve failed: {stderr}")
    result = json.loads(response)
    if result.get("error"):
        raise RuntimeError(f"Swift server solve error: {result['error']}")
    result["_wall_seconds"] = float(wall)
    return result


def summarize_result(engine: str, case_raw: str, result: dict[str, Any], seconds: float) -> dict[str, Any]:
    summary = result.get("summary", {})
    eval_stats = summary.get("eval_stats", {})
    residuals = [float(r.get("residual", float("nan"))) for r in result.get("roots", [])]
    finite_residuals = [r for r in residuals if math.isfinite(r)]
    return {
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
        "refine_seconds": float(summary.get("refine_seconds", summary.get("extract_seconds", 0.0))),
        "metal_select_seconds": float(summary.get("metal_select_process_seconds", 0.0)),
        "metal_polish2_seconds": float(summary.get("metal_polish2_seconds", 0.0)),
        "eval_used": int(eval_stats.get("eval_count", 0)) if eval_stats.get("eval_count") is not None else None,
        "slope_calls": int(eval_stats.get("slope_count", 0)) if eval_stats.get("slope_count") is not None else None,
        "seconds_eval": float(eval_stats.get("seconds_eval", 0.0)),
        "seconds_slope": float(eval_stats.get("seconds_slope", 0.0)),
        "residual_min": float(min(finite_residuals)) if finite_residuals else None,
        "residual_mean": float(sum(finite_residuals) / len(finite_residuals)) if finite_residuals else None,
        "residual_max": float(max(finite_residuals)) if finite_residuals else None,
        "error": None,
    }


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
        "refine_seconds": None,
        "metal_select_seconds": None,
        "metal_polish2_seconds": None,
        "eval_used": None,
        "slope_calls": None,
        "seconds_eval": None,
        "seconds_slope": None,
        "residual_min": None,
        "residual_mean": None,
        "residual_max": None,
        "error": f"{type(exc).__name__}: {exc}",
    }


def compare_pair(row118: dict[str, Any], row128: dict[str, Any]) -> dict[str, Any]:
    seconds118 = float(row118.get("seconds") or 0.0)
    seconds128 = float(row128.get("seconds") or 0.0)
    roots118 = int(row118.get("unique_roots") or 0)
    roots128 = int(row128.get("unique_roots") or 0)
    if roots128 > roots118:
        winner = "128(swift)"
    elif roots118 > roots128:
        winner = "118"
    elif bool(row128.get("success")) and not bool(row118.get("success")):
        winner = "128(swift)"
    elif bool(row118.get("success")) and not bool(row128.get("success")):
        winner = "118"
    elif seconds118 > 0 and seconds128 > 0:
        winner = "128(swift)" if seconds128 < seconds118 else "118"
    else:
        winner = "tie"
    return {
        "case": row118.get("case"),
        "roots_118": roots118,
        "roots_128": roots128,
        "root_delta_128_minus_118": roots128 - roots118,
        "success_118": bool(row118.get("success")),
        "success_128": bool(row128.get("success")),
        "seconds_118": seconds118,
        "seconds_128": seconds128,
        "speedup_128_vs_118": seconds118 / seconds128 if seconds128 > 0 else None,
        "evals_118": row118.get("eval_used"),
        "evals_128": row128.get("eval_used"),
        "slopes_118": row118.get("slope_calls"),
        "slopes_128": row128.get("slope_calls"),
        "winner": winner,
    }


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
        "refine_seconds",
        "metal_select_seconds",
        "metal_polish2_seconds",
        "eval_used",
        "slope_calls",
        "seconds_eval",
        "seconds_slope",
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


def write_pair_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    fields = [
        "case",
        "roots_118",
        "roots_128",
        "root_delta_128_minus_118",
        "success_118",
        "success_128",
        "seconds_118",
        "seconds_128",
        "speedup_128_vs_118",
        "evals_118",
        "evals_128",
        "slopes_118",
        "slopes_128",
        "winner",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    s = payload["summary"]
    p = payload["parameters"]
    lines = [
        "# 128 Swift Metal Port vs 118 Benchmark",
        "",
        f"- cases: `{', '.join(payload['cases'])}`",
        f"- metal candidates/selector top/refine top: `{p['metal_candidates']}/{p['selector_top']}/{p['refine_top']}`",
        f"- count/pool/epochs: `{p['count']}/{p['pool']}/{p['epochs']}`",
        "- scope: same 118 Kostlan system; Metal selector; Metal n=2 polish; Swift Double final Pandrosion refine",
        "",
        "## Summary",
        "",
        f"- 118 roots: `{s['roots_118']}/{s['requested_118']}`; complete cases: `{s['successes_118']}/{s['pairs']}`",
        f"- 128 roots: `{s['roots_128']}/{s['requested_128']}`; complete cases: `{s['successes_128']}/{s['pairs']}`",
        f"- wall seconds: 118=`{s['seconds_118']:.4f}`, 128=`{s['seconds_128']:.4f}`, total=`{s['seconds']:.4f}`",
        f"- root delta 128-118: `{s['root_delta_128_minus_118']}`",
        "",
        "## Pair Results",
        "",
        "| case | 118 roots | 128 roots | delta | 118 sec | 128 sec | 128 speedup | 118 evals | 128 evals | winner |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in payload["pairs"]:
        lines.append(
            f"| `{row.get('case')}` | {row.get('roots_118')} | {row.get('roots_128')} | "
            f"{row.get('root_delta_128_minus_118')} | {fmt_float(row.get('seconds_118'), 4)} | "
            f"{fmt_float(row.get('seconds_128'), 4)} | {fmt_float(row.get('speedup_128_vs_118'), 3)} | "
            f"{row.get('evals_118') or ''} | {row.get('evals_128') or ''} | `{row.get('winner')}` |"
        )
    lines.extend(
        [
            "",
            "## Timing Breakdown",
            "",
            "| case | engine | metal select | metal polish2 | refine/extract | eval sec | slope sec | max residual |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in payload["runs"]:
        lines.append(
            f"| `{row.get('case')}` | `{row.get('engine')}` | {fmt_float(row.get('metal_select_seconds'), 4)} | "
            f"{fmt_float(row.get('metal_polish2_seconds'), 4)} | {fmt_float(row.get('refine_seconds'), 4)} | "
            f"{fmt_float(row.get('seconds_eval'), 4)} | {fmt_float(row.get('seconds_slope'), 4)} | "
            f"{fmt_sci(row.get('residual_max'))} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    p = ENGINE118.build_parser()
    p.description = "128 benchmark: Swift+Metal selector/polish and Swift Double 118-style refine vs Python 118."
    p.set_defaults(
        cases="2,34",
        count=8,
        pool=512,
        epochs=24,
        outdir="verification/128_swift_metal_port_vs_118_benchmark",
    )
    p.add_argument("--metal-candidates", type=int, default=2048)
    p.add_argument("--selector-top", type=int, default=0, help="0 means 8 * --refine-top")
    p.add_argument("--refine-top", type=int, default=96)
    p.add_argument("--diversity-sep", type=float, default=0.45)
    p.add_argument("--metal-polish2", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--metal-polish-epochs", type=int, default=4)
    p.add_argument("--metal-polish-probes", type=int, default=8)
    p.add_argument("--metal-polish-line-search", type=int, default=6)
    p.add_argument("--metal-polish-probe-scale", type=float, default=0.035)
    p.add_argument("--swift-source", default=str(FLOW / "125_swift_metal_start_select_server.swift"))
    p.add_argument("--swift-bin", default=None)
    p.add_argument("--rebuild-swift", action="store_true")
    p.add_argument("--run-order", choices=["118-first", "128-first"], default="118-first")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ENGINE118.ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)
    cases = parse_cases(args.cases)
    out_json = Path(args.out) if args.out else Path(args.outdir) / "128_swift_metal_port_vs_118_benchmark.json"
    swift_source = Path(args.swift_source)
    binary = Path(args.swift_bin) if args.swift_bin else Path(args.outdir) / "bin" / "125_swift_metal_start_select_server"
    helper_build = build_swift_helper(swift_source, binary, force=bool(args.rebuild_swift))

    server_proc = subprocess.Popen(
        [str(binary)],
        cwd=str(ROOT),
        text=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    rows: list[dict[str, Any]] = []
    pairs: list[dict[str, Any]] = []
    full_results: dict[str, Any] = {}
    t_all = now()

    print("=" * 120, flush=True)
    print("128 benchmark: Swift+Metal selector/polish + Swift Double refine vs 118", flush=True)
    print(
        f"cases={cases} count={args.count} pool={args.pool} epochs={args.epochs} "
        f"metal_candidates={args.metal_candidates} selector_top={args.selector_top or int(args.refine_top) * 8} "
        f"refine_top={args.refine_top}",
        flush=True,
    )
    print("=" * 120, flush=True)

    try:
        for case_raw in cases:
            n, d = ENGINE118.parse_case(case_raw)
            if d <= 33:
                raise ValueError(f"case {case_raw!r} is not degree > 33")
            print(f"case={case_raw}", flush=True)
            system = ENGINE118.DenseKostlanSystem.make(
                n,
                d,
                seed_index=int(args.seed_index),
                equation_normalize=bool(args.equation_normalize),
            )
            load_system(server_proc, system)
            payload = solve_payload(args, system)

            if args.run_order == "128-first":
                t0 = now()
                result128 = run_swift_solve(server_proc, payload)
                row128 = summarize_result("128(swift)", case_raw, result128, now() - t0)
                t0 = now()
                result118 = ENGINE118.run_case(args, case_raw)
                row118 = summarize_result("118", case_raw, result118, now() - t0)
            else:
                t0 = now()
                result118 = ENGINE118.run_case(args, case_raw)
                row118 = summarize_result("118", case_raw, result118, now() - t0)
                t0 = now()
                result128 = run_swift_solve(server_proc, payload)
                row128 = summarize_result("128(swift)", case_raw, result128, now() - t0)

            rows.extend([row118, row128])
            print(
                f"  118: status={row118.get('status')} roots={row118.get('unique_roots')}/{row118.get('requested_roots')} "
                f"evals={row118.get('eval_used')} slopes={row118.get('slope_calls')} seconds={float(row118.get('seconds') or 0.0):.4f}",
                flush=True,
            )
            print(
                f"  128(swift): status={row128.get('status')} roots={row128.get('unique_roots')}/{row128.get('requested_roots')} "
                f"evals={row128.get('eval_used')} slopes={row128.get('slope_calls')} seconds={float(row128.get('seconds') or 0.0):.4f}",
                flush=True,
            )
            pair = compare_pair(row118, row128)
            pairs.append(pair)
            print(
                f"  winner={pair['winner']} delta128-118={pair['root_delta_128_minus_118']} "
                f"speedup128={fmt_float(pair.get('speedup_128_vs_118'), 3)}",
                flush=True,
            )
            full_results[case_raw.replace(",", "x")] = {"118": result118, "128(swift)": result128}
    finally:
        try:
            if server_proc.stdin is not None:
                server_proc.stdin.write("__quit__\n")
                server_proc.stdin.flush()
        except Exception:
            pass
        try:
            server_proc.wait(timeout=5)
        except Exception:
            server_proc.kill()

    seconds_all = now() - t_all
    summary = {
        "pairs": len(pairs),
        "runs": len(rows),
        "successes_118": int(sum(1 for row in rows if row.get("engine") == "118" and row.get("success") is True)),
        "successes_128": int(sum(1 for row in rows if row.get("engine") == "128(swift)" and row.get("success") is True)),
        "roots_118": int(sum(int(row.get("unique_roots") or 0) for row in rows if row.get("engine") == "118")),
        "roots_128": int(sum(int(row.get("unique_roots") or 0) for row in rows if row.get("engine") == "128(swift)")),
        "requested_118": int(sum(int(row.get("requested_roots") or 0) for row in rows if row.get("engine") == "118")),
        "requested_128": int(sum(int(row.get("requested_roots") or 0) for row in rows if row.get("engine") == "128(swift)")),
        "seconds_118": float(sum(float(row.get("seconds") or 0.0) for row in rows if row.get("engine") == "118")),
        "seconds_128": float(sum(float(row.get("seconds") or 0.0) for row in rows if row.get("engine") == "128(swift)")),
        "root_delta_128_minus_118": int(sum(int(pair.get("root_delta_128_minus_118") or 0) for pair in pairs)),
        "wins_118": int(sum(1 for pair in pairs if pair.get("winner") == "118")),
        "wins_128": int(sum(1 for pair in pairs if pair.get("winner") == "128(swift)")),
        "ties": int(sum(1 for pair in pairs if pair.get("winner") == "tie")),
        "seconds": float(seconds_all),
    }
    payload_out = {
        "script": "128_pandrosion_swift_metal_port_vs_118_benchmark.py",
        "engines": {
            "118": "flow/118_pandrosion_probe_aware_pure_thales_engine.py",
            "128(swift)": "flow/125_swift_metal_start_select_server.swift op=solve2",
        },
        "helper_build": helper_build,
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
            "seed_index": int(args.seed_index),
            "metal_candidates": int(args.metal_candidates),
            "selector_top": int(args.selector_top) if int(args.selector_top) > 0 else int(args.refine_top) * 8,
            "refine_top": int(args.refine_top),
            "diversity_sep": float(args.diversity_sep),
            "metal_polish2": bool(args.metal_polish2),
            "metal_polish_epochs": int(args.metal_polish_epochs),
            "metal_polish_probes": int(args.metal_polish_probes),
            "metal_polish_line_search": int(args.metal_polish_line_search),
            "metal_polish_probe_scale": float(args.metal_polish_probe_scale),
            "run_order": str(args.run_order),
        },
        "summary": summary,
        "pairs": pairs,
        "runs": rows,
        "full_results": full_results if bool(args.verbose_trials) else {},
    }
    out_csv = out_json.with_suffix(".csv")
    out_pairs_csv = out_json.with_name(out_json.stem + "_pairs.csv")
    out_md = out_json.with_suffix(".md")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload_out, indent=2), encoding="utf-8")
    write_csv(out_csv, rows)
    write_pair_csv(out_pairs_csv, pairs)
    write_markdown(out_md, payload_out)

    print("=" * 120, flush=True)
    print(f"summary: 118 roots={summary['roots_118']}/{summary['requested_118']} successes={summary['successes_118']}/{summary['pairs']} seconds={summary['seconds_118']:.4f}", flush=True)
    print(f"summary: 128 roots={summary['roots_128']}/{summary['requested_128']} successes={summary['successes_128']}/{summary['pairs']} seconds={summary['seconds_128']:.4f}", flush=True)
    print(f"root_delta_128_minus_118={summary['root_delta_128_minus_118']}", flush=True)
    print(f"json={out_json}", flush=True)
    print(f"csv={out_csv}", flush=True)
    print(f"pairs_csv={out_pairs_csv}", flush=True)
    print(f"report={out_md}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
