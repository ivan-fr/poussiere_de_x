#!/usr/bin/env python3
"""
120_pandrosion_probe_vs_vectorized_benchmark.py

Benchmark harness for comparing:
  - 116: multifamily harness over the 115 vectorized PURE Pandrosion engine
  - 119: multifamily harness over the 118 probe-aware PURE Thales engine

The comparison uses the same cases, families, root target, trial pool, epoch
budget, tolerances, and generator options.  By default it also forces 119 to use
the same multifamily seeds as 116, so both engines see the same generated
polynomial systems for every (case, family).
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


def _load_module(filename: str, module_name: str) -> Any:
    path = Path(__file__).resolve().with_name(filename)
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


ENGINE116 = _load_module("116_pandrosion_multifamily_vectorized_pure_pandrosion.py", "pandrosion_116_for_120")
ENGINE119 = _load_module("119_pandrosion_multifamily_probe_aware_pure_thales_engine.py", "pandrosion_119_for_120")
np = ENGINE116.np


def now() -> float:
    return time.perf_counter()


def ensure_numpy() -> None:
    ENGINE116.ensure_numpy()
    ENGINE119.ensure_numpy()


def parse_cases(raw: str) -> list[str]:
    return [c.strip() for c in str(raw).replace("|", ";").split(";") if c.strip()]


def residual_stats(values: Sequence[float]) -> dict[str, Optional[float]]:
    finite = [float(v) for v in values if math.isfinite(float(v))]
    if not finite:
        return {"residual_min": None, "residual_mean": None, "residual_max": None}
    return {
        "residual_min": float(min(finite)),
        "residual_mean": float(sum(finite) / len(finite)),
        "residual_max": float(max(finite)),
    }


def make_engine_args(module: Any, args: argparse.Namespace) -> argparse.Namespace:
    ns = module.build_parser().parse_args([])
    common_names = [
        "seed_index",
        "equation_normalize",
        "linear_scale",
        "count",
        "pool",
        "epochs",
        "tol",
        "accept",
        "cluster_sep",
        "trial_timeout",
        "line_search",
        "probe_scale",
        "powers",
        "power_cap",
        "angles",
        "rays",
        "startopt_steps",
        "startopt_candidates",
        "startopt_gains",
        "startopt_micro_epochs",
        "families",
        "sparse_terms",
        "sparse_frac",
    ]
    for name in common_names:
        if hasattr(ns, name):
            setattr(ns, name, getattr(args, name))
    if hasattr(ns, "probe_candidates"):
        ns.probe_candidates = int(args.probe_candidates)
    if hasattr(ns, "probe_radii"):
        ns.probe_radii = str(args.probe_radii)
    if hasattr(ns, "probe_self"):
        ns.probe_self = bool(args.probe_self)
    if hasattr(ns, "keep_trials"):
        keep = int(args.keep_trials) if int(args.keep_trials) > 0 else int(args.pool)
        ns.keep_trials = max(keep, int(args.pool))
    if hasattr(ns, "verbose_trials"):
        ns.verbose_trials = bool(args.verbose_trials)
    ns.out = None
    ns.outdir = str(args.outdir)
    ns.cases = str(args.cases)
    return ns


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


def summarize_result(engine: str, case_raw: str, family: str, result: dict[str, Any], seconds: float) -> dict[str, Any]:
    summary = result.get("summary", {})
    eval_stats = summary.get("eval_stats", {})
    residuals = [float(r.get("residual", float("nan"))) for r in result.get("roots", [])]
    out = {
        "case": case_raw,
        "family": family,
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
        "probe_eval_sum": probe_eval_sum(result),
        "seed": int(result.get("seed", 0)),
        "terms": int(result.get("terms", 0)),
        "active_terms": int(result.get("active_terms", result.get("terms", 0))),
        "bezout": int(result.get("bezout", 0)),
        "error": None,
    }
    out.update(residual_stats(residuals))
    return out


def error_row(engine: str, case_raw: str, family: str, exc: BaseException, seconds: float) -> dict[str, Any]:
    return {
        "case": case_raw,
        "family": family,
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
        "probe_eval_sum": None,
        "seed": None,
        "terms": None,
        "active_terms": None,
        "bezout": None,
        "residual_min": None,
        "residual_mean": None,
        "residual_max": None,
        "error": f"{type(exc).__name__}: {exc}",
    }


def run_engine(module: Any, engine: str, engine_args: argparse.Namespace, case_raw: str, family: str) -> tuple[dict[str, Any], dict[str, Any]]:
    t0 = now()
    try:
        result = module.run_case(engine_args, case_raw, family)
    except Exception as exc:
        return error_row(engine, case_raw, family, exc, now() - t0), {}
    return summarize_result(engine, case_raw, family, result, now() - t0), result


def compare_pair(row116: dict[str, Any], row119: dict[str, Any]) -> dict[str, Any]:
    roots116 = int(row116.get("unique_roots") or 0)
    roots119 = int(row119.get("unique_roots") or 0)
    seconds116 = float(row116.get("seconds") or 0.0)
    seconds119 = float(row119.get("seconds") or 0.0)
    if roots119 > roots116:
        winner = "119(118)"
    elif roots116 > roots119:
        winner = "116(115)"
    elif bool(row119.get("success")) and not bool(row116.get("success")):
        winner = "119(118)"
    elif bool(row116.get("success")) and not bool(row119.get("success")):
        winner = "116(115)"
    elif seconds119 > 0 and seconds116 > 0:
        winner = "119(118)" if seconds119 < seconds116 else "116(115)"
    else:
        winner = "tie"
    return {
        "case": row116.get("case"),
        "family": row116.get("family"),
        "roots_116": roots116,
        "roots_119": roots119,
        "root_delta_119_minus_116": roots119 - roots116,
        "success_116": bool(row116.get("success")),
        "success_119": bool(row119.get("success")),
        "seconds_116": seconds116,
        "seconds_119": seconds119,
        "speedup_119_vs_116": (seconds116 / seconds119 if seconds119 > 0 else None),
        "evals_116": row116.get("eval_used"),
        "evals_119": row119.get("eval_used"),
        "slopes_116": row116.get("slope_calls"),
        "slopes_119": row119.get("slope_calls"),
        "residual_min_116": row116.get("residual_min"),
        "residual_min_119": row119.get("residual_min"),
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
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "case",
        "family",
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
        "probe_eval_sum",
        "seed",
        "terms",
        "active_terms",
        "bezout",
        "residual_min",
        "residual_mean",
        "residual_max",
        "error",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_pair_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "case",
        "family",
        "roots_116",
        "roots_119",
        "root_delta_119_minus_116",
        "success_116",
        "success_119",
        "seconds_116",
        "seconds_119",
        "speedup_119_vs_116",
        "evals_116",
        "evals_119",
        "slopes_116",
        "slopes_119",
        "residual_min_116",
        "residual_min_119",
        "winner",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    p = payload["parameters"]
    lines = [
        "# 120 Pandrosion Probe-Aware vs Vectorized Benchmark",
        "",
        f"- cases: `{', '.join(payload['cases'])}`",
        f"- families: `{', '.join(payload['families'])}`",
        f"- count/pool/epochs: `{p['count']}/{p['pool']}/{p['epochs']}`",
        f"- same generated systems: `{p['same_systems']}`",
        "",
        "## Summary",
        "",
        f"- 116(115) roots: `{payload['summary']['roots_116']}/{payload['summary']['requested_116']}`; complete families: `{payload['summary']['successes_116']}/{payload['summary']['pairs']}`",
        f"- 119(118) roots: `{payload['summary']['roots_119']}/{payload['summary']['requested_119']}`; complete families: `{payload['summary']['successes_119']}/{payload['summary']['pairs']}`",
        f"- wall seconds: 116(115)=`{payload['summary']['seconds_116']:.4f}`, 119(118)=`{payload['summary']['seconds_119']:.4f}`, total=`{payload['summary']['seconds']:.4f}`",
        f"- root delta 119-116: `{payload['summary']['root_delta_119_minus_116']}`",
        "",
        "## Pair Results",
        "",
        "| case | family | 116 roots | 119 roots | delta | 116 sec | 119 sec | 119 speedup | 116 evals | 119 evals | winner |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in payload["pairs"]:
        speed = fmt_float(row.get("speedup_119_vs_116"), 3)
        lines.append(
            f"| `{row.get('case')}` | `{row.get('family')}` | "
            f"{row.get('roots_116')} | {row.get('roots_119')} | {row.get('root_delta_119_minus_116')} | "
            f"{fmt_float(row.get('seconds_116'), 4)} | {fmt_float(row.get('seconds_119'), 4)} | {speed} | "
            f"{row.get('evals_116') or ''} | {row.get('evals_119') or ''} | `{row.get('winner')}` |"
        )
    lines.extend([
        "",
        "## Residual Minima",
        "",
        "| case | family | 116 min residual | 119 min residual |",
        "|---|---|---:|---:|",
    ])
    for row in payload["pairs"]:
        lines.append(
            f"| `{row.get('case')}` | `{row.get('family')}` | "
            f"{fmt_sci(row.get('residual_min_116'))} | {fmt_sci(row.get('residual_min_119'))} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="120 benchmark: compare 119(118 probe-aware) against 116(115 vectorized pure).")
    p.add_argument("--cases", default="3,4", help="case n,d; multiple cases separated by ';'")
    p.add_argument("--families", default="all", help="family list or group inherited from 116/119")
    p.add_argument("--list-families", action="store_true")
    p.add_argument("--seed-index", type=int, default=0)
    p.add_argument("--same-systems", action=argparse.BooleanOptionalAction, default=True, help="force 119 to use the same multifamily seeds as 116")
    p.add_argument("--equation-normalize", action="store_true", default=False)
    p.add_argument("--no-equation-normalize", dest="equation_normalize", action="store_false")
    p.add_argument("--linear-scale", type=float, default=1.0)
    p.add_argument("--count", type=int, default=4)
    p.add_argument("--pool", type=int, default=512)
    p.add_argument("--epochs", type=int, default=24)
    p.add_argument("--tol", type=float, default=1e-12)
    p.add_argument("--accept", type=float, default=1e-8)
    p.add_argument("--cluster-sep", type=float, default=1e-8)
    p.add_argument("--trial-timeout", type=float, default=0.0)
    p.add_argument("--line-search", type=int, default=12)
    p.add_argument("--probe-scale", type=float, default=0.035)
    p.add_argument("--probe-candidates", type=int, default=8)
    p.add_argument("--probe-radii", default="0,0.35,0.7,1,1.6,2.6,4.2")
    p.add_argument("--probe-self", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--powers", default=None)
    p.add_argument("--power-cap", type=float, default=1048576.0)
    p.add_argument("--angles", default=None)
    p.add_argument("--rays", default=None)
    p.add_argument("--startopt-steps", type=int, default=1)
    p.add_argument("--startopt-candidates", type=int, default=12)
    p.add_argument("--startopt-gains", default=None)
    p.add_argument("--startopt-micro-epochs", type=int, default=0)
    p.add_argument("--sparse-terms", type=int, default=0)
    p.add_argument("--sparse-frac", type=float, default=0.18)
    p.add_argument("--keep-trials", type=int, default=0, help="0 means keep at least --pool trial summaries")
    p.add_argument("--verbose-trials", action="store_true")
    p.add_argument("--outdir", default="verification/120_probe_vs_vectorized_benchmark")
    p.add_argument("--out", default=None)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ensure_numpy()
    parser = build_parser()
    args = parser.parse_args(argv)
    if bool(args.list_families):
        print("120 benchmark families")
        for name in ENGINE116.DEFAULT_FAMILIES:
            print(f"{name}: {ENGINE116.FAMILY_DESCRIPTIONS[name]}")
        print("groups: " + ", ".join(f"{k}={','.join(v)}" for k, v in ENGINE116.FAMILY_GROUPS.items()))
        return 0

    if bool(args.same_systems):
        ENGINE119._family_seed = ENGINE116._family_seed

    cases = parse_cases(args.cases)
    families = ENGINE116.parse_families(args.families)
    args116 = make_engine_args(ENGINE116, args)
    args119 = make_engine_args(ENGINE119, args)

    rows: list[dict[str, Any]] = []
    pairs: list[dict[str, Any]] = []
    full_results: dict[str, Any] = {}
    t_all = now()

    print("=" * 120, flush=True)
    print("120 benchmark: 119(118 probe-aware) vs 116(115 vectorized pure)", flush=True)
    print(f"cases={cases} families={families} count={args.count} pool={args.pool} epochs={args.epochs} same_systems={args.same_systems}", flush=True)
    print("=" * 120, flush=True)

    for case_raw in cases:
        for family in families:
            print(f"case={case_raw} family={family}", flush=True)
            row116, result116 = run_engine(ENGINE116, "116(115)", args116, case_raw, family)
            rows.append(row116)
            roots116 = row116.get("unique_roots")
            req116 = row116.get("requested_roots")
            print(
                f"  116(115): status={row116.get('status')} roots={roots116}/{req116} "
                f"evals={row116.get('eval_used')} slopes={row116.get('slope_calls')} seconds={float(row116.get('seconds') or 0.0):.4f}",
                flush=True,
            )

            row119, result119 = run_engine(ENGINE119, "119(118)", args119, case_raw, family)
            rows.append(row119)
            roots119 = row119.get("unique_roots")
            req119 = row119.get("requested_roots")
            print(
                f"  119(118): status={row119.get('status')} roots={roots119}/{req119} "
                f"evals={row119.get('eval_used')} slopes={row119.get('slope_calls')} seconds={float(row119.get('seconds') or 0.0):.4f}",
                flush=True,
            )

            pair = compare_pair(row116, row119)
            pairs.append(pair)
            print(
                f"  winner={pair['winner']} delta119-116={pair['root_delta_119_minus_116']} "
                f"speedup119={fmt_float(pair.get('speedup_119_vs_116'), 3)}",
                flush=True,
            )
            key = f"{case_raw.replace(',', 'x')}_{family}"
            full_results[key] = {"116(115)": result116, "119(118)": result119}

    seconds_all = now() - t_all
    summary = {
        "pairs": len(pairs),
        "runs": len(rows),
        "successes_116": int(sum(1 for row in rows if row.get("engine") == "116(115)" and row.get("success") is True)),
        "successes_119": int(sum(1 for row in rows if row.get("engine") == "119(118)" and row.get("success") is True)),
        "roots_116": int(sum(int(row.get("unique_roots") or 0) for row in rows if row.get("engine") == "116(115)")),
        "roots_119": int(sum(int(row.get("unique_roots") or 0) for row in rows if row.get("engine") == "119(118)")),
        "requested_116": int(sum(int(row.get("requested_roots") or 0) for row in rows if row.get("engine") == "116(115)")),
        "requested_119": int(sum(int(row.get("requested_roots") or 0) for row in rows if row.get("engine") == "119(118)")),
        "seconds_116": float(sum(float(row.get("seconds") or 0.0) for row in rows if row.get("engine") == "116(115)")),
        "seconds_119": float(sum(float(row.get("seconds") or 0.0) for row in rows if row.get("engine") == "119(118)")),
        "root_delta_119_minus_116": int(sum(int(pair.get("root_delta_119_minus_116") or 0) for pair in pairs)),
        "wins_116": int(sum(1 for pair in pairs if pair.get("winner") == "116(115)")),
        "wins_119": int(sum(1 for pair in pairs if pair.get("winner") == "119(118)")),
        "ties": int(sum(1 for pair in pairs if pair.get("winner") == "tie")),
        "seconds": float(seconds_all),
    }

    payload = {
        "script": "120_pandrosion_probe_vs_vectorized_benchmark.py",
        "engines": {
            "116(115)": "flow/116_pandrosion_multifamily_vectorized_pure_pandrosion.py over flow/115_pandrosion_vectorized_pure_pandrosion.py",
            "119(118)": "flow/119_pandrosion_multifamily_probe_aware_pure_thales_engine.py over flow/118_pandrosion_probe_aware_pure_thales_engine.py",
        },
        "cases": cases,
        "families": families,
        "parameters": {
            "count": int(args.count),
            "pool": int(args.pool),
            "epochs": int(args.epochs),
            "accept": float(args.accept),
            "tol": float(args.tol),
            "cluster_sep": float(args.cluster_sep),
            "line_search": int(args.line_search),
            "probe_scale": float(args.probe_scale),
            "probe_candidates_119": int(args.probe_candidates),
            "probe_radii_119": str(args.probe_radii),
            "probe_self_119": bool(args.probe_self),
            "same_systems": bool(args.same_systems),
            "seed_index": int(args.seed_index),
            "sparse_terms": int(args.sparse_terms),
            "sparse_frac": float(args.sparse_frac),
        },
        "summary": summary,
        "pairs": pairs,
        "runs": rows,
        "full_results": full_results if bool(args.verbose_trials) else {},
    }

    out_json = Path(args.out) if args.out else Path(args.outdir) / "120_probe_vs_vectorized_benchmark.json"
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
        f"summary: 116 roots={summary['roots_116']}/{summary['requested_116']} "
        f"successes={summary['successes_116']}/{summary['pairs']} seconds={summary['seconds_116']:.4f}",
        flush=True,
    )
    print(
        f"summary: 119 roots={summary['roots_119']}/{summary['requested_119']} "
        f"successes={summary['successes_119']}/{summary['pairs']} seconds={summary['seconds_119']:.4f}",
        flush=True,
    )
    print(f"root_delta_119_minus_116={summary['root_delta_119_minus_116']}", flush=True)
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
