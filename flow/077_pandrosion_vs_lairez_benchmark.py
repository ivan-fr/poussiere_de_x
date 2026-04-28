"""
FLOW 077 -- broad benchmark for 076 homothetic Pandrosion vs Lairez-style.

This benchmark is deliberately an orchestration layer.  It treats
flow/076_pandrosion_homothetic_geometry.py as the candidate algorithm and runs
it in isolated subprocesses, then compares against the local Lairez-style
baseline from flow/069 through the dependency module imported by 076.

The default run is modest.  Use --suite high or --suite full with a larger
--max-bezout for high-degree KS stress tests.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import statistics
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FLOW076_PATH = HERE / "076_pandrosion_homothetic_geometry.py"


def load_076():
    spec = importlib.util.spec_from_file_location("flow076_for_077", str(FLOW076_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import 076 from {FLOW076_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow076_for_077"] = mod
    spec.loader.exec_module(mod)
    return mod


f076 = load_076()
m = f076.m
m069 = m.m069


@dataclass(frozen=True)
class Config076:
    name: str
    homothety: str
    equation_normalize: bool
    window_fallback: bool = True
    dt0: float = 0.0
    dtmax: float = 0.0


@dataclass
class BenchRow:
    family: str
    n: int
    d: int
    seed_index: int
    seed: int
    terms: int
    bezout: int
    alg: str
    config: str
    roots: int
    coverage: float
    paths: int
    candidates: int
    batches: int
    retries_used: int
    steps_or_work: int
    max_residual: float
    seconds_observed: float
    wall_seconds: float
    status: str
    notes: str
    artifact_csv: str
    artifact_json: str
    artifact_md: str


def parse_case(s: str) -> Tuple[int, int]:
    if "," in s:
        a, b = s.split(",", 1)
    elif "x" in s:
        a, b = s.split("x", 1)
    else:
        raise argparse.ArgumentTypeError("case must look like 2,8 or 2x8")
    return int(a), int(b)


def encode_case(n: int, d: int) -> str:
    return f"{n},{d}"


def case_slug(family: str, n: int, d: int, seed_index: int, config: str) -> str:
    return f"{family}_{n}x{d}_seed{seed_index}_{config}".replace("-", "_")


def suite_cases(name: str) -> List[Tuple[str, Tuple[int, int]]]:
    name = name.lower()
    if name == "smoke":
        return [("ks", (2, 3)), ("ks", (2, 5))]
    if name == "quick":
        return [
            ("dense064", (2, 3)),
            ("ks", (2, 4)),
            ("ks", (2, 8)),
            ("ks", (3, 3)),
            ("sparse_ks", (2, 10)),
        ]
    if name in {"high", "high_ks"}:
        return [
            ("ks", (2, 8)),
            ("ks", (2, 10)),
            ("ks", (2, 12)),
            ("ks", (2, 14)),
            ("ks", (3, 4)),
            ("ks", (3, 5)),
            ("ks", (4, 3)),
            ("sparse_ks", (2, 16)),
            ("sparse_ks", (3, 5)),
        ]
    if name in {"full", "complete"}:
        return [
            ("dense064", (2, 2)),
            ("dense064", (2, 3)),
            ("dense064", (2, 4)),
            ("dense064", (3, 2)),
            ("dense064", (3, 3)),
            ("dense064", (4, 2)),
            ("ks", (2, 4)),
            ("ks", (2, 6)),
            ("ks", (2, 8)),
            ("ks", (2, 10)),
            ("ks", (2, 12)),
            ("ks", (3, 3)),
            ("ks", (3, 4)),
            ("ks", (3, 5)),
            ("ks", (4, 2)),
            ("ks", (4, 3)),
            ("sparse_ks", (2, 8)),
            ("sparse_ks", (2, 12)),
            ("sparse_ks", (2, 16)),
            ("sparse_ks", (3, 4)),
            ("sparse_ks", (3, 5)),
        ]
    raise ValueError(f"unknown suite: {name}")


def config_set(name: str) -> List[Config076]:
    name = name.lower()
    if name in {"production", "default"}:
        return [Config076("system_eqnorm", "system", True)]
    if name in {"ablation", "ablations"}:
        return [
            Config076("system_eqnorm", "system", True),
            Config076("hybrid_eqnorm", "hybrid", True),
            Config076("system_noeq", "system", False),
            Config076("none_noeq", "none", False),
        ]
    if name in {"full", "complete"}:
        return [
            Config076("system_eqnorm", "system", True),
            Config076("hybrid_eqnorm", "hybrid", True),
            Config076("roots_eqnorm", "roots", True),
            Config076("system_noeq", "system", False),
            Config076("hybrid_noeq", "hybrid", False),
            Config076("none_noeq", "none", False),
        ]
    raise ValueError(f"unknown config set: {name}")


def read_summary_csv(path: Path) -> Dict[str, str]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"empty 076 summary: {path}")
    # 076 may optionally include a Lairez reference.  Keep the candidate row.
    for row in rows:
        if row.get("alg") == "076-homothetic-system-geometry":
            return row
    return rows[0]


def safe_float(value: object, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def safe_int(value: object, default: int = 0) -> int:
    try:
        return int(float(value))
    except Exception:
        return default


def text_payload(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def row_from_076_summary(
    family: str,
    n: int,
    d: int,
    seed_index: int,
    seed: int,
    terms: int,
    bezout: int,
    cfg: Config076,
    summary: Dict[str, str],
    wall_seconds: float,
    artifact_csv: Path,
    artifact_json: Path,
    artifact_md: Path,
) -> BenchRow:
    return BenchRow(
        family=family,
        n=n,
        d=d,
        seed_index=seed_index,
        seed=seed,
        terms=terms,
        bezout=bezout,
        alg="076-homothetic",
        config=cfg.name,
        roots=safe_int(summary.get("roots")),
        coverage=safe_float(summary.get("coverage")),
        paths=safe_int(summary.get("path_rows")),
        candidates=safe_int(summary.get("candidates")),
        batches=safe_int(summary.get("batches")),
        retries_used=safe_int(summary.get("retries_used")),
        steps_or_work=0,
        max_residual=safe_float(summary.get("max_residual"), float("inf")),
        seconds_observed=safe_float(summary.get("seconds_observed"), wall_seconds),
        wall_seconds=wall_seconds,
        status=str(summary.get("status", "unknown")),
        notes=str(summary.get("notes", "")),
        artifact_csv=str(artifact_csv),
        artifact_json=str(artifact_json),
        artifact_md=str(artifact_md),
    )


def timeout_row(
    family: str,
    n: int,
    d: int,
    seed_index: int,
    seed: int,
    terms: int,
    bezout: int,
    cfg: Config076,
    status: str,
    wall_seconds: float,
    artifact_csv: Path,
    artifact_json: Path,
    artifact_md: Path,
    notes: str,
) -> BenchRow:
    return BenchRow(
        family=family,
        n=n,
        d=d,
        seed_index=seed_index,
        seed=seed,
        terms=terms,
        bezout=bezout,
        alg="076-homothetic",
        config=cfg.name,
        roots=0,
        coverage=0.0,
        paths=0,
        candidates=0,
        batches=0,
        retries_used=0,
        steps_or_work=0,
        max_residual=float("inf"),
        seconds_observed=wall_seconds,
        wall_seconds=wall_seconds,
        status=status,
        notes=notes,
        artifact_csv=str(artifact_csv),
        artifact_json=str(artifact_json),
        artifact_md=str(artifact_md),
    )


def run_076_case(
    args: argparse.Namespace,
    family: str,
    n: int,
    d: int,
    seed_index: int,
    seed: int,
    terms: int,
    bezout: int,
    cfg: Config076,
    outdir: Path,
) -> BenchRow:
    slug = case_slug(family, n, d, seed_index, cfg.name)
    case_dir = outdir / "076_batches" / slug
    case_dir.mkdir(parents=True, exist_ok=True)
    artifact_csv = outdir / f"{slug}.csv"
    artifact_batch_csv = outdir / f"{slug}_batches.csv"
    artifact_json = outdir / f"{slug}_roots.json"
    artifact_md = outdir / f"{slug}.md"
    log_path = outdir / f"{slug}.log"
    cmd = [
        sys.executable,
        "-S",
        str(FLOW076_PATH),
        "--family",
        family,
        "--case",
        encode_case(n, d),
        "--seed-index",
        str(seed_index),
        "--retries",
        str(args.retries),
        "--base-chunk-size",
        str(args.base_chunk_size),
        "--parallel-base",
        str(args.parallel_base),
        "--micro-batch",
        str(args.micro_batch),
        "--micro-limit",
        str(args.micro_limit),
        "--cluster-sep",
        str(args.cluster_sep),
        "--tol",
        str(args.tol),
        "--max-steps",
        str(args.max_steps),
        "--max-epochs",
        str(args.max_epochs),
        "--quad-cap",
        str(args.quad_cap),
        "--batch-timeout",
        str(args.batch_timeout),
        "--time-budget",
        str(args.timeout_076),
        "--homothety",
        cfg.homothety,
        "--scale-min",
        str(args.scale_min),
        "--scale-max",
        str(args.scale_max),
        "--system-scale-strength",
        str(args.system_scale_strength),
        "--root-scale-strength",
        str(args.root_scale_strength),
        "--root-scale-quantile",
        str(args.root_scale_quantile),
        "--root-scale-trigger",
        str(args.root_scale_trigger),
        "--modes",
        args.modes,
        "--rescue-modes",
        args.rescue_modes,
        "--outdir",
        str(case_dir),
        "--csv",
        str(artifact_csv),
        "--batch-csv",
        str(artifact_batch_csv),
        "--md",
        str(artifact_md),
        "--roots-json",
        str(artifact_json),
    ]
    if args.stop_at_bezout:
        cmd.append("--stop-at-bezout")
    if cfg.equation_normalize:
        cmd.append("--equation-normalize")
    if cfg.dt0 > 0:
        cmd.extend(["--dt0", str(cfg.dt0)])
    if cfg.dtmax > 0:
        cmd.extend(["--dtmax", str(cfg.dtmax)])
    if not cfg.window_fallback:
        cmd.append("--no-window-fallback")

    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(ROOT),
            timeout=(args.timeout_076 + args.timeout_grace) if args.timeout_076 > 0 else None,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        wall = time.time() - t0
        log_path.write_text(text_payload(proc.stdout))
        if proc.returncode != 0:
            return timeout_row(
                family,
                n,
                d,
                seed_index,
                seed,
                terms,
                bezout,
                cfg,
                f"error:{proc.returncode}",
                wall,
                artifact_csv,
                artifact_json,
                artifact_md,
                f"076 subprocess failed; see {log_path}",
            )
        summary = read_summary_csv(artifact_csv)
        return row_from_076_summary(
            family,
            n,
            d,
            seed_index,
            seed,
            terms,
            bezout,
            cfg,
            summary,
            wall,
            artifact_csv,
            artifact_json,
            artifact_md,
        )
    except subprocess.TimeoutExpired as exc:
        wall = time.time() - t0
        log_path.write_text(text_payload(exc.stdout))
        return timeout_row(
            family,
            n,
            d,
            seed_index,
            seed,
            terms,
            bezout,
            cfg,
            "timeout",
            wall,
            artifact_csv,
            artifact_json,
            artifact_md,
            f"076 subprocess timed out after {args.timeout_076:.2f}s; see {log_path}",
        )


def lairez_row_from_result(
    family: str,
    n: int,
    d: int,
    seed_index: int,
    seed: int,
    terms: int,
    bezout: int,
    result: object,
    seconds: float,
    status: str,
    retries: int,
) -> BenchRow:
    roots = list(getattr(result, "roots", [])) if result is not None else []
    coverage = float(getattr(result, "coverage", len(roots) / max(1, bezout))) if result is not None else 0.0
    return BenchRow(
        family=family,
        n=n,
        d=d,
        seed_index=seed_index,
        seed=seed,
        terms=terms,
        bezout=bezout,
        alg="lairez-style",
        config=f"retries{retries}",
        roots=len(roots),
        coverage=coverage,
        paths=int(getattr(result, "paths", 0)) if result is not None else 0,
        candidates=len(roots),
        batches=0,
        retries_used=retries,
        steps_or_work=int(getattr(result, "steps", 0)) + int(getattr(result, "newton_iters", 0)) if result is not None else 0,
        max_residual=float(getattr(result, "max_residual", float("inf"))) if result is not None else float("inf"),
        seconds_observed=seconds,
        wall_seconds=seconds,
        status=status if status != "ok" else ("ok" if coverage >= 0.999999 else "partial"),
        notes="local gamma total-degree homotopy + analytic Newton corrector",
        artifact_csv="",
        artifact_json="",
        artifact_md="",
    )


def run_lairez_case(
    args: argparse.Namespace,
    family: str,
    n: int,
    d: int,
    seed_index: int,
    seed: int,
    terms: int,
    bezout: int,
    target,
) -> BenchRow:
    kwargs = {
        "seed": 91000 + seed,
        "tol": args.tol,
        "max_steps": args.lairez_max_steps,
        "max_newton_iter": args.lairez_newton_iters,
        "retries": args.lairez_retries,
    }
    timeout = args.timeout_lairez
    try:
        if hasattr(m069, "run_method_with_timeout"):
            result, seconds, status = m069.run_method_with_timeout("lairez", target, kwargs, timeout)
            return lairez_row_from_result(family, n, d, seed_index, seed, terms, bezout, result, seconds, status, args.lairez_retries)
        t0 = time.time()
        result, seconds = m069.run_lairez(target, **kwargs)
        return lairez_row_from_result(family, n, d, seed_index, seed, terms, bezout, result, seconds or (time.time() - t0), "ok", args.lairez_retries)
    except Exception as exc:
        return lairez_row_from_result(family, n, d, seed_index, seed, terms, bezout, None, 0.0, f"error:{exc}", args.lairez_retries)


def print_row(row: BenchRow) -> None:
    print(
        f"{row.family:>9} ({row.n:>2},{row.d:<2}) seed{row.seed_index:<2} "
        f"{row.alg:>15} {row.config:>14} | Bez={row.bezout:<4} terms={row.terms:<4} "
        f"roots={row.roots:<4} cov={100*row.coverage:6.1f}% "
        f"paths={row.paths:<4} cand={row.candidates:<4} batches={row.batches:<3} "
        f"res={row.max_residual:.1e} time={row.wall_seconds:7.2f}s {row.status}",
        flush=True,
    )


def summarize(rows: Sequence[BenchRow]) -> str:
    lines: List[str] = []
    lines.append("=" * 140)
    lines.append("077 SUMMARY BY ALGORITHM / CONFIG")
    lines.append("=" * 140)
    groups: Dict[Tuple[str, str], List[BenchRow]] = {}
    for r in rows:
        groups.setdefault((r.alg, r.config), []).append(r)
    lines.append(f"{'alg':>15} {'config':>14} {'rows':>5} {'avg cov':>9} {'full%':>8} {'avg time':>10} {'avg paths/Bez':>14}")
    lines.append("-" * 140)
    for (alg, cfg), rs in sorted(groups.items()):
        avg_cov = statistics.mean(r.coverage for r in rs)
        full = sum(1 for r in rs if r.coverage >= 0.999999) / len(rs)
        avg_time = statistics.mean(r.wall_seconds for r in rs)
        avg_paths = statistics.mean(r.paths / max(1, r.bezout) for r in rs)
        lines.append(f"{alg:>15} {cfg:>14} {len(rs):>5} {100*avg_cov:8.2f}% {100*full:7.1f}% {avg_time:9.2f}s {avg_paths:13.2f}")

    lines.append("")
    lines.append("=" * 140)
    lines.append("HEAD TO HEAD AGAINST LAIREZ")
    lines.append("=" * 140)
    by_key: Dict[Tuple[str, int, int, int], Dict[str, List[BenchRow]]] = {}
    for r in rows:
        by_key.setdefault((r.family, r.n, r.d, r.seed_index), {}).setdefault(r.alg, []).append(r)

    compared = 0
    cov_wins = 0
    cov_ties = 0
    time_wins = 0
    for key, dct in sorted(by_key.items()):
        lairez_rows = dct.get("lairez-style", [])
        cand_rows = dct.get("076-homothetic", [])
        if not lairez_rows or not cand_rows:
            continue
        lrow = lairez_rows[0]
        for crow in cand_rows:
            compared += 1
            dcov = crow.coverage - lrow.coverage
            ratio = crow.wall_seconds / lrow.wall_seconds if lrow.wall_seconds > 0 else float("inf")
            pratio = crow.paths / max(1, lrow.paths)
            if dcov > 1e-12:
                cov_wins += 1
            elif abs(dcov) <= 1e-12:
                cov_ties += 1
            if crow.wall_seconds < lrow.wall_seconds:
                time_wins += 1
            family, n, d, seed_index = key
            lines.append(
                f"{family:>9} ({n},{d}) seed{seed_index:<2} {crow.config:>14}: "
                f"dcov={100*dcov:+7.1f}% time_ratio={ratio:7.2f} "
                f"paths_ratio={pratio:6.2f} status={crow.status}/{lrow.status}"
            )
    if compared:
        lines.append("-" * 140)
        lines.append(f"coverage: 076 wins {cov_wins}/{compared}, ties {cov_ties}/{compared}; time: 076 faster {time_wins}/{compared}")
    return "\n".join(lines)


def write_csv(rows: Sequence[BenchRow], path: Path) -> None:
    if not rows:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(rows: Sequence[BenchRow], path: Path) -> None:
    path.write_text(json.dumps([asdict(row) for row in rows], indent=2))


def write_md(rows: Sequence[BenchRow], summary: str, path: Path) -> None:
    with path.open("w") as f:
        f.write("# Flow 077 benchmark\n\n")
        f.write("076 homothetic system-generated Pandrosion geometry is compared against the local Lairez-style gamma total-degree continuation baseline.\n\n")
        f.write("| family | n | d | seed | Bezout | terms | algorithm | config | roots | coverage | paths | candidates | batches | max residual | seconds | status |\n")
        f.write("|---|---:|---:|---:|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|\n")
        for r in rows:
            f.write(
                f"| {r.family} | {r.n} | {r.d} | {r.seed_index} | {r.bezout} | {r.terms} | "
                f"{r.alg} | {r.config} | {r.roots} | {100*r.coverage:.1f}% | {r.paths} | "
                f"{r.candidates} | {r.batches} | {r.max_residual:.2e} | {r.wall_seconds:.2f}s | {r.status} |\n"
            )
        f.write("\n```text\n")
        f.write(summary)
        f.write("\n```\n")


def selected_tasks(args: argparse.Namespace) -> List[Tuple[str, Tuple[int, int]]]:
    if args.family:
        if not args.cases:
            raise SystemExit("--family requires --cases")
        return [(args.family, case) for case in args.cases]
    return suite_cases(args.suite)


def main() -> None:
    parser = argparse.ArgumentParser(description="077 broad benchmark: 076 homothetic Pandrosion vs Lairez-style")
    parser.add_argument("--suite", choices=["smoke", "quick", "high", "high_ks", "full", "complete"], default="quick")
    parser.add_argument("--family", default=None, help="Override suite with one family: dense064, ks, sparse_ks")
    parser.add_argument("--cases", nargs="*", type=parse_case, default=None, help="Cases like 2,8 3,4")
    parser.add_argument("--seeds", type=int, default=1)
    parser.add_argument("--only", choices=["both", "076", "lairez"], default="both")
    parser.add_argument("--config-set", choices=["production", "default", "ablation", "ablations", "full", "complete"], default="production")
    parser.add_argument("--max-bezout", type=int, default=256, help="Skip cases with Bezout above this; 0 disables")
    parser.add_argument("--dry-run", action="store_true")

    parser.add_argument("--retries", type=int, default=1)
    parser.add_argument("--base-chunk-size", type=int, default=0)
    parser.add_argument("--parallel-base", type=int, default=0)
    parser.add_argument("--micro-batch", type=int, default=2)
    parser.add_argument("--micro-limit", type=int, default=0)
    parser.add_argument("--batch-timeout", type=float, default=20.0)
    parser.add_argument("--timeout-076", type=float, default=20.0, help="Wall timeout per 076 case/config; 0 disables")
    parser.add_argument("--timeout-grace", type=float, default=5.0, help="extra seconds for 076 to write partial outputs after --timeout-076")
    parser.add_argument("--stop-at-bezout", action="store_true", default=True)
    parser.add_argument("--no-stop-at-bezout", dest="stop_at_bezout", action="store_false")
    parser.add_argument("--cluster-sep", type=float, default=1e-6)
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--max-steps", type=int, default=120)
    parser.add_argument("--max-epochs", type=int, default=4)
    parser.add_argument("--quad-cap", type=int, default=12)
    parser.add_argument("--modes", default="system,integral_gl,blend")
    parser.add_argument("--rescue-modes", default="")
    parser.add_argument("--scale-min", type=float, default=0.25)
    parser.add_argument("--scale-max", type=float, default=4.0)
    parser.add_argument("--system-scale-strength", type=float, default=1.0)
    parser.add_argument("--root-scale-strength", type=float, default=0.60)
    parser.add_argument("--root-scale-quantile", type=float, default=0.75)
    parser.add_argument("--root-scale-trigger", type=float, default=1.75)

    parser.add_argument("--lairez-max-steps", type=int, default=420)
    parser.add_argument("--lairez-newton-iters", type=int, default=12)
    parser.add_argument("--lairez-retries", type=int, default=2)
    parser.add_argument("--timeout-lairez", type=float, default=20.0, help="Wall timeout per Lairez case; 0 disables")

    parser.add_argument("--outdir", default="077_benchmark_artifacts")
    parser.add_argument("--csv", default="077_benchmark.csv")
    parser.add_argument("--json", default="077_benchmark.json")
    parser.add_argument("--md", default="077_benchmark.md")
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    tasks = selected_tasks(args)
    configs = config_set(args.config_set)

    print("=" * 140)
    print("077 -- 076 homothetic Pandrosion vs Lairez-style")
    print("=" * 140)
    print(
        f"suite={args.suite}, seeds={args.seeds}, config_set={args.config_set}, "
        f"max_bezout={args.max_bezout or 'unbounded'}, only={args.only}"
    )
    print(
        f"076: retries={args.retries}, base_chunk={args.base_chunk_size}, parallel_base={args.parallel_base}, "
        f"modes={args.modes}, rescue={args.rescue_modes or '-'}"
    )
    print(
        f"Lairez: retries={args.lairez_retries}, max_steps={args.lairez_max_steps}, "
        f"newton_iters={args.lairez_newton_iters}"
    )
    print("=" * 140, flush=True)

    rows: List[BenchRow] = []
    for family, (n, d) in tasks:
        bezout_est = d ** n
        if args.max_bezout and bezout_est > args.max_bezout:
            print(f"SKIP {family:>9} ({n},{d}) Bezout={bezout_est} > --max-bezout={args.max_bezout}", flush=True)
            continue
        for seed_index in range(max(1, args.seeds)):
            seed = m.seed_for(family, n, d, seed_index)
            target = m.gen_system(family, n, d, seed)
            terms = m.term_count(target)
            bezout = m.bezout(target)
            print(f"CASE {family:>9} ({n},{d}) seed{seed_index} seed={seed} terms={terms} Bezout={bezout}", flush=True)
            if args.dry_run:
                continue
            if args.only in {"both", "076"}:
                for cfg in configs:
                    row = run_076_case(args, family, n, d, seed_index, seed, terms, bezout, cfg, outdir)
                    rows.append(row)
                    print_row(row)
            if args.only in {"both", "lairez"}:
                row = run_lairez_case(args, family, n, d, seed_index, seed, terms, bezout, target)
                rows.append(row)
                print_row(row)

    if rows:
        summary = summarize(rows)
        print(summary, flush=True)
        if args.csv:
            write_csv(rows, Path(args.csv))
            print(f"CSV written to {args.csv}", flush=True)
        if args.json:
            write_json(rows, Path(args.json))
            print(f"JSON written to {args.json}", flush=True)
        if args.md:
            write_md(rows, summary, Path(args.md))
            print(f"Markdown written to {args.md}", flush=True)
    print("=" * 140, flush=True)


if __name__ == "__main__":
    main()
