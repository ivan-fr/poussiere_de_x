#!/usr/bin/env python3
"""Large diagnostic benchmark for 304 Universal Atlas Pandrosion.

No Docker. Runs the standalone 304 engine over a Smale-17-style Kostlan suite,
exports per-run JSON, summary.csv, and report.md.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "flow" / "304_pandrosion_universal_atlas_full_cubic_halley_numpy_engine.py"
OUTDIR = ROOT / "benchmarks" / "304_smale17_large"
RUNS_DIR = OUTDIR / "runs"
SUMMARY_CSV = OUTDIR / "summary.csv"
REPORT = OUTDIR / "report.md"

CASES = (
    [(2, d) for d in range(2, 16)]
    + [(3, d) for d in range(2, 9)]
    + [(4, d) for d in range(2, 6)]
    + [(5, d) for d in range(2, 5)]
    + [(6, d) for d in range(2, 4)]
    + [(7, 2), (8, 2), (10, 2)]
)
SEEDS = list(range(16))
BASE_ARGS = [
    "--pool", "8192",
    "--epochs", "40",
    "--accept", "1e-8",
    "--tol", "1e-12",
    "--line-search", "14",
    "--probe-candidates", "12",
    "--trial-timeout", "0",
    "--keep-trials", "40",
    "--universal-cells", "16",
    "--universal-shells", "5",
]


def load_engine():
    spec = importlib.util.spec_from_file_location("pandrosion304", ENGINE)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def z_from_json(z):
    import numpy as np
    return np.asarray([complex(a, b) for a, b in z], dtype=np.complex128)


def jacobian(system, z):
    import numpy as np
    z = np.asarray(z, dtype=np.complex128)
    n = system.n
    J = np.empty((n, n), dtype=np.complex128)
    pows = system._powers(z)
    for j in range(n):
        term = np.ones(system.terms_per_poly, dtype=np.complex128)
        for k in range(n):
            expk = system.exps[:, k]
            if k == j:
                vals = np.zeros(system.terms_per_poly, dtype=np.complex128)
                mask = expk > 0
                vals[mask] = expk[mask] * pows[k][expk[mask] - 1]
                term *= vals
            else:
                term *= pows[k][expk]
        J[:, j] = system.coeff @ term
    return J


def certify(mod, n, d, seed_index, roots):
    import numpy as np
    system = mod.DenseKostlanSystem.make(n, d, seed_index=seed_index, equation_normalize=False)
    cert = 0
    strong = 0
    residuals = []
    polished = []
    conds = []
    for r in roots:
        z = z_from_json(r["z"])
        F = system.eval(z)
        res = float(np.linalg.norm(F))
        residuals.append(res)
        try:
            J = jacobian(system, z)
            cond = float(np.linalg.cond(J))
            conds.append(cond)
            step = np.linalg.solve(J, -F)
            beta = float(np.linalg.norm(step))
            zp = z + step
            pres = float(np.linalg.norm(system.eval(zp)))
            polished.append(pres)
            # Conservative alpha-lite surrogate, intentionally empirical.
            if math.isfinite(cond) and res < 1e-8 and beta < 1e-5 and pres < 1e-10:
                cert += 1
            if math.isfinite(cond) and res < 1e-10 and beta < 1e-7 and pres < 1e-12:
                strong += 1
        except Exception:
            polished.append(float("inf"))
    return {
        "certified": cert,
        "strong": strong,
        "median_residual": median(residuals) if residuals else float("inf"),
        "median_polished_residual": median(polished) if polished else float("inf"),
        "median_cond": median(conds) if conds else float("inf"),
    }


def q50(xs):
    ys = sorted(float(x) for x in xs if x is not None and math.isfinite(float(x)))
    return ys[len(ys)//2] if ys else float("nan")


def q90(xs):
    ys = sorted(float(x) for x in xs if x is not None and math.isfinite(float(x)))
    if not ys:
        return float("nan")
    return ys[min(len(ys)-1, int(math.ceil(0.9 * len(ys))) - 1)]


def main() -> int:
    mod = load_engine()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    t0 = time.time()
    total = len(CASES) * len(SEEDS)
    run_no = 0
    for n, d in CASES:
        count = min(d ** n, 50)
        for seed in SEEDS:
            run_no += 1
            out = RUNS_DIR / f"304_ks{n}x{d}_seed{seed}.json"
            cmd = [sys.executable, str(ENGINE), "--cases", f"{n},{d}", "--seed-index", str(seed), "--count", str(count), "--out", str(out), *BASE_ARGS]
            print(f"RUN {run_no}/{total} 304 ks({n},{d}) seed={seed} count={count}", flush=True)
            t = time.time()
            proc = subprocess.run(cmd, cwd=str(ROOT), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            seconds = time.time() - t
            roots = 0
            success = False
            trials = None
            failures = None
            duplicates = None
            sec_per_root = float("inf")
            c = {"certified": 0, "strong": 0, "median_residual": float("inf"), "median_polished_residual": float("inf"), "median_cond": float("inf")}
            if out.exists():
                data = json.loads(out.read_text())
                s = data.get("summary", {})
                roots = int(s.get("unique_roots", 0))
                success = bool(s.get("success", False))
                trials = int(s.get("trials_used", 0))
                failures = int(s.get("failures", 0))
                duplicates = int(s.get("duplicates", 0))
                sec_per_root = float(s.get("extract_seconds", seconds)) / max(1, roots)
                c = certify(mod, n, d, seed, data.get("roots", []))
            row = {
                "engine": 304, "n": n, "d": d, "seed": seed, "count": count,
                "rc": proc.returncode, "success": success, "roots": roots,
                "certified": c["certified"], "strong": c["strong"],
                "trials": trials, "failures": failures, "duplicates": duplicates,
                "seconds": seconds, "sec_per_root": sec_per_root,
                "median_residual": c["median_residual"],
                "median_polished_residual": c["median_polished_residual"],
                "median_cond": c["median_cond"],
                "out": str(out.relative_to(ROOT)),
            }
            rows.append(row)
            print(f"  rc={proc.returncode} roots={roots}/{count} cert={c['certified']} strong={c['strong']} sec/root={sec_per_root} medres={c['median_residual']}", flush=True)
            if proc.returncode != 0:
                (RUNS_DIR / f"304_ks{n}x{d}_seed{seed}.log").write_text(proc.stdout[-20000:])
            write_outputs(rows, total, time.time() - t0)
    write_outputs(rows, total, time.time() - t0)
    print(f"WROTE {SUMMARY_CSV}")
    print(f"WROTE {REPORT}")
    return 0


def write_outputs(rows, total_runs, elapsed):
    fields = ["engine","n","d","seed","count","rc","success","roots","certified","strong","trials","failures","duplicates","seconds","sec_per_root","median_residual","median_polished_residual","median_cond","out"]
    with SUMMARY_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    completed = len(rows)
    ok_runs = sum(1 for r in rows if int(r["rc"]) == 0)
    roots = sum(int(r["roots"]) for r in rows)
    requested = sum(int(r["count"]) for r in rows)
    cert = sum(int(r["certified"]) for r in rows)
    strong = sum(int(r["strong"]) for r in rows)
    report = []
    report.append("# Large 304 universal-atlas Smale-17 diagnostic benchmark\n")
    report.append("Empirical stress benchmark for the fixed universal atlas idea. This is not a proof and not a homotopy baseline.\n")
    report.append(f"Cases: {CASES}\n")
    report.append(f"Seeds: {SEEDS}\n")
    report.append(f"Completed runs: {completed}/{total_runs}\n")
    report.append(f"Elapsed wall time: {elapsed:.1f}s\n")
    report.append(f"Base args: `{' '.join(BASE_ARGS)}`\n")
    report.append("Requested roots: `min(Bezout, 50)`\n")
    report.append("\n## Global aggregate\n")
    report.append(f"- Runs OK: {ok_runs}/{completed}\n")
    report.append(f"- Roots: {roots}/{requested}\n")
    report.append(f"- Certified alpha-lite roots: {cert}/{roots}\n")
    report.append(f"- Strong roots: {strong}/{roots}\n")
    report.append(f"- Median seconds/root: {q50([r['sec_per_root'] for r in rows])}\n")
    report.append(f"- P90 seconds/root: {q90([r['sec_per_root'] for r in rows])}\n")
    report.append(f"- Median trials/run: {q50([r['trials'] for r in rows if r['trials'] is not None])}\n")
    report.append(f"- P90 trials/run: {q90([r['trials'] for r in rows if r['trials'] is not None])}\n")
    report.append(f"- Median residual: {q50([r['median_residual'] for r in rows])}\n")
    report.append(f"- P90 run-median residual: {q90([r['median_residual'] for r in rows])}\n")
    report.append(f"- Median polished residual: {q50([r['median_polished_residual'] for r in rows])}\n")
    report.append(f"- Median cond(J): {q50([r['median_cond'] for r in rows])}\n")
    report.append("\n## Per-case aggregate\n\n")
    report.append("| case | runs | roots | certified | strong | median sec/root | p90 sec/root | median trials | median residual | max run-median residual |\n")
    report.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    for n,d in CASES:
        rs = [r for r in rows if r['n'] == n and r['d'] == d]
        if not rs:
            continue
        report.append(f"| ks({n},{d}) | {sum(1 for r in rs if int(r['rc'])==0)}/{len(rs)} | {sum(int(r['roots']) for r in rs)}/{sum(int(r['count']) for r in rs)} | {sum(int(r['certified']) for r in rs)}/{sum(int(r['roots']) for r in rs)} | {sum(int(r['strong']) for r in rs)}/{sum(int(r['roots']) for r in rs)} | {q50([r['sec_per_root'] for r in rs])} | {q90([r['sec_per_root'] for r in rs])} | {q50([r['trials'] for r in rs if r['trials'] is not None])} | {q50([r['median_residual'] for r in rs])} | {max(float(r['median_residual']) for r in rs if math.isfinite(float(r['median_residual'])))} |\n")
    REPORT.write_text("".join(report))

if __name__ == "__main__":
    raise SystemExit(main())
