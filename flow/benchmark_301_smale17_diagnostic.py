#!/usr/bin/env python3
"""Smale-17 diagnostic benchmark for 301.

This is not a proof and not a homotopy comparison. It measures whether the
standalone full-cubic finite-slope Halley-Pandrosion engine behaves like a
plausible research candidate on random dense Kostlan systems as n,d scale.

Stdlib-only harness; launches flow/301 and summarizes coverage, time/root,
trials/root, failures, residual quality, and eval counts.
"""
from __future__ import annotations

import csv, json, math, statistics, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FLOW = ROOT / "flow"
SCRIPT = FLOW / "301_pandrosion_anchored_full_cubic_halley_numpy_engine.py"
OUT = ROOT / "benchmarks" / "301_smale17_diagnostic"
OUT.mkdir(parents=True, exist_ok=True)

# Scaling suite: low-dimensional degree scaling + dimension scaling.
CASES = [
    (2,2), (2,3), (2,4), (2,5), (2,6), (2,8), (2,10),
    (3,2), (3,3), (3,4), (3,5),
    (4,2), (4,3),
    (5,2), (6,2),
]
SEEDS = list(range(8))
BASE = [
    "--pool", "4096",
    "--epochs", "32",
    "--accept", "1e-8",
    "--tol", "1e-12",
    "--line-search", "12",
    "--probe-candidates", "10",
    "--trial-timeout", "0",
    "--keep-trials", "160",
]


def count_for(n: int, d: int) -> int:
    # Full Bezout for small systems; capped sample for larger systems.
    return min(int(d ** n), 20)


def run_one(n: int, d: int, seed: int) -> dict:
    case = f"{n},{d}"
    cnt = count_for(n, d)
    out = OUT / f"301_ks{n}x{d}_seed{seed}.json"
    cmd = [sys.executable, str(SCRIPT), "--cases", case, "--seed-index", str(seed), "--count", str(cnt), "--out", str(out), *BASE]
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    wall = time.time() - t0
    rec = {
        "engine": "301",
        "case": case,
        "n": n,
        "degree": d,
        "seed_index": seed,
        "requested_count": cnt,
        "returncode": proc.returncode,
        "wall_seconds": wall,
        "json": str(out),
        "stderr_tail": proc.stderr[-1200:],
        "stdout_tail": proc.stdout[-1200:],
    }
    if proc.returncode == 0 and out.exists():
        data = json.loads(out.read_text())
        s = data.get("summary", {})
        roots = data.get("roots", [])
        residuals = [float(r.get("residual", float("inf"))) for r in roots]
        eval_stats = s.get("eval_stats", {}) or {}
        unique = int(s.get("unique_roots") or 0)
        trials = int(s.get("trials_used") or 0)
        total_seconds = float(s.get("total_seconds") or 0.0)
        eval_count = int(eval_stats.get("eval_count") or 0)
        slope_count = int(eval_stats.get("slope_count") or 0)
        rec.update({
            "bezout": data.get("bezout"),
            "terms": data.get("terms"),
            "terms_per_poly": data.get("terms_per_poly"),
            "unique_roots": unique,
            "requested_roots": int(s.get("requested_roots") or cnt),
            "coverage": unique / max(1, int(s.get("requested_roots") or cnt)),
            "success": s.get("success"),
            "trials_used": trials,
            "trials_per_root": trials / max(1, unique),
            "duplicates": int(s.get("duplicates") or 0),
            "failures": int(s.get("failures") or 0),
            "failures_per_root": int(s.get("failures") or 0) / max(1, unique),
            "extract_seconds": float(s.get("extract_seconds") or 0.0),
            "total_seconds": total_seconds,
            "seconds_per_root": total_seconds / max(1, unique),
            "best_residual": min(residuals) if residuals else None,
            "median_residual": sorted(residuals)[len(residuals)//2] if residuals else None,
            "worst_residual": max(residuals) if residuals else None,
            "eval_count": eval_count,
            "slope_count": slope_count,
            "evals_per_root": eval_count / max(1, unique),
            "slopes_per_root": slope_count / max(1, unique),
        })
    return rec


def median(xs):
    xs = [x for x in xs if x is not None and math.isfinite(float(x))]
    return statistics.median(xs) if xs else None


def main() -> int:
    rows=[]
    for n,d in CASES:
        for seed in SEEDS:
            print(f"RUN 301 ks({n},{d}) seed={seed} count={count_for(n,d)}", flush=True)
            r=run_one(n,d,seed); rows.append(r)
            print(f"  rc={r.get('returncode')} roots={r.get('unique_roots')}/{r.get('requested_roots')} t/root={r.get('seconds_per_root')} trials/root={r.get('trials_per_root')} medres={r.get('median_residual')}", flush=True)

    csv_path = OUT / "summary.csv"
    keys = sorted({k for r in rows for k in r})
    with csv_path.open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)

    md=[]
    md += ["# 301 Smale-17 diagnostic benchmark", "", "This is empirical only: no proof of average polynomial complexity, no external homotopy baseline.", "", f"Cases: {CASES}", f"Seeds: {SEEDS}", f"Base args: `{' '.join(BASE)}`", "Requested roots: `min(Bezout, 20)`", ""]
    md += ["| case | runs OK | roots | success runs | median sec/root | median trials/root | median failures/root | median residual | worst median residual |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for n,d in CASES:
        sub=[r for r in rows if r.get('n')==n and r.get('degree')==d and r.get('returncode')==0]
        roots=sum(int(r.get('unique_roots') or 0) for r in sub)
        req=sum(int(r.get('requested_roots') or 0) for r in sub)
        succ=sum(1 for r in sub if r.get('success') is True)
        med_sec=median([r.get('seconds_per_root') for r in sub])
        med_trials=median([r.get('trials_per_root') for r in sub])
        med_fail=median([r.get('failures_per_root') for r in sub])
        med_res=median([r.get('median_residual') for r in sub])
        worst_med=max([float(r.get('median_residual') or float('inf')) for r in sub], default=None)
        md.append(f"| ks({n},{d}) | {len(sub)}/{len(SEEDS)} | {roots}/{req} | {succ}/{len(sub)} | {med_sec} | {med_trials} | {med_fail} | {med_res} | {worst_med} |")

    ok=[r for r in rows if r.get('returncode')==0]
    md.append("\n## Global aggregate\n")
    md.append(f"- Runs OK: {len(ok)}/{len(rows)}")
    md.append(f"- Roots: {sum(int(r.get('unique_roots') or 0) for r in ok)}/{sum(int(r.get('requested_roots') or 0) for r in ok)}")
    md.append(f"- Successful runs: {sum(1 for r in ok if r.get('success') is True)}/{len(ok)}")
    md.append(f"- Median seconds/root: {median([r.get('seconds_per_root') for r in ok])}")
    md.append(f"- Median trials/root: {median([r.get('trials_per_root') for r in ok])}")
    md.append(f"- Median failures/root: {median([r.get('failures_per_root') for r in ok])}")
    md.append(f"- Median of run median residuals: {median([r.get('median_residual') for r in ok])}")
    md.append(f"- Median evals/root: {median([r.get('evals_per_root') for r in ok])}")
    md.append(f"- Median slopes/root: {median([r.get('slopes_per_root') for r in ok])}")

    md.append("\n## Smale-17 candidacy notes\n")
    md.append("- Positive empirical sign: full coverage on a scaling Kostlan suite supports local robustness as a research prototype.")
    md.append("- Negative empirical sign: this benchmark only asks for a capped number of roots, not certified approximate zeros for arbitrary input nor all roots.")
    md.append("- Missing theory: no average polynomial complexity proof, no alpha/gamma/mu condition analysis, no probability bound for reaching a local cubic basin, no external homotopy comparison.")

    report = OUT / "report.md"
    report.write_text('\n'.join(md))
    print(f"WROTE {csv_path}")
    print(f"WROTE {report}")
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
