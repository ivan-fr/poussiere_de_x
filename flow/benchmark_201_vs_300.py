#!/usr/bin/env python3
"""Benchmark 201 vs 300 on multivariate Kostlan systems.

Standalone harness using only stdlib. It launches the standalone NumPy engines
with identical cases/seeds/budgets, parses their JSON outputs, and writes a
compact CSV + Markdown report.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FLOW = ROOT / "flow"
OUT = ROOT / "benchmarks" / "201_vs_300"
OUT.mkdir(parents=True, exist_ok=True)

SCRIPTS = {
    "201": FLOW / "201_pandrosion_logstable_optprobe_tensor_halley_numpy_engine.py",
    "300": FLOW / "300_pandrosion_atlas_geodesic_inversejet_numpy_engine.py",
}

# Multivariate suite: enough to see behavior without external solvers.
CASES = ["2,2", "2,3", "2,4", "3,2", "3,3"]
SEEDS = [0, 1, 2]

COMMON = [
    "--pool", "1024",
    "--epochs", "24",
    "--count", "8",
    "--accept", "1e-8",
    "--tol", "1e-12",
    "--line-search", "12",
    "--probe-candidates", "8",
    "--trial-timeout", "0",
    "--keep-trials", "80",
]


def run_one(tag: str, case: str, seed: int) -> dict:
    script = SCRIPTS[tag]
    out = OUT / f"{tag}_ks{case.replace(',', 'x')}_seed{seed}.json"
    cmd = [sys.executable, str(script), "--cases", case, "--seed-index", str(seed), "--out", str(out), *COMMON]
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    wall = time.time() - t0
    rec = {
        "engine": tag,
        "case": case,
        "seed_index": seed,
        "returncode": proc.returncode,
        "wall_seconds": wall,
        "json": str(out),
        "stdout_tail": proc.stdout[-1200:],
        "stderr_tail": proc.stderr[-1200:],
    }
    if proc.returncode == 0 and out.exists():
        data = json.loads(out.read_text())
        summary = data.get("summary", {})
        roots = data.get("roots", [])
        best = roots[0] if roots else {}
        rec.update({
            "n": data.get("n"),
            "degree": data.get("degree"),
            "bezout": data.get("bezout"),
            "unique_roots": summary.get("unique_roots"),
            "requested_roots": summary.get("requested_roots"),
            "success": summary.get("success"),
            "trials_used": summary.get("trials_used"),
            "duplicates": summary.get("duplicates"),
            "failures": summary.get("failures"),
            "extract_seconds": summary.get("extract_seconds"),
            "total_seconds": summary.get("total_seconds"),
            "best_residual": best.get("residual"),
            "best_epochs": best.get("epochs"),
            "best_cond": best.get("cond") or best.get("slope_cond"),
            "halley_used_count_best": best.get("halley_used_count"),
        })
    return rec


def main() -> int:
    rows = []
    for case in CASES:
        for seed in SEEDS:
            for tag in ("201", "300"):
                print(f"RUN {tag} ks({case}) seed={seed}", flush=True)
                rec = run_one(tag, case, seed)
                rows.append(rec)
                if rec.get("returncode") != 0:
                    print(f"  FAIL rc={rec.get('returncode')} stderr={rec.get('stderr_tail')[-300:]}", flush=True)
                else:
                    print(f"  roots={rec.get('unique_roots')}/{rec.get('requested_roots')} best={rec.get('best_residual')} t={rec.get('total_seconds')}", flush=True)

    csv_path = OUT / "summary.csv"
    keys = sorted({k for r in rows for k in r})
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader(); w.writerows(rows)

    # Pairwise Markdown report.
    md = ["# Benchmark 201 vs 300", "", f"Cases: {', '.join(CASES)}", f"Seeds: {SEEDS}", f"Common args: `{' '.join(COMMON)}`", ""]
    md.append("| case | seed | engine | roots | success | best residual | trials | failures | total s |")
    md.append("|---|---:|---|---:|---|---:|---:|---:|---:|")
    for r in rows:
        md.append(f"| ks({r.get('case')}) | {r.get('seed_index')} | {r.get('engine')} | {r.get('unique_roots')}/{r.get('requested_roots')} | {r.get('success')} | {r.get('best_residual')} | {r.get('trials_used')} | {r.get('failures')} | {r.get('total_seconds')} |")

    # Aggregates.
    md.append("\n## Aggregates\n")
    for tag in ("201", "300"):
        sub = [r for r in rows if r.get("engine") == tag and r.get("returncode") == 0]
        if not sub:
            continue
        roots = sum(int(r.get("unique_roots") or 0) for r in sub)
        req = sum(int(r.get("requested_roots") or 0) for r in sub)
        succ = sum(1 for r in sub if r.get("success") is True)
        total_t = sum(float(r.get("total_seconds") or 0) for r in sub)
        fails = sum(int(r.get("failures") or 0) for r in sub)
        md.append(f"- **{tag}**: roots {roots}/{req}, successful runs {succ}/{len(sub)}, failures {fails}, total engine time {total_t:.3f}s")

    report = OUT / "report.md"
    report.write_text("\n".join(md))
    print(f"WROTE {csv_path}")
    print(f"WROTE {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
