#!/usr/bin/env python3
"""Expanded multivariate benchmark 201 vs 300.

Stdlib-only harness. Adaptive requested roots: min(Bezout, 12).
"""
from __future__ import annotations

import csv, json, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FLOW = ROOT / "flow"
OUT = ROOT / "benchmarks" / "201_vs_300_complete"
OUT.mkdir(parents=True, exist_ok=True)

SCRIPTS = {
    "201": FLOW / "201_pandrosion_logstable_optprobe_tensor_halley_numpy_engine.py",
    "300": FLOW / "300_pandrosion_atlas_geodesic_inversejet_numpy_engine.py",
}
CASES = [(2,2), (2,3), (2,4), (2,5), (3,2), (3,3), (3,4), (4,2)]
SEEDS = [0, 1, 2, 3, 4]
BASE = [
    "--pool", "2048",
    "--epochs", "28",
    "--accept", "1e-8",
    "--tol", "1e-12",
    "--line-search", "12",
    "--probe-candidates", "10",
    "--trial-timeout", "0",
    "--keep-trials", "120",
]

def count_for(n: int, d: int) -> int:
    return min(int(d ** n), 12)

def run_one(tag: str, n: int, d: int, seed: int) -> dict:
    case = f"{n},{d}"
    cnt = count_for(n,d)
    out = OUT / f"{tag}_ks{n}x{d}_seed{seed}.json"
    cmd = [sys.executable, str(SCRIPTS[tag]), "--cases", case, "--seed-index", str(seed), "--count", str(cnt), "--out", str(out), *BASE]
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    wall = time.time() - t0
    rec = {"engine": tag, "case": case, "n": n, "degree": d, "seed_index": seed, "requested_count": cnt, "returncode": proc.returncode, "wall_seconds": wall, "json": str(out), "stderr_tail": proc.stderr[-1000:], "stdout_tail": proc.stdout[-1000:]}
    if proc.returncode == 0 and out.exists():
        data = json.loads(out.read_text())
        s = data.get("summary", {})
        roots = data.get("roots", [])
        residuals = [float(r.get("residual", float("inf"))) for r in roots]
        rec.update({
            "bezout": data.get("bezout"),
            "unique_roots": s.get("unique_roots"),
            "requested_roots": s.get("requested_roots"),
            "success": s.get("success"),
            "trials_used": s.get("trials_used"),
            "duplicates": s.get("duplicates"),
            "failures": s.get("failures"),
            "extract_seconds": s.get("extract_seconds"),
            "total_seconds": s.get("total_seconds"),
            "best_residual": min(residuals) if residuals else None,
            "median_residual": sorted(residuals)[len(residuals)//2] if residuals else None,
            "worst_residual": max(residuals) if residuals else None,
        })
    return rec

def main() -> int:
    rows=[]
    for n,d in CASES:
        for seed in SEEDS:
            for tag in ("201","300"):
                print(f"RUN {tag} ks({n},{d}) seed={seed} count={count_for(n,d)}", flush=True)
                r=run_one(tag,n,d,seed); rows.append(r)
                print(f"  rc={r.get('returncode')} roots={r.get('unique_roots')}/{r.get('requested_roots')} best={r.get('best_residual')} t={r.get('total_seconds')}", flush=True)
    keys=sorted({k for r in rows for k in r})
    csv_path=OUT/'summary.csv'
    with csv_path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(rows)
    md=[]
    md += ["# Complete multivariate benchmark 201 vs 300", "", f"Cases: {CASES}", f"Seeds: {SEEDS}", f"Base args: `{' '.join(BASE)}`", "Requested roots: `min(Bezout, 12)`", ""]
    md += ["| case | seed | engine | roots | success | best residual | median residual | trials | failures | total s |", "|---|---:|---|---:|---|---:|---:|---:|---:|---:|"]
    for r in rows:
        md.append(f"| ks({r.get('case')}) | {r.get('seed_index')} | {r.get('engine')} | {r.get('unique_roots')}/{r.get('requested_roots')} | {r.get('success')} | {r.get('best_residual')} | {r.get('median_residual')} | {r.get('trials_used')} | {r.get('failures')} | {r.get('total_seconds')} |")
    md.append("\n## Aggregates\n")
    for tag in ("201","300"):
        sub=[r for r in rows if r.get('engine')==tag and r.get('returncode')==0]
        roots=sum(int(r.get('unique_roots') or 0) for r in sub); req=sum(int(r.get('requested_roots') or 0) for r in sub)
        succ=sum(1 for r in sub if r.get('success') is True); fail=sum(int(r.get('failures') or 0) for r in sub)
        t=sum(float(r.get('total_seconds') or 0) for r in sub)
        trials=sum(int(r.get('trials_used') or 0) for r in sub)
        md.append(f"- **{tag}**: roots {roots}/{req}; successful runs {succ}/{len(sub)}; trials {trials}; failures {fail}; total engine time {t:.3f}s")
    # Pairwise wins
    md.append("\n## Pairwise wins\n")
    pairs={}
    for r in rows:
        pairs.setdefault((r['case'],r['seed_index']),{})[r['engine']]=r
    win_roots={"201":0,"300":0,"tie":0}; win_time={"201":0,"300":0,"tie":0}; win_best={"201":0,"300":0,"tie":0}
    for pair in pairs.values():
        if '201' not in pair or '300' not in pair: continue
        a,b=pair['201'],pair['300']
        ar,br=int(a.get('unique_roots') or 0),int(b.get('unique_roots') or 0)
        win_roots['201' if ar>br else '300' if br>ar else 'tie']+=1
        at,bt=float(a.get('total_seconds') or 1e99),float(b.get('total_seconds') or 1e99)
        win_time['201' if at<bt else '300' if bt<at else 'tie']+=1
        av,bv=a.get('best_residual'),b.get('best_residual')
        av=float(av) if av is not None else 1e99; bv=float(bv) if bv is not None else 1e99
        win_best['201' if av<bv else '300' if bv<av else 'tie']+=1
    md.append(f"- Root-count wins: {win_roots}")
    md.append(f"- Speed wins: {win_time}")
    md.append(f"- Best-residual wins: {win_best}")
    report=OUT/'report.md'; report.write_text('\n'.join(md))
    print(f"WROTE {csv_path}"); print(f"WROTE {report}")
    return 0
if __name__=='__main__': raise SystemExit(main())
