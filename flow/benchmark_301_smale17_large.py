#!/usr/bin/env python3
"""Large Smale-17 diagnostic benchmark for engine 301.

Reviewer-facing empirical stress suite:
- 301 only, many more cases/seeds than the paper's first diagnostic.
- Random dense Kostlan square systems.
- Requested roots: min(Bezout, 50).
- Post-hoc analytic-Jacobian Newton alpha-lite certification.

This is still not a proof and not a replacement for homotopy baselines, but it
is large enough to test whether 301's good behavior survives scaling and many
random instances.
"""
from __future__ import annotations

import csv, json, math, statistics, subprocess, sys, time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
FLOW = ROOT / "flow"
SCRIPT = FLOW / "301_pandrosion_anchored_full_cubic_halley_numpy_engine.py"
OUT = ROOT / "benchmarks" / "301_smale17_large"
OUT.mkdir(parents=True, exist_ok=True)

CASES = [
    # degree scaling in two variables
    *[(2, d) for d in range(2, 16)],
    # degree scaling in three variables
    *[(3, d) for d in range(2, 9)],
    # moderate dimension/degree grid
    (4,2), (4,3), (4,4), (4,5),
    (5,2), (5,3), (5,4),
    (6,2), (6,3),
    # higher-dimensional quadratics
    (7,2), (8,2), (10,2),
]
SEEDS = list(range(16))
BASE = [
    "--pool", "8192",
    "--epochs", "40",
    "--accept", "1e-8",
    "--tol", "1e-12",
    "--line-search", "14",
    "--probe-candidates", "12",
    "--trial-timeout", "0",
    "--keep-trials", "40",
]
ACCEPT = 1e-8
CERT_BETA = 1e-6
CERT_POLISHED_RESIDUAL = 1e-12
STRONG_BETA = 1e-10
STRONG_POLISHED_RESIDUAL = 1e-14


def count_for(n: int, d: int) -> int:
    return min(int(d ** n), 50)


def splitmix64(x: int) -> int:
    x = (int(x) + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return (x ^ (x >> 31)) & 0xFFFFFFFFFFFFFFFF


def stable_seed(n: int, d: int, seed_index: int = 0, salt: int = 0) -> int:
    return int(splitmix64(0x50414E44524F5349 + 1000003 * n + 9176 * d + 97 * seed_index + salt) & 0x7FFFFFFF)


def compositions_leq(d: int, n: int):
    import numpy as np
    out=[]
    def rec(pos, remaining, cur):
        if pos == n-1:
            for k in range(remaining+1): out.append(tuple(cur+[k]))
            return
        for k in range(remaining+1):
            cur.append(k); rec(pos+1, remaining-k, cur); cur.pop()
    rec(0,d,[])
    return np.asarray(out, dtype=np.int16 if d < 32767 else np.int32)


def kostlan_system(n: int, d: int, seed_index: int):
    import numpy as np
    exps = compositions_leq(d,n)
    totals = np.sum(exps, axis=1).astype(np.int64)
    logfac = np.zeros(d+1, dtype=np.float64)
    acc=0.0
    for k in range(1,d+1): acc += math.log(k); logfac[k]=acc
    logs = logfac[d] - logfac[d-totals]
    for j in range(n): logs -= logfac[exps[:,j].astype(np.int64)]
    weights = np.exp(0.5*logs)
    rng = np.random.default_rng(stable_seed(n,d,seed_index))
    coeff = (rng.standard_normal((n, exps.shape[0])) + 1j*rng.standard_normal((n, exps.shape[0]))) / math.sqrt(2.0)
    coeff = (coeff * weights[None,:]).astype(np.complex128)
    return exps, coeff


def eval_and_jac(n: int, d: int, exps: Any, coeff: Any, z: Any):
    import numpy as np
    zz = np.asarray(z, dtype=np.complex128)
    powers=[]
    for zj in zz:
        p=np.empty(d+1, dtype=np.complex128); p[0]=1
        for k in range(1,d+1): p[k]=p[k-1]*zj
        powers.append(p)
    mon=np.ones(exps.shape[0], dtype=np.complex128)
    for j in range(n): mon *= powers[j][exps[:,j]]
    F = coeff @ mon
    J = np.empty((n,n), dtype=np.complex128)
    for j in range(n):
        term=np.ones(exps.shape[0], dtype=np.complex128)
        aj = exps[:,j].astype(np.float64)
        for k in range(n):
            if k == j:
                idx = exps[:,k].astype(np.int64) - 1
                fac = np.zeros(exps.shape[0], dtype=np.complex128)
                mask = idx >= 0
                fac[mask] = powers[k][idx[mask]]
                term *= fac
            else:
                term *= powers[k][exps[:,k]]
        term *= aj
        J[:,j] = coeff @ term
    return F, J


def parse_root_z(root: dict[str, Any]):
    import numpy as np
    return np.asarray([complex(float(a), float(b)) for a,b in root.get("z", [])], dtype=np.complex128)


def finite_median(xs):
    vals=[float(x) for x in xs if x is not None and math.isfinite(float(x))]
    return statistics.median(vals) if vals else None


def quantile(xs, q):
    vals=sorted(float(x) for x in xs if x is not None and math.isfinite(float(x)))
    if not vals: return None
    idx=min(len(vals)-1, max(0, int(round(q*(len(vals)-1)))))
    return vals[idx]


def certify(path: Path, n: int, d: int, seed: int) -> dict[str, Any]:
    import numpy as np
    data=json.loads(path.read_text())
    exps, coeff = kostlan_system(n,d,seed)
    rows=[]
    for r in data.get("roots", []):
        z=parse_root_z(r)
        F,J=eval_and_jac(n,d,exps,coeff,z)
        residual=float(np.linalg.norm(F))
        try:
            dx=np.linalg.solve(J, -F)
            beta=float(np.linalg.norm(dx))
            cond=float(np.linalg.cond(J))
            z1=z+dx
            F1,J1=eval_and_jac(n,d,exps,coeff,z1)
            polished=float(np.linalg.norm(F1))
            beta1=float(np.linalg.norm(np.linalg.solve(J1, -F1)))
        except Exception:
            beta=float('inf'); cond=float('inf'); polished=float('inf'); beta1=float('inf')
        rows.append({
            "residual_check": residual,
            "beta": beta,
            "cond_jac": cond,
            "polished_residual": polished,
            "beta_after_one_newton": beta1,
            "certified": bool(residual <= ACCEPT and beta <= CERT_BETA and polished <= CERT_POLISHED_RESIDUAL),
            "strong": bool(residual <= ACCEPT and beta <= STRONG_BETA and polished <= STRONG_POLISHED_RESIDUAL),
        })
    return {
        "certified_roots": sum(1 for c in rows if c["certified"]),
        "strong_roots": sum(1 for c in rows if c["strong"]),
        "median_beta": finite_median(c["beta"] for c in rows),
        "p90_beta": quantile([c["beta"] for c in rows], 0.9),
        "median_polished_residual": finite_median(c["polished_residual"] for c in rows),
        "p90_polished_residual": quantile([c["polished_residual"] for c in rows], 0.9),
        "max_polished_residual": max([c["polished_residual"] for c in rows], default=None),
        "median_cond_jac": finite_median(c["cond_jac"] for c in rows),
        "p90_cond_jac": quantile([c["cond_jac"] for c in rows], 0.9),
    }


def run_one(n: int, d: int, seed: int) -> dict[str, Any]:
    cnt=count_for(n,d)
    out = OUT / f"301_ks{n}x{d}_seed{seed}.json"
    cmd=[sys.executable, str(SCRIPT), "--cases", f"{n},{d}", "--seed-index", str(seed), "--count", str(cnt), "--out", str(out), *BASE]
    t0=time.time()
    proc=subprocess.run(cmd, cwd=str(ROOT), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    wall=time.time()-t0
    rec={"engine":"301","case":f"{n},{d}","n":n,"degree":d,"seed_index":seed,"requested_count":cnt,"returncode":proc.returncode,"wall_seconds":wall,"json":str(out),"stderr_tail":proc.stderr[-1200:],"stdout_tail":proc.stdout[-1200:]}
    if proc.returncode == 0 and out.exists():
        data=json.loads(out.read_text())
        s=data.get("summary",{})
        roots=data.get("roots",[])
        residuals=[float(r.get("residual", float("inf"))) for r in roots]
        unique=int(s.get("unique_roots") or 0)
        req=int(s.get("requested_roots") or cnt)
        total=float(s.get("total_seconds") or 0.0)
        trials=int(s.get("trials_used") or 0)
        eval_stats=s.get("eval_stats",{}) or {}
        cert=certify(out,n,d,seed)
        rec.update({
            "bezout": data.get("bezout"), "terms": data.get("terms"), "terms_per_poly": data.get("terms_per_poly"),
            "unique_roots": unique, "requested_roots": req, "success": s.get("success"), "coverage": unique/max(1,req),
            "trials_used": trials, "trials_per_root": trials/max(1,unique),
            "duplicates": int(s.get("duplicates") or 0), "failures": int(s.get("failures") or 0), "failures_per_root": int(s.get("failures") or 0)/max(1,unique),
            "total_seconds": total, "seconds_per_root": total/max(1,unique),
            "best_residual": min(residuals) if residuals else None, "median_residual": finite_median(residuals), "p90_residual": quantile(residuals,0.9), "worst_residual": max(residuals) if residuals else None,
            "eval_count": int(eval_stats.get("eval_count") or 0), "slope_count": int(eval_stats.get("slope_count") or 0),
            "evals_per_root": int(eval_stats.get("eval_count") or 0)/max(1,unique), "slopes_per_root": int(eval_stats.get("slope_count") or 0)/max(1,unique),
            **cert,
        })
    return rec


def main() -> int:
    rows=[]
    started=time.time()
    total_runs=len(CASES)*len(SEEDS)
    k=0
    for n,d in CASES:
        for seed in SEEDS:
            k += 1
            print(f"RUN {k}/{total_runs} 301 ks({n},{d}) seed={seed} count={count_for(n,d)}", flush=True)
            r=run_one(n,d,seed); rows.append(r)
            print(f"  rc={r.get('returncode')} roots={r.get('unique_roots')}/{r.get('requested_roots')} cert={r.get('certified_roots')} strong={r.get('strong_roots')} sec/root={r.get('seconds_per_root')} medres={r.get('median_residual')}", flush=True)
            # incrementally save so long runs are still useful if interrupted
            write_outputs(rows, partial=True, elapsed=time.time()-started)
    write_outputs(rows, partial=False, elapsed=time.time()-started)
    return 0


def write_outputs(rows: list[dict[str, Any]], partial: bool, elapsed: float) -> None:
    csv_path=OUT/("summary_partial.csv" if partial else "summary.csv")
    keys=sorted({k for r in rows for k in r})
    with csv_path.open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    ok=[r for r in rows if r.get('returncode')==0]
    md=[]
    md += ["# Large 301 Smale-17 diagnostic benchmark", "", "Empirical stress benchmark for reviewer scale concerns. This is not a proof and not a homotopy baseline.", "", f"Cases: {CASES}", f"Seeds: {SEEDS}", f"Completed runs: {len(rows)}/{len(CASES)*len(SEEDS)}", f"Elapsed wall time: {elapsed:.1f}s", f"Base args: `{' '.join(BASE)}`", "Requested roots: `min(Bezout, 50)`", ""]
    roots=sum(int(r.get('unique_roots') or 0) for r in ok); req=sum(int(r.get('requested_roots') or 0) for r in ok)
    cert=sum(int(r.get('certified_roots') or 0) for r in ok); strong=sum(int(r.get('strong_roots') or 0) for r in ok)
    md += ["## Global aggregate", "", f"- Runs OK: {len(ok)}/{len(rows)}", f"- Roots: {roots}/{req}", f"- Certified alpha-lite roots: {cert}/{roots if roots else 0}", f"- Strong roots: {strong}/{roots if roots else 0}", f"- Median seconds/root: {finite_median([r.get('seconds_per_root') for r in ok])}", f"- P90 seconds/root: {quantile([r.get('seconds_per_root') for r in ok],0.9)}", f"- Median trials/root: {finite_median([r.get('trials_per_root') for r in ok])}", f"- P90 trials/root: {quantile([r.get('trials_per_root') for r in ok],0.9)}", f"- Median residual: {finite_median([r.get('median_residual') for r in ok])}", f"- P90 run-median residual: {quantile([r.get('median_residual') for r in ok],0.9)}", f"- Median polished residual: {finite_median([r.get('median_polished_residual') for r in ok])}", f"- Median cond(J): {finite_median([r.get('median_cond_jac') for r in ok])}", ""]
    md += ["## Per-case aggregate", "", "| case | runs | roots | certified | strong | median sec/root | p90 sec/root | median trials/root | median residual | max run-median residual |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for n,d in CASES:
        sub=[r for r in ok if r.get('n')==n and r.get('degree')==d]
        if not sub: continue
        roots=sum(int(r.get('unique_roots') or 0) for r in sub); req=sum(int(r.get('requested_roots') or 0) for r in sub)
        cert=sum(int(r.get('certified_roots') or 0) for r in sub); strong=sum(int(r.get('strong_roots') or 0) for r in sub)
        medres=[float(r.get('median_residual')) for r in sub if r.get('median_residual') not in (None,"")]
        md.append(f"| ks({n},{d}) | {len(sub)}/{len(SEEDS)} | {roots}/{req} | {cert}/{roots} | {strong}/{roots} | {finite_median([r.get('seconds_per_root') for r in sub])} | {quantile([r.get('seconds_per_root') for r in sub],0.9)} | {finite_median([r.get('trials_per_root') for r in sub])} | {finite_median(medres)} | {max(medres) if medres else None} |")
    report=OUT/("report_partial.md" if partial else "report.md")
    report.write_text('\n'.join(md))
    if not partial:
        print(f"WROTE {csv_path}")
        print(f"WROTE {report}")

if __name__ == '__main__':
    raise SystemExit(main())
