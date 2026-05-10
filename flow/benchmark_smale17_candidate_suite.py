#!/usr/bin/env python3
"""Smale-17 candidate benchmark suite: 200 vs 201 vs 301 + alpha-lite certificates.

Empirical only. This script deliberately does NOT use Docker and does NOT depend
on SciPy/homotopy packages. It runs the three NumPy engines on the same random
Kostlan systems, then reconstructs those systems and applies an analytic
Jacobian Newton-certificate check to each reported root.

Certificate is intentionally conservative but not a formal Smale alpha proof:
  - residual <= accept threshold
  - Newton correction beta = ||J(z)^-1 F(z)|| is small
  - one analytic Newton polish reaches tiny residual
This tells us whether the output is a credible approximate zero candidate.
"""
from __future__ import annotations

import csv, json, math, statistics, subprocess, sys, time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
FLOW = ROOT / "flow"
OUT = ROOT / "benchmarks" / "smale17_candidate_suite"
OUT.mkdir(parents=True, exist_ok=True)

ENGINES = {
    "200": FLOW / "200_pandrosion_118_backbone_gated_halley_numpy_engine.py",
    "201": FLOW / "201_pandrosion_logstable_optprobe_tensor_halley_numpy_engine.py",
    "301": FLOW / "301_pandrosion_anchored_full_cubic_halley_numpy_engine.py",
}
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
    "--keep-trials", "80",
]
ACCEPT = 1e-8
CERT_BETA = 1e-6
CERT_POLISHED_RESIDUAL = 1e-12
STRONG_BETA = 1e-10
STRONG_POLISHED_RESIDUAL = 1e-14


def count_for(n: int, d: int) -> int:
    return min(int(d ** n), 20)


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


def certify_json(path: Path, n: int, d: int, seed: int) -> dict[str, Any]:
    import numpy as np
    data=json.loads(path.read_text())
    exps, coeff = kostlan_system(n,d,seed)
    certs=[]
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
            beta=float("inf"); cond=float("inf"); polished=float("inf"); beta1=float("inf")
        certified = bool(residual <= ACCEPT and beta <= CERT_BETA and polished <= CERT_POLISHED_RESIDUAL)
        strong = bool(residual <= ACCEPT and beta <= STRONG_BETA and polished <= STRONG_POLISHED_RESIDUAL)
        certs.append({"residual_check": residual, "beta": beta, "cond_jac": cond, "polished_residual": polished, "beta_after_one_newton": beta1, "certified": certified, "strong": strong})
    return {
        "certificates": certs,
        "certified_roots": sum(1 for c in certs if c["certified"]),
        "strong_roots": sum(1 for c in certs if c["strong"]),
        "median_beta": median([c["beta"] for c in certs]),
        "median_polished_residual": median([c["polished_residual"] for c in certs]),
        "max_polished_residual": max([c["polished_residual"] for c in certs], default=None),
        "median_cond_jac": median([c["cond_jac"] for c in certs]),
    }


def median(xs):
    vals=[float(x) for x in xs if x is not None and math.isfinite(float(x))]
    return statistics.median(vals) if vals else None


def run_one(engine: str, n: int, d: int, seed: int) -> dict[str, Any]:
    cnt=count_for(n,d)
    out = OUT / f"{engine}_ks{n}x{d}_seed{seed}.json"
    cmd=[sys.executable, str(ENGINES[engine]), "--cases", f"{n},{d}", "--seed-index", str(seed), "--count", str(cnt), "--out", str(out), *BASE]
    t0=time.time()
    proc=subprocess.run(cmd, cwd=str(ROOT), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    wall=time.time()-t0
    rec={"engine":engine,"case":f"{n},{d}","n":n,"degree":d,"seed_index":seed,"requested_count":cnt,"returncode":proc.returncode,"wall_seconds":wall,"json":str(out),"stderr_tail":proc.stderr[-1000:],"stdout_tail":proc.stdout[-1000:]}
    if proc.returncode == 0 and out.exists():
        data=json.loads(out.read_text())
        s=data.get("summary",{})
        roots=data.get("roots",[])
        residuals=[float(r.get("residual", float("inf"))) for r in roots]
        unique=int(s.get("unique_roots") or 0)
        total=float(s.get("total_seconds") or 0.0)
        trials=int(s.get("trials_used") or 0)
        eval_stats=s.get("eval_stats",{}) or {}
        cert=certify_json(out,n,d,seed)
        rec.update({
            "bezout": data.get("bezout"), "terms": data.get("terms"), "terms_per_poly": data.get("terms_per_poly"),
            "unique_roots": unique, "requested_roots": int(s.get("requested_roots") or cnt), "success": s.get("success"),
            "coverage": unique/max(1,int(s.get("requested_roots") or cnt)), "trials_used": trials, "trials_per_root": trials/max(1,unique),
            "duplicates": int(s.get("duplicates") or 0), "failures": int(s.get("failures") or 0), "failures_per_root": int(s.get("failures") or 0)/max(1,unique),
            "total_seconds": total, "seconds_per_root": total/max(1,unique),
            "best_residual": min(residuals) if residuals else None, "median_residual": sorted(residuals)[len(residuals)//2] if residuals else None, "worst_residual": max(residuals) if residuals else None,
            "eval_count": int(eval_stats.get("eval_count") or 0), "slope_count": int(eval_stats.get("slope_count") or 0),
            "evals_per_root": int(eval_stats.get("eval_count") or 0)/max(1,unique), "slopes_per_root": int(eval_stats.get("slope_count") or 0)/max(1,unique),
            **{k:v for k,v in cert.items() if k != "certificates"},
        })
        (OUT / f"{engine}_ks{n}x{d}_seed{seed}_certificates.json").write_text(json.dumps(cert, indent=2))
    return rec


def main() -> int:
    rows=[]
    for n,d in CASES:
        for seed in SEEDS:
            for eng in ("200","201","301"):
                print(f"RUN {eng} ks({n},{d}) seed={seed} count={count_for(n,d)}", flush=True)
                r=run_one(eng,n,d,seed); rows.append(r)
                print(f"  rc={r.get('returncode')} roots={r.get('unique_roots')}/{r.get('requested_roots')} cert={r.get('certified_roots')} strong={r.get('strong_roots')} sec/root={r.get('seconds_per_root')} medres={r.get('median_residual')}", flush=True)
    keys=sorted({k for r in rows for k in r})
    csv_path=OUT/"summary.csv"
    with csv_path.open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)

    md=[]
    md += ["# Smale-17 candidate suite: 200 vs 201 vs 301", "", "Empirical benchmark + analytic-Jacobian Newton certificate check. This is not a Smale alpha proof and not a homotopy baseline.", "", f"Cases: {CASES}", f"Seeds: {SEEDS}", f"Base args: `{' '.join(BASE)}`", "Requested roots: `min(Bezout, 20)`", ""]
    md += ["## Global aggregate", "", "| engine | runs OK | roots | certified | strong | median sec/root | median trials/root | median failures/root | median residual | median polished residual | median cond(J) |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for eng in ("200","201","301"):
        sub=[r for r in rows if r.get('engine')==eng and r.get('returncode')==0]
        roots=sum(int(r.get('unique_roots') or 0) for r in sub); req=sum(int(r.get('requested_roots') or 0) for r in sub)
        cert=sum(int(r.get('certified_roots') or 0) for r in sub); strong=sum(int(r.get('strong_roots') or 0) for r in sub)
        md.append(f"| {eng} | {len(sub)}/{len(CASES)*len(SEEDS)} | {roots}/{req} | {cert}/{roots} | {strong}/{roots} | {median([r.get('seconds_per_root') for r in sub])} | {median([r.get('trials_per_root') for r in sub])} | {median([r.get('failures_per_root') for r in sub])} | {median([r.get('median_residual') for r in sub])} | {median([r.get('median_polished_residual') for r in sub])} | {median([r.get('median_cond_jac') for r in sub])} |")

    md += ["", "## Pairwise case/seed wins", ""]
    groups={}
    for r in rows: groups.setdefault((r.get('case'),r.get('seed_index')),{})[r.get('engine')]=r
    for metric, lower in [("unique_roots", False),("certified_roots", False),("seconds_per_root", True),("median_residual", True),("median_polished_residual", True),("trials_per_root", True)]:
        wins={"200":0,"201":0,"301":0,"tie":0}
        for g in groups.values():
            if not all(e in g and g[e].get(metric) not in (None,"") for e in ENGINES): continue
            vals={e:float(g[e].get(metric)) for e in ENGINES}
            best=min(vals.values()) if lower else max(vals.values())
            got=[e for e,v in vals.items() if abs(v-best) <= max(1e-15, abs(best)*1e-12)]
            if len(got)==1: wins[got[0]]+=1
            else: wins['tie']+=1
        md.append(f"- {metric} ({'lower' if lower else 'higher'} better): {wins}")

    md += ["", "## Per-case aggregate", "", "| case | engine | roots | certified | median sec/root | median trials/root | median residual | max median residual |", "|---|---|---:|---:|---:|---:|---:|---:|"]
    for n,d in CASES:
        for eng in ("200","201","301"):
            sub=[r for r in rows if r.get('engine')==eng and r.get('n')==n and r.get('degree')==d and r.get('returncode')==0]
            roots=sum(int(r.get('unique_roots') or 0) for r in sub); req=sum(int(r.get('requested_roots') or 0) for r in sub)
            cert=sum(int(r.get('certified_roots') or 0) for r in sub)
            medres=[float(r.get('median_residual')) for r in sub if r.get('median_residual') not in (None,"")]
            md.append(f"| ks({n},{d}) | {eng} | {roots}/{req} | {cert}/{roots} | {median([r.get('seconds_per_root') for r in sub])} | {median([r.get('trials_per_root') for r in sub])} | {median(medres)} | {max(medres) if medres else None} |")

    md += ["", "## Interpretation guardrails", "", "- Passing this benchmark means: good empirical approximate-zero production on capped Kostlan samples.", "- It does not mean: Smale 17 solved, average-polynomial proof, certified alpha theory, or superiority to homotopy continuation.", "- External homotopy baseline was not run here because no local `julia`, `phc`, or `bertini` executable is available in this environment."]
    report=OUT/"report.md"; report.write_text('\n'.join(md))
    print(f"WROTE {csv_path}"); print(f"WROTE {report}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
