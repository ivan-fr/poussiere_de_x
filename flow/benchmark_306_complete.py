#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import pathlib
import statistics
import subprocess
import sys
import time
from typing import Any

ROOT = pathlib.Path(__file__).resolve().parents[1]
ENGINE = ROOT / "flow" / "306_pandrosion_universal_atlas_hypercube_inversejet_numpy_engine.py"
OUTDIR = ROOT / "benchmarks" / "306_complete"

CASES = [(2, 2), (2, 3), (2, 4), (2, 5), (3, 2), (3, 3), (3, 4), (4, 2)]
SEEDS = list(range(5))
POOL = 2048
EPOCHS = 28
ACCEPT = 1e-8
MAX_COUNT = 12


def bezout(n: int, d: int) -> int:
    return int(d) ** int(n)


def load_case(path: pathlib.Path) -> dict[str, Any]:
    obj = json.loads(path.read_text())
    if "cases" in obj:
        return obj["cases"][0]
    return obj


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    t0 = time.time()
    for n, d in CASES:
        count = min(bezout(n, d), MAX_COUNT)
        case = f"{n},{d}"
        for seed in SEEDS:
            tag = f"ks{n}x{d}_seed{seed}"
            out = OUTDIR / f"306_{tag}.json"
            cmd = [
                sys.executable,
                str(ENGINE),
                "--cases", case,
                "--count", str(count),
                "--pool", str(POOL),
                "--epochs", str(EPOCHS),
                "--accept", str(ACCEPT),
                "--seed-index", str(seed),
                "--out", str(out),
            ]
            print(f"RUN 306 ks({n},{d}) seed={seed} count={count} pool={POOL} epochs={EPOCHS}", flush=True)
            rt0 = time.time()
            cp = subprocess.run(cmd, cwd=str(ROOT), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            elapsed = time.time() - rt0
            log = OUTDIR / f"306_{tag}.log"
            log.write_text(cp.stdout)
            if cp.returncode != 0:
                print(cp.stdout[-4000:], flush=True)
                row = {
                    "engine": "306",
                    "case": case,
                    "n": n,
                    "degree": d,
                    "seed_index": seed,
                    "requested_roots": count,
                    "unique_roots": 0,
                    "success": False,
                    "trials_used": POOL,
                    "duplicates": None,
                    "failures": None,
                    "eval_count": None,
                    "slope_count": None,
                    "seconds": elapsed,
                    "returncode": cp.returncode,
                    "out": str(out),
                    "log": str(log),
                }
            else:
                r = load_case(out)
                s = r["summary"]
                estats = s.get("eval_stats", {})
                row = {
                    "engine": "306",
                    "case": case,
                    "n": int(r.get("n", n)),
                    "degree": int(r.get("degree", d)),
                    "seed_index": seed,
                    "seed": int(r.get("seed", 0)),
                    "requested_roots": int(s["requested_roots"]),
                    "unique_roots": int(s["unique_roots"]),
                    "success": bool(s["success"]),
                    "trials_used": int(s["trials_used"]),
                    "duplicates": int(s["duplicates"]),
                    "failures": int(s["failures"]),
                    "eval_count": int(estats.get("eval_count", r.get("eval_count", 0) or 0)),
                    "slope_count": int(estats.get("slope_count", 0)),
                    "seconds": float(s["total_seconds"]),
                    "extract_seconds": float(s.get("extract_seconds", 0.0)),
                    "returncode": 0,
                    "out": str(out),
                    "log": str(log),
                }
            rows.append(row)
            print(
                f"DONE 306 ks({n},{d}) seed={seed}: roots={row['unique_roots']}/{row['requested_roots']} "
                f"success={row['success']} trials={row['trials_used']} evals={row['eval_count']} sec={row['seconds']:.3f}",
                flush=True,
            )
            (OUTDIR / "rows.json").write_text(json.dumps(rows, indent=2))

    aggregates: list[dict[str, Any]] = []
    for n, d in CASES:
        case = f"{n},{d}"
        xs = [r for r in rows if r["case"] == case]
        ok = [r for r in xs if r.get("returncode") == 0]
        def mean(key: str) -> float | None:
            vals = [r[key] for r in ok if r.get(key) is not None]
            return float(statistics.mean(vals)) if vals else None
        def med(key: str) -> float | None:
            vals = [r[key] for r in ok if r.get(key) is not None]
            return float(statistics.median(vals)) if vals else None
        aggregates.append({
            "case": case,
            "requested_total": sum(int(r["requested_roots"]) for r in xs),
            "roots_total": sum(int(r["unique_roots"]) for r in xs),
            "successes": sum(1 for r in xs if r.get("success")),
            "runs": len(xs),
            "trials_avg": mean("trials_used"),
            "trials_med": med("trials_used"),
            "evals_avg": mean("eval_count"),
            "evals_med": med("eval_count"),
            "seconds_avg": mean("seconds"),
            "seconds_med": med("seconds"),
            "failures_avg": mean("failures"),
            "duplicates_avg": mean("duplicates"),
        })

    final = {
        "engine": "306",
        "script": str(ENGINE),
        "cases": [f"{n},{d}" for n, d in CASES],
        "seeds": SEEDS,
        "parameters": {"pool": POOL, "epochs": EPOCHS, "accept": ACCEPT, "count": "min(Bezout,12)"},
        "total_wall_seconds": time.time() - t0,
        "rows": rows,
        "aggregates": aggregates,
    }
    (OUTDIR / "summary.json").write_text(json.dumps(final, indent=2))
    print("\nAGGREGATES")
    for a in aggregates:
        print(
            f"{a['case']}: success {a['successes']}/{a['runs']} roots {a['roots_total']}/{a['requested_total']} "
            f"trials_avg {a['trials_avg']:.1f} evals_avg {a['evals_avg']:.1f} sec_avg {a['seconds_avg']:.3f}",
            flush=True,
        )
    print(f"summary={OUTDIR / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
