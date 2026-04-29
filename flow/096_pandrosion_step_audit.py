"""
FLOW 096 -- step audit for the 092 Pandrosion solver.

This script intentionally does not change the 092 algorithm.  It runs the same
stack and records the worker JSON artifacts so we can inspect, per system and
per path:

  * homotopy steps,
  * Pandrosion corrector epochs,
  * wall time,
  * final t, residual, status and candidate presence.

The goal is diagnostic: find whether hard cases are expensive because every path
takes many steps, because a few paths dominate, or because recovery stages retry
the same indices too broadly.
"""
from __future__ import annotations

import csv
import inspect
import importlib.util
import json
import math
import os
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Any, Sequence

HERE = Path(__file__).resolve().parent
FLOW092_PATH = HERE / "092_pandrosion_resultant_pairing.py"


def _load_092():
    spec = importlib.util.spec_from_file_location("flow092_for_096_audit", str(FLOW092_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW092_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow092_for_096_audit"] = mod
    spec.loader.exec_module(mod)
    return mod


def _consume_096_args(argv: list[str]) -> list[str]:
    out = [argv[0]]
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--steps-csv" and i + 1 < len(argv):
            os.environ["PANDROSION_096_STEPS_CSV"] = argv[i + 1]
            i += 2
            continue
        if arg == "--steps-summary-csv" and i + 1 < len(argv):
            os.environ["PANDROSION_096_STEPS_SUMMARY_CSV"] = argv[i + 1]
            i += 2
            continue
        if arg == "--steps-stage-csv" and i + 1 < len(argv):
            os.environ["PANDROSION_096_STEPS_STAGE_CSV"] = argv[i + 1]
            i += 2
            continue
        if arg == "--steps-json" and i + 1 < len(argv):
            os.environ["PANDROSION_096_STEPS_JSON"] = argv[i + 1]
            i += 2
            continue
        if arg == "--steps-md" and i + 1 < len(argv):
            os.environ["PANDROSION_096_STEPS_MD"] = argv[i + 1]
            i += 2
            continue
        if arg == "--steps-top" and i + 1 < len(argv):
            os.environ["PANDROSION_096_STEPS_TOP"] = argv[i + 1]
            i += 2
            continue
        out.append(arg)
        i += 1
    return out


def _env_int(name: str, default: int) -> int:
    try:
        return int(float(os.environ.get(name, str(default))))
    except Exception:
        return default


def _float(v: Any, default: float = 0.0) -> float:
    try:
        out = float(v)
        return out if math.isfinite(out) else default
    except Exception:
        return default


def _int(v: Any, default: int = 0) -> int:
    try:
        return int(float(v))
    except Exception:
        return default


def _q(vals: Sequence[float], q: float) -> float:
    clean = sorted(v for v in vals if math.isfinite(v))
    if not clean:
        return 0.0
    pos = max(0.0, min(1.0, float(q))) * (len(clean) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return clean[lo]
    t = pos - lo
    return (1.0 - t) * clean[lo] + t * clean[hi]


def _default_paths(args) -> dict[str, Path]:
    outdir = Path(getattr(args, "outdir", "") or "/tmp/096_step_audit")
    outdir.mkdir(parents=True, exist_ok=True)
    return {
        "path_csv": Path(os.environ.get("PANDROSION_096_STEPS_CSV", str(outdir / "096_steps_paths.csv"))),
        "summary_csv": Path(os.environ.get("PANDROSION_096_STEPS_SUMMARY_CSV", str(outdir / "096_steps_summary.csv"))),
        "stage_csv": Path(os.environ.get("PANDROSION_096_STEPS_STAGE_CSV", str(outdir / "096_steps_stage_summary.csv"))),
        "json": Path(os.environ.get("PANDROSION_096_STEPS_JSON", str(outdir / "096_steps.json"))),
        "md": Path(os.environ.get("PANDROSION_096_STEPS_MD", str(outdir / "096_steps.md"))),
    }


def _dedupe_batches(batches: Sequence[Any]) -> list[Any]:
    seen: set[str] = set()
    out: list[Any] = []
    for b in batches:
        path = str(getattr(b, "path_json", "") or "")
        key = path or repr(b)
        if key in seen:
            continue
        seen.add(key)
        out.append(b)
    return out


def _load_rows_from_batches(batches: Sequence[Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for b in _dedupe_batches(batches):
        path = str(getattr(b, "path_json", "") or "")
        if not path:
            continue
        try:
            payload = json.loads(Path(path).read_text())
        except Exception:
            continue
        batch_meta = {
            "stage": getattr(b, "stage", ""),
            "batch_retry": getattr(b, "retry", ""),
            "batch_policy": getattr(b, "policy", ""),
            "batch_indices": getattr(b, "indices", ""),
            "batch_status": getattr(b, "status", ""),
            "batch_seconds": getattr(b, "seconds", 0.0),
            "path_json": path,
        }
        for row in payload.get("rows", []):
            if not isinstance(row, dict):
                continue
            item = dict(row)
            item.update(batch_meta)
            rows.append(item)
    return rows


def _case_key(row: dict[str, Any]) -> str:
    if row.get("case") and not row.get("family"):
        return str(row.get("case"))
    return f"{row.get('family', '')} ({_int(row.get('n'))},{_int(row.get('d'))}) seed{_int(row.get('seed', 0))}"


def _summarize_case(rows: Sequence[dict[str, Any]], summary_by_case: dict[str, Any]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[_case_key(row)].append(row)

    summaries: list[dict[str, Any]] = []
    for case, rs in sorted(grouped.items()):
        steps = [_float(r.get("steps")) for r in rs]
        epochs = [_float(r.get("epochs")) for r in rs]
        seconds = [_float(r.get("seconds")) for r in rs]
        residuals = [_float(r.get("residual"), float("inf")) for r in rs if math.isfinite(_float(r.get("residual"), float("inf")))]
        statuses = Counter(str(r.get("status", "")) for r in rs)
        stages = Counter(str(r.get("stage", "")) for r in rs)
        kinds = Counter(str(r.get("kind", "homotopy")) for r in rs)
        candidates = sum(1 for r in rs if r.get("z") is not None or r.get("has_root"))
        ok_paths = sum(1 for r in rs if bool(r.get("ok")))
        final = summary_by_case.get(case, {})
        roots = _int(final.get("roots", 0))
        bezout = _int(final.get("bezout", max(len(rs), 1)))
        total_steps = sum(steps)
        total_epochs = sum(epochs)
        summaries.append({
            "case": case,
            "bezout": bezout,
            "final_roots": roots,
            "coverage": _float(final.get("coverage", roots / max(1, bezout))),
            "tracked_rows": len(rs),
            "homotopy_rows": int(kinds.get("homotopy", 0)),
            "polish_rows": int(kinds.get("polish", 0)),
            "candidate_rows": candidates,
            "ok_path_rows": ok_paths,
            "total_steps": int(total_steps),
            "total_epochs": int(total_epochs),
            "mean_steps": total_steps / max(1, len(rs)),
            "median_steps": _q(steps, 0.50),
            "p90_steps": _q(steps, 0.90),
            "p95_steps": _q(steps, 0.95),
            "p99_steps": _q(steps, 0.99),
            "max_steps": max(steps) if steps else 0.0,
            "mean_epochs": total_epochs / max(1, len(rs)),
            "median_epochs": _q(epochs, 0.50),
            "p90_epochs": _q(epochs, 0.90),
            "max_epochs": max(epochs) if epochs else 0.0,
            "total_path_seconds": sum(seconds),
            "mean_path_seconds": sum(seconds) / max(1, len(rs)),
            "p95_path_seconds": _q(seconds, 0.95),
            "max_path_seconds": max(seconds) if seconds else 0.0,
            "steps_per_final_root": total_steps / max(1, roots),
            "epochs_per_final_root": total_epochs / max(1, roots),
            "max_residual_row": max(residuals) if residuals else float("inf"),
            "statuses": json.dumps(dict(statuses), sort_keys=True),
            "stages": json.dumps(dict(stages), sort_keys=True),
            "kinds": json.dumps(dict(kinds), sort_keys=True),
        })
    return summaries


def _summarize_stage(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(_case_key(row), str(row.get("kind", "")), str(row.get("stage", "")))].append(row)

    out: list[dict[str, Any]] = []
    for (case, kind, stage), rs in sorted(grouped.items()):
        steps = [_float(r.get("steps")) for r in rs]
        epochs = [_float(r.get("epochs")) for r in rs]
        seconds = [_float(r.get("seconds")) for r in rs]
        statuses = Counter(str(r.get("status", "")) for r in rs)
        out.append({
            "case": case,
            "kind": kind,
            "stage": stage,
            "rows": len(rs),
            "ok_rows": sum(1 for r in rs if bool(r.get("ok"))),
            "candidate_rows": sum(1 for r in rs if r.get("has_root") or r.get("z") is not None),
            "total_steps": int(sum(steps)),
            "mean_steps": sum(steps) / max(1, len(rs)),
            "p95_steps": _q(steps, 0.95),
            "max_steps": max(steps) if steps else 0.0,
            "total_epochs": int(sum(epochs)),
            "mean_epochs": sum(epochs) / max(1, len(rs)),
            "p95_epochs": _q(epochs, 0.95),
            "max_epochs": max(epochs) if epochs else 0.0,
            "total_seconds": sum(seconds),
            "mean_seconds": sum(seconds) / max(1, len(rs)),
            "max_seconds": max(seconds) if seconds else 0.0,
            "statuses": json.dumps(dict(statuses), sort_keys=True),
        })
    return out


def _path_csv_rows(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for r in rows:
        steps = _int(r.get("steps"))
        track_epochs = _int(r.get("epochs"))
        corrector_epochs = _int(r.get("corrector_epochs_total", track_epochs))
        out.append({
            "case": _case_key(r),
            "kind": r.get("kind", "homotopy" if steps else "polish"),
            "stage": r.get("stage", ""),
            "retry": r.get("retry", r.get("batch_retry", "")),
            "policy": r.get("policy", r.get("batch_policy", "")),
            "idx": r.get("idx", ""),
            "has_root": bool(r.get("z") is not None or r.get("has_root")),
            "ok": bool(r.get("ok")),
            "status": r.get("status", ""),
            "t": _float(r.get("t")),
            "steps": steps,
            "epochs": corrector_epochs,
            "track_epochs": track_epochs,
            "post_track_epochs": max(0, corrector_epochs - track_epochs),
            "corrector_calls": _int(r.get("corrector_calls_total", r.get("corrector_calls", 0))),
            "seconds": _float(r.get("seconds")),
            "residual": _float(r.get("residual"), float("inf")),
            "residual_scaled": _float(r.get("residual_scaled"), float("inf")),
            "path_json": r.get("path_json", ""),
        })
    return out


def _summary_map(summaries: Sequence[Any]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for s in summaries:
        d = asdict(s) if hasattr(s, "__dataclass_fields__") else dict(s)
        key = f"{d.get('family', '')} ({_int(d.get('n'))},{_int(d.get('d'))}) seed{_int(d.get('seed', 0))}"
        out[key] = d
    return out


def _write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _write_md(path: Path, case_summaries: Sequence[dict[str, Any]], stage_summaries: Sequence[dict[str, Any]], path_rows: Sequence[dict[str, Any]], top_n: int) -> None:
    lines: list[str] = []
    lines.append("# 096 Pandrosion Step Audit\n\n")
    lines.append("## System Summary\n\n")
    lines.append("| case | Bezout | roots | rows | homotopy | polish | candidates | total steps | mean steps | p95 steps | max steps | total epochs | mean sec/path | max sec/path | status counts |\n")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
    for s in case_summaries:
        lines.append(
            f"| {s['case']} | {s['bezout']} | {s['final_roots']} | {s['tracked_rows']} | "
            f"{s['homotopy_rows']} | {s['polish_rows']} | {s['candidate_rows']} | {s['total_steps']} | {s['mean_steps']:.1f} | "
            f"{s['p95_steps']:.1f} | {s['max_steps']:.0f} | {s['total_epochs']} | "
            f"{s['mean_path_seconds']:.3f} | {s['max_path_seconds']:.3f} | `{s['statuses']}` |\n"
        )
    lines.append("\n## Stage Summary\n\n")
    lines.append("| case | kind | stage | rows | ok | candidates | total steps | mean steps | p95 steps | max steps | total epochs | mean sec | max sec | status counts |\n")
    lines.append("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
    for s in stage_summaries:
        lines.append(
            f"| {s['case']} | {s['kind']} | {s['stage']} | {s['rows']} | {s['ok_rows']} | "
            f"{s['candidate_rows']} | {s['total_steps']} | {s['mean_steps']:.1f} | "
            f"{s['p95_steps']:.1f} | {s['max_steps']:.0f} | {s['total_epochs']} | "
            f"{s['mean_seconds']:.3f} | {s['max_seconds']:.3f} | `{s['statuses']}` |\n"
        )
    lines.append("\n## Slowest / Longest Paths\n\n")
    top = sorted(path_rows, key=lambda r: (_int(r.get("steps")), _float(r.get("seconds"))), reverse=True)[:top_n]
    lines.append("| case | kind | stage | idx | status | has root | t | steps | epochs | seconds | residual |\n")
    lines.append("|---|---|---|---:|---|---:|---:|---:|---:|---:|---:|\n")
    for r in top:
        lines.append(
            f"| {r['case']} | {r['kind']} | {r['stage']} | {r['idx']} | {r['status']} | {int(bool(r['has_root']))} | "
            f"{_float(r['t']):.6f} | {_int(r['steps'])} | {_int(r['epochs'])} | "
            f"{_float(r['seconds']):.4f} | {_float(r['residual'], float('inf')):.2e} |\n"
        )
    lines.append("\nNote: homotopy rows have path steps. Polish rows are recovery attempts; they have zero homotopy steps but still count Pandrosion corrector epochs.\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines))


def write_step_audit(mod, summaries, batches, args) -> None:
    captured = getattr(mod, "_096_captured_batches", [])
    all_batches = _dedupe_batches([*captured, *list(batches)])
    rows = _load_rows_from_batches(all_batches)
    rows.extend(getattr(mod, "_096_polish_rows", []))
    path_rows = _path_csv_rows(rows)
    summary_by_case = _summary_map(summaries)
    case_summaries = _summarize_case(path_rows, summary_by_case)
    stage_summaries = _summarize_stage(path_rows)
    paths = _default_paths(args)
    top_n = max(1, _env_int("PANDROSION_096_STEPS_TOP", 20))

    _write_csv(paths["path_csv"], path_rows)
    _write_csv(paths["summary_csv"], case_summaries)
    _write_csv(paths["stage_csv"], stage_summaries)
    payload = {
        "case_summaries": case_summaries,
        "stage_summaries": stage_summaries,
        "path_rows": path_rows,
        "batch_count": len(all_batches),
    }
    paths["json"].parent.mkdir(parents=True, exist_ok=True)
    paths["json"].write_text(json.dumps(payload, indent=2, allow_nan=True))
    _write_md(paths["md"], case_summaries, stage_summaries, path_rows, top_n)

    for s in case_summaries:
        print(
            "096-step-audit "
            f"{s['case']}: rows={s['tracked_rows']} roots={s['final_roots']}/{s['bezout']} "
            f"total_steps={s['total_steps']} mean={s['mean_steps']:.1f} "
            f"p95={s['p95_steps']:.1f} max={s['max_steps']:.0f} "
            f"total_epochs={s['total_epochs']} mean_sec={s['mean_path_seconds']:.3f}",
            flush=True,
        )
    print(f"096-step-audit wrote {paths['summary_csv']}, {paths['stage_csv']} and {paths['md']}", flush=True)


def install_audit_hooks(mod) -> None:
    mod._096_captured_batches = []
    mod._096_polish_rows = []
    mod._096_current_case = {}
    mod._096_corrector_contexts = []

    original_corrector = mod.m.corrector

    def corrector_audit(*args, **kwargs):
        z, ok, epochs = original_corrector(*args, **kwargs)
        for ctx in getattr(mod, "_096_corrector_contexts", []):
            ctx["corrector_calls_total"] = int(ctx.get("corrector_calls_total", 0)) + 1
            ctx["corrector_epochs_total"] = int(ctx.get("corrector_epochs_total", 0)) + _int(epochs)
        return z, ok, epochs

    mod.m.corrector = corrector_audit

    original_track_one_index = mod.track_one_index

    def track_one_index_audit(*args, **kwargs):
        ctx = {"corrector_calls_total": 0, "corrector_epochs_total": 0}
        mod._096_corrector_contexts.append(ctx)
        try:
            row = original_track_one_index(*args, **kwargs)
        finally:
            mod._096_corrector_contexts.pop()
        row = dict(row)
        row["kind"] = "homotopy"
        row["corrector_calls_total"] = int(ctx.get("corrector_calls_total", 0))
        row["corrector_epochs_total"] = int(ctx.get("corrector_epochs_total", row.get("epochs", 0)))
        return row

    mod.track_one_index = track_one_index_audit

    original_polish = mod.m.polish_070

    def infer_polish_stage() -> str:
        frame = inspect.currentframe()
        if frame is not None:
            frame = frame.f_back
        while frame is not None:
            filename = frame.f_code.co_filename
            if filename.endswith("092_pandrosion_resultant_pairing.py"):
                return "092-pair-polish"
            if filename.endswith("091_pandrosion_projective_riemann.py"):
                return "091-projective-polish"
            if filename.endswith("090_pandrosion_multivariate_completion.py"):
                return "090-multivariate-polish"
            if filename.endswith("089_pandrosion_resultant_recovery.py"):
                return "089-resultant-polish"
            frame = frame.f_back
        return "polish"

    def polish_070_audit(target, z, tol, quad_cap, modes, rescue_modes, deadline):
        p0 = time.time()
        stage = infer_polish_stage()
        ctx = {"corrector_calls_total": 0, "corrector_epochs_total": 0}
        mod._096_corrector_contexts.append(ctx)
        try:
            polished = original_polish(target, z, tol, quad_cap, modes, rescue_modes, deadline)
        finally:
            mod._096_corrector_contexts.pop()
        residual = mod.m.residual_norm(target, polished) if polished is not None else float("inf")
        case_ctx = dict(getattr(mod, "_096_current_case", {}) or {})
        if case_ctx:
            row = {
                **case_ctx,
                "kind": "polish",
                "stage": stage,
                "retry": "",
                "policy": "",
                "idx": len(mod._096_polish_rows),
                "ok": bool(polished is not None and math.isfinite(residual) and residual < 1e-7),
                "has_root": polished is not None,
                "t": 1.0,
                "residual_scaled": float(residual),
                "residual": float(residual),
                "steps": 0,
                "epochs": int(ctx.get("corrector_epochs_total", 0)),
                "corrector_calls_total": int(ctx.get("corrector_calls_total", 0)),
                "corrector_epochs_total": int(ctx.get("corrector_epochs_total", 0)),
                "seconds": float(time.time() - p0),
                "status": "ok" if polished is not None and math.isfinite(residual) and residual < 1e-7 else "no-root",
                "z": mod.root_to_json(polished) if polished is not None else None,
                "path_json": "",
            }
            mod._096_polish_rows.append(row)
        return polished

    mod.m.polish_070 = polish_070_audit

    original_run_case = mod.run_case

    def run_case_audit(args, case: str):
        n, d = mod.parse_case(case)
        seed = mod.m.seed_for(args.family, n, d, args.seed_index)
        prev = getattr(mod, "_096_current_case", {})
        mod._096_current_case = {"family": args.family, "n": n, "d": d, "seed": seed}
        try:
            return original_run_case(args, case)
        finally:
            mod._096_current_case = prev

    mod.run_case = run_case_audit

    original_launch_batch = mod.launch_batch

    def launch_batch_audit(*args, **kwargs):
        br, rows = original_launch_batch(*args, **kwargs)
        mod._096_captured_batches.append(br)
        return br, rows

    mod.launch_batch = launch_batch_audit

    original_write_outputs = mod.write_outputs

    def write_outputs_audit(summaries, scores, batches, roots_by_case, args):
        old_fast = os.environ.get("PANDROSION_NO_FAST_EXIT")
        os.environ["PANDROSION_NO_FAST_EXIT"] = "1"
        try:
            original_write_outputs(summaries, scores, batches, roots_by_case, args)
            write_step_audit(mod, summaries, batches, args)
        finally:
            if old_fast is None:
                os.environ.pop("PANDROSION_NO_FAST_EXIT", None)
            else:
                os.environ["PANDROSION_NO_FAST_EXIT"] = old_fast

    mod.write_outputs = write_outputs_audit


def main() -> None:
    sys.argv = _consume_096_args(sys.argv)
    f092 = _load_092()
    sys.argv = f092._consume_092_args(sys.argv)
    sys.argv = f092._apply_best_defaults(sys.argv)
    f091 = f092._load_091()
    sys.argv = f091._consume_091_args(sys.argv)
    f090 = f091._load_090()
    sys.argv = f090._consume_090_args(sys.argv)
    f089 = f090._load_089()
    sys.argv = f089._consume_089_args(sys.argv)
    w087 = f089._load_087()
    sys.argv = w087._consume_branch_args(sys.argv)
    mod = w087.load_081()
    w087.install_branch_safe_tracker(mod)
    w087.install_compat_and_single_policy(mod)
    mod.__file__ = str(Path(__file__).resolve())
    f089.install_resultant_recovery(mod)
    f090.install_multivariate_completion(mod)
    f091.install_projective_completion(mod)
    f092.install_pair_recovery(mod)
    install_audit_hooks(mod)
    sys.argv = f090._rewrite_argv_for_090(w087, sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
