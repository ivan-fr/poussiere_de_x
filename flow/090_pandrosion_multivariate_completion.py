"""
FLOW 090 -- multivariate Pandrosion completion sweeps.

090 keeps the branch-safe Pandrosion tracker from 087 and the bivariate
coordinate-resultant recovery from 089, then adds a dedicated multivariate
completion phase for n >= 3.

The multivariate phase is intentionally Pandrosion-only: it does not call the
Lairez/Newton baseline and does not use a Newton-ELS fallback.  Missing sheets
are recovered by deterministic gamma sweeps over the same total-degree starts,
run in parallel waves and accepted only after Pandrosion polish on the original
system.
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Sequence

Complex = complex
Vector = list[Complex]

HERE = Path(__file__).resolve().parent
FLOW089_PATH = HERE / "089_pandrosion_resultant_recovery.py"


def _load_089():
    spec = importlib.util.spec_from_file_location("flow089_for_090", str(FLOW089_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW089_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow089_for_090"] = mod
    spec.loader.exec_module(mod)
    return mod


def _consume_090_args(argv: list[str]) -> list[str]:
    out = [argv[0]]
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--sweep-retries" and i + 1 < len(argv):
            os.environ["PANDROSION_090_SWEEP_RETRIES"] = argv[i + 1]
            i += 2
            continue
        if arg == "--sweep-limit" and i + 1 < len(argv):
            os.environ["PANDROSION_090_SWEEP_LIMIT"] = argv[i + 1]
            i += 2
            continue
        if arg == "--sweep-batch" and i + 1 < len(argv):
            os.environ["PANDROSION_090_SWEEP_BATCH"] = argv[i + 1]
            i += 2
            continue
        if arg == "--no-multivar-sweep":
            os.environ["PANDROSION_090_SWEEP_OFF"] = "1"
            i += 1
            continue
        out.append(arg)
        i += 1
    return out


def _env_int(name: str, default: int) -> int:
    try:
        return int(float(os.environ.get(name, str(default))))
    except Exception:
        return default


def _finite_vec(z: Sequence[Complex]) -> bool:
    try:
        return all(math.isfinite(complex(v).real) and math.isfinite(complex(v).imag) for v in z)
    except Exception:
        return False


def _unique_append(roots: list[Vector], z: Sequence[Complex] | None, sep: float) -> bool:
    if z is None or not _finite_vec(z):
        return False
    n = len(z)
    for root in roots:
        scale = max(1.0, max(abs(z[i]) for i in range(n)), max(abs(root[i]) for i in range(n)))
        if max(abs(z[i] - root[i]) for i in range(n)) <= sep * scale:
            return False
    roots.append(list(z))
    return True


def _valid_cluster(mod, target, roots: Sequence[Sequence[Complex]], args) -> list[Vector]:
    out: list[Vector] = []
    sep = float(getattr(args, "cluster_sep", 1e-6))
    accept = float(getattr(args, "residual_accept", 1e-7))
    for z in roots:
        if z is None or not _finite_vec(z):
            continue
        res = mod.m.residual_norm(target, z)
        if math.isfinite(res) and res < accept:
            _unique_append(out, z, sep)
    return out


def _valid_roots_from_rows(mod, target, rows: Sequence[dict], args) -> list[Vector]:
    roots: list[Vector] = []
    for row in rows:
        z = mod.root_from_json(row.get("z")) if row.get("z") is not None else None
        roots.extend(_valid_cluster(mod, target, [z], args))
    return roots


def _chunks(items: Sequence[int], size: int) -> list[list[int]]:
    size = max(1, int(size))
    return [list(items[i:i + size]) for i in range(0, len(items), size)]


def _select_policy(mod, target, args, selected_name: str):
    policies, _note = mod.build_policies(target, args)
    for policy in policies:
        if policy.name == selected_name:
            return policy
    return policies[0] if policies else None


def _multivariate_completion(mod, target, roots: list[Vector], args, selected_policy: str):
    if os.environ.get("PANDROSION_090_SWEEP_OFF") == "1" or len(target) < 3:
        return 0, 0, 0, 0, 0.0, roots

    B = mod.m.bezout(target)
    roots = _valid_cluster(mod, target, roots, args)
    if len(roots) >= B:
        return 0, 0, 0, 0, 0.0, roots

    sweep_retries = max(0, _env_int("PANDROSION_090_SWEEP_RETRIES", 2))
    if sweep_retries <= 0:
        return 0, 0, 0, 0, 0.0, roots

    policy = _select_policy(mod, target, args, selected_policy)
    if policy is None:
        return 0, 0, 0, 0, 0.0, roots

    limit = max(1, min(B, _env_int("PANDROSION_090_SWEEP_LIMIT", B)))
    batch = max(1, _env_int("PANDROSION_090_SWEEP_BATCH", max(4, int(getattr(args, "micro_batch", 4)))))
    parallel = max(1, int(getattr(args, "parallel_batches", 1)))
    order = mod.golden_order(B)[:limit]
    outdir = Path(args.outdir) / f"{args.family}_{args.case.replace(',', 'x')}_seed{args.seed_index}" / "090_multivar"
    attempts = hits = batches = 0
    retry_used = 0
    t0 = time.time()

    print(
        f"090-multivar-sweep policy={policy.name} retries={sweep_retries} "
        f"limit={limit} batch={batch} parallel={parallel}",
        flush=True,
    )

    for retry in range(1, sweep_retries + 1):
        retry_chunks = _chunks(order, batch)
        for pos in range(0, len(retry_chunks), parallel):
            retry_used = max(retry_used, retry + 1)
            wave = retry_chunks[pos:pos + parallel]
            brs, rows = mod.run_batches_parallel(args, policy, retry, wave, outdir, f"mv{retry}")
            batches += len(brs)
            attempts += sum(int(b.paths) for b in brs)
            cand_roots = _valid_roots_from_rows(mod, target, rows, args)
            before = len(roots)
            roots = _valid_cluster(mod, target, [*roots, *cand_roots], args)
            hits += max(0, len(roots) - before)
            print(
                f"090 retry={retry} wave={pos // parallel + 1} "
                f"paths={sum(int(b.paths) for b in brs)} roots={len(roots)}/{B} "
                f"new={len(roots) - before}",
                flush=True,
            )
            if getattr(args, "stop_at_bezout", False) and len(roots) >= B:
                return attempts, hits, batches, retry_used, time.time() - t0, roots
        if getattr(args, "stop_at_bezout", False) and len(roots) >= B:
            break
    return attempts, hits, batches, retry_used, time.time() - t0, roots


def install_multivariate_completion(mod) -> None:
    original_run_case = mod.run_case

    def run_case_090(args, case: str):
        summaries, scores, batches, roots = original_run_case(args, case)
        n, d = mod.parse_case(case)
        seed = mod.m.seed_for(args.family, n, d, args.seed_index)
        target = mod.m.gen_system(args.family, n, d, seed)
        B = mod.m.bezout(target)
        primary = summaries[0] if summaries else None
        selected = primary.selected_policy if primary is not None else ""
        attempts, hits, extra_batches, retry_used, sec, roots2 = (
            _multivariate_completion(mod, target, roots, args, selected)
            if len(roots) < B else (0, 0, 0, 0, 0.0, roots)
        )
        roots[:] = roots2

        if primary is not None:
            residuals = [mod.m.residual_norm(target, z) for z in roots]
            maxres = max(residuals) if residuals else float("inf")
            status = "ok" if len(roots) >= B and maxres < args.residual_accept else "partial"
            summaries[0] = replace(
                primary,
                alg="090-pandrosion-multivariate-completion",
                roots=len(roots),
                coverage=len(roots) / max(1, B),
                path_rows=primary.path_rows + attempts,
                candidates=primary.candidates + hits,
                batches=primary.batches + extra_batches,
                retries_used=max(primary.retries_used, retry_used),
                max_residual=maxres,
                seconds_observed=primary.seconds_observed + sec,
                status=status,
                notes=(
                    f"{primary.notes}; 090_multivar_attempts={attempts}; "
                    f"090_multivar_hits={hits}; 090_multivar_batches={extra_batches}; "
                    f"090_multivar_seconds={sec:.2f}; Pandrosion-only gamma completion"
                ),
            )
            print(
                f"    090-pandrosion-multivariate-completion roots={len(roots)}/{B} "
                f"cov={100.0 * len(roots) / max(1, B):.1f}% "
                f"paths={primary.path_rows + attempts} maxres={maxres:.2e} "
                f"sec={primary.seconds_observed + sec:.2f} status={status}",
                flush=True,
            )
        return summaries, scores, batches, roots

    mod.run_case = run_case_090


def _rewrite_argv_for_090(w087, argv: list[str]) -> list[str]:
    out = list(argv)
    if "--batch-worker" not in out and "--outdir" not in out:
        out.extend(["--outdir", "090_batches"])
    if "--batch-worker" not in out and "--force-policy" not in out:
        out.extend(["--force-policy", "polygonS-b0.5-hsoft"])
    return w087.rewrite_argv(out)


def main() -> None:
    sys.argv = _consume_090_args(sys.argv)
    f089 = _load_089()
    sys.argv = f089._consume_089_args(sys.argv)
    w087 = f089._load_087()
    sys.argv = w087._consume_branch_args(sys.argv)
    mod = w087.load_081()
    w087.install_branch_safe_tracker(mod)
    w087.install_compat_and_single_policy(mod)
    f089.install_resultant_recovery(mod)
    install_multivariate_completion(mod)
    sys.argv = _rewrite_argv_for_090(w087, sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
