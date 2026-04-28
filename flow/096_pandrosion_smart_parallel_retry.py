"""
FLOW 096 -- smart parallel Pandrosion retry.

096 keeps 092 as the stable base and changes only the recovery schedule around
the batch tracker.  The first pass is the usual 092/081 run.  If roots are
missing, 096 reads the worker JSON rows, detects suspect sheets, then launches
parallel retry waves under a new gamma:

  * failed paths first,
  * duplicated/collided branches next,
  * then a golden-order gamma sweep until the smart retry limit.

The later 089/090/091/092 recovery layers still run afterwards.  Projective
Riemann recovery is enabled for bivariate systems by default, but it is reached
only if the smart retry has not already closed the Bezout count.
"""
from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Sequence

HERE = Path(__file__).resolve().parent
FLOW092_PATH = HERE / "092_pandrosion_resultant_pairing.py"


def _load_092():
    spec = importlib.util.spec_from_file_location("flow092_for_096", str(FLOW092_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW092_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow092_for_096"] = mod
    spec.loader.exec_module(mod)
    return mod


def _has_option(argv: Sequence[str], opt: str) -> bool:
    prefix = opt + "="
    return any(arg == opt or arg.startswith(prefix) for arg in argv[1:])


def _consume_096_args(argv: list[str]) -> list[str]:
    out = [argv[0]]
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--no-smart-retry":
            os.environ["PANDROSION_096_SMART_OFF"] = "1"
            i += 1
            continue
        if arg == "--smart-retry-limit" and i + 1 < len(argv):
            os.environ["PANDROSION_096_SMART_LIMIT"] = argv[i + 1]
            i += 2
            continue
        if arg == "--smart-retry-batch" and i + 1 < len(argv):
            os.environ["PANDROSION_096_SMART_BATCH"] = argv[i + 1]
            i += 2
            continue
        if arg == "--smart-retry-rounds" and i + 1 < len(argv):
            os.environ["PANDROSION_096_SMART_ROUNDS"] = argv[i + 1]
            i += 2
            continue
        if arg == "--smart-retry-budget" and i + 1 < len(argv):
            os.environ["PANDROSION_096_SMART_BUDGET"] = argv[i + 1]
            i += 2
            continue
        if arg == "--smart-residual-cap" and i + 1 < len(argv):
            os.environ["PANDROSION_096_RESIDUAL_CAP"] = argv[i + 1]
            i += 2
            continue
        out.append(arg)
        i += 1
    return out


def _apply_096_defaults(argv: list[str]) -> list[str]:
    out = list(argv)
    defaults = {
        "--parallel-batches": "8",
        "--batch-size": "16",
        "--micro-batch": "8",
        "--retries": "1",
        "--batch-timeout": "20",
        "--sentinel-count": "1",
    }
    for opt, value in defaults.items():
        if not _has_option(out, opt):
            out.extend([opt, value])
    if not _has_option(out, "--projective-bivariate") and not _has_option(out, "--no-projective"):
        out.append("--projective-bivariate")
    if not _has_option(out, "--projective-budget") and not _has_option(out, "--no-projective"):
        out.extend(["--projective-budget", "4"])
    if not _has_option(out, "--projective-limit") and not _has_option(out, "--no-projective"):
        out.extend(["--projective-limit", "32"])
    if not _has_option(out, "--projective-retries") and not _has_option(out, "--no-projective"):
        out.extend(["--projective-retries", "1"])
    return out


def _env_int(name: str, default: int) -> int:
    try:
        return int(float(os.environ.get(name, str(default))))
    except Exception:
        return default


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def _root_from_json(zjson):
    if zjson is None:
        return None
    return [complex(float(a), float(b)) for a, b in zjson]


def _finite_vec(z) -> bool:
    try:
        return all(math.isfinite(complex(v).real) and math.isfinite(complex(v).imag) for v in z)
    except Exception:
        return False


def _unique_append(roots: list[list[complex]], z, sep: float) -> bool:
    if z is None or not _finite_vec(z):
        return False
    n = len(z)
    for root in roots:
        scale = max(1.0, max(abs(z[i]) for i in range(n)), max(abs(root[i]) for i in range(n)))
        if max(abs(z[i] - root[i]) for i in range(n)) <= sep * scale:
            return False
    roots.append(list(z))
    return True


def _cluster_roots(rows_or_roots: Sequence, sep: float, residual_cap: float | None = None) -> list[list[complex]]:
    roots: list[list[complex]] = []
    for item in rows_or_roots:
        if isinstance(item, dict):
            if residual_cap is not None:
                residual = float(item.get("residual", float("inf")) or float("inf"))
                if not math.isfinite(residual) or residual >= residual_cap:
                    continue
            z = _root_from_json(item.get("z")) if item.get("z") is not None else None
        else:
            z = item
        _unique_append(roots, z, sep)
    return roots


def _adaptive_cluster_roots(rows: Sequence[dict], sep: float, B: int,
                            preferred_cap: float, residual_accept: float) -> tuple[list[list[complex]], float]:
    caps: list[float] = []
    for cap in (preferred_cap, 2e-8, 3e-8, 5e-8, residual_accept):
        cap = max(1e-12, min(float(residual_accept), float(cap)))
        if all(abs(cap - old) > 1e-15 for old in caps):
            caps.append(cap)

    seps: list[float] = []
    for factor in (1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0):
        cur = float(sep) * factor
        if all(abs(cur - old) > 1e-18 for old in seps):
            seps.append(cur)

    best_under: tuple[list[list[complex]], float] | None = None
    best_over: tuple[list[list[complex]], float] | None = None
    for cap in caps:
        for sep_try in seps:
            roots = _cluster_roots(rows, sep_try, cap)
            if len(roots) == B:
                return roots, cap
            if len(roots) < B and (best_under is None or len(roots) > len(best_under[0])):
                best_under = (roots, cap)
            if len(roots) > B and (best_over is None or len(roots) < len(best_over[0])):
                best_over = (roots, cap)
    return best_under or best_over or ([], preferred_cap)


def _read_batch_rows(batches) -> list[dict]:
    rows: list[dict] = []
    for b in batches:
        path = getattr(b, "path_json", "")
        if not path:
            continue
        try:
            rows.extend(json.loads(Path(path).read_text()).get("rows", []))
        except Exception:
            pass
    return rows


def _close(z, w, sep: float) -> bool:
    n = len(z)
    scale = max(1.0, max(abs(z[i]) for i in range(n)), max(abs(w[i]) for i in range(n)))
    return max(abs(z[i] - w[i]) for i in range(n)) <= sep * scale


def _suspect_indices(rows: Sequence[dict], sep: float, residual_accept: float) -> list[int]:
    failed: list[int] = []
    clusters: list[list[tuple[list[complex], int]]] = []
    for row in sorted(rows, key=lambda r: int(r.get("idx", -1))):
        idx = int(row.get("idx", -1))
        z = _root_from_json(row.get("z")) if row.get("z") is not None else None
        residual = float(row.get("residual", float("inf")) or float("inf"))
        if idx < 0 or z is None or not _finite_vec(z) or not math.isfinite(residual) or residual >= residual_accept:
            if idx >= 0:
                failed.append(idx)
            continue
        for cluster in clusters:
            if _close(z, cluster[0][0], sep):
                cluster.append((z, idx))
                break
        else:
            clusters.append([(z, idx)])

    dup: list[int] = []
    for cluster in clusters:
        if len(cluster) > 1:
            dup.extend(idx for _z, idx in cluster)

    seen = set()
    out: list[int] = []
    for idx in [*failed, *dup]:
        if idx not in seen:
            seen.add(idx)
            out.append(idx)
    return out


def _chunks(items: Sequence[int], size: int) -> list[list[int]]:
    size = max(1, int(size))
    return [list(items[i:i + size]) for i in range(0, len(items), size)]


def _select_policy(mod, target, args, selected_name: str):
    policies, _note = mod.build_policies(target, args)
    for policy in policies:
        if policy.name == selected_name:
            return policy
    return policies[0] if policies else None


def install_smart_parallel_retry(mod) -> None:
    original_run_case = mod.run_case

    def run_case_096(args, case: str):
        summaries, scores, batches, roots = original_run_case(args, case)
        if os.environ.get("PANDROSION_096_SMART_OFF") == "1":
            return summaries, scores, batches, roots

        n, d = mod.parse_case(case)
        seed = mod.m.seed_for(args.family, n, d, args.seed_index)
        target = mod.m.gen_system(args.family, n, d, seed)
        B = mod.m.bezout(target)
        primary = summaries[0] if summaries else None
        if primary is None or len(roots) >= B:
            return summaries, scores, batches, roots

        sep = float(getattr(args, "cluster_sep", 1e-6))
        tight_sep = min(sep, 1e-10) if n == 2 else sep
        residual_cap = min(
            float(getattr(args, "residual_accept", 1e-7)),
            max(1e-12, _env_float("PANDROSION_096_RESIDUAL_CAP", 2e-8)),
        )
        rows0 = _read_batch_rows(batches)
        suspects = _suspect_indices(rows0, tight_sep, residual_cap)
        golden = mod.golden_order(B)
        limit = max(0, min(B, _env_int("PANDROSION_096_SMART_LIMIT", min(B, 112))))
        retry_order = list(dict.fromkeys([*suspects, *golden]))[:limit]
        if not retry_order:
            return summaries, scores, batches, roots

        selected = _select_policy(mod, target, args, getattr(primary, "selected_policy", ""))
        if selected is None:
            return summaries, scores, batches, roots

        batch = max(1, _env_int("PANDROSION_096_SMART_BATCH", max(8, int(getattr(args, "micro_batch", 8)))))
        parallel = max(1, int(getattr(args, "parallel_batches", 1)))
        rounds = max(1, _env_int("PANDROSION_096_SMART_ROUNDS", 1))
        budget = max(0.0, _env_float("PANDROSION_096_SMART_BUDGET", 0.0))
        deadline = time.time() + budget if budget > 0.0 else 1e100
        outdir = Path(args.outdir) / f"{args.family}_{n}x{d}_seed{args.seed_index}" / "096_smart"

        all_rows = list(rows0)
        roots, effective_cap = _adaptive_cluster_roots(
            all_rows,
            tight_sep,
            B,
            residual_cap,
            float(getattr(args, "residual_accept", 1e-7)),
        )
        attempts = hits = 0
        extra_batches = 0
        retry_used = int(getattr(primary, "retries_used", 1))
        t0 = time.time()
        print(
            f"096-smart-retry suspects={len(suspects)} limit={len(retry_order)} "
            f"batch={batch} parallel={parallel} rounds={rounds} residual_cap={residual_cap:g}",
            flush=True,
        )

        for retry in range(1, rounds + 1):
            if len(roots) >= B or time.time() >= deadline:
                break
            chunks = _chunks(retry_order, batch)
            for pos in range(0, len(chunks), parallel):
                if len(roots) >= B or time.time() >= deadline:
                    break
                wave = chunks[pos:pos + parallel]
                brs, rows = mod.run_batches_parallel(args, selected, retry, wave, outdir, f"smart{retry}")
                batches.extend(brs)
                extra_batches += len(brs)
                attempts += sum(int(getattr(b, "paths", 0)) for b in brs)
                before = len(roots)
                all_rows.extend(rows)
                roots, effective_cap = _adaptive_cluster_roots(
                    all_rows,
                    tight_sep,
                    B,
                    residual_cap,
                    float(getattr(args, "residual_accept", 1e-7)),
                )
                gained = max(0, len(roots) - before)
                hits += gained
                retry_used = max(retry_used, retry + 1)
                print(
                    f"096 retry={retry} wave={pos // parallel + 1} "
                    f"paths={sum(int(getattr(b, 'paths', 0)) for b in brs)} "
                    f"roots={len(roots)}/{B} new={gained} cap={effective_cap:g}",
                    flush=True,
                )
                if getattr(args, "stop_at_bezout", False) and len(roots) >= B:
                    break

        residuals = [mod.m.residual_norm(target, z) for z in roots]
        maxres = max(residuals) if residuals else float("inf")
        status = "ok" if len(roots) >= B and maxres < args.residual_accept else "partial"
        sec = time.time() - t0
        summaries[0] = replace(
            primary,
            alg="096-pandrosion-smart-parallel-retry",
            roots=len(roots),
            coverage=len(roots) / max(1, B),
            path_rows=primary.path_rows + attempts,
            candidates=primary.candidates + hits,
            batches=primary.batches + extra_batches,
            retries_used=retry_used,
            max_residual=maxres,
            seconds_observed=primary.seconds_observed + sec,
            status=status,
            notes=(
                f"{primary.notes}; 096_smart_attempts={attempts}; "
                f"096_smart_hits={hits}; 096_smart_suspects={len(suspects)}; "
                f"096_smart_seconds={sec:.2f}; smart parallel gamma retry"
            ),
        )
        roots[:] = roots
        print(
            f"    096-pandrosion-smart-parallel-retry roots={len(roots)}/{B} "
            f"cov={100.0 * len(roots) / max(1, B):.1f}% paths={primary.path_rows + attempts} "
            f"maxres={maxres:.2e} sec={primary.seconds_observed + sec:.2f} status={status}",
            flush=True,
        )
        return summaries, scores, batches, roots

    mod.run_case = run_case_096


def main() -> None:
    sys.argv = _consume_096_args(sys.argv)
    sys.argv = _apply_096_defaults(sys.argv)
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
    install_smart_parallel_retry(mod)
    f089.install_resultant_recovery(mod)
    f090.install_multivariate_completion(mod)
    f091.install_projective_completion(mod)
    f092.install_pair_recovery(mod)
    sys.argv = f090._rewrite_argv_for_090(w087, sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
