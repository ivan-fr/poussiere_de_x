"""
FLOW 088 -- branch-safe Pandrosion with deterministic shell recovery.

088 keeps the branch-safe residual-gated tracker from 087 and adds a small,
upfront deterministic recovery phase on the target system itself.  The recovery
starts Pandrosion polishing from logarithmic coordinate shells motivated by the
Jensen/resultant shell used in 081.  This fills roots lost to branch collisions
without switching to Lairez/Newton-ELS and without adding unknowns.
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

HERE = Path(__file__).resolve().parent
FLOW087_PATH = HERE / "087_pandrosion_branch_safe_tracker.py"


def _load_087():
    spec = importlib.util.spec_from_file_location("flow087_for_088", str(FLOW087_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW087_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow087_for_088"] = mod
    spec.loader.exec_module(mod)
    return mod


def _consume_088_args(argv: list[str]) -> list[str]:
    out = [argv[0]]
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--direct-shell-limit" and i + 1 < len(argv):
            os.environ["PANDROSION_088_DIRECT_LIMIT"] = argv[i + 1]
            i += 2
            continue
        if arg == "--direct-shell-timeout" and i + 1 < len(argv):
            os.environ["PANDROSION_088_DIRECT_TIMEOUT"] = argv[i + 1]
            i += 2
            continue
        if arg == "--direct-shell-radii" and i + 1 < len(argv):
            os.environ["PANDROSION_088_DIRECT_RADII"] = argv[i + 1]
            i += 2
            continue
        if arg == "--no-direct-shell":
            os.environ["PANDROSION_088_DIRECT_LIMIT"] = "0"
            i += 1
            continue
        out.append(arg)
        i += 1
    return out


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def _env_int(name: str, default: int) -> int:
    try:
        return int(float(os.environ.get(name, str(default))))
    except Exception:
        return default


def _parse_radii() -> list[float]:
    raw = os.environ.get("PANDROSION_088_DIRECT_RADII", "")
    if raw:
        vals = []
        for part in raw.replace(";", ",").split(","):
            part = part.strip()
            if part:
                vals.append(max(1e-6, float(part)))
        if vals:
            return vals
    return [0.18, 0.28, 0.42, 0.65, 1.0, 1.55, 2.4, 3.7, 5.7, 8.8]


def _unique_append(m, roots: list[list[complex]], z: Sequence[complex], sep: float) -> bool:
    if z is None or not m.safe_z(z):
        return False
    for r in roots:
        scale = max(1.0, max(abs(v) for v in z), max(abs(v) for v in r))
        if max(abs(z[i] - r[i]) for i in range(len(z))) <= sep * scale:
            return False
    roots.append(list(z))
    return True


def _direct_shell_recovery(mod, target, roots: list[list[complex]], args, deadline: float | None) -> tuple[int, int, float]:
    m = mod.m
    f076 = mod.f076
    n = len(target)
    limit = _env_int("PANDROSION_088_DIRECT_LIMIT", 400)
    if limit <= 0 or n != 2:
        return 0, 0, 0.0

    radii = _parse_radii()
    attempts = 0
    hits = 0
    t0 = time.time()
    sep = float(getattr(args, "cluster_sep", 1e-6))
    timeout = _env_float("PANDROSION_088_DIRECT_TIMEOUT", 6.0)
    local_deadline = t0 + timeout if timeout > 0 else None

    for ia, r1 in enumerate(radii):
        for ib, r2 in enumerate(radii):
            for k in range(4):
                if attempts >= limit:
                    return attempts, hits, time.time() - t0
                if deadline is not None and time.time() >= deadline:
                    return attempts, hits, time.time() - t0
                if local_deadline is not None and time.time() >= local_deadline:
                    return attempts, hits, time.time() - t0
                th1 = 2.0 * math.pi * ((k * 0.38196601125 + ia * 0.137 + ib * 0.071) % 1.0)
                th2 = 2.0 * math.pi * ((k * 0.61803398875 + ia * 0.053 + ib * 0.193) % 1.0)
                z0 = [
                    r1 * complex(math.cos(th1), math.sin(th1)),
                    r2 * complex(math.cos(th2), math.sin(th2)),
                ]
                attempts += 1
                zp = m.polish_070(
                    target,
                    z0,
                    args.tol,
                    max(16, int(args.quad_cap)),
                    f076.DEFAULT_MODES,
                    f076.DEFAULT_RESCUE_MODES,
                    None,
                )
                if zp is not None and m.residual_norm(target, zp) < args.residual_accept:
                    if _unique_append(m, roots, zp, sep):
                        hits += 1
                        if len(roots) >= m.bezout(target):
                            return attempts, hits, time.time() - t0
    return attempts, hits, time.time() - t0


def install_shell_recovery(mod) -> None:
    original_run_case = mod.run_case

    def run_case_088(args, case: str):
        summaries, scores, batches, roots = original_run_case(args, case)
        n, d = mod.parse_case(case)
        seed = mod.m.seed_for(args.family, n, d, args.seed_index)
        target = mod.m.gen_system(args.family, n, d, seed)
        B = mod.m.bezout(target)
        if len(roots) < B:
            attempts, hits, sec = _direct_shell_recovery(mod, target, roots, args, None)
        else:
            attempts, hits, sec = 0, 0, 0.0

        if summaries:
            primary = summaries[0]
            residuals = [mod.m.residual_norm(target, z) for z in roots]
            maxres = max(residuals) if residuals else float("inf")
            status = "ok" if len(roots) >= B and maxres < args.residual_accept else "partial"
            summaries[0] = replace(
                primary,
                alg="088-branch-safe-shell-recovery",
                roots=len(roots),
                coverage=len(roots) / max(1, B),
                path_rows=primary.path_rows + attempts,
                candidates=primary.candidates + hits,
                max_residual=maxres,
                seconds_observed=primary.seconds_observed + sec,
                status=status,
                notes=(
                    f"{primary.notes}; direct_shell_attempts={attempts}; "
                    f"direct_shell_hits={hits}; direct_shell_seconds={sec:.2f}; "
                    "Pandrosion-only shell recovery"
                ),
            )
            print(
                f"direct-shell attempts={attempts} hits={hits} roots={len(roots)}/{B} "
                f"sec={sec:.2f}",
                flush=True,
            )
            print(
                f"    088-branch-safe-shell-recovery roots={len(roots)}/{B} "
                f"cov={100.0 * len(roots) / max(1, B):.1f}% "
                f"paths={primary.path_rows + attempts} maxres={maxres:.2e} "
                f"sec={primary.seconds_observed + sec:.2f} status={status}",
                flush=True,
            )
        return summaries, scores, batches, roots

    mod.run_case = run_case_088


def _rewrite_argv_for_088(w087, argv: list[str]) -> list[str]:
    out = list(argv)
    if "--batch-worker" not in out and "--outdir" not in out:
        out.extend(["--outdir", "088_batches"])
    return w087.rewrite_argv(out)


def main() -> None:
    sys.argv = _consume_088_args(sys.argv)
    w087 = _load_087()
    sys.argv = w087._consume_branch_args(sys.argv)
    mod = w087.load_081()
    w087.install_branch_safe_tracker(mod)
    w087.install_compat_and_single_policy(mod)
    install_shell_recovery(mod)
    sys.argv = _rewrite_argv_for_088(w087, sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
