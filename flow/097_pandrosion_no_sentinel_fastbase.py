"""
FLOW 097 -- no-sentinel, fail-fast base tracking for the 092 solver.

This keeps the 092 Pandrosion/resultant-pair stack intact, but removes the
sentinel policy probe when the policy is already forced by the 087/090 defaults.
It also stops late base paths that are stuck in repeated corrective failures.

The intent is not to add a fallback: expensive near-t=1 misses are left for the
existing Pandrosion resultant/pair recovery instead of burning hundreds of extra
corrector epochs on a path that is no longer making useful progress.
"""
from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Sequence

HERE = Path(__file__).resolve().parent
FLOW092_PATH = HERE / "092_pandrosion_resultant_pairing.py"


def _load_092():
    spec = importlib.util.spec_from_file_location("flow092_for_097", str(FLOW092_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW092_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow092_for_097"] = mod
    spec.loader.exec_module(mod)
    return mod


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


def _has_option(argv: Sequence[str], opt: str) -> bool:
    prefix = opt + "="
    return any(arg == opt or arg.startswith(prefix) for arg in argv[1:])


def _consume_097_args(argv: list[str]) -> list[str]:
    out = [argv[0]]
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg in {"--keep-sentinel", "--with-sentinel"}:
            os.environ["PANDROSION_097_SKIP_SENTINEL"] = "0"
            i += 1
            continue
        if arg in {"--skip-sentinel", "--no-sentinel"}:
            os.environ["PANDROSION_097_SKIP_SENTINEL"] = "1"
            i += 1
            continue
        if arg == "--no-base-failfast":
            os.environ["PANDROSION_097_FAILFAST_OFF"] = "1"
            i += 1
            continue
        if arg == "--failfast-t" and i + 1 < len(argv):
            os.environ["PANDROSION_097_FAILFAST_T"] = argv[i + 1]
            i += 2
            continue
        if arg == "--failfast-fails" and i + 1 < len(argv):
            os.environ["PANDROSION_097_FAILFAST_FAILS"] = argv[i + 1]
            i += 2
            continue
        if arg == "--failfast-epochs" and i + 1 < len(argv):
            os.environ["PANDROSION_097_FAILFAST_EPOCHS"] = argv[i + 1]
            i += 2
            continue
        if arg == "--failfast-hard-epochs" and i + 1 < len(argv):
            os.environ["PANDROSION_097_FAILFAST_HARD_EPOCHS"] = argv[i + 1]
            i += 2
            continue
        if arg == "--failfast-min-steps" and i + 1 < len(argv):
            os.environ["PANDROSION_097_FAILFAST_MIN_STEPS"] = argv[i + 1]
            i += 2
            continue
        out.append(arg)
        i += 1
    return out


def _apply_097_defaults(argv: list[str]) -> list[str]:
    out = list(argv)
    if "--batch-worker" not in out and not _has_option(out, "--outdir"):
        out.extend(["--outdir", "097_batches"])
    os.environ.setdefault("PANDROSION_097_SKIP_SENTINEL", "1")
    os.environ.setdefault("PANDROSION_097_FAILFAST_T", "0.94")
    os.environ.setdefault("PANDROSION_097_FAILFAST_FAILS", "12")
    os.environ.setdefault("PANDROSION_097_FAILFAST_EPOCHS", "320")
    os.environ.setdefault("PANDROSION_097_FAILFAST_HARD_EPOCHS", "0")
    os.environ.setdefault("PANDROSION_097_FAILFAST_MIN_STEPS", "48")
    os.environ.setdefault("PANDROSION_097_AGGRESSIVE_MAX_DEGREE", "12")
    return out


def install_no_sentinel(mod) -> None:
    original_launch_batch = mod.launch_batch

    def launch_batch_097(args, policy, retry, indices, outdir, stage):
        if (
            os.environ.get("PANDROSION_097_SKIP_SENTINEL", "1") != "0"
            and stage == "sentinel"
            and getattr(args, "force_policy", "")
        ):
            br = mod.BatchResult080(
                stage,
                retry,
                policy.name,
                mod.encode_indices(indices),
                0,
                0,
                0,
                float("inf"),
                0.0,
                "skipped",
                "",
            )
            return br, []
        return original_launch_batch(args, policy, retry, indices, outdir, stage)

    mod.launch_batch = launch_batch_097

    original_run_case = mod.run_case

    def run_case_no_sentinel(args, case: str):
        if os.environ.get("PANDROSION_097_SKIP_SENTINEL", "1") == "0":
            return original_run_case(args, case)

        args.case = case
        n, d = mod.parse_case(case)
        seed = mod.m.seed_for(args.family, n, d, args.seed_index)
        target = mod.m.gen_system(args.family, n, d, seed)
        B = mod.m.bezout(target)
        terms = mod.m.term_count(target)
        policies, polygon_note = mod.build_policies(target, args)
        outdir = Path(args.outdir) / f"{args.family}_{n}x{d}_seed{args.seed_index}"
        outdir.mkdir(parents=True, exist_ok=True)

        print("=" * 124, flush=True)
        print("097 -- no-sentinel fast-base Pandrosion/resultant pairing", flush=True)
        print("=" * 124, flush=True)
        print(f"family={args.family}, case=({n},{d}), seed={seed}, terms={terms}, Bezout={B}", flush=True)
        print(
            f"policies={len(policies)}, sentinel=off, batch={args.batch_size}, "
            f"parallel_batches={args.parallel_batches}, timeout={args.batch_timeout:g}s, "
            f"eq_norm={args.equation_normalize}; no Newton-ELS",
            flush=True,
        )
        print(f"polygon_note={polygon_note}", flush=True)
        for p in policies:
            print(f"  policy {p.name}: S={mod.encode_floats(p.scales)} start={mod.encode_floats(p.start_radii)}", flush=True)

        t0 = time.time()
        order = mod.golden_order(B)
        selected = next((p for p in policies if getattr(args, "force_policy", "") and p.name == args.force_policy), policies[0])
        print(f"selected_policy={selected.name}: S={mod.encode_floats(selected.scales)} start={mod.encode_floats(selected.start_radii)}", flush=True)
        scores = [
            mod.PolicyScore080(
                selected.name,
                0,
                0,
                0,
                0.0,
                float("inf"),
                0.0,
                "skipped",
                f"{selected.note}; sentinel skipped by 097",
            )
        ]

        all_batches = []
        all_rows: list[dict] = []
        chunks = [list(range(i, min(B, i + args.batch_size))) for i in range(0, B, args.batch_size)]
        base_batches, base_rows = mod.run_batches_parallel(args, selected, 0, chunks, outdir, "base")
        all_batches.extend(base_batches)
        all_rows.extend(base_rows)
        roots = mod.cluster_roots(all_rows, sep=args.cluster_sep)
        candidates = sum(1 for r in all_rows if r.get("z") is not None)
        print(f"base selected={selected.name} candidates={candidates}/{len(all_rows)} roots={len(roots)}/{B}", flush=True)

        retries_used = 1
        if len(roots) < B and args.retries > 1:
            failed = [int(r.get("idx", -1)) for r in base_rows if r.get("z") is None]
            retry_order = list(dict.fromkeys([*failed, *order]))[:max(0, int(args.retry_limit))]
            for retry in range(1, int(args.retries)):
                retries_used = retry + 1
                if args.stop_at_bezout and len(roots) >= B:
                    break
                chunks_r = [retry_order[i:i + args.micro_batch] for i in range(0, len(retry_order), args.micro_batch)]
                for ch in chunks_r:
                    br, rows = mod.launch_batch(args, selected, retry, ch, outdir, f"retry{retry}")
                    all_batches.append(br)
                    all_rows.extend(rows)
                    roots = mod.cluster_roots(all_rows, sep=args.cluster_sep)
                    candidates = sum(1 for r in all_rows if r.get("z") is not None)
                    print(f"retry={retry} policy={selected.name} indices={mod.encode_indices(ch)} roots={len(roots)}/{B} candidates={candidates}", flush=True)
                    if args.stop_at_bezout and len(roots) >= B:
                        break

        roots = mod.cluster_roots(all_rows, sep=args.cluster_sep)
        residuals = [mod.m.residual_norm(target, z) for z in roots]
        maxres = max(residuals) if residuals else float("inf")
        candidates = sum(1 for r in all_rows if r.get("z") is not None)
        status = "ok" if len(roots) >= B and maxres < args.residual_accept else "partial"
        sec = time.time() - t0
        notes = (
            f"097 no-sentinel fast-base; selected={selected.name}; S={mod.encode_floats(selected.scales)}; "
            f"start={mod.encode_floats(selected.start_radii)}; {polygon_note}; sentinel=skipped; "
            f"failfast_t={_env_float('PANDROSION_097_FAILFAST_T', 0.94):.4g}; "
            f"failfast_fails={_env_int('PANDROSION_097_FAILFAST_FAILS', 12)}; "
            f"failfast_epochs={_env_int('PANDROSION_097_FAILFAST_EPOCHS', 320)}; "
            f"aggressive_max_degree={_env_int('PANDROSION_097_AGGRESSIVE_MAX_DEGREE', 12)}; "
            f"guard_scores={[asdict(s) for s in scores]}; no Newton-ELS"
        )
        summaries = [mod.Summary080(args.family, n, d, seed, terms, B,
            "097-no-sentinel-fast-base", selected.name, len(roots), len(roots) / max(1, B),
            len(all_rows), candidates, len(all_batches), retries_used, maxres, sec, status, notes)]

        if args.include_lairez_reference:
            ref = mod.f076.load_lairez_reference(args.family, n, d, seed, B)
            if ref is not None:
                summaries.append(mod.Summary080(ref.family, ref.n, ref.d, ref.seed, ref.terms, ref.bezout,
                    "lairez-style-reference", "n/a", ref.roots, ref.coverage, ref.path_rows,
                    ref.candidates, 0, ref.retries_used, ref.max_residual, ref.seconds_observed,
                    ref.status, ref.notes))
        if args.run_lairez:
            import argparse
            ns = argparse.Namespace(**vars(args))
            ns.case = f"{n},{d}"
            lr = mod.f076.run_lairez_now(ns, target, seed, terms, B)
            summaries.append(mod.Summary080(lr.family, lr.n, lr.d, lr.seed, lr.terms, lr.bezout,
                "lairez-style-run", "n/a", lr.roots, lr.coverage, lr.path_rows,
                lr.candidates, 0, lr.retries_used, lr.max_residual, lr.seconds_observed,
                lr.status, lr.notes))

        for r in summaries:
            print(
                f"{r.alg:>38} roots={r.roots}/{r.bezout} cov={100*r.coverage:.1f}% "
                f"paths={r.path_rows} batches={r.batches} maxres={r.max_residual:.2e} "
                f"sec={r.seconds_observed:.2f} status={r.status}",
                flush=True,
            )
        return summaries, scores, all_batches, roots

    mod.run_case = run_case_no_sentinel


def install_failfast_tracker(mod) -> None:
    if os.environ.get("PANDROSION_097_FAILFAST_OFF") == "1":
        return
    m = mod.m

    def correction_is_ambiguous(z, pred, zc) -> bool:
        corr = max((abs(zc[i] - pred[i]) for i in range(len(z))), default=0.0)
        move = max((abs(pred[i] - z[i]) for i in range(len(z))), default=0.0)
        scale = max(
            1.0,
            max((abs(v) for v in z), default=0.0),
            max((abs(v) for v in pred), default=0.0),
            max((abs(v) for v in zc), default=0.0),
        )
        scale_gate = max(0.02, _env_float("PANDROSION_BRANCH_SCALE", 0.10))
        move_gate = max(0.50, _env_float("PANDROSION_BRANCH_MOVE", 1.25))
        return corr > max(scale_gate * scale, move_gate * max(move, 1e-12))

    def should_stop_late(t: float, steps: int, epochs: int, fails: int, degree: int) -> bool:
        min_steps = max(1, _env_int("PANDROSION_097_FAILFAST_MIN_STEPS", 48))
        late_t = max(0.0, min(0.999999, _env_float("PANDROSION_097_FAILFAST_T", 0.94)))
        fail_cap = max(1, _env_int("PANDROSION_097_FAILFAST_FAILS", 12))
        epoch_cap = max(1, _env_int("PANDROSION_097_FAILFAST_EPOCHS", 320))
        hard_epoch_cap = max(0, _env_int("PANDROSION_097_FAILFAST_HARD_EPOCHS", 0))
        aggressive_max_degree = max(1, _env_int("PANDROSION_097_AGGRESSIVE_MAX_DEGREE", 12))
        if degree > aggressive_max_degree:
            fail_cap = max(fail_cap, 2 * degree)
            epoch_cap = max(epoch_cap, 80 * degree)
        if steps < min_steps:
            return False
        if hard_epoch_cap and epochs >= hard_epoch_cap:
            return True
        return t >= late_t and (fails >= fail_cap or epochs >= epoch_cap)

    def track_one_097(target, start, target_gamma, z0,
                      tol: float = 1e-9, max_steps: int = 420,
                      max_epochs: int = 12, quad_cap: int = 16,
                      modes: Sequence[str] = ("system", "integral_gl", "blend"),
                      rescue_modes: Sequence[str] = (),
                      deadline: float | None = None):
        z = list(z0)
        D = max(m.degrees(target))
        dt = min(0.010, max(0.0015, 0.060 / max(2, D * D)))
        dtmax = min(0.030, max(0.006, 0.18 / max(2, D)))
        t = 0.0
        steps = epochs = fails = guarded = 0
        prev_z = None
        prev_t = None
        status = "ok"
        accept_residual = 80.0 * tol

        while t < 1.0 - 1e-15 and steps < max_steps:
            if m.timed_out(deadline):
                status = "budget"
                break
            if should_stop_late(t, steps, epochs, fails, D):
                status = "early-stop"
                break
            steps += 1
            tnext = min(1.0, t + dt)
            pred = m.tangent_predictor(target_gamma, start, z, t, tnext - t, prev_z=prev_z, prev_t=prev_t)
            Hnext = m.homotopy_polys(target_gamma, start, tnext)
            zc, _ok, ep = m.corrector(Hnext, pred, tol=tol, max_epochs=max_epochs,
                                      quad_cap=quad_cap, modes=modes,
                                      rescue_modes=rescue_modes, deadline=deadline)
            epochs += ep
            rh = m.residual_norm(Hnext, zc)
            ambiguous = correction_is_ambiguous(z, pred, zc)
            if math.isfinite(rh) and rh < accept_residual and not ambiguous:
                prev_z, prev_t = list(z), t
                z, t = list(zc), tnext
                dt = min(dtmax, dt * (1.20 if ep <= max(3, max_epochs // 3) else 1.06))
                fails = max(0, fails - 1)
                continue
            if math.isfinite(rh) and rh < 1e-6 and not ambiguous and not m.timed_out(deadline):
                zr, _okr, epr = m.corrector(Hnext, zc, tol=tol, max_epochs=max_epochs + 6,
                                            quad_cap=quad_cap,
                                            modes=("blend", "system", "integral_gl"),
                                            rescue_modes=rescue_modes,
                                            deadline=deadline)
                epochs += epr
                rr = m.residual_norm(Hnext, zr)
                if math.isfinite(rr) and rr < accept_residual and not correction_is_ambiguous(z, pred, zr):
                    prev_z, prev_t = list(z), t
                    z, t = list(zr), tnext
                    dt = min(dtmax, dt * 1.03)
                    fails = max(0, fails - 1)
                    continue
            if ambiguous:
                guarded += 1
            fails += 1
            dt *= 0.5
            if should_stop_late(t, steps, epochs, fails, D):
                status = "early-stop"
                break
            if dt < 5e-7 or fails > 80:
                status = "guarded" if guarded else "fail"
                break

        if t >= 1.0 - 1e-12 and not m.timed_out(deadline):
            zp, _okp, ep = m.corrector(target, z, tol=tol, max_epochs=max_epochs + 10,
                                       quad_cap=quad_cap,
                                       modes=("system", "integral_gl", "blend"),
                                       rescue_modes=rescue_modes,
                                       deadline=deadline)
            epochs += ep
            rp = m.residual_norm(target, zp)
            if math.isfinite(rp) and rp < 1e-7:
                z = zp
        res = m.residual_norm(target, z)
        ok_final = t >= 1.0 - 1e-12 and res < 1e-7
        return m.PathResult(z=z, ok=ok_final, steps=steps, epochs=epochs, t=t,
                            residual=res, status=status if not ok_final else "ok")

    m.track_one_070 = track_one_097


def main() -> None:
    sys.argv = _consume_097_args(sys.argv)
    sys.argv = _apply_097_defaults(sys.argv)
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
    install_no_sentinel(mod)
    install_failfast_tracker(mod)
    f089.install_resultant_recovery(mod)
    f090.install_multivariate_completion(mod)
    f091.install_projective_completion(mod)
    f092.install_pair_recovery(mod)
    sys.argv = f090._rewrite_argv_for_090(w087, sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
