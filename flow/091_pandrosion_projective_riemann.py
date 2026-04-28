"""
FLOW 091 -- Pandrosion projective/Riemann-sphere completion.

091 builds on 090 and adds a last-resort projective chart portfolio.  Each
coordinate is treated as a Riemann sphere with two affine charts:

    north: z_j = y_j,
    south: z_j = s_j / y_j.

For a south-chart coordinate, denominators are cleared equation-wise so the
tracked system remains polynomial.  The chart may introduce roots at infinity;
therefore every candidate is mapped back to the original z-coordinates and
accepted only if the original residual is small after Pandrosion polish.

No Lairez/Newton-ELS fallback is used.
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Iterable, Sequence

Complex = complex
Vector = list[Complex]
Poly = dict[tuple[int, ...], Complex]
System = list[Poly]

HERE = Path(__file__).resolve().parent
FLOW090_PATH = HERE / "090_pandrosion_multivariate_completion.py"


def _load_090():
    spec = importlib.util.spec_from_file_location("flow090_for_091", str(FLOW090_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW090_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow090_for_091"] = mod
    spec.loader.exec_module(mod)
    return mod


def _consume_091_args(argv: list[str]) -> list[str]:
    out = [argv[0]]
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--projective-retries" and i + 1 < len(argv):
            os.environ["PANDROSION_091_PROJECTIVE_RETRIES"] = argv[i + 1]
            i += 2
            continue
        if arg == "--projective-limit" and i + 1 < len(argv):
            os.environ["PANDROSION_091_PROJECTIVE_LIMIT"] = argv[i + 1]
            i += 2
            continue
        if arg == "--projective-max-charts" and i + 1 < len(argv):
            os.environ["PANDROSION_091_PROJECTIVE_MAX_CHARTS"] = argv[i + 1]
            i += 2
            continue
        if arg == "--projective-budget" and i + 1 < len(argv):
            os.environ["PANDROSION_091_PROJECTIVE_BUDGET"] = argv[i + 1]
            i += 2
            continue
        if arg == "--projective-bivariate":
            os.environ["PANDROSION_091_PROJECTIVE_BIVARIATE"] = "1"
            i += 1
            continue
        if arg == "--projective-pairs":
            os.environ["PANDROSION_091_PROJECTIVE_PAIRS"] = "1"
            i += 1
            continue
        if arg == "--no-projective":
            os.environ["PANDROSION_091_PROJECTIVE_OFF"] = "1"
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


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def _env_bool(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() not in {"0", "false", "no", "off", ""}


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


def _valid_cluster(mod, target: System, roots: Sequence[Sequence[Complex]], args) -> list[Vector]:
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


def _quantile(vals: Sequence[float], q: float, default: float = 1.0) -> float:
    clean = sorted(v for v in vals if math.isfinite(v) and v > 0.0)
    if not clean:
        return default
    idx = int(round(max(0.0, min(1.0, q)) * (len(clean) - 1)))
    return clean[idx]


def _riemann_scales(roots: Sequence[Sequence[Complex]], n: int) -> list[float]:
    out: list[float] = []
    for j in range(n):
        vals = [abs(z[j]) for z in roots if len(z) > j and math.isfinite(abs(z[j]))]
        med = _quantile(vals, 0.50, 1.0)
        q90 = _quantile(vals, 0.90, med)
        s = math.sqrt(max(1e-12, med) * max(1e-12, q90))
        # Keep the south chart near the observed root cloud, without allowing a
        # single extreme root to make the chart numerically wild.
        out.append(max(0.25, min(8.0, s)))
    return out


def _chart_scores(roots: Sequence[Sequence[Complex]], n: int) -> list[tuple[float, int]]:
    scores: list[tuple[float, int]] = []
    for j in range(n):
        vals = [abs(z[j]) for z in roots if len(z) > j and math.isfinite(abs(z[j]))]
        q10 = _quantile(vals, 0.10, 1.0)
        q50 = _quantile(vals, 0.50, 1.0)
        q90 = _quantile(vals, 0.90, q50)
        spread = math.log(max(q90, 1e-12) / max(q10, 1e-12))
        off_unit = abs(math.log(max(q50, 1e-12)))
        scores.append((spread + 0.35 * off_unit, j))
    scores.sort(reverse=True)
    return scores


def _chart_masks(roots: Sequence[Sequence[Complex]], n: int) -> list[tuple[bool, ...]]:
    scored = _chart_scores(roots, n)
    masks: list[tuple[bool, ...]] = []
    for _score, j in scored:
        mask = tuple(k == j for k in range(n))
        masks.append(mask)
    if _env_bool("PANDROSION_091_PROJECTIVE_PAIRS", False):
        top = [j for _score, j in scored[:min(n, 4)]]
        for a_pos, a in enumerate(top):
            for b in top[a_pos + 1:]:
                masks.append(tuple(k == a or k == b for k in range(n)))
    max_charts = max(1, _env_int("PANDROSION_091_PROJECTIVE_MAX_CHARTS", n))
    seen = set()
    out: list[tuple[bool, ...]] = []
    for mask in masks:
        if mask in seen:
            continue
        seen.add(mask)
        out.append(mask)
        if len(out) >= max_charts:
            break
    return out


def _projective_invert_system(polys: System, mask: Sequence[bool], scales: Sequence[float]) -> System:
    out: System = []
    n = len(mask)
    for poly in polys:
        caps = [0 for _ in range(n)]
        for alpha in poly:
            for j, inv in enumerate(mask):
                if inv:
                    caps[j] = max(caps[j], int(alpha[j]))
        q: Poly = {}
        for alpha, coeff in poly.items():
            beta = list(alpha)
            c = complex(coeff)
            for j, inv in enumerate(mask):
                if inv:
                    c *= float(scales[j]) ** int(alpha[j])
                    beta[j] = caps[j] - int(alpha[j])
            e = tuple(int(v) for v in beta)
            q[e] = q.get(e, 0.0 + 0.0j) + c
        out.append({e: c for e, c in q.items() if abs(c) > 1e-14})
    return out


def _from_projective_chart(y: Sequence[Complex], mask: Sequence[bool], scales: Sequence[float]) -> Vector | None:
    z: Vector = []
    for yj, inv, sj in zip(y, mask, scales):
        yj = complex(yj)
        if inv:
            if abs(yj) < 1e-12:
                return None
            z.append(float(sj) / yj)
        else:
            z.append(yj)
    return z if _finite_vec(z) else None


def _mask_key(mask: Sequence[bool]) -> int:
    out = 0
    for i, bit in enumerate(mask):
        if bit:
            out += 1 << i
    return out


def _projective_track_one(mod, target: System, chart_target: System, mask: Sequence[bool],
                          scales: Sequence[float], seed: int, idx: int, retry: int,
                          args, deadline: float | None) -> Vector | None:
    m = mod.m
    degs = m.degrees(chart_target)
    Bc = math.prod(degs)
    if Bc <= 0:
        return None
    idx = int(idx) % Bc
    key = _mask_key(mask)
    phases = m.m068.deterministic_phases(len(target), 0, seed + 17 * Bc + 10007 * retry + 7919 * key)
    start = m.m068.phase_start_system(degs, len(target), phases)
    roots0 = m.m068.phase_start_roots(degs, phases)
    gammas = m.deterministic_gamma_vector(len(target), seed + 991 * Bc + 1291 * key, retry)
    target_gamma = m.scale_system(chart_target, gammas)
    tr = m.track_one_070(
        chart_target,
        start,
        target_gamma,
        roots0[idx],
        tol=args.tol,
        max_steps=args.max_steps,
        max_epochs=args.max_epochs,
        quad_cap=args.quad_cap,
        deadline=deadline,
    )
    if not (tr.ok or tr.t > getattr(args, "accept_t", 0.90)):
        return None
    y = mod.f076.polish_070_compat(chart_target, tr.z, args.tol, max(16, int(args.quad_cap)))
    if y is None:
        return None
    z = _from_projective_chart(y, mask, scales)
    if z is None:
        return None
    res = m.residual_norm(target, z)
    if not (math.isfinite(res) and res < args.residual_accept):
        zp = mod.f076.polish_070_compat(target, z, args.tol, max(16, int(args.quad_cap)))
        if zp is not None:
            z = zp
    res = m.residual_norm(target, z)
    return list(z) if math.isfinite(res) and res < args.residual_accept else None


def _projective_completion(mod, target: System, roots: list[Vector], args):
    n = len(target)
    if os.environ.get("PANDROSION_091_PROJECTIVE_OFF") == "1":
        return 0, 0, 0, 0.0, roots
    if n < 3 and not _env_bool("PANDROSION_091_PROJECTIVE_BIVARIATE", False):
        return 0, 0, 0, 0.0, roots

    B = mod.m.bezout(target)
    roots = _valid_cluster(mod, target, roots, args)
    if len(roots) >= B:
        return 0, 0, 0, 0.0, roots

    retries = max(0, _env_int("PANDROSION_091_PROJECTIVE_RETRIES", 1))
    limit = max(1, _env_int("PANDROSION_091_PROJECTIVE_LIMIT", min(B, 48)))
    budget = max(0.0, _env_float("PANDROSION_091_PROJECTIVE_BUDGET", 8.0))
    if retries <= 0 or budget <= 0.0:
        return 0, 0, 0, 0.0, roots

    scales = _riemann_scales(roots, n)
    masks = _chart_masks(roots, n)
    seed = mod.m.seed_for(args.family, n, mod.parse_case(args.case)[1], args.seed_index)
    sep = float(getattr(args, "cluster_sep", 1e-6))
    attempts = hits = charts_used = 0
    t0 = time.time()
    deadline = t0 + budget

    print(
        f"091-projective-riemann charts={len(masks)} retries={retries} "
        f"limit={limit} budget={budget:.1f}s scales={','.join(f'{s:.4g}' for s in scales)}",
        flush=True,
    )

    for ci, mask in enumerate(masks):
        if time.time() >= deadline:
            break
        chart_target = _projective_invert_system(target, mask, scales)
        degs = mod.m.degrees(chart_target)
        Bc = math.prod(degs)
        order = mod.golden_order(Bc)[:min(limit, Bc)]
        charts_used += 1
        print(
            f"  chart={ci + 1} mask={''.join('S' if b else 'N' for b in mask)} "
            f"degs={','.join(str(d) for d in degs)} paths<= {len(order)}",
            flush=True,
        )
        for retry in range(1, retries + 1):
            for idx in order:
                if time.time() >= deadline:
                    return attempts, hits, charts_used, time.time() - t0, roots
                attempts += 1
                try:
                    z = _projective_track_one(mod, target, chart_target, mask, scales, seed, idx, retry, args, deadline)
                except BaseException:
                    z = None
                if _unique_append(roots, z, sep):
                    hits += 1
                    print(
                        f"    projective hit chart={ci + 1} retry={retry} idx={idx} roots={len(roots)}/{B}",
                        flush=True,
                    )
                    if getattr(args, "stop_at_bezout", False) and len(roots) >= B:
                        return attempts, hits, charts_used, time.time() - t0, roots
    return attempts, hits, charts_used, time.time() - t0, roots


def install_projective_completion(mod) -> None:
    original_run_case = mod.run_case

    def run_case_091(args, case: str):
        summaries, scores, batches, roots = original_run_case(args, case)
        n, d = mod.parse_case(case)
        seed = mod.m.seed_for(args.family, n, d, args.seed_index)
        target = mod.m.gen_system(args.family, n, d, seed)
        B = mod.m.bezout(target)
        primary = summaries[0] if summaries else None
        attempts, hits, charts_used, sec, roots2 = (
            _projective_completion(mod, target, roots, args)
            if len(roots) < B else (0, 0, 0, 0.0, roots)
        )
        roots[:] = roots2

        if primary is not None:
            residuals = [mod.m.residual_norm(target, z) for z in roots]
            maxres = max(residuals) if residuals else float("inf")
            status = "ok" if len(roots) >= B and maxres < args.residual_accept else "partial"
            summaries[0] = replace(
                primary,
                alg="091-pandrosion-projective-riemann",
                roots=len(roots),
                coverage=len(roots) / max(1, B),
                path_rows=primary.path_rows + attempts,
                candidates=primary.candidates + hits,
                batches=primary.batches + charts_used,
                max_residual=maxres,
                seconds_observed=primary.seconds_observed + sec,
                status=status,
                notes=(
                    f"{primary.notes}; 091_projective_attempts={attempts}; "
                    f"091_projective_hits={hits}; 091_projective_charts={charts_used}; "
                    f"091_projective_seconds={sec:.2f}; Riemann south-chart inversion"
                ),
            )
            print(
                f"    091-pandrosion-projective-riemann roots={len(roots)}/{B} "
                f"cov={100.0 * len(roots) / max(1, B):.1f}% "
                f"paths={primary.path_rows + attempts} maxres={maxres:.2e} "
                f"sec={primary.seconds_observed + sec:.2f} status={status}",
                flush=True,
            )
        return summaries, scores, batches, roots

    mod.run_case = run_case_091


def main() -> None:
    sys.argv = _consume_091_args(sys.argv)
    f090 = _load_090()
    sys.argv = f090._consume_090_args(sys.argv)
    f089 = f090._load_089()
    sys.argv = f089._consume_089_args(sys.argv)
    w087 = f089._load_087()
    sys.argv = w087._consume_branch_args(sys.argv)
    mod = w087.load_081()
    w087.install_branch_safe_tracker(mod)
    w087.install_compat_and_single_policy(mod)
    f089.install_resultant_recovery(mod)
    f090.install_multivariate_completion(mod)
    install_projective_completion(mod)
    sys.argv = f090._rewrite_argv_for_090(w087, sys.argv)
    mod.main()


if __name__ == "__main__":
    main()
