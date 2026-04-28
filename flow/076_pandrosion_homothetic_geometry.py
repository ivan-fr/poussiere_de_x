"""
FLOW 076 -- homothetic/scaled Pandrosion geometry.

076 transports the scaling principle of 0pandrosion_pth.tex to dD polynomial
systems.  In the scalar paper one writes x^(1/p)=A^(1/p)(x/A)^(1/p), applies
Pandrosion to the reduced ratio x/A close to 1, then reconstructs by
multiplication.  Here the analogue is a diagonal homothety

    z = S y,      G(y) = D F(S y),

where S is generated from the polynomial system by coefficient/support balancing
and D optionally normalizes equations.  Paths are tracked in y-space and roots
are reconstructed as z=S y.  The Pandrosion telescope identity is preserved: if
Q_G(a,y)(y-a)=G(y)-G(a), then Q_z=D^{-1}Q_G S^{-1} is a valid z-space slope.

The recovery plan is inherited from 075: deterministic low-discrepancy
micro-order plus early stopping.  On KS(2,8), system homothety plus equation
normalization validates all 64 roots in the base pass, without Newton-ELS.
"""
from __future__ import annotations

import argparse
import csv
import concurrent.futures
import importlib.util
import json
import math
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

Complex = complex
Vector = List[Complex]

HERE = Path(__file__).resolve().parent
ADAPTIVE_PATH = HERE / "070_pandrosion_system_adaptive_ks.py"
LAIREZ_REF_CANDIDATES = [
    HERE / "074_benchmark_consolidated.csv",
    HERE / "071_ks28.csv",
    HERE / "071_benchmark_consolidated.csv",
    HERE / "074_ks28_lairez.csv",
]


def load_adaptive():
    spec = importlib.util.spec_from_file_location("flow070_adaptive_for_076", str(ADAPTIVE_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {ADAPTIVE_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow070_adaptive_for_076"] = mod
    spec.loader.exec_module(mod)
    return mod


m = load_adaptive()


@dataclass
class BatchRow:
    stage: str
    retry: int
    indices: str
    paths: int
    candidates: int
    roots_before: int
    roots_after: int
    new_roots: int
    seconds: float
    status: str
    path_json: str


@dataclass
class SummaryRow:
    family: str
    n: int
    d: int
    seed: int
    terms: int
    bezout: int
    alg: str
    roots: int
    coverage: float
    path_rows: int
    candidates: int
    batches: int
    retries_used: int
    max_residual: float
    seconds_observed: float
    status: str
    notes: str


# -----------------------------------------------------------------------------
# small helpers
# -----------------------------------------------------------------------------


def parse_case(s: str) -> Tuple[int, int]:
    a, b = s.split(",")
    return int(a), int(b)


def encode_indices(indices: Sequence[int]) -> str:
    return ",".join(str(int(i)) for i in indices)


def parse_indices(s: str) -> List[int]:
    out: List[int] = []
    for part in str(s).replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            a, b = part.split(":", 1)
            out.extend(range(int(a), int(b)))
        else:
            out.append(int(part))
    # preserve order, remove duplicates
    seen = set()
    uniq = []
    for i in out:
        if i not in seen:
            seen.add(i)
            uniq.append(i)
    return uniq


def root_from_json(zjson) -> Vector | None:
    if zjson is None:
        return None
    return [complex(float(a), float(b)) for a, b in zjson]


def root_to_json(z: Sequence[Complex]):
    return [[float(v.real), float(v.imag)] for v in z]


# -----------------------------------------------------------------------------
# Homothetic scaling: dD analogue of the scaling principle in 0pandrosion.tex.
# -----------------------------------------------------------------------------


def parse_scales(s: str | None, n: int) -> List[float]:
    if not s:
        return [1.0 for _ in range(n)]
    vals: List[float] = []
    for part in str(s).replace(';', ',').split(','):
        part = part.strip()
        if part:
            vals.append(float(part))
    if len(vals) == 1 and n > 1:
        vals = vals * n
    if len(vals) != n:
        raise ValueError(f"expected {n} scales, got {len(vals)} from {s!r}")
    return [max(1e-12, float(v)) for v in vals]


def encode_scales(scales: Sequence[float]) -> str:
    return ','.join(f"{float(v):.12g}" for v in scales)


def _solve_real(A: List[List[float]], b: List[float], ridge: float = 0.0) -> List[float] | None:
    n = len(b)
    if n == 0:
        return []
    M = [list(row) for row in A]
    x = list(b)
    if ridge:
        for i in range(n):
            M[i][i] += ridge
    for k in range(n):
        piv = max(range(k, n), key=lambda i: abs(M[i][k]))
        if abs(M[piv][k]) < 1e-14:
            return None
        if piv != k:
            M[k], M[piv] = M[piv], M[k]
            x[k], x[piv] = x[piv], x[k]
        pv = M[k][k]
        for j in range(k, n):
            M[k][j] /= pv
        x[k] /= pv
        for i in range(n):
            if i == k:
                continue
            f = M[i][k]
            if f == 0.0:
                continue
            for j in range(k, n):
                M[i][j] -= f * M[k][j]
            x[i] -= f * x[k]
    return x


def system_diagonal_scales(polys: m.System, min_scale: float = 0.25, max_scale: float = 4.0,
                           strength: float = 1.0) -> List[float]:
    """Coefficient/support homothety chosen from the system only.

    We fit log-scales x_j and equation offsets b_i so that the magnitudes of
    coefficients of F(S y) are more balanced:

        log |c_{i,alpha}| + alpha.x + b_i \approx 0.

    This is a tropical/Cauchy-style proxy for placing the unknown roots in a
    near-unit box, the dD analogue of choosing A=floor(x^{1/p})^p in the scalar
    Pandrosion scaling principle.  The fit is intentionally conservative.
    """
    n = len(polys)
    if n == 0:
        return []
    unk = n + len(polys)
    A = [[0.0 for _ in range(unk)] for _ in range(unk)]
    b = [0.0 for _ in range(unk)]
    count = 0
    for i, poly in enumerate(polys):
        for alpha, coeff in poly.items():
            c = abs(coeff)
            if c <= 0.0:
                continue
            # Constant terms help equation normalization but not variable scale.
            row = [0.0 for _ in range(unk)]
            deg = sum(alpha)
            if deg == 0:
                continue
            for j, aj in enumerate(alpha):
                row[j] = float(aj)
            row[n + i] = 1.0
            rhs = -math.log(max(c, 1e-300))
            # Small terms are noisy; log compression prevents a single monomial
            # from dictating the homothety.
            w = 0.25 + min(2.0, math.sqrt(c) / (1.0 + math.sqrt(c)))
            for r in range(unk):
                b[r] += w * row[r] * rhs
                rr = row[r]
                if rr == 0.0:
                    continue
                for q in range(unk):
                    if row[q] != 0.0:
                        A[r][q] += w * rr * row[q]
            count += 1
    if count < n:
        return [1.0 for _ in range(n)]
    # Ridge on variable scales; stronger ridge on equation offsets only fixes
    # numerical gauge without forcing a scale.
    for j in range(n):
        A[j][j] += 1e-2
    for i in range(len(polys)):
        A[n + i][n + i] += 1e-6
    sol = _solve_real(A, b)
    if sol is None:
        return [1.0 for _ in range(n)]
    lo, hi = math.log(min_scale), math.log(max_scale)
    out = []
    for j in range(n):
        x = max(lo, min(hi, float(sol[j]) * float(strength)))
        # Avoid meaningless tiny perturbations around 1.
        if abs(x) < 0.06:
            x = 0.0
        out.append(math.exp(x))
    return out


def root_stat_scales(roots: Sequence[Vector], n: int, q: float = 0.75,
                     min_scale: float = 0.25, max_scale: float = 4.0,
                     strength: float = 1.0) -> List[float]:
    """Homothety estimated from already discovered sheets.

    This is the multistart analogue of picking beta close to the target root in
    the scalar paper.  Since the missing roots are unknown, use a robust quantile
    of the roots already obtained by the base pass.
    """
    if not roots:
        return [1.0 for _ in range(n)]
    out: List[float] = []
    qq = max(0.0, min(1.0, float(q)))
    for j in range(n):
        vals = sorted(max(1e-12, abs(z[j])) for z in roots if len(z) > j and math.isfinite(abs(z[j])))
        if not vals:
            out.append(1.0)
            continue
        idx = int(round(qq * (len(vals) - 1)))
        s = vals[idx]
        # Blend toward 1: strength=0 means no scaling, strength=1 full scale.
        log_s = float(strength) * math.log(max(1e-12, s))
        s = math.exp(max(math.log(min_scale), min(math.log(max_scale), log_s)))
        if abs(math.log(s)) < 0.06:
            s = 1.0
        out.append(s)
    return out


def combine_scales(system_s: Sequence[float], root_s: Sequence[float], mode: str,
                   min_scale: float, max_scale: float) -> List[float]:
    n = max(len(system_s), len(root_s))
    ss = list(system_s) if system_s else [1.0 for _ in range(n)]
    rr = list(root_s) if root_s else [1.0 for _ in range(n)]
    if len(ss) < n:
        ss += [1.0] * (n - len(ss))
    if len(rr) < n:
        rr += [1.0] * (n - len(rr))
    lo, hi = math.log(min_scale), math.log(max_scale)
    mode = mode.lower()
    out: List[float] = []
    for a, b0 in zip(ss, rr):
        if mode == 'none':
            x = 0.0
        elif mode == 'system':
            x = math.log(max(1e-12, a))
        elif mode == 'roots':
            x = math.log(max(1e-12, b0))
        else:  # hybrid: coefficient homothety first, root-stat correction second.
            x = 0.50 * math.log(max(1e-12, a)) + 0.75 * math.log(max(1e-12, b0))
        x = max(lo, min(hi, x))
        if abs(x) < 0.06:
            x = 0.0
        out.append(math.exp(x))
    return out


def variable_scale_system(polys: m.System, scales: Sequence[float], normalize_equations: bool = False) -> m.System:
    """Return G(y)=D F(S y), with optional equation normalization D.

    Roots transport by z=S y.  For every exact Pandrosion geometry computed in
    y-space, the z-space slope is Q_y S^{-1} and still satisfies the telescope
    identity for F.  The implementation tracks y directly and reconstructs z.
    """
    n = len(scales)
    out: m.System = []
    for poly in polys:
        q: m.Poly = {}
        for alpha, coeff in poly.items():
            factor = 1.0 + 0.0j
            for j, aj in enumerate(alpha):
                if aj:
                    factor *= float(scales[j]) ** aj
            q[alpha] = coeff * factor
        if normalize_equations:
            norm = max((abs(c) for c in q.values()), default=1.0)
            if norm > 0 and math.isfinite(norm):
                q = {a: c / norm for a, c in q.items()}
        out.append(q)
    return out


def vec_to_scaled(z: Sequence[Complex], scales: Sequence[float]) -> Vector:
    return [complex(z[j]) / float(scales[j]) for j in range(len(scales))]


def vec_from_scaled(y: Sequence[Complex], scales: Sequence[float]) -> Vector:
    return [complex(y[j]) * float(scales[j]) for j in range(len(scales))]


def cluster_roots(roots: Sequence[Vector], sep: float = 1e-6) -> List[Vector]:
    clusters: List[Vector] = []
    for z in roots:
        n = len(z)
        dup = False
        for r in clusters:
            scale = max(1.0, max(abs(z[i]) for i in range(n)), max(abs(r[i]) for i in range(n)))
            if max(abs(z[i] - r[i]) for i in range(n)) <= sep * scale:
                dup = True
                break
        if not dup:
            clusters.append(list(z))
    return clusters


def count_new_roots(current: Sequence[Vector], candidates: Sequence[Vector], sep: float) -> Tuple[List[Vector], int]:
    before = len(current)
    roots = list(current)
    for z in candidates:
        roots = cluster_roots([*roots, z], sep=sep)
    return roots, len(roots) - before


def read_batch_candidates(batch_json: str) -> Tuple[List[Vector], int, int, List[dict]]:
    p = Path(batch_json)
    if not p.exists():
        return [], 0, 0, []
    data = json.loads(p.read_text())
    rows = list(data.get("rows", []))
    roots: List[Vector] = []
    candidates = 0
    for r in rows:
        z = root_from_json(r.get("z"))
        if z is not None:
            candidates += 1
            roots.append(z)
    return roots, len(rows), candidates, rows


def read_all_roots(batch_rows: Sequence[BatchRow]) -> Tuple[List[Vector], int, int]:
    roots: List[Vector] = []
    path_rows = 0
    candidates = 0
    for b in batch_rows:
        rs, pr, cand, _ = read_batch_candidates(b.path_json)
        roots.extend(rs)
        path_rows += pr
        candidates += cand
    return roots, path_rows, candidates


# -----------------------------------------------------------------------------
# deterministic predictive recovery orders
# -----------------------------------------------------------------------------


def golden_order(B: int) -> List[int]:
    """Low-discrepancy order on the start-index torus.

    For B=64 this starts [0, 39, 15, 54, ...].  The first two entries are the
    two useful retry-1 sentinels in the 074 KS(2,8) validation, but the rule is
    deterministic and size-scaled, not hard-coded to KS(2,8).
    """
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    out: List[int] = []
    seen = set()
    k = 0
    while len(out) < B and k < 4 * B + 10:
        if k == 0:
            idx = 0
        else:
            idx = int(math.floor(B * ((k * phi) % 1.0)))
        idx = max(0, min(B - 1, idx))
        if idx not in seen:
            seen.add(idx)
            out.append(idx)
        k += 1
    # Fill any holes caused by flooring collisions.
    for idx in range(B):
        if idx not in seen:
            out.append(idx)
    return out


def chunk_ranges(B: int, chunk_size: int) -> List[List[int]]:
    return [list(range(s, min(B, s + chunk_size))) for s in range(0, B, chunk_size)]


def stratified_window_order(B: int, block: int) -> List[List[int]]:
    block = max(1, min(B, int(block)))
    starts = [0, B // 2, B // 4, (3 * B) // 4]
    starts.extend(range(0, B, block))
    out: List[List[int]] = []
    seen = set()
    for s in starts:
        s = max(0, min(B, int(s)))
        e = min(B, s + block)
        if s >= e:
            continue
        key = (s, e)
        if key in seen:
            continue
        seen.add(key)
        out.append(list(range(s, e)))
    return out


def anomaly_order_from_rows(base_rows: Sequence[dict], B: int) -> List[int]:
    """System/path diagnostics from the base pass.

    Failures, very slow paths, and high final residuals are not guaranteed to
    identify the missing sheets under a new gamma, so this is appended after the
    golden sentinels rather than used alone.
    """
    scored = []
    for r in base_rows:
        idx = int(r.get("idx", -1))
        if idx < 0 or idx >= B:
            continue
        ok = bool(r.get("z") is not None)
        t = float(r.get("t", 0.0) or 0.0)
        res = float(r.get("residual", float("inf")) or float("inf"))
        sec = float(r.get("seconds", 0.0) or 0.0)
        score = 0.0
        if not ok:
            score += 1000.0
        score += max(0.0, 1.0 - min(1.0, t)) * 100.0
        score += min(50.0, math.log10(1.0 + max(0.0, res)))
        score += sec
        scored.append((score, idx))
    scored.sort(reverse=True)
    out: List[int] = []
    seen = set()
    for _, idx in scored:
        for j in (idx, idx - 1, idx + 1, idx - 2, idx + 2):
            if 0 <= j < B and j not in seen:
                seen.add(j)
                out.append(j)
    return out


def combined_micro_order(base_rows: Sequence[dict], B: int) -> List[int]:
    pieces = [golden_order(B), anomaly_order_from_rows(base_rows, B)]
    out: List[int] = []
    seen = set()
    for piece in pieces:
        for idx in piece:
            if idx not in seen:
                seen.add(idx)
                out.append(idx)
    return out


# -----------------------------------------------------------------------------
# worker + orchestration
# -----------------------------------------------------------------------------


def worker_run(args: argparse.Namespace) -> None:
    n, d = parse_case(args.case)
    seed = m.seed_for(args.family, n, d, args.seed_index)
    target = m.gen_system(args.family, n, d, seed)
    degs = m.degrees(target)
    B = math.prod(degs)
    retry = int(args.retry)
    indices = [i for i in parse_indices(args.indices) if 0 <= i < B]

    scales = parse_scales(getattr(args, "scales", ""), n)
    use_scaled = any(abs(math.log(max(1e-12, s))) > 1e-12 for s in scales)

    # Match the 074 chunked deterministic structure, then transport the whole
    # homotopy by z=S y.  This is the dD homothety analogue of x=A x'.
    phases = m.m068.deterministic_phases(n, 0, seed + 17 * B + 10007 * retry)
    start_z = m.m068.phase_start_system(degs, n, phases)
    roots0_z = m.m068.phase_start_roots(degs, phases)
    if use_scaled:
        target_track = variable_scale_system(target, scales, normalize_equations=args.equation_normalize)
        start_track = variable_scale_system(start_z, scales, normalize_equations=args.equation_normalize)
        roots0 = [vec_to_scaled(r, scales) for r in roots0_z]
    else:
        target_track = target
        start_track = start_z
        roots0 = roots0_z
    gammas = m.deterministic_gamma_vector(n, seed + 991 * B, retry)
    target_gamma = m.scale_system(target_track, gammas)

    rows = []
    t0 = time.time()
    for idx in indices:
        pt0 = time.time()
        z_polished = None
        y_polished = None
        tr_payload: Dict[str, object] = {}
        try:
            track_kwargs = {}
            if args.dt0 and args.dt0 > 0:
                track_kwargs["dt0"] = float(args.dt0)
            if args.dtmax and args.dtmax > 0:
                track_kwargs["dtmax"] = float(args.dtmax)
            tr = m.track_one_070(target_track, start_track, target_gamma, roots0[idx],
                                 tol=args.tol, max_steps=args.max_steps,
                                 max_epochs=args.max_epochs, quad_cap=args.quad_cap,
                                 **track_kwargs)
            tr_payload = {
                "ok": bool(tr.ok), "t": float(tr.t), "residual_scaled": float(tr.residual),
                "steps": int(tr.steps), "epochs": int(tr.epochs),
            }
            if tr.ok or tr.t > 0.90:
                y_polished = m.polish_070(target_track, tr.z, args.tol, args.quad_cap)
                if y_polished is not None:
                    z_polished = vec_from_scaled(y_polished, scales) if use_scaled else list(y_polished)
                    # Verify in the original coordinates.  A short original
                    # Pandrosion polish is allowed: it is still derivative-free
                    # and prevents equation-normalization artifacts.
                    rz = m.residual_norm(target, z_polished)
                    if not (math.isfinite(rz) and rz < 1e-7):
                        zp2 = m.polish_070(target, z_polished, args.tol, args.quad_cap)
                        if zp2 is not None:
                            z_polished = zp2
                    tr_payload["residual"] = float(m.residual_norm(target, z_polished)) if z_polished is not None else float("inf")
                else:
                    tr_payload["residual"] = float("inf")
            else:
                tr_payload["residual"] = float("inf")
            rows.append({
                "idx": int(idx),
                **tr_payload,
                "seconds": float(time.time() - pt0),
                "scale": encode_scales(scales),
                "z": (root_to_json(z_polished) if z_polished is not None else None),
            })
        except BaseException as exc:
            rows.append({
                "idx": int(idx), "ok": False, "t": 0.0,
                "residual": float("inf"), "residual_scaled": float("inf"),
                "steps": 0, "epochs": 0,
                "seconds": float(time.time() - pt0), "scale": encode_scales(scales), "z": None,
                "error": repr(exc),
            })
        if args.verbose:
            r = rows[-1]
            print(f"retry={retry} idx={idx} scale={r.get('scale')} ok={r.get('ok')} t={r.get('t'):.3f} "
                  f"res={r.get('residual'):.2e} z={r.get('z') is not None} sec={r.get('seconds'):.2f}", flush=True)
    out = {
        "family": args.family, "case": args.case, "seed_index": args.seed_index,
        "seed": seed, "retry": retry, "indices": indices, "rows": rows,
        "scales": scales, "equation_normalize": bool(args.equation_normalize),
        "seconds": float(time.time() - t0),
    }
    Path(args.out).write_text(json.dumps(out, indent=2))


def run_batch(args: argparse.Namespace, stage: str, retry: int, indices: Sequence[int],
              outdir: Path, roots_before: Sequence[Vector], sep: float,
              scales: Sequence[float] | None = None) -> Tuple[BatchRow, List[dict]]:
    outdir.mkdir(parents=True, exist_ok=True)
    Btag = f"{args.family}_{args.case.replace(',', 'x')}_seed{args.seed_index}"
    if len(indices) <= 8:
        itag = "i" + "_".join(str(i) for i in indices)
    else:
        itag = f"{min(indices)}_{max(indices)+1}_n{len(indices)}"
    outfile = outdir / f"076_{Btag}_{stage}_r{retry}_{itag}.json"
    cmd = [
        sys.executable, "-S", str(Path(__file__).resolve()),
        "--worker", "--family", args.family, "--case", args.case,
        "--seed-index", str(args.seed_index), "--retry", str(retry),
        "--indices", encode_indices(indices), "--out", str(outfile),
        "--tol", str(args.tol), "--max-steps", str(args.max_steps),
        "--max-epochs", str(args.max_epochs), "--quad-cap", str(args.quad_cap),
        "--dt0", str(args.dt0), "--dtmax", str(args.dtmax),
        "--scales", encode_scales(scales if scales is not None else getattr(args, "active_scales", [])),
    ]
    if args.equation_normalize:
        cmd.append("--equation-normalize")
    if args.verbose:
        cmd.append("--verbose")
    t0 = time.time()
    status = "ok"
    try:
        subprocess.run(cmd, timeout=args.batch_timeout if args.batch_timeout > 0 else None,
                       check=True, stdout=(None if args.verbose else subprocess.DEVNULL),
                       stderr=(None if args.verbose else subprocess.DEVNULL))
    except subprocess.TimeoutExpired:
        status = "timeout"
    except subprocess.CalledProcessError as exc:
        status = f"error:{exc.returncode}"
    seconds = time.time() - t0
    candidates_roots, paths, candidates, rows = read_batch_candidates(str(outfile))
    if outfile.exists():
        try:
            data = json.loads(outfile.read_text())
            seconds = float(data.get("seconds", seconds))
        except Exception:
            pass
    after_roots, new = count_new_roots(roots_before, candidates_roots, sep=sep)
    br = BatchRow(stage=stage, retry=retry, indices=encode_indices(indices),
                  paths=paths if paths else len(indices), candidates=candidates,
                  roots_before=len(roots_before), roots_after=len(after_roots),
                  new_roots=new, seconds=float(seconds), status=status,
                  path_json=str(outfile))
    return br, rows




def recompute_batch_progress(batch_rows: Sequence[BatchRow], sep: float) -> Tuple[List[BatchRow], List[Vector], int, int]:
    """Recompute before/after/new fields in deterministic batch order.

    Parallel base batches finish out of order and cannot know their contribution
    while running.  This pass rebuilds the cumulative clustering exactly as the
    sequential orchestrator would report it.
    """
    roots: List[Vector] = []
    fixed: List[BatchRow] = []
    total_paths = 0
    total_candidates = 0
    for b in batch_rows:
        cand_roots, paths, candidates, _ = read_batch_candidates(b.path_json)
        before = len(roots)
        roots, new = count_new_roots(roots, cand_roots, sep=sep)
        fixed.append(BatchRow(stage=b.stage, retry=b.retry, indices=b.indices,
                              paths=paths if paths else b.paths, candidates=candidates,
                              roots_before=before, roots_after=len(roots),
                              new_roots=new, seconds=b.seconds, status=b.status,
                              path_json=b.path_json))
        total_paths += paths if paths else b.paths
        total_candidates += candidates
    return fixed, roots, total_paths, total_candidates


def run_base_batches(args: argparse.Namespace, B: int, outdir: Path, base_scales: Sequence[float]) -> Tuple[List[BatchRow], List[dict]]:
    chunks = chunk_ranges(B, args.base_chunk_size)
    trace_rows: List[dict] = []

    # Sequential baseline, unless the user explicitly requests concurrency.
    if args.parallel_base <= 1 or len(chunks) <= 1:
        base_rows: List[BatchRow] = []
        roots_seq: List[Vector] = []
        for inds in chunks:
            before = list(roots_seq)
            br, rows = run_batch(args, "base", 0, inds, outdir, before, args.cluster_sep, scales=base_scales)
            base_rows.append(br)
            trace_rows.extend(rows)
            roots_seq, _, _ = read_all_roots(base_rows)
            roots_seq = cluster_roots(roots_seq, sep=args.cluster_sep)
            print(f"base retry=0 indices={inds[0]}:{inds[-1]+1} candidates={br.candidates}/{br.paths} roots={len(roots_seq)}/{B} sec={br.seconds:.2f}", flush=True)
            if args.stop_at_bezout and len(roots_seq) >= B:
                break
        return base_rows, trace_rows

    # Concurrent base chunks.  This does not change the mathematics; it only
    # runs independent path-index chunks in parallel subprocesses.
    jobs = max(1, int(args.parallel_base))
    spec_indices: List[int] = []
    if args.speculative_micro and args.retries > 1:
        spec_indices = golden_order(B)[:max(1, int(args.micro_batch))]
        # One extra worker may be used for the tiny speculative retry batch.
        jobs += 1

    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
        future_items = []
        for pos, inds in enumerate(chunks):
            fut = ex.submit(run_batch, args, "base", 0, inds, outdir, [], args.cluster_sep, base_scales)
            future_items.append(("base", pos, inds, fut))
        if spec_indices:
            fut = ex.submit(run_batch, args, "micro-spec", 1, spec_indices, outdir, [], args.cluster_sep, base_scales)
            future_items.append(("micro-spec", len(chunks), spec_indices, fut))

        tmp: Dict[Tuple[str, int], Tuple[BatchRow, List[dict], List[int]]] = {}
        for stage, pos, inds, fut in future_items:
            br, rows = fut.result()
            tmp[(stage, pos)] = (br, rows, inds)

    ordered_base = [tmp[("base", i)] for i in range(len(chunks))]
    raw_base_rows = [x[0] for x in ordered_base]
    for _, rows, _ in ordered_base:
        trace_rows.extend(rows)
    fixed_base_rows, roots_after_base, _, _ = recompute_batch_progress(raw_base_rows, args.cluster_sep)
    for br, (_, _, inds) in zip(fixed_base_rows, ordered_base):
        print(f"base-parallel retry=0 indices={inds[0]}:{inds[-1]+1} candidates={br.candidates}/{br.paths} roots={br.roots_after}/{B} sec={br.seconds:.2f}", flush=True)

    all_initial = list(fixed_base_rows)
    if spec_indices:
        spec_br, _spec_rows, _ = tmp[("micro-spec", len(chunks))]
        # Recompute the speculative batch contribution after the base rows.
        recomputed, roots_after_spec, _, _ = recompute_batch_progress([*all_initial, spec_br], args.cluster_sep)
        spec_fixed = recomputed[-1]
        all_initial.append(spec_fixed)
        print(f"micro-spec retry=1 indices={encode_indices(spec_indices)} candidates={spec_fixed.candidates}/{spec_fixed.paths} +{spec_fixed.new_roots} roots={spec_fixed.roots_after}/{B} sec={spec_fixed.seconds:.2f}", flush=True)
    return all_initial, trace_rows


def load_lairez_reference(family: str, n: int, d: int, seed: int, B: int) -> SummaryRow | None:
    """Load the strongest matching local Lairez-style reference row.

    Several earlier experiments saved Lairez rows with different retry budgets.
    For reporting 075 head-to-head, prefer the best coverage, then lower time.
    """
    candidates: List[SummaryRow] = []
    for path in LAIREZ_REF_CANDIDATES:
        if not path.exists():
            continue
        try:
            with path.open() as f:
                rdr = csv.DictReader(f)
                for r in rdr:
                    alg = r.get("alg", "")
                    if "lairez" not in alg:
                        continue
                    if r.get("family") and r.get("family") != family:
                        continue
                    if r.get("n") and int(float(r.get("n"))) != n:
                        continue
                    if r.get("d") and int(float(r.get("d"))) != d:
                        continue
                    roots = int(float(r.get("roots", 0) or 0))
                    cov = float(r.get("coverage", 0.0) or 0.0)
                    row = SummaryRow(
                        family=family, n=n, d=d, seed=seed,
                        terms=int(float(r.get("terms", 0) or 0)), bezout=B,
                        alg="lairez-style-reference", roots=roots, coverage=cov,
                        path_rows=int(float(r.get("paths", r.get("path_rows", 0)) or 0)),
                        candidates=roots, batches=0,
                        retries_used=int(float(r.get("gamma_retries", r.get("retries_used", 0)) or 0)),
                        max_residual=float(r.get("max_residual", float("inf")) or float("inf")),
                        seconds_observed=float(r.get("seconds", r.get("seconds_observed", 0.0)) or 0.0),
                        status=str(r.get("status", "ok" if cov >= .999999 else "partial")),
                        notes=f"loaded best local reference from {path.name}: gamma total-degree + analytic Newton",
                    )
                    candidates.append(row)
        except Exception:
            continue
    if not candidates:
        return None
    candidates.sort(key=lambda r: (-r.coverage, r.seconds_observed))
    return candidates[0]


def run_lairez_now(args: argparse.Namespace, target, seed: int, terms: int, B: int) -> SummaryRow:
    t0 = time.time()
    res = m.m069.run_lairez(target, seed=91000 + seed, tol=args.tol,
                            max_steps=args.lairez_max_steps,
                            max_newton_iter=args.lairez_newton_iters,
                            retries=args.lairez_retries)
    if isinstance(res, tuple):
        resl, seconds = res
    else:
        resl, seconds = res, time.time() - t0
    if not isinstance(resl, dict):
        roots = getattr(resl, "roots", [])
        coverage = float(getattr(resl, "coverage", len(roots) / max(1, B)))
        row = SummaryRow(args.family, *parse_case(args.case), seed, terms, B,
                         "lairez-style", len(roots), coverage,
                         int(getattr(resl, "paths", 0)), len(roots), 0, args.lairez_retries,
                         float(getattr(resl, "max_residual", float("inf"))),
                         float(seconds), "ok" if coverage >= .999999 else "partial",
                         "fresh run: gamma total-degree + analytic Newton")
        return row
    roots = resl.get("roots", [])
    return SummaryRow(
        family=args.family, n=parse_case(args.case)[0], d=parse_case(args.case)[1],
        seed=seed, terms=terms, bezout=B, alg="lairez-style",
        roots=len(roots), coverage=float(resl.get("coverage", len(roots) / max(1, B))),
        path_rows=int(resl.get("paths", 0)), candidates=len(roots), batches=0,
        retries_used=args.lairez_retries,
        max_residual=float(resl.get("max_residual", float("inf"))),
        seconds_observed=float(seconds),
        status=str(resl.get("status", "ok" if len(roots) >= B else "partial")),
        notes="fresh run: gamma total-degree + analytic Newton",
    )


def orchestrate(args: argparse.Namespace) -> None:
    n, d = parse_case(args.case)
    seed = m.seed_for(args.family, n, d, args.seed_index)
    target = m.gen_system(args.family, n, d, seed)
    B = m.bezout(target)
    terms = m.term_count(target)
    outdir = Path(args.outdir or "/mnt/data/076_batches")
    batch_rows: List[BatchRow] = []
    base_trace_rows: List[dict] = []

    system_scales = system_diagonal_scales(target, min_scale=args.scale_min, max_scale=args.scale_max,
                                           strength=args.system_scale_strength)
    base_root_scales = [1.0 for _ in range(n)]
    base_scales = combine_scales(system_scales, base_root_scales, args.homothety, args.scale_min, args.scale_max)
    # For the pure roots mode, do not scale the first pass because no roots have
    # been discovered yet.
    if args.homothety.lower() == "roots":
        base_scales = [1.0 for _ in range(n)]

    print("=" * 128)
    print("076 -- homothetic system-generated Pandrosion geometry")
    print("=" * 128)
    print(f"family={args.family}, case=({n},{d}), seed={seed}, terms={terms}, Bezout={B}")
    print(f"base_chunk={args.base_chunk_size}, parallel_base={args.parallel_base}, speculative_micro={args.speculative_micro}, micro_batch={args.micro_batch}, micro_limit={args.micro_limit}, retries={args.retries}, dt0={args.dt0 or 'default'}, dtmax={args.dtmax or 'default'}")
    print(f"homothety={args.homothety}, system_scales={encode_scales(system_scales)}, base_scales={encode_scales(base_scales)}, eq_norm={args.equation_normalize}")
    print("076 uses the 070/074 system-generated Pandrosion tracker on homothetically scaled systems; no Newton-ELS.", flush=True)

    t0 = time.time()
    roots: List[Vector] = []

    # Base retry 0: complete coverage attempt, chunked for isolation.
    base_batches, base_trace_rows = run_base_batches(args, B, outdir, base_scales)
    batch_rows.extend(base_batches)
    all_roots, _, _ = read_all_roots(batch_rows)
    roots = cluster_roots(all_roots, sep=args.cluster_sep)
    root_scales_raw = root_stat_scales(roots, n, q=args.root_scale_quantile,
                                       min_scale=args.scale_min, max_scale=args.scale_max,
                                       strength=args.root_scale_strength)
    trig = max(1.0, float(args.root_scale_trigger))
    root_scales = [1.0 if (1.0 / trig <= s <= trig) else s for s in root_scales_raw]
    retry_scales = combine_scales(system_scales, root_scales, args.homothety, args.scale_min, args.scale_max)
    if args.homothety.lower() == "none":
        retry_scales = [1.0 for _ in range(n)]
    print(f"root_scales_raw={encode_scales(root_scales_raw)}, root_scales={encode_scales(root_scales)}, retry_scales={encode_scales(retry_scales)}", flush=True)

    # Micro-recovery: tiny batches in generated order. Stop as soon as full.
    retries_used = 1
    for retry in range(1, max(1, args.retries)):
        if args.stop_at_bezout and len(roots) >= B:
            break
        retries_used = retry + 1
        order = combined_micro_order(base_trace_rows, B)
        limit = min(B, max(0, args.micro_limit))
        tried_retry = set()
        for b in batch_rows:
            if b.retry == retry:
                tried_retry.update(parse_indices(b.indices))
        micro_order = [i for i in (order[:limit] if limit > 0 else []) if i not in tried_retry]
        if micro_order:
            for s in range(0, len(micro_order), args.micro_batch):
                if args.stop_at_bezout and len(roots) >= B:
                    break
                inds = micro_order[s:s + args.micro_batch]
                before = list(roots)
                br, _rows = run_batch(args, "micro", retry, inds, outdir, before, args.cluster_sep, scales=retry_scales)
                batch_rows.append(br)
                all_roots, _, _ = read_all_roots(batch_rows)
                roots = cluster_roots(all_roots, sep=args.cluster_sep)
                print(f"micro retry={retry} indices={encode_indices(inds)} candidates={br.candidates}/{br.paths} +{br.new_roots} roots={len(roots)}/{B} sec={br.seconds:.2f}", flush=True)
                if args.stop_at_bezout and len(roots) >= B:
                    break

        # Fallback to 074-style stratified windows if micro-order is insufficient.
        if len(roots) < B and args.window_fallback:
            for inds in stratified_window_order(B, args.recovery_block):
                # Avoid exact repeats of indices already tried in this retry.
                tried = set()
                for b in batch_rows:
                    if b.retry == retry:
                        tried.update(parse_indices(b.indices))
                inds2 = [i for i in inds if i not in tried]
                if not inds2:
                    continue
                before = list(roots)
                br, _rows = run_batch(args, "window", retry, inds2, outdir, before, args.cluster_sep, scales=retry_scales)
                batch_rows.append(br)
                all_roots, _, _ = read_all_roots(batch_rows)
                roots = cluster_roots(all_roots, sep=args.cluster_sep)
                print(f"window retry={retry} indices={inds2[0]}:{inds2[-1]+1} n={len(inds2)} candidates={br.candidates}/{br.paths} +{br.new_roots} roots={len(roots)}/{B} sec={br.seconds:.2f}", flush=True)
                if args.stop_at_bezout and len(roots) >= B:
                    break

    all_roots, path_rows, candidates = read_all_roots(batch_rows)
    roots = cluster_roots(all_roots, sep=args.cluster_sep)
    residuals = [m.residual_norm(target, z) for z in roots]
    max_res = max(residuals) if residuals else float("inf")
    status = "ok" if len(roots) >= B and max_res < 1e-7 else "partial"
    sec_total = time.time() - t0
    sec_observed = sum(b.seconds for b in batch_rows)
    # observed path seconds excludes subprocess startup; wall seconds includes it.
    notes = (f"golden micro-order first; observed_path_seconds={sec_observed:.2f}; "
             f"base_chunks={args.base_chunk_size}; parallel_base={args.parallel_base}; speculative_micro={args.speculative_micro}; micro_limit={args.micro_limit}; "
             f"fallback={'on' if args.window_fallback else 'off'}; homothety={args.homothety}; "
             f"base_scales={encode_scales(base_scales)}; root_scales_raw={encode_scales(root_scales_raw)}; root_scales={encode_scales(root_scales)}; retry_scales={encode_scales(retry_scales)}; "
             f"eq_norm={args.equation_normalize}; dt0={args.dt0 or 'default'}; dtmax={args.dtmax or 'default'}; no Newton-ELS")
    row075 = SummaryRow(args.family, n, d, seed, terms, B, "076-homothetic-system-geometry",
                        len(roots), len(roots) / max(1, B), path_rows, candidates,
                        len(batch_rows), retries_used, max_res, sec_total, status, notes)
    summary_rows = [row075]
    if args.include_lairez_reference:
        ref = load_lairez_reference(args.family, n, d, seed, B)
        if ref is not None:
            summary_rows.append(ref)
    if args.run_lairez:
        summary_rows.append(run_lairez_now(args, target, seed, terms, B))

    print("-" * 128)
    for r in summary_rows:
        print(f"{r.alg:>32} roots={r.roots}/{r.bezout} cov={100*r.coverage:.1f}% paths={r.path_rows} batches={r.batches} maxres={r.max_residual:.2e} sec={r.seconds_observed:.2f} status={r.status}", flush=True)

    write_outputs(summary_rows, batch_rows, roots, args)


# -----------------------------------------------------------------------------
# outputs
# -----------------------------------------------------------------------------


def write_outputs(summary_rows: Sequence[SummaryRow], batch_rows: Sequence[BatchRow], roots: Sequence[Vector], args: argparse.Namespace) -> None:
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(asdict(summary_rows[0]).keys()))
            w.writeheader()
            for r in summary_rows:
                w.writerow(asdict(r))
    if args.batch_csv and batch_rows:
        with open(args.batch_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(asdict(batch_rows[0]).keys()))
            w.writeheader()
            for r in batch_rows:
                w.writerow(asdict(r))
    if args.roots_json:
        payload = {
            "family": args.family, "case": args.case, "roots": len(roots),
            "clusters": [root_to_json(z) for z in roots],
        }
        Path(args.roots_json).write_text(json.dumps(payload, indent=2))
    if args.md:
        with open(args.md, "w") as f:
            f.write("# Flow 076 benchmark\n\n")
            f.write("076 keeps the system-generated Pandrosion geometry from 074/075, and adds a diagonal homothety z=S y inspired by the scaling principle of 0pandrosion. No Newton-ELS is used by 076.\n\n")
            f.write("## Summary\n\n")
            f.write("| alg | case | Bezout | roots | coverage | paths | batches | retries | max residual | seconds | status | notes |\n")
            f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|\n")
            for r in summary_rows:
                note = r.notes.replace("|", "/")[:180]
                f.write(f"| {r.alg} | ({r.n},{r.d}) | {r.bezout} | {r.roots} | {100*r.coverage:.1f}% | {r.path_rows} | {r.batches} | {r.retries_used} | {r.max_residual:.2e} | {r.seconds_observed:.2f} | {r.status} | {note} |\n")
            if batch_rows:
                f.write("\n## Batch details\n\n")
                f.write("| stage | retry | indices | paths | candidates | roots before | roots after | new roots | seconds | status |\n")
                f.write("|---|---:|---|---:|---:|---:|---:|---:|---:|---|\n")
                for b in batch_rows:
                    inds = b.indices if len(b.indices) <= 70 else b.indices[:67] + "..."
                    f.write(f"| {b.stage} | {b.retry} | `{inds}` | {b.paths} | {b.candidates} | {b.roots_before} | {b.roots_after} | {b.new_roots} | {b.seconds:.2f} | {b.status} |\n")
            f.write("\n## Geometry note\n\n")
            f.write("Each tracked path uses exact Pandrosion slopes generated from the polynomial system profile, preserving `Q(a,z)(z-a)=F(z)-F(a)` up to roundoff. The homothety changes the coordinate geometry used by the paths; micro-order only changes which start paths are retried.\n")


def main() -> None:
    ap = argparse.ArgumentParser(description="076 homothetic predictive micro-recovery for generated Pandrosion geometry")
    ap.add_argument("--family", default="ks")
    ap.add_argument("--case", default="2,8")
    ap.add_argument("--seed-index", type=int, default=0)
    ap.add_argument("--retries", type=int, default=2)
    ap.add_argument("--base-chunk-size", type=int, default=32)
    ap.add_argument("--parallel-base", type=int, default=1, help="run base chunks concurrently in isolated subprocesses")
    ap.add_argument("--speculative-micro", action="store_true", help="run the first retry micro-batch concurrently with the base chunks")
    ap.add_argument("--micro-batch", type=int, default=2)
    ap.add_argument("--micro-limit", type=int, default=16)
    ap.add_argument("--window-fallback", action="store_true", default=True)
    ap.add_argument("--no-window-fallback", dest="window_fallback", action="store_false")
    ap.add_argument("--recovery-block", type=int, default=16)
    ap.add_argument("--batch-timeout", type=float, default=120.0)
    ap.add_argument("--stop-at-bezout", action="store_true")
    ap.add_argument("--cluster-sep", type=float, default=1e-6)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--max-steps", type=int, default=120)
    ap.add_argument("--max-epochs", type=int, default=4)
    ap.add_argument("--quad-cap", type=int, default=12)
    ap.add_argument("--dt0", type=float, default=0.0, help="optional initial homotopy step for track_one_070; 0 uses tracker default")
    ap.add_argument("--dtmax", type=float, default=0.0, help="optional maximum homotopy step for track_one_070; 0 uses tracker default")
    ap.add_argument("--homothety", choices=["none", "system", "roots", "hybrid"], default="hybrid",
                    help="diagonal scaling policy z=S y: system coefficients, root quantiles after base, or hybrid")
    ap.add_argument("--scale-min", type=float, default=0.25)
    ap.add_argument("--scale-max", type=float, default=4.0)
    ap.add_argument("--system-scale-strength", type=float, default=1.0)
    ap.add_argument("--root-scale-strength", type=float, default=0.60)
    ap.add_argument("--root-scale-quantile", type=float, default=0.75)
    ap.add_argument("--root-scale-trigger", type=float, default=1.75,
                    help="ignore root-stat homothety unless the robust root scale is outside [1/T,T]")
    ap.add_argument("--equation-normalize", action="store_true", help="normalize equations after variable homothety")
    ap.add_argument("--outdir", default="/mnt/data/076_batches")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--batch-csv", default=None)
    ap.add_argument("--md", default=None)
    ap.add_argument("--roots-json", default=None)
    ap.add_argument("--include-lairez-reference", action="store_true")
    ap.add_argument("--run-lairez", action="store_true")
    ap.add_argument("--lairez-max-steps", type=int, default=420)
    ap.add_argument("--lairez-newton-iters", type=int, default=12)
    ap.add_argument("--lairez-retries", type=int, default=2)
    ap.add_argument("--verbose", action="store_true")
    # worker-only
    ap.add_argument("--worker", action="store_true")
    ap.add_argument("--retry", type=int, default=0)
    ap.add_argument("--indices", default="")
    ap.add_argument("--scales", default="")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    args.base_chunk_size = max(1, int(args.base_chunk_size))
    args.micro_batch = max(1, int(args.micro_batch))
    args.scale_min = max(1e-6, float(args.scale_min))
    args.scale_max = max(args.scale_min, float(args.scale_max))
    if args.worker:
        if not args.out:
            raise SystemExit("--worker requires --out")
        worker_run(args)
    else:
        orchestrate(args)


if __name__ == "__main__":
    main()
