"""
FLOW 078 -- Pandrosion-polygon homothetic geometry for dD systems.

078 keeps the tracker, batching, and output contract of flow/076.  Scaling is
the 076 diagonal homothety.  In pandrosion_vs_lairez.py, k_optimal_dyadic(coeffs)
uses the Pandrosion/Jensen polygon of a univariate polynomial to choose a
dyadic shell factor k, then runs the scaled problem Q(w)=P((R_geo/k)w).

For a dD system we can do the same conservatively by projecting each equation
to univariate coordinate slices.  The resulting per-coordinate dyadic factors
k_j can move the 076 system homothety inward:

    z = S y,      S_j = S_system,j / k_j^strength.

The default strength is 0, so 078 inherits the 076 homothety unless the caller
explicitly enables polygonal scale pressure with --polygon-k-strength.

078 also adapts the "start optimize point" to the machinery already present in
070: the first Pandrosion corrector anchor is selected from a tiny deterministic
portfolio when the default anchor does not improve the residual enough.  This
reuses pandrosion_h + T2_step; it does not introduce a T3 method or a Lairez
fallback.
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple

Complex = complex

HERE = Path(__file__).resolve().parent
FLOW076_PATH = HERE / "076_pandrosion_homothetic_geometry.py"
DEFAULT_OUTDIR_078 = "/tmp/pandrosion_078_batches"


def load_076():
    spec = importlib.util.spec_from_file_location("flow076_for_078", str(FLOW076_PATH))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {FLOW076_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["flow076_for_078"] = mod
    spec.loader.exec_module(mod)
    return mod


f076 = load_076()
m = f076.m


@dataclass
class PolygonConfig:
    enabled: bool = True
    projections: int = 5
    n_radii: int = 32
    n_phases: int = 8
    k_strength: float = 0.0
    max_k: float = 8.0
    quantile: float = 0.75
    force_k: str = ""
    portfolio_base_passes: int = 1
    portfolio_missing_threshold: int = 4
    start_optimize: bool = True
    start_strength: float = 1.0
    start_anchor_trials: int = 2
    start_trigger: float = 0.985


LAST_POLYGON_NOTES: List[str] = []
LAST_K_VECTOR: List[float] = []


# -----------------------------------------------------------------------------
# Univariate Pandrosion polygon, adapted from pandrosion_vs_lairez.py.
# -----------------------------------------------------------------------------


def horner1(coeffs: Sequence[Complex], z: Complex) -> Complex:
    p = 0.0 + 0.0j
    for c in coeffs:
        p = p * z + c
    return p


def trim_coeffs(coeffs: Sequence[Complex], eps: float = 1e-14) -> List[Complex]:
    out = list(coeffs)
    while len(out) > 1 and abs(out[0]) <= eps:
        out.pop(0)
    return out


def geometric_root_scale(coeffs: Sequence[Complex]) -> float:
    coeffs = trim_coeffs(coeffs)
    d = len(coeffs) - 1
    if d <= 0:
        return 1.0
    a0 = abs(coeffs[-1])
    ad = abs(coeffs[0])
    if a0 < 1e-300 or ad < 1e-300:
        return 1.0
    return max(1e-12, min(1e12, (a0 / ad) ** (1.0 / d)))


def pandrosion_polygon_clusters(
    coeffs: Sequence[Complex],
    n_radii: int = 32,
    n_phases: int = 8,
) -> List[Tuple[float, int]]:
    coeffs = trim_coeffs(coeffs)
    d = len(coeffs) - 1
    if d <= 0:
        return []

    R_geo = geometric_root_scale(coeffs)
    ad = abs(coeffs[0])
    if ad > 1e-300:
        upper = []
        for k in range(1, len(coeffs)):
            v = abs(coeffs[k]) / ad
            if v > 1e-300:
                upper.append(v ** (1.0 / k))
        rho = 2.0 * max(upper) if upper else 2.0
    else:
        rho = 2.0

    rev = list(reversed(coeffs))
    if abs(rev[0]) > 1e-300:
        lower = []
        for k in range(1, len(rev)):
            v = abs(rev[k]) / abs(rev[0])
            if v > 1e-300:
                lower.append(v ** (1.0 / k))
        rho_lo = 1.0 / (2.0 * max(lower)) if lower else 1.0
    else:
        rho_lo = 1e-3

    n_radii = max(8, int(n_radii))
    n_phases = max(4, int(n_phases))
    r_min = max(min(rho_lo * 0.5, R_geo / 100.0), 1e-12)
    r_max = max(rho * 1.5, R_geo * 100.0, r_min * 1.01)
    log_min, log_max = math.log(r_min), math.log(r_max)
    log_rs = [log_min + (log_max - log_min) * i / (n_radii - 1) for i in range(n_radii)]
    phases = [2.0 * math.pi * k / n_phases for k in range(n_phases)]

    log_avg: List[float] = []
    for log_r in log_rs:
        r = math.exp(log_r)
        total = 0.0
        for phi in phases:
            z = complex(r * math.cos(phi), r * math.sin(phi))
            val = abs(horner1(coeffs, z))
            total += math.log(val) if val > 1e-300 else -700.0
        avg = total / n_phases
        log_avg.append(max(-700.0, min(700.0, avg)) if math.isfinite(avg) else -700.0)

    counts: List[int] = []
    for i in range(n_radii - 1):
        dr = log_rs[i + 1] - log_rs[i]
        slope = (log_avg[i + 1] - log_avg[i]) / dr if dr else 0.0
        if not math.isfinite(slope):
            count = d
        else:
            count = max(0, min(d, int(round(slope))))
        counts.append(count)
    for i in range(1, len(counts)):
        if counts[i] < counts[i - 1]:
            counts[i] = counts[i - 1]

    clusters: List[Tuple[float, int]] = []
    prev = 0
    i = 0
    while i < len(counts):
        if counts[i] > prev:
            j = i
            while j + 1 < len(counts) and counts[j + 1] == counts[i]:
                j += 1
            multiplicity = counts[i] - prev
            radius = math.exp((log_rs[i] + log_rs[j + 1]) / 2.0)
            clusters.append((radius, multiplicity))
            prev = counts[i]
            i = j + 1
        else:
            i += 1
    total = sum(mult for _, mult in clusters)
    if total < d:
        clusters.append((math.exp(log_rs[-1]), d - total))
    return clusters


def k_optimal_dyadic(
    coeffs: Sequence[Complex],
    n_radii: int = 32,
    n_phases: int = 8,
    max_k: float = 8.0,
) -> float:
    clusters = list(pandrosion_polygon_clusters(coeffs, n_radii=n_radii, n_phases=n_phases)) + [(1.0, 1)]
    clusters = [(max(float(r), 1e-300), max(1, int(mul))) for r, mul in clusters if math.isfinite(r)]
    if not clusters:
        return 1.0
    clusters.sort(key=lambda item: item[0])
    total = max(1, sum(mul for _, mul in clusters))
    lo_target = 0.10 * total
    hi_target = 0.90 * total
    cumul = 0
    r_p10 = clusters[0][0]
    r_p90 = clusters[-1][0]
    for radius, mul in clusters:
        cumul += mul
        if cumul >= lo_target:
            r_p10 = radius
            break
    cumul = 0
    for radius, mul in clusters:
        cumul += mul
        if cumul >= hi_target:
            r_p90 = radius
            break
    k_raw = max(1.0, r_p90 / max(r_p10, 1e-300))
    k = 2.0 ** round(math.log2(k_raw)) if k_raw > 0 else 1.0
    return max(1.0, min(float(max_k), float(k)))


# -----------------------------------------------------------------------------
# dD projection adapter.
# -----------------------------------------------------------------------------


def parse_force_k(value: str, n: int) -> List[float] | None:
    if not value:
        return None
    vals = [float(x.strip()) for x in value.replace(";", ",").split(",") if x.strip()]
    if len(vals) == 1 and n > 1:
        vals *= n
    if len(vals) != n:
        raise ValueError(f"--polygon-force-k expected 1 or {n} values, got {len(vals)}")
    return [max(1.0, v) for v in vals]


def projection_base(n: int, projection_index: int) -> List[Complex]:
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    out: List[Complex] = []
    for j in range(n):
        theta = 2.0 * math.pi * (((projection_index + 1) * (j + 1) * phi) % 1.0)
        out.append(complex(math.cos(theta), math.sin(theta)))
    return out


def coordinate_projection_coeffs(poly, coord: int, base: Sequence[Complex]) -> List[Complex]:
    max_power = max((alpha[coord] for alpha in poly.keys()), default=0)
    coeffs = [0.0 + 0.0j for _ in range(max_power + 1)]
    for alpha, coeff in poly.items():
        term = complex(coeff)
        for j, power in enumerate(alpha):
            if j == coord:
                continue
            if power:
                term *= base[j] ** power
        coeffs[alpha[coord]] += term
    return trim_coeffs([coeffs[p] for p in range(max_power, -1, -1)])


def dyadic_k_vector(polys, cfg: PolygonConfig) -> List[float]:
    n = len(polys)
    forced = parse_force_k(cfg.force_k, n)
    if forced is not None:
        return [min(float(cfg.max_k), k) for k in forced]

    out: List[float] = []
    for coord in range(n):
        ks: List[float] = []
        for pidx in range(max(1, int(cfg.projections))):
            base = projection_base(n, pidx)
            for poly in polys:
                coeffs = coordinate_projection_coeffs(poly, coord, base)
                if len(coeffs) < 3:
                    continue
                k = k_optimal_dyadic(
                    coeffs,
                    n_radii=cfg.n_radii,
                    n_phases=cfg.n_phases,
                    max_k=cfg.max_k,
                )
                if math.isfinite(k):
                    ks.append(k)
        if not ks:
            out.append(1.0)
            continue
        logs = sorted(math.log(max(1.0, k)) for k in ks)
        q = max(0.0, min(1.0, float(cfg.quantile)))
        idx = min(len(logs) - 1, int(round(q * (len(logs) - 1))))
        kval = 2.0 ** round(logs[idx] / math.log(2.0))
        out.append(max(1.0, min(float(cfg.max_k), float(kval))))
    return out


def start_anchor_candidates(
    z: Sequence[Complex],
    epoch: int,
    cfg: PolygonConfig,
) -> List[List[Complex]]:
    """Extra near-iterate anchors for the existing 070 Pandrosion corrector."""
    if not cfg.start_optimize or cfg.start_anchor_trials <= 0:
        return []
    n = len(z)
    radius0 = max(1.0, max(abs(v) for v in z))
    radius0 *= max(0.0, float(cfg.start_strength)) / (1.0 + 0.13 * max(0, epoch))
    if radius0 <= 0.0:
        return []

    anchors: List[List[Complex]] = []
    golden = 2.399963229728653
    strengths = (0.028, 0.050, 0.090, 0.140, 0.210)
    for trial in range(max(0, int(cfg.start_anchor_trials))):
        r = strengths[trial % len(strengths)] * radius0
        vec: List[Complex] = []
        for j, zj in enumerate(z):
            theta = golden * (trial + 1) * (j + 1) + 0.37 * epoch
            skew = 1.0 + 0.11 * ((trial + j) % 3)
            vec.append(
                complex(zj)
                + r * skew * complex(math.cos(theta), math.sin(theta))
                + complex(0.0017 * (j + 1), -0.0013 * (trial + 1))
            )
        anchors.append(vec)
    return anchors


def accept_from_anchor(
    polys,
    z: Sequence[Complex],
    anchor: Sequence[Complex],
    rz: float,
    best_r: float,
    best_z: List[Complex],
    modes: Sequence[str],
    quad_cap: int,
    deadline: float | None,
) -> Tuple[float, List[Complex]]:
    if m.timed_out(deadline) or not m.safe_z(anchor):
        return best_r, best_z
    Fa = m.F_eval(polys, anchor)
    if not m.finite_vec(Fa):
        return best_r, best_z
    for mode in modes:
        if m.timed_out(deadline):
            break
        local_before = best_r
        cand, ok = m.pandrosion_h(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
        if ok:
            best_r, best_z = m.accept_damped(polys, z, cand, best_r, best_z)
        if best_r < 0.92 * local_before:
            cand, ok = m.T2_step(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
            if ok:
                best_r, best_z = m.accept_damped(polys, z, cand, best_r, best_z)
        if best_r < 0.28 * rz:
            break
    return best_r, best_z


def install_start_anchor_optimizer(cfg: PolygonConfig) -> None:
    """Patch 070.corrector with a small anchor-start optimizer, no new method."""
    if not cfg.start_optimize:
        return

    def corrector_078(polys, z_init: Sequence[Complex], tol: float = 1e-9,
                      max_epochs: int = 12, quad_cap: int = 16,
                      modes: Sequence[str] = ("system", "integral_gl", "blend"),
                      rescue_modes: Sequence[str] = (),
                      deadline: float | None = None):
        z = list(z_init)
        epochs_used = 0
        anchor = m.deterministic_anchor(z, 0)
        for epoch in range(max_epochs):
            if m.timed_out(deadline):
                return z, False, epochs_used
            epochs_used = epoch + 1
            rz = m.residual_norm(polys, z)
            if rz < tol:
                return z, True, epoch
            if not math.isfinite(rz):
                return z, False, epoch
            Fa = m.F_eval(polys, anchor)
            if m.finite_vec(Fa) and max(abs(v) for v in Fa) < tol:
                return list(anchor), True, epoch
            best_r, best_z = rz, list(z)

            tried = set()
            for mode in modes:
                if m.timed_out(deadline):
                    break
                tried.add(mode)
                local_before = best_r
                cand, ok = m.pandrosion_h(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
                if ok:
                    best_r, best_z = m.accept_damped(polys, z, cand, best_r, best_z)
                if best_r < 0.92 * local_before:
                    cand, ok = m.T2_step(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
                    if ok:
                        best_r, best_z = m.accept_damped(polys, z, cand, best_r, best_z)
                if best_r < 0.22 * rz:
                    break

            trigger = max(0.05, min(0.999, float(cfg.start_trigger)))
            if epoch == 0 and best_r > trigger * rz and not m.timed_out(deadline):
                probe_modes = tuple(modes[:1]) if len(modes) > 1 else tuple(modes)
                for alt_anchor in start_anchor_candidates(z, epoch, cfg):
                    best_r, best_z = accept_from_anchor(
                        polys, z, alt_anchor, rz, best_r, best_z,
                        probe_modes, quad_cap, deadline,
                    )
                    if best_r < 0.30 * rz:
                        break

            if best_r > 0.60 * rz:
                for mode in rescue_modes:
                    if mode in tried or m.timed_out(deadline):
                        continue
                    local_before = best_r
                    cand, ok = m.pandrosion_h(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
                    if ok:
                        best_r, best_z = m.accept_damped(polys, z, cand, best_r, best_z)
                    if best_r < 0.85 * local_before:
                        cand, ok = m.T2_step(polys, anchor, z, Fa, mode, quad_cap=quad_cap)
                        if ok:
                            best_r, best_z = m.accept_damped(polys, z, cand, best_r, best_z)
                    if best_r < 0.38 * rz:
                        break

            if best_r > 0.985 * rz and not m.timed_out(deadline):
                cand, ok = m.m068.armijo_fallback(polys, z)
                if ok:
                    best_r, best_z = m.accept_damped(polys, z, cand, best_r, best_z)

            if best_r >= rz:
                anchor = m.deterministic_anchor(z, epoch + 1, strength=0.050)
                continue
            z = list(best_z)
            anchor = list(z)
        return z, m.residual_norm(polys, z) < tol, epochs_used

    m.corrector = corrector_078


def install_polygon_scaling(cfg: PolygonConfig) -> None:
    original_scale_fn = f076.system_diagonal_scales
    original_write_outputs = f076.write_outputs
    original_run_base_batches = f076.run_base_batches

    def polygon_system_scales(polys, min_scale: float = 0.25, max_scale: float = 4.0, strength: float = 1.0):
        base = original_scale_fn(polys, min_scale=min_scale, max_scale=max_scale, strength=strength)
        LAST_POLYGON_NOTES.clear()
        if not cfg.enabled:
            LAST_POLYGON_NOTES.append("078 polygon scaling disabled")
            return base
        ks = dyadic_k_vector(polys, cfg)
        LAST_K_VECTOR[:] = list(ks)
        lo, hi = math.log(min_scale), math.log(max_scale)
        adjusted: List[float] = []
        for scale, kval in zip(base, ks):
            x = math.log(max(1e-12, float(scale))) - float(cfg.k_strength) * math.log(max(1.0, kval))
            adjusted.append(math.exp(max(lo, min(hi, x))))
        note = (
            "078 k_optimal_dyadic projections="
            f"{cfg.projections}, radii={cfg.n_radii}, phases={cfg.n_phases}, "
            f"k={f076.encode_scales(ks)}, strength={cfg.k_strength}, "
            f"base={f076.encode_scales(base)}, adjusted={f076.encode_scales(adjusted)}"
        )
        LAST_POLYGON_NOTES.append(note)
        print(note, flush=True)
        return adjusted

    def write_outputs_078(summary_rows, batch_rows, roots, args):
        start_note = (
            "078 start-anchor optimize "
            f"trials={cfg.start_anchor_trials}, trigger={cfg.start_trigger}, "
            f"strength={cfg.start_strength}"
            if cfg.start_optimize
            else "078 start-anchor optimize disabled"
        )
        prefix = "; ".join([*LAST_POLYGON_NOTES, start_note])
        for row in summary_rows:
            if row.alg == "076-homothetic-system-geometry":
                row.alg = "078-pandrosion-polygon-homothetic"
                row.notes = f"{prefix}; {row.notes}" if prefix else row.notes
        return original_write_outputs(summary_rows, batch_rows, roots, args)

    def run_batch_078(args, stage, retry, indices, outdir, roots_before, sep, scales=None):
        outdir.mkdir(parents=True, exist_ok=True)
        Btag = f"{args.family}_{args.case.replace(',', 'x')}_seed{args.seed_index}"
        if len(indices) <= 8:
            itag = "i" + "_".join(str(i) for i in indices)
        else:
            itag = f"{min(indices)}_{max(indices)+1}_n{len(indices)}"
        outfile = outdir / f"078_{Btag}_{stage}_r{retry}_{itag}.json"
        cmd = [
            sys.executable, "-S", str(Path(__file__).resolve()),
            "--start-strength", str(cfg.start_strength),
            "--start-anchor-trials", str(cfg.start_anchor_trials),
            "--start-trigger", str(cfg.start_trigger),
            "--worker", "--family", args.family, "--case", args.case,
            "--seed-index", str(args.seed_index), "--retry", str(retry),
            "--indices", f076.encode_indices(indices), "--out", str(outfile),
            "--tol", str(args.tol), "--max-steps", str(args.max_steps),
            "--max-epochs", str(args.max_epochs), "--quad-cap", str(args.quad_cap),
            "--dt0", str(args.dt0), "--dtmax", str(args.dtmax),
            "--scales", f076.encode_scales(scales if scales is not None else getattr(args, "active_scales", [])),
            "--modes", str(getattr(args, "modes", ",".join(f076.DEFAULT_MODES))),
            "--rescue-modes", str(getattr(args, "rescue_modes", ",".join(f076.DEFAULT_RESCUE_MODES))),
        ]
        if not cfg.start_optimize:
            cmd.insert(3, "--no-start-optimize")
        if args.equation_normalize:
            cmd.append("--equation-normalize")
        if args.verbose:
            cmd.append("--verbose")
        t0 = time.time()
        status = "ok"
        try:
            f076.subprocess.run(
                cmd,
                timeout=f076.batch_timeout_for(args),
                check=True,
                stdout=(None if args.verbose else f076.subprocess.DEVNULL),
                stderr=(None if args.verbose else f076.subprocess.DEVNULL),
            )
        except f076.subprocess.TimeoutExpired:
            status = "timeout"
        except f076.subprocess.CalledProcessError as exc:
            status = f"error:{exc.returncode}"
        seconds = time.time() - t0
        candidates_roots, paths, candidates, rows = f076.read_batch_candidates(str(outfile))
        if outfile.exists():
            try:
                data = f076.json.loads(outfile.read_text())
                seconds = float(data.get("seconds", seconds))
            except Exception:
                pass
        after_roots, new = f076.count_new_roots(roots_before, candidates_roots, sep=sep)
        br = f076.BatchRow(
            stage=stage,
            retry=retry,
            indices=f076.encode_indices(indices),
            paths=paths if paths else len(indices),
            candidates=candidates,
            roots_before=len(roots_before),
            roots_after=len(after_roots),
            new_roots=new,
            seconds=float(seconds),
            status=status,
            path_json=str(outfile),
        )
        return br, rows

    def run_base_batches_078(args, B: int, outdir: Path, base_scales: Sequence[float]):
        base_rows, trace_rows = original_run_base_batches(args, B, outdir, base_scales)
        roots, _, _ = f076.read_all_roots(base_rows)
        roots = f076.cluster_roots(roots, sep=args.cluster_sep)
        missing = B - len(roots)
        if cfg.portfolio_base_passes <= 1 or missing <= 0:
            return base_rows, trace_rows
        if missing < max(1, int(cfg.portfolio_missing_threshold)):
            return base_rows, trace_rows
        if getattr(args, "_deadline", None) is not None and time.time() >= args._deadline:
            return base_rows, trace_rows

        chunks = f076.chunk_ranges(B, args.base_chunk_size)
        max_workers = max(1, int(args.parallel_base))
        extra_rows = []
        extra_trace = []
        for retry in range(1, max(1, int(cfg.portfolio_base_passes))):
            if getattr(args, "_deadline", None) is not None and time.time() >= args._deadline:
                break
            print(
                f"078 portfolio-base retry={retry} missing={missing}/{B} "
                f"chunks={len(chunks)} workers={max_workers}",
                flush=True,
            )
            with f076.concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
                future_items = []
                for pos, inds in enumerate(chunks):
                    fut = ex.submit(
                        f076.run_batch,
                        args,
                        "portfolio-base",
                        retry,
                        inds,
                        outdir,
                        [],
                        args.cluster_sep,
                        base_scales,
                    )
                    future_items.append((pos, inds, fut))
                raw = []
                for pos, inds, fut in future_items:
                    br, rows = fut.result()
                    raw.append((pos, inds, br, rows))
            raw.sort(key=lambda item: item[0])
            retry_rows = [item[2] for item in raw]
            fixed_rows, roots_after, _, _ = f076.recompute_batch_progress(
                [*base_rows, *extra_rows, *retry_rows],
                args.cluster_sep,
            )
            fixed_retry_rows = fixed_rows[-len(retry_rows):]
            for fixed, (_, inds, _br, rows) in zip(fixed_retry_rows, raw):
                extra_trace.extend(rows)
                print(
                    f"portfolio-base retry={retry} indices={inds[0]}:{inds[-1]+1} "
                    f"candidates={fixed.candidates}/{fixed.paths} roots={fixed.roots_after}/{B} "
                    f"sec={fixed.seconds:.2f}",
                    flush=True,
                )
            extra_rows.extend(fixed_retry_rows)
            roots = roots_after
            missing = B - len(roots)
            if args.stop_at_bezout and len(roots) >= B:
                break

        if extra_rows:
            all_fixed, _, _, _ = f076.recompute_batch_progress([*base_rows, *extra_rows], args.cluster_sep)
            return all_fixed, [*trace_rows, *extra_trace]
        return base_rows, trace_rows

    f076.system_diagonal_scales = polygon_system_scales
    f076.write_outputs = write_outputs_078
    f076.run_batch = run_batch_078
    f076.run_base_batches = run_base_batches_078


def has_option(argv: Sequence[str], name: str) -> bool:
    return any(arg == name or arg.startswith(name + "=") for arg in argv)


def add_078_defaults(argv: List[str], use_defaults: bool) -> List[str]:
    if not use_defaults:
        return argv
    out = list(argv)
    if not has_option(out, "--outdir"):
        out.extend(["--outdir", DEFAULT_OUTDIR_078])
    if not has_option(out, "--time-budget"):
        out.extend(["--time-budget", "20"])
    if not has_option(out, "--batch-timeout"):
        out.extend(["--batch-timeout", "20"])
    if not has_option(out, "--base-chunk-size"):
        out.extend(["--base-chunk-size", "0"])
    if not has_option(out, "--parallel-base"):
        out.extend(["--parallel-base", "0"])
    if not has_option(out, "--homothety") and not has_option(out, "--fast"):
        out.extend(["--homothety", "system"])
    if not has_option(out, "--equation-normalize") and not has_option(out, "--fast"):
        out.append("--equation-normalize")
    if (
        not has_option(out, "--window-fallback")
        and not has_option(out, "--no-window-fallback")
        and not has_option(out, "--fast")
    ):
        out.append("--no-window-fallback")
    return out


def main() -> None:
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--no-polygon-scaling", dest="polygon_scaling", action="store_false", default=True)
    pre.add_argument("--polygon-projections", type=int, default=5)
    pre.add_argument("--polygon-radii", type=int, default=32)
    pre.add_argument("--polygon-phases", type=int, default=8)
    pre.add_argument("--polygon-k-strength", type=float, default=0.0)
    pre.add_argument("--polygon-max-k", type=float, default=8.0)
    pre.add_argument("--polygon-quantile", type=float, default=0.75)
    pre.add_argument("--polygon-force-k", default="")
    pre.add_argument("--portfolio-base-passes", type=int, default=1)
    pre.add_argument("--portfolio-missing-threshold", type=int, default=4)
    pre.add_argument("--start-optimize", dest="start_optimize", action="store_true", default=True)
    pre.add_argument("--no-start-optimize", dest="start_optimize", action="store_false")
    pre.add_argument("--start-strength", type=float, default=1.0)
    pre.add_argument("--start-anchor-trials", type=int, default=2)
    pre.add_argument("--start-trigger", type=float, default=0.985)
    pre.add_argument("--no-078-defaults", dest="defaults_078", action="store_false", default=True)
    cfg_args, rest = pre.parse_known_args()

    cfg = PolygonConfig(
        enabled=bool(cfg_args.polygon_scaling),
        projections=max(1, int(cfg_args.polygon_projections)),
        n_radii=max(8, int(cfg_args.polygon_radii)),
        n_phases=max(4, int(cfg_args.polygon_phases)),
        k_strength=max(0.0, float(cfg_args.polygon_k_strength)),
        max_k=max(1.0, float(cfg_args.polygon_max_k)),
        quantile=max(0.0, min(1.0, float(cfg_args.polygon_quantile))),
        force_k=str(cfg_args.polygon_force_k or ""),
        portfolio_base_passes=max(1, int(cfg_args.portfolio_base_passes)),
        portfolio_missing_threshold=max(1, int(cfg_args.portfolio_missing_threshold)),
        start_optimize=bool(cfg_args.start_optimize),
        start_strength=max(0.0, float(cfg_args.start_strength)),
        start_anchor_trials=max(0, int(cfg_args.start_anchor_trials)),
        start_trigger=max(0.05, min(0.999, float(cfg_args.start_trigger))),
    )
    install_start_anchor_optimizer(cfg)
    install_polygon_scaling(cfg)
    sys.argv = [sys.argv[0], *add_078_defaults(rest, bool(cfg_args.defaults_078))]
    f076.main()


if __name__ == "__main__":
    main()
