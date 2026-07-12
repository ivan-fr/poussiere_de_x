"""Pandrosion 321: compact standalone NumPy root-harvesting engine.

Core geometry: logarithmic atlas starts, vectorised swarm, paired local jets,
immediate Broyden updates, overflow-safe root deflation, adaptive LM damping,
and direct/reciprocal IRP rescue.  The base system is always queried as a
black-box oracle.  Exact expression polynomials and dense Kostlan systems are
supported; large Kostlan cases use an explicitly labelled random-feature
surrogate.  No local project imports are used.
"""
from __future__ import annotations

import argparse
import ast
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np


# ---------- small utilities -------------------------------------------------

def clock() -> float:
    return time.perf_counter()


def norm(v: Any) -> float:
    try:
        x = float(np.linalg.norm(v))
        return x if math.isfinite(x) else float("inf")
    except Exception:
        return float("inf")


def norms(F: Any) -> np.ndarray:
    r = np.asarray(np.linalg.norm(np.asarray(F, np.complex128), axis=1), float)
    r[~np.isfinite(r)] = np.inf
    return r


def strict_json(x: Any) -> Any:
    if isinstance(x, dict):
        return {str(k): strict_json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [strict_json(v) for v in x]
    if isinstance(x, np.ndarray):
        return strict_json(x.tolist())
    if isinstance(x, np.generic):
        return strict_json(x.item())
    if isinstance(x, complex):
        return [float(x.real), float(x.imag)]
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def parse_case(raw: str) -> tuple[int, int]:
    p = [s.strip() for s in str(raw).lower().replace("x", ",").replace(":", ",").split(",") if s.strip()]
    if len(p) != 2 or int(p[0]) <= 0 or int(p[1]) <= 0:
        raise ValueError(f"invalid case {raw!r}; expected positive n,d")
    return int(p[0]), int(p[1])


def floats(raw: Optional[str], default: Sequence[float]) -> list[float]:
    if raw is None or not str(raw).strip():
        return list(default)
    out = []
    for token in str(raw).replace(";", ",").split(","):
        try:
            x = float(token)
        except ValueError:
            continue
        if math.isfinite(x) and x > 0:
            out.append(x)
    return out or list(default)


def starts_from_text(raw: Optional[str], n: int) -> list[np.ndarray]:
    if not raw or not str(raw).strip():
        return []
    text = str(raw).replace("i", "j")
    if n == 1 and ";" not in text and "|" not in text:
        return [np.asarray([complex(x)], np.complex128) for x in text.split(",") if x.strip()]
    out = []
    for point in text.replace("|", ";").split(";"):
        q = [x.strip() for x in point.split(",") if x.strip()]
        if not q:
            continue
        if len(q) != n:
            raise ValueError(f"start {point!r} has {len(q)} coordinates; expected {n}")
        out.append(np.asarray([complex(x) for x in q], np.complex128))
    return out


def mix64(x: int) -> int:
    x = (int(x) + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return (x ^ (x >> 31)) & 0xFFFFFFFFFFFFFFFF


def stable_seed(n: int, d: int, index: int) -> int:
    return int(mix64(0x50414E44524F5349 + 1000003 * n + 9176 * d + 97 * index) & 0x7FFFFFFF)


def direction(n: int, index: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(mix64(seed + 104729 * index))
    v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    return np.asarray(v / max(norm(v), 1e-300) * math.sqrt(n), np.complex128)


# ---------- black-box oracles ----------------------------------------------

class PolyExpr:
    allowed = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name,
               ast.Add, ast.Sub, ast.Mult, ast.Pow, ast.USub, ast.UAdd, ast.Load)

    def __init__(self, raw: str) -> None:
        self.raw = raw.strip()
        self.tree = ast.parse(self.raw.replace("^", "**"), mode="eval")
        for node in ast.walk(self.tree):
            if not isinstance(node, self.allowed):
                raise ValueError(f"unsupported polynomial syntax: {type(node).__name__}")
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
                if not isinstance(node.right, ast.Constant):
                    raise ValueError("polynomial exponents must be constants")
                e = node.right.value
                if isinstance(e, complex) or float(e) != int(e) or int(e) < 0:
                    raise ValueError("polynomial exponents must be non-negative integers")

    def eval(self, env: dict[str, Any]) -> Any:
        def go(node: ast.AST) -> Any:
            if isinstance(node, ast.Constant):
                return node.value
            if isinstance(node, ast.Name):
                if node.id not in env:
                    raise ValueError(f"unknown variable {node.id!r}")
                return env[node.id]
            if isinstance(node, ast.UnaryOp):
                v = go(node.operand)
                return -v if isinstance(node.op, ast.USub) else v
            if isinstance(node, ast.BinOp):
                a, b = go(node.left), go(node.right)
                if isinstance(node.op, ast.Add): return a + b
                if isinstance(node.op, ast.Sub): return a - b
                if isinstance(node.op, ast.Mult): return a * b
                if isinstance(node.op, ast.Pow): return a ** int(b)
            raise ValueError("invalid polynomial expression")
        return go(self.tree.body)


@dataclass
class Oracle:
    n: int
    d: int
    seed: int
    kind: str
    eval_count: int = 0
    seconds: float = 0.0

    def _eval_batch(self, Z: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def eval_batch(self, Z: Any) -> np.ndarray:
        t0 = clock()
        zz = np.asarray(Z, np.complex128)
        if zz.ndim == 1:
            zz = zz[None, :]
        with np.errstate(all="ignore"):
            out = np.asarray(self._eval_batch(zz), np.complex128)
        out[~np.isfinite(out)] = complex(1e300, 0)
        self.eval_count += len(zz)
        self.seconds += clock() - t0
        return out

    def eval(self, z: Any) -> np.ndarray:
        return self.eval_batch(np.asarray(z, np.complex128)[None, :])[0]


class ExpressionOracle(Oracle):
    def __init__(self, n: int, d: int, seed: int, raw: str, names: Optional[str]) -> None:
        parts = [x.strip() for x in raw.replace("|", ";").split(";") if x.strip()]
        variables = [x.strip() for x in str(names or "").replace(";", ",").split(",") if x.strip()]
        variables = variables or (["x"] if n == 1 else [f"x{i+1}" for i in range(n)])
        if len(parts) != n or len(variables) != n:
            raise ValueError(f"expected {n} polynomials and {n} variable names")
        super().__init__(n, d, seed, "polynomial-exact")
        self.expressions, self.names = [PolyExpr(x) for x in parts], variables

    def _eval_batch(self, Z: np.ndarray) -> np.ndarray:
        env: dict[str, Any] = {"I": 1j, "j": 1j}
        for i, name in enumerate(self.names): env[name] = Z[:, i]
        for i in range(self.n):
            env.setdefault(f"x{i+1}", Z[:, i]); env.setdefault(f"z{i+1}", Z[:, i])
        if self.n == 1:
            env.setdefault("x", Z[:, 0]); env.setdefault("z", Z[:, 0])
        cols = []
        for expr in self.expressions:
            v = np.asarray(expr.eval(env), np.complex128)
            cols.append(np.full(len(Z), complex(v), np.complex128) if v.ndim == 0 else v.reshape(len(Z)))
        return np.stack(cols, axis=1)


def compositions(d: int, n: int) -> np.ndarray:
    out: list[tuple[int, ...]] = []
    def rec(pos: int, rem: int, cur: list[int]) -> None:
        if pos == n:
            out.append(tuple(cur)); return
        for k in range(rem + 1):
            rec(pos + 1, rem - k, cur + [k])
    rec(0, d, [])
    return np.asarray(out, np.int16 if d < 32767 else np.int32)


class KostlanOracle(Oracle):
    def __init__(self, n: int, d: int, seed: int, dense_cap: int, features: int, normalize: bool) -> None:
        terms = math.comb(n + d, d)
        rng = np.random.default_rng(seed)
        if terms <= dense_cap:
            exps = compositions(d, n)
            totals = exps.sum(1)
            logs = np.asarray([math.lgamma(d + 1) - math.lgamma(d - int(t) + 1)
                               - sum(math.lgamma(int(a) + 1) for a in row)
                               for row, t in zip(exps, totals)])
            log_scales = 0.5*logs; scales = np.exp(log_scales - float(np.max(log_scales)))
            kind = "kostlan-dense"
        else:
            m = max(n + 8, features)
            draws = rng.multinomial(d, np.full(n + 1, 1 / (n + 1)), size=m)
            exps = draws[:, 1:]
            scales = np.ones(m)
            kind = "kostlan-random-feature-surrogate"
        coeff = (rng.standard_normal((n, len(exps))) + 1j * rng.standard_normal((n, len(exps)))) / math.sqrt(2)
        coeff *= scales[None, :]
        if normalize:
            coeff /= np.maximum(np.linalg.norm(coeff, axis=1), 1e-300)[:, None]
        super().__init__(n, d, seed, kind)
        self.exps, self.degrees, self.coeff = exps, exps.sum(1), np.asarray(coeff, np.complex128)

    def _eval_batch(self, Z: np.ndarray) -> np.ndarray:
        scale = np.maximum(1, np.max(np.abs(Z), axis=1)); W = Z/scale[:, None]
        mon = np.ones((len(Z), len(self.exps)), np.complex128)
        for j in range(self.n):
            powers = np.ones((len(Z), self.d + 1), np.complex128)
            for k in range(1, self.d + 1): powers[:, k] = powers[:, k - 1] * W[:, j]
            mon *= powers[:, self.exps[:, j]]
        mon *= np.exp(np.clip(np.log(scale)[:, None]*(self.degrees[None, :]-self.d), -745, 0))
        return mon @ self.coeff.T


def make_oracle(args: argparse.Namespace, n: int, d: int) -> Oracle:
    seed = stable_seed(n, d, args.seed_index)
    if args.system_source in {"poly", "polynomial"}:
        if not args.polys:
            raise ValueError("--polys is required for polynomial systems")
        return ExpressionOracle(n, d, seed, args.polys, args.variables)
    return KostlanOracle(n, d, seed, args.dense_max_terms, args.features, args.equation_normalize)


# ---------- targets, jets and linear algebra --------------------------------

class Target:
    def __init__(self, oracle: Oracle) -> None:
        self.oracle = oracle

    def to_base(self, y: Any) -> np.ndarray:
        return np.asarray(y, np.complex128)

    def from_base(self, y: Any) -> np.ndarray:
        return np.asarray(y, np.complex128)

    def eval_batch(self, Y: Any) -> np.ndarray:
        return self.oracle.eval_batch(self.to_base(Y))

    def eval(self, y: Any) -> np.ndarray:
        return self.eval_batch(np.asarray(y)[None, :])[0]


class IRPTarget(Target):
    def __init__(self, base: Target, scale: complex, reciprocal: bool) -> None:
        self.base, self.oracle = base, base.oracle
        self.scale, self.reciprocal = complex(scale), bool(reciprocal)

    def to_base(self, u: Any) -> np.ndarray:
        x = np.asarray(u, np.complex128)
        y = self.scale * x if not self.reciprocal else 1 / np.where(np.abs(self.scale * x) < 1e-14, self.scale * x + 1e-14, self.scale * x)
        return self.base.to_base(y)

    def from_base(self, y: Any) -> np.ndarray:
        x = self.base.from_base(y)
        return x / self.scale if not self.reciprocal else 1 / np.where(np.abs(self.scale * x) < 1e-14, self.scale * x + 1e-14, self.scale * x)


def paired_jets(target: Target, Y: Any, F0: Optional[Any], hrel: float) -> tuple[np.ndarray, np.ndarray, int]:
    yy = np.asarray(Y, np.complex128)
    if yy.ndim == 1: yy = yy[None, :]
    b, n = yy.shape
    f0 = target.eval_batch(yy) if F0 is None else np.asarray(F0, np.complex128)
    h = hrel * np.maximum(1, np.abs(yy))
    plus = np.repeat(yy[:, None, :], n, axis=1)
    minus = plus.copy()
    idx = np.arange(n)
    plus[:, idx, idx] += h; minus[:, idx, idx] -= h
    f = target.eval_batch(np.vstack((plus.reshape(b*n, n), minus.reshape(b*n, n))))
    fp, fm = f[:b*n].reshape(b, n, n), f[b*n:].reshape(b, n, n)
    J = np.transpose(fp - fm, (0, 2, 1)) / (2 * h[:, None, :])
    return f0, J, (2*n*b + (b if F0 is None else 0))


def linear_step(J: np.ndarray, f: np.ndarray, mu: float, trust: float, y: np.ndarray) -> tuple[np.ndarray, str]:
    try:
        if mu > 0:
            h = J.conj().T
            scale = float(np.trace(h @ J).real) / max(1, len(y)) + 1e-300
            delta = np.linalg.solve(h @ J + mu * scale * np.eye(len(y)), h @ (-f))
            method = "lm"
        else:
            delta = np.linalg.solve(J, -f); method = "solve"
    except Exception:
        try:
            delta = np.linalg.pinv(J, rcond=1e-12) @ (-f); method = "pinv"
        except Exception:
            delta = np.zeros_like(y); method = "failed"
    limit = trust * max(1, norm(y)) if trust > 0 else 10 * max(1, norm(y))
    nd = norm(delta)
    if math.isfinite(nd) and nd > limit: delta *= limit / nd
    return np.asarray(delta, np.complex128), method


def deflation_logs(Y: Any, roots: Sequence[np.ndarray], alpha: float) -> np.ndarray:
    yy = np.asarray(Y, np.complex128)
    if yy.ndim == 1: yy = yy[None, :]
    out = np.zeros(len(yy))
    for root in roots:
        out += np.log1p(alpha / np.maximum(np.linalg.norm(yy - root[None, :], axis=1), 1e-12))
    return out


def line_choice(R: Any, Y: Any, y: Any, current: float, best: float,
                roots: Sequence[np.ndarray], alpha: float) -> tuple[Optional[int], np.ndarray]:
    rr = np.asarray(R, float)
    logw = deflation_logs(Y, roots, alpha) if roots else np.zeros(len(rr))
    logw0 = float(deflation_logs(np.asarray(y)[None, :], roots, alpha)[0]) if roots else 0
    merit = np.log(np.maximum(rr, 1e-300)) + logw
    admissible = np.isfinite(rr) & ((merit < math.log(max(current, 1e-300)) + logw0) | (rr < best))
    return (int(np.argmin(np.where(admissible, merit, np.inf))), merit) if np.any(admissible) else (None, merit)


# ---------- compact Pandrosion corrector ------------------------------------

@dataclass
class CorrectResult:
    y: np.ndarray
    f: np.ndarray
    residual: float
    ok: bool
    status: str
    epochs: int
    rebuilds: int
    broyden: int
    line_evals: int
    jet_evals: int
    parabolic_evals: int
    mu: float


def correct(target: Target, y0: Any, args: argparse.Namespace, known: Sequence[np.ndarray],
            epochs: int, deadline: Optional[float]) -> CorrectResult:
    y = np.asarray(y0, np.complex128).copy(); f = target.eval(y); r = norm(f)
    best_y, best_f, best_r = y.copy(), f.copy(), r
    A: Optional[np.ndarray] = None; age = rebuilds = updates = line_evals = jet_evals = para = done = 0
    mu = max(0, args.lm_damping); status = "max-epochs"
    full_line = np.asarray(floats(args.line_grid, [1, .75, .5, .25, .125, .0625]))
    short_line = full_line[:min(4, len(full_line))]
    for ep in range(max(1, epochs)):
        if deadline is not None and clock() >= deadline: status = "timeout"; break
        if best_r < args.accept: status = "converged"; break
        force_rebuild = A is None or not args.broyden or age >= args.jet_refresh
        moved = False
        for attempt in range(3):
            rebuilt = force_rebuild or A is None
            if rebuilt:
                _, batch_J, used = paired_jets(target, y[None, :], f[None, :], args.jet_radius)
                A = batch_J[0]; age = 0; rebuilds += 1; jet_evals += used; force_rebuild = False
            delta, _ = linear_step(A, f, mu if args.adaptive_lm else args.lm_damping, args.trust_radius, y)
            if not np.all(np.isfinite(delta)) or norm(delta) <= 1e-15:
                if args.adaptive_lm and attempt < 2:
                    mu = min(1, max(1e-4, mu*5)); continue
                status = "invalid-step"; break
            L = full_line if rebuilt else short_line
            Y = y[None, :] + L[:, None] * delta[None, :]
            F = target.eval_batch(Y); R = norms(F); line_evals += len(Y)
            i, merit = line_choice(R, Y, y, r, best_r, known, args.deflation_alpha)
            if i is None:
                if not rebuilt:
                    force_rebuild = True; continue
                if args.adaptive_lm and attempt < 2:
                    mu = min(1, max(1e-4, mu * 5)); continue
                status = "no-decrease"; break
            lam, yn, fn, rn = float(L[i]), Y[i].copy(), F[i].copy(), float(R[i])
            if args.parabolic and 0 < i < len(L) - 1:
                x1, x2, x3 = float(L[i-1]), float(L[i]), float(L[i+1]); q1, q2, q3 = R[i-1]**2, R[i]**2, R[i+1]**2
                den = (x1-x2)*(x1-x3)*(x2-x3)
                a = (x3*(q2-q1)+x2*(q1-q3)+x1*(q3-q2))/den if abs(den) > 1e-300 else 0
                b = (x3*x3*(q1-q2)+x2*x2*(q3-q1)+x1*x1*(q2-q3))/den if abs(den) > 1e-300 else 0
                star = float(-b/(2*a)) if a > 0 else -1
                if a > 0 and min(x1, x2, x3) < star < max(x1, x2, x3) and abs(star-lam) > .03*max(lam, 1e-12):
                    yp = y + star*delta; fp = target.eval(yp); rp = norm(fp); para += 1; line_evals += 1
                    wp = float(deflation_logs(yp[None, :], known, args.deflation_alpha)[0]) if known else 0
                    if math.isfinite(rp) and math.log(max(rp, 1e-300)) + wp < merit[i]:
                        lam, yn, fn, rn = star, yp, fp, rp
            predicted = norm(f + A @ (lam*delta)) if A is not None else r
            dy, df = yn-y, fn-f; denom = complex(np.vdot(dy, dy))
            if args.broyden and abs(denom) > 1e-300:
                A = A + np.outer(df - A@dy, dy.conj()) / denom; age += 1; updates += 1
            old_r = r; y, f, r = yn, fn, rn; done = ep + 1; moved = True
            if r < best_r: best_y, best_f, best_r = y.copy(), f.copy(), r
            if args.adaptive_lm:
                if r < .5*old_r or r <= predicted: mu *= .25
                elif r > .9*old_r: mu = min(1, max(1e-4, mu*2))
            break
        if not moved: break
    ok = best_r < args.accept
    if ok: status = "converged"
    return CorrectResult(best_y, best_f, best_r, ok, status, done, rebuilds, updates,
                         line_evals, jet_evals, para, mu)


def irp_rescue(base: Target, result: CorrectResult, args: argparse.Namespace,
               known: Sequence[np.ndarray], deadline: Optional[float]) -> CorrectResult:
    if result.ok or not args.irp: return result
    best = result
    gains = floats(args.irp_scales, [1, 2**(1/3), 2**(-1/3), 2, .5])
    candidates = []
    for gain in gains:
        for reciprocal in (False, True):
            chart = IRPTarget(base, complex(gain), reciprocal)
            u = chart.from_base(best.y)
            energy = float(np.mean(np.abs(np.log(np.maximum(np.abs(u), 1e-300)))))
            candidates.append((energy, chart, u))
    for _, chart, u in sorted(candidates, key=lambda x: x[0])[:args.irp_top]:
        if deadline is not None and clock() >= deadline: break
        loc = correct(chart, u, args, [], args.irp_epochs, deadline)
        y = chart.to_base(loc.y); f = loc.f.copy(); r = norm(f)
        if r < best.residual:
            best = CorrectResult(y, f, r, r < args.accept, "irp-converged" if r < args.accept else "irp-improved",
                                 result.epochs + loc.epochs, result.rebuilds + loc.rebuilds,
                                 result.broyden + loc.broyden, result.line_evals + loc.line_evals,
                                 result.jet_evals + loc.jet_evals, result.parabolic_evals + loc.parabolic_evals, loc.mu)
        if best.ok: break
    return best


# ---------- starts, swarm and root finishing --------------------------------

def atlas_start(n: int, trial: int, seed: int) -> np.ndarray:
    radii = [2**(k/3) for k in range(-15, 16)]
    return radii[(trial*17) % len(radii)] * direction(n, trial, seed)


def swarm(base: Target, args: argparse.Namespace, n: int, seed: int) -> tuple[list[np.ndarray], dict[str, Any]]:
    size = args.swarm_size or min(args.pool, max(32, 8*args.count))
    keep = args.swarm_keep or max(8, 3*args.count)
    Y = np.stack([atlas_start(n, i, seed + 31337) for i in range(max(1, size))])
    F = base.eval_batch(Y); R = norms(F); alive = np.isfinite(R); jet_samples = line_samples = 0
    L = np.asarray([1, .5, .25])
    for _ in range(args.swarm_iters):
        ids = np.flatnonzero(alive)
        if not len(ids): break
        _, J, used = paired_jets(base, Y[ids], F[ids], args.jet_radius); jet_samples += used
        try:
            D = np.linalg.solve(J, -F[ids, :, None])[:, :, 0]
        except Exception:
            rows = []
            for j, f in zip(J, F[ids]):
                try: rows.append(np.linalg.pinv(j, rcond=1e-12) @ (-f))
                except Exception: rows.append(np.zeros(n, np.complex128))
            D = np.asarray(rows)
        bad = ~np.all(np.isfinite(D), axis=1); D[bad] = 0
        dn = np.linalg.norm(D, axis=1); cap = 10*np.maximum(1, np.linalg.norm(Y[ids], axis=1))
        D *= np.minimum(1, cap/np.maximum(dn, 1e-300))[:, None]
        cand = Y[ids, None, :] + L[None, :, None]*D[:, None, :]
        fc = base.eval_batch(cand.reshape(-1, n)).reshape(len(ids), len(L), n); line_samples += len(ids)*len(L)
        rc = norms(fc.reshape(-1, n)).reshape(len(ids), len(L)); pick = np.argmin(rc, axis=1); rows = np.arange(len(ids))
        improved = np.isfinite(rc[rows, pick]) & (rc[rows, pick] < R[ids])
        chosen = ids[improved]; Y[chosen] = cand[rows, pick][improved]; F[chosen] = fc[rows, pick][improved]; R[chosen] = rc[rows, pick][improved]
        alive[ids[~improved]] = False
    selected: list[int] = []
    sep = args.swarm_sep or .15*math.sqrt(n)
    for i in np.argsort(R):
        if not math.isfinite(float(R[i])): break
        if all(norm(Y[i]-Y[j]) > sep for j in selected): selected.append(int(i))
        if len(selected) >= keep: break
    return [Y[i].copy() for i in selected], {"size": size, "alive": int(alive.sum()), "kept": len(selected),
        "best_residual": float(R[selected[0]]) if selected else None, "jet_samples": jet_samples, "line_samples": line_samples}


def polish(base: Target, y0: Any, args: argparse.Namespace, deadline: Optional[float]) -> tuple[np.ndarray, np.ndarray, float, dict[str, Any]]:
    y = np.asarray(y0, np.complex128).copy(); f = base.eval(y); r = norm(f); Jfinal = None
    for _ in range(args.polish_steps):
        if r < args.tol or (deadline is not None and clock() >= deadline): break
        _, JJ, _ = paired_jets(base, y[None, :], f[None, :], args.jet_radius); J = JJ[0]; Jfinal = J.copy()
        d, _ = linear_step(J, f, 0, args.trust_radius, y)
        Y = y[None, :] + np.asarray([1, .5, .25, .125])[:, None]*d[None, :]
        F = base.eval_batch(Y); R = norms(F); i = int(np.argmin(R))
        if R[i] >= r: break
        y, f, r, Jfinal = Y[i].copy(), F[i].copy(), float(R[i]), None
    if Jfinal is None and (deadline is None or clock() < deadline):
        _, JJ, _ = paired_jets(base, y[None, :], f[None, :], args.jet_radius); Jfinal = JJ[0]
    if Jfinal is None:
        return y, f, r, {"smin": None, "smax": None, "cond": None, "near_multiple": None, "singular": None, "status": "deadline"}
    try:
        s = np.linalg.svd(Jfinal, compute_uv=False); smin, smax = float(s[-1]), float(s[0])
    except Exception:
        return y, f, r, {"smin": None, "smax": None, "cond": None, "near_multiple": None, "singular": None, "status": "unavailable"}
    cond = None if smin <= 0 else float(smax/smin)
    meta = {"smin": smin, "smax": smax, "cond": cond if cond is not None and math.isfinite(cond) else None,
            "near_multiple": bool(smin <= 1e-8*max(smax, 1e-300)), "singular": bool(smin <= np.finfo(float).eps*max(smax, 1e-300))}
    return y, f, r, meta


# ---------- case driver ------------------------------------------------------

def run_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    started = clock(); n, d = parse_case(case_raw); oracle = make_oracle(args, n, d); base = Target(oracle)
    explicit = starts_from_text(args.starts, n)
    swarm_points, swarm_meta = ([], {"enabled": False})
    if args.swarm and not explicit and args.count:
        swarm_points, swarm_meta = swarm(base, args, n, oracle.seed); swarm_meta["enabled"] = True
    priority = explicit + swarm_points; roots: list[dict[str, Any]] = []; trials = []; duplicates = failures = 0
    for trial in range(args.pool):
        if len(roots) >= args.count: break
        y0 = priority[trial].copy() if trial < len(priority) else atlas_start(n, trial-len(priority), oracle.seed)
        known = [np.asarray(r["z_complex"], np.complex128) for r in roots]
        if known and min(norm(y0-r) for r in known) <= args.early_dup_sep:
            duplicates += 1; trials.append({"trial": trial, "status": "start-near-root", "accepted": False}); continue
        deadline = clock()+args.trial_timeout if args.trial_timeout > 0 else None
        before = oracle.eval_count
        loc = correct(base, y0, args, known, args.epochs, deadline)
        loc = irp_rescue(base, loc, args, known, deadline)
        if loc.residual <= max(args.polish_gate, 100*args.accept):
            z, f, residual, conditioning = polish(base, loc.y, args, deadline)
        else:
            z, f, residual = loc.y.copy(), loc.f.copy(), loc.residual
            conditioning = {"smin": None, "smax": None, "cond": None, "near_multiple": None, "singular": None, "status": "polish-gate"}
        accepted = bool(math.isfinite(residual) and residual < args.accept)
        rec = {"trial": trial, "status": loc.status, "accepted": accepted, "residual": residual,
               "epochs": loc.epochs, "jet_rebuilds": loc.rebuilds, "broyden_updates": loc.broyden,
               "jet_evals": loc.jet_evals, "line_evals": loc.line_evals, "parabolic_evals": loc.parabolic_evals,
               "oracle_samples": oracle.eval_count-before, "conditioning": conditioning}
        if args.verbose_trials: rec.update({"start": y0, "z": z})
        if not accepted: failures += 1; trials.append(rec); continue
        dup = next((i for i, root in enumerate(roots) if norm(z-root["z_complex"]) <= args.cluster_sep), None)
        if dup is not None:
            duplicates += 1; rec.update({"status": "duplicate", "cluster": dup}); trials.append(rec); continue
        root = {"id": len(roots), "trial": trial, "z_complex": z.copy(), "residual": residual,
                "realness": norm(z.imag)/max(norm(z), 1e-300), **conditioning}
        roots.append(root); rec.update({"status": "new-root", "root_id": root["id"]}); trials.append(rec)
    encoded = []
    for root in roots:
        r = dict(root); r["z"] = [[float(x.real), float(x.imag)] for x in r.pop("z_complex")]; encoded.append(r)
    return {"script": Path(__file__).name, "version": 321, "standalone": True, "case": f"{n},{d}", "n": n, "degree": d,
            "oracle": {"kind": oracle.kind, "seed": oracle.seed, "samples": oracle.eval_count, "seconds": oracle.seconds,
                       "residual_mode": "raw" if isinstance(oracle, ExpressionOracle) else "root-equivalent-degree-normalized"},
            "parameters": {"count": args.count, "pool": args.pool, "epochs": args.epochs, "accept": args.accept,
                "broyden": args.broyden, "deflation_alpha": args.deflation_alpha, "irp": args.irp, "swarm": args.swarm},
            "swarm": swarm_meta, "roots": encoded, "trials": trials if args.verbose_trials else trials[:args.keep_trials],
            "summary": {"requested": args.count, "unique": len(roots), "success": len(roots) >= args.count,
                "trials": len(trials), "duplicates": duplicates, "failures": failures,
                "oracle_samples": oracle.eval_count, "seconds": clock()-started}}


# ---------- CLI and deterministic regressions -------------------------------

def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Pandrosion 321 compact standalone NumPy engine")
    p.add_argument("--cases", default="2,4"); p.add_argument("--seed-index", type=int, default=0)
    p.add_argument("--system-source", choices=["ks", "kostlan", "poly", "polynomial"], default="ks")
    p.add_argument("--polys"); p.add_argument("--variables"); p.add_argument("--starts")
    p.add_argument("--dense-max-terms", type=int, default=250000); p.add_argument("--features", type=int, default=2048)
    p.add_argument("--equation-normalize", action="store_true", default=False)
    p.add_argument("--count", type=int, default=8); p.add_argument("--pool", type=int, default=1024)
    p.add_argument("--epochs", type=int, default=18); p.add_argument("--accept", type=float, default=1e-8)
    p.add_argument("--tol", type=float, default=1e-12); p.add_argument("--cluster-sep", type=float, default=1e-7)
    p.add_argument("--early-dup-sep", type=float, default=1e-4); p.add_argument("--trial-timeout", type=float, default=0)
    p.add_argument("--jet-radius", type=float, default=1e-5); p.add_argument("--jet-refresh", type=int, default=4)
    p.add_argument("--broyden", action="store_true", default=True); p.add_argument("--no-broyden", dest="broyden", action="store_false")
    p.add_argument("--adaptive-lm", action="store_true", default=True); p.add_argument("--no-adaptive-lm", dest="adaptive_lm", action="store_false")
    p.add_argument("--lm-damping", type=float, default=0); p.add_argument("--trust-radius", type=float, default=10)
    p.add_argument("--deflation-alpha", type=float, default=.15); p.add_argument("--line-grid", default="1,.75,.5,.25,.125,.0625")
    p.add_argument("--parabolic", action="store_true", default=True); p.add_argument("--no-parabolic", dest="parabolic", action="store_false")
    p.add_argument("--irp", action="store_true", default=True); p.add_argument("--no-irp", dest="irp", action="store_false")
    p.add_argument("--irp-epochs", type=int, default=4); p.add_argument("--irp-top", type=int, default=4)
    p.add_argument("--irp-scales", default="1,1.2599210499,.793700526,2,.5")
    p.add_argument("--swarm", action="store_true", default=True); p.add_argument("--no-swarm", dest="swarm", action="store_false")
    p.add_argument("--swarm-size", type=int, default=0); p.add_argument("--swarm-keep", type=int, default=0)
    p.add_argument("--swarm-iters", type=int, default=2); p.add_argument("--swarm-sep", type=float, default=0)
    p.add_argument("--polish-steps", type=int, default=4); p.add_argument("--polish-gate", type=float, default=1e-2)
    p.add_argument("--keep-trials", type=int, default=100)
    p.add_argument("--verbose-trials", action="store_true"); p.add_argument("--self-test", action="store_true")
    p.add_argument("--out"); p.add_argument("--outdir", default="/tmp/321_pandrosion")
    return p


def validate(args: argparse.Namespace) -> None:
    for raw in str(args.cases).replace("|", ";").split(";"):
        if raw.strip(): parse_case(raw)
    if args.count < 0 or args.pool < 0 or args.epochs <= 0: raise ValueError("invalid count/pool/epochs")
    if args.accept <= 0 or args.tol < 0 or args.jet_radius <= 0 or args.polish_gate < 0: raise ValueError("invalid numerical tolerance")
    if args.cluster_sep <= 0 or args.early_dup_sep < 0 or args.deflation_alpha < 0: raise ValueError("invalid separation/deflation")
    if args.dense_max_terms <= 0 or args.features <= 0 or args.jet_refresh <= 0 or args.swarm_iters < 0: raise ValueError("invalid backend/iteration budget")


def self_test(args: argparse.Namespace) -> dict[str, Any]:
    specs = [("two-roots", "1,2", "x^2-3*x-10", "-8,4", 2),
             ("four-roots", "2,2", "x1^2-1;x2^2-1", "-2,-2;-2,2;2,-2;2,2", 4),
             ("multiple-root", "1,2", "x^2", "0", 1)]
    checks = []
    for name, case, poly, starts, count in specs:
        a = argparse.Namespace(**vars(args)); a.self_test = False; a.system_source = "poly"; a.cases = case
        a.polys = poly; a.variables = None; a.starts = starts; a.count = count; a.pool = max(8, count)
        a.swarm = False; a.epochs = 20; result = run_case(a, case)
        passed = bool(result["summary"]["success"])
        if name == "multiple-root": passed = bool(passed and result["roots"][0]["singular"] and result["roots"][0]["cond"] is None)
        checks.append({"name": name, "passed": passed, "result": result})
    idx, _ = line_choice(np.asarray([11., 9.]), np.asarray([[1+0j], [.5+0j]]), np.asarray([0j]),
                         10., 10., [np.asarray([.5+0j])], 1.)
    checks.append({"name": "acceptance-safe-deflation", "passed": idx == 1, "selected": idx})
    a = argparse.Namespace(**vars(args)); a.system_source = "poly"; a.polys = "x^2-1"; a.variables = None
    oracle = make_oracle(a, 1, 2); target = Target(oracle)
    loc = correct(target, np.asarray([2+0j]), a, [], 2, None)
    checks.append({"name": "immediate-broyden", "passed": loc.rebuilds == 1 and loc.broyden >= 1,
                   "rebuilds": loc.rebuilds, "updates": loc.broyden})
    return {"script": Path(__file__).name, "self_test": True, "passed": all(c["passed"] for c in checks), "checks": checks}


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parser().parse_args(argv); validate(args)
    if args.self_test:
        final = self_test(args); cases = ["self-test"]
    else:
        cases = [x.strip() for x in str(args.cases).replace("|", ";").split(";") if x.strip()]
        results = [run_case(args, case) for case in cases]
        final = results[0] if len(results) == 1 else {"script": Path(__file__).name, "cases": results}
    out = Path(args.out) if args.out else Path(args.outdir) / ("self_test.json" if args.self_test else f"321_{cases[0].replace(',', 'x')}.json")
    out.parent.mkdir(parents=True, exist_ok=True); out.write_text(json.dumps(strict_json(final), indent=2, allow_nan=False), encoding="utf-8")
    if args.self_test:
        print(f"321 self-test: {'PASS' if final['passed'] else 'FAIL'} ({sum(int(c['passed']) for c in final['checks'])}/{len(final['checks'])})")
        print(f"out={out}"); return 0 if final["passed"] else 1
    for result in (results if len(cases) > 1 else [final]):
        s = result["summary"]; print(f"321 case={result['case']} oracle={result['oracle']['kind']} roots={s['unique']}/{s['requested']} trials={s['trials']} samples={s['oracle_samples']} seconds={s['seconds']:.3f}")
    print(f"out={out}"); return 0


if __name__ == "__main__":
    raise SystemExit(main())
