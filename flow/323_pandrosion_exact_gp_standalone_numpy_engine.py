"""Pandrosion 323: compact standalone NumPy root-harvesting engine.

Core geometry: logarithmic atlas starts, vectorised swarm, paired local jets,
immediate Broyden updates, overflow-safe root deflation, augmented-system LM,
adaptive second-order inverse jets, and direct/reciprocal IRP rescue.  Exact
expression polynomials and dense Kostlan systems are supported.  Large KS cases
default to an exact-in-law Gaussian-process oracle sampled conditionally through
an incremental rank-adaptive covariance factor, so no monomial or feature
truncation is needed.
The validated stratified surrogate remains available as an explicit fallback.
No local project imports are used.
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


def stable_seed(n: int, d: int, index: int, salt: int = 0) -> int:
    return int(mix64(0x50414E44524F5349 + 1000003 * n + 9176 * d + 97 * index + salt) & 0x7FFFFFFF)


def direction(n: int, index: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(mix64(seed + 104729*index))
    v = rng.standard_normal(n) + 1j*rng.standard_normal(n)
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

    def backward_error(self, z: Any) -> float:
        return norm(self.eval(z))


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
    def __init__(self, n: int, d: int, seed: int, dense_cap: int, features: int, normalize: bool, eval_block: int) -> None:
        terms = math.comb(n + d, d)
        rng = np.random.default_rng(seed)
        if terms <= dense_cap:
            exps = compositions(d, n)
            totals = exps.sum(1)
            logs = np.asarray([math.lgamma(d + 1) - math.lgamma(d - int(t) + 1)
                               - sum(math.lgamma(int(a) + 1) for a in row)
                               for row, t in zip(exps, totals)])
            log_scales = 0.5*logs; scales = np.exp(log_scales - float(np.max(log_scales)))
            kind, feature_mode, feature_logs = "kostlan-dense", False, None
        else:
            m = max(n + 8, features)
            exps = np.zeros((m, n), np.int16 if d < 32767 else np.int32); degrees = np.zeros(m, int)
            idx = 1
            for j in range(n):
                if idx >= m: break
                exps[idx, j] = 1; degrees[idx] = 1; idx += 1
            while idx < m:
                q = (idx-(n+1)+.5)/max(1, m-(n+1)); k = min(d, max(0, int(q*(d+1))))
                if idx % 7 == 0: k = int(rng.integers(0, d+1))
                if k: exps[idx] = rng.multinomial(k, np.full(n, 1/n))
                degrees[idx] = k; idx += 1
            feature_logs = np.asarray([.5*(math.lgamma(d+1)-math.lgamma(k+1)-math.lgamma(d-k+1)
                                    + math.log(d+1)+k*math.log(max(1, n))-math.log(m)) for k in degrees])
            scales = np.ones(m); kind, feature_mode = "kostlan-stratified-feature-surrogate", True
        coeff = (rng.standard_normal((n, len(exps))) + 1j * rng.standard_normal((n, len(exps)))) / math.sqrt(2)
        coeff *= scales[None, :]
        if normalize:
            coeff /= np.maximum(np.linalg.norm(coeff, axis=1), 1e-300)[:, None]
        super().__init__(n, d, seed, kind)
        self.exps, self.degrees, self.coeff = exps, exps.sum(1), np.asarray(coeff, np.complex128)
        self.feature_mode, self.feature_logs, self.eval_block = feature_mode, feature_logs, max(1, eval_block)

    def _basis(self, Z: np.ndarray) -> np.ndarray:
        if self.feature_mode:
            az = np.abs(Z); logamp = np.log(np.maximum(az, 1e-300)) @ self.exps.T + self.feature_logs[None, :]
            logamp -= .5*self.d*np.log1p(np.sum(az*az, axis=1))[:, None]
            shift = np.max(np.where(np.isfinite(logamp), logamp, -np.inf), axis=1)
            logamp -= np.where(np.isfinite(shift), shift, 0)[:, None]
            return np.exp(np.clip(logamp, -745, 0))*np.exp(1j*(np.angle(Z) @ self.exps.T))
        scale = np.maximum(1, np.max(np.abs(Z), axis=1)); W = Z/scale[:, None]
        mon = np.ones((len(Z), len(self.exps)), np.complex128)
        for j in range(self.n):
            powers = np.ones((len(Z), self.d + 1), np.complex128)
            for k in range(1, self.d + 1): powers[:, k] = powers[:, k - 1] * W[:, j]
            mon *= powers[:, self.exps[:, j]]
        mon *= np.exp(np.clip(np.log(scale)[:, None]*(self.degrees[None, :]-self.d), -745, 0))
        return mon

    def _eval_batch(self, Z: np.ndarray) -> np.ndarray:
        return np.vstack([self._basis(Z[i:i+self.eval_block]) @ self.coeff.T for i in range(0, len(Z), self.eval_block)])

    def backward_error(self, z: Any) -> float:
        t0 = clock(); Z = np.asarray(z, np.complex128)[None, :]
        with np.errstate(all="ignore"):
            phi = self._basis(Z); f = phi @ self.coeff.T
            value = norm(f[0])/max(norm(phi[0])*norm(self.coeff), 1e-300)
        self.eval_count += 1; self.seconds += clock()-t0
        return value if math.isfinite(value) else float("inf")


class ExactKSGPOracle(Oracle):
    """Lazy exact-in-law normalized Kostlan Gaussian-process oracle.

    Every adaptively requested finite set has the exact joint KS distribution.
    The incremental covariance factor avoids monomial enumeration.  It stays
    block-triangular while full rank and drops only numerically null conditional
    modes; every such stabilization is reported.
    """
    def __init__(self, n: int, d: int, seed: int, jitter: float, max_points: int, min_radius: float) -> None:
        super().__init__(n, d, seed, "kostlan-exact-gp-lazy")
        self.rng = np.random.default_rng(seed); self.jitter = max(1e-16, jitter)
        self.max_points, self.recommended_radius = max_points, max(1e-6, min_radius)
        self.points = np.empty((0, n), np.complex128); self.L = np.empty((0, 0), np.complex128)
        self.xi = np.empty((0, n), np.complex128); self.cache: dict[bytes, np.ndarray] = {}
        self.stabilizations = 0; self.max_correction = 0.0
        self.min_conditional_eigenvalue = float("inf")

    def _key(self, z: np.ndarray) -> bytes:
        return np.ascontiguousarray(z, np.complex128).view(np.float64).tobytes()

    def kernel(self, A: Any, B: Any) -> np.ndarray:
        a, b = np.asarray(A, np.complex128), np.asarray(B, np.complex128)
        if a.ndim == 1: a = a[None, :]
        if b.ndim == 1: b = b[None, :]
        den = np.sqrt(1+np.sum(np.abs(a)**2, axis=1))[:, None]*np.sqrt(1+np.sum(np.abs(b)**2, axis=1))[None, :]
        base = (1+a@b.conj().T)/np.maximum(den, 1e-300)
        mag = np.abs(base); base = np.where(mag > 1, base/np.maximum(mag, 1e-300), base)
        return base**self.d

    def _solve_factor(self, B: np.ndarray) -> np.ndarray:
        X = np.empty_like(B)
        for i in range(len(self.L)):
            X[i] = (B[i]-self.L[i, :i]@X[:i])/self.L[i, i]
        return X

    def _append(self, Z: np.ndarray) -> None:
        for z in Z:
            key = self._key(z)
            if key in self.cache: continue
            if len(self.cache)+1 > self.max_points:
                raise RuntimeError(f"GP point cap exceeded: {len(self.cache)+1}>{self.max_points}")
            old = len(self.points)
            if old:
                k = self.kernel(self.points, z)[:, 0]; v = self._solve_factor(k[:, None])[:, 0]
                conditional = float((self.kernel(z, z)[0, 0]-np.vdot(v, v)).real)
                self.min_conditional_eigenvalue = min(self.min_conditional_eigenvalue, conditional)
            else:
                v = np.empty(0, np.complex128); conditional = float(self.kernel(z, z)[0, 0].real)
            floor = self.jitter
            if conditional > floor:
                sigma = math.sqrt(conditional)
                fresh = (self.rng.standard_normal(self.n)+1j*self.rng.standard_normal(self.n))/math.sqrt(2)
                fnew = v.conj()@self.xi + sigma*fresh
                enlarged = np.zeros((old+1, old+1), np.complex128)
                enlarged[:old, :old] = self.L; enlarged[old, :old] = v.conj(); enlarged[old, old] = sigma
                self.L = enlarged; self.xi = np.vstack((self.xi, fresh)); self.points = np.vstack((self.points, z))
            else:
                self.stabilizations += 1; self.max_correction = max(self.max_correction, abs(conditional))
                fnew = v.conj()@self.xi
            self.cache[key] = np.asarray(fnew, np.complex128)

    def _eval_batch(self, Z: np.ndarray) -> np.ndarray:
        self._append(Z)
        return np.asarray([self.cache[self._key(z)] for z in Z], np.complex128)

    def backward_error(self, z: Any) -> float:
        return norm(self.eval(z))/math.sqrt(max(1, self.n))


def make_oracle(args: argparse.Namespace, n: int, d: int) -> Oracle:
    if args.system_source in {"poly", "polynomial"}:
        seed = stable_seed(n, d, args.seed_index)
        if not args.polys:
            raise ValueError("--polys is required for polynomial systems")
        return ExpressionOracle(n, d, seed, args.polys, args.variables)
    terms = math.comb(n+d, d); backend = args.ks_backend
    if backend == "auto": backend = "dense" if terms <= args.dense_max_terms else "gp"
    if backend == "dense":
        if terms > args.dense_max_terms: raise ValueError(f"dense KS has {terms} terms/poly; raise --dense-max-terms explicitly")
        return KostlanOracle(n, d, stable_seed(n, d, args.seed_index), terms, args.features, args.equation_normalize, args.eval_block)
    if backend == "feature":
        return KostlanOracle(n, d, stable_seed(n, d, args.seed_index, 0x314314), 0, args.features, args.equation_normalize, args.eval_block)
    return ExactKSGPOracle(n, d, stable_seed(n, d, args.seed_index, 0x323600), args.gp_jitter, args.gp_max_points, args.gp_jet_radius)


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
    hrel = max(hrel, float(getattr(target.oracle, "recommended_radius", 0)))
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
    jscale = float(np.max(np.abs(J))) if J.size else 0
    if not math.isfinite(jscale) or jscale <= 1e-300:
        return np.zeros_like(y), "failed"
    with np.errstate(all="ignore"):
        jn, fn = J/jscale, f/jscale
    if not np.all(np.isfinite(jn)) or not np.all(np.isfinite(fn)):
        return np.zeros_like(y), "failed"
    try:
        if mu > 0:
            aug = np.vstack((jn, math.sqrt(mu)*np.eye(len(y), dtype=np.complex128)))
            rhs = np.concatenate((-fn, np.zeros(len(y), np.complex128)))
            delta, _, _, _ = np.linalg.lstsq(aug, rhs, rcond=1e-12); method = "lm-augmented"
        else:
            delta = np.linalg.solve(jn, -fn); method = "solve"
    except Exception:
        try:
            delta, _, _, _ = np.linalg.lstsq(jn, -fn, rcond=1e-12); method = "svd-lstsq"
        except Exception:
            delta = np.zeros_like(y); method = "failed"
    limit = trust * max(1, norm(y)) if trust > 0 else 10 * max(1, norm(y))
    nd = norm(delta)
    if math.isfinite(nd) and nd > limit: delta *= limit / nd
    return np.asarray(delta, np.complex128), method


def higher_order_step(target: Target, y: np.ndarray, f: np.ndarray, J: np.ndarray,
                      d1: np.ndarray, r: float, args: argparse.Namespace) -> tuple[np.ndarray, int, bool]:
    if not args.higher_order or r > args.higher_order_gate or norm(f+J@d1) > args.higher_order_trigger*r:
        return d1, 0, False
    try:
        cond = float(np.linalg.cond(J/np.max(np.abs(J))))
    except Exception:
        return d1, 0, False
    nd = norm(d1)
    if not math.isfinite(cond) or cond > args.higher_order_cond or nd <= 1e-14:
        return d1, 0, False
    h = max(1e-4, min(1e-2, .05*max(1, norm(y))/nd))
    fp, fm = target.eval(y+h*d1), target.eval(y-h*d1)
    curvature = (fp-2*f+fm)/(h*h)
    d2, _ = linear_step(J, .5*curvature, 0, args.trust_radius, y)
    if np.all(np.isfinite(d2)) and norm(d2) <= args.higher_order_ratio*nd:
        return d1+d2, 2, True
    return d1, 2, False


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
    higher_order_evals: int
    higher_order_used: int
    mu: float


def correct(target: Target, y0: Any, args: argparse.Namespace, known: Sequence[np.ndarray],
            epochs: int, deadline: Optional[float]) -> CorrectResult:
    y = np.asarray(y0, np.complex128).copy(); f = target.eval(y); r = norm(f)
    best_y, best_f, best_r = y.copy(), f.copy(), r
    A: Optional[np.ndarray] = None; age = rebuilds = updates = line_evals = jet_evals = para = high_evals = high_used = done = 0
    mu = max(0, args.lm_damping); status = "max-epochs"
    full_line = np.asarray(floats(args.line_grid, [1, .75, .5, .25, .125, .0625]))
    short_line = full_line[:min(4, len(full_line))]
    for ep in range(max(1, epochs)):
        if deadline is not None and clock() >= deadline: status = "timeout"; break
        if best_r < args.accept: status = "converged"; break
        force_rebuild = A is None or not args.broyden or age >= args.jet_refresh
        moved = False
        for attempt in range(6):
            rebuilt = force_rebuild or A is None
            if rebuilt:
                _, batch_J, used = paired_jets(target, y[None, :], f[None, :], args.jet_radius)
                A = batch_J[0]; age = 0; rebuilds += 1; jet_evals += used; force_rebuild = False
            delta, _ = linear_step(A, f, mu if args.adaptive_lm else args.lm_damping, args.trust_radius, y)
            if rebuilt:
                delta, extra, used_high = higher_order_step(target, y, f, A, delta, r, args)
                jet_evals += extra; high_evals += extra; high_used += int(used_high)
            if not np.all(np.isfinite(delta)) or norm(delta) <= 1e-15:
                if args.adaptive_lm and attempt < 5:
                    mu = min(1, max(1e-3, mu*10)); continue
                status = "invalid-step"; break
            L = full_line if rebuilt else short_line
            Y = y[None, :] + L[:, None] * delta[None, :]
            F = target.eval_batch(Y); R = norms(F); line_evals += len(Y)
            i, merit = line_choice(R, Y, y, r, best_r, known, args.deflation_alpha)
            if i is None:
                if not rebuilt and age > 0:
                    force_rebuild = True; continue
                if args.adaptive_lm and attempt < 5:
                    mu = min(1, max(1e-3, mu*10)); continue
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
            with np.errstate(all="ignore"):
                predicted = norm(f + A @ (lam*delta)) if A is not None else r
            dy, df = yn-y, fn-f; denom = complex(np.vdot(dy, dy))
            if args.broyden and abs(denom) > 1e-300:
                with np.errstate(all="ignore"):
                    candidate_A = A + np.outer(df-A@dy, dy.conj())/denom
                if np.all(np.isfinite(candidate_A)):
                    A = candidate_A; age += 1; updates += 1
                else:
                    A = None; age = args.jet_refresh
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
                         line_evals, jet_evals, para, high_evals, high_used, mu)


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
                                 result.jet_evals + loc.jet_evals, result.parabolic_evals + loc.parabolic_evals,
                                 result.higher_order_evals + loc.higher_order_evals,
                                 result.higher_order_used + loc.higher_order_used, loc.mu)
        if best.ok: break
    return best


# ---------- starts, swarm and root finishing --------------------------------

def atlas_start(n: int, trial: int, seed: int) -> np.ndarray:
    radii = [2**(k/3) for k in range(-15, 16)]
    return radii[(trial*17) % len(radii)] * direction(n, trial, seed)


def swarm(base: Target, args: argparse.Namespace, n: int, seed: int) -> tuple[list[np.ndarray], dict[str, Any]]:
    size = args.swarm_size or min(args.pool, max(32, 8*args.count))
    if isinstance(base.oracle, ExactKSGPOracle): size = min(size, args.gp_swarm_cap)
    keep = args.swarm_keep or max(8, 3*args.count)
    Y = np.stack([atlas_start(n, i, seed+31337) for i in range(max(1, size))])
    F = base.eval_batch(Y); R = norms(F); alive = np.isfinite(R); jet_samples = line_samples = 0
    L = np.asarray([1, .5, .25, .1])
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
    origin_error = oracle.backward_error(np.zeros(n, np.complex128))
    explicit = starts_from_text(args.starts, n)
    swarm_points, swarm_meta, cap_error = ([], {"enabled": False}, None)
    if args.swarm and not explicit and args.count:
        try:
            swarm_points, swarm_meta = swarm(base, args, n, oracle.seed); swarm_meta["enabled"] = True
        except RuntimeError as exc:
            cap_error = str(exc); swarm_meta = {"enabled": True, "status": "oracle-cap", "error": cap_error}
    priority = explicit + swarm_points; roots: list[dict[str, Any]] = []; trials = []; duplicates = failures = 0
    for trial in range(args.pool):
        if len(roots) >= args.count: break
        y0 = priority[trial].copy() if trial < len(priority) else atlas_start(n, trial-len(priority), oracle.seed)
        known = [np.asarray(r["z_complex"], np.complex128) for r in roots]
        if known and min(norm(y0-r) for r in known) <= args.early_dup_sep:
            duplicates += 1; trials.append({"trial": trial, "status": "start-near-root", "accepted": False}); continue
        deadline = clock()+args.trial_timeout if args.trial_timeout > 0 else None
        before = oracle.eval_count
        try:
            loc = correct(base, y0, args, known, args.epochs, deadline)
            loc = irp_rescue(base, loc, args, known, deadline)
            if loc.residual <= max(args.polish_gate, 100*args.accept):
                z, f, residual, conditioning = polish(base, loc.y, args, deadline)
            else:
                z, f, residual = loc.y.copy(), loc.f.copy(), loc.residual
                conditioning = {"smin": None, "smax": None, "cond": None, "near_multiple": None, "singular": None, "status": "polish-gate"}
            validation_error = oracle.backward_error(z)
        except RuntimeError as exc:
            cap_error = str(exc); failures += 1
            trials.append({"trial": trial, "status": "oracle-cap", "accepted": False,
                           "error": cap_error, "oracle_samples": oracle.eval_count-before})
            break
        accepted = bool(math.isfinite(residual) and residual < args.accept and validation_error < args.validation_accept)
        rec = {"trial": trial, "status": loc.status, "accepted": accepted, "residual": residual,
               "validation_error": validation_error,
               "epochs": loc.epochs, "jet_rebuilds": loc.rebuilds, "broyden_updates": loc.broyden,
               "jet_evals": loc.jet_evals, "line_evals": loc.line_evals, "parabolic_evals": loc.parabolic_evals,
               "higher_order_evals": loc.higher_order_evals, "higher_order_used": loc.higher_order_used,
               "oracle_samples": oracle.eval_count-before, "conditioning": conditioning}
        if args.verbose_trials: rec.update({"start": y0, "z": z})
        if not accepted:
            if residual < args.accept: rec["status"] = "validation-failed"
            failures += 1; trials.append(rec); continue
        dup = next((i for i, root in enumerate(roots) if norm(z-root["z_complex"]) <= args.cluster_sep), None)
        if dup is not None:
            duplicates += 1; rec.update({"status": "duplicate", "cluster": dup}); trials.append(rec); continue
        root = {"id": len(roots), "trial": trial, "z_complex": z.copy(), "residual": residual, "validation_error": validation_error,
                "realness": norm(z.imag)/max(norm(z), 1e-300), **conditioning}
        roots.append(root); rec.update({"status": "new-root", "root_id": root["id"]}); trials.append(rec)
    encoded = []
    for root in roots:
        r = dict(root); r["z"] = [[float(x.real), float(x.imag)] for x in r.pop("z_complex")]; encoded.append(r)
    return {"script": Path(__file__).name, "version": 323, "standalone": True, "case": f"{n},{d}", "n": n, "degree": d,
            "oracle": {"kind": oracle.kind, "seed": oracle.seed, "samples": oracle.eval_count, "seconds": oracle.seconds,
                       "residual_mode": "raw" if isinstance(oracle, ExpressionOracle) else "root-equivalent-degree-normalized",
                       "validation_mode": ("raw-residual" if isinstance(oracle, ExpressionOracle) else
                                           "normalized-gp-residual" if isinstance(oracle, ExactKSGPOracle) else "scale-invariant-backward-error"),
                       "origin_backward_error": origin_error,
                       "exact_terms_per_polynomial": math.comb(n+d, d),
                       "feature_count": (len(oracle.exps) if isinstance(oracle, KostlanOracle) else None),
                       "multiresolution_stability_established": (None if isinstance(oracle, ExactKSGPOracle) else
                                                                  not isinstance(oracle, KostlanOracle) or not oracle.feature_mode),
                       "surrogate_warning": ("root is validated only for this fixed feature bank; nested-bank convergence is not established"
                                             if isinstance(oracle, KostlanOracle) and oracle.feature_mode else None),
                       "feature_degree_range": ([int(np.min(oracle.degrees)), int(np.max(oracle.degrees))] if isinstance(oracle, KostlanOracle) and oracle.feature_mode else None),
                       "gp_unique_points": (len(oracle.cache) if isinstance(oracle, ExactKSGPOracle) else None),
                       "gp_active_covariance_rank": (len(oracle.points) if isinstance(oracle, ExactKSGPOracle) else None),
                       "gp_stabilizations": (oracle.stabilizations if isinstance(oracle, ExactKSGPOracle) else None),
                       "gp_max_covariance_correction": (oracle.max_correction if isinstance(oracle, ExactKSGPOracle) else None),
                       "gp_min_conditional_eigenvalue": (oracle.min_conditional_eigenvalue if isinstance(oracle, ExactKSGPOracle) else None),
                       "gp_rank_tolerance": (oracle.jitter if isinstance(oracle, ExactKSGPOracle) else None),
                       "gp_exact_finite_dimensional_law": ("up to floating-point arithmetic and the reported rank tolerance"
                                                           if isinstance(oracle, ExactKSGPOracle) else None),
                       "gp_reproducibility_scope": ("exact in law for this adaptive query sequence; another query order samples another realization"
                                                    if isinstance(oracle, ExactKSGPOracle) else None)},
            "parameters": {"count": args.count, "pool": args.pool, "epochs": args.epochs, "accept": args.accept,
                "ks_backend": args.ks_backend, "broyden": args.broyden, "higher_order": args.higher_order, "validation_accept": args.validation_accept,
                "deflation_alpha": args.deflation_alpha, "irp": args.irp, "swarm": args.swarm},
            "swarm": swarm_meta, "roots": encoded, "trials": trials if args.verbose_trials else trials[:args.keep_trials],
            "summary": {"requested": args.count, "unique": len(roots), "success": len(roots) >= args.count,
                "trials": len(trials), "duplicates": duplicates, "failures": failures,
                "oracle_samples": oracle.eval_count, "seconds": clock()-started, "oracle_cap_error": cap_error}}


# ---------- CLI and deterministic regressions -------------------------------

def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Pandrosion 323 compact standalone NumPy engine")
    p.add_argument("--cases", default="2,4"); p.add_argument("--seed-index", type=int, default=0)
    p.add_argument("--system-source", choices=["ks", "kostlan", "poly", "polynomial"], default="ks")
    p.add_argument("--ks-backend", choices=["auto", "dense", "gp", "feature"], default="auto")
    p.add_argument("--polys"); p.add_argument("--variables"); p.add_argument("--starts")
    p.add_argument("--dense-max-terms", type=int, default=250000); p.add_argument("--features", type=int, default=2048)
    p.add_argument("--eval-block", type=int, default=128)
    p.add_argument("--gp-jitter", type=float, default=1e-9); p.add_argument("--gp-max-points", type=int, default=10000)
    p.add_argument("--gp-jet-radius", type=float, default=1e-3); p.add_argument("--gp-swarm-cap", type=int, default=8)
    p.add_argument("--equation-normalize", action="store_true", default=False)
    p.add_argument("--count", type=int, default=8); p.add_argument("--pool", type=int, default=1024)
    p.add_argument("--epochs", type=int, default=18); p.add_argument("--accept", type=float, default=1e-8)
    p.add_argument("--tol", type=float, default=1e-12); p.add_argument("--cluster-sep", type=float, default=1e-7)
    p.add_argument("--early-dup-sep", type=float, default=1e-4); p.add_argument("--trial-timeout", type=float, default=0)
    p.add_argument("--jet-radius", type=float, default=1e-5); p.add_argument("--jet-refresh", type=int, default=4)
    p.add_argument("--broyden", action="store_true", default=True); p.add_argument("--no-broyden", dest="broyden", action="store_false")
    p.add_argument("--adaptive-lm", action="store_true", default=True); p.add_argument("--no-adaptive-lm", dest="adaptive_lm", action="store_false")
    p.add_argument("--lm-damping", type=float, default=0); p.add_argument("--trust-radius", type=float, default=10)
    p.add_argument("--higher-order", action="store_true", default=True); p.add_argument("--no-higher-order", dest="higher_order", action="store_false")
    p.add_argument("--higher-order-gate", type=float, default=.1); p.add_argument("--higher-order-trigger", type=float, default=.25)
    p.add_argument("--higher-order-cond", type=float, default=1e8); p.add_argument("--higher-order-ratio", type=float, default=.75)
    p.add_argument("--deflation-alpha", type=float, default=.15); p.add_argument("--line-grid", default="1,.75,.5,.25,.125,.0625")
    p.add_argument("--parabolic", action="store_true", default=True); p.add_argument("--no-parabolic", dest="parabolic", action="store_false")
    p.add_argument("--irp", action="store_true", default=True); p.add_argument("--no-irp", dest="irp", action="store_false")
    p.add_argument("--irp-epochs", type=int, default=4); p.add_argument("--irp-top", type=int, default=4)
    p.add_argument("--irp-scales", default="1,1.2599210499,.793700526,2,.5")
    p.add_argument("--swarm", action="store_true", default=True); p.add_argument("--no-swarm", dest="swarm", action="store_false")
    p.add_argument("--swarm-size", type=int, default=0); p.add_argument("--swarm-keep", type=int, default=0)
    p.add_argument("--swarm-iters", type=int, default=2); p.add_argument("--swarm-sep", type=float, default=0)
    p.add_argument("--polish-steps", type=int, default=4); p.add_argument("--polish-gate", type=float, default=1e-2)
    p.add_argument("--validation-accept", type=float, default=1e-8)
    p.add_argument("--keep-trials", type=int, default=100)
    p.add_argument("--verbose-trials", action="store_true"); p.add_argument("--self-test", action="store_true")
    p.add_argument("--out"); p.add_argument("--outdir", default="/tmp/323_pandrosion")
    return p


def validate(args: argparse.Namespace) -> None:
    for raw in str(args.cases).replace("|", ";").split(";"):
        if raw.strip(): parse_case(raw)
    if args.count < 0 or args.pool < 0 or args.epochs <= 0: raise ValueError("invalid count/pool/epochs")
    if args.accept <= 0 or args.tol < 0 or args.jet_radius <= 0 or args.polish_gate < 0: raise ValueError("invalid numerical tolerance")
    if args.cluster_sep <= 0 or args.early_dup_sep < 0 or args.deflation_alpha < 0: raise ValueError("invalid separation/deflation")
    if args.dense_max_terms <= 0 or args.features <= 0 or args.eval_block <= 0 or args.jet_refresh <= 0 or args.swarm_iters < 0: raise ValueError("invalid backend/iteration budget")
    if args.gp_jitter <= 0 or args.gp_max_points <= 0 or args.gp_jet_radius <= 0 or args.gp_swarm_cap <= 0: raise ValueError("invalid GP controls")
    if args.validation_accept <= 0 or args.higher_order_gate < 0 or not 0 < args.higher_order_trigger <= 1 or args.higher_order_cond <= 0 or args.higher_order_ratio <= 0: raise ValueError("invalid validation/higher-order controls")


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
    a.higher_order_gate = 1.; loc = correct(target, np.asarray([1.2+0j]), a, [], 1, None)
    checks.append({"name": "adaptive-higher-order", "passed": loc.higher_order_used >= 1, "used": loc.higher_order_used})
    huge, method = linear_step(np.eye(2, dtype=np.complex128)*1e200, np.ones(2, np.complex128)*1e200, .1, 10, np.zeros(2))
    checks.append({"name": "scaled-lm", "passed": bool(np.all(np.isfinite(huge)) and norm(huge) > 0), "method": method})
    a = argparse.Namespace(**vars(args)); a.self_test = False; a.system_source = "ks"; a.cases = "2,3"
    a.starts = None; a.count = 2; a.pool = 32; a.epochs = 14; a.swarm = True; a.swarm_size = 32
    result = run_case(a, "2,3")
    checks.append({"name": "kostlan-smoke", "passed": bool(result["summary"]["success"]), "result": result})
    a = argparse.Namespace(**vars(args)); a.system_source = "ks"; a.ks_backend = "feature"; a.dense_max_terms = 1; a.features = 512
    feature = make_oracle(a, 20, 20); origin = feature.backward_error(np.zeros(20, np.complex128))
    axes = sum(bool(np.sum(row) == 1) for row in feature.exps)
    checks.append({"name": "feature-origin-guard", "passed": bool(feature.feature_mode and np.min(feature.degrees) == 0 and axes >= 20 and origin > a.validation_accept),
                   "origin_backward_error": origin, "degree_range": [int(np.min(feature.degrees)), int(np.max(feature.degrees))], "linear_features": axes})
    a.features = 2048; a.count = 1; a.pool = 8; a.epochs = 40; a.swarm = True; a.swarm_size = 16
    a.swarm_iters = 3; a.swarm_keep = 16; a.irp_epochs = 8; result = run_case(a, "20,20")
    root = result["roots"][0] if result["roots"] else {}
    root_norm = norm(np.asarray([complex(*q) for q in root.get("z", [])]))
    checks.append({"name": "feature-20x20-nondegenerate", "passed": bool(result["summary"]["success"] and root_norm > 1e-3 and root.get("validation_error", 1) < a.validation_accept),
                   "root_norm": root_norm, "validation_error": root.get("validation_error"), "result": result})
    a = argparse.Namespace(**vars(args)); a.system_source = "ks"; a.ks_backend = "gp"; a.gp_max_points = 256
    gp = make_oracle(a, 2, 3); Z = np.asarray([[0, 0], [.2+.1j, -.3j], [1, -.5j]], np.complex128)
    first = gp.eval_batch(Z); cached = gp.eval_batch(Z); covariance_error = norm(gp.L@gp.L.conj().T-gp.kernel(Z, Z))
    checks.append({"name": "exact-gp-covariance-cache", "passed": bool(np.array_equal(first, cached) and len(gp.cache) == 3 and covariance_error < 1e-9),
                   "covariance_error": covariance_error, "unique_points": len(gp.cache), "stabilizations": gp.stabilizations})
    dense = KostlanOracle(2, 3, 7, 10, 32, False, 32); phi = dense._basis(Z)
    dense_cov = phi@phi.conj().T; dense_cov /= np.sqrt(np.diag(dense_cov))[:, None]*np.sqrt(np.diag(dense_cov))[None, :]
    kernel_error = norm(dense_cov-gp.kernel(Z, Z))
    checks.append({"name": "gp-kernel-dense-equivalence", "passed": kernel_error < 1e-12, "kernel_error": kernel_error})
    finite = ExactKSGPOracle(1, 2, 9, 1e-9, 64, 1e-3)
    finite.eval_batch(np.linspace(-2, 2, 17, dtype=np.complex128)[:, None])
    checks.append({"name": "gp-finite-rank-stability", "passed": len(finite.cache) == 17 and len(finite.points) <= 3,
                   "queries": len(finite.cache), "active_rank": len(finite.points), "stabilizations": finite.stabilizations})
    return {"script": Path(__file__).name, "self_test": True, "passed": all(c["passed"] for c in checks), "checks": checks}


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parser().parse_args(argv); validate(args)
    if args.self_test:
        final = self_test(args); cases = ["self-test"]
    else:
        cases = [x.strip() for x in str(args.cases).replace("|", ";").split(";") if x.strip()]
        results = [run_case(args, case) for case in cases]
        final = results[0] if len(results) == 1 else {"script": Path(__file__).name, "cases": results}
    out = Path(args.out) if args.out else Path(args.outdir) / ("self_test.json" if args.self_test else f"323_{cases[0].replace(',', 'x')}.json")
    out.parent.mkdir(parents=True, exist_ok=True); out.write_text(json.dumps(strict_json(final), indent=2, allow_nan=False), encoding="utf-8")
    if args.self_test:
        print(f"323 self-test: {'PASS' if final['passed'] else 'FAIL'} ({sum(int(c['passed']) for c in final['checks'])}/{len(final['checks'])})")
        print(f"out={out}"); return 0 if final["passed"] else 1
    for result in (results if len(cases) > 1 else [final]):
        s = result["summary"]; print(f"323 case={result['case']} oracle={result['oracle']['kind']} roots={s['unique']}/{s['requested']} trials={s['trials']} samples={s['oracle_samples']} seconds={s['seconds']:.3f}")
    print(f"out={out}"); return 0


if __name__ == "__main__":
    raise SystemExit(main())
