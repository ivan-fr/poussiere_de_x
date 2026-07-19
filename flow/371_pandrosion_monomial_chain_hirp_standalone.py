"""Pandrosion 371: monomial-chain chart plus residual HIRP finite-slope solver.

Version 371 extends 370 and adapts the residual-homothety results of the accompanying
Pandrosion/Halley/IRP paper to general square black-box systems.  The raw
Pandrosion layer is still the exact coordinate-telescopic finite-slope solve
``Q p = -F(x)``.  The five oracle samples available after each telescope are
then used in one of three intrinsic residual programs:

* safeguarded raw dyadic globalization;
* true two-layer IRP, with the residual reevaluated and the finite-slope chart
  reclosed between layers; or
* HIRP-1: one directional Halley homothety followed by one reclosed residual
  Pandrosion layer.

Unlike 369's frozen one-shot analogue, the sequential IRP/HIRP branches use
actual intermediate residuals, matching the paper's defining requirement that
the residual and chart be recomputed after every layer.  A frozen diagonal
equation chart is also selected when equation scales span many decades; this
is a root-equivalent homothety and makes stopping/selection invariant to
arbitrary equation units.  Residual corrections are optionally snapped to a
dyadic logarithmic palette before return, mirroring the paper's chart-scaling
principle.

Every epoch still uses exactly ``n+5`` oracle evaluations.  There are no
analytic derivatives, Jacobians, quasi-Newton histories, homotopies, external
rescue solvers, restarts, or SciPy calls.  The implementation is standalone.
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

VERSION = "371"
BASE_TRAJECTORY = "370"


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
    if backend == "auto": backend = "dense" if terms <= args.dense_max_terms else ("gaussian" if n >= 12 and d >= 12 else "gp")
    if backend == "dense":
        if terms > args.dense_max_terms: raise ValueError(f"dense KS has {terms} terms/poly; raise --dense-max-terms explicitly")
        return KostlanOracle(n, d, stable_seed(n, d, args.seed_index), terms, args.features, args.equation_normalize, args.eval_block)
    if backend in {"feature", "gaussian"}:
        o = KostlanOracle(n, d, stable_seed(n, d, args.seed_index, 0x314314), 0, args.features, args.equation_normalize, args.eval_block)
        if backend == "gaussian":
            o.kind = "kostlan-gaussian-stratified-fixed-bank"
        return o
    return ExactKSGPOracle(n, d, stable_seed(n, d, args.seed_index, 0x325600), args.gp_jitter, args.gp_max_points, args.gp_jet_radius)


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


# ---------- exact finite Pandrosion slopes ---------------------------------

def finite_slope(target: Target, a: Any, b: Any) -> tuple[np.ndarray, float, int]:
    """Coordinate-telescopic Q with F(b)-F(a)=Q(a,b)(b-a), no derivatives."""
    aa, bb = np.asarray(a, np.complex128), np.asarray(b, np.complex128)
    n = len(aa); path = np.repeat(aa[None, :], n+1, axis=0)
    for j in range(n): path[j+1:, j] = bb[j]
    F = target.eval_batch(path); dz = bb-aa
    if np.any(np.abs(dz) <= 1e-14*np.maximum(1, np.maximum(np.abs(aa), np.abs(bb)))):
        raise ValueError("degenerate finite-slope endpoint")
    Q = (F[1:]-F[:-1]).T/dz[None, :]
    telescoped = Q@dz
    defect = norm((F[-1]-F[0])-telescoped)/max(norm(F[-1]-F[0]), norm(telescoped), 1e-300)
    return Q, defect, len(path)


def close_finite_slope(Q: np.ndarray, b: np.ndarray, fb: np.ndarray, y: np.ndarray,
                       f: np.ndarray, limit: float) -> tuple[np.ndarray,float,float,bool]:
    """Reclose one column so the transported matrix satisfies the finite identity exactly."""
    d=b-y; j=int(np.argmax(np.abs(d)))
    if abs(d[j])<=1e-14*max(1,norm(b),norm(y)): return Q,float("inf"),float("inf"),False
    mismatch=(fb-f)-Q@d; correction=mismatch/d[j]
    ratio=norm(correction)/max(norm(Q[:,j]),1e-300)
    if not math.isfinite(ratio) or ratio>limit: return Q,float("inf"),ratio,False
    out=Q.copy(); out[:,j]+=correction
    defect=norm((fb-f)-out@d)/max(norm(fb-f),norm(out@d),1e-300)
    return out,defect,ratio,bool(math.isfinite(defect))


def deflation_logs(Y: Any, roots: Sequence[np.ndarray], alpha: float) -> np.ndarray:
    yy = np.asarray(Y, np.complex128)
    if yy.ndim == 1: yy = yy[None, :]
    out = np.zeros(len(yy))
    for root in roots:
        out += np.log1p(alpha/np.maximum(np.linalg.norm(yy-root[None, :], axis=1), 1e-12))
    return out


def line_choice(R: Any, Y: Any, y: Any, current: float, best: float,
                roots: Sequence[np.ndarray], alpha: float) -> tuple[Optional[int], np.ndarray]:
    rr = np.asarray(R, float); logw = deflation_logs(Y, roots, alpha) if roots else np.zeros(len(rr))
    logw0 = float(deflation_logs(np.asarray(y)[None, :], roots, alpha)[0]) if roots else 0
    merit = np.log(np.maximum(rr, 1e-300))+logw
    admissible = np.isfinite(rr) & ((merit < math.log(max(current, 1e-300))+logw0) | (rr < best))
    return (int(np.argmin(np.where(admissible, merit, np.inf))), merit) if np.any(admissible) else (None, merit)


def probe_endpoint(target: Target, y: np.ndarray, f: np.ndarray, prev: Optional[np.ndarray],
                   args: argparse.Namespace, seed: int, ep: int) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    n = len(y); yn = max(1, norm(y)); radii = floats(args.probe_radii, [.5, 1, 2])
    gp = isinstance(target.oracle, ExactKSGPOracle); budget = min(args.probes, args.gp_probe_cap) if gp else args.probes
    base_radius = args.probe_scale*yn
    min_radius = yn*max(1e-6, math.sqrt(target.oracle.jitter/max(1, target.oracle.d))) if gp else 1e-10*yn
    local_radius = base_radius if prev is None else min(base_radius, max(2*norm(prev), min_radius))
    candidates: list[np.ndarray] = []
    if prev is not None and norm(prev) > 1e-14:
        b = y + prev
        for j in range(n):
            if abs(b[j]-y[j]) < min_radius/max(1, math.sqrt(n)):
                b[j] += min_radius/max(1, math.sqrt(n))*np.exp(2j*math.pi*(j+1)*.61803398875)
        candidates.append(b)
    k = 0
    while len(candidates) < budget:
        v = direction(n, ep*budget+k, seed+104729*(ep+1)+7919*(k+1))
        rad = local_radius*radii[k % len(radii)]
        v = v/max(norm(v), 1e-300)*math.sqrt(max(1, n))
        b = y + rad*v
        tiny = 1e-10*yn
        for j in range(n):
            if abs(b[j]-y[j]) < tiny:
                b[j] += tiny*np.exp(2j*math.pi*((j+1)*.61803398875+(k+1)*.41421356237))
        candidates.append(b); k += 1
    B = np.asarray(candidates[:budget]); FB = target.eval_batch(B); R = norms(FB)
    gap = np.linalg.norm(FB-f[None, :], axis=1)/np.maximum(R+norm(f), 1e-300)
    scores = np.log(np.maximum(R, 1e-300)) + args.equal_value_weight*np.log1p(1/np.maximum(gap, 1e-14))
    top = 1 if gp else args.probe_top
    order = np.argsort(scores)[:max(1, min(top, len(B)))]
    best = None
    for idx in order:
        try:
            Q, defect, used = finite_slope(target, y, B[int(idx)])
            s = np.linalg.svd(Q, compute_uv=False); cond = float(s[0]/s[-1]) if s[-1] > 0 else float("inf")
            total = float(scores[int(idx)]) + args.probe_cond_weight*math.log1p(cond if math.isfinite(cond) else 1e300)
            item = (total, B[int(idx)].copy(), FB[int(idx)].copy(), Q, defect, cond, float(s[-1]), float(s[0]), used)
            if best is None or item[0] < best[0]: best = item
        except Exception:
            continue
    if best is None: raise RuntimeError("no usable finite-slope endpoint")
    _, b, fb, Q, defect, cond, smin, smax, used = best
    return b, Q, {"probe_evals": len(B), "slope_evals": used, "telescope_defect": defect,
                   "cond": cond, "smin": smin, "smax": smax, "endpoint_residual": norm(fb), "endpoint_value":fb}


@dataclass
class CorrectResult:
    y: np.ndarray
    f: np.ndarray
    residual: float
    ok: bool
    status: str
    epochs: int
    slopes: int
    transported_steps: int
    probe_evals: int
    slope_evals: int
    line_evals: int
    parabolic_evals: int
    max_telescope_defect: float
    max_transport_correction: float
    smin: Optional[float]
    smax: Optional[float]
    cond: Optional[float]
    budget: int = 0
    budget_used: int = 0
    budget_exhausted: bool = False
    halley_proposals: int = 0
    halley_accepted: int = 0
    irp_proposals: int = 0
    irp_accepted: int = 0
    irp2_proposals: int = 0
    irp2_accepted: int = 0
    hirp_proposals: int = 0
    hirp_accepted: int = 0
    monomial_chain_proposals: int = 0
    monomial_chain_accepted: int = 0
    monomial_chain_degree_max: int = 0
    micro_steps_accepted: int = 0
    max_halley_factor: float = 0.0
    equation_chart_active: bool = False
    equation_chart_spread: float = 0.0
    equation_chart_center: float = 0.0


@dataclass
class CachedStart:
    """Already observed endpoint; retained only to avoid a duplicate sample."""
    y: np.ndarray
    f: np.ndarray


def _chart_log_scale(target: Target, Y: Any) -> np.ndarray:
    """Common oracle normalization removed inside one Pandrosion chart.

    Multiplying every equation value at a point by the same positive scalar
    does not change its zero set.  Freezing that scalar across a telescopic
    path restores a coherent root-equivalent target without coefficients,
    derivatives, or a second solver.
    """
    yy = np.asarray(Y, np.complex128)
    if yy.ndim == 1:
        yy = yy[None, :]
    z = np.asarray(target.to_base(yy), np.complex128)
    oracle = target.oracle
    if not isinstance(oracle, KostlanOracle):
        return np.zeros(len(z), float)
    if not oracle.feature_mode:
        return oracle.d*np.log(np.maximum(1.0, np.max(np.abs(z), axis=1)))
    with np.errstate(all="ignore"):
        amp = np.log(np.maximum(np.abs(z), 1e-300)) @ oracle.exps.T
        amp += oracle.feature_logs[None, :]
        out = np.max(np.where(np.isfinite(amp), amp, -np.inf), axis=1)
    return np.where(np.isfinite(out), out, 0.0)


def _values_in_chart(target: Target, Y: Any, reference: float) -> np.ndarray:
    yy = np.asarray(Y, np.complex128)
    if yy.ndim == 1:
        yy = yy[None, :]
    values = target.eval_batch(yy)
    exponent = np.clip(_chart_log_scale(target, yy)-reference, -700.0, 700.0)
    with np.errstate(all="ignore"):
        return np.asarray(values*np.exp(exponent)[:, None], np.complex128)


def _telescope_delta(x: np.ndarray, radius: np.ndarray,
                     epoch: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Return the deterministic endpoint and coordinate increments of a telescope.

    The helper is shared by the ordinary finite-slope layer and the singular
    monomial-chain chart.  Keeping one implementation prevents a subtle phase
    mismatch when the latter decodes degree information from the already
    observed telescope columns.
    """
    n = len(x)
    b = x + np.asarray(radius, float)*direction(n, epoch, seed)/math.sqrt(max(1, n))
    dz = b-x
    tiny = 1e-14*np.maximum(1.0, np.maximum(np.abs(x), np.abs(b)))
    for j in range(n):
        if abs(dz[j]) <= tiny[j]:
            dz[j] = tiny[j]*np.exp(2j*math.pi*((j+1)*.61803398875+(epoch+1)*.41421356237))
            b[j] = x[j]+dz[j]
    return np.asarray(b, np.complex128), np.asarray(dz, np.complex128)


def _pure_telescope(target: Target, x: np.ndarray, fx: np.ndarray, radius: np.ndarray,
                    epoch: int, seed: int) -> tuple[np.ndarray, float, float, int]:
    """One exact coordinate-telescopic slope in a frozen root-equivalent chart."""
    n = len(x)
    b, dz = _telescope_delta(x, radius, epoch, seed)
    path = np.repeat(x[None, :], n+1, axis=0)
    for j in range(n):
        path[j+1:, j] = b[j]
    reference = float(_chart_log_scale(target, x)[0])
    F = np.vstack((fx[None, :], _values_in_chart(target, path[1:], reference)))
    # Columns are independent finite slopes.  The vectorized assignment keeps
    # exactly the same adjacent differences and scalar divisions as 367.
    Q = ((F[1:]-F[:-1])/dz[:, None]).T.copy()
    with np.errstate(all="ignore"):
        projected = Q@dz
    defect = norm((F[-1]-F[0])-projected)/max(
        norm(F[-1]-F[0]), norm(projected), 1e-300
    )
    return np.asarray(Q, np.complex128), reference, defect, n


def _lu_factor(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Small standalone complex LU factorization with partial pivoting."""
    LU = np.asarray(A, np.complex128).copy()
    n = LU.shape[0]
    piv = np.arange(n)
    for k in range(n-1):
        p = k + int(np.argmax(np.abs(LU[k:, k])))
        if abs(LU[p, k]) <= 1e-15*max(1.0, float(np.max(np.abs(LU)))):
            raise np.linalg.LinAlgError("singular equilibrated telescope")
        if p != k:
            LU[[k, p]] = LU[[p, k]]
            piv[[k, p]] = piv[[p, k]]
        LU[k+1:, k] /= LU[k, k]
        LU[k+1:, k+1:] -= LU[k+1:, k, None]*LU[k, k+1:][None, :]
    if n and abs(LU[-1, -1]) <= 1e-15*max(1.0, float(np.max(np.abs(LU)))):
        raise np.linalg.LinAlgError("singular equilibrated telescope")
    return LU, piv


def _lu_solve(LU: np.ndarray, piv: np.ndarray, b: np.ndarray) -> np.ndarray:
    y = np.asarray(b, np.complex128)[piv].copy()
    n = len(y)
    for i in range(n):
        if i: y[i] -= LU[i, :i]@y[:i]
    for i in range(n-1, -1, -1):
        if i+1 < n: y[i] -= LU[i, i+1:]@y[i+1:]
        y[i] /= LU[i, i]
    return y


@dataclass
class _SlopeFactor:
    LU: np.ndarray
    piv: np.ndarray
    row: np.ndarray
    col: np.ndarray
    # Reused triangular-solve workspace.  It removes one complex allocation per
    # slope solve without changing pivoting, arithmetic order, or candidates.
    work: np.ndarray


def _factor_slope(Q: np.ndarray) -> Optional[_SlopeFactor]:
    """Two-sided diagonal equilibration followed by one reusable LU."""
    A = np.asarray(Q, np.complex128)
    if not A.size or not np.all(np.isfinite(A)):
        return None
    row_norm = np.max(np.abs(A), axis=1)
    row = 1.0/np.maximum(row_norm, 1e-300)
    Ar = row[:, None]*A
    col_norm = np.max(np.abs(Ar), axis=0)
    col = 1.0/np.maximum(col_norm, 1e-300)
    Ae = Ar*col[None, :]
    try:
        LU, piv = _lu_factor(Ae)
    except Exception:
        return None
    return _SlopeFactor(LU, piv, row, col, np.empty(A.shape[0], np.complex128))


def _factored_slope_step(factor: Optional[_SlopeFactor], f: np.ndarray,
                         x: np.ndarray, cap_factor: float) -> tuple[np.ndarray, str]:
    if factor is None:
        return np.zeros_like(x), "failed"
    try:
        # Same forward/back substitution as 367, but in a factor-owned buffer.
        # The resulting ``step`` is a fresh array, so the workspace can safely
        # be reused by the later IRP residual solve in the same epoch.
        rhs = -factor.row*np.asarray(f, np.complex128)
        y = factor.work
        np.take(rhs, factor.piv, out=y)
        n = len(y)
        for i in range(n):
            if i:
                y[i] -= factor.LU[i, :i]@y[:i]
        for i in range(n-1, -1, -1):
            if i+1 < n:
                y[i] -= factor.LU[i, i+1:]@y[i+1:]
            y[i] /= factor.LU[i, i]
        step = factor.col*y
    except Exception:
        return np.zeros_like(x), "failed"
    cap = cap_factor*max(1.0, norm(x))
    size = norm(step)
    if not math.isfinite(size):
        return np.zeros_like(x), "failed"
    if size > cap:
        step = step*(cap/size)
    return np.asarray(step, np.complex128), "equilibrated-reused-lu"


def _adaptive_pandrosion_radius(target: Target, x: np.ndarray,
                                previous_step: Optional[np.ndarray],
                                singular_rebuilds: int,
                                args: argparse.Namespace) -> np.ndarray:
    """Componentwise scale-invariant telescope radii.

    Each coordinate receives a radius proportional to its own current scale and
    recent displacement.  The geometric mean is normalized back to the global
    requested radius, preventing one large-unit coordinate from dominating the
    entire telescope.
    """
    n = len(x)
    growth = 10.0**min(singular_rebuilds, 4)
    coord = np.maximum(1.0, np.abs(x))
    if previous_step is not None:
        ps = np.abs(np.asarray(previous_step, np.complex128))
        coord = np.maximum(coord, ps/max(float(args.pandrosion_radius_step_ratio), 1e-12))
    # Keep coordinate scaling relative, not explosively absolute.
    gm = math.exp(float(np.mean(np.log(np.maximum(coord, 1e-300)))))
    relative = np.clip(coord/max(gm, 1e-300), 1e-4, 1e4)
    global_scale = 1.0 if args.pandrosion_radius_mode == "fixed" else max(1.0, norm(x)/math.sqrt(max(1, n)))
    base = float(args.pandrosion_radius)*growth*global_scale
    floor = math.sqrt(np.finfo(float).eps)*coord
    if singular_rebuilds > 0 and args.degree_radius_chart and int(target.oracle.d) > 1:
        # A monomial response of degree d has finite slope O(r^(d-1)) near a
        # flat chart.  Choose r so this response is at least sqrt(eps), then
        # cap the chart width to avoid a violent high-degree excursion.
        d = int(target.oracle.d)
        degree_cell = min(float(args.degree_radius_cap),
                          np.finfo(float).eps**(1.0/(2.0*max(1, d-1))))
        floor = np.maximum(floor, degree_cell*coord)
    if isinstance(target.oracle, ExactKSGPOracle):
        floor = np.maximum(floor, float(target.oracle.recommended_radius))
    radii = np.maximum(floor, base*relative)
    if previous_step is not None:
        ceiling = np.maximum(floor, float(args.pandrosion_radius_step_ratio)*np.maximum(np.abs(previous_step), floor))
        radii = np.minimum(radii, ceiling)
    return np.asarray(radii, float)

def _cap_step(step: np.ndarray, x: np.ndarray, cap_factor: float) -> np.ndarray:
    out = np.asarray(step, np.complex128).copy()
    size = norm(out); cap = cap_factor*max(1.0, norm(x))
    if math.isfinite(size) and size > cap:
        out *= cap/size
    return out


def _directional_halley_factor(f0: np.ndarray, fhalf: np.ndarray,
                               ffull: np.ndarray,
                               args: argparse.Namespace) -> tuple[Optional[float], dict[str, float]]:
    """Derivative-free Halley factor on the projected Pandrosion ray.

    For g(t)=Re <F(x),F(x+t p)>/||F(x)||^2, use the observations at
    t=0,1/2,1 to form second-order forward differences at zero, then apply the
    scalar Halley formula.  The proposal is only a globalization candidate;
    the full vector residual decides whether it is accepted.
    """
    scale = float(np.max(np.abs(f0))) if len(f0) else 0.0
    if not math.isfinite(scale) or scale <= 1e-300:
        return None, {"g_half": float("nan"), "g_full": float("nan"),
                      "g_prime": float("nan"), "g_second": float("nan")}
    a = np.asarray(f0/scale, np.complex128)
    b = np.asarray(fhalf/scale, np.complex128)
    c = np.asarray(ffull/scale, np.complex128)
    den = float(np.vdot(a, a).real)
    if not math.isfinite(den) or den <= 1e-300:
        return None, {"g_half": float("nan"), "g_full": float("nan"),
                      "g_prime": float("nan"), "g_second": float("nan")}
    g0 = 1.0
    gh = float(np.vdot(a, b).real/den)
    g1 = float(np.vdot(a, c).real/den)
    gp = -3.0*g0 + 4.0*gh - g1
    gpp = 4.0*(g0 - 2.0*gh + g1)
    hden = 2.0*gp*gp - g0*gpp
    diag = {"g_half": gh, "g_full": g1, "g_prime": gp, "g_second": gpp}
    if (not args.directional_halley or not all(math.isfinite(v) for v in (gh, g1, gp, gpp, hden))
            or gp >= -float(args.halley_min_projected_slope)
            or abs(hden) <= 1e-14*max(1.0, gp*gp, abs(gpp))):
        return None, diag
    factor = -2.0*g0*gp/hden
    if not math.isfinite(factor):
        return None, diag
    factor = float(np.clip(factor, args.halley_min_factor, args.halley_max_factor))
    return factor, diag


def _distinct_step(step: np.ndarray, existing: Sequence[np.ndarray], tol: float = 1e-10) -> bool:
    """Scale-aware step distinctness.

    The previous absolute unit floor treated all corrections below ``1e-10``
    as identical.  That is inappropriate in reciprocal/small-root charts, for
    example when the root itself is of size ``1e-8``.  Distinctness is now
    measured relative to the corrections being compared.  The zero/zero case
    remains non-distinct.
    """
    a = norm(step)
    for other in existing:
        b = norm(other)
        scale = max(a, b, np.finfo(float).tiny)
        if norm(step-other) <= tol*scale:
            return False
    return True


def _select_equation_chart(Q: np.ndarray, f: np.ndarray, x: np.ndarray,
                           radius: np.ndarray, args: argparse.Namespace
                           ) -> tuple[np.ndarray, dict[str, float | bool]]:
    """Select one frozen diagonal equation homothety.

    The chart is inferred primarily from the *observed equation values at the
    common start*.  This is deliberate: arbitrary equation-unit multipliers
    appear directly in those values, whereas extreme variable units need not.
    A tiny relative floor handles components that happen to vanish at the
    start without letting a large telescope cell masquerade as equation scale.
    """
    del Q, x, radius
    ff = np.asarray(f, np.complex128)
    magnitude = np.abs(ff)
    positive = magnitude[np.isfinite(magnitude) & (magnitude > 0)]
    reference = float(np.median(positive)) if len(positive) else 1.0
    floor = max(1e-3*reference, 1e-300)
    scale = np.maximum(magnitude, floor)
    logs = np.log10(scale)
    finite = np.isfinite(logs)
    if not np.any(finite):
        return np.ones(len(ff), float), {"active": False, "spread": 0.0, "center": 0.0}
    good = logs[finite]
    spread = float(np.max(good)-np.min(good)) if len(good) else 0.0
    center = float(np.median(good)) if len(good) else 0.0
    active = bool(args.equation_chart and (
        spread >= float(args.equation_chart_spread)
        or (len(ff) == 1 and abs(center) >= float(args.equation_chart_center))
    ))
    if not active:
        return np.ones(len(ff), float), {"active": False, "spread": spread, "center": center}
    weights = np.clip(1.0/scale, 1e-150, 1e150)
    return np.asarray(weights, float), {"active": True, "spread": spread, "center": center}


def _apply_equation_chart(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    v = np.asarray(values, np.complex128)
    return np.asarray(v*np.asarray(weights, float), np.complex128)


def _palette_cap_residual_layer(layer: np.ndarray, anchor: np.ndarray,
                                x: np.ndarray, args: argparse.Namespace
                                ) -> tuple[np.ndarray, float]:
    """Dyadic logarithmic palette for an additive residual layer.

    The exact frozen-chart residual solve is returned unchanged in the local
    cell.  If it is too large relative to the preceding layer, it is snapped
    down by a power of the requested palette base.  This is the additive
    finite-slope analogue of the paper's logarithmic residual-scale palette.
    """
    q = np.asarray(layer, np.complex128).copy()
    qn = norm(q)
    an = max(norm(anchor), 1e-14*max(1.0, norm(x)))
    ratio = qn/max(an, 1e-300)
    cap = float(args.residual_layer_ratio_cap)
    if not math.isfinite(ratio) or ratio <= cap or qn <= 0:
        return q, 1.0
    base = max(float(args.residual_palette_base), 1.0000001)
    exponent = math.floor(math.log(max(cap/ratio, 1e-300), base))
    factor = min(1.0, base**exponent)
    factor = max(factor, float(args.residual_palette_min_factor))
    return np.asarray(factor*q, np.complex128), float(factor)


def _reclosed_residual_layer(Q: np.ndarray, x0: np.ndarray, f0: np.ndarray,
                             x1: np.ndarray, f1: np.ndarray,
                             anchor: np.ndarray, args: argparse.Namespace
                             ) -> tuple[Optional[np.ndarray], Optional[np.ndarray], dict[str, float | bool]]:
    """Reclose a finite-slope chart and solve one residual Pandrosion layer."""
    try:
        Qt, defect, correction, ok = close_finite_slope(
            np.asarray(Q, np.complex128), np.asarray(x1, np.complex128),
            np.asarray(f1, np.complex128), np.asarray(x0, np.complex128),
            np.asarray(f0, np.complex128), float(args.residual_transport_limit),
        )
    except Exception:
        return None, None, {"ok": False, "defect": float("inf"),
                            "correction": float("inf"), "palette": 0.0}
    if not ok:
        return None, None, {"ok": False, "defect": defect,
                            "correction": correction, "palette": 0.0}
    factor = _factor_slope(Qt)
    q, mode = _factored_slope_step(
        factor, np.asarray(f1, np.complex128), np.asarray(x1, np.complex128),
        float(args.pandrosion_step_cap),
    )
    if mode == "failed" or not np.all(np.isfinite(q)):
        return None, None, {"ok": False, "defect": defect,
                            "correction": correction, "palette": 0.0}
    q, palette = _palette_cap_residual_layer(q, anchor, x1, args)
    return q, Qt, {"ok": True, "defect": defect,
                   "correction": correction, "palette": palette}


def _raw_chart_values(target: Target, Y: np.ndarray, reference: float,
                      equation_weights: np.ndarray
                      ) -> tuple[np.ndarray, np.ndarray]:
    raw = target.eval_batch(Y)
    factors = np.exp(np.clip(_chart_log_scale(target, Y)-reference, -700.0, 700.0))
    with np.errstate(all="ignore"):
        common = np.asarray(raw*factors[:, None], np.complex128)
    balanced = _apply_equation_chart(common, equation_weights[None, :])
    return np.asarray(raw, np.complex128), np.asarray(balanced, np.complex128)


def _monomial_chain_hint(x: np.ndarray, f: np.ndarray,
                         degree: int, args: argparse.Namespace) -> bool:
    """Detect a flat, residual-concentrated polynomial chart.

    The intended case is a chain of proportional monomial closures started on
    a singular affine chart, such as

    ``z_j**p - c_j z_{j+1} = 0`` and ``z_n**p-c_n = 0`` at the origin.

    The gate is deliberately narrow.  It does not inspect family names and it
    remains inactive on ordinary dense systems, including the zero starts of
    the Kostlan suite, whose residual energy is spread across many equations.
    """
    if (not getattr(args, "monomial_chain_chart", True)
            or len(x) < 2 or int(degree) < 2):
        return False
    ff = np.asarray(f, np.complex128)
    metric = norm(ff)
    if (not math.isfinite(metric)
            or metric < float(args.monomial_chain_residual_min)):
        return False
    concentration = float(np.max(np.abs(ff))/max(metric, 1e-300))
    rms_x = norm(x)/math.sqrt(max(1, len(x)))
    flat_scale = max(1.0, rms_x)
    flat_fraction = float(np.mean(
        np.abs(np.asarray(x, np.complex128))
        <= float(args.monomial_chain_flat_tolerance)*flat_scale
    ))
    return bool(
        concentration >= float(args.monomial_chain_concentration)
        and flat_fraction >= float(args.monomial_chain_flat_fraction)
    )


def _band_orientation(Q: np.ndarray) -> tuple[str, float]:
    """Return the closest diagonal-plus-one-neighbour chain orientation."""
    A = np.asarray(Q, np.complex128)
    n = A.shape[0]
    den = max(norm(A), 1e-300)

    upper = A.copy()
    lower = A.copy()
    for j in range(n):
        upper[j, j] = 0.0
        lower[j, j] = 0.0
    for j in range(n-1):
        upper[j, j+1] = 0.0
        lower[j+1, j] = 0.0
    eu = norm(upper)/den
    el = norm(lower)/den
    return ("upper", eu) if eu <= el else ("lower", el)


def _complex_pth_roots(value: complex, degree: int) -> list[complex]:
    p = max(1, int(degree))
    if abs(value) == 0.0:
        return [0j]*p
    mag = abs(value)**(1.0/p)
    angle = math.atan2(value.imag, value.real)
    return [mag*np.exp(1j*(angle+2.0*math.pi*k)/p) for k in range(p)]


def _infer_probe_degree(xj: complex, f0: complex,
                        deltas: Sequence[complex], values: Sequence[complex],
                        max_degree: int, args: argparse.Namespace
                        ) -> tuple[Optional[int], Optional[complex], dict[str, float]]:
    """Infer the effective monomial degree from two homothetic probes.

    For every candidate ``p`` we compare the two coefficients obtained from

    ``(F(x+delta)-F(x))/((x+delta)^p-x^p)``.

    This is more reliable than trusting the declared maximum polynomial degree:
    the benchmark's cyclic family, for example, advertises degree five while
    individual instances use powers two, three, or five.
    """
    best: Optional[tuple[float, int, complex, float]] = None
    scale = max(abs(f0), 1e-300)
    for p in range(2, max(2, int(max_degree))+1):
        phi = [((xj+d)**p-xj**p) for d in deltas]
        if min(abs(z) for z in phi) <= 1e-300:
            continue
        coeff = [(values[k]-f0)/phi[k] for k in range(2)]
        if not all(np.isfinite(z.real) and np.isfinite(z.imag) for z in coeff):
            continue
        mismatch = abs(coeff[0]-coeff[1])/max(
            abs(coeff[0]), abs(coeff[1]), 1e-300
        )
        signal = min(abs(values[k]-f0) for k in range(2))/scale
        score = mismatch
        if best is None or score < best[0]:
            # Weight the larger probe slightly more; it has the better
            # floating-point signal for high powers.
            a = (2.0*coeff[0]+coeff[1])/3.0
            best = (score, p, a, signal)
    if best is None:
        return None, None, {"mismatch": float("inf"), "signal": 0.0}
    mismatch, p, a, signal = best
    ok = bool(
        mismatch <= float(args.monomial_chain_degree_mismatch)
        and signal >= float(args.monomial_chain_probe_signal)
        and abs(a) > 1e-300
    )
    return (p if ok else None), (a if ok else None), {
        "mismatch": float(mismatch), "signal": float(signal),
        "degree": float(p),
    }


def _monomial_chain_portfolio(
    target: Target,
    x: np.ndarray,
    fx: np.ndarray,
    fx_bal: np.ndarray,
    Q_bal: np.ndarray,
    radius: np.ndarray,
    epoch: int,
    seed: int,
    reference: float,
    equation_weights: np.ndarray,
    args: argparse.Namespace,
) -> tuple[Optional[np.ndarray], Optional[np.ndarray], float, str, dict[str, float | bool]]:
    """Five-sample singular-chart portfolio for monomial chains.

    Two homothetic coordinate probes infer the effective degree and the
    terminal monomial coefficient.  If the already observed telescope is close
    to a diagonal-plus-super/subdiagonal chain, three model roots are then
    constructed recursively and evaluated.  The routine never accepts its
    model prediction without checking the full black-box residual.

    The five observations replace, rather than augment, the normal five line
    samples, so the strict ``n+5`` epoch budget is preserved.
    """
    n = len(x)
    orientation, band_error = _band_orientation(Q_bal)
    diag = {
        "reliable": False,
        "samples": 0.0,
        "band_error": float(band_error),
        "degree": 0.0,
        "degree_mismatch": float("inf"),
        "probe_signal": 0.0,
        "terminal_coefficient_mismatch": float("inf"),
    }
    terminal = n-1 if orientation == "upper" else 0
    max_component = float(np.max(np.abs(fx_bal))) if n else 0.0
    if (band_error > float(args.monomial_chain_band_tolerance)
            or abs(fx_bal[terminal])
            < float(args.monomial_chain_terminal_fraction)*max(max_component, 1e-300)):
        return None, None, float("inf"), "chain-unreliable", diag

    probe_radius = float(args.monomial_chain_probe_radius)*max(1.0, abs(x[terminal]))
    phase = np.exp(2j*math.pi*((terminal+1)*.61803398875+(epoch+1)*.41421356237))
    deltas = [probe_radius*phase, 0.5*probe_radius*phase]
    probe_points = []
    for delta in deltas:
        y = np.asarray(x, np.complex128).copy()
        y[terminal] += delta
        probe_points.append(y)
    probe_Y = np.vstack(probe_points)
    probe_raw, probe_bal = _raw_chart_values(
        target, probe_Y, reference, equation_weights
    )
    diag["samples"] = 2.0

    pmax = min(int(target.oracle.d), int(args.monomial_chain_degree_max))
    degree, probe_coefficient, pdiag = _infer_probe_degree(
        complex(x[terminal]), complex(fx_bal[terminal]), deltas,
        [complex(probe_bal[0, terminal]), complex(probe_bal[1, terminal])],
        pmax, args,
    )
    diag["degree_mismatch"] = float(pdiag["mismatch"])
    diag["probe_signal"] = float(pdiag["signal"])
    diag["degree"] = float(pdiag.get("degree", 0.0))

    # If the degree fit fails, spend the remaining three samples on phase-
    # separated chart probes.  This keeps the epoch budget exact and can still
    # escape a pathological cancellation, but it is not counted as a decoded
    # monomial-chain proposal.
    if degree is None or probe_coefficient is None:
        extra_points = []
        for k in (1, 2, 3):
            y = np.asarray(x, np.complex128).copy()
            y[terminal] += probe_radius*phase*np.exp(0.5j*math.pi*k)
            extra_points.append(y)
        extra_raw, extra_bal = _raw_chart_values(
            target, np.vstack(extra_points), reference, equation_weights
        )
        diag["samples"] = 5.0
        all_points = np.vstack((probe_Y, np.vstack(extra_points)))
        all_raw = np.vstack((probe_raw, extra_raw))
        all_bal = np.vstack((probe_bal, extra_bal))
        metrics = norms(all_bal)
        k = int(np.argmin(metrics))
        return all_points[k], all_raw[k], float(metrics[k]), "chain-probe", diag

    _, dz = _telescope_delta(x, radius, epoch, seed)
    p = int(degree)
    coefficients = np.zeros(n, np.complex128)
    for j in range(n):
        phi = ((x[j]+dz[j])**p-x[j]**p)/dz[j]
        if abs(phi) <= 1e-300:
            return None, None, float("inf"), "chain-unreliable", diag
        coefficients[j] = Q_bal[j, j]/phi
    if (not np.all(np.isfinite(coefficients))
            or np.min(np.abs(coefficients)) <= 1e-300):
        return None, None, float("inf"), "chain-unreliable", diag

    terminal_mismatch = abs(coefficients[terminal]-probe_coefficient)/max(
        abs(coefficients[terminal]), abs(probe_coefficient), 1e-300
    )
    diag["terminal_coefficient_mismatch"] = float(terminal_mismatch)
    if terminal_mismatch > float(args.monomial_chain_coefficient_mismatch):
        return None, None, float("inf"), "chain-unreliable", diag

    constants = np.asarray(fx_bal, np.complex128).copy()
    constants -= coefficients*np.asarray(x, np.complex128)**p
    if orientation == "upper":
        coupling = np.asarray([Q_bal[j, j+1] for j in range(n-1)], np.complex128)
        constants[:-1] -= coupling*np.asarray(x[1:], np.complex128)
    else:
        coupling = np.asarray([Q_bal[j+1, j] for j in range(n-1)], np.complex128)
        constants[1:] -= coupling*np.asarray(x[:-1], np.complex128)

    # Three deterministic branch policies.  Every policy solves the decoded
    # triangular model exactly; the full oracle selects the useful branch.
    policies: list[tuple[str, int]] = [
        ("fixed", 0),
        ("fixed", max(1, int(round(p/3.0))) % p),
        ("nearest", 0),
        ("alternating", 1),
        ("fixed", max(1, int(round(2.0*p/3.0))) % p),
    ]
    chain_points: list[np.ndarray] = []

    def choose_root(value: complex, coordinate: int,
                    policy: tuple[str, int]) -> complex:
        roots = _complex_pth_roots(value, p)
        kind, offset = policy
        if kind == "fixed":
            return roots[offset % p]
        if kind == "alternating":
            return roots[(offset+coordinate) % p]
        target_value = x[coordinate] if abs(x[coordinate]) > 1e-14 else x[coordinate]+dz[coordinate]
        return min(roots, key=lambda z: abs(z-target_value))

    for policy in policies:
        z = np.asarray(x, np.complex128).copy()
        try:
            if orientation == "upper":
                value = -constants[-1]/coefficients[-1]
                z[-1] = choose_root(value, n-1, policy)
                for j in range(n-2, -1, -1):
                    value = -(coupling[j]*z[j+1]+constants[j])/coefficients[j]
                    z[j] = choose_root(value, j, policy)
            else:
                value = -constants[0]/coefficients[0]
                z[0] = choose_root(value, 0, policy)
                for j in range(1, n):
                    value = -(coupling[j-1]*z[j-1]+constants[j])/coefficients[j]
                    z[j] = choose_root(value, j, policy)
        except Exception:
            continue
        step = _cap_step(z-x, x, float(args.pandrosion_step_cap))
        z = np.asarray(x+step, np.complex128)
        if all(norm(z-old) > 1e-11*max(norm(z), norm(old), 1e-300)
               for old in chain_points):
            chain_points.append(z)
        if len(chain_points) == 3:
            break
    while len(chain_points) < 3:
        # Rare duplicate-branch case (notably p=2): use an alternating branch
        # with a shifted parity.  It still solves the decoded chain model.
        policy = ("alternating", len(chain_points)+1)
        z = np.asarray(x, np.complex128).copy()
        if orientation == "upper":
            z[-1] = choose_root(-constants[-1]/coefficients[-1], n-1, policy)
            for j in range(n-2, -1, -1):
                z[j] = choose_root(
                    -(coupling[j]*z[j+1]+constants[j])/coefficients[j], j, policy
                )
        else:
            z[0] = choose_root(-constants[0]/coefficients[0], 0, policy)
            for j in range(1, n):
                z[j] = choose_root(
                    -(coupling[j-1]*z[j-1]+constants[j])/coefficients[j], j, policy
                )
        chain_points.append(np.asarray(x+_cap_step(
            z-x, x, float(args.pandrosion_step_cap)
        ), np.complex128))

    chain_Y = np.vstack(chain_points[:3])
    chain_raw, chain_bal = _raw_chart_values(
        target, chain_Y, reference, equation_weights
    )
    diag["samples"] = 5.0
    all_points = np.vstack((probe_Y, chain_Y))
    all_raw = np.vstack((probe_raw, chain_raw))
    all_bal = np.vstack((probe_bal, chain_bal))
    metrics = norms(all_bal)
    k = int(np.argmin(metrics))
    diag["reliable"] = True
    label = "monomial-chain" if k >= 2 else "chain-probe"
    return all_points[k], all_raw[k], float(metrics[k]), label, diag


def correct(target: Target, y0: Any, args: argparse.Namespace, known: Sequence[np.ndarray],
            epochs: int, deadline: Optional[float], cached: Optional[CachedStart]=None) -> CorrectResult:
    """Residual-chart Pandrosion with sequential IRP-2 and HIRP-1.

    One exact coordinate telescope supplies the raw layer.  The two fixed pilot
    samples at half/full step diagnose the local residual chart.  In a stable
    cell, the remaining three samples execute either a genuine Halley-then-IRP
    composition or two genuinely renormalized raw residual layers.  Outside
    that cell the robust 369-style dyadic portfolio is retained.
    """
    del known
    n = target.oracle.n
    call_start = int(target.oracle.eval_count)
    budget = int(args.strict_budget_factor*n+args.strict_budget_base)
    if cached is None:
        x = np.asarray(y0, np.complex128).copy()
        fx = target.eval(x)
    else:
        x = np.asarray(cached.y, np.complex128).copy()
        fx = np.asarray(cached.f, np.complex128).copy()

    if norm(fx) == 0.0:
        return CorrectResult(
            y=x.copy(), f=fx.copy(), residual=0.0, ok=True,
            status="converged", epochs=0, slopes=0, transported_steps=0,
            probe_evals=0, slope_evals=0, line_evals=0,
            parabolic_evals=0, max_telescope_defect=0.0,
            max_transport_correction=0.0, smin=None, smax=None, cond=None,
            budget=budget, budget_used=int(target.oracle.eval_count)-call_start,
        )

    best_x, best_f, best_r = x.copy(), fx.copy(), float("inf")
    raw_grid = floats(args.line_grid, [1.0, .5, .25, .125, .0625])
    raw_fill = [t for t in raw_grid if abs(t-1.0) > 1e-14 and abs(t-.5) > 1e-14]
    raw_fill.extend([.25, .125, .0625, 1.5, 2.0])
    seed = mix64(int(target.oracle.seed)+0x371A71A5)
    slopes = slope_evals = line_evals = done = 0
    singular_rebuilds = 0
    transported_steps = 0
    max_defect = 0.0
    max_transport_correction = 0.0
    status = "max-epochs"
    last_Q: Optional[np.ndarray] = None
    max_norm = float(args.max_solution_norm)
    previous_step: Optional[np.ndarray] = None
    equation_weights: Optional[np.ndarray] = None
    equation_chart_active = False
    equation_chart_spread = 0.0
    equation_chart_center = 0.0
    halley_proposals = halley_accepted = 0
    irp_proposals = irp_accepted = 0
    irp2_proposals = irp2_accepted = 0
    hirp_proposals = hirp_accepted = 0
    monomial_chain_proposals = monomial_chain_accepted = 0
    monomial_chain_degree_max = 0
    micro_steps_accepted = 0
    monomial_chain_disabled = False
    max_halley_factor = 0.0

    # The solver oracle may be expressed in arbitrary equation units while an
    # independent validator uses the unscaled equations.  Keep a modest
    # internal safety margin before terminating, but retain ``args.accept`` as
    # the final mathematical acceptance threshold.  If the strict budget is
    # reached, the best already admissible point is still returned.
    stop_threshold = float(args.accept)*float(args.internal_accept_factor)

    for ep in range(max(1, int(epochs))):
        if equation_weights is not None and best_r < stop_threshold:
            status = "converged"
            break
        if deadline is not None and clock() >= deadline:
            status = "timeout"
            break
        used = int(target.oracle.eval_count)-call_start
        if used+n+5 > budget:
            status = "strict-budget-exhausted"
            break

        chain_hint = bool(
            (not monomial_chain_disabled)
            and _monomial_chain_hint(x, fx, int(target.oracle.d), args)
        )
        radius = _adaptive_pandrosion_radius(
            target, x, previous_step, singular_rebuilds, args
        )
        if chain_hint:
            # The ordinary tiny finite-slope cell erases a p-th power response
            # at a flat start.  Move only the telescope chart—not the iterate—
            # to a constructible monomial scale where the diagonal powers are
            # visible above floating-point cancellation.
            radius = np.maximum(
                radius,
                float(args.monomial_chain_telescope_radius)
                *np.maximum(1.0, np.abs(np.asarray(x, np.complex128))),
            )
        try:
            Q, reference, defect, slope_used = _pure_telescope(
                target, x, fx, radius, ep, seed,
            )
        except Exception as exc:
            status = f"telescope-error:{type(exc).__name__}"
            break
        last_Q = Q
        slopes += 1
        slope_evals += slope_used
        max_defect = max(max_defect, defect)

        if equation_weights is None:
            equation_weights, eqdiag = _select_equation_chart(
                Q, fx, x, radius, args
            )
            equation_chart_active = bool(eqdiag["active"])
            equation_chart_spread = float(eqdiag["spread"])
            equation_chart_center = float(eqdiag["center"])
            current_metric = norm(_apply_equation_chart(fx, equation_weights))
            best_x, best_f, best_r = x.copy(), fx.copy(), current_metric
            if current_metric < stop_threshold:
                status = "converged"
                done = ep
                break

        assert equation_weights is not None
        fx_bal = _apply_equation_chart(fx, equation_weights)
        Q_bal = np.asarray(Q*equation_weights[:, None], np.complex128)

        if chain_hint:
            chain_x, chain_f, chain_r, chain_label, chaindiag = _monomial_chain_portfolio(
                target, x, fx, fx_bal, Q_bal, radius, ep, seed,
                reference, equation_weights, args,
            )
            chain_samples = int(round(float(chaindiag.get("samples", 0.0))))
            line_evals += chain_samples
            if bool(chaindiag.get("reliable", False)):
                monomial_chain_proposals += 1
                monomial_chain_degree_max = max(
                    monomial_chain_degree_max,
                    int(round(float(chaindiag.get("degree", 0.0)))),
                )
            current_chain_metric = max(norm(fx_bal), 1e-300)
            improved = bool(
                chain_x is not None and chain_f is not None
                and math.isfinite(chain_r)
                and chain_r <= current_chain_metric
                    *(1.0-float(args.monomial_chain_min_improvement))
            )
            if improved:
                previous_x = x.copy()
                x = np.asarray(chain_x, np.complex128).copy()
                fx = np.asarray(chain_f, np.complex128).copy()
                previous_step = np.asarray(x-previous_x, np.complex128)
                singular_rebuilds = 0
                if bool(chaindiag.get("reliable", False)):
                    monomial_chain_accepted += 1
                r = float(chain_r)
                if r < best_r and norm(x) <= max_norm:
                    best_x, best_f, best_r = x.copy(), fx.copy(), r
                done = ep+1
                continue
            if chain_samples:
                # A failed decoded chart is not retried indefinitely.  The
                # next epoch returns to the ordinary HIRP/dyadic mechanism.
                monomial_chain_disabled = True
                singular_rebuilds += 1
                previous_step = None
                done = ep+1
                continue

        slope_factor = _factor_slope(Q_bal)
        step, mode = _factored_slope_step(
            slope_factor, fx_bal, x, float(args.pandrosion_step_cap)
        )
        step_size = norm(step)
        step_resolution = 1e-15*max(norm(x), norm(radius), 1e-12)
        if mode == "failed" or not np.all(np.isfinite(step)) or step_size <= step_resolution:
            status = "singular-telescope-rebuild"
            singular_rebuilds += 1
            previous_step = None
            done = ep+1
            continue
        singular_rebuilds = 0

        # Fixed pilot observations: two raw candidates and the complete local
        # data needed for the directional Halley/linearity diagnostics.
        pilot_steps = [0.5*step, step]
        pilot_Y = np.vstack([x+s for s in pilot_steps])
        pilot_raw, pilot_bal = _raw_chart_values(
            target, pilot_Y, reference, equation_weights
        )
        denom = max(norm(fx_bal), norm(pilot_bal[0]), norm(pilot_bal[1]), 1e-300)
        linearity_defect = norm(fx_bal-2.0*pilot_bal[0]+pilot_bal[1])/denom
        current_metric = max(norm(fx_bal), 1e-300)
        full_ratio = norm(pilot_bal[1])/current_metric
        half_ratio = norm(pilot_bal[0])/current_metric

        halley_factor, hdiag = _directional_halley_factor(
            fx_bal, pilot_bal[0], pilot_bal[1], args
        )
        halley_step: Optional[np.ndarray] = None
        halley_predicted = float("inf")
        if (halley_factor is not None
                and float(args.linearity_halley_min) <= linearity_defect
                <= float(args.linearity_halley_max)):
            t = float(halley_factor)
            predicted = 1.0+hdiag["g_prime"]*t+0.5*hdiag["g_second"]*t*t
            pilot_best = min(abs(hdiag["g_half"]), abs(hdiag["g_full"]), 1.0)
            if (math.isfinite(predicted)
                    and abs(predicted) <= float(args.halley_prediction_ratio)*pilot_best):
                proposal = _cap_step(t*step, x, float(args.pandrosion_step_cap))
                if _distinct_step(proposal, pilot_steps):
                    halley_step = proposal
                    halley_predicted = abs(float(predicted))

        # Decide whether the remaining three samples should execute a true
        # sequential residual composition.  These gates are intentionally
        # stricter than the one-shot 369 accelerator gates.
        use_hirp = bool(
            args.hirp_layer and halley_step is not None
            and linearity_defect <= float(args.hirp_linearity_max)
            and full_ratio <= float(args.hirp_full_ratio)
            and half_ratio <= float(args.hirp_half_ratio)
            and abs(float(halley_factor)-1.0) <= float(args.hirp_factor_window)
            and (n == 1 or current_metric <= float(args.hirp_residual_trigger))
        )
        use_irp2 = bool(
            (not use_hirp) and args.irp2_layer
            and linearity_defect <= float(args.irp2_linearity_max)
            and full_ratio <= float(args.irp2_full_ratio)
            and half_ratio < 1.0
        )

        extra_steps: list[np.ndarray] = []
        extra_labels: list[str] = []
        extra_raw_rows: list[np.ndarray] = []
        extra_bal_rows: list[np.ndarray] = []
        occupied: list[np.ndarray] = [pilot_steps[0], pilot_steps[1]]

        def append_dyadics(count: int) -> None:
            for factor_value in raw_fill:
                if len(extra_steps) >= count:
                    break
                proposal = _cap_step(
                    float(factor_value)*step, x, float(args.pandrosion_step_cap)
                )
                if _distinct_step(proposal, occupied):
                    extra_steps.append(proposal)
                    extra_labels.append("dyadic")
                    occupied.append(proposal)
            while len(extra_steps) < count:
                factor_value = 2.0**(-len(extra_steps)-3)
                proposal = _cap_step(
                    factor_value*step, x, float(args.pandrosion_step_cap)
                )
                extra_steps.append(proposal)
                extra_labels.append("dyadic")
                occupied.append(proposal)

        sequential = False
        if use_hirp and halley_step is not None:
            sequential = True
            # Layer 1: actual Halley homothety.
            hY = np.asarray([x+halley_step], np.complex128)
            hraw, hbal = _raw_chart_values(target, hY, reference, equation_weights)
            extra_steps.append(halley_step)
            extra_labels.append("halley")
            extra_raw_rows.append(hraw[0])
            extra_bal_rows.append(hbal[0])
            occupied.append(halley_step)
            halley_proposals += 1
            max_halley_factor = max(max_halley_factor, abs(float(halley_factor)))

            # Layer 2: recompute the residual, reclose the chart, and apply one
            # raw Pandrosion residual layer.  This is HIRP-1, not a prediction.
            qh, _, qdiag = _reclosed_residual_layer(
                Q_bal, x, fx_bal, hY[0], hbal[0], halley_step, args
            )
            if qh is not None:
                transported_steps += 1
                max_transport_correction = max(
                    max_transport_correction, float(qdiag["correction"])
                )
                hirstep = _cap_step(
                    halley_step+qh, x, float(args.pandrosion_step_cap)
                )
                if _distinct_step(hirstep, occupied):
                    hirY = np.asarray([x+hirstep], np.complex128)
                    hirraw, hirbal = _raw_chart_values(
                        target, hirY, reference, equation_weights
                    )
                    extra_steps.append(hirstep)
                    extra_labels.append("hirp")
                    extra_raw_rows.append(hirraw[0])
                    extra_bal_rows.append(hirbal[0])
                    occupied.append(hirstep)
                    hirp_proposals += 1

            append_dyadics(3)
            # Evaluate only the still-unevaluated dyadic tail.
            missing = len(extra_steps)-len(extra_raw_rows)
            if missing > 0:
                tail_steps = extra_steps[len(extra_raw_rows):]
                tail_Y = np.vstack([x+s for s in tail_steps])
                tail_raw, tail_bal = _raw_chart_values(
                    target, tail_Y, reference, equation_weights
                )
                extra_raw_rows.extend(list(tail_raw))
                extra_bal_rows.extend(list(tail_bal))

        elif use_irp2:
            # First true residual layer from the already observed full step.
            q1, Q1, q1diag = _reclosed_residual_layer(
                Q_bal, x, fx_bal, pilot_Y[1], pilot_bal[1], step, args
            )
            if q1 is not None and Q1 is not None:
                sequential = True
                transported_steps += 1
                max_transport_correction = max(
                    max_transport_correction, float(q1diag["correction"])
                )
                irp1step = _cap_step(step+q1, x, float(args.pandrosion_step_cap))
                if _distinct_step(irp1step, occupied):
                    i1Y = np.asarray([x+irp1step], np.complex128)
                    i1raw, i1bal = _raw_chart_values(
                        target, i1Y, reference, equation_weights
                    )
                    extra_steps.append(irp1step)
                    extra_labels.append("irp")
                    extra_raw_rows.append(i1raw[0])
                    extra_bal_rows.append(i1bal[0])
                    occupied.append(irp1step)
                    irp_proposals += 1

                    # Recompute and reclose once more: true IRP-2.
                    q2, _, q2diag = _reclosed_residual_layer(
                        Q1, pilot_Y[1], pilot_bal[1], i1Y[0], i1bal[0], q1, args
                    )
                    if q2 is not None:
                        transported_steps += 1
                        max_transport_correction = max(
                            max_transport_correction, float(q2diag["correction"])
                        )
                        irp2step = _cap_step(
                            step+q1+q2, x, float(args.pandrosion_step_cap)
                        )
                        if _distinct_step(irp2step, occupied):
                            i2Y = np.asarray([x+irp2step], np.complex128)
                            i2raw, i2bal = _raw_chart_values(
                                target, i2Y, reference, equation_weights
                            )
                            extra_steps.append(irp2step)
                            extra_labels.append("irp2")
                            extra_raw_rows.append(i2raw[0])
                            extra_bal_rows.append(i2bal[0])
                            occupied.append(irp2step)
                            irp2_proposals += 1

                append_dyadics(3)
                missing = len(extra_steps)-len(extra_raw_rows)
                if missing > 0:
                    tail_steps = extra_steps[len(extra_raw_rows):]
                    tail_Y = np.vstack([x+s for s in tail_steps])
                    tail_raw, tail_bal = _raw_chart_values(
                        target, tail_Y, reference, equation_weights
                    )
                    extra_raw_rows.extend(list(tail_raw))
                    extra_bal_rows.extend(list(tail_bal))

        if not sequential:
            # Robust 369-style branch: at most one acceleration and two raw
            # dyadic safeguards.  The IRP proposal now uses a reclosed chart.
            accelerator_candidates: list[tuple[float, np.ndarray, str]] = []
            if halley_step is not None:
                accelerator_candidates.append(
                    (halley_predicted, halley_step, "halley")
                )
            if (args.irp_layer
                    and linearity_defect <= float(args.linearity_irp_max)
                    and full_ratio <= float(args.linearity_irp_full_ratio)
                    and half_ratio < 1.0):
                # Outside a true sequential IRP/HIRP cell, preserve 369's
                # frozen residual layer.  Reclosure is reserved for branches
                # that actually reevaluate the intermediate residual.
                q, qmode = _factored_slope_step(
                    slope_factor, pilot_bal[1], x,
                    float(args.pandrosion_step_cap),
                )
                if qmode != "failed" and np.all(np.isfinite(q)):
                    q, _ = _palette_cap_residual_layer(q, step, x, args)
                    proposal = _cap_step(step+q, x, float(args.pandrosion_step_cap))
                    if _distinct_step(proposal, occupied):
                        score = max(0.0, full_ratio*linearity_defect)
                        accelerator_candidates.append((score, proposal, "irp"))
            if accelerator_candidates:
                accelerator_candidates.sort(key=lambda z: z[0])
                _, proposal, label = accelerator_candidates[0]
                extra_steps.append(proposal)
                extra_labels.append(label)
                occupied.append(proposal)
                if label == "halley":
                    halley_proposals += 1
                    max_halley_factor = max(
                        max_halley_factor, abs(float(halley_factor))
                    )
                else:
                    irp_proposals += 1
            append_dyadics(3)
            extra_Y = np.vstack([x+s for s in extra_steps[:3]])
            extra_raw, extra_bal = _raw_chart_values(
                target, extra_Y, reference, equation_weights
            )
            extra_raw_rows = list(extra_raw)
            extra_bal_rows = list(extra_bal)

        # Exactly three extra samples have now been consumed.
        extra_steps = extra_steps[:3]
        extra_labels = extra_labels[:3]
        extra_raw = np.asarray(extra_raw_rows[:3], np.complex128)
        extra_bal = np.asarray(extra_bal_rows[:3], np.complex128)
        line_evals += 5
        Y = np.vstack((pilot_Y, np.vstack([x+s for s in extra_steps])))
        raw_values = np.vstack((pilot_raw, extra_raw))
        balanced_values = np.vstack((pilot_bal, extra_bal))
        chart_residuals = norms(balanced_values)
        finite = np.isfinite(chart_residuals)
        if not np.any(finite):
            status = "invalid-chart-residual"
            break
        labels = ["raw-half", "raw-full"]+extra_labels
        i = int(np.argmin(np.where(finite, chart_residuals, np.inf)))
        previous_x = x.copy()
        x, fx = Y[i].copy(), raw_values[i].copy()
        previous_step = np.asarray(x-previous_x, np.complex128)
        if 0.0 < norm(previous_step) <= 1e-15:
            micro_steps_accepted += 1
        chosen = labels[i]
        if chosen == "halley":
            halley_accepted += 1
        elif chosen == "irp":
            irp_accepted += 1
        elif chosen == "irp2":
            irp2_accepted += 1
        elif chosen == "hirp":
            hirp_accepted += 1

        # Persistent equation chart: comparable root-equivalent backward error
        # across epochs, independent of the temporary common overflow chart.
        r = norm(_apply_equation_chart(fx, equation_weights))
        if r < best_r and norm(x) <= max_norm:
            best_x, best_f, best_r = x.copy(), fx.copy(), r
        done = ep+1

    ok = bool(best_r < args.accept and norm(best_x) <= max_norm)
    if ok:
        status = "converged"
    budget_used = int(target.oracle.eval_count)-call_start
    smin = smax = cond = None
    if args.diagnostic_svd and last_Q is not None:
        singular = np.linalg.svd(last_Q, compute_uv=False)
        if len(singular):
            smax, smin = float(singular[0]), float(singular[-1])
            cond = float(smax/smin) if smin > 0 else float("inf")
    return CorrectResult(
        y=best_x, f=best_f, residual=best_r, ok=ok, status=status,
        epochs=done, slopes=slopes, transported_steps=transported_steps,
        probe_evals=0, slope_evals=slope_evals, line_evals=line_evals,
        parabolic_evals=0, max_telescope_defect=max_defect,
        max_transport_correction=max_transport_correction,
        smin=smin, smax=smax, cond=cond, budget=budget,
        budget_used=budget_used,
        budget_exhausted=(status == "strict-budget-exhausted"),
        halley_proposals=halley_proposals,
        halley_accepted=halley_accepted,
        irp_proposals=irp_proposals,
        irp_accepted=irp_accepted,
        irp2_proposals=irp2_proposals,
        irp2_accepted=irp2_accepted,
        hirp_proposals=hirp_proposals,
        hirp_accepted=hirp_accepted,
        monomial_chain_proposals=monomial_chain_proposals,
        monomial_chain_accepted=monomial_chain_accepted,
        monomial_chain_degree_max=monomial_chain_degree_max,
        micro_steps_accepted=micro_steps_accepted,
        max_halley_factor=max_halley_factor,
        equation_chart_active=equation_chart_active,
        equation_chart_spread=equation_chart_spread,
        equation_chart_center=equation_chart_center,
    )


def add_results(base: CorrectResult, loc: CorrectResult, y: np.ndarray, f: np.ndarray,
                residual: float, status: str) -> CorrectResult:
    """Audit-only aggregation; it never launches another numerical method."""
    return CorrectResult(
        y=np.asarray(y, np.complex128), f=np.asarray(f, np.complex128),
        residual=residual, ok=bool(base.ok or loc.ok), status=status,
        epochs=base.epochs+loc.epochs, slopes=base.slopes+loc.slopes,
        transported_steps=base.transported_steps+loc.transported_steps,
        probe_evals=0, slope_evals=base.slope_evals+loc.slope_evals,
        line_evals=base.line_evals+loc.line_evals, parabolic_evals=0,
        max_telescope_defect=max(base.max_telescope_defect, loc.max_telescope_defect),
        max_transport_correction=max(base.max_transport_correction, loc.max_transport_correction),
        smin=loc.smin if loc.smin is not None else base.smin,
        smax=loc.smax if loc.smax is not None else base.smax,
        cond=loc.cond if loc.cond is not None else base.cond,
        budget=base.budget+loc.budget,
        budget_used=base.budget_used+loc.budget_used,
        budget_exhausted=base.budget_exhausted or loc.budget_exhausted,
        halley_proposals=base.halley_proposals+loc.halley_proposals,
        halley_accepted=base.halley_accepted+loc.halley_accepted,
        irp_proposals=base.irp_proposals+loc.irp_proposals,
        irp_accepted=base.irp_accepted+loc.irp_accepted,
        irp2_proposals=base.irp2_proposals+loc.irp2_proposals,
        irp2_accepted=base.irp2_accepted+loc.irp2_accepted,
        hirp_proposals=base.hirp_proposals+loc.hirp_proposals,
        hirp_accepted=base.hirp_accepted+loc.hirp_accepted,
        monomial_chain_proposals=(base.monomial_chain_proposals
                                  +loc.monomial_chain_proposals),
        monomial_chain_accepted=(base.monomial_chain_accepted
                                 +loc.monomial_chain_accepted),
        monomial_chain_degree_max=max(base.monomial_chain_degree_max,
                                      loc.monomial_chain_degree_max),
        micro_steps_accepted=base.micro_steps_accepted+loc.micro_steps_accepted,
        max_halley_factor=max(base.max_halley_factor, loc.max_halley_factor),
        equation_chart_active=base.equation_chart_active or loc.equation_chart_active,
        equation_chart_spread=max(base.equation_chart_spread, loc.equation_chart_spread),
        equation_chart_center=(loc.equation_chart_center if loc.equation_chart_active
                               else base.equation_chart_center),
    )


# ---------- starts, swarm, polish -------------------------------------------

def _portfolio_trial_index(trial: int) -> int:
    """Small-chart-first permutation preserving the original deterministic atlas.

    The successful small-radius directions of 328 are moved forward without
    changing their actual coordinates. Large charts remain available later.
    """
    head = [0, 11, 2, 13, 4, 15, 6, 17, 8, 19, 10, 1, 12, 3, 14, 5,
            16, 7, 18, 9, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
    if trial < len(head):
        return head[trial]
    return trial


def atlas_start(n: int, trial: int, seed: int) -> np.ndarray:
    original = _portfolio_trial_index(trial)
    radii = [2**(k/3) for k in range(-15,16)]
    return radii[(original*17)%len(radii)]*direction(n, original, seed)



def screened_start_portfolio(base: Target, args: argparse.Namespace, n: int, seed: int) -> tuple[list[np.ndarray], dict[str, Any]]:
    """Cheap automatic-start wrapper for 330.

    SciPy ``root`` itself does not choose x0; it accepts the user's x0.  For a
    fair standalone harvester, 330 builds a deterministic origin-first bank,
    evaluates it in two vectorized radial charts, and ranks candidates by both
    residual and observed radial decrease.  This costs two batched oracle calls
    and no derivatives/Jacobian.
    """
    wanted=max(args.pool, args.start_keep)
    raw=max(wanted, args.start_candidates)
    radii=floats(args.start_radii, [.125,.25,.5,1.,2.])
    C=[np.zeros(n,np.complex128)]
    # low-discrepancy deterministic directions; interleave radii to avoid a
    # large-radius block dominating the top of the list.
    for i in range(raw-1):
        r=radii[(i*7)%len(radii)]
        C.append(r*direction(n,i,seed+0x330330))
    Y=np.stack(C)
    F0=base.eval_batch(Y); R0=norms(F0)
    shrink=float(args.start_shrink)
    Ys=shrink*Y
    Fs=base.eval_batch(Ys); Rs=norms(Fs)
    # prefer low residuals, but reward a clear inward contraction; logarithms
    # make the score scale robust across KS realizations.
    eps=1e-300
    contraction=np.log((R0+eps)/(Rs+eps))
    score=np.log(R0+eps)-args.start_contraction_weight*np.maximum(contraction,0.0)
    # Origin is always tried first, matching the common neutral x0 used in
    # local-solver benchmarks; remaining starts are residual-ranked/diverse.
    selected=[0]
    sep=args.start_sep or .08*math.sqrt(n)
    for i in np.argsort(score):
        i=int(i)
        if i==0: continue
        # choose whichever of y and shrink*y empirically has lower residual
        cand=Ys[i] if Rs[i] < R0[i] else Y[i]
        if all(norm(cand-(Ys[j] if Rs[j] < R0[j] else Y[j]))>sep for j in selected[1:]):
            selected.append(i)
        if len(selected)>=wanted: break
    out=[]
    for i in selected:
        if i==0: out.append(Y[0].copy())
        else: out.append((Ys[i] if Rs[i] < R0[i] else Y[i]).copy())
    return out,{"enabled":True,"mode":"origin-first-two-chart-screening","candidates":len(Y),
        "kept":len(out),"batch_evals":2*len(Y),"best_initial_residual":float(min(np.min(R0),np.min(Rs))),
        "origin_residual":float(R0[0]),"shrink":shrink}


def pilot_race_portfolio(base: Target, starts: list[np.ndarray], args: argparse.Namespace) -> tuple[list[np.ndarray], list[Optional[CachedStart]], dict[str, Any]]:
    """Race candidates and transfer the winning endpoint value without reevaluation.

    The local model is deliberately rebuilt by the same telescopic Pandrosion rule;
    only the already-known endpoint and residual are cached.  This preserves 335's
    robust dynamics while removing a redundant oracle call at every continuation.
    """
    if not starts or args.pilot_candidates <= 0 or args.pilot_epochs <= 0:
        return starts,[None for _ in starts],{"enabled":False}
    k=min(len(starts),max(1,args.pilot_candidates)); raced=[]; t0=clock()
    for i in range(k):
        y0=np.asarray(starts[i],np.complex128); r0=norm(base.eval(y0)); before=base.oracle.eval_count
        loc=correct(base,y0,args,[],args.pilot_epochs,None)
        used=max(1,base.oracle.eval_count-before); gain=math.log(max(r0,1e-300)/max(loc.residual,1e-300)); score=gain/used
        cache=CachedStart(loc.y.copy(),loc.f.copy())
        raced.append((0 if loc.ok else 1,-score,loc.residual,i,loc.y.copy(),cache,used,r0,gain))
    raced.sort(key=lambda x:(x[0],x[1],x[2]))
    order=[r[4] for r in raced]; caches=[r[5] for r in raced]
    order.extend(np.asarray(x,np.complex128).copy() for x in starts[k:]); caches.extend(None for _ in starts[k:])
    return order,caches,{"enabled":True,"candidates":k,"epochs":args.pilot_epochs,"seconds":clock()-t0,
        "oracle_samples":sum(r[6] for r in raced),"winner_original_index":int(raced[0][3]),
        "winner_initial_residual":float(raced[0][7]),"winner_pilot_residual":float(raced[0][2]),
        "winner_log_gain":float(raced[0][8]),"cached_endpoint_handoff":True,
        "mode":"356-single-origin-three-step-cached-base-slope"}


def swarm(base: Target, args: argparse.Namespace, n: int, seed: int) -> tuple[list[np.ndarray], dict[str, Any]]:
    """Oversampled, diversity-preserving Pandrosion atlas.

    329 first evaluates a radial/phase portfolio in one batch, retains a
    separated beam across residual quantiles, and only then spends finite-slope
    layers on the beam.  This prevents the swarm from collapsing into one basin
    before the expensive corrector starts.
    """
    size = args.swarm_size or min(args.pool, max(32, 8*args.count))
    if isinstance(base.oracle, ExactKSGPOracle): size = min(size, args.gp_swarm_cap)
    keep = min(size, args.swarm_keep or max(8, 3*args.count))
    oversample = 1 if isinstance(base.oracle, ExactKSGPOracle) else 4
    total = max(size, oversample*size)
    Y0 = np.stack([atlas_start(n, i, seed+31337) for i in range(total)])
    R0 = norms(base.eval_batch(Y0))
    sep = args.swarm_sep or .12*math.sqrt(n)

    # Residual-ranked but spatially separated pre-beam.
    pre=[]
    for i in np.argsort(R0):
        if not math.isfinite(float(R0[i])): continue
        if all(norm(Y0[i]-Y0[j]) > sep for j in pre): pre.append(int(i))
        if len(pre)>=size: break
    if len(pre)<size:
        for i in np.argsort(R0):
            if int(i) not in pre: pre.append(int(i))
            if len(pre)>=size: break
    Y=Y0[pre].copy(); R=R0[pre].copy()

    slope_layers=transports=0
    deadline=clock()+args.swarm_timeout if args.swarm_timeout>0 else None
    # Correct only the best half plus evenly spaced explorers.
    order=list(np.argsort(R))
    budget=max(keep, min(size, 2*keep))
    chosen=order[:budget]
    for rank,i in enumerate(chosen):
        if args.swarm_epochs<=0 or (deadline is not None and clock()>=deadline): break
        loc=correct(base,Y[int(i)],args,[],args.swarm_epochs,deadline)
        Y[int(i)],R[int(i)]=loc.y.copy(),loc.residual
        slope_layers+=loc.slopes; transports+=loc.transported_steps

    selected=[]
    # Strong diversity after local preconditioning.
    for i in np.argsort(R):
        if not math.isfinite(float(R[i])): continue
        local_sep=sep*(.65 if len(selected)>=keep//2 else 1.0)
        if all(norm(Y[i]-Y[j]) > local_sep for j in selected): selected.append(int(i))
        if len(selected)>=keep: break
    return [Y[i].copy() for i in selected], {"size":size,"oversampled":total,"kept":len(selected),
        "best_residual":float(R[selected[0]]) if selected else None,"evals":len(Y0),
        "finite_slope_layers":slope_layers,"transported_steps":transports,
        "mode":"327-adaptive-competitive-full-pandrosion-atlas"}


def conditioning(loc: CorrectResult) -> dict[str, Any]:
    if loc.smin is None or loc.smax is None:
        return {"smin":None,"smax":None,"cond":None,"near_multiple":None,"singular":None,"mode":"unavailable"}
    return {"smin":loc.smin,"smax":loc.smax,"cond":loc.cond,
        "near_multiple":bool(loc.smin<=1e-8*max(loc.smax,1e-300)),
        "singular":bool(loc.smin<=np.finfo(float).eps*max(loc.smax,1e-300)),
        "mode":"finite-slope-endpoint"}


# ---------- case driver ------------------------------------------------------

def run_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    started=clock(); n,d=parse_case(case_raw); oracle=make_oracle(args,n,d); base=Target(oracle)
    origin_error=oracle.backward_error(np.zeros(n,np.complex128)); explicit=starts_from_text(args.starts,n)
    swarm_points,swarm_meta,cap_error=[],{"enabled":False},None
    start_points,start_meta=[],{"enabled":False}
    if not explicit and args.count:
        start_points=[np.zeros(n,np.complex128)]
        start_meta={"enabled":True,"mode":"single-origin-no-screen","candidates":1,"evals":0}
    pilot_meta={"enabled":False}; start_caches=[None for _ in start_points]
    if not explicit and start_points and args.count:
        start_points,start_caches,pilot_meta=pilot_race_portfolio(base,start_points,args)
    if args.swarm and not explicit and args.count:
        try: swarm_points,swarm_meta=swarm(base,args,n,oracle.seed); swarm_meta["enabled"]=True
        except RuntimeError as exc: cap_error=str(exc); swarm_meta={"enabled":True,"status":"oracle-cap","error":cap_error}
    legacy_points=[]
    priority=explicit+start_points+legacy_points+swarm_points
    priority_caches=[None for _ in explicit]+start_caches+[None for _ in legacy_points]+[None for _ in swarm_points]
    roots=[]; trials=[]; duplicates=failures=0
    for trial in range(args.pool):
        if len(roots)>=args.count: break
        cached0=None
        if trial < len(priority):
            y0 = priority[trial].copy(); cached0=priority_caches[trial]
        else:
            y0 = atlas_start(n, trial-len(priority), oracle.seed)
        known=[np.asarray(r["z_complex"],np.complex128) for r in roots]
        if known and min(norm(y0-r) for r in known)<=args.early_dup_sep:
            duplicates+=1; trials.append({"trial":trial,"status":"start-near-root","accepted":False}); continue
        deadline=clock()+args.trial_timeout if args.trial_timeout>0 else None; before=oracle.eval_count
        try:
            loc=correct(base,y0,args,known,args.epochs,deadline,cached0)
            z,f,residual=loc.y.copy(),loc.f.copy(),loc.residual; condmeta=conditioning(loc)
            validation_error=oracle.backward_error(z)
        except RuntimeError as exc:
            cap_error=str(exc); failures+=1; trials.append({"trial":trial,"status":"oracle-cap","accepted":False,
                "error":cap_error,"oracle_samples":oracle.eval_count-before}); break
        accepted=bool(math.isfinite(residual) and residual<args.accept and validation_error<args.validation_accept)
        rec={"trial":trial,"status":loc.status,"accepted":accepted,"residual":residual,
            "validation_error":validation_error,"epochs":loc.epochs,"finite_slopes":loc.slopes,
            "transported_finite_slope_steps":loc.transported_steps,
            "probe_evals":loc.probe_evals,"slope_path_evals":loc.slope_evals,"line_evals":loc.line_evals,
            "parabolic_evals":loc.parabolic_evals,"max_telescope_defect":loc.max_telescope_defect,
            "max_transport_correction":loc.max_transport_correction,
            "strict_budget":loc.budget,"strict_budget_used":loc.budget_used,
            "strict_budget_exhausted":loc.budget_exhausted,
            "halley_proposals":loc.halley_proposals,"halley_accepted":loc.halley_accepted,
            "irp_proposals":loc.irp_proposals,"irp_accepted":loc.irp_accepted,
            "irp2_proposals":loc.irp2_proposals,"irp2_accepted":loc.irp2_accepted,
            "hirp_proposals":loc.hirp_proposals,"hirp_accepted":loc.hirp_accepted,
            "monomial_chain_proposals":loc.monomial_chain_proposals,
            "monomial_chain_accepted":loc.monomial_chain_accepted,
            "monomial_chain_degree_max":loc.monomial_chain_degree_max,
            "micro_steps_accepted":loc.micro_steps_accepted,
            "max_halley_factor":loc.max_halley_factor,
            "equation_chart_active":loc.equation_chart_active,
            "equation_chart_spread":loc.equation_chart_spread,
            "equation_chart_center":loc.equation_chart_center,
            "oracle_samples":oracle.eval_count-before,"conditioning":condmeta}
        if args.verbose_trials: rec.update({"start":y0,"z":z})
        if not accepted:
            if residual<args.accept: rec["status"]="validation-failed"
            failures+=1; trials.append(rec); continue
        dup=next((i for i,r in enumerate(roots) if norm(z-r["z_complex"])<=args.cluster_sep),None)
        if dup is not None:
            duplicates+=1; rec.update({"status":"duplicate","cluster":dup}); trials.append(rec); continue
        root={"id":len(roots),"trial":trial,"z_complex":z.copy(),"residual":residual,
            "validation_error":validation_error,"realness":norm(z.imag)/max(norm(z),1e-300),**condmeta}
        roots.append(root); rec.update({"status":"new-root","root_id":root["id"]}); trials.append(rec)
    encoded=[]
    for root in roots:
        rr=dict(root); rr["z"]=[[float(x.real),float(x.imag)] for x in rr.pop("z_complex")]; encoded.append(rr)
    gp=isinstance(oracle,ExactKSGPOracle); feature=isinstance(oracle,KostlanOracle) and oracle.feature_mode
    return {"script":Path(__file__).name,"version":371,"standalone":True,"case":f"{n},{d}","n":n,"degree":d,
      "method":{"family":"pandrosion-371-monomial-chain-hirp","derivatives":False,"jacobian":False,"newton":False,
        "halley":True,"halley_kind":"projected derivative-free residual homothety from t=0,1/2,1",
        "monomial_irp":False,"irp_analogue":True,"true_residual_recompute":True,
        "irp_kind":"reclosed finite-slope IRP-1/IRP-2 with actual residual recomputation; HIRP-1 is Halley followed by one reclosed raw layer",
        "singular_chart":"degree-inferred upper/lower triangular monomial-chain closure with full residual validation",
        "residual_palette":"dyadic additive-layer chart scaling",
        "quasi_newton":False,"homotopy":False,"fallback":False,"restarts":False,
        "levenberg_marquardt":False,"scipy":False,
        "identity":"F(b)-F(a)=Q(a,b)(b-a) by coordinate telescoping and exact one-column chart reclosure",
        "chart":"common overflow chart plus optional frozen diagonal equation homothety",
        "globalization":"two raw pilots plus either HIRP-1, IRP-2, or a one-accelerator/two-dyadic portfolio; five samples total"},
      "oracle":{"kind":oracle.kind,"seed":oracle.seed,"samples":oracle.eval_count,"seconds":oracle.seconds,
        "residual_mode":"raw" if isinstance(oracle,ExpressionOracle) else "root-equivalent-degree-normalized",
        "validation_mode":"raw-residual" if isinstance(oracle,ExpressionOracle) else "normalized-gp-residual" if gp else "scale-invariant-backward-error",
        "origin_backward_error":origin_error,"exact_terms_per_polynomial":math.comb(n+d,d),
        "feature_count":len(oracle.exps) if isinstance(oracle,KostlanOracle) else None,
        "surrogate_warning":"root is validated only for this fixed feature bank" if feature else None,
        "gp_unique_points":len(oracle.cache) if gp else None,"gp_active_covariance_rank":len(oracle.points) if gp else None,
        "gp_stabilizations":oracle.stabilizations if gp else None,
        "gp_max_covariance_correction":oracle.max_correction if gp else None,
        "gp_rank_tolerance":oracle.jitter if gp else None,
        "gp_exact_finite_dimensional_law":"up to floating-point arithmetic and reported rank tolerance" if gp else None,
        "gp_reproducibility_scope":"exact in law for this adaptive query sequence; query order selects the realization" if gp else None},
      "parameters":{"count":args.count,"pool":args.pool,"epochs":args.epochs,"accept":args.accept,
        "ks_backend":args.ks_backend,"pandrosion_radius":args.pandrosion_radius,
        "pandrosion_radius_mode":args.pandrosion_radius_mode,
        "pandrosion_radius_step_ratio":args.pandrosion_radius_step_ratio,
        "pandrosion_step_cap":args.pandrosion_step_cap,"line_grid":args.line_grid,
        "directional_halley":args.directional_halley,"irp_layer":args.irp_layer,
        "irp2_layer":args.irp2_layer,"hirp_layer":args.hirp_layer,
        "monomial_chain_chart":args.monomial_chain_chart,
        "equation_chart":args.equation_chart,
        "residual_palette_base":args.residual_palette_base,
        "halley_factor_bounds":[args.halley_min_factor,args.halley_max_factor],
        "strict_budget":f"{args.strict_budget_factor}*n+{args.strict_budget_base}",
        "swarm":args.swarm},
      "portfolio":{"start_meta":start_meta,"pilot_meta":pilot_meta},"swarm":swarm_meta,"roots":encoded,"trials":trials if args.verbose_trials else trials[:args.keep_trials],
      "summary":{"requested":args.count,"unique":len(roots),"success":len(roots)>=args.count,
        "trials":len(trials),"duplicates":duplicates,"failures":failures,
        "oracle_samples":oracle.eval_count,"seconds":clock()-started,"oracle_cap_error":cap_error}}


# ---------- CLI and regressions ---------------------------------------------

def parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser(description="Pandrosion 371 monomial-chain chart plus residual IRP-2/HIRP-1 standalone engine")
    p.add_argument("--cases",default="2,4"); p.add_argument("--seed-index",type=int,default=0)
    p.add_argument("--system-source",choices=["ks","kostlan","poly","polynomial"],default="ks")
    p.add_argument("--ks-backend",choices=["auto","dense","gp","feature","gaussian"],default="auto")
    p.add_argument("--polys"); p.add_argument("--variables"); p.add_argument("--starts")
    p.add_argument("--dense-max-terms",type=int,default=250000); p.add_argument("--features",type=int,default=3072)
    p.add_argument("--eval-block",type=int,default=128); p.add_argument("--equation-normalize",action="store_true",default=False)
    p.add_argument("--gp-jitter",type=float,default=1e-9); p.add_argument("--gp-max-points",type=int,default=10000)
    p.add_argument("--gp-jet-radius",type=float,default=1e-3); p.add_argument("--gp-swarm-cap",type=int,default=8)
    p.add_argument("--gp-probe-cap",type=int,default=4)
    p.add_argument("--count",type=int,default=1); p.add_argument("--pool",type=int,default=1)
    p.add_argument("--epochs",type=int,default=96); p.add_argument("--accept",type=float,default=1e-8)
    p.add_argument("--internal-accept-factor",type=float,default=.5,
                   help="termination safety factor for oracle equation units")
    p.add_argument("--screen-starts",action="store_true",default=True)
    p.add_argument("--no-screen-starts",dest="screen_starts",action="store_false")
    p.add_argument("--start-candidates",type=int,default=24)
    p.add_argument("--start-keep",type=int,default=8)
    p.add_argument("--start-radii",default=".125,.25,.5,1,2")
    p.add_argument("--start-shrink",type=float,default=.5)
    p.add_argument("--start-contraction-weight",type=float,default=.75)
    p.add_argument("--start-sep",type=float,default=0)
    p.add_argument("--pilot-candidates",type=int,default=0)
    p.add_argument("--pilot-epochs",type=int,default=3)
    p.add_argument("--legacy-head",type=int,default=2)
    p.add_argument("--tol",type=float,default=1e-12); p.add_argument("--cluster-sep",type=float,default=1e-7)
    p.add_argument("--early-dup-sep",type=float,default=1e-4); p.add_argument("--trial-timeout",type=float,default=0)
    p.add_argument("--probes",type=int,default=1); p.add_argument("--probe-top",type=int,default=1)
    p.add_argument("--probe-scale",type=float,default=.03); p.add_argument("--probe-radii",default=".5,1,2")
    p.add_argument("--probe-cond-weight",type=float,default=.02); p.add_argument("--equal-value-weight",type=float,default=.02)
    p.add_argument("--trust-radius",type=float,default=10); p.add_argument("--line-grid",default="1,.5,.25,.125,.0625")
    p.add_argument("--transport-steps",type=int,default=8); p.add_argument("--transport-ratio",type=float,default=.85)
    p.add_argument("--transport-trust-radius",type=float,default=5); p.add_argument("--transport-line-grid",default="1,.5,.25")
    p.add_argument("--transport-closure-limit",type=float,default=2.)
    p.add_argument("--slope-retries",type=int,default=2); p.add_argument("--slope-retry-factor",type=float,default=.25)
    p.add_argument("--deflation-alpha",type=float,default=.15)
    p.add_argument("--parabolic",action="store_true",default=True); p.add_argument("--no-parabolic",dest="parabolic",action="store_false")
    p.add_argument("--swarm",action="store_true",default=False); p.add_argument("--no-swarm",dest="swarm",action="store_false")
    p.add_argument("--swarm-size",type=int,default=48); p.add_argument("--swarm-keep",type=int,default=12)
    p.add_argument("--swarm-sep",type=float,default=0); p.add_argument("--swarm-epochs",type=int,default=1)
    p.add_argument("--swarm-timeout",type=float,default=0)
    p.add_argument("--validation-accept",type=float,default=1e-8); p.add_argument("--keep-trials",type=int,default=100)
    p.add_argument("--pandrosion-radius",type=float,default=1e-5)
    p.add_argument("--pandrosion-radius-mode",choices=["fixed","adaptive"],default="adaptive")
    p.add_argument("--pandrosion-radius-step-ratio",type=float,default=.25)
    p.add_argument("--pandrosion-step-cap",type=float,default=10.0)
    p.add_argument("--degree-radius-chart",action="store_true",default=True)
    p.add_argument("--no-degree-radius-chart",dest="degree_radius_chart",action="store_false")
    p.add_argument("--degree-radius-cap",type=float,default=.05)
    p.add_argument("--directional-halley",action="store_true",default=True)
    p.add_argument("--no-directional-halley",dest="directional_halley",action="store_false")
    p.add_argument("--irp-layer",action="store_true",default=True)
    p.add_argument("--no-irp-layer",dest="irp_layer",action="store_false")
    p.add_argument("--irp2-layer",action="store_true",default=True)
    p.add_argument("--no-irp2-layer",dest="irp2_layer",action="store_false")
    p.add_argument("--hirp-layer",action="store_true",default=True)
    p.add_argument("--no-hirp-layer",dest="hirp_layer",action="store_false")
    p.add_argument("--halley-min-factor",type=float,default=.0625)
    p.add_argument("--halley-max-factor",type=float,default=2.0)
    p.add_argument("--halley-prediction-ratio",type=float,default=.92)
    p.add_argument("--linearity-irp-max",type=float,default=.12)
    p.add_argument("--linearity-irp-full-ratio",type=float,default=.92)
    p.add_argument("--linearity-halley-min",type=float,default=1e-4)
    p.add_argument("--linearity-halley-max",type=float,default=.45)
    p.add_argument("--halley-min-projected-slope",type=float,default=.05)
    p.add_argument("--hirp-linearity-max",type=float,default=.06)
    p.add_argument("--hirp-full-ratio",type=float,default=.80)
    p.add_argument("--hirp-half-ratio",type=float,default=.95)
    p.add_argument("--hirp-factor-window",type=float,default=.075)
    p.add_argument("--hirp-residual-trigger",type=float,default=.003)
    p.add_argument("--irp2-linearity-max",type=float,default=.035)
    p.add_argument("--irp2-full-ratio",type=float,default=.82)
    p.add_argument("--residual-transport-limit",type=float,default=1.0)
    p.add_argument("--residual-layer-ratio-cap",type=float,default=1.0)
    p.add_argument("--residual-palette-base",type=float,default=2.0)
    p.add_argument("--residual-palette-min-factor",type=float,default=.0625)
    p.add_argument("--equation-chart",action="store_true",default=True)
    p.add_argument("--no-equation-chart",dest="equation_chart",action="store_false")
    p.add_argument("--equation-chart-spread",type=float,default=6.0)
    p.add_argument("--equation-chart-center",type=float,default=3.5)
    p.add_argument("--monomial-chain-chart",action="store_true",default=True)
    p.add_argument("--no-monomial-chain-chart",dest="monomial_chain_chart",action="store_false")
    p.add_argument("--monomial-chain-telescope-radius",type=float,default=.25)
    p.add_argument("--monomial-chain-probe-radius",type=float,default=.25)
    p.add_argument("--monomial-chain-residual-min",type=float,default=1e-4)
    p.add_argument("--monomial-chain-concentration",type=float,default=.85)
    p.add_argument("--monomial-chain-flat-fraction",type=float,default=.80)
    p.add_argument("--monomial-chain-flat-tolerance",type=float,default=1e-10)
    p.add_argument("--monomial-chain-band-tolerance",type=float,default=1e-8)
    p.add_argument("--monomial-chain-terminal-fraction",type=float,default=.80)
    p.add_argument("--monomial-chain-degree-max",type=int,default=32)
    p.add_argument("--monomial-chain-degree-mismatch",type=float,default=1e-5)
    p.add_argument("--monomial-chain-probe-signal",type=float,default=1e-10)
    p.add_argument("--monomial-chain-coefficient-mismatch",type=float,default=1e-3)
    p.add_argument("--monomial-chain-min-improvement",type=float,default=1e-10)
    p.add_argument("--strict-budget-factor",type=int,default=40)
    p.add_argument("--strict-budget-base",type=int,default=100)
    p.add_argument("--max-solution-norm",type=float,default=1e8)
    p.add_argument("--diagnostic-svd",action="store_true",default=False,help="compute singular-value diagnostics on every rebuilt slope")
    p.add_argument("--verbose-trials",action="store_true"); p.add_argument("--self-test",action="store_true")
    p.add_argument("--out"); p.add_argument("--outdir",default="/tmp/371_pandrosion")
    return p


def validate(args: argparse.Namespace) -> None:
    for raw in str(args.cases).replace("|",";").split(";"):
        if raw.strip(): parse_case(raw)
    if args.count<0 or args.pool<0 or args.epochs<=0 or args.probes<=0 or args.probe_top<=0: raise ValueError("invalid budgets")
    if args.pilot_candidates<0 or args.pilot_epochs<0: raise ValueError("invalid pilot controls")
    if args.legacy_head<0 or args.start_candidates<=0 or args.start_keep<=0 or not 0<args.start_shrink<1 or args.start_contraction_weight<0: raise ValueError("invalid start-screen controls")
    if args.swarm_epochs<0 or args.swarm_timeout<0: raise ValueError("invalid swarm controls")
    if args.accept<=0 or args.validation_accept<=0 or args.trust_radius<=0: raise ValueError("invalid tolerances")
    if not (0 < args.internal_accept_factor <= 1):
        raise ValueError("internal_accept_factor must lie in (0,1]")
    if args.gp_jitter<=0 or args.gp_max_points<=0 or args.gp_jet_radius<=0 or args.gp_swarm_cap<=0 or args.gp_probe_cap<=0: raise ValueError("invalid GP controls")
    if args.probe_scale<=0 or args.probe_cond_weight<0 or args.equal_value_weight<0: raise ValueError("invalid probe controls")
    if args.transport_steps<0 or not 0<args.transport_ratio<=1 or args.transport_trust_radius<=0 or args.transport_closure_limit<=0: raise ValueError("invalid finite-slope transport controls")
    if args.slope_retries<0 or not 0<args.slope_retry_factor<1: raise ValueError("invalid finite-slope retry controls")
    if args.pandrosion_radius<=0 or args.pandrosion_step_cap<=0 or args.max_solution_norm<=0 or args.degree_radius_cap<=0: raise ValueError("invalid pure Pandrosion controls")
    if not 0<args.pandrosion_radius_step_ratio<=1: raise ValueError("invalid adaptive-radius ratio")
    if not 0<args.halley_min_factor<=args.halley_max_factor or args.halley_min_projected_slope<0: raise ValueError("invalid Halley controls")
    if not (0<=args.hirp_linearity_max and 0<args.hirp_full_ratio<=1 and 0<args.hirp_half_ratio<=1 and args.hirp_factor_window>0 and args.hirp_residual_trigger>0): raise ValueError("invalid HIRP controls")
    if not (0<=args.irp2_linearity_max and 0<args.irp2_full_ratio<=1): raise ValueError("invalid IRP-2 controls")
    if args.residual_transport_limit<=0 or args.residual_layer_ratio_cap<=0: raise ValueError("invalid residual layer controls")
    if args.residual_palette_base<=1 or not 0<args.residual_palette_min_factor<=1: raise ValueError("invalid residual palette")
    if args.equation_chart_spread<0 or args.equation_chart_center<0: raise ValueError("invalid equation chart")
    if (args.monomial_chain_telescope_radius<=0 or args.monomial_chain_probe_radius<=0
            or args.monomial_chain_residual_min<0 or args.monomial_chain_flat_tolerance<0
            or args.monomial_chain_degree_max<2):
        raise ValueError("invalid monomial-chain chart")
    if not (0<args.monomial_chain_concentration<=1
            and 0<args.monomial_chain_flat_fraction<=1
            and args.monomial_chain_band_tolerance>=0
            and 0<args.monomial_chain_terminal_fraction<=1
            and args.monomial_chain_degree_mismatch>=0
            and args.monomial_chain_probe_signal>=0
            and args.monomial_chain_coefficient_mismatch>=0
            and 0<=args.monomial_chain_min_improvement<1):
        raise ValueError("invalid monomial-chain thresholds")
    if args.strict_budget_factor<=0 or args.strict_budget_base<0: raise ValueError("invalid strict budget")


def self_test(args: argparse.Namespace) -> dict[str,Any]:
    checks=[]

    # Exact black-box telescope.
    poly=ExpressionOracle(2,2,7,"x1^2+x1*x2-1;x2^2-x1",None)
    target=Target(poly)
    aa=np.asarray([.2+.1j,-.3j]); bb=np.asarray([1-.2j,.4+.1j])
    Q,defect,_=finite_slope(target,aa,bb)
    identity=norm((target.eval(bb)-target.eval(aa))-Q@(bb-aa))
    checks.append({"name":"exact-coordinate-telescope","passed":identity<1e-12 and defect<1e-12,
        "identity_error":identity,"defect":defect})

    # Two distinct roots, each reached by the same corrector.
    a=argparse.Namespace(**vars(args)); a.system_source="poly"
    a.polys="x^2-3*x-10"; a.variables=None; a.starts="-8,4"
    a.count=2; a.pool=2; a.swarm=False; a.pilot_candidates=0; a.epochs=96
    roots=run_case(a,"1,2")
    checks.append({"name":"two-roots-irp-halley-pandrosion","passed":roots["summary"]["success"],
        "unique":roots["summary"]["unique"]})

    # The scalar cube-root closure uses the adaptive Halley/IRP portfolio and
    # reaches high accuracy under the same hard sample cap.
    cube_args=argparse.Namespace(**vars(args)); cube_args.accept=1e-12
    cube=ExpressionOracle(1,3,stable_seed(1,3,0),"x^3-2",None)
    cube_loc=correct(Target(cube),np.asarray([1+0j]),cube_args,[],96,None)
    checks.append({"name":"cube-root-hirp-residual-composition","passed":cube_loc.ok
        and cube.eval_count<=25 and cube_loc.halley_proposals>0
        and cube_loc.hirp_proposals>0 and cube_loc.hirp_accepted>0,
        "residual":cube_loc.residual,"samples":cube.eval_count,
        "epochs":cube_loc.epochs,"halley_accepted":cube_loc.halley_accepted,
        "hirp_accepted":cube_loc.hirp_accepted,
        "transported_layers":cube_loc.transported_steps})

    # Disable Halley to audit true residual recomputation across two raw layers.
    irp2_args=argparse.Namespace(**vars(args)); irp2_args.accept=1e-12
    irp2_args.directional_halley=False; irp2_args.hirp_layer=False
    irp2=ExpressionOracle(1,3,stable_seed(1,3,1),"x^3-2",None)
    irp2_loc=correct(Target(irp2),np.asarray([1+0j]),irp2_args,[],96,None)
    checks.append({"name":"two-layer-irp-recomputes-residual","passed":irp2_loc.ok
        and irp2_loc.irp2_proposals>0 and irp2_loc.irp2_accepted>0
        and irp2_loc.transported_steps>=2,
        "residual":irp2_loc.residual,"samples":irp2.eval_count,
        "irp2_accepted":irp2_loc.irp2_accepted,
        "transported_layers":irp2_loc.transported_steps})

    # A globally tiny equation multiplier must not create a false initial root.
    class TinyEquationOracle(Oracle):
        def __init__(self) -> None:
            super().__init__(1,3,stable_seed(1,3,2),"tiny-equation-scale")
        def _eval_batch(self,Z:np.ndarray)->np.ndarray:
            return 1e-12*(Z[:,[0]]**3-2.0)
    tiny=TinyEquationOracle()
    tiny_loc=correct(Target(tiny),np.asarray([1+0j]),cube_args,[],96,None)
    tiny_unscaled=abs(complex(tiny_loc.y[0])**3-2.0)
    checks.append({"name":"frozen-equation-homothety","passed":tiny_loc.ok
        and tiny_loc.equation_chart_active and tiny_unscaled<1e-10,
        "internal_residual":tiny_loc.residual,"unscaled_residual":tiny_unscaled,
        "samples":tiny.eval_count,"chart_center":tiny_loc.equation_chart_center})

    # Projective dense normalization is frozen, preventing roots at infinity.
    dense_seed=stable_seed(1,3,4)
    dense=KostlanOracle(1,3,dense_seed,math.comb(4,3),512,False,128)
    dense_start=atlas_start(1,0,dense_seed)
    dense_loc=correct(Target(dense),dense_start,args,[],96,None)
    checks.append({"name":"dense-frozen-chart","passed":dense_loc.ok and norm(dense_loc.y)<1e8
        and dense.eval_count<=40*1+100,"residual":dense_loc.residual,
        "solution_norm":norm(dense_loc.y),"samples":dense.eval_count})

    # Fixed-feature normalization uses the same root-equivalent chart rule.
    n,d=10,20; feature_seed=stable_seed(n,d,0,0x3595C1)
    feature=KostlanOracle(n,d,feature_seed,0,512,False,128)
    feature_loc=correct(Target(feature),np.zeros(n,np.complex128),args,[],96,None)
    checks.append({"name":"feature-frozen-chart","passed":feature_loc.ok and norm(feature_loc.y)<1e8
        and feature.eval_count<=40*n+100,"residual":feature_loc.residual,
        "solution_norm":norm(feature_loc.y),"samples":feature.eval_count})

    # Stationary coupled chain: no escape and no alternate solver.
    n,d=10,2; rng=np.random.default_rng(stable_seed(n,d,0,0x360C0A9))
    phases=rng.uniform(-math.pi,math.pi,n); magnitudes=rng.uniform(.7,1.3,n)
    coefficients=magnitudes*np.exp(1j*phases)
    literal=lambda z:f"({z.real:.17g}{z.imag:+.17g}j)"
    equations=[f"x{i+1}^2-{literal(complex(coefficients[i]))}*x{i+2}" for i in range(n-1)]
    equations.append(f"x{n}^2-{literal(complex(coefficients[-1]))}")
    coupled=ExpressionOracle(n,d,stable_seed(n,d,0),";".join(equations),None)
    coupled_loc=correct(Target(coupled),np.zeros(n,np.complex128),args,[],96,None)
    checks.append({"name":"stationary-monomial-chain-chart","passed":coupled_loc.ok
        and coupled.eval_count<=40*n+100
        and coupled_loc.monomial_chain_accepted>0
        and coupled_loc.monomial_chain_degree_max==2,
        "residual":coupled_loc.residual,"samples":coupled.eval_count,
        "chain_accepted":coupled_loc.monomial_chain_accepted,
        "inferred_degree":coupled_loc.monomial_chain_degree_max})

    # Reciprocal/small-root chart: corrections below the former absolute
    # 1e-15 step cutoff must remain representable relative to a 1e-8 root.
    class TinyRootOracle(Oracle):
        def __init__(self) -> None:
            super().__init__(1,2,stable_seed(1,2,371),"tiny-root-chart")
        def _eval_batch(self,Z:np.ndarray)->np.ndarray:
            y=Z[:,[0]]/1e-8
            return y+0.08*y*y-1.08
    tiny_root=TinyRootOracle()
    tiny_root_loc=correct(Target(tiny_root),np.asarray([6.99e-9+0j]),args,[],96,None)
    checks.append({"name":"small-root-scale-aware-corrections","passed":tiny_root_loc.ok
        and abs(complex(tiny_root_loc.y[0])-1e-8)<1e-18,
        "residual":tiny_root_loc.residual,"samples":tiny_root.eval_count,
        "root":complex(tiny_root_loc.y[0])})

    # Singular double root under the same finite-slope dynamics.
    multiple=ExpressionOracle(1,2,19,"(x-(.7+.4j))^2",None)
    multiple_loc=correct(Target(multiple),np.zeros(1,np.complex128),args,[],96,None)
    checks.append({"name":"multiple-root-irp-halley","passed":multiple_loc.ok,
        "residual":multiple_loc.residual,"samples":multiple.eval_count})

    # An impossible system must stop without crossing the cap.
    impossible=ExpressionOracle(1,1,23,"1",None)
    impossible_loc=correct(Target(impossible),np.zeros(1,np.complex128),args,[],1000,None)
    checks.append({"name":"hard-cap","passed":not impossible_loc.ok
        and impossible_loc.budget_used<=impossible_loc.budget,
        "status":impossible_loc.status,"used":impossible_loc.budget_used,
        "budget":impossible_loc.budget})

    flags=roots["method"]
    forbidden=("derivatives","jacobian","newton","quasi_newton","homotopy",
               "fallback","restarts","levenberg_marquardt","scipy")
    checks.append({"name":"method-audit","passed":all(not flags[k] for k in forbidden)
        and flags.get("halley") is True and flags.get("irp_analogue") is True
        and flags.get("true_residual_recompute") is True
        and flags.get("monomial_irp") is False,
        "method":flags})
    return {"script":Path(__file__).name,"self_test":True,
        "passed":all(c["passed"] for c in checks),"checks":checks}


def main(argv:Optional[Sequence[str]]=None)->int:
    args=parser().parse_args(argv); validate(args)
    if args.self_test: final=self_test(args); cases=["self-test"]
    else:
        cases=[x.strip() for x in str(args.cases).replace("|",";").split(";") if x.strip()]
        results=[run_case(args,c) for c in cases]; final=results[0] if len(results)==1 else {"script":Path(__file__).name,"cases":results}
    out=Path(args.out) if args.out else Path(args.outdir)/("self_test.json" if args.self_test else f"371_{cases[0].replace(',','x')}.json")
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(strict_json(final),indent=2,allow_nan=False),encoding="utf-8")
    if args.self_test:
        print(f"371 self-test: {'PASS' if final['passed'] else 'FAIL'} ({sum(int(c['passed']) for c in final['checks'])}/{len(final['checks'])})")
        print(f"out={out}"); return 0 if final["passed"] else 1
    for result in (results if len(cases)>1 else [final]):
        s=result["summary"]; print(f"371 case={result['case']} oracle={result['oracle']['kind']} roots={s['unique']}/{s['requested']} trials={s['trials']} samples={s['oracle_samples']} seconds={s['seconds']:.3f}")
    print(f"out={out}"); return 0


if __name__=="__main__": raise SystemExit(main())
