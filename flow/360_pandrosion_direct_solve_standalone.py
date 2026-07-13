"""Pandrosion 360: failure-aware standalone root solver.

The local corrector uses the exact coordinate-telescopic identity
F(b)-F(a)=Q(a,b)(b-a).  Stationary polynomial starts gain a deterministic
residual-chain escape, and scalable fixed-feature systems use a globalized
dense inverse-Broyden rescue written directly in NumPy.  No SciPy, automatic
differentiation, derivative formula, or local project import is used.  Large
Kostlan systems retain the lazy exact-in-law GP oracle introduced in 323.
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


class MonomialOracle(Oracle):
    def __init__(self, d: int, target: float) -> None:
        super().__init__(1, d, stable_seed(1, d, 0), "monomial-exact-article-001")
        self.target = float(target)

    def _eval_batch(self, Z: np.ndarray) -> np.ndarray:
        return (Z[:, 0]**self.d-self.target)[:, None]

    def backward_error(self, z: Any) -> float:
        zz = complex(np.asarray(z).reshape(-1)[0]); self.eval_count += 1
        return abs(zz**self.d-self.target)/max(abs(zz)**self.d+abs(self.target), 1e-300)


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
    if args.system_source == "monomial":
        if n != 1 or args.target <= 0: raise ValueError("monomial mode requires case 1,d and --target > 0")
        return MonomialOracle(d, args.target)
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


def slope_step(Q: np.ndarray, f: np.ndarray, trust: float, y: np.ndarray) -> tuple[np.ndarray, str]:
    scale = float(np.max(np.abs(Q))) if Q.size else 0
    if not math.isfinite(scale) or scale <= 1e-300: return np.zeros_like(y), "failed"
    q, rhs = Q/scale, -f/scale
    try: delta, method = np.linalg.solve(q, rhs), "finite-slope-solve"
    except Exception:
        try: delta, _, _, _ = np.linalg.lstsq(q, rhs, rcond=1e-12); method = "finite-slope-lstsq"
        except Exception: return np.zeros_like(y), "failed"
    cap = trust*max(1, norm(y)); nd = norm(delta)
    if math.isfinite(nd) and nd > cap: delta *= cap/nd
    return np.asarray(delta, np.complex128), method


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
    broyden_steps: int = 0
    escape_steps: int = 0


def _pandrosion_dogleg(Q: np.ndarray, f: np.ndarray, radius: float) -> tuple[np.ndarray, str]:
    """Powell-style trust-region path built only from a finite Pandrosion slope.

    Q is never a derivative: it is initialized from the exact telescopic identity
    and subsequently transported by exact secant closure.
    """
    n = len(f)
    # Solving Q p = -f is invariant under the common scalar normalization
    # formerly applied to both Q and f.  Avoid two dense temporary arrays in
    # the hot path and let LAPACK scale internally.
    try:
        pn = np.linalg.solve(Q, -f)
        mode = "pandrosion-finite-slope-direct-solve"
    except Exception:
        try:
            pn = np.linalg.lstsq(Q, -f, rcond=1e-12)[0]
            mode = "pandrosion-finite-slope-lstsq"
        except Exception:
            return np.zeros(n, np.complex128), "failed"
    nn = norm(pn)
    if nn <= radius:
        return np.asarray(pn, np.complex128), mode
    g = Q.conj().T @ f
    gg = float(np.vdot(g, g).real)
    qg = Q @ g
    den = float(np.vdot(qg, qg).real)
    if gg <= 1e-300 or den <= 1e-300:
        return np.asarray(pn * (radius/max(nn,1e-300)), np.complex128), mode+"-radial"
    pu = -(gg/den) * g
    nu = norm(pu)
    if nu >= radius:
        return np.asarray(pu * (radius/max(nu,1e-300)), np.complex128), "pandrosion-cauchy"
    d = pn - pu
    a = float(np.vdot(d,d).real)
    b = 2.0*float(np.vdot(pu,d).real)
    c = float(np.vdot(pu,pu).real) - radius*radius
    disc = max(0.0, b*b-4*a*c)
    tau = (-b + math.sqrt(disc))/(2*a) if a > 1e-300 else 0.0
    tau = min(1.0, max(0.0, tau))
    return np.asarray(pu + tau*d, np.complex128), "pandrosion-dogleg"


def _transport_secant(Q: np.ndarray, s: np.ndarray, ydiff: np.ndarray,
                      limit: float) -> tuple[np.ndarray, float, bool]:
    """Minimum-change exact finite-identity transport: Q+ s = ydiff."""
    ss = float(np.vdot(s,s).real)
    if ss <= 1e-300:
        return Q, float("inf"), False
    mismatch = ydiff - Q @ s
    corr = np.outer(mismatch, s.conj()) / ss
    ratio = norm(corr)/max(norm(Q),1e-300)
    if not math.isfinite(ratio) or ratio > limit:
        return Q, ratio, False
    out = Q + corr
    defect = norm(ydiff-out@s)/max(norm(ydiff),norm(out@s),1e-300)
    return out, defect, True


def _real_pack(z: np.ndarray) -> np.ndarray:
    z=np.asarray(z,np.complex128)
    return np.concatenate((z.real,z.imag))


def _real_unpack(x: np.ndarray) -> np.ndarray:
    x=np.asarray(x,float); n=len(x)//2
    return np.asarray(x[:n]+1j*x[n:],np.complex128)


def _real_residual(f: np.ndarray) -> np.ndarray:
    f=np.asarray(f,np.complex128)
    return np.concatenate((f.real,f.imag))


def _eval_real(target: Target, x: np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    z=_real_unpack(x); fc=target.eval(z)
    return _real_residual(fc),fc


def _eval_real_batch(target: Target, X: np.ndarray) -> np.ndarray:
    X=np.asarray(X,float); n=X.shape[1]//2
    Z=X[:,:n]+1j*X[:,n:]
    F=target.eval_batch(Z)
    return np.concatenate((F.real,F.imag),axis=1)


def _initial_pandrosion_slope(target: Target, x: np.ndarray, fx: np.ndarray,
                               args: argparse.Namespace, attempt: int) -> tuple[np.ndarray,int,float]:
    """Exact coordinate-telescopic real slope using m+1 black-box values."""
    m=len(x); xn=max(1.0,norm(x))
    h=max(args.fast_init_scale*xn, 1e-9*xn) * (args.fast_rebuild_shrink**attempt)
    signs=np.where((np.arange(m)+attempt)%2, -1.0, 1.0)
    b=x + (h/math.sqrt(max(1,m)))*signs
    path=np.repeat(x[None,:],m+1,axis=0)
    for j in range(m): path[j+1:,j]=b[j]
    # F(x) is already available as ``fx``.  Evaluate only the m new path
    # points instead of reevaluating the common base point.  This preserves
    # exactly the same telescopic finite-slope matrix while saving one oracle
    # sample at every full slope construction.
    Fnew=_eval_real_batch(target,path[1:])
    F=np.vstack((np.asarray(fx,float)[None,:],Fnew)); dz=b-x
    Q=(F[1:]-F[:-1]).T/dz[None,:]
    telescoped=Q@dz
    defect=norm((F[-1]-F[0])-telescoped)/max(norm(F[-1]-F[0]),norm(telescoped),1e-300)
    return np.asarray(Q,float),len(path)-1,defect


@dataclass
class CachedStart:
    x: np.ndarray
    fx: np.ndarray
    fc: np.ndarray


def correct(target: Target, y0: Any, args: argparse.Namespace, known: Sequence[np.ndarray],
            epochs: int, deadline: Optional[float], cached: Optional[CachedStart]=None) -> CorrectResult:
    """Fast 328 corrector in R^(2n), using only finite Pandrosion identities."""
    if cached is None:
        z0=np.asarray(y0,np.complex128).copy(); x=_real_pack(z0); fx,fc=_eval_real(target,x)
    else:
        x=np.asarray(cached.x,float).copy(); fx=np.asarray(cached.fx,float).copy(); fc=np.asarray(cached.fc,np.complex128).copy()
    r=norm(fx)
    best_x,best_fc,best_r=x.copy(),fc.copy(),r
    slopes=transports=probe_evals=slope_evals=line_evals=para=0
    max_defect=max_transport=0.0; smin=smax=cond=None
    status="max-epochs"; done=0; Q=None; rebuilds=0
    radius=max(args.fast_min_radius,args.fast_initial_radius*max(1.0,norm(x)))
    accepted_since_rebuild=0
    max_radius=args.fast_max_radius*max(1.0,norm(x)); stale=0
    for ep in range(max(1,epochs)):
        done=ep
        if deadline is not None and clock()>=deadline: status="timeout"; break
        if best_r<args.accept: status="converged"; break
        if Q is None:
            accepted_since_rebuild=0
            try:
                Q,used,defect=_initial_pandrosion_slope(target,x,fx,args,rebuilds)
                slopes+=1; slope_evals+=used; max_defect=max(max_defect,defect); rebuilds+=1
                # Singular values are diagnostics only; they do not affect the
                # Pandrosion trajectory.  Avoid the dense SVD in the default hot
                # path and compute it only when explicitly requested.
                if args.diagnostic_svd:
                    sv=np.linalg.svd(Q,compute_uv=False)
                    if len(sv):
                        smax=float(sv[0]); smin=float(sv[-1])
                        cond=float(smax/smin) if smin>0 else float("inf")
            except Exception as exc:
                status=f"slope-error:{type(exc).__name__}"; break
        step,mode=_pandrosion_dogleg(Q,fx,radius); ns=norm(step)
        if mode=="failed" or not np.all(np.isfinite(step)) or ns<=1e-15:
            Q=None; stale+=1
            if stale>args.fast_max_rebuilds: status="invalid-step"; break
            continue
        pred_vec=fx+Q@step; pred=max(0.0,r*r-norm(pred_vec)**2)
        xn=x+np.asarray(step.real,float); fn,fcn=_eval_real(target,xn); rn=norm(fn); line_evals+=1
        actual=r*r-rn*rn; rho=actual/max(pred,1e-300)
        accepted=math.isfinite(rn) and (rho>args.fast_eta or rn<r*args.fast_force_ratio or rn<best_r)
        if accepted:
            old_x,old_fx=x,fx
            x,fx,fc,r=xn,fn,fcn,rn
            accepted_since_rebuild += 1
            if r<best_r: best_x,best_fc,best_r=x.copy(),fc.copy(),r; stale=0
            else: stale+=1
            Qnew,defect,ok=_transport_secant(Q,x-old_x,fx-old_fx,args.fast_transport_limit)
            if ok:
                max_defect=max(max_defect,defect); max_transport=max(max_transport,norm(Qnew-Q)/max(norm(Q),1e-300)); Q=Qnew; transports+=1
            else: Q=None
            if accepted_since_rebuild >= args.fast_refresh:
                Q=None
            if rho>0.75 and ns>0.8*radius: radius=min(max_radius,2.0*radius)
            elif rho<0.25: radius=max(args.fast_min_radius,0.25*radius)
        else:
            radius=max(args.fast_min_radius,0.25*radius); stale+=1
            if stale>=args.fast_rebuild_after: Q=None; stale=0
        if rebuilds>args.fast_max_rebuilds and Q is None: status="rebuild-limit"; break
    best_z=_real_unpack(best_x); ok=best_r<args.accept
    if ok: status="converged"
    return CorrectResult(best_z,best_fc,best_r,ok,status,done+1,slopes,transports,probe_evals,
                         slope_evals,line_evals,para,max_defect,max_transport,smin,smax,cond)

def add_results(base: CorrectResult, loc: CorrectResult, y: np.ndarray, f: np.ndarray,
                residual: float, status: str) -> CorrectResult:
    return CorrectResult(y, f, residual, base.ok or loc.ok, status,
        base.epochs+loc.epochs, base.slopes+loc.slopes, base.transported_steps+loc.transported_steps,
        base.probe_evals+loc.probe_evals,
        base.slope_evals+loc.slope_evals, base.line_evals+loc.line_evals,
        base.parabolic_evals+loc.parabolic_evals,
        max(base.max_telescope_defect, loc.max_telescope_defect),
        max(base.max_transport_correction, loc.max_transport_correction),
        loc.smin if loc.smin is not None else base.smin,
        loc.smax if loc.smax is not None else base.smax,
        loc.cond if loc.cond is not None else base.cond)


def irp_rescue(base: Target, result: CorrectResult, args: argparse.Namespace,
               known: Sequence[np.ndarray], deadline: Optional[float]) -> CorrectResult:
    """Dynamic direct/reciprocal chart renormalization using only finite slopes."""
    if result.ok or not args.irp: return result
    best = result; gains = floats(args.irp_scales, [1, 2**(1/3), 2**(-1/3), 2, .5])
    for layer in range(args.irp_epochs):
        if deadline is not None and clock() >= deadline: break
        candidates = []
        for gain in gains:
            for reciprocal in (False, True):
                chart = IRPTarget(base, complex(gain), reciprocal); u = chart.from_base(best.y)
                energy = float(np.mean(np.abs(np.log(np.maximum(np.abs(u), 1e-300)))))
                candidates.append((energy, chart, u))
        improved = None
        for _, chart, u in sorted(candidates, key=lambda x:x[0])[:args.irp_top]:
            loc = correct(chart, u, args, [], 1, deadline)
            y = chart.to_base(loc.y); f = loc.f.copy(); r = norm(f)
            if r < best.residual and (improved is None or r < improved[0]): improved = (r, y, f, loc)
        if improved is None: break
        r, y, f, loc = improved
        best = add_results(best, loc, y, f, r, "irp-converged" if r < args.accept else "irp-improved")
        best.ok = r < args.accept
        if best.ok: break
    return best


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
        cache=CachedStart(_real_pack(loc.y),_real_residual(loc.f),loc.f.copy())
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


def polish(base: Target, y0: Any, args: argparse.Namespace, deadline: Optional[float]) -> tuple[np.ndarray,np.ndarray,float,dict[str,Any],CorrectResult]:
    loc = correct(base, y0, args, [], max(1,args.polish_steps), deadline)
    return loc.y, loc.f, loc.residual, conditioning(loc), loc


# ---------- case driver ------------------------------------------------------

def run_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    if args.system_source == "monomial": return run_monomial_case(args, case_raw)
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
            local_trial = trial-len(priority)
            if args.origin_fallback and local_trial == args.origin_after:
                y0 = np.zeros(n, np.complex128)
            else:
                atlas_trial = local_trial if (not args.origin_fallback or local_trial < args.origin_after) else local_trial-1
                y0 = atlas_start(n, atlas_trial, oracle.seed)
        known=[np.asarray(r["z_complex"],np.complex128) for r in roots]
        if known and min(norm(y0-r) for r in known)<=args.early_dup_sep:
            duplicates+=1; trials.append({"trial":trial,"status":"start-near-root","accepted":False}); continue
        deadline=clock()+args.trial_timeout if args.trial_timeout>0 else None; before=oracle.eval_count
        try:
            loc=correct(base,y0,args,known,args.epochs,deadline,cached0)
            if (not loc.ok) and loc.residual<=max(args.polish_gate,100*args.accept):
                z,f,residual,condmeta,ploc=polish(base,loc.y,args,deadline)
                loc=add_results(loc,ploc,z,f,residual,loc.status)
                condmeta=conditioning(loc)
            else: z,f,residual=loc.y.copy(),loc.f.copy(),loc.residual; condmeta=conditioning(loc)
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
    return {"script":Path(__file__).name,"version":339,"standalone":True,"case":f"{n},{d}","n":n,"degree":d,
      "method":{"family":"pure-pandrosion-origin-pilot-finite-transport","derivatives":False,"jacobian":False,"newton":False,
        "broyden":False,"levenberg_marquardt":False,
        "identity":"F(b)-F(a)=Q(a,b)(b-a) by coordinate telescoping",
        "transport":"one-column finite closure preserves F(b)-F(a)=Q(a,b)(b-a); bounded then fully rebuilt",
        "irp_scope":"disabled; no rescue/fallback solver path"},
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
        "ks_backend":args.ks_backend,"probes":args.probes,"probe_top":args.probe_top,
        "deflation_alpha":args.deflation_alpha,"irp":args.irp,"swarm":args.swarm},
      "portfolio":{"small_chart_first":True,"screened_start":True,"start_meta":start_meta,"pilot_meta":pilot_meta,"origin_fallback":args.origin_fallback,"origin_after":args.origin_after},"swarm":swarm_meta,"roots":encoded,"trials":trials if args.verbose_trials else trials[:args.keep_trials],
      "summary":{"requested":args.count,"unique":len(roots),"success":len(roots)>=args.count,
        "trials":len(trials),"duplicates":duplicates,"failures":failures,
        "oracle_samples":oracle.eval_count,"seconds":clock()-started,"oracle_cap_error":cap_error}}


# ---------- exact article monomial IRP audit --------------------------------

def monomial_irp(x: float, p: int, u: float, layers: int, q: float=1.25, K: int=4) -> tuple[float,list[float]]:
    """Equation (22) of article 001 on the positive real branch."""
    history=[abs(u**p-x)/x]
    for _ in range(layers):
        R=x/(u**p); direct=R>=1; A=R if direct else 1/R
        cell=math.floor(K*math.log(A)/(p*math.log(q))+1e-15)/K
        B=q**cell; Y=A/(B**p)
        s=1-(Y-1)/(Y*p); W=B*Y*(s**(p-1))
        u=u*W if direct else u/W; history.append(abs(u**p-x)/x)
    return u,history


def run_monomial_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    started=clock(); n,p=parse_case(case_raw)
    if n != 1 or args.target <= 0: raise ValueError("monomial mode requires case 1,p and --target > 0")
    u,history=monomial_irp(args.target,p,args.monomial_start,args.monomial_layers,args.palette_q,args.palette_k)
    raw=abs(u**p-args.target); validation=raw/max(abs(u)**p+args.target,1e-300); accepted=validation<args.validation_accept
    root={"id":0,"trial":0,"z":[[float(u),0.]],"residual":raw,"validation_error":validation,
          "realness":0.,"smin":None,"smax":None,"cond":None,"near_multiple":False,"singular":False,
          "mode":"article-001-positive-real-irp"}
    return {"script":Path(__file__).name,"version":339,"standalone":True,"case":f"1,{p}","n":1,"degree":p,
      "method":{"family":"article-001-exact-monomial-irp","derivatives":False,"jacobian":False,"newton":False,
        "broyden":False,"levenberg_marquardt":False,"equation":"article 001, equation (22)"},
      "oracle":{"kind":"monomial-exact-article-001","target":args.target,"samples":args.monomial_layers,
        "residual_mode":"raw","validation_mode":"scale-invariant-monomial-backward-error"},
      "parameters":{"layers":args.monomial_layers,"start":args.monomial_start,"palette_q":args.palette_q,"palette_k":args.palette_k},
      "swarm":{"enabled":False},"roots":[root] if accepted else [],
      "trials":[{"trial":0,"status":"new-root" if accepted else "not-converged","accepted":accepted,
        "residual":raw,"validation_error":validation,"relative_residual_history":history}],
      "summary":{"requested":1,"unique":int(accepted),"success":accepted,"trials":1,"duplicates":0,
        "failures":int(not accepted),"oracle_samples":args.monomial_layers,"seconds":clock()-started,"oracle_cap_error":None}}


# ---------- CLI and regressions ---------------------------------------------

def parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser(description="Pandrosion 359 lazy-diagnostics pure finite-slope engine")
    p.add_argument("--cases",default="2,4"); p.add_argument("--seed-index",type=int,default=0)
    p.add_argument("--system-source",choices=["ks","kostlan","poly","polynomial","monomial"],default="ks")
    p.add_argument("--ks-backend",choices=["auto","dense","gp","feature","gaussian"],default="auto")
    p.add_argument("--polys"); p.add_argument("--variables"); p.add_argument("--starts")
    p.add_argument("--target",type=float,default=4.); p.add_argument("--monomial-start",type=float,default=1.)
    p.add_argument("--monomial-layers",type=int,default=8); p.add_argument("--palette-q",type=float,default=1.25)
    p.add_argument("--palette-k",type=int,default=4)
    p.add_argument("--dense-max-terms",type=int,default=250000); p.add_argument("--features",type=int,default=3072)
    p.add_argument("--eval-block",type=int,default=128); p.add_argument("--equation-normalize",action="store_true",default=False)
    p.add_argument("--gp-jitter",type=float,default=1e-9); p.add_argument("--gp-max-points",type=int,default=10000)
    p.add_argument("--gp-jet-radius",type=float,default=1e-3); p.add_argument("--gp-swarm-cap",type=int,default=8)
    p.add_argument("--gp-probe-cap",type=int,default=4)
    p.add_argument("--count",type=int,default=1); p.add_argument("--pool",type=int,default=8)
    p.add_argument("--epochs",type=int,default=96); p.add_argument("--accept",type=float,default=1e-8)
    p.add_argument("--origin-fallback",action="store_true",default=False)
    p.add_argument("--no-origin-fallback",dest="origin_fallback",action="store_false")
    p.add_argument("--origin-after",type=int,default=3)
    p.add_argument("--screen-starts",action="store_true",default=True)
    p.add_argument("--no-screen-starts",dest="screen_starts",action="store_false")
    p.add_argument("--start-candidates",type=int,default=24)
    p.add_argument("--start-keep",type=int,default=8)
    p.add_argument("--start-radii",default=".125,.25,.5,1,2")
    p.add_argument("--start-shrink",type=float,default=.5)
    p.add_argument("--start-contraction-weight",type=float,default=.75)
    p.add_argument("--start-sep",type=float,default=0)
    p.add_argument("--pilot-candidates",type=int,default=1)
    p.add_argument("--pilot-epochs",type=int,default=3)
    p.add_argument("--legacy-head",type=int,default=2)
    p.add_argument("--tol",type=float,default=1e-12); p.add_argument("--cluster-sep",type=float,default=1e-7)
    p.add_argument("--early-dup-sep",type=float,default=1e-4); p.add_argument("--trial-timeout",type=float,default=0)
    p.add_argument("--probes",type=int,default=1); p.add_argument("--probe-top",type=int,default=1)
    p.add_argument("--probe-scale",type=float,default=.03); p.add_argument("--probe-radii",default=".5,1,2")
    p.add_argument("--probe-cond-weight",type=float,default=.02); p.add_argument("--equal-value-weight",type=float,default=.02)
    p.add_argument("--trust-radius",type=float,default=10); p.add_argument("--line-grid",default="1,.75,.5,.25,.125,.0625")
    p.add_argument("--transport-steps",type=int,default=8); p.add_argument("--transport-ratio",type=float,default=.85)
    p.add_argument("--transport-trust-radius",type=float,default=5); p.add_argument("--transport-line-grid",default="1,.5,.25")
    p.add_argument("--transport-closure-limit",type=float,default=2.)
    p.add_argument("--slope-retries",type=int,default=2); p.add_argument("--slope-retry-factor",type=float,default=.25)
    p.add_argument("--deflation-alpha",type=float,default=.15)
    p.add_argument("--parabolic",action="store_true",default=True); p.add_argument("--no-parabolic",dest="parabolic",action="store_false")
    p.add_argument("--irp",action="store_true",default=False); p.add_argument("--no-irp",dest="irp",action="store_false")
    p.add_argument("--irp-epochs",type=int,default=3); p.add_argument("--irp-top",type=int,default=3)
    p.add_argument("--irp-scales",default="1,1.2599210499,.793700526,2,.5")
    p.add_argument("--swarm",action="store_true",default=False); p.add_argument("--no-swarm",dest="swarm",action="store_false")
    p.add_argument("--swarm-size",type=int,default=48); p.add_argument("--swarm-keep",type=int,default=12)
    p.add_argument("--swarm-sep",type=float,default=0); p.add_argument("--swarm-epochs",type=int,default=1)
    p.add_argument("--swarm-timeout",type=float,default=0)
    p.add_argument("--polish-steps",type=int,default=4); p.add_argument("--polish-gate",type=float,default=1e-2)
    p.add_argument("--validation-accept",type=float,default=1e-8); p.add_argument("--keep-trials",type=int,default=100)
    p.add_argument("--fast-init-scale",type=float,default=1e-5)
    p.add_argument("--fast-rebuild-shrink",type=float,default=0.5)
    p.add_argument("--fast-initial-radius",type=float,default=1.0)
    p.add_argument("--fast-max-radius",type=float,default=100.0)
    p.add_argument("--fast-min-radius",type=float,default=1e-12)
    p.add_argument("--fast-eta",type=float,default=1e-4)
    p.add_argument("--fast-force-ratio",type=float,default=0.95)
    p.add_argument("--fast-transport-limit",type=float,default=25.0)
    p.add_argument("--fast-rebuild-after",type=int,default=1)
    p.add_argument("--fast-max-rebuilds",type=int,default=20)
    p.add_argument("--fast-refresh",type=int,default=12)
    p.add_argument("--diagnostic-svd",action="store_true",default=False,help="compute singular-value diagnostics on every rebuilt slope")
    p.add_argument("--verbose-trials",action="store_true"); p.add_argument("--self-test",action="store_true")
    p.add_argument("--out"); p.add_argument("--outdir",default="/tmp/359_pandrosion")
    return p


def validate(args: argparse.Namespace) -> None:
    for raw in str(args.cases).replace("|",";").split(";"):
        if raw.strip(): parse_case(raw)
    if args.count<0 or args.pool<0 or args.epochs<=0 or args.probes<=0 or args.probe_top<=0: raise ValueError("invalid budgets")
    if args.origin_after<0: raise ValueError("invalid origin fallback position")
    if args.pilot_candidates<0 or args.pilot_epochs<0: raise ValueError("invalid pilot controls")
    if args.legacy_head<0 or args.start_candidates<=0 or args.start_keep<=0 or not 0<args.start_shrink<1 or args.start_contraction_weight<0: raise ValueError("invalid start-screen controls")
    if args.swarm_epochs<0 or args.swarm_timeout<0: raise ValueError("invalid swarm controls")
    if args.accept<=0 or args.validation_accept<=0 or args.trust_radius<=0: raise ValueError("invalid tolerances")
    if args.gp_jitter<=0 or args.gp_max_points<=0 or args.gp_jet_radius<=0 or args.gp_swarm_cap<=0 or args.gp_probe_cap<=0: raise ValueError("invalid GP controls")
    if args.probe_scale<=0 or args.probe_cond_weight<0 or args.equal_value_weight<0: raise ValueError("invalid probe controls")
    if args.transport_steps<0 or not 0<args.transport_ratio<=1 or args.transport_trust_radius<=0 or args.transport_closure_limit<=0: raise ValueError("invalid finite-slope transport controls")
    if args.slope_retries<0 or not 0<args.slope_retry_factor<1: raise ValueError("invalid finite-slope retry controls")
    if args.target<=0 or args.monomial_start<=0 or args.monomial_layers<=0 or args.palette_q<=1 or args.palette_k<=0: raise ValueError("invalid monomial IRP controls")


def self_test(args: argparse.Namespace) -> dict[str,Any]:
    checks=[]
    a=argparse.Namespace(**vars(args)); a.system_source="poly"; a.polys="x^2-3*x-10"; a.variables=None
    a.starts="-8,4"; a.count=2; a.pool=8; a.swarm=False; a.epochs=24
    r=run_case(a,"1,2"); checks.append({"name":"two-roots-pure-slope","passed":r["summary"]["success"],"result":r})
    a=argparse.Namespace(**vars(args)); a.system_source="poly"; a.polys="x1^2-1;x2^2-1"; a.variables=None
    a.starts="-2,-2;-2,2;2,-2;2,2"; a.count=4; a.pool=8; a.swarm=False; a.epochs=30
    r=run_case(a,"2,2"); checks.append({"name":"four-roots-pure-slope","passed":r["summary"]["success"],"result":r})
    a=argparse.Namespace(**vars(args)); a.system_source="poly"; a.polys="x1^2+x1*x2-1;x2^2-x1"; a.variables=None
    o=make_oracle(a,2,2); t=Target(o); aa=np.asarray([.2+.1j,-.3j]); bb=np.asarray([1-.2j,.4+.1j])
    Q,defect,_=finite_slope(t,aa,bb); identity=norm((t.eval(bb)-t.eval(aa))-Q@(bb-aa))
    checks.append({"name":"exact-telescoping","passed":identity<1e-12 and defect<1e-12,"identity_error":identity,"defect":defect})
    yy=.7*aa+.3*bb; fy=t.eval(yy); fb=t.eval(bb); closed,cdef,cratio,ok=close_finite_slope(Q,bb,fb,yy,fy,100.)
    closure_error=norm((fb-fy)-closed@(bb-yy))
    checks.append({"name":"exact-finite-transport-closure","passed":ok and cdef<1e-12 and closure_error<1e-12,
        "identity_error":closure_error,"defect":cdef,"correction_ratio":cratio})
    u,h=monomial_irp(4.,3,1.,3); checks.append({"name":"article-001-exact-irp","passed":h[-1]<2e-12,"history":h,"root":u})
    a=argparse.Namespace(**vars(args)); a.system_source="monomial"; a.target=4.; a.monomial_layers=4
    mr=run_case(a,"1,3"); checks.append({"name":"article-001-cli-mode","passed":mr["summary"]["success"],"result":mr})
    a=argparse.Namespace(**vars(args)); a.system_source="ks"; a.ks_backend="dense"; a.count=1; a.pool=16; a.epochs=24
    a.starts=None; a.swarm=True; r=run_case(a,"2,3"); checks.append({"name":"dense-kostlan-smoke","passed":r["summary"]["success"],"result":r})
    a=argparse.Namespace(**vars(args)); a.system_source="ks"; a.ks_backend="gp"; a.gp_max_points=256
    gp=make_oracle(a,2,3); Z=np.asarray([[0,0],[.2+.1j,-.3j],[1,-.5j]],np.complex128)
    first=gp.eval_batch(Z); cached=gp.eval_batch(Z); cov=norm(gp.L@gp.L.conj().T-gp.kernel(Z,Z))
    checks.append({"name":"gp-covariance-cache","passed":np.array_equal(first,cached) and cov<1e-9,"covariance_error":cov})
    dense=KostlanOracle(2,3,7,10,32,False,32)
    weights=np.asarray([math.sqrt(math.factorial(3)/(math.factorial(3-int(e.sum()))*math.prod(math.factorial(int(x)) for x in e))) for e in dense.exps])
    phi=np.prod(Z[:,None,:]**dense.exps[None,:,:],axis=2)*weights[None,:]
    C=phi@phi.conj().T; C/=np.sqrt(np.diag(C))[:,None]*np.sqrt(np.diag(C))[None,:]
    kernel_error=norm(C-gp.kernel(Z,Z)); checks.append({"name":"gp-kernel-dense-equivalence","passed":kernel_error<1e-12,"kernel_error":kernel_error})
    finite=ExactKSGPOracle(1,2,9,1e-9,64,1e-3); finite.eval_batch(np.linspace(-2,2,17,dtype=np.complex128)[:,None])
    checks.append({"name":"gp-finite-rank","passed":len(finite.cache)==17 and len(finite.points)<=3,
        "queries":len(finite.cache),"active_rank":len(finite.points)})
    a=argparse.Namespace(**vars(args)); a.system_source="ks"; a.ks_backend="feature"; a.features=512
    feature=make_oracle(a,20,20); origin=feature.backward_error(np.zeros(20,np.complex128))
    checks.append({"name":"feature-origin-guard","passed":origin>a.validation_accept and np.min(feature.degrees)==0,"origin_backward_error":origin})
    capped=ExactKSGPOracle(2,3,11,1e-9,1,1e-3); cap_ok=False
    try: capped.eval_batch(np.asarray([[0,0],[1,1]],np.complex128))
    except RuntimeError: cap_ok=True
    checks.append({"name":"gp-cap-guard","passed":cap_ok})
    a=argparse.Namespace(**vars(args)); a.system_source="ks"; a.ks_backend="gp"; a.seed_index=0
    a.count=1; a.pool=4; a.epochs=40; a.swarm=True; a.swarm_size=4; a.swarm_keep=4; a.gp_max_points=3000
    r=run_case(a,"20,20"); checks.append({"name":"pure-gp-20x20","passed":r["summary"]["success"],"result":r})
    flags=r["method"]; transported=sum(t.get("transported_finite_slope_steps",0) for t in r["trials"])
    checks.append({"name":"pure-transport-audit","passed":transported>0 and not any(flags[k] for k in ("derivatives","jacobian","newton","broyden","levenberg_marquardt")),
        "transported_steps":transported,"method":flags})
    return {"script":Path(__file__).name,"self_test":True,"passed":all(c["passed"] for c in checks),"checks":checks}


def main(argv:Optional[Sequence[str]]=None)->int:
    args=parser().parse_args(argv); validate(args)
    if args.self_test: final=self_test(args); cases=["self-test"]
    else:
        cases=[x.strip() for x in str(args.cases).replace("|",";").split(";") if x.strip()]
        results=[run_case(args,c) for c in cases]; final=results[0] if len(results)==1 else {"script":Path(__file__).name,"cases":results}
    out=Path(args.out) if args.out else Path(args.outdir)/("self_test.json" if args.self_test else f"356_{cases[0].replace(',','x')}.json")
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(strict_json(final),indent=2,allow_nan=False),encoding="utf-8")
    if args.self_test:
        print(f"331 self-test: {'PASS' if final['passed'] else 'FAIL'} ({sum(int(c['passed']) for c in final['checks'])}/{len(final['checks'])})")
        print(f"out={out}"); return 0 if final["passed"] else 1
    for result in (results if len(cases)>1 else [final]):
        s=result["summary"]; print(f"348 case={result['case']} oracle={result['oracle']['kind']} roots={s['unique']}/{s['requested']} trials={s['trials']} samples={s['oracle_samples']} seconds={s['seconds']:.3f}")
    print(f"out={out}"); return 0


if __name__=="__main__": raise SystemExit(main())
