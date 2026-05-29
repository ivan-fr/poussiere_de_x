#!/usr/bin/env python3
"""
316_pandrosion_standalone_exact_system_geometry_irp_hypercube_inversejet_numpy_engine.py

Standalone NumPy engine with the same operational surface as 315, but without
loading 314 or 315 by file path.  It supports:

  - KS/Kostlan systems: dense for small cases, lazy-feature or geometry-kernel
    for large cases;
  - exact user polynomial systems with optional fitted geometry-kernel starts;
  - explicit starts, universal/projective starts, StartOpt, 306-style hypercube
    inverse-jet correction, lazy IRP chart rescue, exact residual validation;
  - JSON export compatible with the 315 summary shape.

Dependencies: Python stdlib + NumPy only.  No local flow imports.
"""
from __future__ import annotations

import argparse
import ast
import dataclasses
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def now() -> float:
    return time.time()


def cjson(z: complex) -> list[float]:
    w = complex(z)
    return [float(w.real), float(w.imag)]


def root_to_json(z: Sequence[complex]) -> list[list[float]]:
    return [cjson(v) for v in z]


def parse_case(raw: str) -> tuple[int, int]:
    s = str(raw).strip().lower().replace("x", ",").replace(":", ",")
    p = [q.strip() for q in s.split(",") if q.strip()]
    if len(p) != 2:
        raise ValueError(f"case must be n,d, got {raw!r}")
    n, d = int(p[0]), int(p[1])
    if n <= 0 or d <= 0:
        raise ValueError("n,d must be positive")
    return n, d


def parse_float_list(raw: Optional[str], default: Sequence[float], positive: bool = False) -> list[float]:
    if raw is None or str(raw).strip() == "":
        return list(default)
    out: list[float] = []
    for part in str(raw).replace(";", ",").split(","):
        try:
            x = float(part.strip())
        except Exception:
            continue
        if math.isfinite(x) and ((not positive) or x > 0):
            out.append(x)
    return out or list(default)


def splitmix64(x: int) -> int:
    x = (int(x) + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return (x ^ (x >> 31)) & 0xFFFFFFFFFFFFFFFF


def u01(x: int) -> float:
    return ((splitmix64(x) >> 11) & ((1 << 53) - 1)) / float(1 << 53)


def phase(theta: float) -> complex:
    return complex(math.cos(theta), math.sin(theta))


def stable_seed(n: int, d: int, seed_index: int = 0, salt: int = 0) -> int:
    return int(splitmix64(0x50414E44524F5349 + 1000003 * n + 9176 * d + 97 * seed_index + salt) & 0x7FFFFFFF)


def exact_kostlan_terms(n: int, d: int) -> int:
    return int(math.comb(int(n) + int(d), int(d)))


def auto_lazy_feature_count(n: int, d: int, cap: int) -> int:
    base = int(32 * math.sqrt(max(1, (int(n) + 1) * (int(d) + 1))))
    return int(max(int(n) + 8, min(max(1, int(cap)), max(512, base))))


def auto_geometry_anchor_count(n: int, d: int, cap: int) -> int:
    base = int(24 * math.sqrt(max(1, (int(n) + 1) * (int(d) + 1))))
    return int(max(int(n) + 16, min(max(1, int(cap)), max(512, base))))


def finite_norm(v: Any) -> float:
    try:
        r = float(np.linalg.norm(v))
        return r if math.isfinite(r) else float("inf")
    except Exception:
        return float("inf")


def log_energy(y: Sequence[complex]) -> float:
    yy = np.asarray(y, dtype=np.complex128)
    if yy.size == 0:
        return 0.0
    return float(np.mean(np.abs(np.log(np.maximum(np.abs(yy), 1e-300)))))


def realness(z: Sequence[complex]) -> float:
    zz = np.asarray(z, dtype=np.complex128)
    den = max(1e-300, float(np.linalg.norm(zz)))
    return float(np.linalg.norm(zz.imag) / den)


def score_root(residual: float, realness_value: float, cond: Optional[float]) -> float:
    c = float(cond) if cond is not None and math.isfinite(float(cond)) else 1e300
    return float(math.log1p(max(0.0, residual)) + 0.01 * math.log1p(max(0.0, realness_value)) + 0.001 * math.log1p(max(0.0, c)))


def cluster_index(roots: list[dict[str, Any]], z: Any, sep: float) -> Optional[int]:
    zz = np.asarray(z, dtype=np.complex128)
    for i, r in enumerate(roots):
        if float(np.linalg.norm(zz - np.asarray(r["z_complex"], dtype=np.complex128))) <= float(sep):
            return i
    return None


def parse_complex_token(raw: str) -> complex:
    s = str(raw).strip().replace("i", "j")
    if not s:
        raise ValueError("empty complex token")
    return complex(s)


def parse_start_points(raw: Optional[str], n: int) -> list[Any]:
    if raw is None or str(raw).strip() == "":
        return []
    text = str(raw).strip()
    if int(n) == 1 and ";" not in text and "|" not in text:
        return [np.asarray([parse_complex_token(p)], dtype=np.complex128) for p in text.split(",") if p.strip()]
    pts = []
    for part in text.replace("|", ";").split(";"):
        if not part.strip():
            continue
        coords = [p.strip() for p in part.split(",") if p.strip()]
        if len(coords) != int(n):
            raise ValueError(f"start point {part!r} has {len(coords)} coord(s), expected {n}")
        pts.append(np.asarray([parse_complex_token(p) for p in coords], dtype=np.complex128))
    return pts


# ---------------------------------------------------------------------------
# KS systems
# ---------------------------------------------------------------------------

def compositions_leq(d: int, n: int) -> Any:
    out: list[tuple[int, ...]] = []

    def rec(pos: int, rem: int, cur: list[int]) -> None:
        if pos == n - 1:
            for k in range(rem + 1):
                out.append(tuple(cur + [k]))
            return
        for k in range(rem + 1):
            cur.append(k)
            rec(pos + 1, rem - k, cur)
            cur.pop()

    rec(0, int(d), [])
    return np.asarray(out, dtype=np.int16 if d < 32767 else np.int32)


def multinomial_kostlan_weights(exps: Any, d: int) -> Any:
    totals = np.sum(exps, axis=1).astype(np.int64)
    logfac = np.zeros(int(d) + 1, dtype=np.float64)
    for k in range(1, int(d) + 1):
        logfac[k] = logfac[k - 1] + math.log(k)
    logs = logfac[d] - logfac[d - totals]
    for j in range(exps.shape[1]):
        logs -= logfac[exps[:, j].astype(np.int64)]
    return np.exp(0.5 * logs)


@dataclasses.dataclass
class DenseKostlanSystem:
    n: int
    d: int
    seed: int
    exps: Any
    coeff: Any
    weights: Any
    equation_normalize: bool = True
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, seed_index: int = 0, equation_normalize: bool = True) -> "DenseKostlanSystem":
        t0 = now()
        exps = compositions_leq(d, n)
        weights = multinomial_kostlan_weights(exps, d)
        seed = stable_seed(n, d, seed_index)
        rng = np.random.default_rng(seed)
        coeff = (rng.standard_normal((n, exps.shape[0])) + 1j * rng.standard_normal((n, exps.shape[0]))) / math.sqrt(2.0)
        coeff = coeff * weights[None, :]
        if equation_normalize:
            row = np.linalg.norm(coeff, axis=1)
            coeff = coeff / np.where(row > 0, row, 1.0)[:, None]
        obj = cls(int(n), int(d), seed, exps, coeff.astype(np.complex128), weights, bool(equation_normalize))
        obj._generation_seconds = now() - t0
        return obj

    @property
    def terms_per_poly(self) -> int:
        return int(self.exps.shape[0])

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return float(getattr(self, "_generation_seconds", 0.0))

    def monomials_batch(self, Z: Any) -> Any:
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        B, M = int(ZZ.shape[0]), self.terms_per_poly
        mon = np.ones((B, M), dtype=np.complex128)
        for j in range(self.n):
            p = np.empty((B, self.d + 1), dtype=np.complex128)
            p[:, 0] = 1.0
            if self.d > 0:
                p[:, 1] = ZZ[:, j]
                for k in range(2, self.d + 1):
                    p[:, k] = p[:, k - 1] * ZZ[:, j]
            mon *= p[:, self.exps[:, j]]
        return mon

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = (self.monomials_batch(np.asarray(z, dtype=np.complex128))[0] @ self.coeff.T)
        self.eval_count += 1
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        out = self.monomials_batch(ZZ) @ self.coeff.T
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms}


@dataclasses.dataclass
class LazyFeatureKostlanSystem:
    n: int
    d: int
    seed: int
    feature_exps: Any
    feature_exps_t: Any
    feature_log_scales: Any
    coeff: Any
    lazy_features: int
    projective_normalize: bool = True
    dynamic_normalize: bool = True
    equation_normalize: bool = True
    eval_block: int = 128
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, seed_index: int = 0, equation_normalize: bool = True, lazy_features: int = 0, lazy_feature_cap: int = 8192, projective_normalize: bool = True, dynamic_normalize: bool = True, eval_block: int = 128) -> "LazyFeatureKostlanSystem":
        t0 = now()
        n, d = int(n), int(d)
        m = max(n + 1, int(lazy_features) if int(lazy_features) > 0 else auto_lazy_feature_count(n, d, int(lazy_feature_cap)))
        seed = stable_seed(n, d, seed_index, salt=0x314314)
        rng = np.random.default_rng(seed)
        exps = np.zeros((m, n), dtype=np.int16 if d < 32767 else np.int32)
        degrees = np.zeros(m, dtype=np.int64)
        idx = 1
        for j in range(n):
            if idx >= m:
                break
            exps[idx, j] = 1 if d >= 1 else 0
            degrees[idx] = 1 if d >= 1 else 0
            idx += 1
        probs = np.full(n, 1.0 / max(1, n))
        while idx < m:
            q = (idx - (n + 1) + 0.5) / max(1, m - (n + 1))
            k = int(min(d, max(0, math.floor(q * (d + 1)))))
            if d > 0 and idx % 7 == 0:
                k = int(rng.integers(0, d + 1))
            if k:
                exps[idx, :] = rng.multinomial(k, probs).astype(exps.dtype)
            degrees[idx] = k
            idx += 1
        log_scales = np.empty(m, dtype=np.float64)
        log_m, log_deg, log_n = math.log(max(1, m)), math.log(max(1, d + 1)), math.log(max(1, n))
        for i, k in enumerate(degrees):
            lc = math.lgamma(d + 1) - math.lgamma(int(k) + 1) - math.lgamma(d - int(k) + 1)
            log_scales[i] = 0.5 * (lc + log_deg + int(k) * log_n - log_m)
        coeff = (rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m))) / math.sqrt(2.0)
        if equation_normalize:
            row = np.linalg.norm(coeff, axis=1)
            coeff = coeff / np.where(row > 0, row, 1.0)[:, None] * math.sqrt(float(m))
        obj = cls(n, d, seed, exps, exps.astype(np.float64).T, log_scales, coeff.astype(np.complex128), m, bool(projective_normalize), bool(dynamic_normalize), bool(equation_normalize), max(1, int(eval_block)))
        obj._generation_seconds = now() - t0
        return obj

    @property
    def terms_per_poly(self) -> int:
        return exact_kostlan_terms(self.n, self.d)

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return float(getattr(self, "_generation_seconds", 0.0))

    def _eval_block(self, ZZ: Any) -> Any:
        ZZ = np.asarray(ZZ, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        az = np.abs(ZZ)
        log_amp = np.log(np.maximum(az, 1e-300)) @ self.feature_exps_t + self.feature_log_scales[None, :]
        phase_arg = np.angle(ZZ) @ self.feature_exps_t
        if self.projective_normalize:
            log_amp -= (0.5 * float(self.d) * np.log1p(np.sum(az * az, axis=1)))[:, None]
        if self.dynamic_normalize:
            shift = np.max(np.where(np.isfinite(log_amp), log_amp, -np.inf), axis=1)
            log_amp = log_amp - np.where(np.isfinite(shift), shift, 0.0)[:, None]
            log_amp = np.clip(log_amp, -745.0, 0.0)
        else:
            log_amp = np.clip(log_amp, -745.0, 700.0)
        Phi = np.exp(log_amp) * np.exp(1j * phase_arg)
        out = Phi @ self.coeff.T
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        return out

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = self._eval_block(z)[0]
        self.eval_count += 1
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        chunks = [self._eval_block(ZZ[i:i + self.eval_block]) for i in range(0, int(ZZ.shape[0]), self.eval_block)]
        out = np.vstack(chunks) if chunks else np.empty((0, self.n), dtype=np.complex128)
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms, "lazy_features": int(self.lazy_features)}


@dataclasses.dataclass
class GeometryKernelKostlanSystem:
    n: int
    d: int
    seed: int
    anchors: Any
    anchor_conj_t: Any
    anchor_den: Any
    coeff: Any
    geometry_anchors: int
    anchor_scales: list[float]
    dynamic_normalize: bool = True
    self_normalize: bool = True
    equation_normalize: bool = True
    eval_block: int = 128
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, seed_index: int = 0, equation_normalize: bool = True, geometry_anchors: int = 0, geometry_anchor_cap: int = 4096, geometry_anchor_scales: Optional[Sequence[float]] = None, dynamic_normalize: bool = True, self_normalize: bool = True, eval_block: int = 128) -> "GeometryKernelKostlanSystem":
        t0 = now()
        n, d = int(n), int(d)
        m = max(n + 2, int(geometry_anchors) if int(geometry_anchors) > 0 else auto_geometry_anchor_count(n, d, int(geometry_anchor_cap)))
        seed = stable_seed(n, d, seed_index, salt=0x314C0DE)
        rng = np.random.default_rng(seed)
        scales = [float(x) for x in (geometry_anchor_scales or []) if math.isfinite(float(x)) and float(x) > 0] or [0.25, 0.5, 1.0, 2.0, 4.0]
        anchors = np.zeros((m, n), dtype=np.complex128)
        idx = 1
        axis_scale = scales[min(2, len(scales) - 1)]
        for j in range(n):
            if idx >= m:
                break
            anchors[idx, j] = complex(axis_scale)
            idx += 1
        sqrt_n = math.sqrt(max(1, n))
        while idx < m:
            shell = scales[(idx - 1 - n) % len(scales)]
            v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
            v = v / max(1e-300, float(np.linalg.norm(v))) * (shell * sqrt_n)
            anchors[idx, :] = v
            idx += 1
        coeff = (rng.standard_normal((n, m)) + 1j * rng.standard_normal((n, m))) / math.sqrt(2.0 * max(1, m))
        if equation_normalize:
            row = np.linalg.norm(coeff, axis=1)
            coeff = coeff / np.where(row > 0, row, 1.0)[:, None]
        anchor_norm2 = np.sum(np.abs(anchors) ** 2, axis=1)
        obj = cls(n, d, seed, anchors, np.conjugate(anchors).T, np.sqrt(1.0 + anchor_norm2), coeff.astype(np.complex128), m, [float(x) for x in scales], bool(dynamic_normalize), bool(self_normalize), bool(equation_normalize), max(1, int(eval_block)))
        obj._generation_seconds = now() - t0
        return obj

    @property
    def terms_per_poly(self) -> int:
        return exact_kostlan_terms(self.n, self.d)

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return float(getattr(self, "_generation_seconds", 0.0))

    def _kernel_block(self, ZZ: Any) -> Any:
        ZZ = np.asarray(ZZ, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        dot = ZZ @ self.anchor_conj_t
        zn = np.sqrt(1.0 + np.sum(np.abs(ZZ) ** 2, axis=1))
        base = (1.0 + dot) / np.maximum(1e-300, zn[:, None] * self.anchor_den[None, :])
        mag = np.abs(base)
        base = np.where(mag > 1.0, base / np.maximum(mag, 1e-300), base)
        mag = np.minimum(mag, 1.0)
        log_amp = float(self.d) * np.log(np.maximum(mag, 1e-300))
        phase_arg = float(self.d) * np.angle(base)
        if self.dynamic_normalize:
            shift = np.max(np.where(np.isfinite(log_amp), log_amp, -np.inf), axis=1)
            log_amp = log_amp - np.where(np.isfinite(shift), shift, 0.0)[:, None]
        K = np.exp(np.clip(log_amp, -745.0, 0.0)) * np.exp(1j * phase_arg)
        if self.self_normalize:
            row = np.sqrt(np.mean(np.abs(K) ** 2, axis=1))
            K = K / np.maximum(row[:, None], 1e-300)
        return K

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = (self._kernel_block(z) @ self.coeff.T)[0]
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        self.eval_count += 1
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            return self.eval(ZZ)[None, :]
        chunks = []
        for i in range(0, int(ZZ.shape[0]), self.eval_block):
            F = self._kernel_block(ZZ[i:i + self.eval_block]) @ self.coeff.T
            F[~np.isfinite(F)] = complex(1e300, 0.0)
            chunks.append(F)
        out = np.vstack(chunks) if chunks else np.empty((0, self.n), dtype=np.complex128)
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms, "geometry_anchors": int(self.geometry_anchors), "geometry_anchor_scales": list(self.anchor_scales)}


def safe_norms(F: Any) -> Any:
    R = np.linalg.norm(np.asarray(F, dtype=np.complex128), axis=1)
    R = np.asarray(R, dtype=float)
    R[~np.isfinite(R)] = np.inf
    return R


def finite_slope_matrix(system: Any, a: Sequence[complex], b: Sequence[complex]) -> Any:
    t0 = now()
    aa = np.asarray(a, dtype=np.complex128)
    bb = np.asarray(b, dtype=np.complex128)
    cur = aa.copy()
    f_prev = system.eval(cur)
    Q = np.zeros((int(system.n), int(system.n)), dtype=np.complex128)
    for j in range(int(system.n)):
        old = cur[j]
        cur[j] = bb[j]
        f_next = system.eval(cur)
        dz = bb[j] - old
        if abs(dz) > 1e-300:
            Q[:, j] = (f_next - f_prev) / dz
        else:
            h = 1e-6 * max(1.0, abs(old))
            plus, minus = cur.copy(), cur.copy()
            plus[j] = old + h
            minus[j] = old - h
            Q[:, j] = (system.eval(plus) - system.eval(minus)) / (2.0 * h)
        f_prev = f_next
    system.slope_count = int(getattr(system, "slope_count", 0)) + 1
    system.seconds_slope = float(getattr(system, "seconds_slope", 0.0)) + (now() - t0)
    return Q


def make_kostlan_system(args: argparse.Namespace, n: int, d: int) -> Any:
    mode = str(getattr(args, "system_mode", "auto")).strip().lower().replace("_", "-")
    dense_terms = exact_kostlan_terms(n, d)
    if mode == "dense" or (mode == "auto" and dense_terms <= int(args.dense_max_terms)):
        return DenseKostlanSystem.make(n, d, int(args.seed_index), bool(args.equation_normalize))
    if mode in {"auto", "geometry", "geometry-kernel", "kernel", "projective-kernel"}:
        return GeometryKernelKostlanSystem.make(n, d, int(args.seed_index), bool(args.equation_normalize), int(args.geometry_anchors), int(args.geometry_anchor_cap), parse_float_list(args.geometry_anchor_scales, [0.25, 0.5, 1.0, 2.0, 4.0], positive=True), bool(args.geometry_dynamic_normalize), bool(args.geometry_self_normalize), int(args.geometry_eval_block))
    if mode not in {"lazy", "lazy-feature", "feature", "stream"}:
        raise ValueError(f"unknown system mode {mode!r}")
    return LazyFeatureKostlanSystem.make(n, d, int(args.seed_index), bool(args.equation_normalize), int(args.lazy_features), int(args.lazy_feature_cap), bool(args.lazy_projective_normalize), bool(args.lazy_dynamic_normalize), int(args.lazy_eval_block))


# ---------------------------------------------------------------------------
# Exact expression polynomial systems and geometry guides
# ---------------------------------------------------------------------------

class SafeExpression:
    _allowed = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub, ast.UAdd, ast.Load)

    def __init__(self, raw: str) -> None:
        self.raw = str(raw).strip()
        tree = ast.parse(self.raw.replace("^", "**"), mode="eval")
        self._validate(tree)
        self.tree = tree

    def _validate(self, node: ast.AST) -> None:
        if not isinstance(node, self._allowed):
            raise ValueError(f"unsupported expression node: {type(node).__name__}")
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
            if not isinstance(node.right, ast.Constant):
                raise ValueError("polynomial powers must be numeric constants")
            if isinstance(node.right.value, complex) or float(node.right.value) != int(node.right.value):
                raise ValueError("polynomial powers must be integer constants")
        for child in ast.iter_child_nodes(node):
            self._validate(child)

    def eval(self, env: dict[str, Any]) -> Any:
        return self._eval(self.tree.body, env)

    def _eval(self, node: ast.AST, env: dict[str, Any]) -> Any:
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            if node.id not in env:
                raise ValueError(f"unknown variable {node.id!r}")
            return env[node.id]
        if isinstance(node, ast.UnaryOp):
            v = self._eval(node.operand, env)
            return -v if isinstance(node.op, ast.USub) else v
        if isinstance(node, ast.BinOp):
            a, b = self._eval(node.left, env), self._eval(node.right, env)
            if isinstance(node.op, ast.Add):
                return a + b
            if isinstance(node.op, ast.Sub):
                return a - b
            if isinstance(node.op, ast.Mult):
                return a * b
            if isinstance(node.op, ast.Div):
                return a / b
            if isinstance(node.op, ast.Pow):
                return a ** int(b)
        raise ValueError(f"unsupported expression node: {type(node).__name__}")


@dataclasses.dataclass
class ExpressionPolynomialSystem:
    n: int
    d: int
    seed: int
    expressions_raw: list[str]
    expressions: list[SafeExpression]
    variable_names: list[str]
    eval_count: int = 0
    slope_count: int = 0
    seconds_eval: float = 0.0
    seconds_slope: float = 0.0

    @classmethod
    def make(cls, n: int, d: int, raw: str, variable_names: Optional[Sequence[str]] = None, seed_index: int = 0) -> "ExpressionPolynomialSystem":
        parts = [p.strip() for p in str(raw).replace("|", ";").split(";") if p.strip()]
        if not parts:
            raise ValueError("--polys/--poly is required for polynomial systems")
        names = [str(p).strip() for p in (variable_names or []) if str(p).strip()] or (["x"] if int(n) == 1 else [f"x{i + 1}" for i in range(int(n))])
        if len(names) != int(n) or len(parts) != int(n):
            raise ValueError(f"expected {n} variables and {n} polynomials")
        return cls(int(n), int(d), stable_seed(n, d, seed_index, salt=0x315EAC7), parts, [SafeExpression(p) for p in parts], names)

    @property
    def terms_per_poly(self) -> int:
        return exact_kostlan_terms(self.n, self.d)

    @property
    def total_terms(self) -> int:
        return int(self.n * self.terms_per_poly)

    @property
    def bezout(self) -> int:
        return int(self.d ** self.n)

    @property
    def generation_seconds(self) -> float:
        return 0.0

    def _env(self, ZZ: Any) -> dict[str, Any]:
        env: dict[str, Any] = {"I": 1j, "j": 1j}
        for i, name in enumerate(self.variable_names):
            env[name] = ZZ[:, i]
        if self.n == 1:
            env.setdefault("x", ZZ[:, 0])
            env.setdefault("z", ZZ[:, 0])
        for i in range(self.n):
            env.setdefault(f"x{i + 1}", ZZ[:, i])
            env.setdefault(f"z{i + 1}", ZZ[:, i])
        return env

    def eval(self, z: Sequence[complex]) -> Any:
        t0 = now()
        out = self.eval_batch(np.asarray(z, dtype=np.complex128)[None, :])[0]
        self.seconds_eval += now() - t0
        return out

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        ZZ = np.asarray(Z, dtype=np.complex128)
        if ZZ.ndim == 1:
            ZZ = ZZ[None, :]
        env = self._env(ZZ)
        cols = []
        for expr in self.expressions:
            val = expr.eval(env)
            arr = np.asarray(val, dtype=np.complex128)
            if arr.ndim == 0:
                arr = np.full(int(ZZ.shape[0]), complex(arr), dtype=np.complex128)
            cols.append(arr.reshape(int(ZZ.shape[0])))
        out = np.stack(cols, axis=1)
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        self.eval_count += int(ZZ.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"eval_count": int(self.eval_count), "slope_count": int(self.slope_count), "seconds_eval": float(self.seconds_eval), "seconds_slope": float(self.seconds_slope), "terms_per_poly": self.terms_per_poly, "total_terms": self.total_terms, "expressions": list(self.expressions_raw), "variables": list(self.variable_names)}


@dataclasses.dataclass
class FittedGeometryKernelSystem:
    base: Any
    kernel: GeometryKernelKostlanSystem
    weights: Any
    ridge: float
    fit_seconds: float
    eval_count: int = 0
    seconds_eval: float = 0.0

    @classmethod
    def make(cls, base: Any, n: int, d: int, args: argparse.Namespace) -> "FittedGeometryKernelSystem":
        t0 = now()
        kernel = GeometryKernelKostlanSystem.make(n, d, int(args.seed_index), False, int(args.geometry_anchors), int(args.geometry_anchor_cap), parse_float_list(args.geometry_anchor_scales, [0.25, 0.5, 1.0, 2.0, 4.0], positive=True), bool(args.geometry_dynamic_normalize), bool(args.geometry_self_normalize), int(args.geometry_eval_block))
        A = kernel._kernel_block(kernel.anchors)
        F = base.eval_batch(kernel.anchors)
        ridge = max(0.0, float(args.geometry_fit_ridge))
        eye = np.eye(int(A.shape[0]), dtype=np.complex128)
        try:
            W = np.linalg.solve(A + ridge * eye, F)
            method = "solve"
        except Exception:
            W, _, _, _ = np.linalg.lstsq(A + ridge * eye, F, rcond=None)
            method = "lstsq"
        obj = cls(base, kernel, np.asarray(W, dtype=np.complex128), ridge, now() - t0)
        obj._fit_method = method
        return obj

    @property
    def n(self) -> int:
        return int(self.base.n)

    @property
    def d(self) -> int:
        return int(self.base.d)

    @property
    def seed(self) -> int:
        return int(self.kernel.seed)

    @property
    def terms_per_poly(self) -> int:
        return int(self.base.terms_per_poly)

    @property
    def total_terms(self) -> int:
        return int(self.base.total_terms)

    @property
    def bezout(self) -> int:
        return int(self.base.bezout)

    @property
    def generation_seconds(self) -> float:
        return float(self.fit_seconds)

    def eval(self, z: Sequence[complex]) -> Any:
        return self.eval_batch(np.asarray(z, dtype=np.complex128)[None, :])[0]

    def eval_batch(self, Z: Any) -> Any:
        t0 = now()
        out = self.kernel._kernel_block(Z) @ self.weights
        out[~np.isfinite(out)] = complex(1e300, 0.0)
        self.eval_count += int(out.shape[0])
        self.seconds_eval += now() - t0
        return out

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return finite_slope_matrix(self, a, b)

    def stats(self) -> dict[str, Any]:
        return {"geometry_fit_eval_count": int(self.eval_count), "geometry_fit_seconds_eval": float(self.seconds_eval), "geometry_fit_seconds": float(self.fit_seconds), "geometry_fit_ridge": float(self.ridge), "geometry_fit_method": str(getattr(self, "_fit_method", "unknown")), "geometry_anchors": int(self.kernel.geometry_anchors)}


@dataclasses.dataclass
class GeometryWrappedSystem:
    exact: Any
    geometry_system: FittedGeometryKernelSystem
    geometry_mode: str = "kernel"

    @property
    def n(self) -> int:
        return int(self.exact.n)

    @property
    def d(self) -> int:
        return int(self.exact.d)

    @property
    def seed(self) -> int:
        return int(self.exact.seed)

    @property
    def terms_per_poly(self) -> int:
        return int(self.exact.terms_per_poly)

    @property
    def total_terms(self) -> int:
        return int(self.exact.total_terms)

    @property
    def bezout(self) -> int:
        return int(self.exact.bezout)

    @property
    def generation_seconds(self) -> float:
        return float(self.exact.generation_seconds) + float(self.geometry_system.fit_seconds)

    def eval(self, z: Sequence[complex]) -> Any:
        return self.exact.eval(z)

    def eval_batch(self, Z: Any) -> Any:
        return self.exact.eval_batch(Z)

    def residuals_batch(self, Z: Any) -> Any:
        return safe_norms(self.eval_batch(Z))

    def slope_matrix(self, a: Sequence[complex], b: Sequence[complex]) -> Any:
        return self.exact.slope_matrix(a, b)

    def stats(self) -> dict[str, Any]:
        out = dict(self.exact.stats())
        out.update({"geometry_mode": self.geometry_mode, **self.geometry_system.stats()})
        return out


def make_system_316(args: argparse.Namespace, n: int, d: int) -> tuple[Any, str, str]:
    source = str(args.system_source).strip().lower().replace("_", "-")
    geom = str(args.geometry_mode).strip().lower().replace("_", "-")
    if source in {"poly", "polynomial", "expr", "expression"}:
        vars_raw = str(args.variables or "").strip()
        variables = [p.strip() for p in vars_raw.replace(";", ",").split(",") if p.strip()] if vars_raw else None
        exact = ExpressionPolynomialSystem.make(n, d, str(args.polys or ""), variables, int(args.seed_index))
        mode = str(args.system_mode).strip().lower().replace("_", "-")
        if geom in {"", "none"} and mode in {"geometry", "geometry-kernel", "kernel", "projective-kernel"}:
            geom = "kernel"
        if geom in {"kernel", "geometry", "geometry-kernel", "projective-kernel"}:
            guide = FittedGeometryKernelSystem.make(exact, n, d, args)
            return GeometryWrappedSystem(exact, guide, "kernel"), "polynomial", "exact+geometry-kernel"
        return exact, "polynomial", "exact"
    if source not in {"ks", "kostlan"}:
        raise ValueError(f"unknown --system-source {source!r}")
    system = make_kostlan_system(args, n, d)
    if isinstance(system, GeometryKernelKostlanSystem):
        return system, "ks", "geometry-kernel"
    if isinstance(system, LazyFeatureKostlanSystem):
        return system, "ks", "lazy-feature"
    return system, "ks", "dense"


# ---------------------------------------------------------------------------
# Charts, starts, and correctors
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class LinearChart:
    A: Any
    Ainv: Any

    @classmethod
    def identity(cls, n: int, scale: float = 1.0) -> "LinearChart":
        A = np.eye(int(n), dtype=np.complex128) * complex(scale)
        Ainv = np.eye(int(n), dtype=np.complex128) / complex(scale)
        return cls(A, Ainv)

    def z_from_y(self, y: Sequence[complex]) -> Any:
        return self.A @ np.asarray(y, dtype=np.complex128)

    def y_from_z(self, z: Sequence[complex]) -> Any:
        return self.Ainv @ np.asarray(z, dtype=np.complex128)


@dataclasses.dataclass
class TargetTrack:
    system: Any
    chart: LinearChart

    def eval(self, y: Sequence[complex]) -> Any:
        return self.system.eval(self.chart.z_from_y(y))

    def eval_batch(self, Y: Any) -> Any:
        YY = np.asarray(Y, dtype=np.complex128)
        if YY.ndim == 1:
            return self.eval(YY)[None, :]
        return self.system.eval_batch(YY @ np.asarray(self.chart.A, dtype=np.complex128).T)

    def residual(self, y: Sequence[complex]) -> float:
        return finite_norm(self.eval(y))

    def residuals_batch(self, Y: Any) -> Any:
        try:
            return safe_norms(self.eval_batch(Y))
        except Exception:
            YY = np.asarray(Y, dtype=np.complex128)
            return np.full(1 if YY.ndim == 1 else int(YY.shape[0]), np.inf, dtype=float)


DEFAULT_POWERS = [2.0 ** k for k in range(-20, 21)] + [3.0 * (2.0 ** k) for k in range(-12, 18)] + [10.0 ** k for k in range(-6, 7)]
DEFAULT_ANGLES_DEG = [0, 6, 12, 18, 24, 32, 40, 48, 56, 64, 72, 80, 86, 89, 90, 91, 94, 100, 108, 116, 128, 140, 152, 164, 172]
DEFAULT_RADII = [0.025, 0.04, 0.06, 0.08, 0.12, 0.18, 0.27, 0.4, 0.6, 0.85, 1.15, 1.55, 2.05, 2.75, 3.6, 4.8, 6.4]
DEFAULT_GAINS = [0.035, 0.055, 0.085, 0.12, 0.18, 0.27, 0.4, 0.58, 0.78, 1.0, 1.28, 1.65, 2.2, 3.0, 4.2, 6.0, 8.5, 12.0]


def raw_direction(n: int, trial: int, seed: int) -> Any:
    vals = []
    for j in range(int(n)):
        h1 = splitmix64(seed + 0xD1A5E + 1000003 * trial + 4099 * (j + 1))
        h2 = splitmix64(seed + 0xBADC0DE + 1000033 * trial + 9176 * (j + 1))
        vals.append(math.exp(0.45 * (2.0 * u01(h2) - 1.0)) * phase(2.0 * math.pi * u01(h1)))
    v = np.asarray(vals, dtype=np.complex128)
    return v / max(1e-300, float(np.linalg.norm(v))) * math.sqrt(max(1, int(n)))


def universal_atlas_start(target: TargetTrack, n: int, trial: int, seed: int, powers: Sequence[float], angles: Sequence[float], radii: Sequence[float], cap: float, roots_found: int, duplicates: int, failures: int, target_count: int, universal_cells: int = 16, universal_shells: int = 5, **_: Any) -> tuple[Any, dict[str, Any]]:
    cells = max(1, int(universal_cells))
    shells = max(1, int(universal_shells))
    candidates = []
    metas = []
    pows = [min(max(float(x), 1e-300), float(cap)) for x in powers if float(x) > 0] or [1.0]
    for k in range(cells):
        idx = trial * cells + k
        layer = (idx // cells) % shells
        pwr = pows[(idx * 37 + 11 * layer + roots_found) % len(pows)]
        rad = radii[(idx * 13 + failures) % len(radii)] if radii else 1.0
        theta = angles[(idx * 19 + duplicates) % len(angles)] if angles else 0.0
        thrust = [1.0, 1.6, 2.5, 4.0, 6.5, 10.0, 16.0][idx % 7]
        lam = min(float(cap), pwr * (thrust ** (0.2 + 0.8 * min(1.0, roots_found / max(1.0, float(target_count))))))
        d = raw_direction(n, idx, seed + 65537 * k)
        y = rad * d
        c, s = math.cos(theta), math.sin(theta)
        pole = phase(0.37 * (idx + 1))
        w = y / pole
        den = -s * w + c
        den = np.where(np.abs(den) < 1e-12, den + 1e-12, den)
        y = lam * pole * ((c * w + s) / den)
        candidates.append(np.asarray(y, dtype=np.complex128))
        metas.append({"homothety": float(lam), "base_homothety": float(pwr), "thales_thrust": float(thrust), "theta_deg": float(math.degrees(theta)), "base_radius": float(rad), "atlas_selected_index": int(idx), "atlas_selected_layer": int(layer)})
    Y = np.asarray(candidates, dtype=np.complex128)
    R = target.residuals_batch(Y)
    idx_best = int(np.nanargmin(R)) if np.any(np.isfinite(R)) else 0
    meta = dict(metas[idx_best])
    meta.update({"chart": "316-standalone-universal-mobius-kernel-atlas", "atlas_mode": "316-compact-batched-score", "atlas_cells_tested": int(cells), "atlas_admissible_cells": int(np.sum(np.isfinite(R))), "atlas_cell_residual": float(R[idx_best]) if idx_best < len(R) else float("inf"), "dup_pressure": float((duplicates + 1.0) / (roots_found + 1.0)), "fail_pressure": float((failures + 1.0) / (trial + 1.0)), "progress": float(min(1.0, roots_found / max(1.0, float(target_count))))})
    return Y[idx_best].copy(), meta


def origin_affine_start(target: TargetTrack, n: int, h: float, max_norm: float) -> tuple[Optional[Any], dict[str, Any]]:
    h0 = max(1e-12, float(h))
    y0 = np.zeros(int(n), dtype=np.complex128)
    try:
        f0 = target.eval(y0)
        P = np.eye(int(n), dtype=np.complex128) * h0
        J = ((target.eval_batch(P) - f0[None, :]) / h0).T
        try:
            y = np.linalg.solve(J, -f0)
            method = "solve"
        except Exception:
            y, _, _, _ = np.linalg.lstsq(J, -f0, rcond=None)
            method = "lstsq"
        if not np.all(np.isfinite(y)):
            return None, {"origin_seed_enabled": True, "origin_seed_status": "nonfinite"}
        cap = float(max_norm)
        yn = float(np.linalg.norm(y))
        if math.isfinite(cap) and cap > 0 and yn > cap:
            y = y * (cap / max(yn, 1e-300))
        return y, {"origin_seed_enabled": True, "origin_seed_status": "ok", "origin_seed_method": method, "origin_seed_h": float(h0), "origin_seed_r0": finite_norm(f0), "origin_seed_r1": target.residual(y), "origin_seed_norm": float(np.linalg.norm(y))}
    except Exception as exc:
        return None, {"origin_seed_enabled": True, "origin_seed_status": f"error:{type(exc).__name__}", "origin_seed_h": float(h0)}


def startopt(target: TargetTrack, y0: Any, trial: int, seed: int, steps: int, candidates: int, gains: Sequence[float], micro_epochs: int) -> tuple[Any, dict[str, Any]]:
    best = np.asarray(y0, dtype=np.complex128).copy()
    best_r = target.residual(best)
    initial = best_r
    evals = 1
    chosen_gain = 1.0
    for step in range(max(0, int(steps))):
        pool = [(1.0, best)]
        for c in range(max(0, int(candidates) - 1)):
            gain = float(gains[(trial + 3 * step + c) % len(gains)])
            cand = gain * best
            if c % 3 == 1:
                noise = raw_direction(len(best), trial + 31 * step + c, seed)
                cand = cand + 0.08 * max(1.0, float(np.linalg.norm(best))) * noise
            elif c % 3 == 2:
                cand = cand * phase(0.11 * (c + 1))
            pool.append((gain, np.asarray(cand, dtype=np.complex128)))
        Y = np.asarray([p[1] for p in pool], dtype=np.complex128)
        R = target.residuals_batch(Y)
        evals += len(pool)
        idx = int(np.nanargmin(R)) if np.any(np.isfinite(R)) else 0
        if float(R[idx]) < best_r:
            best, best_r, chosen_gain = Y[idx].copy(), float(R[idx]), float(pool[idx][0])
    return best, {"startopt_enabled": bool(steps > 0), "startopt_r0": float(initial), "startopt_r1": float(best_r), "startopt_ratio": (float(best_r / initial) if math.isfinite(best_r) and math.isfinite(initial) and initial > 0 else None), "startopt_steps": int(max(0, steps)), "startopt_evals": int(evals), "startopt_micro_epochs": int(max(0, micro_epochs)), "startopt_gain": float(chosen_gain), "startopt_batch_numpy": True}


def _line_lambdas(line_search: int, line_grid: Sequence[float]) -> list[float]:
    vals = [float(x) for x in line_grid if math.isfinite(float(x)) and float(x) > 0]
    if vals:
        return vals
    return [1.0, 0.75, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.03125, 0.02][:max(1, int(line_search))]


def hypercube_delta(target: TargetTrack, y: Any, f: Any, ep: int, seed: int, cloud_nodes: int = 0) -> tuple[Any, dict[str, Any]]:
    n = len(y)
    M = max(int(cloud_nodes), 2 * n + 4, 16)
    yn = max(1.0, float(np.linalg.norm(y)))
    h = 1e-5 * yn
    rng = np.random.default_rng(int(seed) + int(ep) * 1337)
    signs = rng.choice([-1.0, 1.0], size=(M, n))
    dY = h * signs
    dF = target.eval_batch(y[None, :] + dY) - f[None, :]
    A, _, _, _ = np.linalg.lstsq(dY, dF, rcond=None)
    A = A.T
    if not np.all(np.isfinite(A)):
        return np.zeros_like(y), {"hypercube_error": "nonfinite_jacobian", "hypercube_order": 0, "hypercube_nodes": M}
    try:
        d1 = np.linalg.solve(A, -f)
    except Exception:
        return np.zeros_like(y), {"hypercube_error": "singular_jacobian", "hypercube_order": 0, "hypercube_nodes": M}
    nd1 = float(np.linalg.norm(d1))
    if not math.isfinite(nd1) or nd1 < 1e-14 or nd1 > 1e4 * yn:
        return d1, {"hypercube_order": 1, "hypercube_delta1_norm": nd1, "hypercube_nodes": M}
    hh = max(1e-5, min(1e-2, 0.05 * yn / nd1))
    p1, m1 = target.eval(y + hh * d1), target.eval(y - hh * d1)
    p2, m2 = target.eval(y + 2.0 * hh * d1), target.eval(y - 2.0 * hh * d1)
    D2 = (-p2 + 16.0 * p1 - 30.0 * f + 16.0 * m1 - m2) / (12.0 * hh**2)
    D3 = (p2 - 2.0 * p1 + 2.0 * m1 - m2) / (2.0 * hh**3)
    try:
        d2 = np.linalg.solve(A, -0.5 * D2)
    except Exception:
        return d1, {"hypercube_order": 1, "hypercube_error": "singular_delta2", "hypercube_nodes": M, "hypercube_delta1_norm": nd1}
    nd2 = float(np.linalg.norm(d2))
    if not math.isfinite(nd2) or nd2 > 2.0 * max(nd1, yn):
        return d1, {"hypercube_order": 1, "hypercube_error": "delta2_rejected", "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2}
    try:
        f2 = target.eval(y + hh * d2)
        f12 = target.eval(y + hh * d1 + hh * d2)
        cross = (f12 - p1 - f2 + f) / (hh**2)
        d3 = np.linalg.solve(A, -(cross + (1.0 / 6.0) * D3))
    except Exception:
        return d1 + d2, {"hypercube_order": 2, "hypercube_error": "singular_delta3", "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2}
    nd3 = float(np.linalg.norm(d3))
    if not math.isfinite(nd3) or nd3 > 2.0 * max(nd1, yn):
        return d1 + d2, {"hypercube_order": 2, "hypercube_error": "delta3_rejected", "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2, "hypercube_delta3_norm": nd3}
    return d1 + d2 + d3, {"hypercube_order": 4, "hypercube_nodes": M, "hypercube_delta1_norm": nd1, "hypercube_delta2_norm": nd2, "hypercube_delta3_norm": nd3}


def hypercube_inversejet_corrector(target: TargetTrack, y0: Sequence[complex], max_epochs: int, tol: float, accept: float, trial_timeout: float, line_search: int = 12, line_grid: Sequence[float] = (), direction_seed: int = 0, cloud_nodes: int = 0) -> dict[str, Any]:
    t0 = now()
    deadline = t0 + float(trial_timeout) if trial_timeout and trial_timeout > 0 else None
    y = np.asarray(y0, dtype=np.complex128).copy()
    best_y, best_r = y.copy(), target.residual(y)
    status, ok, epochs = "started", False, 0
    total_line, total_hyper, used, last = 0, 0, 0, {}
    lambdas = _line_lambdas(line_search, line_grid)
    for ep in range(max(1, int(max_epochs))):
        if deadline is not None and now() > deadline:
            status = "timeout"
            break
        f = target.eval(y)
        r = finite_norm(f)
        if r < best_r:
            best_y, best_r = y.copy(), r
        if r <= max(float(tol), float(accept)) and (accept <= 0 or r < accept):
            ok, status = True, "converged"
            break
        delta, meta = hypercube_delta(target, y, f, ep, int(direction_seed), int(cloud_nodes))
        last = dict(meta)
        total_hyper += int(meta.get("hypercube_nodes", max(2 * len(y) + 4, 16))) + 7
        used += 1
        if not np.all(np.isfinite(delta)):
            status = "nonfinite-hypercube-step"
            break
        L = np.asarray(lambdas, dtype=float)
        Y = y[None, :] + L[:, None] * delta[None, :]
        R = target.residuals_batch(Y)
        total_line += int(len(Y))
        better = np.isfinite(R) & ((R < r) | (R < best_r))
        if not np.any(better):
            status = "no-hypercube-decrease"
            break
        idx = int(np.nanargmin(np.where(better, R, np.inf)))
        y = Y[idx].copy()
        if float(R[idx]) < best_r:
            best_y, best_r = y.copy(), float(R[idx])
        epochs = ep + 1
        last["hypercube_chosen_lambda"] = float(L[idx])
    else:
        status = "max-epochs"
    final_r = target.residual(best_y)
    if final_r <= max(float(tol), float(accept)) and (accept <= 0 or final_r < accept):
        ok, status, best_r = True, "converged", final_r
    return {"accepted": bool(ok if accept <= 0 else (math.isfinite(best_r) and best_r < accept)), "ok": bool(ok), "status": status, "epochs": int(epochs), "residual": float(best_r), "y": best_y, "seconds": float(now() - t0), "slope_cond": None, "corrector": "316-306-hypercube-inversejet", "line_search_evals": int(total_line), "line_lambdas": [float(x) for x in lambdas], "hypercube_total_evals": int(total_hyper), "hypercube_used_count": int(used), "halley_enabled": False, "halley_total_evals": 0, "halley_used_count": 0, **last}


class NonlinearChartTarget:
    def __init__(self, base: TargetTrack, scale: complex, reciprocal: bool) -> None:
        self.base, self.scale, self.reciprocal = base, complex(scale), bool(reciprocal)

    def to_base(self, u: Any) -> Any:
        uu = np.asarray(u, dtype=np.complex128)
        if not self.reciprocal:
            return self.scale * uu
        den = self.scale * uu
        den = np.where(np.abs(den) < 1e-14, den + 1e-14, den)
        return 1.0 / den

    def from_base(self, y: Any) -> Any:
        yy = np.asarray(y, dtype=np.complex128)
        if not self.reciprocal:
            return yy / self.scale
        den = self.scale * yy
        den = np.where(np.abs(den) < 1e-14, den + 1e-14, den)
        return 1.0 / den

    def eval(self, u: Sequence[complex]) -> Any:
        return self.base.eval(self.to_base(u))

    def eval_batch(self, U: Any) -> Any:
        UU = np.asarray(U, dtype=np.complex128)
        return self.base.eval_batch(self.to_base(UU))

    def residual(self, u: Sequence[complex]) -> float:
        return finite_norm(self.eval(u))

    def residuals_batch(self, U: Any) -> Any:
        return safe_norms(self.eval_batch(U))


def complex_scale_palette(gains: Sequence[float], phases: Sequence[float], top: int) -> list[complex]:
    out: list[complex] = []
    for g in gains:
        for ph in phases:
            out.append(float(g) * phase(float(ph)))
            if len(out) >= int(top):
                return out
    return out or [1.0 + 0.0j]


def chart_candidates(base: TargetTrack, y: Any, scales: Sequence[complex], inversion: bool, top: int) -> list[tuple[float, NonlinearChartTarget, Any, dict[str, Any]]]:
    recs = []
    for s in list(scales)[:max(1, int(top))]:
        for recip in ([False, True] if inversion else [False]):
            ct = NonlinearChartTarget(base, complex(s), recip)
            u0 = ct.from_base(y)
            r = ct.residual(u0)
            e = log_energy(u0)
            if math.isfinite(r) and math.isfinite(e):
                recs.append((r * (1.0 + 0.002 * e), ct, u0, {"kind": "reciprocal" if recip else "direct", "scale": s, "score": r, "log_energy": e}))
    recs.sort(key=lambda x: x[0])
    return recs


def lazy_irp_hypercube_inversejet_corrector_312(base: TargetTrack, y0: Sequence[complex], max_epochs: int, tol: float, accept: float, trial_timeout: float, line_search: int = 12, line_grid: Sequence[float] = (), direction_seed: int = 0, cloud_nodes: int = 0, irp_layers: int = 2, irp_inner_epochs: int = 2, irp_scales: Optional[Sequence[complex]] = None, irp_chart_top: int = 2, irp_inversion: bool = True, collapse: bool = True, collapse_residual: float = 1e-4, collapse_drop: float = 0.42, collapse_rel_step: float = 0.35, collapse_after: int = 2, local_inner_epochs: int = 3, lazy_direct_epochs: int = 1, lazy_trigger_drop: float = 0.82, lazy_trigger_after: int = 1, lazy_bad_cond: float = 1e10, lazy_log_energy: float = 8.0, eager_irp: bool = False, rescue_collapsed: bool = False) -> dict[str, Any]:
    del lazy_bad_cond
    t0 = now()
    deadline = t0 + float(trial_timeout) if trial_timeout and trial_timeout > 0 else None
    y = np.asarray(y0, dtype=np.complex128).copy()
    best_y, best_r = y.copy(), base.residual(y)
    status, ok, epochs_done = "started", False, 0
    scales = list(irp_scales or [1.0 + 0.0j])
    total_line = total_hyper = hyper_used = triggers = rescues = direct_steps = skipped = chart_switches = recip_uses = direct_uses = 0
    collapsed, locality_hits, stagnation_hits = False, 0, 0
    collapse_epoch = collapse_reason = None
    last: dict[str, Any] = {}
    last_chart: dict[str, Any] = {"kind": "direct", "scale": 1.0 + 0j, "score": best_r, "log_energy": log_energy(y)}

    def absorb(loc: dict[str, Any]) -> None:
        nonlocal total_line, total_hyper, hyper_used
        total_line += int(loc.get("line_search_evals", 0) or 0)
        total_hyper += int(loc.get("hypercube_total_evals", 0) or 0)
        hyper_used += int(loc.get("hypercube_used_count", 0) or 0)

    for ep in range(max(1, int(max_epochs))):
        if deadline is not None and now() > deadline:
            status = "timeout"
            break
        r0 = base.residual(y)
        if r0 < best_r:
            best_y, best_r = y.copy(), r0
        if r0 <= max(float(tol), float(accept)) and (accept <= 0 or r0 < accept):
            ok, status = True, "converged"
            break
        loc = hypercube_inversejet_corrector(base, y, max(1, int(local_inner_epochs if collapsed else lazy_direct_epochs)), tol, accept, 0.0 if deadline is None else max(0.0, deadline - now()), line_search, line_grid, int(direction_seed) + 1000003 * ep + 17, cloud_nodes)
        absorb(loc)
        last = dict(loc)
        yd = np.asarray(loc.get("y", y), dtype=np.complex128)
        rd = base.residual(yd)
        improved = math.isfinite(rd) and rd <= r0 * (1.0 - 1e-14)
        if improved:
            direct_steps += 1
            stagnation_hits = 0
            if rd < best_r:
                best_y, best_r = yd.copy(), rd
        else:
            stagnation_hits += 1
        strong = bool(improved and rd / max(r0, 1e-300) <= float(lazy_trigger_drop) and log_energy(yd) <= float(lazy_log_energy))
        trigger = bool(eager_irp and not collapsed) or (not improved and stagnation_hits >= max(1, int(lazy_trigger_after))) or (improved and not strong)
        if collapsed and not rescue_collapsed:
            trigger = False
        if not trigger:
            if improved:
                step = float(np.linalg.norm(yd - y)) / max(1.0, float(np.linalg.norm(y)))
                drop = rd / max(r0, 1e-300)
                y = yd.copy()
                direct_uses += 1
                epochs_done = ep + 1
                if collapse and not collapsed:
                    if rd <= float(collapse_residual):
                        collapsed, collapse_epoch, collapse_reason = True, ep, "residual-threshold"
                    elif drop <= float(collapse_drop) and step <= float(collapse_rel_step):
                        locality_hits += 1
                        if locality_hits >= max(1, int(collapse_after)):
                            collapsed, collapse_epoch, collapse_reason = True, ep, "observed-local-contraction"
                status = "local-collapsed" if collapsed else "direct-306-running"
                continue
            status = "no-lazy-direct-decrease"
            break
        triggers += 1
        rescues += 1 if not improved else 0
        base_y, base_r = (yd.copy(), rd) if improved else (y.copy(), r0)
        irp_y, irp_r = base_y.copy(), base_r
        any_irp = False
        for layer in range(max(1, int(irp_layers))):
            cands = chart_candidates(base, irp_y, scales, bool(irp_inversion), int(irp_chart_top))
            if not cands:
                status = "no-admissible-lazy-irp-chart"
                break
            best_cy = None
            best_cr = irp_r
            best_cloc: dict[str, Any] = {}
            best_chart: dict[str, Any] = {}
            for _, ct, u0, cm in cands:
                loc2 = hypercube_inversejet_corrector(ct, u0, max(1, int(irp_inner_epochs)), tol, accept, 0.0, line_search, line_grid, int(direction_seed) + 1000003 * ep + 65537 * layer, cloud_nodes)
                absorb(loc2)
                yb = ct.to_base(loc2.get("y", u0))
                rb = base.residual(yb)
                if math.isfinite(rb) and rb < best_cr:
                    best_cy, best_cr, best_cloc, best_chart = yb.copy(), rb, dict(loc2), dict(cm)
            if best_cy is None:
                status = "lazy-irp-no-extra-decrease"
                break
            irp_y, irp_r, last, last_chart = best_cy, best_cr, best_cloc, best_chart
            any_irp = True
            chart_switches += 1
            recip_uses += 1 if best_chart.get("kind") == "reciprocal" else 0
            direct_uses += 0 if best_chart.get("kind") == "reciprocal" else 1
            if irp_r < best_r:
                best_y, best_r = irp_y.copy(), irp_r
            if irp_r <= max(float(tol), float(accept)) and (accept <= 0 or irp_r < accept):
                ok, status = True, "converged"
                break
        if ok:
            break
        if any_irp and irp_r < r0:
            y = irp_y.copy()
            epochs_done = ep + 1
            status = "lazy-irp-running"
            continue
        if improved:
            y = yd.copy()
            epochs_done = ep + 1
            skipped += 1
            status = "direct-after-weak-irp"
            continue
        status = status if status != "started" else "no-lazy-irp-or-direct-decrease"
        break
    final_r = base.residual(best_y)
    if final_r <= max(float(tol), float(accept)) and (accept <= 0 or final_r < accept):
        ok, status, best_r = True, "converged", final_r
    return {"accepted": bool(ok if accept <= 0 else (math.isfinite(best_r) and best_r < accept)), "ok": bool(ok), "status": status, "epochs": int(epochs_done), "residual": float(best_r), "y": best_y, "seconds": float(now() - t0), "slope_cond": None, "corrector": "316-312-lazy-irp-plus-306-hypercube-inversejet", "line_search_evals": int(total_line), "hypercube_total_evals": int(total_hyper), "hypercube_used_count": int(hyper_used), "irp_layers_completed": int(chart_switches), "irp_chart_switches": int(chart_switches), "irp_collapsed": bool(collapsed), "irp_collapse_epoch": collapse_epoch, "irp_collapse_reason": collapse_reason, "irp_reciprocal_uses": int(recip_uses), "irp_direct_uses": int(direct_uses), "irp_lazy_triggers": int(triggers), "irp_lazy_rescues": int(rescues), "irp_lazy_direct_steps": int(direct_steps), "irp_lazy_direct_good": int(direct_steps - triggers if direct_steps >= triggers else 0), "irp_lazy_direct_weak": int(triggers), "irp_lazy_skipped": int(skipped), "irp_last_chart_kind": last_chart.get("kind"), "irp_last_chart_scale": cjson(last_chart.get("scale", 1.0 + 0j)), "irp_last_chart_log_energy": last_chart.get("log_energy"), "halley_enabled": False, "halley_total_evals": 0, "halley_used_count": 0, **{k: v for k, v in last.items() if k != "y"}}


def polish_exact_root(system: Any, chart: LinearChart, y: Any, accept: float, max_steps: int) -> tuple[Any, Any, float, dict[str, Any]]:
    z = chart.z_from_y(y)
    r = finite_norm(system.eval(z))
    meta = {"exact_polish_enabled": bool(max_steps > 0), "exact_polish_steps": 0, "exact_polish_r0": float(r), "exact_polish_r1": float(r)}
    if max_steps <= 0 or not math.isfinite(r):
        return z, np.asarray(y, dtype=np.complex128), r, meta
    best_z, best_f, best_r = np.asarray(z, dtype=np.complex128).copy(), system.eval(z), r
    for step in range(int(max_steps)):
        if best_r <= float(accept):
            break
        try:
            J = system.slope_matrix(best_z, best_z)
            delta, _, _, _ = np.linalg.lstsq(np.asarray(J, dtype=np.complex128), -best_f, rcond=None)
        except Exception:
            meta["exact_polish_status"] = "linear-solve-failed"
            break
        improved = False
        for gain in (1.0, 0.5, 0.25, 0.125, 0.0625):
            cand_z = best_z + gain * delta
            cand_f = system.eval(cand_z)
            cand_r = finite_norm(cand_f)
            if cand_r < best_r:
                best_z, best_f, best_r = cand_z, cand_f, cand_r
                meta["exact_polish_steps"] = step + 1
                improved = True
                break
        if not improved:
            meta["exact_polish_status"] = "no-decrease"
            break
    meta["exact_polish_r1"] = float(best_r)
    return best_z, chart.y_from_z(best_z), best_r, meta


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_case(args: argparse.Namespace, case_raw: str) -> dict[str, Any]:
    t_case = now()
    n, d = parse_case(case_raw)
    system, source, backend = make_system_316(args, n, d)
    chart = LinearChart.identity(n, float(args.linear_scale))
    target = TargetTrack(system, chart)
    guide_system = getattr(system, "geometry_system", None)
    guide_target = TargetTrack(guide_system, chart) if bool(args.geometry_starts) and guide_system is not None else target
    starts = parse_start_points(args.starts, n)
    powers = sorted(set(round(float(x), 16) for x in parse_float_list(args.powers, DEFAULT_POWERS, positive=True)))
    powers = [min(max(x, 1e-300), float(args.power_cap)) for x in powers]
    angles = [math.radians(x) for x in parse_float_list(args.angles, DEFAULT_ANGLES_DEG)]
    radii = parse_float_list(args.rays, DEFAULT_RADII, positive=True)
    gains = parse_float_list(args.startopt_gains, DEFAULT_GAINS, positive=True)
    roots: list[dict[str, Any]] = []
    trials: list[dict[str, Any]] = []
    failures = duplicates = 0
    t_extract = now()
    for trial in range(int(args.pool)):
        if len(roots) >= int(args.count):
            break
        explicit = trial < len(starts)
        if explicit:
            y_raw = np.asarray(starts[trial], dtype=np.complex128)
            geom = {"chart": "316-explicit-start", "atlas_mode": "explicit-start", "homothety": 1.0, "theta_deg": None}
            start_target = target
        else:
            start_target = guide_target
            y_raw = None
            geom: dict[str, Any] = {}
            if bool(args.origin_seed) and backend in {"lazy-feature", "geometry-kernel"} and (trial == 0 if int(args.origin_seed_period) <= 0 else trial % max(1, int(args.origin_seed_period)) == 0):
                cap = float(args.origin_seed_max_norm) if float(args.origin_seed_max_norm) > 0 else 2.0 * math.sqrt(max(1, n))
                y_origin, ometa = origin_affine_start(start_target, n, float(args.origin_seed_h), cap)
                if y_origin is not None:
                    y_raw = y_origin
                    geom = {"chart": "316-origin-affine-finite-difference-seed", "atlas_mode": "origin-affine-seed", "homothety": 1.0, "theta_deg": None, **ometa}
            if y_raw is None:
                y_raw, geom = universal_atlas_start(start_target, n, trial, system.seed + 0x113000, powers, angles, radii, float(args.power_cap), len(roots), duplicates, failures, int(args.count), int(args.universal_cells), int(args.universal_shells), cell_probe_radius=float(args.cell_probe_radius), cell_descent_min=float(args.cell_descent_min), cell_equal_gap_min=float(args.cell_equal_gap_min), cell_log_max=float(args.cell_log_max), universal_cycle=bool(args.universal_cycle))
        y0, smeta = startopt(start_target, y_raw, trial, system.seed + 0x112555, int(args.startopt_steps), int(args.startopt_candidates), gains, int(args.startopt_micro_epochs))
        loc = lazy_irp_hypercube_inversejet_corrector_312(target, y0, int(args.epochs), float(args.tol), float(args.accept), float(args.trial_timeout), int(args.line_search), parse_float_list(args.line_grid, []), system.seed + 7919 * trial, int(args.hypercube_nodes), int(args.irp_layers), int(args.irp_inner_epochs), complex_scale_palette(parse_float_list(args.irp_gains, [1.0, 0.5, 2.0, 0.25, 4.0, 0.125, 8.0], positive=True), parse_float_list(args.irp_phases, [0.0, 0.08, -0.08, 0.19, -0.19]), int(args.irp_top)), int(args.irp_chart_top), bool(args.irp_inversion), bool(args.collapse), float(args.collapse_residual), float(args.collapse_drop), float(args.collapse_rel_step), int(args.collapse_after), int(args.local_inner_epochs), int(args.lazy_direct_epochs), float(args.lazy_trigger_drop), int(args.lazy_trigger_after), float(args.lazy_bad_cond), float(args.lazy_log_energy), bool(args.eager_irp), bool(args.rescue_collapsed))
        if source == "polynomial":
            z, y_final, residual, polish_meta = polish_exact_root(system, chart, loc["y"], float(args.accept), int(args.exact_polish_steps))
        else:
            y_final = np.asarray(loc["y"], dtype=np.complex128)
            z = chart.z_from_y(y_final)
            residual = finite_norm(system.eval(z))
            polish_meta = {"exact_polish_enabled": False, "exact_polish_steps": 0, "exact_polish_r0": float(residual), "exact_polish_r1": float(residual)}
        accepted = bool(math.isfinite(residual) and residual < float(args.accept))
        rec = {"trial": int(trial), "accepted": accepted, "status": loc.get("status"), "r1": float(residual), "epochs": int(loc.get("epochs", 0)), "seconds": float(loc.get("seconds", 0.0)), **{k: v for k, v in loc.items() if k != "y"}, **geom, **smeta, **polish_meta}
        if bool(args.verbose_trials):
            rec["z"] = root_to_json(z)
            rec["y0"] = root_to_json(y0)
        if not accepted:
            failures += 1
            trials.append(rec)
            continue
        dup = cluster_index(roots, z, float(args.cluster_sep))
        if dup is not None:
            duplicates += 1
            rec["status"] = "duplicate"
            rec["cluster"] = int(dup)
            trials.append(rec)
            continue
        root = {"id": len(roots), "source": "316-standalone-exact-system-geometry-guided", "trial": int(trial), "z_complex": np.asarray(z, dtype=np.complex128).copy(), "y_complex": np.asarray(y_final, dtype=np.complex128).copy(), "residual": float(residual), "realness": realness(z), "cond": None, "epochs": int(loc.get("epochs", 0)), "seconds": float(loc.get("seconds", 0.0)), **{k: v for k, v in loc.items() if k != "y"}, **geom, **smeta, **polish_meta}
        root["score"] = score_root(float(root["residual"]), float(root["realness"]), root["cond"])
        roots.append(root)
        rec["status"] = "new-root"
        rec["root_id"] = int(root["id"])
        trials.append(rec)
    encoded_roots = []
    for root in sorted(roots, key=lambda q: (float(q.get("score", float("inf"))), int(q.get("id", 0)))):
        rr = dict(root)
        rr["z"] = root_to_json(rr.pop("z_complex"))
        rr["y"] = root_to_json(rr.pop("y_complex"))
        encoded_roots.append(rr)
    geom_kernel = getattr(getattr(system, "geometry_system", None), "kernel", system)
    result = {"script": Path(__file__).name, "autonomous": True, "dependencies": {"python_scripts": [], "numpy": True}, "mode": "316-standalone-exact-system-plus-geometry-guided-lazy-irp-hypercube-inversejet", "flow_formula": "standalone exact/polynomial or KS system -> optional fitted geometry starts -> StartOpt -> 312 lazy IRP / 306 hypercube inverse-jet -> exact residual validation", "case": f"{n},{d}", "family": source, "system_source": source, "system_backend": backend, "seed_index": int(args.seed_index), "seed": int(system.seed), "n": int(n), "degree": int(d), "terms_per_poly": int(system.terms_per_poly), "terms": int(system.total_terms), "bezout": int(system.bezout), "equation_normalize": bool(args.equation_normalize), "linear_A": [[cjson(chart.A[i, j]) for j in range(n)] for i in range(n)], "parameters": {"system_source": str(args.system_source), "polys": str(args.polys or ""), "variables": str(args.variables or ""), "geometry_mode": str(args.geometry_mode), "geometry_starts": bool(args.geometry_starts), "geometry_fit_ridge": float(args.geometry_fit_ridge), "exact_polish_steps": int(args.exact_polish_steps), "starts": str(args.starts or ""), "system_mode": str(args.system_mode), "count": int(args.count), "pool": int(args.pool), "accept": float(args.accept), "tol": float(args.tol), "epochs": int(args.epochs), "cluster_sep": float(args.cluster_sep), "line_search": int(args.line_search), "hypercube_nodes": int(args.hypercube_nodes), "startopt_steps": int(args.startopt_steps), "startopt_candidates": int(args.startopt_candidates), "geometry_anchors": int(getattr(geom_kernel, "geometry_anchors", 0) or 0), "geometry_anchor_cap": int(args.geometry_anchor_cap), "geometry_anchor_scales": parse_float_list(args.geometry_anchor_scales, [0.25, 0.5, 1.0, 2.0, 4.0], positive=True)}, "roots": encoded_roots, "trials": trials if bool(args.verbose_trials) else trials[: min(len(trials), int(args.keep_trials))], "summary": {"requested_roots": int(args.count), "unique_roots": len(roots), "success": bool(len(roots) >= int(args.count)), "trials_used": len(trials), "duplicates": int(duplicates), "failures": int(failures), "generation_seconds": float(system.generation_seconds), "extract_seconds": float(now() - t_extract), "total_seconds": float(now() - t_case), "eval_stats": system.stats()}}
    return result


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="316 standalone NumPy exact-system + geometry-guided lazy IRP / hypercube inverse-jet engine.")
    p.add_argument("--cases", default="2,4")
    p.add_argument("--seed-index", type=int, default=0)
    p.add_argument("--equation-normalize", action="store_true", default=False)
    p.add_argument("--no-equation-normalize", dest="equation_normalize", action="store_false")
    p.add_argument("--system-mode", choices=["auto", "dense", "geometry-kernel", "geometry", "kernel", "projective-kernel", "lazy-feature", "lazy", "feature", "stream"], default="auto")
    p.add_argument("--dense-max-terms", type=int, default=250000)
    p.add_argument("--geometry-anchors", type=int, default=0)
    p.add_argument("--geometry-anchor-cap", type=int, default=4096)
    p.add_argument("--geometry-anchor-scales", default="0.25,0.5,1,2,4")
    p.add_argument("--geometry-dynamic-normalize", action="store_true", default=True)
    p.add_argument("--no-geometry-dynamic-normalize", dest="geometry_dynamic_normalize", action="store_false")
    p.add_argument("--geometry-self-normalize", action="store_true", default=True)
    p.add_argument("--no-geometry-self-normalize", dest="geometry_self_normalize", action="store_false")
    p.add_argument("--geometry-eval-block", type=int, default=128)
    p.add_argument("--lazy-features", type=int, default=0)
    p.add_argument("--lazy-feature-cap", type=int, default=8192)
    p.add_argument("--lazy-projective-normalize", action="store_true", default=True)
    p.add_argument("--no-lazy-projective-normalize", dest="lazy_projective_normalize", action="store_false")
    p.add_argument("--lazy-dynamic-normalize", action="store_true", default=True)
    p.add_argument("--no-lazy-dynamic-normalize", dest="lazy_dynamic_normalize", action="store_false")
    p.add_argument("--lazy-eval-block", type=int, default=128)
    p.add_argument("--origin-seed", action="store_true", default=True)
    p.add_argument("--no-origin-seed", dest="origin_seed", action="store_false")
    p.add_argument("--origin-seed-h", type=float, default=1e-5)
    p.add_argument("--origin-seed-max-norm", type=float, default=0.0)
    p.add_argument("--origin-seed-period", type=int, default=0)
    p.add_argument("--linear-scale", type=float, default=1.0)
    p.add_argument("--count", "--thales-count", "--useful-count", type=int, default=8)
    p.add_argument("--pool", "--thales-pool", "--useful-pool", type=int, default=4096)
    p.add_argument("--epochs", "--thales-epochs", "--useful-epochs", type=int, default=24)
    p.add_argument("--tol", type=float, default=1e-12)
    p.add_argument("--accept", "--residual-accept", type=float, default=1e-8)
    p.add_argument("--cluster-sep", type=float, default=1e-8)
    p.add_argument("--trial-timeout", "--thales-trial-timeout", "--useful-trial-timeout", type=float, default=0.0)
    p.add_argument("--line-search", type=int, default=12)
    p.add_argument("--hypercube-nodes", type=int, default=0)
    p.add_argument("--irp-layers", type=int, default=2)
    p.add_argument("--irp-inner-epochs", type=int, default=2)
    p.add_argument("--local-inner-epochs", type=int, default=3)
    p.add_argument("--lazy-direct-epochs", type=int, default=1)
    p.add_argument("--lazy-trigger-drop", type=float, default=0.82)
    p.add_argument("--lazy-trigger-after", type=int, default=1)
    p.add_argument("--lazy-bad-cond", type=float, default=1e10)
    p.add_argument("--lazy-log-energy", type=float, default=8.0)
    p.add_argument("--eager-irp", action="store_true", default=False)
    p.add_argument("--rescue-collapsed", action="store_true", default=False)
    p.add_argument("--irp-chart-top", type=int, default=2)
    p.add_argument("--irp-gains", default="1,0.5,2,0.25,4,0.125,8")
    p.add_argument("--irp-phases", default="0,0.08,-0.08,0.19,-0.19")
    p.add_argument("--irp-top", type=int, default=14)
    p.add_argument("--irp-inversion", action="store_true", default=True)
    p.add_argument("--no-irp-inversion", dest="irp_inversion", action="store_false")
    p.add_argument("--collapse", action="store_true", default=True)
    p.add_argument("--no-collapse", dest="collapse", action="store_false")
    p.add_argument("--collapse-residual", type=float, default=1e-4)
    p.add_argument("--collapse-drop", type=float, default=0.42)
    p.add_argument("--collapse-rel-step", type=float, default=0.35)
    p.add_argument("--collapse-after", type=int, default=2)
    p.add_argument("--probe-scale", type=float, default=0.035)
    p.add_argument("--probe-candidates", type=int, default=8)
    p.add_argument("--probe-radii", default="0,0.35,0.7,1,1.6,2.6,4.2")
    p.add_argument("--probe-self", action="store_true", default=True)
    p.add_argument("--no-probe-self", dest="probe_self", action="store_false")
    p.add_argument("--probe-axis", action="store_true", default=True)
    p.add_argument("--no-probe-axis", dest="probe_axis", action="store_false")
    p.add_argument("--probe-condition-top", type=int, default=2)
    p.add_argument("--probe-condition-weight", type=float, default=0.0025)
    p.add_argument("--probe-equal-value-weight", type=float, default=0.015)
    p.add_argument("--probe-equal-value-eps", type=float, default=1e-12)
    p.add_argument("--probe-curvature-top", type=int, default=2)
    p.add_argument("--probe-curvature-weight", type=float, default=0.003)
    p.add_argument("--probe-curvature-mid", type=float, default=0.5)
    p.add_argument("--trust-cond-weight", type=float, default=0.06)
    p.add_argument("--trust-cond-min", type=float, default=4.0)
    p.add_argument("--halley", action="store_true", default=True)
    p.add_argument("--no-halley", dest="halley", action="store_false")
    p.add_argument("--halley-gate-residual", type=float, default=0.25)
    p.add_argument("--halley-probe-fraction", type=float, default=0.50)
    p.add_argument("--halley-cond-weight", type=float, default=0.025)
    p.add_argument("--halley-min-gate", type=float, default=0.04)
    p.add_argument("--halley-max-correction", type=float, default=1.25)
    p.add_argument("--probe-log-weight", type=float, default=0.0005)
    p.add_argument("--probe-log-delta-weight", type=float, default=0.0015)
    p.add_argument("--halley-log-weight", type=float, default=0.03)
    p.add_argument("--halley-log-scale", type=float, default=1.0)
    p.add_argument("--tensor-extra-directions", type=int, default=0)
    p.add_argument("--line-grid", default="1,0.75,0.5,0.35,0.25,0.18,0.125,0.09,0.0625,0.045,0.03125,0.02")
    p.add_argument("--powers", "--thales-powers", default=None)
    p.add_argument("--power-cap", "--thales-power-cap", type=float, default=1048576.0)
    p.add_argument("--angles", "--thales-angles", default=None)
    p.add_argument("--rays", "--thales-rays", default=None)
    p.add_argument("--startopt-steps", type=int, default=1)
    p.add_argument("--startopt-candidates", type=int, default=12)
    p.add_argument("--startopt-gains", default=None)
    p.add_argument("--startopt-micro-epochs", type=int, default=0)
    p.add_argument("--universal-cells", type=int, default=16)
    p.add_argument("--universal-shells", type=int, default=5)
    p.add_argument("--cell-probe-radius", type=float, default=0.14)
    p.add_argument("--cell-descent-min", type=float, default=1.02)
    p.add_argument("--cell-equal-gap-min", type=float, default=1e-10)
    p.add_argument("--cell-log-max", type=float, default=80.0)
    p.add_argument("--universal-cycle", action="store_true", default=True)
    p.add_argument("--no-universal-cycle", dest="universal_cycle", action="store_false")
    p.add_argument("--out", "--thales-out", "--useful-out", default=None)
    p.add_argument("--outdir", default="/private/tmp/316_standalone_out")
    p.add_argument("--keep-trials", type=int, default=160)
    p.add_argument("--verbose-trials", action="store_true")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--system-source", choices=["ks", "kostlan", "polynomial", "poly", "expr", "expression"], default="ks")
    p.add_argument("--polys", "--poly", default=None)
    p.add_argument("--variables", default=None)
    p.add_argument("--geometry-mode", choices=["none", "kernel", "geometry", "geometry-kernel", "projective-kernel"], default="none")
    p.add_argument("--geometry-starts", action="store_true", default=True)
    p.add_argument("--no-geometry-starts", dest="geometry_starts", action="store_false")
    p.add_argument("--geometry-fit-ridge", type=float, default=1e-8)
    p.add_argument("--exact-polish-steps", type=int, default=6)
    p.add_argument("--starts", default=None)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if bool(args.self_test):
        args.system_source = "polynomial"
        args.system_mode = "geometry-kernel"
        args.geometry_mode = "kernel"
        args.cases = "1,2"
        args.polys = "x^2 - 3*x - 10"
        args.starts = "-8,4"
        args.count = 2
        args.pool = min(int(args.pool), 32)
        args.epochs = min(int(args.epochs), 16)
        args.accept = min(float(args.accept), 1e-8)
        args.out = args.out or "/private/tmp/316_standalone_out/self_test_316.json"
    cases = [c.strip() for c in str(args.cases).replace("|", ";").split(";") if c.strip()]
    outputs = [run_case(args, c) for c in cases]
    final = outputs[0] if len(outputs) == 1 else {"script": Path(__file__).name, "autonomous": True, "cases": outputs}
    out = Path(args.out) if args.out else Path(args.outdir) / f"316_standalone_{cases[0].replace(',', 'x') if cases else 'case'}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(final, indent=2), encoding="utf-8")
    print("=" * 120, flush=True)
    print("316 STANDALONE EXACT SYSTEM + GEOMETRY-GUIDED LAZY IRP / HYPERCUBE INVERSE-JET", flush=True)
    print("NumPy + stdlib only; no local flow imports.", flush=True)
    print("=" * 120, flush=True)
    for r in outputs:
        s = r["summary"]
        geom = r.get("parameters", {}).get("geometry_anchors", 0)
        suffix = f", geometry_anchors={geom}" if geom else ""
        print(f"case={r.get('family')}({r['n']},{r['degree']}), backend={r.get('system_backend')}{suffix}, seed={r['seed']}", flush=True)
        print(f"roots={s['unique_roots']}/{s['requested_roots']} success={s['success']} trials={s['trials_used']} duplicates={s['duplicates']} failures={s['failures']}", flush=True)
        print(f"seconds: generation={s['generation_seconds']:.2f}, extract={s['extract_seconds']:.2f}, total={s['total_seconds']:.2f}", flush=True)
        if r.get("roots"):
            best = r["roots"][0]
            print(f"best_root: residual={float(best.get('residual', float('inf'))):.3e}, trial={best.get('trial')}", flush=True)
    print(f"out={out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
