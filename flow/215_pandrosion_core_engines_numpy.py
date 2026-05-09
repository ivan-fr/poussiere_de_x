#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
215 Pandrosion Core Engines - standalone NumPy

This file contains only generic Pandrosion numerical engines.
No trading logic.
No backtest.
No data download.
No pandas/SciPy/PyTorch/JAX/sklearn.

Core engines:
  1. raw finite-slope Pandrosion;
  2. log-stable optimal-probe Pandrosion;
  3. equal-value pole avoidance;
  4. singularity-aware finite-slope probe scoring;
  5. gated finite-slope Halley correction;
  6. spectral tensor correction filtering;
  7. optional entropic probe blending.

The central operation is always:
    Q_G(a,b) delta = -G(a)

where Q_G(a,b) is a Pandrosion finite-slope matrix built from probes.
No analytic Jacobian or Hessian is used.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Optional, Sequence

import numpy as np


# ---------------------------------------------------------------------------
# Basic utilities
# ---------------------------------------------------------------------------

def now() -> float:
    return time.time()


def phase(theta: float) -> complex:
    return complex(math.cos(float(theta)), math.sin(float(theta)))


def log_energy(x: Any) -> float:
    z = np.asarray(x, dtype=np.complex128).ravel()
    if z.size == 0:
        return 0.0
    v = np.log1p(np.abs(z))
    v[~np.isfinite(v)] = 0.0
    return float(np.linalg.norm(v) / math.sqrt(max(1, z.size)))


def realness(x: Any) -> float:
    z = np.asarray(x, dtype=np.complex128).ravel()
    return float(np.linalg.norm(z.imag) / max(1e-300, np.linalg.norm(z)))


def stable_norm(x: Any) -> float:
    try:
        return float(np.linalg.norm(np.asarray(x, dtype=np.complex128)))
    except Exception:
        return float("inf")


def to_serializable_complex_array(x: Any) -> list[list[float]]:
    z = np.asarray(x, dtype=np.complex128).ravel()
    return [[float(v.real), float(v.imag)] for v in z]


def raw_direction(n: int, seed: int, spiral: bool = True) -> np.ndarray:
    """Deterministic complex pseudo-random direction."""
    rng = np.random.default_rng(int(seed) & 0xFFFFFFFF)
    if spiral:
        k = np.arange(int(n), dtype=float)
        theta = 2.399963229728653 * (k + 1.0) + 0.6180339887498948 * (seed % 10007)
        base = np.cos(theta) + 1j * np.sin(theta)
        noise = 0.35 * (rng.normal(size=n) + 1j * rng.normal(size=n))
        v = base + noise
    else:
        v = rng.normal(size=n) + 1j * rng.normal(size=n)
    nv = max(1e-300, float(np.linalg.norm(v)))
    return np.asarray(v / nv, dtype=np.complex128)


# ---------------------------------------------------------------------------
# Configuration and results
# ---------------------------------------------------------------------------

@dataclass
class PandrosionConfig:
    max_epochs: int = 16
    accept: float = 1e-10
    tol: float = 1e-12

    probe_scale: float = 0.035
    probe_candidates: int = 10
    probe_radii: tuple[float, ...] = (0.0, 0.5, 1.0, 2.0, 4.0)
    include_self_probe: bool = True
    axis_probes: bool = True
    condition_top: int = 3
    curvature_top: int = 2

    equal_value_weight: float = 0.015
    equal_value_eps: float = 1e-12
    condition_weight: float = 0.0025
    curvature_weight: float = 0.003
    curvature_mid: float = 0.5
    log_weight: float = 0.0005
    log_delta_weight: float = 0.0015
    singularity_weight: float = 0.002

    entropic_probe_top: int = 1
    entropic_probe_temperature: float = 0.025
    entropic_probe_blend: float = 0.0

    line_lambdas: tuple[float, ...] = (1.0, 0.75, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045)
    projective_line_weight: float = 0.0005

    trust_cond_weight: float = 0.06
    trust_cond_min: float = 4.0

    halley_enabled: bool = True
    halley_gate_residual: float = 0.25
    halley_probe_fraction: float = 0.50
    halley_cond_weight: float = 0.025
    halley_log_weight: float = 0.03
    halley_log_scale: float = 1.0
    halley_min_gate: float = 0.04
    halley_max_correction: float = 1.25
    tensor_extra_directions: int = 0

    spectral_filter: bool = True
    spectral_min_ratio: float = 1e-10
    spectral_shrink_power: float = 0.35

    seed: int = 215
    trial_timeout: float = 0.0


@dataclass
class SolveResult:
    accepted: bool
    ok: bool
    status: str
    residual: float
    epochs: int
    x: np.ndarray
    seconds: float
    eval_count: int
    slope_count: int
    line_eval_count: int
    probe_eval_count: int
    halley_used_count: int
    last_condition: Optional[float]
    metadata: dict[str, Any]


# ---------------------------------------------------------------------------
# Target wrapper
# ---------------------------------------------------------------------------

class PandrosionTarget:
    def __init__(self, func: Callable[[np.ndarray], np.ndarray], n: int, name: str = "target"):
        self.func = func
        self.n = int(n)
        self.name = str(name)
        self.eval_count = 0
        self.slope_count = 0

    def eval(self, x: Any) -> np.ndarray:
        z = np.asarray(x, dtype=np.complex128)
        self.eval_count += 1
        y = np.asarray(self.func(z), dtype=np.complex128).ravel()
        return y

    def eval_batch(self, X: Any) -> np.ndarray:
        XX = np.asarray(X, dtype=np.complex128)
        if XX.ndim == 1:
            XX = XX[None, :]
        out = []
        for row in XX:
            out.append(self.eval(row))
        return np.asarray(out, dtype=np.complex128)

    def residual(self, x: Any) -> float:
        try:
            return float(np.linalg.norm(self.eval(x)))
        except Exception:
            return float("inf")

    def residuals_batch(self, X: Any) -> np.ndarray:
        XX = np.asarray(X, dtype=np.complex128)
        if XX.ndim == 1:
            XX = XX[None, :]
        out = []
        for row in XX:
            try:
                out.append(float(np.linalg.norm(self.eval(row))))
            except Exception:
                out.append(float("inf"))
        return np.asarray(out, dtype=float)

    def slope_matrix(self, a: Any, b: Any) -> np.ndarray:
        self.slope_count += 1
        return finite_slope_matrix(self, a, b)


# ---------------------------------------------------------------------------
# Finite-slope matrix and diagnostics
# ---------------------------------------------------------------------------

def finite_slope_matrix(target: PandrosionTarget, a: Any, b: Any) -> np.ndarray:
    """Telescopic finite-slope matrix.

    It builds Q such that F(b)-F(a) is represented by coordinate finite slopes
    along a path from a to b. If one coordinate step is zero, a micro-probe is
    used. This is derivative-free and does not compute an analytic Jacobian.
    """
    aa = np.asarray(a, dtype=np.complex128).ravel()
    bb = np.asarray(b, dtype=np.complex128).ravel()
    f0 = target.eval(aa)
    m = len(f0)
    n = len(aa)
    Q = np.zeros((m, n), dtype=np.complex128)

    x_prev = aa.copy()
    f_prev = f0.copy()
    for j in range(n):
        x_next = x_prev.copy()
        step = bb[j] - x_prev[j]
        if abs(step) <= 1e-12 * max(1.0, abs(aa[j]), abs(bb[j])):
            step = (1e-7 + 1e-7j) * max(1.0, abs(aa[j]))
        x_next[j] = x_prev[j] + step
        f_next = target.eval(x_next)
        Q[:, j] = (f_next - f_prev) / step
        x_prev = x_next
        f_prev = f_next
    return Q


def singularity_score(Q: Any) -> tuple[float, float, float]:
    try:
        s = np.linalg.svd(np.asarray(Q, dtype=np.complex128), compute_uv=False)
        if s.size == 0:
            return float("inf"), 0.0, float("inf")
        smin = float(np.min(s))
        smax = float(np.max(s))
        cond = smax / max(1e-14, smin)
        score = math.log1p(1.0 / max(smin, 1e-14)) + 0.15 * math.log1p(max(0.0, cond))
        return float(score), float(smin), float(cond)
    except Exception:
        return float("inf"), 0.0, float("inf")


def softmax_weights(values: Any, temperature: float) -> np.ndarray:
    v = np.asarray(values, dtype=float)
    temp = max(1e-12, float(temperature))
    finite = np.isfinite(v)
    if not np.any(finite):
        return np.ones_like(v) / max(1, v.size)
    m = float(np.nanmin(v[finite]))
    z = np.where(finite, -(v - m) / temp, -1e300)
    z = z - float(np.max(z))
    w = np.exp(z)
    sw = float(np.sum(w))
    if not math.isfinite(sw) or sw <= 0:
        return np.ones_like(v) / max(1, v.size)
    return w / sw


# ---------------------------------------------------------------------------
# Probe selection
# ---------------------------------------------------------------------------

def build_probe_candidates(
    y: np.ndarray,
    prev_delta: Optional[np.ndarray],
    ep: int,
    direction_seed: int,
    cfg: PandrosionConfig,
) -> tuple[list[str], np.ndarray]:
    n = len(y)
    ynorm = max(1.0, float(np.linalg.norm(y)))
    radii = [float(r) for r in cfg.probe_radii if float(r) >= 0]
    if not radii:
        radii = [1.0]

    budget = max(1, int(cfg.probe_candidates))
    candidates: list[tuple[str, np.ndarray]] = []

    if cfg.include_self_probe:
        candidates.append(("self", y.copy()))

    if prev_delta is not None and np.all(np.isfinite(prev_delta)):
        pdn = max(1e-300, float(np.linalg.norm(prev_delta)))
        base = prev_delta / pdn * min(max(pdn, cfg.probe_scale * ynorm), 2.5 * ynorm)
        candidates.append(("inertial", y + base))

    if cfg.axis_probes:
        rad0 = cfg.probe_scale * ynorm * (radii[min(1, len(radii) - 1)] if radii else 1.0)
        for j in range(n):
            if len(candidates) >= budget:
                break
            ej = np.zeros(n, dtype=np.complex128)
            ej[j] = 1.0
            ph = phase(0.173 + 0.6180339887498948 * (ep + 1) + 1.41421356237 * (j + 1))
            candidates.append((f"axis-{j}", y + rad0 * ph * ej))

    k = 0
    while len(candidates) < budget:
        rad = cfg.probe_scale * ynorm * radii[k % len(radii)]
        qdir = raw_direction(
            n,
            direction_seed + 104729 * (ep + 1) + 7919 * (k + 1),
            True,
        )
        qdir = qdir / max(1e-300, float(np.linalg.norm(qdir))) * math.sqrt(max(1, n))
        ph = phase(0.6180339887498948 * (ep + 1) + 2.399963229728653 * (k + 1))
        step = rad * ph * qdir
        if float(np.linalg.norm(y)) > 0:
            step = step + (0.12 * rad) * y / ynorm * phase(0.38196601125 * (k + 1))
        tiny = 1e-12 * ynorm
        for j in range(n):
            if abs(step[j]) < tiny:
                step[j] += tiny * phase(0.17 + j + ep + k)
        candidates.append((f"geom-{k}", y + step))
        k += 1

    cand_slice = candidates[:budget]
    names = [name for name, _ in cand_slice]
    B = np.asarray([np.asarray(b, dtype=np.complex128) for _, b in cand_slice], dtype=np.complex128)
    return names, B


def choose_probe(
    target: PandrosionTarget,
    y: np.ndarray,
    f: np.ndarray,
    residual: float,
    prev_delta: Optional[np.ndarray],
    ep: int,
    direction_seed: int,
    cfg: PandrosionConfig,
) -> tuple[np.ndarray, dict[str, Any]]:
    names, B = build_probe_candidates(y, prev_delta, ep, direction_seed, cfg)
    FB = target.eval_batch(B)
    R = np.linalg.norm(FB, axis=1)
    D = np.linalg.norm(B - y[None, :], axis=1)
    ynorm = max(1.0, float(np.linalg.norm(y)))

    scores = np.asarray([
        math.log1p(max(0.0, float(R[i]))) + 1e-14 * math.log1p(float(D[i]))
        if math.isfinite(float(R[i])) else float("inf")
        for i in range(len(B))
    ], dtype=float)

    # log-stable/projective score
    ly = log_energy(y)
    lB = np.asarray([log_energy(B[i]) for i in range(len(B))], dtype=float)
    lD = np.abs(lB - ly)
    if cfg.log_weight > 0:
        scores += float(cfg.log_weight) * np.log1p(lB)
    if cfg.log_delta_weight > 0:
        scores += float(cfg.log_delta_weight) * np.log1p(lD)

    # equal-value pole avoidance
    Fy_norm = float(np.linalg.norm(f))
    Fdiff = np.linalg.norm(FB - np.asarray(f, dtype=np.complex128)[None, :], axis=1)
    Fscale = R + Fy_norm + 1e-300
    equal_gap = Fdiff / Fscale
    nonself = D > (1e-12 * ynorm)
    pole_penalty = np.zeros_like(scores)
    if cfg.equal_value_weight > 0:
        safe_gap = np.maximum(equal_gap, cfg.equal_value_eps)
        pole_penalty = np.where(nonself, cfg.equal_value_weight * np.log1p(1.0 / safe_gap), 0.0)
        scores += pole_penalty

    cond_best = None
    cond_scored = 0
    singularity_best = None
    singularity_sigma_min_best = None
    singularity_scored = 0
    curvature_best = None
    curvature_scored = 0

    # condition + singularity scoring over top probes
    if cfg.condition_top > 0 and np.any(np.isfinite(scores)):
        order = np.argsort(scores)
        for idx in order[:max(0, min(cfg.condition_top, len(order)))]:
            ii = int(idx)
            try:
                Q = target.slope_matrix(y, B[ii])
                cond = float(np.linalg.cond(Q))
                cond_scored += 1
                if math.isfinite(cond):
                    scores[ii] += cfg.condition_weight * math.log1p(max(0.0, cond))
                else:
                    scores[ii] = float("inf")
                if cfg.singularity_weight > 0:
                    ss, _smin, _cond2 = singularity_score(Q)
                    singularity_scored += 1
                    if math.isfinite(ss):
                        scores[ii] += cfg.singularity_weight * ss
            except Exception:
                scores[ii] = float("inf")

    # curvature proxy over top probes
    if cfg.curvature_top > 0 and cfg.curvature_weight > 0 and np.any(np.isfinite(scores)):
        order = np.argsort(scores)
        curv_vals = []
        for idx in order[:max(0, min(cfg.curvature_top, len(order)))]:
            ii = int(idx)
            if not bool(nonself[ii]):
                continue
            try:
                bvec = B[ii]
                mid = y + cfg.curvature_mid * (bvec - y)
                Q_full = target.slope_matrix(y, bvec)
                Q_mid = target.slope_matrix(y, mid)
                curv = float(np.linalg.norm(Q_full - Q_mid) / max(1e-300, np.linalg.norm(Q_full)))
                curvature_scored += 1
                if math.isfinite(curv):
                    scores[ii] += cfg.curvature_weight * math.log1p(max(0.0, curv))
                    curv_vals.append(curv)
                else:
                    scores[ii] = float("inf")
            except Exception:
                scores[ii] = float("inf")
        if curv_vals:
            curvature_best = float(min(curv_vals))

    # optional entropic blending
    entropic_used = False
    entropic_entropy = None
    if cfg.entropic_probe_top > 1 and cfg.entropic_probe_blend > 0 and len(B) > 1 and np.any(np.isfinite(scores)):
        finite_idx = np.where(np.isfinite(scores))[0]
        top = finite_idx[np.argsort(scores[finite_idx])[:max(1, min(cfg.entropic_probe_top, len(finite_idx)))]]
        if len(top) > 1:
            w = softmax_weights(scores[top], cfg.entropic_probe_temperature)
            entropy = float(-np.sum(w * np.log(np.maximum(w, 1e-300))))
            blended = np.sum(B[top] * w[:, None], axis=0)
            base_best = B[int(np.nanargmin(scores))]
            blend = max(0.0, min(1.0, cfg.entropic_probe_blend))
            candidate = (1.0 - blend) * base_best + blend * blended
            try:
                fc = target.eval(candidate)
                rc = float(np.linalg.norm(fc))
                gap = float(np.linalg.norm(fc - f) / max(1e-300, rc + Fy_norm))
                sc = math.log1p(max(0.0, rc))
                sc += cfg.equal_value_weight * math.log1p(1.0 / max(cfg.equal_value_eps, gap))
                sc += cfg.log_weight * math.log1p(log_energy(candidate))
                if math.isfinite(sc) and sc <= float(np.nanmin(scores)) * (1.0 + 0.15) + 1e-12:
                    B = np.vstack([B, candidate[None, :]])
                    R = np.concatenate([R, [rc]])
                    D = np.concatenate([D, [float(np.linalg.norm(candidate - y))]])
                    names.append("entropic-blend")
                    scores = np.concatenate([scores, [sc]])
                    equal_gap = np.concatenate([equal_gap, [gap]])
                    lB = np.concatenate([lB, [log_energy(candidate)]])
                    lD = np.concatenate([lD, [abs(log_energy(candidate) - ly)]])
                    pole_penalty = np.concatenate([pole_penalty, [0.0]])
                    entropic_used = True
                    entropic_entropy = entropy
            except Exception:
                pass

    if not np.any(np.isfinite(scores)):
        raise RuntimeError("no-finite-probe")

    best_idx = int(np.nanargmin(scores))
    best_b = B[best_idx].copy()
    best_name = names[best_idx]

    try:
        Qbest = target.slope_matrix(y, best_b)
        cond_best = float(np.linalg.cond(Qbest))
        ss, smin, _ = singularity_score(Qbest)
        singularity_best = float(ss) if math.isfinite(ss) else None
        singularity_sigma_min_best = float(smin) if math.isfinite(smin) else None
    except Exception:
        cond_best = None

    meta = {
        "probe_mode": "215-logstable-equalvalue-singularity-aware",
        "probe_name": best_name,
        "probe_candidates": int(len(B)),
        "probe_evals": int(len(B)),
        "probe_residual": float(R[best_idx]) if best_idx < len(R) else None,
        "probe_distance": float(D[best_idx]) if best_idx < len(D) else None,
        "probe_improvement_proxy": (float(residual / R[best_idx]) if best_idx < len(R) and R[best_idx] > 0 and math.isfinite(residual) else None),
        "probe_equal_value_gap": float(equal_gap[best_idx]) if best_idx < len(equal_gap) else None,
        "probe_equal_value_gap_min": float(np.nanmin(equal_gap)) if len(equal_gap) else None,
        "probe_equal_value_penalty": float(pole_penalty[best_idx]) if best_idx < len(pole_penalty) else 0.0,
        "probe_log_energy_current": float(ly),
        "probe_log_energy_best": float(lB[best_idx]) if best_idx < len(lB) else None,
        "probe_log_delta_best": float(lD[best_idx]) if best_idx < len(lD) else None,
        "probe_condition_best": cond_best,
        "probe_condition_scored": int(cond_scored),
        "probe_curvature_best": curvature_best,
        "probe_curvature_scored": int(curvature_scored),
        "probe_singularity_best": singularity_best,
        "probe_singularity_sigma_min_best": singularity_sigma_min_best,
        "probe_singularity_scored": int(singularity_scored),
        "probe_entropic_used": bool(entropic_used),
        "probe_entropic_entropy": entropic_entropy,
    }
    return best_b, meta


# ---------------------------------------------------------------------------
# Gated finite-slope Halley correction
# ---------------------------------------------------------------------------

def halley_gate(residual: float, cond: Optional[float], le: float, cfg: PandrosionConfig) -> float:
    r = max(0.0, float(residual))
    rg = 1.0 / (1.0 + (r / max(1e-300, cfg.halley_gate_residual)) ** 2)
    if cond is None or not math.isfinite(float(cond)):
        cg = 0.0
    else:
        cg = 1.0 / (1.0 + cfg.halley_cond_weight * math.log1p(max(0.0, float(cond))))
    lg = 1.0 / (1.0 + cfg.halley_log_weight * max(0.0, le) / max(1e-12, cfg.halley_log_scale))
    g = max(0.0, min(1.0, rg * cg * lg))
    return float(g if g >= cfg.halley_min_gate else 0.0)


def gated_halley_delta(
    target: PandrosionTarget,
    y: np.ndarray,
    f: np.ndarray,
    Q: np.ndarray,
    delta_raw: np.ndarray,
    residual: float,
    cond: Optional[float],
    cfg: PandrosionConfig,
) -> tuple[np.ndarray, dict[str, Any]]:
    gate = halley_gate(residual, cond, log_energy(y), cfg)
    meta = {
        "halley_gate": float(gate),
        "halley_used": False,
        "halley_delta2_norm": None,
        "halley_raw_norm": float(np.linalg.norm(delta_raw)),
        "halley_tensor_terms": 0,
        "spectral_ratio": None,
        "spectral_shrink": None,
    }
    if gate <= 0.0:
        return delta_raw, meta

    dnorm = float(np.linalg.norm(delta_raw))
    ynorm = max(1.0, float(np.linalg.norm(y)))
    if not math.isfinite(dnorm) or dnorm <= 1e-300:
        return delta_raw, meta

    frac = max(1e-6, min(1.0, cfg.halley_probe_fraction))
    hvec = frac * delta_raw

    try:
        fp = target.eval(y + hvec)
        fm = target.eval(y - hvec)
        sec_terms = [(fp - 2.0 * f + fm) / (frac * frac)]

        # optional extra tensor directions
        extra = max(0, int(cfg.tensor_extra_directions))
        for kk in range(extra):
            if kk == 0:
                qdir = 1j * delta_raw
            else:
                qdir = raw_direction(len(y), cfg.seed + 0x215000 + 7919 * (kk + 1), True)
                qdir = qdir / max(1e-300, float(np.linalg.norm(qdir))) * dnorm
            h2 = frac * qdir
            try:
                fp2 = target.eval(y + h2)
                fm2 = target.eval(y - h2)
                sec_terms.append((fp2 - 2.0 * f + fm2) / (frac * frac))
            except Exception:
                pass

        sec2 = np.mean(np.asarray(sec_terms, dtype=np.complex128), axis=0)
        delta2 = np.linalg.solve(Q, -0.5 * sec2)
        if not np.all(np.isfinite(delta2)):
            return delta_raw, meta

        spectral_ratio = None
        spectral_shrink = 1.0
        if cfg.spectral_filter:
            try:
                sv = np.linalg.svd(np.asarray(Q, dtype=np.complex128), compute_uv=False)
                if sv.size:
                    smax = float(np.max(sv))
                    smin = float(np.min(sv))
                    spectral_ratio = smin / max(smax, 1e-300)
                    if spectral_ratio < cfg.spectral_min_ratio:
                        spectral_shrink = max(0.05, spectral_ratio / max(cfg.spectral_min_ratio, 1e-300))
                        delta2 = delta2 * spectral_shrink
            except Exception:
                pass

        d2norm = float(np.linalg.norm(delta2))
        max_d2 = cfg.halley_max_correction * max(dnorm, 1e-300)
        if d2norm > max_d2:
            delta2 = delta2 * (max_d2 / max(d2norm, 1e-300))
            d2norm = float(np.linalg.norm(delta2))

        if dnorm + d2norm > 20.0 * ynorm:
            scale = (20.0 * ynorm) / max(dnorm + d2norm, 1e-300)
            delta2 *= scale
            d2norm = float(np.linalg.norm(delta2))

        delta_halley = delta_raw + delta2
        delta_mix = (1.0 - gate) * delta_raw + gate * delta_halley
        meta.update({
            "halley_used": True,
            "halley_delta2_norm": float(d2norm),
            "halley_sec2_norm": float(np.linalg.norm(sec2)),
            "halley_tensor_terms": int(len(sec_terms)),
            "halley_mix_norm": float(np.linalg.norm(delta_mix)),
            "spectral_ratio": (float(spectral_ratio) if spectral_ratio is not None and math.isfinite(float(spectral_ratio)) else None),
            "spectral_shrink": float(spectral_shrink),
        })
        return delta_mix, meta
    except Exception as exc:
        meta["halley_error"] = type(exc).__name__
        return delta_raw, meta


# ---------------------------------------------------------------------------
# Main solve engine
# ---------------------------------------------------------------------------

def solve_pandrosion(
    func: Callable[[np.ndarray], np.ndarray],
    x0: Sequence[complex],
    config: Optional[PandrosionConfig] = None,
    name: str = "target",
) -> SolveResult:
    cfg = config or PandrosionConfig()
    target = PandrosionTarget(func, n=len(x0), name=name)

    y = np.asarray(x0, dtype=np.complex128).copy()
    t0 = now()
    deadline = t0 + float(cfg.trial_timeout) if cfg.trial_timeout and cfg.trial_timeout > 0 else None

    best_y = y.copy()
    best_r = target.residual(y)
    ok = False
    status = "started"
    epochs = 0
    prev_delta = None
    last_cond = None
    last_probe_meta: dict[str, Any] = {}
    last_halley_meta: dict[str, Any] = {}
    total_probe_evals = 0
    total_line_evals = 0
    halley_used_count = 0

    for ep in range(max(1, int(cfg.max_epochs))):
        if deadline is not None and now() > deadline:
            status = "timeout"
            break

        try:
            f = target.eval(y)
            r = float(np.linalg.norm(f))
        except Exception as exc:
            status = f"eval-error:{type(exc).__name__}"
            break

        if math.isfinite(r) and r < best_r:
            best_r = r
            best_y = y.copy()

        if r <= max(float(cfg.tol), float(cfg.accept)):
            ok = True
            status = "converged"
            break

        try:
            b, pmeta = choose_probe(
                target=target,
                y=y,
                f=f,
                residual=r,
                prev_delta=prev_delta,
                ep=ep,
                direction_seed=cfg.seed + 7919 * ep,
                cfg=cfg,
            )
            last_probe_meta = pmeta
            total_probe_evals += int(pmeta.get("probe_evals", 0))
        except Exception as exc:
            status = f"probe-error:{type(exc).__name__}"
            break

        try:
            Q = target.slope_matrix(y, b)
            last_cond = float(np.linalg.cond(Q))
            delta_raw = np.linalg.solve(Q, -f)
        except Exception as exc:
            status = f"slope-solve-error:{type(exc).__name__}"
            break

        if not np.all(np.isfinite(delta_raw)):
            status = "nonfinite-step"
            break

        ynorm = max(1.0, float(np.linalg.norm(y)))
        dnorm = float(np.linalg.norm(delta_raw))
        cap_factor = 18.0
        if cfg.trust_cond_weight > 0 and last_cond is not None and math.isfinite(last_cond):
            cap_factor = max(cfg.trust_cond_min, 18.0 / (1.0 + cfg.trust_cond_weight * math.log1p(max(0.0, last_cond))))
        if dnorm > cap_factor * ynorm:
            delta_raw *= (cap_factor * ynorm) / max(dnorm, 1e-300)

        delta = delta_raw
        if cfg.halley_enabled:
            delta, hmeta = gated_halley_delta(target, y, f, Q, delta_raw, r, last_cond, cfg)
            last_halley_meta = hmeta
            if hmeta.get("halley_used"):
                halley_used_count += 1
        else:
            last_halley_meta = {"halley_gate": 0.0, "halley_used": False}

        if not np.all(np.isfinite(delta)):
            status = "nonfinite-halley-step"
            break

        L = np.asarray(cfg.line_lambdas, dtype=float)
        Yraw = y[None, :] + L[:, None] * delta_raw[None, :]
        if cfg.halley_enabled and last_halley_meta.get("halley_used"):
            Yhal = y[None, :] + L[:, None] * delta[None, :]
            Ytry = np.vstack([Yraw, Yhal])
            direction_tags = np.asarray(["raw"] * len(L) + ["halley"] * len(L), dtype=object)
            lambda_tags = np.concatenate([L, L])
        else:
            Ytry = Yraw
            direction_tags = np.asarray(["raw"] * len(L), dtype=object)
            lambda_tags = L

        Rtry = target.residuals_batch(Ytry)
        total_line_evals += int(len(Ytry))
        finite = np.isfinite(Rtry)
        better = finite & ((Rtry < r) | (Rtry < best_r))
        if not np.any(better):
            status = "no-decrease"
            break

        scores = np.where(better, Rtry * (1.0 + 1e-15 / np.maximum(lambda_tags, 1e-15)), np.inf)
        if cfg.projective_line_weight > 0:
            try:
                line_log = np.asarray([log_energy(v) for v in Ytry], dtype=float)
                scores = scores + np.where(better, cfg.projective_line_weight * np.log1p(line_log), 0.0)
            except Exception:
                pass

        idx = int(np.nanargmin(scores))
        lam = float(lambda_tags[idx])
        rr = float(Rtry[idx])
        yy = Ytry[idx].copy()
        chosen_direction = str(direction_tags[idx])

        prev_delta = lam * (delta if chosen_direction == "halley" else delta_raw)
        y = yy
        if rr < best_r:
            best_y = yy.copy()
            best_r = rr
        epochs = ep + 1
        last_halley_meta["halley_chosen_direction"] = chosen_direction
        last_halley_meta["halley_chosen_lambda"] = lam
    else:
        status = "max-epochs"

    final_r = target.residual(best_y)
    if final_r <= max(float(cfg.tol), float(cfg.accept)):
        ok = True
        status = "converged"

    meta = {
        "target": target.name,
        "config": asdict(cfg),
        "corrector": "215-core-logstable-equalvalue-singularity-aware-spectral-halley-pandrosion",
        "realness": realness(best_y),
        "log_energy": log_energy(best_y),
        "line_lambdas": [float(x) for x in cfg.line_lambdas],
        **last_probe_meta,
        **last_halley_meta,
    }

    return SolveResult(
        accepted=bool(final_r <= cfg.accept),
        ok=bool(ok),
        status=status,
        residual=float(final_r),
        epochs=int(epochs),
        x=best_y,
        seconds=float(now() - t0),
        eval_count=int(target.eval_count),
        slope_count=int(target.slope_count),
        line_eval_count=int(total_line_evals),
        probe_eval_count=int(total_probe_evals),
        halley_used_count=int(halley_used_count),
        last_condition=last_cond,
        metadata=meta,
    )


# ---------------------------------------------------------------------------
# Scalar engine
# ---------------------------------------------------------------------------

def solve_scalar_pandrosion(
    f_scalar: Callable[[complex], complex],
    x0: complex,
    config: Optional[PandrosionConfig] = None,
    name: str = "scalar",
) -> tuple[complex, SolveResult]:
    def F(z: np.ndarray) -> np.ndarray:
        return np.asarray([f_scalar(z[0])], dtype=np.complex128)

    res = solve_pandrosion(F, [complex(x0)], config=config, name=name)
    return complex(res.x[0]), res


# ---------------------------------------------------------------------------
# Root collection helper for self-tests
# ---------------------------------------------------------------------------

def cluster_index(roots: list[np.ndarray], x: np.ndarray, sep: float) -> Optional[int]:
    if not roots:
        return None
    xx = np.asarray(x, dtype=np.complex128)
    for i, r in enumerate(roots):
        if float(np.linalg.norm(xx - r)) <= sep * max(1.0, float(np.linalg.norm(r))):
            return i
    return None


def collect_roots(
    func: Callable[[np.ndarray], np.ndarray],
    n: int,
    count: int,
    pool: int,
    cfg: PandrosionConfig,
    cluster_sep: float = 1e-8,
    radius: float = 1.2,
    name: str = "collect",
) -> dict[str, Any]:
    rng = np.random.default_rng(int(cfg.seed))
    roots: list[np.ndarray] = []
    root_records: list[dict[str, Any]] = []
    trials_used = 0
    failures = 0
    duplicates = 0
    t0 = now()

    # deterministic cyclotomic-ish starts first
    starts = []
    for k in range(int(pool)):
        theta = 2.0 * math.pi * k / max(1, int(pool))
        v = np.zeros(n, dtype=np.complex128)
        for j in range(n):
            v[j] = radius * phase(theta * (j + 1) + 0.37 * j)
        starts.append(v)

    for trial, x0 in enumerate(starts):
        if len(roots) >= int(count):
            break
        local_cfg = PandrosionConfig(**asdict(cfg))
        local_cfg.seed = int(cfg.seed) + 104729 * (trial + 1)
        res = solve_pandrosion(func, x0, config=local_cfg, name=f"{name}-trial{trial}")
        trials_used += 1
        if res.accepted and math.isfinite(res.residual):
            idx = cluster_index(roots, res.x, cluster_sep)
            if idx is None:
                roots.append(res.x.copy())
                root_records.append({
                    "root_id": len(roots) - 1,
                    "trial": trial,
                    "residual": float(res.residual),
                    "epochs": int(res.epochs),
                    "seconds": float(res.seconds),
                    "eval_count": int(res.eval_count),
                    "slope_count": int(res.slope_count),
                    "halley_used_count": int(res.halley_used_count),
                    "realness": float(res.metadata.get("realness", 0.0)),
                    "x": to_serializable_complex_array(res.x),
                    "status": res.status,
                })
            else:
                duplicates += 1
        else:
            failures += 1

    return {
        "summary": {
            "requested_roots": int(count),
            "unique_roots": int(len(roots)),
            "success": bool(len(roots) >= int(count)),
            "trials_used": int(trials_used),
            "duplicates": int(duplicates),
            "failures": int(failures),
            "total_seconds": float(now() - t0),
            "config": asdict(cfg),
        },
        "roots": root_records,
    }


# ---------------------------------------------------------------------------
# Self-test systems
# ---------------------------------------------------------------------------

def separated_power_system(n: int, d: int) -> Callable[[np.ndarray], np.ndarray]:
    """Simple n-dimensional system with d^n roots:
        F_j(z)=z_j^d-1.
    """
    n = int(n)
    d = int(d)

    def F(z: np.ndarray) -> np.ndarray:
        z = np.asarray(z, dtype=np.complex128).ravel()
        return np.asarray([z[j] ** d - 1.0 for j in range(n)], dtype=np.complex128)
    return F


def coupled_power_system(n: int, d: int, coupling: float = 0.05) -> Callable[[np.ndarray], np.ndarray]:
    """Weakly coupled polynomial system for engine testing."""
    n = int(n)
    d = int(d)
    c = complex(coupling)

    def F(z: np.ndarray) -> np.ndarray:
        z = np.asarray(z, dtype=np.complex128).ravel()
        out = []
        for j in range(n):
            neighbor = z[(j + 1) % n]
            out.append(z[j] ** d - 1.0 + c * (neighbor - z[j]))
        return np.asarray(out, dtype=np.complex128)
    return F


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_cases(s: str) -> list[tuple[int, int]]:
    out = []
    for block in str(s).replace(";", " ").split():
        if "," in block:
            a, b = block.split(",", 1)
        elif "x" in block:
            a, b = block.split("x", 1)
        else:
            continue
        try:
            n = int(a.strip())
            d = int(b.strip())
            if n > 0 and d > 0:
                out.append((n, d))
        except Exception:
            pass
    return out or [(2, 4)]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="215 standalone NumPy core Pandrosion engines only; no trading/backtest.")
    p.add_argument("--cases", default="2,2 2,3 2,4", help="cases n,d separated by spaces or semicolons")
    p.add_argument("--system", choices=["separated", "coupled"], default="separated")
    p.add_argument("--count", type=int, default=8)
    p.add_argument("--pool", type=int, default=256)
    p.add_argument("--epochs", type=int, default=16)
    p.add_argument("--accept", type=float, default=1e-10)
    p.add_argument("--seed", type=int, default=215)
    p.add_argument("--probe-candidates", type=int, default=10)
    p.add_argument("--probe-scale", type=float, default=0.035)
    p.add_argument("--cluster-sep", type=float, default=1e-8)
    p.add_argument("--no-halley", dest="halley_enabled", action="store_false")
    p.add_argument("--halley", dest="halley_enabled", action="store_true", default=True)
    p.add_argument("--tensor-extra-directions", type=int, default=0)
    p.add_argument("--entropic-probe-top", type=int, default=1)
    p.add_argument("--entropic-probe-blend", type=float, default=0.0)
    p.add_argument("--out", default="/mnt/data/215_pandrosion_core_selftest.json")
    p.add_argument("--self-test", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    cfg = PandrosionConfig(
        max_epochs=int(args.epochs),
        accept=float(args.accept),
        seed=int(args.seed),
        probe_candidates=int(args.probe_candidates),
        probe_scale=float(args.probe_scale),
        halley_enabled=bool(args.halley_enabled),
        tensor_extra_directions=int(args.tensor_extra_directions),
        entropic_probe_top=int(args.entropic_probe_top),
        entropic_probe_blend=float(args.entropic_probe_blend),
    )

    results = []
    for n, d in parse_cases(args.cases):
        F = separated_power_system(n, d) if args.system == "separated" else coupled_power_system(n, d)
        target_count = min(int(args.count), int(d) ** int(n))
        if args.self_test:
            target_count = min(target_count, 8)
        res = collect_roots(
            func=F,
            n=n,
            count=target_count,
            pool=int(args.pool),
            cfg=cfg,
            cluster_sep=float(args.cluster_sep),
            name=f"{args.system}-{n}-{d}",
        )
        res["case"] = {"n": int(n), "d": int(d), "system": args.system}
        results.append(res)
        s = res["summary"]
        print(f"case={args.system}({n},{d}) roots={s['unique_roots']}/{s['requested_roots']} success={s['success']} trials={s['trials_used']} seconds={s['total_seconds']:.3f}")

    payload = {
        "source": "215-pandrosion-core-engines-standalone-numpy",
        "description": "generic Pandrosion engines only; no trading/backtest",
        "results": results,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"out={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
