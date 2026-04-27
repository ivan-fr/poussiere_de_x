"""
PAPER: 011 (research piste: dD Pandrosion as multivariate homotopy corrector)
TITLE: Can dD Pandrosion structure improve multivariate root-finding?
STATUS: research-grade exploration

CORE IDEA
=========
The dD Pandrosion iteration has a STRIKING property (rank-1 collapse,
flow/004 Theorem 4): from ANY initial vector, ONE iteration projects
to the diagonal {s_1 = ... = s_n}.  In a homotopy context, this means
that if we use dD Pandrosion as the CORRECTOR step, paths that diverge
in different components are realigned in a single iteration -- preventing
the path-crossing failures that plague Lairez on Kostlan-Smale adversarial.

PROPOSED ALGORITHM
==================
Multivariate F: C^n -> C^n.  Standard Lairez homotopy:
    G_i(z) = z_i^d - 1,
    H(z, t) = gamma (1-t) G(z) + t F(z).
At each tracking step (predictor + corrector), we use:

    PREDICTOR: standard Euler  z_pred = z + h * dz/dt
    CORRECTOR: dD-PANDROSION step instead of Newton

The dD-Pandrosion corrector at fixed t for the diagonal-symmetric subset:
    s_i_new = 1 - (H_i(z, t) - 1) / (some scalar Delta_d-like quantity)
applied componentwise with rank-1 collapse forcing components to align.

We test this on Kostlan-Smale random systems at small (n, d) and compare
with plain Lairez Newton corrector.

KEY QUESTION
============
Does the rank-1 alignment property of dD Pandrosion provide better
robustness than Newton on adversarial systems?
"""
from __future__ import annotations
import math
import time
import warnings
from itertools import product

import numpy as np


# ---------------------------------------------------------------------------
# Multivariate KS system (paper 9 mv ensemble)
# ---------------------------------------------------------------------------

def make_kostlan_smale_mv(n: int, d: int, rng):
    """Random F : C^n -> C^n with Kostlan-Smale (Bombieri-Weyl) measure."""
    def gen_indices(remaining_d, depth):
        if depth == 0:
            yield (); return
        for k in range(remaining_d + 1):
            for r in gen_indices(remaining_d - k, depth - 1):
                yield (k,) + r
    indices = list(gen_indices(d, n))
    F = []
    for i in range(n):
        f = {}
        for alpha in indices:
            abs_a = sum(alpha)
            if abs_a > d: continue
            mc = math.factorial(d)
            den = math.factorial(d - abs_a)
            for a in alpha:
                den *= math.factorial(a)
            sigma = math.sqrt(mc // den)
            c = (rng.standard_normal() + 1j * rng.standard_normal()) * sigma
            f[alpha] = c
        F.append(f)
    return F


def evaluate_F(F, z):
    n = len(F)
    out = np.zeros(n, dtype=complex)
    for i, f in enumerate(F):
        s = 0.0 + 0j
        for alpha, c in f.items():
            term = c
            for j, a in enumerate(alpha):
                if a:
                    term *= z[j]**a
            s += term
        out[i] = s
    return out


def jacobian_F(F, z, h=1e-7):
    n = len(z)
    F0 = evaluate_F(F, z)
    J = np.zeros((n, n), dtype=complex)
    for j in range(n):
        z_h = z.copy(); z_h[j] += h
        J[:, j] = (evaluate_F(F, z_h) - F0) / h
    return J


# ---------------------------------------------------------------------------
# Standard Lairez homotopy with Newton corrector (baseline)
# ---------------------------------------------------------------------------

def lairez_newton_track(F, z0, gamma, d, dt_init=0.05, dt_min=1e-6,
                        max_steps=400, tol=1e-10):
    z = np.asarray(z0, dtype=complex).copy()
    t = 0.0
    dt = dt_init
    fails = 0
    for _ in range(max_steps):
        if t >= 1.0 - 1e-9:
            return z, True
        h = min(dt, 1.0 - t)
        # Predictor: dz/dt = -DH^{-1} * dH/dt
        try:
            J = gamma * (1 - t) * np.diag(d * z**(d-1)) + t * jacobian_F(F, z)
            G = z**d - 1
            Ht = -gamma * G + evaluate_F(F, z)
            dzdt = -np.linalg.solve(J, Ht)
        except np.linalg.LinAlgError:
            return z, False
        z_pred = z + h * dzdt
        # Corrector: Newton on H(., t+h) = 0
        z_corr = z_pred.copy()
        ok = False
        for _ in range(10):
            G = z_corr**d - 1
            Hv = gamma * (1 - t - h) * G + (t + h) * evaluate_F(F, z_corr)
            if np.linalg.norm(Hv) < tol:
                ok = True; break
            try:
                Jc = gamma * (1 - t - h) * np.diag(d * z_corr**(d-1)) \
                     + (t + h) * jacobian_F(F, z_corr)
                z_corr -= np.linalg.solve(Jc, Hv)
            except np.linalg.LinAlgError:
                break
        if ok:
            z = z_corr
            t = t + h
            dt = min(dt * 1.5, 0.2)
            fails = 0
        else:
            fails += 1
            if fails >= 5 or dt < dt_min:
                return z, False
            dt *= 0.5
    return z, t >= 1.0 - 1e-9


def find_all_lairez(F, n, d, gamma=None, tol=1e-10):
    if gamma is None:
        gamma = complex(0.6, 0.8)
    omega = np.exp(2j * np.pi * np.arange(d) / d)
    found = []
    for k_tuple in product(range(d), repeat=n):
        z0 = np.array([omega[k] for k in k_tuple], dtype=complex)
        z_end, _ = lairez_newton_track(F, z0, gamma, d, tol=tol)
        found.append(z_end)
    return found


# ---------------------------------------------------------------------------
# dD-Pandrosion-style corrector on H
# ---------------------------------------------------------------------------

def pandrosion_dD_corrector_on_H(F, z, t, h, gamma, d, max_iter=15, tol=1e-10):
    """Pandrosion-style corrector: at fixed t' = t + h, iterate
        z_i_new = 1 - (z_i_old - 1) / (some scalar quantity)
    by analogy with paper-0 sec:11.

    For our purposes, we adapt as: at each corrector step, compute the
    Newton step delta = J^{-1} H_v, then APPLY a Pandrosion-flavoured
    update that uses a SCALAR damper (analogous to the rank-1 dD collapse).

    Concretely, let u = J^{-1} H_v (Newton direction) and let
    sigma = mean(u) (a scalar summary).  Then
        z_corr = z_old - sigma * 1   (project Newton step onto diagonal)
                + (u - sigma * 1)    (residual off-diagonal correction)
    The full step equals u (= Newton), but a damped version with
    weight w in [0, 1] interpolates between Newton (w=1) and
    diagonal-only (w=0).  Empirical tuning may help adversarial paths.
    """
    z = np.asarray(z, dtype=complex).copy()
    n = len(z)
    for it in range(max_iter):
        G = z**d - 1
        Hv = gamma * (1 - t - h) * G + (t + h) * evaluate_F(F, z)
        if np.linalg.norm(Hv) < tol:
            return z, True
        try:
            J = gamma * (1 - t - h) * np.diag(d * z**(d-1)) \
                + (t + h) * jacobian_F(F, z)
            u = np.linalg.solve(J, Hv)
        except np.linalg.LinAlgError:
            return z, False
        # dD-style: split into diagonal and off-diagonal modes
        # The diagonal mode = mean(u), off-diagonal = u - mean(u)
        # Damp the off-diagonal slightly (analogous to rank-1 collapse)
        u_mean = np.mean(u)
        u_off = u - u_mean
        # Apply dampening on off-diagonal (factor decreases with iter)
        w_off = 0.7  # tunable: 1 = pure Newton, 0 = diagonal-only
        delta = u_mean * np.ones(n, dtype=complex) + w_off * u_off
        z = z - delta
    return z, np.linalg.norm(Hv) < tol


def pandrosion_dD_track(F, z0, gamma, d, dt_init=0.05, dt_min=1e-6,
                       max_steps=400, tol=1e-10):
    z = np.asarray(z0, dtype=complex).copy()
    t = 0.0
    dt = dt_init
    fails = 0
    for _ in range(max_steps):
        if t >= 1.0 - 1e-9:
            return z, True
        h = min(dt, 1.0 - t)
        try:
            J = gamma * (1 - t) * np.diag(d * z**(d-1)) + t * jacobian_F(F, z)
            G = z**d - 1
            Ht = -gamma * G + evaluate_F(F, z)
            dzdt = -np.linalg.solve(J, Ht)
        except np.linalg.LinAlgError:
            return z, False
        z_pred = z + h * dzdt
        z_corr, ok = pandrosion_dD_corrector_on_H(F, z_pred, t, h, gamma, d, tol=tol)
        if ok:
            z = z_corr
            t = t + h
            dt = min(dt * 1.5, 0.2)
            fails = 0
        else:
            fails += 1
            if fails >= 5 or dt < dt_min:
                return z, False
            dt *= 0.5
    return z, t >= 1.0 - 1e-9


def find_all_pandrosion_dD(F, n, d, gamma=None, tol=1e-10):
    if gamma is None:
        gamma = complex(0.6, 0.8)
    omega = np.exp(2j * np.pi * np.arange(d) / d)
    found = []
    for k_tuple in product(range(d), repeat=n):
        z0 = np.array([omega[k] for k in k_tuple], dtype=complex)
        z_end, _ = pandrosion_dD_track(F, z0, gamma, d, tol=tol)
        found.append(z_end)
    return found


def max_residual(F, roots):
    if not roots:
        return float('inf')
    return max(float(np.linalg.norm(evaluate_F(F, np.asarray(r)))) for r in roots)


def main():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    print("=" * 76)
    print("dD-Pandrosion as multivariate corrector vs Lairez Newton baseline")
    print("=" * 76)
    rng = np.random.default_rng(2026)
    test_systems = []
    for _ in range(2):
        for n, d in [(2, 2), (2, 3), (3, 2)]:
            test_systems.append((n, d, make_kostlan_smale_mv(n, d, rng)))
    print(f"\n  {'(n, d)':>8}  {'method':<26}  {'time':>8}  {'max ||F||':>14}")
    print("-" * 76)
    for n, d, F in test_systems:
        D = d**n
        # Lairez baseline
        t0 = time.perf_counter()
        try:
            roots_L = find_all_lairez(F, n, d, tol=1e-10)
            t_L = time.perf_counter() - t0
            res_L = max_residual(F, roots_L)
        except Exception as e:
            t_L = float('nan'); res_L = float('inf')
        # dD-Pandrosion experiment
        t0 = time.perf_counter()
        try:
            roots_P = find_all_pandrosion_dD(F, n, d, tol=1e-10)
            t_P = time.perf_counter() - t0
            res_P = max_residual(F, roots_P)
        except Exception as e:
            t_P = float('nan'); res_P = float('inf')
        print(f"  ({n}, {d})    Lairez Newton              {t_L*1000:>7.2f}ms  {res_L:>14.2e}")
        print(f"  ({n}, {d})    dD-Pandrosion damped       {t_P*1000:>7.2f}ms  {res_P:>14.2e}")
        print()


if __name__ == "__main__":
    main()
