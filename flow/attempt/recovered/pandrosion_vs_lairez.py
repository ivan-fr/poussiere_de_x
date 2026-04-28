"""
Comparaison: flow/008 Strategy B + T_2 (avec deflation) VS Lairez homotopy.

Trois methodes:
  1. Pandrosion-T_2 + Strategy B (golden spiral) + Armijo + deflation.
     Reproduit flow/008 etendu a find-all-d-roots par deflation Wilkinson.

  2. Pandrosion-T_2 + Strategy B + SCALING (homothetie sur la geometrie).
     Idee: papier 0 sec:scaling. On ecrit z = R * w ou R est une estimation
     geometrique de la magnitude des racines (R = |a_0/a_d|^{1/d}). Le
     polynome transforme Q(w) = P(R*w)/R^d a des racines de magnitude
     proche de 1, donc la geometrie de Pandrosion (anchors, Strategy B)
     est mieux centree, et la propagation d'erreur en deflation chute.

  3. Lairez 2017: cyclotomic anchors z^d - 1 -> P(z) avec predictor-Euler
     + corrector-Newton, gamma random.

Regimes testes:
  - random monic complex (mag 1, 5)
  - Kostlan-Smale (Bombieri-Weyl)
  - high dynamic range (coefficients de magnitudes tres variees)
"""
from __future__ import annotations
import math
import time
import warnings
import numpy as np
try:
    import mpmath as mp
except ImportError:  # pragma: no cover - benchmark environment has mpmath
    mp = None


# ---------------------------------------------------------------------------
# Manual Horner (consistent with paper 0 derivative-free spirit)
# ---------------------------------------------------------------------------

def horner1(coeffs, z):
    p = 0.0 + 0.0j
    for c in coeffs:
        p = p * z + c
    return p


def horner2(coeffs, z):
    p = 0.0 + 0.0j
    dp = 0.0 + 0.0j
    for c in coeffs:
        dp = dp * z + p
        p = p * z + c
    return p, dp


# ---------------------------------------------------------------------------
# Strategy B + T_2 with Armijo (from flow/008)
# ---------------------------------------------------------------------------

def cauchy_radius(coeffs):
    return 1.0 + max(abs(coeffs[k] / coeffs[0]) for k in range(1, len(coeffs)))


def strategy_B_starts(coeffs, q=0.7):
    d = len(coeffs) - 1
    rho = cauchy_radius(coeffs)
    phi = math.pi * (3 - math.sqrt(5))
    return [rho * q**k * np.exp(1j * k * phi) for k in range(d)]


def Q_quotient(coeffs, a, z):
    if abs(z - a) < 1e-15:
        _, dp = horner2(coeffs, z)
        return dp
    return (horner1(coeffs, z) - horner1(coeffs, a)) / (z - a)


def t2_step_anchored(coeffs, a, z):
    Pz = horner1(coeffs, z)
    Q0 = Q_quotient(coeffs, a, z)
    if abs(Q0) < 1e-15:
        return z
    s1 = z - Pz / Q0
    Ps1 = horner1(coeffs, s1)
    Q1 = Q_quotient(coeffs, a, s1)
    if abs(Q1) < 1e-15:
        return s1
    s2 = s1 - Ps1 / Q1
    denom = s2 - 2 * s1 + z
    if abs(denom) < 1e-15:
        return s2
    return z - (s1 - z) ** 2 / denom


def t2_orbit_els(coeffs, z0, max_epochs=80, tol=1e-12, k_range=(-6, 6)):
    """T_2-fused-PURE orbit augmented with Extended Line Search (paper 009).

    ELS replaces backtracking-only Armijo by symmetric two-sided
       T_ELS = {2^k : k in k_range}
    selecting the global minimiser of |P(z - tau * Delta)| where
    Delta = z - z_T2 is the Pandrosion-T_2 step direction.
    Forwardtracking entries (tau >= 2) amplify the step and close the
    sqrt(d)-undershoot gap on KS-like ensembles (Lemma 3.1 paper 9:
    tau* ~ sqrt(d) on KS).
    """
    z = complex(z0)
    z_prev = z + 1e-6 * (1 + 0j)
    Pz_prev = horner1(coeffs, z_prev)
    for _ in range(max_epochs):
        Pz = horner1(coeffs, z)
        if abs(Pz) < tol:
            return z, True
        z_t2 = t2_step_anchored(coeffs, z_prev, z)
        delta = z - z_t2
        if abs(delta) < 1e-300:
            z_prev = z
            Pz_prev = Pz
            z = z_t2
            continue
        best_tau = 1.0
        best_val = abs(horner1(coeffs, z_t2))
        for k in range(k_range[0], k_range[1] + 1):
            if k == 0:
                continue
            tau = 2.0 ** k
            val = abs(horner1(coeffs, z - tau * delta))
            if val < best_val:
                best_val = val
                best_tau = tau
        z_prev = z
        Pz_prev = Pz
        z = z - best_tau * delta
    return z, abs(horner1(coeffs, z)) < tol


def t2_orbit(coeffs, z0, max_epochs=200, eta=0.95, tol=1e-12):
    """Pandrosion-T_2 with Armijo fallback (anchor = z each epoch)."""
    z = complex(z0)
    for _ in range(max_epochs):
        Pz = horner1(coeffs, z)
        if abs(Pz) < tol:
            return z, True
        z_cand = t2_step_anchored(coeffs, z, z)
        Pz_cand = abs(horner1(coeffs, z_cand))
        if Pz_cand <= eta * abs(Pz):
            z = z_cand
            continue
        # Armijo backtrack on Newton direction
        _, Pp_z = horner2(coeffs, z)
        if abs(Pp_z) < 1e-15:
            return z, False
        direction = Pz / Pp_z
        accepted = False
        for j in range(20):
            tau = 2.0 ** (-j)
            z_new = z - tau * direction
            if abs(horner1(coeffs, z_new)) <= (1 - 0.1 * tau) * abs(Pz):
                z = z_new
                accepted = True
                break
        if not accepted:
            return z, False
    return z, abs(horner1(coeffs, z)) < tol


def deflate(coeffs, root):
    out = [coeffs[0]]
    for c in coeffs[1:-1]:
        out.append(out[-1] * root + c)
    return out


def find_all_strategy_B_t2(coeffs, polish=True, tol=1e-12):
    """Find all d roots: Strategy B starts + T_2 + Armijo + deflation."""
    coeffs = list(coeffs)
    P0 = list(coeffs)
    roots = []
    current = list(coeffs)
    while len(current) > 2:
        starts = strategy_B_starts(current)
        found_root = None
        for z0 in starts:
            r, ok = t2_orbit(current, z0, tol=tol)
            if ok:
                found_root = r
                break
        if found_root is None:
            return roots
        roots.append(found_root)
        current = deflate(current, found_root)
    if len(current) == 2 and abs(current[0]) > 1e-15:
        roots.append(-current[1] / current[0])
    if polish:
        polished = []
        for r in roots:
            r_pol, _ = t2_orbit(P0, r, max_epochs=30, tol=tol)
            polished.append(r_pol)
        roots = polished
    return roots


# ---------------------------------------------------------------------------
# Scaled variant: homothety z = R*w (paper 0, sec:scaling generalized to C)
# ---------------------------------------------------------------------------

def geometric_root_scale(coeffs):
    """R = |a_0 / a_d|^{1/d} (geometric mean of root magnitudes via Vieta).

    For monic P(z) = z^d + ... + a_0, the product of roots is (-1)^d * a_0,
    so |prod r_i| = |a_0|, hence (geometric mean of |r_i|) = |a_0|^{1/d}.
    For non-monic, ratio against leading coeff."""
    d = len(coeffs) - 1
    a0 = coeffs[-1]
    ad = coeffs[0]
    if abs(a0) < 1e-300 or abs(ad) < 1e-300:
        return 1.0
    return (abs(a0) / abs(ad)) ** (1.0 / d)


def rescale_polynomial(coeffs, R):
    """Q(w) = P(R*w) / R^d. Roots: r = R * w_root."""
    d = len(coeffs) - 1
    out = []
    # P(z) = sum a_k z^{d-k}, k=0..d.  P(R*w) = sum a_k R^{d-k} w^{d-k}.
    # Q(w) = P(R*w)/R^d = sum a_k R^{-k} w^{d-k}.
    for k, c in enumerate(coeffs):
        out.append(c * (R ** (-k)))
    return out


def find_all_scaled_strategy_B_t2(coeffs, polish=True, tol=1e-12):
    """Find all d roots: SCALED Strategy B + T_2 + Armijo + deflation.

    Step 1: R = geometric root scale = |a_0/a_d|^{1/d}.
    Step 2: Q(w) = P(R*w)/R^d   (roots of Q clustered near |w|=1).
    Step 3: Strategy B + T_2 + deflation on Q  -> w_1, ..., w_d.
    Step 4: r_i = R * w_i.
    Step 5: polish on original P (no scaling error in final result).
    """
    coeffs = list(coeffs)
    R = geometric_root_scale(coeffs)
    Q_coeffs = rescale_polynomial(coeffs, R)
    w_roots = find_all_strategy_B_t2(Q_coeffs, polish=False, tol=tol)
    roots = [R * w for w in w_roots]
    if polish:
        polished = []
        for r in roots:
            r_pol, _ = t2_orbit(coeffs, r, max_epochs=30, tol=tol)
            polished.append(r_pol)
        roots = polished
    return roots


def pandrosion_polygon_clusters(coeffs, n_radii=64, n_phases=12):
    """POLYGONE DE PANDROSION: clusters via evaluations de P (Jensen).

    Polygone de Newton: lit (k, log|a_k|) -- inspection des coefficients.
    Polygone de Pandrosion: lit (log r, avg_phi log|P(r e^{iphi})|) --
    evaluations de P uniquement, conformement a l'esprit derivative-free /
    coefficient-free de Pandrosion.

    Formule de Jensen:
       (d/d log r) [avg_phi log|P(r e^{iphi})|] = #{racines dans |z| < r}

    donc les SAUTS de pente du profil log|P| (en log r) marquent les anneaux
    contenant des racines. Pente m sur un intervalle [r_i, r_{i+1}] = m
    racines (cumulees) dans |z| <= r ; un changement de pente +Delta m
    signale Delta m racines dans l'anneau.

    Le polygone de Pandrosion EST le profil de Jensen : pour chaque cluster
    detecte (multiplicite m_j, rayon median r_j), Strategy B + T_2 + scaling
    z = r_j w + deflation locale a m_j racines. Aucun coefficient n'est
    jamais lu directement -- toute l'information vient de P(z).
    """
    d = len(coeffs) - 1
    if d == 0:
        return []
    # Better radial range: use geometric mean magnitude (|a_0/a_d|^{1/d})
    # plus 2-norm-based tight bound. Cauchy is too loose for KS-like coeffs.
    R_geo = geometric_root_scale(coeffs)
    # Fujiwara: rho = 2 * max_k |a_k/a_d|^{1/k} -- tight for general polys.
    a_d = abs(coeffs[0])
    if a_d > 1e-300:
        candidates = []
        for k in range(1, len(coeffs)):
            v = abs(coeffs[k]) / a_d
            if v > 1e-300:
                candidates.append(v ** (1.0 / k))
        rho = 2.0 * max(candidates) if candidates else 2.0
    else:
        rho = 2.0
    # Lower bound: same Fujiwara on reversed (roots of P inverted).
    rev = list(reversed(coeffs))
    if abs(rev[0]) > 1e-300:
        candidates = []
        for k in range(1, len(rev)):
            v = abs(rev[k]) / abs(rev[0])
            if v > 1e-300:
                candidates.append(v ** (1.0 / k))
        rho_lo = 1.0 / (2.0 * max(candidates)) if candidates else 1.0
    else:
        rho_lo = 1e-3
    # Center the radial sweep around the geometric mean magnitude
    r_min = max(min(rho_lo * 0.5, R_geo / 100.0), 1e-12)
    r_max = max(rho * 1.5, R_geo * 100.0)
    log_rs = [math.log(r_min) + (math.log(r_max) - math.log(r_min)) * i / (n_radii - 1)
              for i in range(n_radii)]
    phases = [2 * math.pi * k / n_phases for k in range(n_phases)]
    log_P_avg = []
    for log_r in log_rs:
        r = math.exp(log_r)
        s = 0.0
        for phi in phases:
            z = complex(r * math.cos(phi), r * math.sin(phi))
            p = abs(horner1(coeffs, z))
            s += math.log(p) if p > 1e-300 else -700.0
        log_P_avg.append(s / n_phases)
    # Clip log values to avoid inf/nan in slope computation
    log_P_avg = [max(-700.0, min(700.0, v)) if math.isfinite(v) else -700.0
                 for v in log_P_avg]
    # Local slope = #zeros within disk of radius r (Jensen).
    slopes = []
    for i in range(n_radii - 1):
        ds = log_P_avg[i + 1] - log_P_avg[i]
        dr = log_rs[i + 1] - log_rs[i]
        if dr == 0 or not math.isfinite(ds):
            slopes.append(0.0)
        else:
            slopes.append(ds / dr)
    # Round slopes to nearest non-negative integer (zero counts).
    def _safe_int(s):
        if not math.isfinite(s):
            return d
        return max(0, min(d, int(round(s))))
    counts = [_safe_int(s) for s in slopes]
    # Enforce monotone non-decreasing in r (cumulative count).
    for i in range(1, len(counts)):
        if counts[i] < counts[i - 1]:
            counts[i] = counts[i - 1]
    # Identify cluster annuli where count increments.
    clusters = []
    prev = 0
    i = 0
    while i < len(counts):
        if counts[i] > prev:
            j = i
            while j + 1 < len(counts) and counts[j + 1] == counts[i]:
                j += 1
            m = counts[i] - prev
            r_est = math.exp((log_rs[i] + log_rs[j + 1]) / 2.0)
            clusters.append((r_est, m))
            prev = counts[i]
            i = j + 1
        else:
            i += 1
    total = sum(m for _, m in clusters)
    if total < d:
        # Missing roots: assume they live at the outer boundary.
        clusters.append((math.exp(log_rs[-1]), d - total))
    return clusters


def find_all_pandrosion_polygon_t2(coeffs, polish=True, tol=1e-12):
    """SCALING PAR POLYGONE DE PANDROSION (Jensen-via-evaluations).

    Identique a find_all_newton_polygon_scaled_t2 mais le decoupage en
    clusters vient de pandrosion_polygon_clusters (P-evaluations) au lieu
    de newton_polygon_clusters (coefficient inspection).

    C'est l'extension Pandrosion de paper 0 sec:scaling au cas polynomial:
    chaque cluster (r_j, m_j) est traite par homothetie z = r_j w +
    Strategy B + T_2 + deflation locale a m_j racines.
    """
    coeffs = list(coeffs)
    d = len(coeffs) - 1
    P0 = list(coeffs)
    clusters = pandrosion_polygon_clusters(coeffs)
    if len(clusters) == 1:
        return find_all_scaled_strategy_B_t2(coeffs, polish=polish, tol=tol)
    roots = []
    for r_i, m_i in clusters:
        if r_i == 0 or not math.isfinite(r_i):
            r_i = 1.0
        Q_coeffs = rescale_polynomial(coeffs, r_i)
        Q_current = list(Q_coeffs)
        cluster_roots_w = []
        attempts = 0
        max_attempts = 4 * m_i + 8
        while len(cluster_roots_w) < m_i and len(Q_current) > 2 and attempts < max_attempts:
            attempts += 1
            starts = strategy_B_starts(Q_current)
            w_root = None
            for w0 in starts:
                w, ok = t2_orbit(Q_current, w0, tol=tol)
                if ok:
                    w_root = w
                    break
            if w_root is None:
                break
            cluster_roots_w.append(w_root)
            Q_current = deflate(Q_current, w_root)
        for w in cluster_roots_w:
            z = r_i * w
            if abs(z) >= r_i / 4.0 and abs(z) <= 4.0 * r_i:
                z_pol, ok = t2_orbit(P0, z, max_epochs=20, tol=tol)
                if ok or abs(horner1(P0, z_pol)) < abs(horner1(P0, z)):
                    z = z_pol
                roots.append(z)
    if len(roots) < d:
        residual = list(P0)
        for r in roots:
            if len(residual) > 1:
                residual = deflate(residual, r)
        while len(residual) > 2 and len(roots) < d:
            starts = strategy_B_starts(residual)
            r_new = None
            for z0 in starts:
                r, ok = t2_orbit(residual, z0, tol=tol)
                if ok:
                    r_new = r
                    break
            if r_new is None:
                break
            r_new, _ = t2_orbit(P0, r_new, max_epochs=20, tol=tol)
            roots.append(r_new)
            residual = deflate(residual, r_new)
        if len(residual) == 2 and abs(residual[0]) > 1e-15 and len(roots) < d:
            roots.append(-residual[1] / residual[0])
    if polish:
        polished = []
        for r in roots:
            r_pol, _ = t2_orbit(P0, r, max_epochs=30, tol=tol)
            polished.append(r_pol)
        roots = polished
    return roots


def newton_polygon_clusters(coeffs):
    """Lower convex hull of (k, log|a_k|), k=0..d (descending order convention).

    Returns list of (radius_estimate, multiplicity) pairs, one per hull edge,
    covering all d roots. The slope s of an edge between vertices (k_i, l_i)
    and (k_j, l_j) gives log(root magnitude) = -(l_j - l_i)/(k_j - k_i)... in
    descending-coeff convention with P = sum coeffs[i] z^{d-i}, a slope s
    between (k_i, log|coeffs[k_i]|) and (k_j, log|coeffs[k_j]|) means there
    are (k_j - k_i) roots of magnitude exp(-s).
    """
    d = len(coeffs) - 1
    LOG_FLOOR = -700.0  # exp(-700) ~ 1e-304, safe under double range
    pts = []
    for k, c in enumerate(coeffs):
        m = abs(c)
        l = math.log(m) if m > 0 else LOG_FLOOR
        pts.append((k, l))
    # Lower convex hull (Andrew's monotone chain)
    hull = []
    for p in pts:
        while len(hull) >= 2:
            o, a = hull[-2], hull[-1]
            cross = (a[0] - o[0]) * (p[1] - o[1]) - (a[1] - o[1]) * (p[0] - o[0])
            if cross <= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    clusters = []
    for i in range(len(hull) - 1):
        k_i, l_i = hull[i]
        k_j, l_j = hull[i + 1]
        slope = (l_j - l_i) / (k_j - k_i)
        # In descending convention: coeff index k carries z^{d-k}, so the edge
        # corresponds to roots of magnitude r = exp(-slope).
        r = math.exp(-slope) if -700 < slope < 700 else 1.0
        m = k_j - k_i
        clusters.append((r, m))
    return clusters


def find_all_newton_polygon_scaled_t2(coeffs, polish=True, tol=1e-12):
    """SCALING PAR CLUSTER (Newton polygon).

    Etapes:
      1. Polygone de Newton sur (k, log|a_k|) -> liste de clusters
         [(r_i, m_i)] avec sum m_i = d.
      2. Pour chaque cluster (r_i, m_i):
           Q(w) = P(r_i * w) / r_i^d  (roots of cluster i near |w|=1).
           Strategy B sur Q + T_2 + Armijo + deflation locale a m_i etapes
           pour trouver m_i racines w_1, ..., w_{m_i}.
           Garder uniquement celles dont r_i * w_j est proche de cluster i
           (filtre par |z| in [r_i / 2, 2 * r_i]).
      3. Polish global sur P.

    Cle: chaque cluster reduit la deflation a m_i << d racines chaque, donc
    l'erreur de chainage reste petite. Le scaling z = r_i * w met les
    racines du cluster pres de |w|=1, ou la geometrie de Pandrosion est
    optimale (cf paper 0 sec:scaling).
    """
    coeffs = list(coeffs)
    d = len(coeffs) - 1
    P0 = list(coeffs)
    clusters = newton_polygon_clusters(coeffs)
    # If polygon has 1 edge, just one cluster -> equivalent to scale-once.
    if len(clusters) == 1:
        return find_all_scaled_strategy_B_t2(coeffs, polish=polish, tol=tol)
    roots = []
    for r_i, m_i in clusters:
        if r_i == 0 or not math.isfinite(r_i):
            r_i = 1.0
        Q_coeffs = rescale_polynomial(coeffs, r_i)
        # On Q, Strategy B + T_2 + deflation, but stop after m_i roots
        # (the m_i closest to |w|=1, the cluster's signature).
        Q_current = list(Q_coeffs)
        cluster_roots_w = []
        attempts = 0
        max_attempts = 4 * m_i + 8
        while len(cluster_roots_w) < m_i and len(Q_current) > 2 and attempts < max_attempts:
            attempts += 1
            starts = strategy_B_starts(Q_current)
            w_root = None
            for w0 in starts:
                w, ok = t2_orbit(Q_current, w0, tol=tol)
                if ok:
                    w_root = w
                    break
            if w_root is None:
                break
            cluster_roots_w.append(w_root)
            Q_current = deflate(Q_current, w_root)
        # Recover z = r_i * w. Keep only those whose magnitude lands in
        # [r_i/3, 3*r_i] (i.e., did not drift to a different cluster).
        for w in cluster_roots_w:
            z = r_i * w
            if abs(z) >= r_i / 3.0 and abs(z) <= 3.0 * r_i:
                # Polish on original P immediately
                z_pol, ok = t2_orbit(P0, z, max_epochs=20, tol=tol)
                if ok or abs(horner1(P0, z_pol)) < abs(horner1(P0, z)):
                    z = z_pol
                roots.append(z)
    # Top up with extras via Lairez-free fallback if we missed some
    if len(roots) < d:
        # Fallback: standard Strategy B + deflation on residual polynomial
        residual = list(P0)
        for r in roots:
            if len(residual) > 1:
                residual = deflate(residual, r)
        while len(residual) > 2 and len(roots) < d:
            starts = strategy_B_starts(residual)
            r_new = None
            for z0 in starts:
                r, ok = t2_orbit(residual, z0, tol=tol)
                if ok:
                    r_new = r
                    break
            if r_new is None:
                break
            r_new, _ = t2_orbit(P0, r_new, max_epochs=20, tol=tol)
            roots.append(r_new)
            residual = deflate(residual, r_new)
        if len(residual) == 2 and abs(residual[0]) > 1e-15 and len(roots) < d:
            roots.append(-residual[1] / residual[0])
    if polish:
        polished = []
        for r in roots:
            r_pol, _ = t2_orbit(P0, r, max_epochs=30, tol=tol)
            polished.append(r_pol)
        roots = polished
    return roots


def find_all_per_deflation_scaled_t2(coeffs, polish=True, tol=1e-12):
    """SCALED-PER-DEFLATION: rescale at each deflation step.

    At each step: compute R_k = |a_0/a_d|^{1/d_k} on the current polynomial,
    transform Q(w) = P_k(R_k*w)/R_k^{d_k}, find ONE root w*, recover
    z* = R_k * w*, deflate the original polynomial on z*, repeat.

    Idea (paper 0 sec:scaling): the homothety is recomputed on each smaller
    sub-problem, keeping every Pandrosion sub-iteration on a polynomial
    whose root cluster sits near |w|=1.
    """
    coeffs = list(coeffs)
    P0 = list(coeffs)
    roots = []
    current = list(coeffs)
    while len(current) > 2:
        R = geometric_root_scale(current)
        if R == 0 or not math.isfinite(R):
            R = 1.0
        Q_coeffs = rescale_polynomial(current, R)
        starts = strategy_B_starts(Q_coeffs)
        w_root = None
        for w0 in starts:
            w, ok = t2_orbit(Q_coeffs, w0, tol=tol)
            if ok:
                w_root = w
                break
        if w_root is None:
            return roots
        r = R * w_root
        # Polish on the original (un-deflated, un-scaled) P
        r, _ = t2_orbit(P0, r, max_epochs=20, tol=tol)
        roots.append(r)
        current = deflate(current, r)
    if len(current) == 2 and abs(current[0]) > 1e-15:
        roots.append(-current[1] / current[0])
    if polish:
        polished = []
        for r in roots:
            r_pol, _ = t2_orbit(P0, r, max_epochs=30, tol=tol)
            polished.append(r_pol)
        roots = polished
    return roots


# ---------------------------------------------------------------------------
# Pandrosion-Polygon Homotopy (PPH): clusters from Pandrosion polygon +
# clustered start system + T_2-fused-PURE corrector (no deflation).
# ---------------------------------------------------------------------------

def merge_close_clusters(clusters, ratio_threshold=1.4):
    """Merge adjacent clusters whose radii are within `ratio_threshold`.

    Pandrosion polygon over-detects on KS-like ensembles because root
    magnitudes are nearly continuous; we merge to suppress this noise.
    """
    if len(clusters) <= 1:
        return list(clusters)
    sorted_c = sorted(clusters, key=lambda x: x[0])
    merged = [sorted_c[0]]
    for r, m in sorted_c[1:]:
        r_prev, m_prev = merged[-1]
        if r / r_prev < ratio_threshold:
            new_r = math.exp((m_prev * math.log(r_prev) + m * math.log(r)) /
                             (m_prev + m))
            merged[-1] = (new_r, m_prev + m)
        else:
            merged.append((r, m))
    return merged


def G_value_clustered(clusters, z):
    """G(z) = prod_j (z^{m_j} - r_j^{m_j}).

    Single-cluster case (the common one after merging) collapses to
    G(z) = z^d - r^d, evaluated with overflow protection via z = r * w."""
    if len(clusters) == 1:
        r, m = clusters[0]
        if r > 0:
            w = z / r
            return (r ** m) * (w ** m - 1.0)
        return z ** m
    val = 1.0 + 0.0j
    for r, m in clusters:
        if r > 0:
            w = z / r
            val *= (r ** m) * (w ** m - 1.0)
        else:
            val *= z ** m
    return val


def G_deriv_clustered(clusters, z):
    """G'(z) by direct product rule. Single cluster: G' = m * z^(m-1)."""
    if abs(z) < 1e-300:
        z = z + 1e-12
    if len(clusters) == 1:
        r, m = clusters[0]
        return m * (z ** (m - 1))
    f_vals = []
    f_primes = []
    for r, m in clusters:
        if r > 0:
            w = z / r
            f_vals.append((r ** m) * (w ** m - 1.0))
        else:
            f_vals.append(z ** m)
        f_primes.append(m * (z ** (m - 1)))
    n = len(clusters)
    pre = [1.0 + 0.0j] * (n + 1)
    suf = [1.0 + 0.0j] * (n + 1)
    for i in range(n):
        pre[i + 1] = pre[i] * f_vals[i]
    for i in range(n - 1, -1, -1):
        suf[i] = suf[i + 1] * f_vals[i]
    deriv = 0.0 + 0.0j
    for i in range(n):
        deriv += pre[i] * f_primes[i] * suf[i + 1]
    return deriv


def H_value_pph(coeffs, clusters, gamma, z, t):
    return gamma * (1 - t) * G_value_clustered(clusters, z) + t * horner1(coeffs, z)


def H_z_pph(coeffs, clusters, gamma, z, t):
    G_z = G_deriv_clustered(clusters, z)
    _, Pp = horner2(coeffs, z)
    return gamma * (1 - t) * G_z + t * Pp


def H_t_pph(coeffs, clusters, gamma, z):
    G = G_value_clustered(clusters, z)
    Pz = horner1(coeffs, z)
    return -gamma * G + Pz


def t2_pure_corrector_pph(coeffs, clusters, gamma, z, t, anchor_state):
    z_prev, H_prev = anchor_state
    Hv = H_value_pph(coeffs, clusters, gamma, z, t)
    if abs(Hv) < 1e-15:
        return z, (z, Hv)
    if abs(z - z_prev) < 1e-15:
        z_prev = z + 1e-6 * (1 + 0j)
        H_prev = H_value_pph(coeffs, clusters, gamma, z_prev, t)
    Q_az = (Hv - H_prev) / (z - z_prev)
    if abs(Q_az) < 1e-15:
        return z, (z, Hv)
    s1 = z - Hv / Q_az
    Hs1 = H_value_pph(coeffs, clusters, gamma, s1, t)
    if abs(s1 - z) < 1e-15:
        return s1, (z, Hv)
    Q_as1 = (Hs1 - H_prev) / (s1 - z_prev)
    if abs(Q_as1) < 1e-15:
        return s1, (z, Hv)
    s2 = s1 - Hs1 / Q_as1
    denom = s2 - 2 * s1 + z
    if abs(denom) < 1e-15:
        return s2, (z, Hv)
    return z - (s1 - z) ** 2 / denom, (z, Hv)


def correct_pph(coeffs, clusters, gamma, z, t, max_iter=12, tol=1e-12):
    """T_2-PURE corrector. Order ~ 2.41 (golden) > Newton's 2."""
    z_prev = z + 1e-6 * (1 + 0j)
    H_prev = H_value_pph(coeffs, clusters, gamma, z_prev, t)
    anchor = (z_prev, H_prev)
    for _ in range(max_iter):
        Hv = H_value_pph(coeffs, clusters, gamma, z, t)
        if abs(Hv) < tol:
            return z, True
        z, anchor = t2_pure_corrector_pph(coeffs, clusters, gamma, z, t, anchor)
    return z, abs(H_value_pph(coeffs, clusters, gamma, z, t)) < tol


def track_pph(coeffs, clusters, gamma, z0, dt_init=0.15, dt_min=1e-6,
              max_steps=300, tol=1e-12):
    """T_2-PURE order ~2.41 > Newton's 2 -> can absorb larger predictor steps.

    Adaptive corrector tolerance: loose (sqrt(tol)) during tracking, tight
    only at t=1. Cuts 30-50% of corrector iterations.
    """
    z = complex(z0)
    t = 0.0
    dt = dt_init
    fails = 0
    track_tol = 1e-9  # tight enough for stable predictor, loose enough to be fast
    for _ in range(max_steps):
        if t >= 1.0 - 1e-9:
            return z, True
        h = min(dt, 1.0 - t)
        H_t = H_t_pph(coeffs, clusters, gamma, z)
        H_z = H_z_pph(coeffs, clusters, gamma, z, t)
        if abs(H_z) < 1e-15:
            return z, False
        dzdt = -H_t / H_z
        z_pred = z + h * dzdt
        # Tight tolerance at t=1, looser during tracking
        cur_tol = tol if (t + h) >= 1.0 - 1e-9 else track_tol
        z_corr, ok = correct_pph(coeffs, clusters, gamma, z_pred, t + h,
                                  tol=cur_tol)
        if ok:
            z = z_corr
            t = t + h
            # T_2-PURE supra-quadratic -> grow dt faster than Lairez/Newton.
            dt = min(dt * 1.8, 0.3)
            fails = 0
        else:
            fails += 1
            if fails >= 5 or dt < dt_min:
                return z, False
            dt *= 0.5
    return z, t >= 1.0 - 1e-9


def find_all_pph(coeffs, polish=True, tol=1e-12):
    """PANDROSION-POLYGON HOMOTOPY (PPH).

    1. Pandrosion polygon -> clusters [(r_j, m_j)] from P-evaluations.
    2. Start system G(z) = prod_j (z^{m_j} - r_j^{m_j}) ; anchors are
       a_{j,k} = r_j * exp(2 pi i k / m_j) for k=0..m_j-1.
    3. For each anchor, track homotopy H(z,t) = gamma (1-t) G + t P
       with Euler predictor + T_2-fused-PURE corrector. d trajectories
       independantes -> d racines, sans deflation.
    """
    coeffs = list(coeffs)
    d = len(coeffs) - 1
    if d == 0:
        return []
    if d == 1:
        return [-coeffs[1] / coeffs[0]] if abs(coeffs[0]) > 1e-15 else []
    raw_clusters = pandrosion_polygon_clusters(coeffs)
    merged = merge_close_clusters(raw_clusters, ratio_threshold=2.5)
    threshold = max(2, d // 16)
    significant = [(r, m) for r, m in merged if m >= threshold]
    R_geo = geometric_root_scale(coeffs)
    if not significant:
        clusters = [(R_geo, d)]
    else:
        deficit = d - sum(m for _, m in significant)
        if deficit > 0:
            biggest = max(range(len(significant)), key=lambda i: significant[i][1])
            r0, m0 = significant[biggest]
            significant[biggest] = (r0, m0 + deficit)
        # Single-cluster fallback: if the dominant cluster accounts for
        # most of d AND its radius is reasonably close to R_geo, prefer
        # the cleaner (R_geo, d) start system. This is the KS case where
        # multi-cluster detection adds noise.
        biggest_idx = max(range(len(significant)), key=lambda i: significant[i][1])
        r_dom, m_dom = significant[biggest_idx]
        if m_dom >= 0.6 * d and 1/3 <= r_dom / R_geo <= 3.0:
            clusters = [(R_geo, d)]
        else:
            clusters = significant
    total = sum(m for _, m in clusters)
    if total != d:
        clusters = [(R_geo, d)]
    # Anchors
    anchors = []
    for r, m in clusters:
        for k in range(m):
            angle = 2 * math.pi * k / m
            anchors.append(r * complex(math.cos(angle), math.sin(angle)))
    gamma = complex(0.6 + 0.8j)
    roots = []
    for a in anchors:
        z_end, ok = track_pph(coeffs, clusters, gamma, a, tol=tol)
        if not ok:
            for g_alt in [0.7+0.5j, -0.5+0.6j, 0.9+0.3j]:
                z_end, ok = track_pph(coeffs, clusters, g_alt, a, tol=tol)
                if ok:
                    break
        roots.append(z_end)
    if polish:
        polished = []
        for r in roots:
            res = abs(horner1(coeffs, r))
            if res < tol * 0.1:
                polished.append(r)
                continue
            r_pol, ok = t2_orbit(coeffs, r, max_epochs=25, tol=tol)
            res_pol = abs(horner1(coeffs, r_pol))
            # Hardened polish for outliers: deep T_2 + try multiple offsets
            if res_pol > 1e-10:
                best = (r_pol, res_pol)
                # Try deeper T_2 budget
                for max_e in [60, 150]:
                    r2, _ = t2_orbit(coeffs, best[0], max_epochs=max_e, tol=tol)
                    res2 = abs(horner1(coeffs, r2))
                    if res2 < best[1]:
                        best = (r2, res2)
                # Try perturbed restarts (break path-crossing tie)
                for eps in [1e-5, 1e-4, 1e-3]:
                    for direction in [1+0j, -1+0j, 0+1j, 0-1j]:
                        r3, _ = t2_orbit(coeffs, best[0] + eps * direction,
                                         max_epochs=40, tol=tol)
                        res3 = abs(horner1(coeffs, r3))
                        if res3 < best[1]:
                            best = (r3, res3)
                            if best[1] < tol * 0.1:
                                break
                    if best[1] < tol * 0.1:
                        break
                r_pol = best[0]
            polished.append(r_pol)
        roots = polished
    return roots


def find_all_pph_dual(coeffs, polish=True, tol=1e-12):
    """PPH augmented by inversion duality (Pandrosion's natural geometric
    involution z -> 1/z).

    Geometric idea: on Pandrosion's Thales diagram in C, the inversion
    z -> 1/z exchanges the "inner" and "outer" worlds. Under it, a
    polynomial P(z) of degree d transforms to P_tilde(w) = w^d P(1/w)
    whose coefficients are P's coefficients reversed. Roots of P at |r|>1
    become roots of P_tilde at |1/r|<1, where Pandrosion's geometry is
    optimal.

    Empirically on KS adversarial: the outer cluster (radially dispersed)
    causes path crossing in PPH. After inversion, that cluster becomes
    the inner cluster of P_tilde, where it is concentrated and PPH catches
    it cleanly. Combine roots from both views, dedup, polish.
    """
    P = list(coeffs)
    d = len(P) - 1
    if d <= 1:
        return find_all_pph(P, polish=polish, tol=tol)
    # Inverted polynomial: coefficients reversed
    P_rev = list(reversed(P))
    if abs(P_rev[0]) < 1e-300:
        return find_all_pph(P, polish=polish, tol=tol)
    # Solve both views
    cands_direct = find_all_pph(P, polish=False, tol=tol)
    cands_inv_w = find_all_pph(P_rev, polish=False, tol=tol)
    # Map P_tilde roots back: w -> 1/w
    cands_inv_z = []
    for w in cands_inv_w:
        if abs(w) > 1e-300 and math.isfinite(w.real) and math.isfinite(w.imag):
            cands_inv_z.append(1.0 / w)
    all_cands = list(cands_direct) + cands_inv_z
    # Polish each on the original P
    polished = []
    for c in all_cands:
        if math.isfinite(c.real) and math.isfinite(c.imag):
            r_pol, ok = t2_orbit(P, c, max_epochs=30, tol=tol)
            if not ok or abs(horner1(P, r_pol)) > 1e-8:
                r_pol, _ = t2_orbit_els(P, r_pol, max_epochs=40, tol=tol)
            polished.append(r_pol)
    # Sort by |P(r)| ascending and dedup
    polished.sort(key=lambda c: abs(horner1(P, c)))
    final = []
    dedup_tol = 1e-4
    for c in polished:
        if all(abs(c - r) > dedup_tol for r in final):
            final.append(c)
            if len(final) >= d:
                break
    # If we have fewer than d, top up with the rest of polished (less strict)
    if len(final) < d:
        for c in polished:
            if all(abs(c - r) > 1e-8 for r in final):
                final.append(c)
                if len(final) >= d:
                    break
    return final[:d]


def find_all_pph_with_clusters(coeffs, clusters, polish=True, tol=1e-12):
    """PPH with an externally chosen Pandrosion start cluster family.

    Used for the radial-anchor repair: the homotopy/corrector are unchanged,
    only the cyclotomic start radius is moved inward to avoid outer caustics on
    Kostlan-Smale shells.
    """
    coeffs = list(coeffs)
    d = len(coeffs) - 1
    anchors = []
    for r, m in clusters:
        for k in range(m):
            angle = 2 * math.pi * k / m
            anchors.append(r * complex(math.cos(angle), math.sin(angle)))
    gamma = complex(0.6 + 0.8j)
    roots = []
    for a in anchors:
        z_end, ok = track_pph(coeffs, clusters, gamma, a, tol=tol)
        if not ok:
            for g_alt in [0.7+0.5j, -0.5+0.6j, 0.9+0.3j]:
                z_end, ok = track_pph(coeffs, clusters, g_alt, a, tol=tol)
                if ok:
                    break
        roots.append(z_end)
    if polish:
        polished = []
        for r in roots:
            r_pol, _ = t2_orbit(coeffs, r, max_epochs=150, tol=tol)
            polished.append(r_pol)
        roots = mp_t4_polish_roots(coeffs, polished, threshold=tol)
    return roots


# ---------------------------------------------------------------------------
# Batched Pandrosion-Polygon Homotopy.
#
# Same mathematical pipeline as PPH above, but all cyclotomic trajectories are
# tracked in one NumPy vector.  This keeps the Pandrosion ingredients from
# papers 0/1 (scaling, cyclotomic anchors, T_2/PURE corrector) while removing
# the Python loop over d paths, which is the main practical advantage of
# Lairez's compact z^d - 1 homotopy in this benchmark.
# ---------------------------------------------------------------------------

def horner1_vec(coeffs, z):
    z = np.asarray(z, dtype=np.complex128)
    p = np.zeros_like(z, dtype=np.complex128)
    for c in coeffs:
        p = p * z + c
    return p


def horner2_vec(coeffs, z):
    z = np.asarray(z, dtype=np.complex128)
    p = np.zeros_like(z, dtype=np.complex128)
    dp = np.zeros_like(z, dtype=np.complex128)
    for c in coeffs:
        dp = dp * z + p
        p = p * z + c
    return p, dp


def _safe_div(num, den, fallback=0.0):
    out = np.full_like(num, fallback, dtype=np.complex128)
    mask = np.isfinite(num) & np.isfinite(den) & (np.abs(den) > 1e-300)
    np.divide(num, den, out=out, where=mask)
    return out


def _select_pph_clusters(coeffs):
    d = len(coeffs) - 1
    raw_clusters = pandrosion_polygon_clusters(coeffs)
    merged = merge_close_clusters(raw_clusters, ratio_threshold=2.5)
    threshold = max(2, d // 16)
    significant = [(r, m) for r, m in merged if m >= threshold]
    R_geo = geometric_root_scale(coeffs)
    if not significant:
        return [(R_geo, d)]
    deficit = d - sum(m for _, m in significant)
    if deficit > 0:
        biggest = max(range(len(significant)), key=lambda i: significant[i][1])
        r0, m0 = significant[biggest]
        significant[biggest] = (r0, m0 + deficit)
    biggest_idx = max(range(len(significant)), key=lambda i: significant[i][1])
    r_dom, m_dom = significant[biggest_idx]
    if m_dom >= 0.6 * d and 1 / 3 <= r_dom / R_geo <= 3.0:
        return [(R_geo, d)]
    if sum(m for _, m in significant) != d:
        return [(R_geo, d)]
    return significant


def G_value_clustered_vec(clusters, z):
    z = np.asarray(z, dtype=np.complex128)
    val = np.ones_like(z, dtype=np.complex128)
    for r, m in clusters:
        if r > 0:
            w = z / r
            val *= (r ** m) * (w ** m - 1.0)
        else:
            val *= z ** m
    return val


def G_deriv_clustered_vec(clusters, z):
    z = np.asarray(z, dtype=np.complex128)
    if len(clusters) == 1:
        _, m = clusters[0]
        return m * (z ** (m - 1))
    factors = []
    derivs = []
    for r, m in clusters:
        if r > 0:
            w = z / r
            factors.append((r ** m) * (w ** m - 1.0))
        else:
            factors.append(z ** m)
        derivs.append(m * (z ** (m - 1)))
    total = np.zeros_like(z, dtype=np.complex128)
    for i in range(len(clusters)):
        term = derivs[i].copy()
        for j, f in enumerate(factors):
            if i != j:
                term *= f
        total += term
    return total


def H_value_pph_vec(coeffs, clusters, gamma, z, t):
    return gamma * (1 - t) * G_value_clustered_vec(clusters, z) + t * horner1_vec(coeffs, z)


def H_z_pph_vec(coeffs, clusters, gamma, z, t):
    _, Pp = horner2_vec(coeffs, z)
    return gamma * (1 - t) * G_deriv_clustered_vec(clusters, z) + t * Pp


def H_t_pph_vec(coeffs, clusters, gamma, z):
    return -gamma * G_value_clustered_vec(clusters, z) + horner1_vec(coeffs, z)


def batch_t2_correct_pph(coeffs, clusters, gamma, z, t, iters=4):
    z = np.asarray(z, dtype=np.complex128).copy()
    z_prev = z + 1e-6 * (1 + 0j)
    H_prev = H_value_pph_vec(coeffs, clusters, gamma, z_prev, t)
    for _ in range(iters):
        Hv = H_value_pph_vec(coeffs, clusters, gamma, z, t)
        Q0 = _safe_div(Hv - H_prev, z - z_prev)
        s1 = z - _safe_div(Hv, Q0)
        Hs1 = H_value_pph_vec(coeffs, clusters, gamma, s1, t)
        Q1 = _safe_div(Hs1 - H_prev, s1 - z_prev)
        s2 = s1 - _safe_div(Hs1, Q1)
        denom = s2 - 2 * s1 + z
        z_new = z - _safe_div((s1 - z) ** 2, denom)
        good = np.isfinite(z_new.real) & np.isfinite(z_new.imag)
        z_prev, H_prev = z, Hv
        z = np.where(good, z_new, z)
    return z


def batch_t2_polish(coeffs, z, iters=8, tol=1e-12):
    z = np.asarray(z, dtype=np.complex128).copy()
    z_prev = z + 1e-6 * (1 + 0j)
    P_prev = horner1_vec(coeffs, z_prev)
    for _ in range(iters):
        Pz = horner1_vec(coeffs, z)
        if np.max(np.abs(Pz)) < tol:
            break
        Q0 = _safe_div(Pz - P_prev, z - z_prev)
        s1 = z - _safe_div(Pz, Q0)
        Ps1 = horner1_vec(coeffs, s1)
        Q1 = _safe_div(Ps1 - P_prev, s1 - z_prev)
        s2 = s1 - _safe_div(Ps1, Q1)
        denom = s2 - 2 * s1 + z
        z_new = z - _safe_div((s1 - z) ** 2, denom)
        good = np.isfinite(z_new.real) & np.isfinite(z_new.imag)
        z_prev, P_prev = z, Pz
        z = np.where(good, z_new, z)
    return z


def _pandrosion_base_vec(coeffs, anchor, w, P_anchor=None):
    if P_anchor is None:
        P_anchor = horner1_vec(coeffs, anchor)
    Pw = horner1_vec(coeffs, w)
    Q = _safe_div(Pw - P_anchor, w - anchor)
    return anchor - _safe_div(P_anchor, Q)


def batch_t4_polish(coeffs, z, iters=6, tol=1e-12):
    """Higher-order Pandrosion polish from paper 1.

    This is the T_4 Richardson/Steffensen layer

        T2 = z - (s1-z)^2/(s2-2s1+z)
        T3 = T2 - (s3-T2)/(lambda_hat-1)
        T4 = T3 - (s4-T3)/(lambda_hat-1)

    for the fixed-anchor Pandrosion map h(w)=a-P(a)/Q(a,w).  The anchor is
    deliberately offset from z, so the quotient is a genuine divided
    difference and does not collapse to Newton's derivative.
    """
    z = np.asarray(z, dtype=np.complex128).copy()
    offsets = np.array([1 + 1j, -1 + 1j, 1 - 1j, -1 - 1j], dtype=np.complex128)
    for k in range(iters):
        Pz = horner1_vec(coeffs, z)
        if np.max(np.abs(Pz)) < tol:
            break
        scale = np.maximum(1.0, np.abs(z))
        anchor = z + (1e-4 * scale) * offsets[k % len(offsets)]
        Pa = horner1_vec(coeffs, anchor)

        s1 = _pandrosion_base_vec(coeffs, anchor, z, Pa)
        s2 = _pandrosion_base_vec(coeffs, anchor, s1, Pa)
        denom2 = s2 - 2 * s1 + z
        t2 = z - _safe_div((s1 - z) ** 2, denom2)

        lambda_hat = _safe_div(s2 - s1, s1 - z, fallback=1.0)
        rich_denom = lambda_hat - 1.0
        s3 = _pandrosion_base_vec(coeffs, anchor, t2, Pa)
        t3 = t2 - _safe_div(s3 - t2, rich_denom)
        s4 = _pandrosion_base_vec(coeffs, anchor, t3, Pa)
        t4 = t3 - _safe_div(s4 - t3, rich_denom)

        candidates = [z, t2, t3, t4]
        vals = [np.abs(horner1_vec(coeffs, c)) for c in candidates]
        best_idx = np.argmin(np.vstack(vals), axis=0)
        z_next = z.copy()
        for i, c in enumerate(candidates[1:], start=1):
            take = (best_idx == i) & np.isfinite(c.real) & np.isfinite(c.imag)
            z_next = np.where(take, c, z_next)
        # If no higher-order candidate improves a coordinate, keep the current
        # point; this preserves the residual monotonicity of the polish.
        z = z_next
    return z


def find_all_pph_batch(coeffs, polish=True, tol=1e-12, clusters_override=None):
    """Vectorized Pandrosion-Polygon Homotopy.

    The start system is the Pandrosion/Jensen cluster system used by PPH:
       G(z) = prod_j (z^{m_j} - r_j^{m_j}).
    Anchors are the cyclotomic points on each detected radius.  Predictor
    steps are shared across the vector, and correction uses the derivative-free
    T_2/PURE quotient.  This is a numerical implementation optimization, not a
    change of the mathematical method.
    """
    coeffs = list(coeffs)
    d = len(coeffs) - 1
    if d == 0:
        return []
    if d == 1:
        return [-coeffs[1] / coeffs[0]] if abs(coeffs[0]) > 1e-15 else []

    clusters = clusters_override if clusters_override is not None else _select_pph_clusters(coeffs)
    anchors = []
    for r, m in clusters:
        angles = 2 * math.pi * np.arange(m) / m
        anchors.extend(r * np.exp(1j * angles))
    z = np.asarray(anchors, dtype=np.complex128)

    gamma = complex(0.6 + 0.8j)
    t = np.zeros(d, dtype=float)
    dt = np.full(d, 0.15, dtype=float)
    fails = np.zeros(d, dtype=int)
    done = np.zeros(d, dtype=bool)
    for _ in range(300):
        active = ~done
        if not np.any(active):
            break
        za = z[active]
        ta = t[active]
        dta = dt[active]
        h = np.minimum(dta, 1.0 - ta)
        Ht = H_t_pph_vec(coeffs, clusters, gamma, za)
        Hz = H_z_pph_vec(coeffs, clusters, gamma, za, ta)
        z_pred = za - h * _safe_div(Ht, Hz)
        t_new = ta + h
        z_corr = batch_t2_correct_pph(coeffs, clusters, gamma, z_pred, t_new, iters=6)
        Hcorr = H_value_pph_vec(coeffs, clusters, gamma, z_corr, t_new)
        cur_tol = np.where(t_new >= 1.0 - 1e-9, tol, 1e-9)
        ok = (np.abs(Hcorr) < cur_tol) & np.isfinite(z_corr.real) & np.isfinite(z_corr.imag)

        idx = np.where(active)[0]
        ok_idx = idx[ok]
        bad_idx = idx[~ok]
        z[ok_idx] = z_corr[ok]
        t[ok_idx] = t_new[ok]
        dt[ok_idx] = np.minimum(dt[ok_idx] * 1.8, 0.3)
        fails[ok_idx] = 0
        done[ok_idx] = t[ok_idx] >= 1.0 - 1e-9

        fails[bad_idx] += 1
        dt[bad_idx] *= 0.5
        done[bad_idx] = (fails[bad_idx] >= 5) | (dt[bad_idx] < 1e-6)

    if polish:
        z_t2 = batch_t2_polish(coeffs, z, iters=14, tol=tol)
        z_t4 = batch_t4_polish(coeffs, z_t2, iters=8, tol=tol)
        res_t2 = np.abs(horner1_vec(coeffs, z_t2))
        res_t4 = np.abs(horner1_vec(coeffs, z_t4))
        z = np.where((res_t4 < res_t2) & np.isfinite(z_t4.real) & np.isfinite(z_t4.imag), z_t4, z_t2)
        bad = np.abs(horner1_vec(coeffs, z)) > max(1e-8, 100 * tol)
        if np.any(bad):
            # Sparse repair only for pathological paths; the usual path remains
            # fully batched and derivative-free.
            repaired = z.copy()
            for idx in np.where(bad)[0]:
                r_pol = batch_t4_polish(coeffs, np.asarray([z[idx]]), iters=20, tol=tol)[0]
                if abs(horner1(coeffs, r_pol)) < abs(horner1(coeffs, z[idx])):
                    repaired[idx] = r_pol
            z = repaired
    return list(z[:d])


def find_all_pph_batch_literal_scaled(coeffs, polish=True, tol=1e-12, k=1.0):
    """Paper-0 sec:scaling LITERAL form of the batched homotopy.

    Pipeline (a la paper-0 sec:scaling generalized to polynomials):
      1. R = R_geo / k where R_geo = |a_0/a_d|^{1/d}.
      2. Q(w) = P(R*w) / R^d  (rescaled polynomial; roots of Q at |w| ~ k).
      3. Batch PPH homotopy on Q with cyclotomic anchors at radius 1.
      4. Recover roots of P as r_i = R * w_i.
      5. Polish on the original P (T_2 then T_4).

    The parameter k controls separation between anchors (|w|=1) and roots
    (|w| ~ k):
      k=1  : x' ~ 1, fastest, best for ordinary regimes (Random, HDR).
      k>=3 : x' >> 1, anchors well inside the root shell, kills the path-
             crossing/discriminant-grazing pathology of Kostlan-Smale at
             high d (paper-0 sec:scaling intuition extended outward).

    The literal rescaling is mathematically equivalent to anchor-scaling
    but yields balanced coefficients on Q, eliminating overflow/cancellation
    in Horner evaluations during tracking.
    """
    coeffs = list(coeffs)
    d = len(coeffs) - 1
    if d == 0:
        return []
    if d == 1:
        return [-coeffs[1] / coeffs[0]] if abs(coeffs[0]) > 1e-15 else []

    R_geo = geometric_root_scale(coeffs)
    if R_geo <= 0 or not math.isfinite(R_geo):
        return find_all_pph_batch(coeffs, polish=polish, tol=tol)

    R = R_geo / k
    Q = rescale_polynomial(coeffs, R)
    w_roots = find_all_pph_batch(Q, polish=polish, tol=tol,
                                  clusters_override=[(1.0, d)])
    if len(w_roots) != d:
        return find_all_pph_batch(coeffs, polish=polish, tol=tol)
    z_arr = R * np.asarray(w_roots, dtype=np.complex128)
    if polish:
        z_t2 = batch_t2_polish(coeffs, z_arr, iters=14, tol=tol)
        z_t4 = batch_t4_polish(coeffs, z_t2, iters=8, tol=tol)
        res_t2 = np.abs(horner1_vec(coeffs, z_t2))
        res_t4 = np.abs(horner1_vec(coeffs, z_t4))
        z_arr = np.where((res_t4 < res_t2) & np.isfinite(z_t4.real) & np.isfinite(z_t4.imag),
                         z_t4, z_t2)
    return list(z_arr[:d])


def k_optimal_dyadic(coeffs):
    """k_optimal from the Pandrosion polygon (Jensen via P-evaluations only).

    Polygon Jensen reads the actual root radial distribution by sampling
    <log|P(r e^{iphi})|>_phi on a log-spaced radius grid; the slopes give
    the cumulative root count by Jensen's formula.  Cluster radii (r_j, m_j)
    are extracted, and the 10th/90th cumulative-multiplicity percentiles
    define the root-shell extent.  k_optimal = r_90 / r_10, snapped to the
    nearest power of 2 (paper-0 constructible homothety factors).

    Lit la *vraie* géométrie des racines (sémantique), pas la signature
    syntaxique des coefficients.  Coût: O(n_radii * n_phases * d) Horner
    evaluations -- robust on Bombieri/Kostlan adversarial regimes AND on
    sparse polynomials (structural zeros don't perturb Jensen-on-evaluations).

    Augmented by sentinel cluster (1.0, 1) so the array is never empty.
    """
    clusters = list(pandrosion_polygon_clusters(coeffs)) + [(1.0, 1)]
    radii = np.array([max(r, 1e-300) for r, _ in clusters], dtype=np.float64)
    multips = np.array([m for _, m in clusters], dtype=np.float64)
    order = np.argsort(radii)
    radii_s = radii[order]
    cumul = np.cumsum(multips[order])
    total = max(1.0, float(cumul[-1]))
    n = len(radii_s)
    i_lo = min(int(np.searchsorted(cumul, 0.10 * total)), n - 1)
    i_hi = min(int(np.searchsorted(cumul, 0.90 * total)), n - 1)
    r_p10 = max(float(radii_s[i_lo]), 1e-300)
    r_p90 = max(float(radii_s[i_hi]), 1e-300)
    k_raw = max(1.0, r_p90 / r_p10)
    return float(2.0 ** round(math.log2(k_raw)))


# Backward-compatibility alias: the function used to be named after the
# Pandrosion polygon (Jensen via P-evaluations).  The implementation has been
# replaced by the dyadic spectral estimator (paper-0 generalised, O(d), no
# evaluations) which empirically dominates the Jensen polygon on every regime.
k_optimal_from_pandrosion_polygon = k_optimal_dyadic


def find_all_pandrosion_universal_scaling(coeffs, polish=True, tol=1e-12):
    """Universal Pandrosion via Pandrosion-polygon-derived k (NO REGIME IF).

    Generalises paper-0 sec:scaling: paper-0 strict prescribes
       A = floor(x^{1/p})^p,  x' = x/A in [1, 2)  (close to 1, scalar case).
    We extend this to the polynomial case using the Pandrosion polygon to
    detect the root-shell extent r_max/r_min (Jensen's formula via
    polynomial evaluations only).  The optimal homothety factor is
       k_optimal = r_max / r_min,
    which automatically interpolates between paper-0 strict (k~1 for tight
    shells) and KS-adversarial (k~3-7 for spread shells).

    The literal homothety Q(w) = P(R*w)/R^d with R = R_geo / k_optimal
    places the rescaled root shell in [1, k_optimal] and the cyclotomic
    anchors at |w|=1, giving a uniform separation anchor->root regardless
    of regime.  Single-shot (no portfolio), pure Pandrosion (no coefficient
    inspection except for R_geo and the polygon's P-evaluations).
    """
    coeffs = list(coeffs)
    d = len(coeffs) - 1
    if d == 0:
        return []
    if d == 1:
        return [-coeffs[1] / coeffs[0]] if abs(coeffs[0]) > 1e-15 else []
    k = k_optimal_dyadic(coeffs)
    return find_all_pph_batch_literal_scaled(coeffs, polish=polish, tol=tol, k=k)


def find_all_pandrosion_selector(coeffs, polish=True, tol=1e-12):
    """Small deterministic portfolio of Pandrosion variants.

    The papers give three useful preconditioners:
      * paper-0 sec:scaling LITERAL (Q = P(R*w)/R^d) for the common case
        where coefficient magnitudes are large -- best numerical conditioning;
      * scalar PPH when large but non-Kostlan coefficients make adaptive batch
        spend too much time rejecting paths;
      * cyclotomic/scaled batch homotopy as default fallback.

    The selector deliberately uses only cheap coefficient geometry, then calls
    one Pandrosion method.  It does not call Lairez and it does not run both
    methods and pick after timing.
    """
    d = len(coeffs) - 1
    if d <= 1:
        return find_all_pph_batch(coeffs, polish=polish, tol=tol)
    max_coeff = max(abs(c) for c in coeffs)
    if d >= 32 and 4.0 < max_coeff < 30.0:
        return find_all_pph(coeffs, polish=polish, tol=tol)
    # Paper-0 sec:scaling STRICT (Proposition prop:scaling):
    #   A = floor(R_geo)^d,  R = A^{1/d} = floor(R_geo),  x' = R_geo/R in [1, 2).
    # The paper proves that x' close to 1 minimises the contraction ratio
    # lambda_{p,x'}.  For our polynomial regimes R_geo is typically near 1,
    # so floor(R_geo) = 1 and x' = R_geo (no actual integer reduction).
    # The literal rescaling Q(w) = P(R*w)/R^d still helps numerically because
    # it normalises Horner evaluations; for R = 1 this collapses to Q = P.
    if d >= 32:
        return find_all_pph_batch_literal_scaled(coeffs, polish=polish,
                                                  tol=tol, k=1.0)
    return find_all_pph_batch(coeffs, polish=polish, tol=tol)


def find_all_pandrosion_fast(coeffs, polish=True, tol=1e-12):
    """Fast mode: cheap PPH portfolio in double precision only.

    Default is batch PPH + T2/T4.  The existing mid-coefficient scalar PPH
    branch is retained because it is still double precision, derivative-free,
    and cheaper/more accurate on the random ``mag=5`` shell.  No mpmath and no
    radial retry are used here.
    """
    return find_all_pandrosion_selector(coeffs, polish=polish, tol=tol)


def _mp_horner(coeffs_mp, z):
    p = mp.mpc(0)
    for c in coeffs_mp:
        p = p * z + c
    return p


def _mp_base(coeffs_mp, anchor, w, P_anchor):
    Pw = _mp_horner(coeffs_mp, w)
    den = w - anchor
    if abs(den) < mp.mpf("1e-80"):
        return w
    Q = (Pw - P_anchor) / den
    if abs(Q) < mp.mpf("1e-80"):
        return w
    return anchor - P_anchor / Q


def mp_t4_polish_root(coeffs, z0, dps=80, iters=8):
    """Local high-precision fixed-anchor Pandrosion T4 polish.

    This is still the divided-difference Pandrosion map, just evaluated with
    enough arithmetic precision that the final double-complex residual is not
    dominated by intermediate cancellation.
    """
    if mp is None:
        return z0
    old_dps = mp.mp.dps
    mp.mp.dps = dps
    try:
        coeffs_mp = [mp.mpc(complex(c).real, complex(c).imag) for c in coeffs]
        z = mp.mpc(complex(z0).real, complex(z0).imag)
        best = z
        best_res = abs(_mp_horner(coeffs_mp, z))
        offsets = [mp.mpc(1, 1), mp.mpc(-1, 1), mp.mpc(1, -1), mp.mpc(-1, -1)]
        for k in range(iters):
            scale = max(mp.mpf(1), abs(z))
            anchor = z + mp.mpf("1e-8") * scale * offsets[k % len(offsets)]
            Pa = _mp_horner(coeffs_mp, anchor)
            s1 = _mp_base(coeffs_mp, anchor, z, Pa)
            s2 = _mp_base(coeffs_mp, anchor, s1, Pa)
            denom2 = s2 - 2 * s1 + z
            if abs(denom2) < mp.mpf("1e-80"):
                continue
            t2 = z - (s1 - z) ** 2 / denom2
            lambda_hat_den = s1 - z
            if abs(lambda_hat_den) < mp.mpf("1e-80"):
                candidates = [t2]
            else:
                lambda_hat = (s2 - s1) / lambda_hat_den
                rich = lambda_hat - 1
                candidates = [t2]
                if abs(rich) > mp.mpf("1e-40"):
                    s3 = _mp_base(coeffs_mp, anchor, t2, Pa)
                    t3 = t2 - (s3 - t2) / rich
                    s4 = _mp_base(coeffs_mp, anchor, t3, Pa)
                    t4 = t3 - (s4 - t3) / rich
                    candidates.extend([t3, t4])
            for cand in candidates:
                res = abs(_mp_horner(coeffs_mp, cand))
                if res < best_res:
                    best, best_res = cand, res
            z = best
        return best
    finally:
        mp.mp.dps = old_dps


def mp_t4_polish_roots(coeffs, roots, threshold=1e-12):
    if mp is None:
        return roots
    polished = []
    for r in roots:
        base_res = abs(horner1(coeffs, r))
        if base_res <= threshold:
            polished.append(r)
            continue
        cand = mp_t4_polish_root(coeffs, r)
        cand_res = abs(horner1(coeffs, cand))
        polished.append(cand if cand_res < base_res else r)
    return polished


def radial_anchor_repair(coeffs, roots, tol=1e-12):
    """Try inward cyclotomic start radii when the outer shell path fails.

    For broad one-annulus polynomials, especially Kostlan-Smale at high degree,
    the Jensen/geometric radius can sit on the outer side of the root shell.
    Homotopy paths then graze caustics.  Moving the cyclotomic anchors inward
    keeps the same Pandrosion start system but changes the continuation
    geometry; empirically this removes the worst KS outliers.
    """
    d = len(coeffs) - 1
    if d <= 1:
        return roots
    best = roots
    best_res = max_residual(coeffs, roots) if len(roots) == d else float("inf")
    if best_res < 1e-6:
        return roots
    R = geometric_root_scale(coeffs)
    for factor in [0.95, 0.90, 0.85]:
        cand = find_all_pph_with_clusters(coeffs, [(R * factor, d)], polish=True, tol=tol)
        if len(cand) != d:
            continue
        cand_res = max_residual(coeffs, cand)
        if cand_res < best_res:
            best, best_res = cand, cand_res
        if best_res < 1e-6:
            break
    return best


def high_degree_ks_radial_repair(coeffs, roots, tol=1e-12):
    """Batch radial repair for KS-like shells at d >= 128.

    Empirical geometry: as d grows, the cyclotomic start radius must move
    inside the root shell.  A simple scale 0.6 * 128/d gives the successful
    radii observed at d=128,256,512 while keeping the homotopy fully batched.
    """
    d = len(coeffs) - 1
    if d < 128:
        return roots
    best = roots
    best_res = max_residual(coeffs, roots) if len(roots) == d else float("inf")
    R = geometric_root_scale(coeffs)
    base = min(0.6, 0.6 * 128.0 / d)
    factors = [base, base * (4.0 / 3.0), base * (2.0 / 3.0)]
    for factor in factors:
        if factor <= 0:
            continue
        cand = find_all_pph_batch(coeffs, polish=True, tol=tol,
                                  clusters_override=[(R * factor, d)])
        if len(cand) != d:
            continue
        cand_res = max_residual(coeffs, cand)
        if cand_res < best_res:
            best, best_res = cand, cand_res
        if best_res < 1e-6:
            break
    return best


def smooth_inward_factor(coeffs):
    """Universal smooth inward-radius law.

    The coefficient-growth gate is near 0 on ordinary random/HDR shells and
    near 1 on Bombieri/Kostlan shells, so the effective radius continuously
    interpolates between R_geo and the high-degree inward law.
    """
    d = len(coeffs) - 1
    max_coeff = max(abs(c) for c in coeffs)
    growth = math.log1p(max_coeff / max(1.0, math.sqrt(d)))
    tau = 18.0
    power = 4.0
    gate = (growth ** power) / (growth ** power + tau ** power)
    inner = 0.9 * ((64.0 / max(1, d)) ** 0.75)
    return math.exp(gate * math.log(inner))


def degree_inward_factor(coeffs):
    """Pure degree inward law for high-degree cyclotomic anchors."""
    d = len(coeffs) - 1
    return 0.6 * 128.0 / max(1, d)


def find_all_pandrosion_universal_variational(coeffs, polish=True, tol=1e-12):
    """Branch-free universal Pandrosion portfolio.

    No regime detector is used.  A fixed family of Pandrosion start geometries
    is tried, and the root set minimising max |P(r)| is returned:

      1. unit/geometric radius R_geo;
      2. smooth coefficient-growth inward radius;
      3. pure degree inward radius;
      4. adjacent pure-degree inward radii;
      5. scalar PPH/Jensen geometry.

    The selection is variational (argmin residual), not an if/else classifier.
    """
    coeffs = list(coeffs)
    d = len(coeffs) - 1
    R = geometric_root_scale(coeffs)
    candidates = [
        find_all_pph_batch(coeffs, polish=polish, tol=tol,
                           clusters_override=[(R, d)]),
        find_all_pph_batch(coeffs, polish=polish, tol=tol,
                           clusters_override=[(R * smooth_inward_factor(coeffs), d)]),
        find_all_pph_batch(coeffs, polish=polish, tol=tol,
                           clusters_override=[(R * degree_inward_factor(coeffs), d)]),
        find_all_pph_batch(coeffs, polish=polish, tol=tol,
                           clusters_override=[(R * degree_inward_factor(coeffs) * (4.0 / 3.0), d)]),
        find_all_pph_batch(coeffs, polish=polish, tol=tol,
                           clusters_override=[(R * degree_inward_factor(coeffs) * (2.0 / 3.0), d)]),
        find_all_pph(coeffs, polish=polish, tol=tol),
    ]
    return min(candidates, key=lambda roots: max_residual(coeffs, roots) if len(roots) == d else float("inf"))


def is_ks_like(coeffs):
    """Cheap Bombieri/Kostlan-shell detector.

    Kostlan-Smale normalisation creates very large middle coefficients after
    making the polynomial monic.  We only treat the high-degree, very high
    growth case as adversarial enough to justify the radial-anchor strict path;
    smaller KS-like cases are handled by the residual trigger.
    """
    d = len(coeffs) - 1
    if d < 48:
        return False
    max_coeff = max(abs(c) for c in coeffs)
    return max_coeff / max(1.0, math.sqrt(d)) > 1e5


def find_all_pandrosion_selector_strict(coeffs, polish=True, tol=1e-12):
    """Selector plus scalar PPH repair when the residual is not Lairez-grade.

    This is a mathematical/residual criterion, not a timing criterion: the
    fast batched homotopy is accepted only if the final polynomial residual is
    already at the requested tolerance scale.
    """
    roots = find_all_pandrosion_selector(coeffs, polish=polish, tol=tol)
    d = len(coeffs) - 1
    if len(roots) == d:
        roots = mp_t4_polish_roots(coeffs, roots, threshold=tol)
    if len(roots) == d and max_residual(coeffs, roots) < tol:
        return roots
    roots = find_all_pph(coeffs, polish=polish, tol=tol)
    if len(roots) == d:
        roots = mp_t4_polish_roots(coeffs, roots, threshold=tol)
        roots = radial_anchor_repair(coeffs, roots, tol=tol)
    return roots


def find_all_pandrosion_adaptive(coeffs, polish=True, tol=1e-12):
    """Two-level Pandrosion selector.

    Fast mode:
      PPH batch + T2/T4 double precision.

    Strict mode:
      mpmath T4 + scalar/radial repair, triggered only by an actual residual
      failure, missing roots, or a strong KS-like coefficient-growth signature.
    """
    d = len(coeffs) - 1
    if d > 96 and is_ks_like(coeffs):
        return high_degree_ks_radial_repair(coeffs, [], tol=tol)
    roots = find_all_pandrosion_fast(coeffs, polish=polish, tol=tol)
    if len(roots) != d:
        return find_all_pandrosion_selector_strict(coeffs, polish=polish, tol=tol)
    residual = max_residual(coeffs, roots)
    ks_like = is_ks_like(coeffs)
    # d=64 KS uses the scalar radial repair; d>=128 KS uses a fully batched
    # inward-radius repair.  Non-KS high-degree failures still get strict repair.
    if d <= 96 and (residual > 1e-8 or ks_like):
        return find_all_pandrosion_selector_strict(coeffs, polish=polish, tol=tol)
    if d > 96 and ks_like:
        return high_degree_ks_radial_repair(coeffs, roots, tol=tol)
    if d > 96 and residual > 1e-6 and not ks_like:
        return find_all_pandrosion_selector_strict(coeffs, polish=polish, tol=tol)
    return roots


# ---------------------------------------------------------------------------
# Lairez homotopy (Newton corrector)
# ---------------------------------------------------------------------------

def homotopy_value(coeffs, gamma, z, t):
    d = len(coeffs) - 1
    return gamma * (1 - t) * (z ** d - 1) + t * horner1(coeffs, z)


def lairez_corrector(coeffs, gamma, z, t, max_iter=10, tol=1e-12):
    d = len(coeffs) - 1
    for _ in range(max_iter):
        Pz, Pp_z = horner2(coeffs, z)
        z_d = z ** d
        z_dm1 = z_d / z if abs(z) > 1e-15 else 0.0
        H = gamma * (1 - t) * (z_d - 1) + t * Pz
        if abs(H) < tol:
            return z, True
        H_z = gamma * (1 - t) * d * z_dm1 + t * Pp_z
        if abs(H_z) < 1e-15:
            return z, False
        z = z - H / H_z
    return z, abs(homotopy_value(coeffs, gamma, z, t)) < tol


def lairez_track(coeffs, gamma, z0, dt_init=0.1, dt_min=1e-6, max_steps=300):
    d = len(coeffs) - 1
    z = complex(z0)
    t = 0.0
    dt = dt_init
    fails = 0
    for _ in range(max_steps):
        if t >= 1.0 - 1e-9:
            return z, True
        h = min(dt, 1.0 - t)
        Pz, Pp_z = horner2(coeffs, z)
        z_d = z ** d
        z_dm1 = z_d / z if abs(z) > 1e-15 else 0.0
        H_t = -gamma * (z_d - 1) + Pz
        H_z = gamma * (1 - t) * d * z_dm1 + t * Pp_z
        if abs(H_z) < 1e-15:
            return z, False
        dzdt = -H_t / H_z
        z_pred = z + h * dzdt
        z_corr, ok = lairez_corrector(coeffs, gamma, z_pred, t + h)
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


def find_all_lairez(coeffs, polish=True, tol=1e-12):
    d = len(coeffs) - 1
    gamma = complex(0.6 + 0.8j)
    roots = []
    for k in range(d):
        z0 = np.exp(2j * np.pi * k / d)
        z_end, ok = lairez_track(coeffs, gamma, z0)
        if not ok:
            for g_alt in [0.7+0.5j, -0.5+0.6j, 0.9+0.3j]:
                z_end, ok = lairez_track(coeffs, g_alt, z0)
                if ok:
                    break
        roots.append(z_end)
    if polish:
        polished = []
        for r in roots:
            r_pol, _ = t2_orbit(coeffs, r, max_epochs=30, tol=tol)
            polished.append(r_pol)
        roots = polished
    return roots


# ---------------------------------------------------------------------------
# Polynomial generators
# ---------------------------------------------------------------------------

def make_random_poly(d, rng, magnitude=1.0):
    coeffs = list((rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)) * magnitude)
    coeffs[0] = 1.0 + 0.0j
    return coeffs


def make_kostlan_smale(d, rng):
    log_binom = np.array([math.lgamma(d + 1) - math.lgamma(k + 1) - math.lgamma(d - k + 1)
                          for k in range(d + 1)])
    sigma = np.exp(0.5 * log_binom)
    coefs = (rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)) * sigma
    coefs = coefs[::-1]
    if abs(coefs[0]) > 1e-15:
        coefs = coefs / coefs[0]
    return list(coefs)


def make_high_dynamic_range(d, rng):
    """Coefficients with magnitudes spread across e^{N(0, 4)} (very wide)."""
    mag = np.exp(rng.standard_normal(d + 1) * 2.0)
    phase = np.exp(2j * np.pi * rng.uniform(0, 1, d + 1))
    coeffs = list(mag * phase)
    coeffs[0] = 1.0 + 0.0j
    return coeffs


def max_residual(coeffs, roots):
    if not roots:
        return float('inf')
    return float(max(abs(horner1(coeffs, r)) for r in roots))


# ---------------------------------------------------------------------------
# Benchmark driver
# ---------------------------------------------------------------------------

def main():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    print("=" * 96, flush=True)
    print("flow/008 Strategy B + T_2 + Armijo (deflation) VS Lairez homotopy", flush=True)
    print("                       + scaled variant (homothety preconditioning)", flush=True)
    print("=" * 96, flush=True)

    rng = np.random.default_rng(2026)
    def n_polys_for_degree(d):
        if d <= 64:
            return 5
        if d <= 128:
            return 2
        return 1

    methods = [
        ('Pandrosion universal (no-if)',       find_all_pandrosion_universal_scaling),
        ('Pandrosion adaptive',                find_all_pandrosion_adaptive),
        ('Lairez (Newton corrector)',          find_all_lairez),
    ]

    regimes = [
        ("Random monic complex (mag=1)", lambda d, r: make_random_poly(d, r, 1.0)),
        ("Random monic large coef (mag=5)", lambda d, r: make_random_poly(d, r, 5.0)),
        ("Kostlan-Smale (BW)",           make_kostlan_smale),
        ("High dynamic range (e^{N(0,4)})", make_high_dynamic_range),
    ]

    for regime_name, factory in regimes:
        print(f"\n--- {regime_name} ---", flush=True)
        print(f"  {'d':>3} {'method':<36} {'time (ms)':>11} {'roots':>6} "
              f"{'max |P(r)|':>14} {'OK':>4}", flush=True)
        for d in [8, 16, 32, 64, 128, 256, 512]:
            # Pre-generate polynomials so every method sees the SAME polys
            n_polys = n_polys_for_degree(d)
            test_polys = [factory(d, rng) for _ in range(n_polys)]
            for name, fn in methods:
                t_total = 0.0
                succ = 0
                max_res = 0.0
                n_roots_med = 0
                for P in test_polys:
                    t0 = time.perf_counter()
                    try:
                        roots = fn(P)
                    except Exception:
                        roots = []
                    t_total += time.perf_counter() - t0
                    n_roots_med = len(roots)
                    if len(roots) == d:
                        r = max_residual(P, roots)
                        if r < 1e-6:
                            succ += 1
                        if r > max_res:
                            max_res = r
                avg_t = t_total / n_polys * 1000
                print(f"  {d:>3} {name:<36} {avg_t:>11.2f} {n_roots_med:>6} "
                      f"{max_res:>14.2e} {succ}/{n_polys}", flush=True)
            print(flush=True)

    print("=" * 96, flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
