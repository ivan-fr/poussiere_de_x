"""
SU(2)-equivariant Atiyah-Sutcliffe toolkit
=========================================

A clean, exploratory Python implementation of all prerequisites needed to attack
Atiyah-Sutcliffe via the Pandrosion-SU(2) framework. Each section builds on the
previous, and the final section runs exploratory experiments.

Sections:
  1. Spinor calculus on CP^1 (Hopf, skew product, Vandermonde)
  2. Atiyah-Berry polynomials in homogeneous form
  3. SU(2) representation theory (Wigner D, Clebsch-Gordan formula)
  4. Closed-form formula for Gram off-diagonals via spinors
  5. Hadamard-Fischer / Schur complement / Bhatia bounds
  6. Real-stability tests and mixed characteristic polynomials
  7. Platonic-configuration exact computations (tetrahedron, cube, icosahedron)
  8. DPP / reproducing kernel structure
  9. Exploratory experiments: hunting for a structural breakthrough
"""

from __future__ import annotations
import math
import itertools
import numpy as np
from typing import Callable, Tuple, List


# =====================================================================
# Section 1: Spinor calculus on CP^1
# =====================================================================

def hopf_lift(x: np.ndarray) -> np.ndarray:
    """Lift x in S^2 to a unit spinor u in C^2.

    For x = (sin(t) cos(p), sin(t) sin(p), cos(t)):
        u = (cos(t/2), sin(t/2) e^{i p})
    Convention: u up to U(1) phase. Hopf-projection: x = (u^* sigma u) where
    sigma are Pauli matrices.
    """
    z = float(x[2])
    if z > 1 - 1e-13:
        return np.array([1.0 + 0j, 0.0 + 0j])
    if z < -1 + 1e-13:
        return np.array([0.0 + 0j, 1.0 + 0j])
    cos_half = math.sqrt((1 + z) / 2)
    sin_half = math.sqrt((1 - z) / 2)
    phi = math.atan2(float(x[1]), float(x[0]))
    return np.array([cos_half + 0j,
                     sin_half * complex(math.cos(phi), math.sin(phi))])


def spinor_skew(u: np.ndarray, v: np.ndarray) -> complex:
    """[u, v] = u^1 v^2 - u^2 v^1.  SU(2)-invariant skew bilinear form."""
    return u[0] * v[1] - u[1] * v[0]


def spinor_inner(u: np.ndarray, v: np.ndarray) -> complex:
    """<u, v> = bar(u^1) v^1 + bar(u^2) v^2.  Standard Hermitian on C^2."""
    return np.conj(u[0]) * v[0] + np.conj(u[1]) * v[1]


def spinor_vandermonde(spinors: List[np.ndarray]) -> np.ndarray:
    """V_{k, j} := (u_k^1)^j * (u_k^2)^{n-1-j}.  Determinant = prod_{i<j} [u_i, u_j]."""
    n = len(spinors)
    V = np.zeros((n, n), dtype=complex)
    for k, u in enumerate(spinors):
        for j in range(n):
            V[k, j] = u[0]**j * u[1]**(n - 1 - j)
    return V


def hopf_project(u: np.ndarray) -> np.ndarray:
    """Inverse Hopf: u in C^2 -> x in S^2.  x = (u^* sigma u)."""
    a, b = u[0], u[1]
    x1 = 2 * np.real(np.conj(a) * b)
    x2 = 2 * np.imag(np.conj(a) * b)
    x3 = abs(a)**2 - abs(b)**2
    return np.array([x1, x2, x3])


# =====================================================================
# Section 2: Atiyah-Berry polynomials
# =====================================================================

def stereo_proj(v: np.ndarray) -> complex:
    """Stereographic projection from north pole."""
    z = float(v[2])
    if z >= 1.0 - 1e-13:
        return complex('inf')
    return complex(float(v[0]), float(v[1])) / (1.0 - z)


def atiyah_berry_polys(points: List[np.ndarray]) -> List[np.ndarray]:
    """Build n Atiyah-Berry polynomials.  Returns list of coefficient arrays."""
    n = len(points)
    polys = []
    for i in range(n):
        finite = []
        n_inf = 0
        for j in range(n):
            if i == j:
                continue
            v = points[j] - points[i]
            v = v / np.linalg.norm(v)
            a = stereo_proj(v)
            if math.isinf(a.real) or math.isinf(a.imag):
                n_inf += 1
            else:
                finite.append(a)
        c = np.array([1.0 + 0j])
        for a in finite:
            c = np.convolve(c, np.array([-a, 1.0 + 0j]))
        if n_inf > 0:
            cc = np.zeros(len(c) + n_inf, dtype=complex)
            cc[n_inf:] = c
            c = cc
        polys.append(c)
    return polys


def atiyah_M_matrix(points: List[np.ndarray]) -> np.ndarray:
    """Coefficient matrix M, rows = Atiyah polynomials in monomial basis."""
    n = len(points)
    polys = atiyah_berry_polys(points)
    M = np.zeros((n, n), dtype=complex)
    for i, p in enumerate(polys):
        if len(p) < n:
            pp = np.zeros(n, dtype=complex)
            pp[:len(p)] = p
            p = pp
        elif len(p) > n:
            p = p[:n]
        M[i] = p
    return M


def su2_norm_sq(coeffs: np.ndarray, n: int) -> float:
    """SU(2)-invariant norm squared.  ||p||^2 = sum |c_k|^2 / C(n-1, k)."""
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    return float(np.sum(np.abs(coeffs)**2 / binom))


def gram_matrix(points: List[np.ndarray]) -> np.ndarray:
    """SU(2)-Gram matrix G = M W M^*, W = diag(1/C(n-1, k))."""
    n = len(points)
    M = atiyah_M_matrix(points)
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    W = np.diag(1.0 / binom)
    return M @ W @ M.conj().T


def gram_normalized(points: List[np.ndarray]) -> np.ndarray:
    """G_norm with unit diagonal."""
    G = gram_matrix(points)
    d = np.sqrt(np.maximum(np.real(np.diag(G)), 1e-300))
    return G / np.outer(d, d)


def atiyah_D_squared(points: List[np.ndarray]) -> Tuple[float, float]:
    """Returns (log|D|^2, L_n) where L_n = sum log C(n-1, k).

    Atiyah-Sutcliffe: log|D|^2 >= 0.
    Identity: log|D|^2 = log det G_norm + L_n.
    """
    n = len(points)
    M = atiyah_M_matrix(points)
    binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
    log_norms_2 = np.log(np.maximum(
        np.sum(np.abs(M)**2 / binom, axis=1), 1e-300))
    sg, logabs = np.linalg.slogdet(M)
    log_D_sq = 2 * logabs - log_norms_2.sum()
    L_n = float(np.sum(np.log(binom)))
    return float(log_D_sq), L_n


# =====================================================================
# Section 3: SU(2) representation theory
# =====================================================================

def wigner_D_matrix(j: float, R: np.ndarray) -> np.ndarray:
    """Compute Wigner D-matrix for spin j given a 3x3 rotation R.

    Strategy: factor R = R_z(alpha) R_y(beta) R_z(gamma) (Euler angles ZYZ),
    then use closed-form formula for D^j_{m,m'}(alpha, beta, gamma).
    """
    # Extract Euler angles ZYZ from R
    # R = R_z(alpha) R_y(beta) R_z(gamma)
    # cos(beta) = R[2,2]
    cos_b = float(np.clip(R[2, 2], -1.0, 1.0))
    beta = math.acos(cos_b)
    if abs(math.sin(beta)) > 1e-10:
        alpha = math.atan2(R[1, 2], R[0, 2])
        gamma = math.atan2(R[2, 1], -R[2, 0])
    else:
        alpha = math.atan2(R[1, 0], R[0, 0])
        gamma = 0.0

    # Build D^j matrix
    dim = int(round(2*j + 1))
    D = np.zeros((dim, dim), dtype=complex)
    ms = [j - i for i in range(dim)]  # m = j, j-1, ..., -j
    for i_m, m in enumerate(ms):
        for i_mp, mp in enumerate(ms):
            d_small = wigner_d_small(j, m, mp, beta)
            D[i_m, i_mp] = (np.exp(-1j * m * alpha) * d_small *
                            np.exp(-1j * mp * gamma))
    return D


def wigner_d_small(j: float, m: float, mp: float, beta: float) -> float:
    """Wigner small d-function d^j_{m,mp}(beta). Closed-form via factorial sum."""
    cos_h = math.cos(beta / 2)
    sin_h = math.sin(beta / 2)
    if abs(sin_h) < 1e-15:
        return 1.0 if abs(m - mp) < 1e-9 else 0.0
    # factor outside sum
    pre_num = (math.factorial(int(j+m)) * math.factorial(int(j-m)) *
               math.factorial(int(j+mp)) * math.factorial(int(j-mp)))
    pre = math.sqrt(pre_num)
    # sum over s
    total = 0.0
    s_min = int(max(0, mp - m))
    s_max = int(min(j + mp, j - m))
    for s in range(s_min, s_max + 1):
        denom = (math.factorial(int(j + mp - s)) *
                 math.factorial(s) *
                 math.factorial(int(m - mp + s)) *
                 math.factorial(int(j - m - s)))
        sign = (-1) ** (m - mp + s)
        c_pow = 2*j + mp - m - 2*s
        s_pow = m - mp + 2*s
        total += sign * (cos_h ** c_pow) * (sin_h ** s_pow) / denom
    return pre * total


def clebsch_gordan(j1: float, m1: float, j2: float, m2: float,
                   J: float, M: float) -> float:
    """Clebsch-Gordan coefficient <j1 m1; j2 m2 | J M>. Closed-form formula."""
    if abs(m1 + m2 - M) > 1e-9:
        return 0.0
    if J < abs(j1 - j2) - 1e-9 or J > j1 + j2 + 1e-9:
        return 0.0
    pre = math.sqrt(
        (2*J + 1) * math.factorial(int(J + j1 - j2)) *
        math.factorial(int(J - j1 + j2)) *
        math.factorial(int(j1 + j2 - J)) /
        math.factorial(int(J + j1 + j2 + 1))
    )
    pre *= math.sqrt(
        math.factorial(int(J + M)) * math.factorial(int(J - M)) *
        math.factorial(int(j1 - m1)) * math.factorial(int(j1 + m1)) *
        math.factorial(int(j2 - m2)) * math.factorial(int(j2 + m2))
    )
    total = 0.0
    k_min = int(max(0, j2 - J - m1, j1 - J + m2))
    k_max = int(min(j1 + j2 - J, j1 - m1, j2 + m2))
    for k in range(k_min, k_max + 1):
        denom = (math.factorial(k) * math.factorial(int(j1 + j2 - J - k)) *
                 math.factorial(int(j1 - m1 - k)) *
                 math.factorial(int(j2 + m2 - k)) *
                 math.factorial(int(J - j2 + m1 + k)) *
                 math.factorial(int(J - j1 - m2 + k)))
        total += (-1) ** k / denom
    return pre * total


# =====================================================================
# Section 4: Closed-form spinor formula for Gram off-diagonals
# =====================================================================

def gram_via_spinors(points: List[np.ndarray]) -> np.ndarray:
    """Compute the Gram matrix using spinor lifts directly.

    For Atiyah-Berry: each p_i is associated with spinor u_i. The off-diagonal
    <p_i, p_j>_SU(2) has a closed form via spinor invariants. We compute it
    as a verification / alternative pathway.

    For the LAGRANGE version (P_i normalized to P_i(u_i) = 1), the Gram
    matrix off-diagonals are related to the spinor inner products as:
        <P_i^Lag, P_j^Lag> = (overlap of n-th symmetric powers).
    """
    n = len(points)
    spinors = [hopf_lift(x) for x in points]
    # Build Lagrange basis polynomials in homogeneous form
    G_lag = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            # <P_i^Lag, P_j^Lag>_SU(2) = ?
            # Using P_i(z, w) = prod_{k != i} (u_k^1 w - u_k^2 z) / [u_i, u_k]
            # The SU(2) inner product on Sym^{n-1}(C^2) of products of linear
            # forms has the permanent-like formula:
            #   <prod L_k, prod M_k> = (1/(n-1)!) sum_sigma prod <L_k, M_sigma(k)>
            # where <l, m>_C2 = bar(l^1) m^1 + bar(l^2) m^2 for l = (l^1)z + (l^2)w.
            # For our convention: linear form (u^1 w - u^2 z) = -u^2 z + u^1 w
            # has C^2-coefficients (-u^2, u^1). Inner product:
            # <(-u_i^2, u_i^1), (-u_j^2, u_j^1)> = bar(-u_i^2)(-u_j^2) + bar(u_i^1)(u_j^1)
            #                                   = bar(u_i^2) u_j^2 + bar(u_i^1) u_j^1
            #                                   = <u_i, u_j>_C2 (Hermitian inner product)
            num = 0.0 + 0j
            # Lists of "other" indices
            others_i = [k for k in range(n) if k != i]
            others_j = [k for k in range(n) if k != j]
            # sum over permutations sigma: others_i -> others_j
            for perm in itertools.permutations(range(len(others_j))):
                prod = 1.0 + 0j
                for a, k_i in enumerate(others_i):
                    k_j = others_j[perm[a]]
                    prod *= spinor_inner(spinors[k_i], spinors[k_j])
                num += prod
            # Scale by symmetric-power normalization 1/(n-1)!
            num /= math.factorial(n - 1)
            # Lagrange normalization factor: divide by [u_i, u_*] [u_j, u_*]^* products
            denom_i = 1.0 + 0j
            denom_j = 1.0 + 0j
            for k in others_i:
                denom_i *= spinor_skew(spinors[i], spinors[k])
            for k in others_j:
                denom_j *= spinor_skew(spinors[j], spinors[k])
            G_lag[i, j] = num / (denom_i * np.conj(denom_j))
    return G_lag


# =====================================================================
# Section 5: Hadamard-Fischer, Schur complement, Bhatia
# =====================================================================

def hadamard_upper(G: np.ndarray) -> float:
    """Hadamard upper bound: det G <= prod G_ii."""
    return float(np.real(np.prod(np.diag(G))))


def schur_residuals(G: np.ndarray) -> np.ndarray:
    """r_i^2 = g_i^* G_<i^{-1} g_i. det G = prod (G_ii - r_i^2 G_<i diag), or for unit diag,
    det G = prod (1 - r_i^2). Sum of -log(1-r^2) = -log det.
    """
    n = G.shape[0]
    r2 = np.zeros(n)
    for i in range(n):
        if i == 0:
            r2[i] = 0.0
            continue
        G_lt = G[:i, :i]
        g_i = G[:i, i]
        try:
            sol = np.linalg.solve(G_lt, g_i)
            r2[i] = float(np.real(np.dot(g_i.conj(), sol)))
        except np.linalg.LinAlgError:
            r2[i] = 1.0
    return r2


def bhatia_perturbation_bound(G: np.ndarray) -> float:
    """If ||G - I||_op < 1: -log det G <= ||R||_F^2 / (1 - ||R||_op).
    Returns +inf if precondition fails."""
    n = G.shape[0]
    R = G - np.eye(n)
    op = np.linalg.norm(R, 2)
    if op >= 1.0:
        return float('inf')
    fro_sq = np.linalg.norm(R, 'fro')**2
    return float(fro_sq / (1 - op))


def hadamard_fischer_iterated(G: np.ndarray) -> float:
    """det G via iterated Schur complement. Equivalent to det G = prod (1-r_i^2)
    for unit-diagonal G. Return -log det."""
    r2 = schur_residuals(G)
    r2 = np.clip(r2, 0.0, 1.0 - 1e-15)
    return float(-np.sum(np.log(1 - r2)))


# =====================================================================
# Section 6: Real-stability tests, mixed characteristic polynomial
# =====================================================================

def is_real_rooted(coeffs: np.ndarray, tol: float = 1e-8) -> bool:
    """Check whether univariate poly has all real roots.

    coeffs: high-to-low order (numpy convention).
    """
    roots = np.roots(coeffs)
    return all(abs(r.imag) <= tol * (1 + abs(r.real)) for r in roots)


def mixed_char_poly(matrices: List[np.ndarray]) -> np.ndarray:
    """Mixed characteristic polynomial:
        mu(t) = E_{eps in {0,1}^m}[det(tI - sum_i eps_i M_i)]
    Computed by enumeration (exponential in m, only feasible for small m).
    """
    m = len(matrices)
    if m == 0:
        return np.array([1.0])
    n = matrices[0].shape[0]
    avg_coeffs = np.zeros(n + 1)
    for s in range(2**m):
        eps = [(s >> i) & 1 for i in range(m)]
        S = sum(eps[i] * matrices[i] for i in range(m))
        if not isinstance(S, np.ndarray):
            S = np.zeros((n, n))
        coeffs = np.poly(np.real(np.linalg.eigvalsh(S)))
        avg_coeffs += np.real(coeffs)
    avg_coeffs /= 2**m
    return avg_coeffs


def mss_barrier(coeffs: np.ndarray, t: float) -> float:
    """MSS barrier function Phi(t) = -p'(t)/p(t) for t > max root."""
    p = coeffs
    pp = np.polyder(p)
    p_t = np.polyval(p, t)
    pp_t = np.polyval(pp, t)
    return -pp_t / p_t


# =====================================================================
# Section 7: Platonic configurations (exact Atiyah-Sutcliffe verification)
# =====================================================================

def tetrahedron_S2() -> List[np.ndarray]:
    """4 vertices of a regular tetrahedron inscribed in S^2."""
    s = 1.0 / math.sqrt(3)
    return [np.array([s, s, s]),
            np.array([s, -s, -s]),
            np.array([-s, s, -s]),
            np.array([-s, -s, s])]


def octahedron_S2() -> List[np.ndarray]:
    """6 vertices of regular octahedron."""
    return [np.array([1, 0, 0]), np.array([-1, 0, 0]),
            np.array([0, 1, 0]), np.array([0, -1, 0]),
            np.array([0, 0, 1]), np.array([0, 0, -1])]


def cube_S2() -> List[np.ndarray]:
    """8 vertices of cube inscribed in S^2."""
    s = 1.0 / math.sqrt(3)
    pts = []
    for sx in [-1, 1]:
        for sy in [-1, 1]:
            for sz in [-1, 1]:
                pts.append(np.array([sx*s, sy*s, sz*s]))
    return pts


def icosahedron_S2() -> List[np.ndarray]:
    """12 vertices of regular icosahedron."""
    phi = (1 + math.sqrt(5)) / 2
    nrm = math.sqrt(1 + phi**2)
    pts = []
    for s1 in [-1, 1]:
        for s2 in [-1, 1]:
            pts.append(np.array([0, s1/nrm, s2*phi/nrm]))
            pts.append(np.array([s1/nrm, s2*phi/nrm, 0]))
            pts.append(np.array([s2*phi/nrm, 0, s1/nrm]))
    return pts


def dodecahedron_S2() -> List[np.ndarray]:
    """20 vertices of regular dodecahedron."""
    phi = (1 + math.sqrt(5)) / 2
    inv = 1.0 / phi
    nrm = math.sqrt(3)
    pts = []
    for s1 in [-1, 1]:
        for s2 in [-1, 1]:
            for s3 in [-1, 1]:
                pts.append(np.array([s1, s2, s3]) / nrm)
    for s1 in [-1, 1]:
        for s2 in [-1, 1]:
            pts.append(np.array([0, s1*inv, s2*phi]) / nrm)
            pts.append(np.array([s1*inv, s2*phi, 0]) / nrm)
            pts.append(np.array([s2*phi, 0, s1*inv]) / nrm)
    return pts


def equispaced_circle(n: int, eps: float = 0.0) -> List[np.ndarray]:
    """n points equispaced on a circle at height z = eps."""
    pts = []
    r = math.sqrt(1 - eps**2)
    for k in range(n):
        theta = 2 * math.pi * k / n
        pts.append(np.array([r*math.cos(theta), r*math.sin(theta), eps]))
    return pts


# =====================================================================
# Section 8: DPP / reproducing kernel structure
# =====================================================================

def reproducing_kernel_matrix(spinors: List[np.ndarray], n: int) -> np.ndarray:
    """K_{ij} := <u_i, u_j>^{n-1}, the (n-1)-th symmetric power of the standard
    Hermitian inner product on C^2. This is the reproducing kernel of
    Sym^{n-1}(C^2) in coherent-state language.

    For Atiyah: this kernel is the natural SU(2)-equivariant matrix
    associated to the configuration. Its determinant is a SU(2) invariant
    of the configuration.
    """
    n_pts = len(spinors)
    K = np.zeros((n_pts, n_pts), dtype=complex)
    for i in range(n_pts):
        for j in range(n_pts):
            ip = spinor_inner(spinors[i], spinors[j])
            K[i, j] = ip ** (n - 1)
    return K


def dpp_kernel_atiyah_relate(points: List[np.ndarray]) -> dict:
    """Compute the relationship between log|D|^2 (Atiyah-Berry) and
    log det K (DPP-like coherent-state kernel)."""
    n = len(points)
    spinors = [hopf_lift(x) for x in points]
    K = reproducing_kernel_matrix(spinors, n)
    sg, logabs = np.linalg.slogdet(K)
    logD_sq, L_n = atiyah_D_squared(points)
    return dict(
        log_det_K=float(np.real(logabs)),
        log_D_sq=logD_sq,
        L_n=L_n,
        diff=logD_sq - float(np.real(logabs)),
    )


# =====================================================================
# Section 9: Exploratory experiments — hunting for breakthrough
# =====================================================================

def random_S2_points(rng, n: int) -> List[np.ndarray]:
    pts = []
    for _ in range(n):
        v = rng.standard_normal(3)
        pts.append(v / np.linalg.norm(v))
    return pts


def antipodal_split(rng, n: int) -> List[np.ndarray]:
    pts = []
    for _ in range(n // 2):
        v = rng.standard_normal(3) * 0.1
        v[2] += 1.0
        pts.append(v / np.linalg.norm(v))
    for _ in range(n - n // 2):
        v = rng.standard_normal(3) * 0.1
        v[2] -= 1.0
        pts.append(v / np.linalg.norm(v))
    return pts


# =====================================================================
# Verification / exploration runner
# =====================================================================

def run_all():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 100, flush=True)
    print("SU(2)-EQUIVARIANT ATIYAH-SUTCLIFFE TOOLKIT: SELF-TEST + EXPLORATION",
          flush=True)
    print("=" * 100, flush=True)

    # --- Test 1: spinor calculus
    print("\n[1] Spinor calculus self-test:")
    rng = np.random.default_rng(2026)
    pts = random_S2_points(rng, 5)
    spinors = [hopf_lift(x) for x in pts]
    # Hopf inverse: u -> x, then check
    for i, (x, u) in enumerate(zip(pts, spinors)):
        x_back = hopf_project(u)
        err = np.linalg.norm(x - x_back)
        print(f"  point {i}: hopf round-trip error = {err:.2e}")
    # Spinor Vandermonde identity
    V = spinor_vandermonde(spinors)
    sg, ld = np.linalg.slogdet(V)
    log_pred = sum(2 * math.log(abs(spinor_skew(spinors[i], spinors[j])))
                  for i in range(5) for j in range(i+1, 5))
    print(f"  log|det V|^2 = {2*ld:.4f}, predicted = {log_pred:.4f}, "
          f"diff = {2*ld - log_pred:.2e}")

    # --- Test 2: Atiyah-Berry vs Lagrange via spinors
    print("\n[2] Atiyah-Berry |D|^2 via direct M-construction (n=5):")
    logD_sq, L_n = atiyah_D_squared(pts)
    print(f"  log|D|^2 = {logD_sq:.4f}, L_n = {L_n:.4f}")
    print(f"  Atiyah-Sutcliffe: log|D|^2 >= 0?  {logD_sq >= 0}  "
          f"(margin {logD_sq:.4f})")

    # --- Test 3: Wigner D-matrix self-test
    print("\n[3] Wigner D-matrix self-test (j = 1):")
    # Identity rotation
    D_id = wigner_D_matrix(1.0, np.eye(3))
    print(f"  D^1(I): max off-diag = {np.max(np.abs(D_id - np.eye(3))):.2e}")
    # Clebsch-Gordan: <1/2, 1/2; 1/2, -1/2 | 1, 0> = 1/sqrt(2)
    cg = clebsch_gordan(0.5, 0.5, 0.5, -0.5, 1.0, 0.0)
    print(f"  CG(1/2,1/2; 1/2,-1/2 | 1,0) = {cg:.6f}, expected = {1/math.sqrt(2):.6f}")

    # --- Test 4: Closed-form Gram via spinors (small n only)
    print("\n[4] Spinor formula for Lagrange Gram (n=4):")
    pts4 = random_S2_points(rng, 4)
    spinors4 = [hopf_lift(x) for x in pts4]
    G_lag = gram_via_spinors(pts4)
    print(f"  Lagrange Gram (size {G_lag.shape}):")
    print(f"  Hermitian? max|G_lag - G_lag^*| = {np.max(np.abs(G_lag - G_lag.conj().T)):.2e}")
    eigs_lag = np.real(np.linalg.eigvalsh((G_lag + G_lag.conj().T)/2))
    print(f"  eigenvalues: {eigs_lag}")
    print(f"  PSD? min eigenvalue = {eigs_lag.min():.4f}")

    # --- Test 5: Reproducing kernel matrix (DPP)
    print("\n[5] Reproducing-kernel structure K_{ij} = <u_i, u_j>^{n-1}:")
    for n in [4, 6, 8, 10]:
        rng2 = np.random.default_rng(42 + n)
        pts_n = random_S2_points(rng2, n)
        result = dpp_kernel_atiyah_relate(pts_n)
        print(f"  n = {n}: log|D|^2 = {result['log_D_sq']:.4f}, "
              f"log det K = {result['log_det_K']:.4f}, "
              f"diff = {result['diff']:.4f}, L_n = {result['L_n']:.4f}")

    # --- Test 6: Schur complement decomposition exactness
    print("\n[6] Schur complement exactness (sum -log(1-r^2) = -log det):")
    for n in [4, 6, 10]:
        rng3 = np.random.default_rng(100 + n)
        pts_n = antipodal_split(rng3, n)
        Gn = gram_normalized(pts_n)
        r2 = schur_residuals(Gn)
        r2 = np.clip(r2, 0.0, 1.0 - 1e-15)
        sum_neg = -np.sum(np.log(1 - r2))
        sg, ld = np.linalg.slogdet(Gn)
        print(f"  n = {n}: schur sum = {sum_neg:.4f}, -log det = {-ld:.4f}, "
              f"diff = {sum_neg - (-ld):.2e}")

    # --- Test 7: Bhatia bound applicability
    print("\n[7] Bhatia perturbation bound:")
    for n in [4, 6, 10, 15, 20]:
        rng4 = np.random.default_rng(200 + n)
        pts_n = random_S2_points(rng4, n)
        Gn = gram_normalized(pts_n)
        b = bhatia_perturbation_bound(Gn)
        sg, ld = np.linalg.slogdet(Gn)
        print(f"  n = {n}: -log det = {-ld:.4f}, Bhatia bound = "
              f"{'inf (||R||_op>=1)' if b == float('inf') else f'{b:.4f}'}")

    # --- Test 8: Mixed char poly real-stability
    print("\n[8] Mixed char poly under reflection randomization (small n):")
    n = 4
    pts_n = antipodal_split(np.random.default_rng(1234), n)
    G_base = gram_normalized(pts_n)
    # Build "rank-1 like" decomposition: G = sum_k row_k row_k^* (Cholesky-like)
    # Use eigendecomposition: G = U L U^* = sum eig_k * u_k u_k^*
    eigs, U = np.linalg.eigh(G_base)
    matrices = [eigs[k] * np.outer(U[:, k], U[:, k].conj()) for k in range(n)]
    mu = mixed_char_poly(matrices)
    is_rs = is_real_rooted(mu)
    roots = np.roots(mu)
    print(f"  base G eigenvalues: {eigs}")
    print(f"  mixed char poly real-rooted? {is_rs}")
    print(f"  mu roots: {sorted([(r.real, r.imag) for r in roots])}")

    # --- Test 9: Platonic configurations (exact!)
    print("\n[9] Atiyah-Sutcliffe on Platonic solids (high-symmetry exact cases):")
    plat = [
        ("tetrahedron", tetrahedron_S2(), 4),
        ("octahedron", octahedron_S2(), 6),
        ("cube", cube_S2(), 8),
        ("icosahedron", icosahedron_S2(), 12),
        ("dodecahedron", dodecahedron_S2(), 20),
    ]
    print(f"  {'config':>15} {'n':>4} {'log|D|^2':>10} {'L_n':>10} {'log det G':>11} "
          f"{'beta':>8}")
    for name, pts_p, n in plat:
        try:
            logD_sq, L_n = atiyah_D_squared(pts_p)
            Gn = gram_normalized(pts_p)
            sg, ld = np.linalg.slogdet(Gn)
            beta = -float(ld) / L_n if L_n > 0 else 0.0
            print(f"  {name:>15} {n:>4} {logD_sq:>10.4f} {L_n:>10.4f} "
                  f"{ld:>11.4f} {beta:>8.4f}")
        except Exception as e:
            print(f"  {name}: error {e}")

    # --- Test 10: HUNT — does log det K (DPP) bound log|D|^2?
    print("\n[10] HUNT — relationship between DPP kernel and Atiyah |D|:")
    print(f"  {'n':>4} {'log|D|^2':>10} {'log det K':>11} {'diff':>10} {'L_n':>10}")
    for n in [4, 6, 8, 10, 12, 15]:
        rng_h = np.random.default_rng(2026 + n)
        diffs = []
        for _ in range(20):
            pts_h = random_S2_points(rng_h, n)
            res = dpp_kernel_atiyah_relate(pts_h)
            diffs.append(res['diff'])
        diffs = np.array(diffs)
        # Test: is the diff a function of L_n only?
        binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        print(f"  {n:>4} {'':>10} {'':>11} mean diff = {diffs.mean():>+.4f}, "
              f"std = {diffs.std():>.4f}, L_n = {L_n:.2f}")

    # --- Test 11: HUNT — is log|D|^2 + log det K constant (per n)?
    print("\n[11] Test: log|D|^2 + log det K = constant of n only?")
    for n in [4, 5, 6, 8, 10, 12]:
        rng_h = np.random.default_rng(3026 + n)
        sums = []
        for _ in range(30):
            pts_h = random_S2_points(rng_h, n)
            res = dpp_kernel_atiyah_relate(pts_h)
            sums.append(res['log_D_sq'] + res['log_det_K'])
        sums = np.array(sums)
        binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        # log K + log|D|^2 ?? compare with various candidates
        print(f"  n = {n}: mean = {sums.mean():+.4f}, std = {sums.std():.4f}, "
              f"L_n = {L_n:.2f}, ratio mean/L_n = {sums.mean()/L_n:.4f}")

    # --- Test 12: HUNT — what about log|D|^2 - log det K?
    print("\n[12] Test: log|D|^2 - log det K (signed difference)?")
    for n in [4, 5, 6, 8, 10, 12, 15]:
        rng_h = np.random.default_rng(4026 + n)
        diffs = []
        for _ in range(30):
            pts_h = random_S2_points(rng_h, n)
            res = dpp_kernel_atiyah_relate(pts_h)
            diffs.append(res['log_D_sq'] - res['log_det_K'])
        diffs = np.array(diffs)
        binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        print(f"  n = {n}: log|D|^2 - log det K: mean = {diffs.mean():+.4f}, "
              f"std = {diffs.std():.4f}, ratio/L_n = {diffs.mean()/L_n:.4f}")

    # --- Test 13: HUNT — DPP kernel positivity / properties
    print("\n[13] DPP kernel K_{ij} = <u_i, u_j>^{n-1}: spectral properties")
    print(f"  {'n':>4} {'min eig K':>11} {'max eig K':>11} {'cond K':>11} "
          f"{'min eig G_norm':>15}")
    for n in [4, 6, 8, 10, 12]:
        rng_h = np.random.default_rng(5026 + n)
        pts_h = antipodal_split(rng_h, n)
        spinors_h = [hopf_lift(x) for x in pts_h]
        K = reproducing_kernel_matrix(spinors_h, n)
        K_herm = (K + K.conj().T) / 2
        eigs_K = np.real(np.linalg.eigvalsh(K_herm))
        Gn = gram_normalized(pts_h)
        eigs_G = np.real(np.linalg.eigvalsh(Gn))
        cond = float(eigs_K.max() / max(abs(eigs_K.min()), 1e-300))
        print(f"  {n:>4} {eigs_K.min():>11.4e} {eigs_K.max():>11.4e} "
              f"{cond:>11.4e} {eigs_G.min():>15.4e}")


if __name__ == "__main__":
    run_all()
