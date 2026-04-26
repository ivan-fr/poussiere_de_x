"""
KEY HYPOTHESIS:  det G_norm  >=  det K   for ALL configurations.

If TRUE, then  -log det G_norm <= -log det K, and we just need to bound
-log det K analytically. The DPP kernel K has KNOWN structure:

  K_{ij} = <u_i, u_j>^{n-1}  with |<u_i, u_j>|^2 = (1 + cos(theta_ij))/2

So K = (M_inner)^{*(n-1)} where M_inner_{ij} = <u_i, u_j> is the unit-diag
PSD Gram of the spinors themselves (which is easier to analyze).

By Schur product (Hadamard product preserving PSD), and the fact that
M_inner is PSD with unit diag, K is also PSD with unit diag.

Test 1: verify det G_norm >= det K on every individual configuration.
Test 2: bound -log det K via the spinor inner product matrix.
"""
from __future__ import annotations
import math
import numpy as np
from su2_atiyah_toolkit import (
    hopf_lift, atiyah_D_squared, gram_normalized,
    reproducing_kernel_matrix, random_S2_points, antipodal_split,
    spinor_inner, equispaced_circle,
)


def spinor_inner_matrix(spinors):
    """M_{ij} = <u_i, u_j>. Hermitian, PSD, unit diag."""
    n = len(spinors)
    M = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            M[i, j] = spinor_inner(spinors[i], spinors[j])
    return M


def main():
    np.seterr(over='ignore', invalid='ignore', divide='ignore')

    print("=" * 100, flush=True)
    print("TEST 1: Is det G_norm >= det K for EVERY configuration?", flush=True)
    print("=" * 100, flush=True)

    violations = 0
    total = 0
    worst_violation = 0.0
    worst_excess = 0.0

    for family_name, gen in [
        ('random', lambda rng, n: random_S2_points(rng, n)),
        ('antipodal', lambda rng, n: antipodal_split(rng, n)),
        ('equispaced', lambda rng, n: equispaced_circle(n, 0.01)),
        ('equispaced-tiny', lambda rng, n: equispaced_circle(n, 1e-5)),
    ]:
        print(f"\n{family_name}:")
        for n in [4, 5, 6, 8, 10, 12, 15, 20]:
            rng = np.random.default_rng(42 + n)
            n_cfg = 30 if family_name not in ['equispaced', 'equispaced-tiny'] else 1
            cfg_violations = 0
            cfg_total = 0
            for _ in range(n_cfg):
                pts = gen(rng, n)
                spinors = [hopf_lift(x) for x in pts]
                K = reproducing_kernel_matrix(spinors, n)
                K_h = (K + K.conj().T) / 2
                Gn = gram_normalized(pts)
                sg1, ld_K = np.linalg.slogdet(K_h)
                sg2, ld_G = np.linalg.slogdet(Gn)
                if np.isfinite(ld_K) and np.isfinite(ld_G):
                    cfg_total += 1
                    total += 1
                    diff = float(ld_G) - float(np.real(ld_K))
                    if diff < -1e-10:  # G < K
                        cfg_violations += 1
                        violations += 1
                        if diff < worst_violation:
                            worst_violation = diff
                    if diff > worst_excess:
                        worst_excess = diff
            if cfg_total:
                print(f"  n={n}: {cfg_total} configs, "
                      f"{cfg_violations} violations of det G >= det K")
    print(f"\nGRAND TOTAL: {violations}/{total} violations")
    print(f"  worst violation (most negative diff): {worst_violation:.4e}")
    print(f"  largest excess (det G/det K ratio): exp({worst_excess:.4f})")

    if violations == 0:
        print("\n*** HYPOTHESIS CONFIRMED: det G_norm >= det K UNIFORMLY ***")
    else:
        print(f"\n*** HYPOTHESIS REFUTED: {violations} counterexamples ***")

    # Test 2: bound -log det K
    print("\n" + "=" * 100, flush=True)
    print("TEST 2: Can we bound -log det K analytically?", flush=True)
    print("=" * 100, flush=True)
    print("\nKey identity: K = M_inner^{*(n-1)} (Schur power)")
    print("where M_inner_{ij} = <u_i, u_j>, and Schur power preserves PSD.")
    print("\nClassical DPP fact: -log det K <= -(n-1) log det M_inner ?")
    print("Test:")
    print(f"\n{'n':>4} {'-log det K':>13} {'-(n-1) log det M_inner':>23} "
          f"{'L_n':>10} {'-log det K / L_n':>20}",
          flush=True)
    for n in [4, 5, 6, 8, 10, 12, 15, 20, 30, 50]:
        rng = np.random.default_rng(2026 + n)
        # Worst-case over multiple families
        worst_K = 0.0
        worst_M = 0.0
        for trial_family in [random_S2_points, antipodal_split]:
            for _ in range(20):
                pts = trial_family(rng, n)
                spinors = [hopf_lift(x) for x in pts]
                K = reproducing_kernel_matrix(spinors, n)
                K_h = (K + K.conj().T) / 2
                M_inner = spinor_inner_matrix(spinors)
                M_h = (M_inner + M_inner.conj().T) / 2
                sg1, ld_K = np.linalg.slogdet(K_h)
                sg2, ld_M = np.linalg.slogdet(M_h)
                if np.isfinite(ld_K) and np.isfinite(ld_M):
                    if -float(np.real(ld_K)) > worst_K:
                        worst_K = -float(np.real(ld_K))
                        worst_M = -float(np.real(ld_M)) * (n-1)
        binom = np.array([math.comb(n - 1, k) for k in range(n)], dtype=float)
        L_n = float(np.sum(np.log(binom)))
        print(f"{n:>4} {worst_K:>13.4f} {worst_M:>23.4f} {L_n:>10.2f} "
              f"{worst_K/L_n:>20.4f}",
              flush=True)

    # Test 3: M_inner has very specific structure - it's a Gram of spinors
    # which is a SU(2)-invariant function of (x_i)
    print("\n" + "=" * 100, flush=True)
    print("TEST 3: Structure of M_inner (Gram of spinors)", flush=True)
    print("=" * 100, flush=True)
    n = 6
    pts = random_S2_points(np.random.default_rng(12345), n)
    spinors = [hopf_lift(x) for x in pts]
    M_inner = spinor_inner_matrix(spinors)
    print(f"M_inner (n=6 random):")
    eigs_M = np.real(np.linalg.eigvalsh((M_inner + M_inner.conj().T)/2))
    print(f"  eigenvalues: {eigs_M}")
    print(f"  rank: {np.sum(eigs_M > 1e-10)}")
    print(f"  M_inner has rank 2 (because spinors live in C^2 = 2-dim)!")

    print("\nKey observation: M_inner = U^* U where U is the 2 x n matrix of spinors.")
    print("So M_inner has rank exactly 2, and det M_inner = 0 for n >= 3.")
    print("Therefore -log det M_inner = +inf, useless as bound.")
    print("\nBUT: K = M_inner^{(n-1)} (Schur/Hadamard power) has FULL rank!")
    print("Schur-product theorem: A o B is PSD if A, B are. Iterated power preserves PSD.")
    print("Schur powers can be analyzed via Hadamard's inequality on Schur products...")

    # Let's see what structure K has
    print("\nDirect computation of K eigenvalues (n=6 random):")
    K_test = reproducing_kernel_matrix(spinors, n)
    K_h = (K_test + K_test.conj().T) / 2
    eigs_K = np.real(np.linalg.eigvalsh(K_h))
    print(f"  K eigenvalues: {eigs_K}")
    print(f"  K is full-rank PSD with unit diagonal — like a correlation matrix.")

    # Test 4: if det G_norm >= det K, can we approximate?
    print("\n" + "=" * 100, flush=True)
    print("TEST 4: Can we PROVE det G_norm >= det K via Schur structure?", flush=True)
    print("=" * 100, flush=True)
    print("Conjecture: G_norm = K + R where R is also PSD?")
    print("This would imply det G_norm = det(K + R) >= det K (since K PSD, R PSD).")
    print("\nNumerical test: is G_norm - K PSD?")
    for n in [4, 6, 8, 10]:
        rng = np.random.default_rng(99 + n)
        violations = 0
        for _ in range(30):
            pts = random_S2_points(rng, n)
            spinors = [hopf_lift(x) for x in pts]
            K = reproducing_kernel_matrix(spinors, n)
            K_h = (K + K.conj().T) / 2
            Gn = gram_normalized(pts)
            R = Gn - K_h
            R_h = (R + R.conj().T) / 2
            eigs_R = np.real(np.linalg.eigvalsh(R_h))
            if eigs_R.min() < -1e-10:
                violations += 1
        print(f"  n={n}: G_norm - K is PSD? "
              f"violations: {violations}/30")


if __name__ == "__main__":
    main()
