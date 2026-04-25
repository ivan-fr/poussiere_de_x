import numpy as np

def eval_poly(roots, z):
    return np.prod(z - roots)

def Q_diff(roots, z0, z):
    if abs(z - z0) < 1e-18:
        mat = z0 - roots
        prod_tot = np.prod(mat)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.sum(np.where(mat != 0, prod_tot / mat, 0))
    return (eval_poly(roots, z) - eval_poly(roots, z0)) / (z - z0)

def pandrosion_T3_step(roots, z0, z):
    P0 = eval_poly(roots, z0)
    Q = Q_diff(roots, z0, z)
    if abs(Q) < 1e-60: return None
    s1 = z0 - P0 / Q
    Q2 = Q_diff(roots, z0, s1)
    if abs(Q2) < 1e-60: return s1
    s2 = z0 - P0 / Q2
    den = s2 - 2 * s1 + z
    if abs(den) < 1e-60: return s1
    return z - (s1 - z) ** 2 / den

roots = np.array([-0.96731984-0.06459744j,  0.22987706-0.45925592j,  0.06489314-0.22842954j])
z0 = -0.6496697856287815-0.28032195127130505j

print(f"Original z0: {z0}")
print(f"|P(z0)|: {abs(eval_poly(roots, z0))}")

# 1. Try just taking a second T3 step (K=2) without changing anchor
z1 = pandrosion_T3_step(roots, z0, z0)
z2 = pandrosion_T3_step(roots, z0, z1)
print(f"|P(z1)|: {abs(eval_poly(roots, z1))} (K=1)")
print(f"|P(z2)|: {abs(eval_poly(roots, z2))} (K=2)")

# 2. Try re-anchoring at z1
z2_reanchored = pandrosion_T3_step(roots, z1, z1)
print(f"|P(z2_reanchored)|: {abs(eval_poly(roots, z2_reanchored))}")

# 3. Try standard Newton step from z0
P0 = eval_poly(roots, z0)
Pprime0 = Q_diff(roots, z0, z0)
z_newton = z0 - P0 / Pprime0
print(f"|P(z_newton)|: {abs(eval_poly(roots, z_newton))}")

# 4. Try multi-start from Cauchy circle R = 2
print("\nMulti-start from Cauchy circle R = 2:")
R = 2.0
for s in range(3):
    theta = 2 * np.pi * s / 3 + np.pi/6
    z_start = R * np.exp(1j * theta)
    z = z_start
    print(f" Start {s}: {z_start}, |P| = {abs(eval_poly(roots, z))}")
    for i in range(5): # 5 epochs
        z_next = pandrosion_T3_step(roots, z, z)
        z = z_next
        print(f"   Epoch {i+1} |P|: {abs(eval_poly(roots, z))}")
